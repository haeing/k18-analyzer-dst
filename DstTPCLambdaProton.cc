// -*- C++ -*-

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "Kinematics.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCVertex.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>
#include <UnpackerManager.hh>

#define RawHit 0
#define RawCluster 0
#define TrackSearch 1
#define TrackCluster 1
#define TruncatedMean 0
#define CalibHist 0

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;
  const auto& gUnpacker = GUnpacker::get_instance();
  auto&       gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const Double_t TRUNCATED_MEAN_RATIO = 0.8; // keep the lowest 80% of per-hit dE/dx
  constexpr Double_t kInvalidBeamMomentum = -9999.;
  constexpr Double_t kBeamMomentumLoss = 0.010; // GeV/c
  const std::vector<TString> kUserParamKeys = {};

  constexpr Double_t kTightCloseDist = 5.;
  constexpr Double_t kTightChiSqr = 2.;
  constexpr Double_t kPiMinusMomCut = 0.3;
  constexpr Int_t kSigmaPion = 0x1;
  constexpr Int_t kSigmaProton = 0x4;
  constexpr Double_t kMassPi = 0.13957039;
  constexpr Double_t kMassP = 0.93827208;
  constexpr Double_t kMassLambda = 1.115683;

  Bool_t HasPIDBit(Int_t pid, Int_t bit) { return (pid & bit) != 0; }

  Bool_t IsProtonCandidate(Int_t charge, Double_t mom0, Double_t dedx)
  {
    const Int_t pid = Kinematics::HypTPCdEdxPID(dedx, charge * mom0);
    return HasPIDBit(pid, kSigmaProton) && charge > 0 && dedx >= 25.;
  }

  Bool_t IsPionCandidate(Int_t charge, Double_t mom0, Double_t dedx)
  {
    const Int_t pid = Kinematics::HypTPCdEdxPID(dedx, charge * mom0);
    return HasPIDBit(pid, kSigmaPion) && charge < 0;
  }

  Bool_t IsTightLambdaCandidate(Double_t close_dist,
                                Double_t chisqr_p,
                                Double_t chisqr_pi,
                                Double_t mom0_pi)
  {
    return std::isfinite(close_dist) &&
      close_dist < kTightCloseDist &&
      chisqr_p <= kTightChiSqr &&
      chisqr_pi <= kTightChiSqr &&
      mom0_pi < kPiMinusMomCut;
  }

  TVector3 ClosestHelixPoint(const TVector3& point, Double_t par[5],
                           Double_t tmin, Double_t tmax)
  {
    constexpr Int_t kScanSteps = 720;
    TVector3 closest(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
    Double_t best_dist2 = std::numeric_limits<Double_t>::max();
    for (Int_t i=0; i<=kScanSteps; ++i) {
      const Double_t t = tmin + (tmax-tmin)*i/kScanSteps;
      const TVector3 pos = Kinematics::CalcHelixPosition(par, t);
      const Double_t dist2 = (pos-point).Mag2();
      if (dist2 < best_dist2) { best_dist2 = dist2; closest = pos; }
    }
    return closest;
  }

  TLorentzVector
  MissingFromLambda(const TVector3& lambda_mom, Double_t lambda_mass,
                    Double_t beam_px, Double_t beam_py, Double_t beam_pz)
  {
    TLorentzVector beam;
    beam.SetXYZM(beam_px, beam_py, beam_pz, pdg::KaonMass());

    TLorentzVector target;
    target.SetXYZM(0., 0., 0., pdg::ProtonMass());

    TLorentzVector lambda;
    lambda.SetXYZM(lambda_mom.X(), lambda_mom.Y(), lambda_mom.Z(), kMassLambda);

    return beam + target - lambda;
  }

  void SanitizeBeamMomentum(Double_t& px, Double_t& py, Double_t& pz)
  {
    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        px <= kInvalidBeamMomentum + 0.5 ||
        py <= kInvalidBeamMomentum + 0.5 ||
        pz <= kInvalidBeamMomentum + 0.5) {
      px = 0.;
      py = 0.;
      pz = 0.;
    }
  }

  void ApplyBeamMaterialLoss(Double_t& px, Double_t& py, Double_t& pz)
  {
    const Double_t p = TMath::Sqrt(px*px + py*py + pz*pz);
    if (!std::isfinite(p) || p <= kBeamMomentumLoss) {
      px = py = pz = 0.;
      return;
    }
    const Double_t scale = (p - kBeamMomentumLoss)/p;
    px *= scale;
    py *= scale;
    pz *= scale;
  }

  Bool_t CopyBestK18BeamMomentum(const std::vector<Double_t>& px,
                                 const std::vector<Double_t>& py,
                                 const std::vector<Double_t>& pz,
                                 Double_t& beam_px,
                                 Double_t& beam_py,
                                 Double_t& beam_pz)
  {
    if (px.empty() || py.empty() || pz.empty()) return false;
    if (px.size() != py.size() || px.size() != pz.size()) return false;

    beam_px = px.at(0);
    beam_py = py.at(0);
    beam_pz = pz.at(0);
    return true;
  }

  // The K18 tree is ordered by its best matrix-fit candidate.  Pick the RK
  // momentum of that same candidate at the VP nearest to the Lambda vertex z.
  Bool_t CopyK18BeamMomentumAtVertexZ(
    const std::vector<std::vector<Double_t>>& zvp,
    const std::vector<std::vector<Double_t>>& pxvp,
    const std::vector<std::vector<Double_t>>& pyvp,
    const std::vector<std::vector<Double_t>>& pzvp,
    Double_t vertex_z,
    Double_t& beam_px, Double_t& beam_py, Double_t& beam_pz,
    Double_t& selected_vp_z, Int_t& selected_vp)
  {
    selected_vp_z = TMath::QuietNaN();
    selected_vp = -1;
    if (!std::isfinite(vertex_z) || zvp.empty() || pxvp.empty() ||
        pyvp.empty() || pzvp.empty()) return false;

    const auto& zs = zvp.front();
    const auto& pxs = pxvp.front();
    const auto& pys = pyvp.front();
    const auto& pzs = pzvp.front();
    if (zs.empty() || zs.size() != pxs.size() || zs.size() != pys.size() ||
        zs.size() != pzs.size()) return false;

    Double_t min_dz = std::numeric_limits<Double_t>::max();
    for (std::size_t ivp = 0; ivp < zs.size(); ++ivp) {
      if (!std::isfinite(zs[ivp]) || !std::isfinite(pxs[ivp]) ||
          !std::isfinite(pys[ivp]) || !std::isfinite(pzs[ivp])) continue;
      const Double_t dz = TMath::Abs(zs[ivp] - vertex_z);
      if (dz >= min_dz) continue;
      min_dz = dz;
      beam_px = pxs[ivp];
      beam_py = pys[ivp];
      beam_pz = pzs[ivp];
      selected_vp_z = zs[ivp];
      selected_vp = static_cast<Int_t>(ivp);
    }
    return selected_vp >= 0;
  }

}


namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpcHit, kBlc, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[DstTPCK18HelixTracking]", "[K18Tracking]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "k18", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
  Bool_t SetupReaders();
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;
  Int_t beamflag;
  Double_t beam_px;
  Double_t beam_py;
  Double_t beam_pz;
  std::vector<Double_t> clkTpc;

  Int_t ntTpc; // Number of Tracks
  Int_t effective_ntTpc; // Number of tracks with no beam/accidental/k18 flag
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> is_beam; // 1 = beam-tagged, 0 = not
  std::vector<Int_t> is_k18; // 1 = K18-tagged, 0 = not
  std::vector<Int_t> is_accidental; // 1 = accidental-tagged, 0 = not
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  // Track-wise helix theta range (TPCLocalTrackHelix::GetMint / GetMaxt); per-hit angles stay in helix_t.
  std::vector<Double_t> helix_theta_min;
  std::vector<Double_t> helix_theta_max;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx

  std::vector<Double_t> dz_factor;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Int_t> charge;//Helix charge
  std::vector<Int_t> pid;//HypTPC dE/dx PID bit pattern (same convention as TPCLocalTrackHelix::GetPid)
  std::vector<Double_t> path;//Helix path
  std::vector<std::vector<Double_t>> combi_id;
  std::vector<std::vector<Double_t>> closeDistTpc;
  std::vector<std::vector<Double_t>> vtxTpc;
  std::vector<std::vector<Double_t>> vtyTpc;
  std::vector<std::vector<Double_t>> vtzTpc;
  std::vector<std::vector<Double_t>> mom_vtx;
  std::vector<std::vector<Double_t>> mom_vty;
  std::vector<std::vector<Double_t>> mom_vtz;

  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> pathhit_cor;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  std::vector<Double_t> lambda_mass;
  std::vector<Double_t> lambda_close_dist;
  std::vector<Double_t> lambda_vtx_x;
  std::vector<Double_t> lambda_vtx_y;
  std::vector<Double_t> lambda_vtx_z;
  std::vector<Double_t> lambda_mom_x;
  std::vector<Double_t> lambda_mom_y;
  std::vector<Double_t> lambda_mom_z;
  std::vector<Double_t> lambda_proton_mom_x;
  std::vector<Double_t> lambda_proton_mom_y;
  std::vector<Double_t> lambda_proton_mom_z;
  std::vector<Double_t> lambda_pion_mom_x;
  std::vector<Double_t> lambda_pion_mom_y;
  std::vector<Double_t> lambda_pion_mom_z;

  std::vector<Int_t> lambdap_lambda_id, lambdap_proton_id;
  std::vector<Double_t> lambdap_mass, lambdap_close_dist;
  std::vector<Double_t> lambdap_vtx_x, lambdap_vtx_y, lambdap_vtx_z;
  std::vector<Double_t> lambdap_proton_vtx_x, lambdap_proton_vtx_y, lambdap_proton_vtx_z;
  std::vector<Double_t> lambdap_lambda_mom_x, lambdap_lambda_mom_y, lambdap_lambda_mom_z;
  std::vector<Double_t> lambdap_proton_mom_x, lambdap_proton_mom_y, lambdap_proton_mom_z;
  
  std::vector<Double_t> X_mass;
  std::vector<Double_t> X_mass2;
  std::vector<Double_t> X_vtx_x;
  std::vector<Double_t> X_vtx_y;
  std::vector<Double_t> X_vtx_z;
  std::vector<Double_t> X_mom_x;
  std::vector<Double_t> X_mom_y;
  std::vector<Double_t> X_mom_z;
  std::vector<Double_t> X_beam_px;
  std::vector<Double_t> X_beam_py;
  std::vector<Double_t> X_beam_pz;
  std::vector<Double_t> X_beam_vp_z;
  std::vector<Int_t> X_beam_vp_id;
  
  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    beam_px = 0.;
    beam_py = 0.;
    beam_pz = 0.;
    dst::clear_all(trigpat, trigflag, clkTpc);
  }


  void clearHelixTracks() {
    ntTpc = 0;
    effective_ntTpc = 0;
    dst::clear_all(
      nhtrack, is_beam, is_k18, is_accidental, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz,
      helix_theta_min, helix_theta_max, dE, dEdx,
      
#if TruncatedMean
      dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, 
      dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60,
#endif
      dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, pid, path,
      combi_id, closeDistTpc, vtxTpc, vtyTpc, vtzTpc, mom_vtx, mom_vty, mom_vtz,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center
    );
  }

  void clearSelectedLambdaProton() {
    dst::clear_all(
      lambda_mass, lambda_close_dist,
      lambda_vtx_x, lambda_vtx_y, lambda_vtx_z,
      lambda_mom_x, lambda_mom_y, lambda_mom_z,
      lambda_proton_mom_x, lambda_proton_mom_y, lambda_proton_mom_z,
      lambda_pion_mom_x, lambda_pion_mom_y, lambda_pion_mom_z,
      lambdap_lambda_id, lambdap_proton_id, lambdap_mass, lambdap_close_dist,
      lambdap_vtx_x, lambdap_vtx_y, lambdap_vtx_z,
      lambdap_proton_vtx_x, lambdap_proton_vtx_y, lambdap_proton_vtx_z,
      lambdap_lambda_mom_x, lambdap_lambda_mom_y, lambdap_lambda_mom_z,
      lambdap_proton_mom_x, lambdap_proton_mom_y, lambdap_proton_mom_z,
      X_mass, X_mass2,
      X_vtx_x, X_vtx_y, X_vtx_z, X_mom_x, X_mom_y, X_mom_z,
      X_beam_px, X_beam_py, X_beam_pz, X_beam_vp_z, X_beam_vp_id
    );
  }

  
  void clear()
  {
    clearBasicInfo();
#if RawHit
    clearRawHits();
#endif
#if RawCluster
    clearClusters();
#endif
    clearHelixTracks();
    clearSelectedLambdaProton();
  }

  void resizeTracks(Int_t nTracks) {
    dst::resize_all(nTracks,
      nhtrack, is_beam, is_k18, is_accidental, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz,
      helix_theta_min, helix_theta_max, dE, dEdx,

#if TruncatedMean
      dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, 
      dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60,
#endif
      dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, pid, path,
      combi_id, closeDistTpc, vtxTpc, vtyTpc, vtzTpc, mom_vtx, mom_vty, mom_vtz,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center
    );
  }

  void resizeTrackHits(Int_t it, Int_t nh) {
    dst::resize_all(nh,
      hitlayer[it], hitpos_x[it], hitpos_y[it], hitpos_z[it], calpos_x[it], calpos_y[it], calpos_z[it], residual[it], residual_x[it], residual_y[it], residual_z[it], helix_t[it],

      pathhit[it], pathhit_cor[it], theta_diff[it], track_cluster_de[it], track_cluster_size[it], track_cluster_mrow[it], track_cluster_de_center[it], track_cluster_x_center[it], track_cluster_y_center[it], track_cluster_z_center[it], track_cluster_row_center[it]
    );
  }

  void resizeTrackCombi(Int_t it, Int_t n) {
    dst::resize_all(n,
      combi_id[it], closeDistTpc[it], 
      vtxTpc[it], vtyTpc[it], vtzTpc[it], 
      mom_vtx[it], mom_vty[it], mom_vtz[it]
    );
  }

  
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;
  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<Int_t>* effective_ntTpc;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack;
  TTreeReaderValue<std::vector<Int_t>>* is_beam;
  TTreeReaderValue<std::vector<Int_t>>* is_k18;
  TTreeReaderValue<std::vector<Int_t>>* is_accidental;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* helix_theta_min;
  TTreeReaderValue<std::vector<Double_t>>* helix_theta_max;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;
  TTreeReaderValue<std::vector<Double_t>>* dz_factor;
  TTreeReaderValue<std::vector<Double_t>>* mom0_x;
  TTreeReaderValue<std::vector<Double_t>>* mom0_y;
  TTreeReaderValue<std::vector<Double_t>>* mom0_z;
  TTreeReaderValue<std::vector<Double_t>>* mom0;
  TTreeReaderValue<std::vector<Int_t>>* charge;
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* path;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* combi_id;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* closeDistTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtyTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtzTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_vtx;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_vty;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_vtz;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* helix_t;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* theta_diff;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit_cor;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<UInt_t>* blc_runnum;
  TTreeReaderValue<UInt_t>* blc_evnum;
  TTreeReaderValue<std::vector<Double_t>>* blc_px;
  TTreeReaderValue<std::vector<Double_t>>* blc_py;
  TTreeReaderValue<std::vector<Double_t>>* blc_pz;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_zvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_pxvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_pyvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_pzvpHS;
};

namespace root
{
  Event  event;
  Src    src;
  TTree *tree;
}

//_____________________________________________________________________________
namespace
{
  
#if RawHit
  //_____________________________________________________________________________  
  void FillRawHits(TPCAnalyzer& tpc_ana)
  {
    Int_t nh_tpc = 0;
    for (Int_t layer=0; layer<NumOfLayersTPC; ++layer) {
      auto hc = tpc_ana.GetTPCHC(layer);
      for (const auto& hit : hc) {
        if (!hit || !hit->IsGood())
          continue;
        const auto& pos = hit->GetPosition();
        event.raw_hitpos_x.push_back(pos.X());
        event.raw_hitpos_y.push_back(pos.Y());
        event.raw_hitpos_z.push_back(pos.Z());
        event.raw_de.push_back(hit->GetCDe());
        event.raw_padid.push_back(hit->GetPad());
        event.raw_layer.push_back(layer);
        event.raw_row.push_back(hit->GetRow());
        ++nh_tpc;
      }
    }
    event.nhTpc = nh_tpc;
  }
#endif

#if RawCluster
  //_____________________________________________________________________________  
  void FillClusters(TPCAnalyzer& tpc_ana)
  {
    Int_t n_cl_tpc = 0;
    for (Int_t layer=0; layer<NumOfLayersTPC; ++layer) {
      auto hc = tpc_ana.GetTPCClCont(layer);
      for (const auto& cl : hc) {
        if (!cl || !cl->IsGood())
          continue;
        TPCHit* center_hit = cl->GetCenterHit();
        event.cluster_x.push_back(cl->GetX());
        event.cluster_y.push_back(cl->GetY());
        event.cluster_z.push_back(cl->GetZ());
        event.cluster_de.push_back(cl->GetDe());
        event.cluster_size.push_back(cl->GetClusterSize());
        event.cluster_layer.push_back(layer);
        event.cluster_mrow.push_back(cl->MeanRow());
        event.cluster_houghflag.push_back(cl->GetHoughFlag());
        if (center_hit) {
          const TVector3& center_pos = center_hit->GetPosition();
          event.cluster_de_center.push_back(center_hit->GetCDe());
          event.cluster_x_center.push_back(center_pos.X());
          event.cluster_y_center.push_back(center_pos.Y());
          event.cluster_z_center.push_back(center_pos.Z());
          event.cluster_row_center.push_back(center_hit->GetRow());
        } else {
          event.cluster_de_center.push_back(TMath::QuietNaN());
          event.cluster_x_center.push_back(TMath::QuietNaN());
          event.cluster_y_center.push_back(TMath::QuietNaN());
          event.cluster_z_center.push_back(TMath::QuietNaN());
          event.cluster_row_center.push_back(-1);
        }
        ++n_cl_tpc;
      }
    }
    event.nclTpc = n_cl_tpc;
  }
#endif

  //_____________________________________________________________________________  
  void FillHelixPairKinematics(Int_t it, Int_t nt_tpc, TPCAnalyzer& tpc_ana)
  {
    const Double_t qnan = TMath::QuietNaN();
    for (Int_t it_pair = 0; it_pair < nt_tpc; ++it_pair) {
      event.combi_id[it][it_pair] = static_cast<Double_t>(it_pair);
      event.closeDistTpc[it][it_pair] = qnan;
      event.vtxTpc[it][it_pair] = qnan;
      event.vtyTpc[it][it_pair] = qnan;
      event.vtzTpc[it][it_pair] = qnan;
      event.mom_vtx[it][it_pair] = qnan;
      event.mom_vty[it][it_pair] = qnan;
      event.mom_vtz[it][it_pair] = qnan;
    }

    for (Int_t it_pair=0; it_pair<nt_tpc; ++it_pair) {
      if (it_pair==it) continue;
      TPCVertex* vertex = tpc_ana.FindVertexTPC(it, it_pair);
      if (!vertex || !vertex->IsCalculated())
        continue;

      Int_t track_index = -1;
      for (Int_t i = 0; i < vertex->GetNTracks(); ++i) {
        if (vertex->GetTrackId(i) == it) {
          track_index = i;
          break;
        }
      }
      if (track_index < 0)
        continue;

      const TVector3 vertex_pos = vertex->GetVertex();
      const TVector3 mom_vtx_vec = vertex->GetTrackMom(track_index);
      event.closeDistTpc[it][it_pair] = vertex->GetClosestDist();
      event.vtxTpc[it][it_pair] = vertex_pos.x();
      event.vtyTpc[it][it_pair] = vertex_pos.y();
      event.vtzTpc[it][it_pair] = vertex_pos.z();
      event.mom_vtx[it][it_pair] = mom_vtx_vec.x();
      event.mom_vty[it][it_pair] = mom_vtx_vec.y();
      event.mom_vtz[it][it_pair] = mom_vtx_vec.z();
    }
  }

  //_____________________________________________________________________________  
  void ProcessOneHelixTrack(Int_t it, Int_t nt_tpc, TPCLocalTrackHelix* helix_track,
                            TPCAnalyzer& tpc_ana,
                            TPCEventAnalyzer& event_ana)
  {
    if (!helix_track) return;

    // Per-track basic parameters and containers.
    Int_t n_hits = helix_track->GetNHit();
    Double_t chi_sqr = helix_track->GetChiSquare();
    Double_t helix_cx = helix_track->Getcx(), helix_cy = helix_track->Getcy();
    Double_t helix_z0 = helix_track->Getz0(), helix_r = helix_track->Getr();
    Double_t helix_dz = helix_track->Getdz();
    TVector3 mom0_vec = helix_track->GetMom0();
    const Int_t is_beam       = helix_track->GetIsBeam();
    const Int_t is_accidental = helix_track->GetIsAccidental();
    const Int_t is_k18        = helix_track->GetIsK18();
    const Bool_t is_no_flag   = (is_beam == 0 && is_accidental == 0 && is_k18 == 0);

    event.nhtrack[it] = n_hits;
    event.is_beam[it] = is_beam;
    event.is_accidental[it] = is_accidental;
    if (is_no_flag) {
      ++event.effective_ntTpc;
    }
    event.chisqr[it] = chi_sqr;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r;
    event.helix_dz[it] = helix_dz;
    event.dz_factor[it] = TMath::Hypot(1., helix_dz);
    if (helix_track->IsThetaCalculated()) {
      event.helix_theta_min[it] = helix_track->GetMint();
      event.helix_theta_max[it] = helix_track->GetMaxt();
    } else {
      const Double_t qnan = TMath::QuietNaN();
      event.helix_theta_min[it] = qnan;
      event.helix_theta_max[it] = qnan;
    }
    event.mom0_x[it] = mom0_vec.x();
    event.mom0_y[it] = mom0_vec.y();
    event.mom0_z[it] = mom0_vec.z();
    event.mom0[it] = mom0_vec.Mag();
    const Int_t track_pid = helix_track->GetPid();
    event.pid[it] = track_pid;
    event.resizeTrackHits(it, n_hits);
    event.resizeTrackCombi(it, nt_tpc);

    // Pair-wise helix kinematics (closest approach / vertex / momentum at vertex).
    FillHelixPairKinematics(it, nt_tpc, tpc_ana);

    // Per-hit observables and cluster-linked quantities.
#if TruncatedMean
    std::vector<Double_t> dedx_vec;
    std::vector<Double_t> dedx_cor_vec;
#endif
    for (Int_t ih = 0; ih < n_hits; ++ih) {
      TPCLTrackHit* hit = helix_track->GetHit(ih);
      if (!hit) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hit_pos = hit->GetLocalHitPos();
      const TVector3& cal_pos = hit->GetLocalCalPosHelix();
      const TVector3& residual_vec = hit->GetResidualVect();
#if TrackCluster
      event_ana.FillHelixHitHist(hit, (nt_tpc == 1 && n_hits >= 20), track_pid);
#else
      event_ana.FillHelixHitHist(hit, false, track_pid);
#endif

      TPCHit* cl_hit = hit->GetHit();
      TPCCluster* cl = cl_hit->GetParentCluster();
      Int_t cl_size = cl->GetClusterSize();
      Double_t cl_de = cl->GetDe();
      Double_t mrow = cl->MeanRow();
      TPCHit* center_hit = cl->GetCenterHit();
      const TVector3 center_pos =
        center_hit ? center_hit->GetPosition() : cl->GetPosition();
      Double_t center_de = center_hit ? center_hit->GetCDe() : TMath::QuietNaN();
      Int_t center_row = center_hit ? center_hit->GetRow() : -1;

      event.track_cluster_de[it][ih] = cl_de;
      event.track_cluster_size[it][ih] = cl_size;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = center_de;
      event.track_cluster_x_center[it][ih] = center_pos.X();
      event.track_cluster_y_center[it][ih] = center_pos.Y();
      event.track_cluster_z_center[it][ih] = center_pos.Z();
      event.track_cluster_row_center[it][ih] = center_row;

      Double_t pad_theta = tpc::GetTheta(layer, mrow)*TMath::DegToRad();
      Double_t t_cal = hit->GetTheta();
      Double_t theta_diff = t_cal - pad_theta;
      event.theta_diff[it][ih] = theta_diff;
      Double_t path_hit = tpc::padParameter[layer][5];
      event.pathhit[it][ih] = path_hit;
      Double_t path_hit_cor = (path_hit/TMath::Abs(TMath::Cos(theta_diff)))*TMath::Hypot(1., helix_dz);
      event.pathhit_cor[it][ih] = path_hit_cor;

#if TruncatedMean
      Double_t dedx_step = cl_de/path_hit;
      Double_t dedx_cor_step = cl_de/path_hit_cor;
      dedx_vec.push_back(dedx_step);
      dedx_cor_vec.push_back(dedx_cor_step);
#endif

      Double_t residual = hit->GetResidual();
      event.hitlayer[it][ih] = static_cast<Double_t>(layer);
      event.hitpos_x[it][ih] = hit_pos.x();
      event.hitpos_y[it][ih] = hit_pos.y();
      event.hitpos_z[it][ih] = hit_pos.z();
      event.calpos_x[it][ih] = cal_pos.x();
      event.calpos_y[it][ih] = cal_pos.y();
      event.calpos_z[it][ih] = cal_pos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = residual_vec.x();
      event.residual_y[it][ih] = residual_vec.y();
      event.residual_z[it][ih] = residual_vec.z();
      event.helix_t[it][ih] = t_cal;
    }

    // Track-level integrals and dE/dx summary.
    event.charge[it] = helix_track->GetCharge();
    event.path[it] = helix_track->GetPath();
    event.dE[it] = helix_track->GetTrackdE();
    event.dEdx[it] = helix_track->GetdEdx(TRUNCATED_MEAN_RATIO);

#if TruncatedMean
    std::sort(dedx_vec.begin(), dedx_vec.end());
    std::sort(dedx_cor_vec.begin(), dedx_cor_vec.end());
    std::vector<Double_t> dedx_cumulative(dedx_vec.size() + 1, 0.);
    std::vector<Double_t> dedx_cor_cumulative(dedx_cor_vec.size() + 1, 0.);
    for (std::size_t i = 0; i < dedx_vec.size(); ++i) {
      dedx_cumulative[i + 1] = dedx_cumulative[i] + dedx_vec[i];
      dedx_cor_cumulative[i + 1] = dedx_cor_cumulative[i] + dedx_cor_vec[i];
    }

    event.dEdx_0[it] = TMath::Mean(dedx_vec.size(), dedx_vec.data());
    event.dEdx_cor_0[it] = TMath::Mean(dedx_cor_vec.size(), dedx_cor_vec.data());
    event.dEdx_10[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.9);
    event.dEdx_20[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.8);
    event.dEdx_30[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.7);
    event.dEdx_40[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.6);
    event.dEdx_50[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.5);
    event.dEdx_60[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cumulative, 0.4);
    event.dEdx_cor_10[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.9);
    event.dEdx_cor_20[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.8);
    event.dEdx_cor_30[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.7);
    event.dEdx_cor_40[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.6);
    event.dEdx_cor_50[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.5);
    event.dEdx_cor_60[it] = TPCEventAnalyzer::CalcTruncatedMean(dedx_cor_cumulative, 0.4);
#endif
    HF1("TPCTrk_Num_TrackHits", n_hits);
    HF1("TPCTrk_Chisqr", chi_sqr);
    event_ana.FillHelixPidHist(helix_track, track_pid, event.dEdx[it]);
  }


} // namespace

//_____________________________________________________________________________
Int_t
main(Int_t argc, char** argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if (!CheckArg(arg))
    return EXIT_FAILURE;
  if (!DstOpen(arg))
    return EXIT_FAILURE;
  if (!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if (!dst::ValidateUserParams(gUser, kUserParamKeys))
    return EXIT_FAILURE;
  if (!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if (!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if (!dst::SetupReaders())
    return EXIT_FAILURE;
  
  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();
  Int_t ievent = skip;
  for (; ievent<nevent && !CatchSignal::Stop(); ++ievent) {
    gCounter.check();
    InitializeEvent();
    if (DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t n_input_files = 0;
  for (const auto& name : TreeName) if (name != "") n_input_files++;

  Int_t open_file = 0, open_tree = 0;
  for (Int_t i=0; i<nArgc; ++i) {
    if (TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if (open_file!=n_input_files || open_tree!=n_input_files) {
    spdlog::error("DstOpen Failed: opened files/trees mismatch based on TreeName definitions."
                  " expected: {}, open_file: {}, open_tree: {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  std::vector<TTree*> aligned_trees = { TTreeCont[kTpcHit] };
  if (!CheckEntries(aligned_trees)) return false;
  if (TTreeCont[kBlc]->BuildIndex("run_number", "event_number") < 0) {
    spdlog::error("failed to build k18 index with run_number/event_number");
    return false;
  }

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if (ievent%1==0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;

  const Long64_t blc_entry = TTreeCont[kBlc]->GetEntryNumberWithIndex(event.runnum, event.evnum);
  if (blc_entry >= 0) {
    TTreeReaderCont[kBlc]->SetEntry(blc_entry);
    if (**src.blc_runnum == event.runnum && **src.blc_evnum == event.evnum) {
      CopyBestK18BeamMomentum(**src.blc_px, **src.blc_py, **src.blc_pz,
                              event.beam_px, event.beam_py, event.beam_pz);
    }
  }
  SanitizeBeamMomentum(event.beam_px, event.beam_py, event.beam_pz);
  HF1("Status", event.status++);

  event.ntTpc = **src.ntTpc;
  event.effective_ntTpc = **src.effective_ntTpc;
  event.nhtrack = **src.nhtrack;
  event.is_beam = **src.is_beam;
  event.is_k18 = **src.is_k18;
  event.is_accidental = **src.is_accidental;
  event.chisqr = **src.chisqr;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;
  event.helix_theta_min = **src.helix_theta_min;
  event.helix_theta_max = **src.helix_theta_max;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.dz_factor = **src.dz_factor;
  event.mom0_x = **src.mom0_x;
  event.mom0_y = **src.mom0_y;
  event.mom0_z = **src.mom0_z;
  event.mom0 = **src.mom0;
  event.charge = **src.charge;
  event.pid = **src.pid;
  event.path = **src.path;
  event.combi_id = **src.combi_id;
  event.closeDistTpc = **src.closeDistTpc;
  event.vtxTpc = **src.vtxTpc;
  event.vtyTpc = **src.vtyTpc;
  event.vtzTpc = **src.vtzTpc;
  event.mom_vtx = **src.mom_vtx;
  event.mom_vty = **src.mom_vty;
  event.mom_vtz = **src.mom_vtz;
  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.helix_t = **src.helix_t;
  event.theta_diff = **src.theta_diff;
  event.pathhit = **src.pathhit;
  event.pathhit_cor = **src.pathhit_cor;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;
  HF1("NTracks_TPC", event.ntTpc);
  HF1("Status", event.status++);

  for (Int_t ip = 0; ip < event.ntTpc; ++ip) {
    if (event.is_beam.at(ip) || event.is_k18.at(ip) || event.is_accidental.at(ip)) continue;
    const Int_t charge_p = event.charge.at(ip);
    const Double_t mom0_p = event.mom0.at(ip);
    const Double_t dedx_p = event.dEdx.at(ip);
    if (!IsProtonCandidate(charge_p, mom0_p, dedx_p)) continue;

    for (Int_t ipi = 0; ipi < event.ntTpc; ++ipi) {
      if (ip == ipi) continue;
      if (event.is_beam.at(ipi) || event.is_k18.at(ipi) || event.is_accidental.at(ipi)) continue;

      const Int_t charge_pi = event.charge.at(ipi);
      const Double_t mom0_pi = event.mom0.at(ipi);
      const Double_t dedx_pi = event.dEdx.at(ipi);
      if (!IsPionCandidate(charge_pi, mom0_pi, dedx_pi)) continue;

      const Double_t close_dist = event.closeDistTpc.at(ip).at(ipi);
      const Double_t vtx_x = event.vtxTpc.at(ip).at(ipi);
      const Double_t vtx_y = event.vtyTpc.at(ip).at(ipi);
      const Double_t vtx_z = event.vtzTpc.at(ip).at(ipi);
      if (!std::isfinite(close_dist) ||
          !std::isfinite(vtx_x) || !std::isfinite(vtx_y) || !std::isfinite(vtx_z)) continue;

      const Double_t p_px = event.mom_vtx.at(ip).at(ipi);
      const Double_t p_py = event.mom_vty.at(ip).at(ipi);
      const Double_t p_pz = event.mom_vtz.at(ip).at(ipi);
      const Double_t pi_px = event.mom_vtx.at(ipi).at(ip);
      const Double_t pi_py = event.mom_vty.at(ipi).at(ip);
      const Double_t pi_pz = event.mom_vtz.at(ipi).at(ip);
      if (!std::isfinite(p_px) || !std::isfinite(p_py) || !std::isfinite(p_pz) ||
          !std::isfinite(pi_px) || !std::isfinite(pi_py) || !std::isfinite(pi_pz)) continue;

      const Double_t chisqr_p = event.chisqr.at(ip);
      const Double_t chisqr_pi = event.chisqr.at(ipi);
      if (!IsTightLambdaCandidate(close_dist, chisqr_p, chisqr_pi, mom0_pi)) continue;

      const TVector3 p_mom(p_px, p_py, p_pz);
      const TVector3 pi_mom(pi_px, pi_py, pi_pz);

      TLorentzVector lv_p;
      TLorentzVector lv_pi;
      lv_p.SetXYZM(p_mom.X(), p_mom.Y(), p_mom.Z(), kMassP);
      lv_pi.SetXYZM(pi_mom.X(), pi_mom.Y(), pi_mom.Z(), kMassPi);
      const TLorentzVector lambda = lv_p + lv_pi;
      const TVector3 lambda_mom = lambda.Vect();
      // The Lambda daughters are evaluated at vtx_z above, so evaluate the
      // beam at the same z with the K18 RK trajectory.  Old K18 files without
      // a valid VP trajectory retain the VO momentum as a fallback.
      Double_t beam_px = event.beam_px;
      Double_t beam_py = event.beam_py;
      Double_t beam_pz = event.beam_pz;
      Double_t beam_vp_z = TMath::QuietNaN();
      Int_t beam_vp_id = -1;
      CopyK18BeamMomentumAtVertexZ(**src.blc_zvpHS,
                                   **src.blc_pxvpHS,
                                   **src.blc_pyvpHS,
                                   **src.blc_pzvpHS,
                                   vtx_z,
                                   beam_px, beam_py, beam_pz,
                                   beam_vp_z, beam_vp_id);
      Double_t beam_px_for_missing = beam_px;
      Double_t beam_py_for_missing = beam_py;
      Double_t beam_pz_for_missing = beam_pz;
      ApplyBeamMaterialLoss(beam_px_for_missing, beam_py_for_missing,
                            beam_pz_for_missing);
      const TLorentzVector missing = MissingFromLambda(lambda_mom, lambda.M(),
                                                     beam_px_for_missing,
                                                     beam_py_for_missing,
                                                     beam_pz_for_missing);
      const TVector3 missing_mom = missing.Vect();

      HF1("LambdaEta_Mass", lambda.M());
      HF1("LambdaEta_CloseDist", close_dist);
      HF1("LambdaEta_VtxX", vtx_x);
      HF1("LambdaEta_VtxY", vtx_y);
      HF1("LambdaEta_VtxZ", vtx_z);
      HF1("LambdaEta_Proton_MomX", p_mom.X());
      HF1("LambdaEta_Proton_MomY", p_mom.Y());
      HF1("LambdaEta_Proton_MomZ", p_mom.Z());
      HF1("LambdaEta_Proton_Mom", lv_p.P());
      HF1("LambdaEta_Pion_MomX", pi_mom.X());
      HF1("LambdaEta_Pion_MomY", pi_mom.Y());
      HF1("LambdaEta_Pion_MomZ", pi_mom.Z());
      HF1("LambdaEta_Pion_Mom", lv_pi.P());
      HF1("LambdaEta_Lambda_MomX", lambda_mom.X());
      HF1("LambdaEta_Lambda_MomY", lambda_mom.Y());
      HF1("LambdaEta_Lambda_MomZ", lambda_mom.Z());
      HF1("LambdaEta_Lambda_Mom", lambda_mom.Mag());

      event.lambda_mass.push_back(lambda.M());
      event.lambda_close_dist.push_back(close_dist);
      event.lambda_vtx_x.push_back(vtx_x);
      event.lambda_vtx_y.push_back(vtx_y);
      event.lambda_vtx_z.push_back(vtx_z);
      event.lambda_mom_x.push_back(lambda_mom.X());
      event.lambda_mom_y.push_back(lambda_mom.Y());
      event.lambda_mom_z.push_back(lambda_mom.Z());
      event.lambda_proton_mom_x.push_back(p_mom.X());
      event.lambda_proton_mom_y.push_back(p_mom.Y());
      event.lambda_proton_mom_z.push_back(p_mom.Z());
      event.lambda_pion_mom_x.push_back(pi_mom.X());
      event.lambda_pion_mom_y.push_back(pi_mom.Y());
      event.lambda_pion_mom_z.push_back(pi_mom.Z());
      event.X_mass.push_back(missing.M());
      event.X_mass2.push_back(missing.M2());
      event.X_vtx_x.push_back(vtx_x);
      event.X_vtx_y.push_back(vtx_y);
      event.X_vtx_z.push_back(vtx_z);
      event.X_mom_x.push_back(missing_mom.X());
      event.X_mom_y.push_back(missing_mom.Y());
      event.X_mom_z.push_back(missing_mom.Z());
      event.X_beam_px.push_back(beam_px);
      event.X_beam_py.push_back(beam_py);
      event.X_beam_pz.push_back(beam_pz);
      event.X_beam_vp_z.push_back(beam_vp_z);
      event.X_beam_vp_id.push_back(beam_vp_id);

      const Int_t lambda_id = static_cast<Int_t>(event.lambda_mass.size()) - 1;
      const TVector3 lambda_vtx(vtx_x, vtx_y, vtx_z);
      for (Int_t ip_extra=0; ip_extra<event.ntTpc; ++ip_extra) {
        if (ip_extra == ip || ip_extra == ipi) continue;
        if (event.is_beam.at(ip_extra) || event.is_k18.at(ip_extra) ||
            event.is_accidental.at(ip_extra)) continue;
        if (!IsProtonCandidate(event.charge.at(ip_extra), event.mom0.at(ip_extra),
                               event.dEdx.at(ip_extra))) continue;

        Double_t helix_par[5] = {event.helix_cx.at(ip_extra),
                                 event.helix_cy.at(ip_extra),
                                 event.helix_z0.at(ip_extra),
                                 event.helix_r.at(ip_extra),
                                 event.helix_dz.at(ip_extra)};
        Double_t tmin = event.helix_theta_min.at(ip_extra);
        Double_t tmax = event.helix_theta_max.at(ip_extra);
        if (tmin > tmax) std::swap(tmin, tmax);
        const Double_t vertex_dist =
          Kinematics::CalcHelixCloseDist(lambda_vtx, helix_par, tmin, tmax);
        if (!std::isfinite(vertex_dist)) continue;
        const TVector3 proton_vtx =
          ClosestHelixPoint(lambda_vtx, helix_par, tmin, tmax);
        if (!std::isfinite(proton_vtx.X()) || !std::isfinite(proton_vtx.Y()) ||
            !std::isfinite(proton_vtx.Z())) continue;

        const TVector3 extra_proton_mom(event.mom0_x.at(ip_extra),
                                        event.mom0_y.at(ip_extra),
                                        event.mom0_z.at(ip_extra));
        if (!std::isfinite(extra_proton_mom.X()) ||
            !std::isfinite(extra_proton_mom.Y()) ||
            !std::isfinite(extra_proton_mom.Z())) continue;
        TLorentzVector extra_proton;
        extra_proton.SetXYZM(extra_proton_mom.X(), extra_proton_mom.Y(),
                             extra_proton_mom.Z(), kMassP);
        const TLorentzVector lambda_proton = lambda + extra_proton;

        event.lambdap_lambda_id.push_back(lambda_id);
        event.lambdap_proton_id.push_back(ip_extra);
        event.lambdap_mass.push_back(lambda_proton.M());
        event.lambdap_close_dist.push_back(vertex_dist);
        event.lambdap_vtx_x.push_back(vtx_x);
        event.lambdap_vtx_y.push_back(vtx_y);
        event.lambdap_vtx_z.push_back(vtx_z);
        event.lambdap_proton_vtx_x.push_back(proton_vtx.X());
        event.lambdap_proton_vtx_y.push_back(proton_vtx.Y());
        event.lambdap_proton_vtx_z.push_back(proton_vtx.Z());
        event.lambdap_lambda_mom_x.push_back(lambda_mom.X());
        event.lambdap_lambda_mom_y.push_back(lambda_mom.Y());
        event.lambdap_lambda_mom_z.push_back(lambda_mom.Z());
        event.lambdap_proton_mom_x.push_back(extra_proton_mom.X());
        event.lambdap_proton_mom_y.push_back(extra_proton_mom.Y());
        event.lambdap_proton_mom_z.push_back(extra_proton_mom.Z());
        HF1("LambdaProton_Mass", lambda_proton.M());
        HF1("LambdaProton_CloseDist", vertex_dist);
        HF2("LambdaProton_LambdaVtxXY", vtx_x, vtx_y);
        HF2("LambdaProton_LambdaVtxYZ", vtx_z, vtx_y);
        HF2("LambdaProton_LambdaVtxZX", vtx_z, vtx_x);
        HF2("LambdaProton_ProtonVtxXY", proton_vtx.X(), proton_vtx.Y());
        HF2("LambdaProton_ProtonVtxYZ", proton_vtx.Z(), proton_vtx.Y());
        HF2("LambdaProton_ProtonVtxZX", proton_vtx.Z(), proton_vtx.X());
      }
    }
  }
  HF1("Status", event.status++);

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  spdlog::info(" Close : {}", TFileCont[kOutFile]->GetName());
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for (Int_t i=0; i<n; ++i) {
    if (TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if (TTreeCont[i]) delete TTreeCont[i];
    if (TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::SetupReaders()
{
  if (!dst::SetupReader(kTpcHit, "kTpcHit")) return false;
  if (!dst::SetupReader(kBlc, "kBlc")) return false;

  const auto& reader = TTreeReaderCont[kTpcHit];

  dst::SetBranch(reader, "run_number", src.runnum);
  dst::SetBranch(reader, "event_number", src.evnum);
  dst::SetBranch(reader, "trig_pat", src.trigpat);
  dst::SetBranch(reader, "trig_flag", src.trigflag);
  dst::SetBranch(reader, "beam_flag", src.beamflag);
  dst::SetBranch(reader, "clkTpc", src.clkTpc);
  dst::SetBranch(reader, "ntTpc", src.ntTpc);
  dst::SetBranch(reader, "effective_ntTpc", src.effective_ntTpc);
  dst::SetBranch(reader, "nhtrack", src.nhtrack);
  dst::SetBranch(reader, "is_beam", src.is_beam);
  dst::SetBranch(reader, "is_k18", src.is_k18);
  dst::SetBranch(reader, "is_accidental", src.is_accidental);
  dst::SetBranch(reader, "chisqr", src.chisqr);
  dst::SetBranch(reader, "helix_cx", src.helix_cx);
  dst::SetBranch(reader, "helix_cy", src.helix_cy);
  dst::SetBranch(reader, "helix_z0", src.helix_z0);
  dst::SetBranch(reader, "helix_r", src.helix_r);
  dst::SetBranch(reader, "helix_dz", src.helix_dz);
  dst::SetBranch(reader, "helix_theta_min", src.helix_theta_min);
  dst::SetBranch(reader, "helix_theta_max", src.helix_theta_max);
  dst::SetBranch(reader, "mom0_x", src.mom0_x);
  dst::SetBranch(reader, "mom0_y", src.mom0_y);
  dst::SetBranch(reader, "mom0_z", src.mom0_z);
  dst::SetBranch(reader, "mom0", src.mom0);
  dst::SetBranch(reader, "dE", src.dE);
  dst::SetBranch(reader, "dEdx", src.dEdx);
  dst::SetBranch(reader, "pid", src.pid);
  dst::SetBranch(reader, "dz_factor", src.dz_factor);
  dst::SetBranch(reader, "charge", src.charge);
  dst::SetBranch(reader, "path", src.path);
  dst::SetBranch(reader, "combi_id", src.combi_id);
  dst::SetBranch(reader, "closeDistTpc", src.closeDistTpc);
  dst::SetBranch(reader, "vtxTpc", src.vtxTpc);
  dst::SetBranch(reader, "vtyTpc", src.vtyTpc);
  dst::SetBranch(reader, "vtzTpc", src.vtzTpc);
  dst::SetBranch(reader, "mom_vtx", src.mom_vtx);
  dst::SetBranch(reader, "mom_vty", src.mom_vty);
  dst::SetBranch(reader, "mom_vtz", src.mom_vtz);
  dst::SetBranch(reader, "hitlayer", src.hitlayer);
  dst::SetBranch(reader, "hitpos_x", src.hitpos_x);
  dst::SetBranch(reader, "hitpos_y", src.hitpos_y);
  dst::SetBranch(reader, "hitpos_z", src.hitpos_z);
  dst::SetBranch(reader, "calpos_x", src.calpos_x);
  dst::SetBranch(reader, "calpos_y", src.calpos_y);
  dst::SetBranch(reader, "calpos_z", src.calpos_z);
  dst::SetBranch(reader, "residual", src.residual);
  dst::SetBranch(reader, "residual_x", src.residual_x);
  dst::SetBranch(reader, "residual_y", src.residual_y);
  dst::SetBranch(reader, "residual_z", src.residual_z);
  dst::SetBranch(reader, "helix_t", src.helix_t);
  dst::SetBranch(reader, "theta_diff", src.theta_diff);
  dst::SetBranch(reader, "pathhit", src.pathhit);
  dst::SetBranch(reader, "pathhit_cor", src.pathhit_cor);
  dst::SetBranch(reader, "track_cluster_de", src.track_cluster_de);
  dst::SetBranch(reader, "track_cluster_size", src.track_cluster_size);
  dst::SetBranch(reader, "track_cluster_mrow", src.track_cluster_mrow);
  dst::SetBranch(reader, "track_cluster_de_center", src.track_cluster_de_center);
  dst::SetBranch(reader, "track_cluster_x_center", src.track_cluster_x_center);
  dst::SetBranch(reader, "track_cluster_y_center", src.track_cluster_y_center);
  dst::SetBranch(reader, "track_cluster_z_center", src.track_cluster_z_center);
  dst::SetBranch(reader, "track_cluster_row_center", src.track_cluster_row_center);

  const auto& blc_reader = TTreeReaderCont[kBlc];
  dst::SetBranch(blc_reader, "run_number", src.blc_runnum);
  dst::SetBranch(blc_reader, "event_number", src.blc_evnum);
  dst::SetBranch(blc_reader, "px", src.blc_px);
  dst::SetBranch(blc_reader, "py", src.blc_py);
  dst::SetBranch(blc_reader, "pz", src.blc_pz);
  dst::SetBranch(blc_reader, "zvpHS", src.blc_zvpHS);
  dst::SetBranch(blc_reader, "pxvpHS", src.blc_pxvpHS);
  dst::SetBranch(blc_reader, "pyvpHS", src.blc_pyvpHS);
  dst::SetBranch(blc_reader, "pzvpHS", src.blc_pzvpHS);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  TPCEventAnalyzer::SetDstCalibFlag(
#if CalibHist
    true
#else
    false
#endif
  );

  hist::BuildStatus();
  hist::BuildTPCHelixTracking(TPCEventAnalyzer::GetDstCalibFlag());
  hist::BuildTPCLambdaEta();
  hist::BuildTPCLambdaProton();

  tree = new TTree("tpc", "tree of DstTPCLambdaProton");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("beam_px", &event.beam_px);
  tree->Branch("beam_py", &event.beam_py);
  tree->Branch("beam_pz", &event.beam_pz);
  tree->Branch("clkTpc", &event.clkTpc);

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "effective_ntTpc", &event.effective_ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "is_beam", &event.is_beam );
  tree->Branch( "is_k18", &event.is_k18 );
  tree->Branch( "is_accidental", &event.is_accidental );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "helix_theta_min", &event.helix_theta_min );
  tree->Branch( "helix_theta_max", &event.helix_theta_max );
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "pid", &event.pid );
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

  tree->Branch( "combi_id", &event.combi_id );
  tree->Branch( "closeDistTpc", &event.closeDistTpc );
  tree->Branch( "vtxTpc", &event.vtxTpc );
  tree->Branch( "vtyTpc", &event.vtyTpc );
  tree->Branch( "vtzTpc", &event.vtzTpc );
  tree->Branch( "mom_vtx", &event.mom_vtx );
  tree->Branch( "mom_vty", &event.mom_vty );
  tree->Branch( "mom_vtz", &event.mom_vtz );

  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "theta_diff", &event.theta_diff);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "pathhit_cor", &event.pathhit_cor);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  tree->Branch("lambda_mass", &event.lambda_mass);
  tree->Branch("lambda_close_dist", &event.lambda_close_dist);
  tree->Branch("lambda_vtx_x", &event.lambda_vtx_x);
  tree->Branch("lambda_vtx_y", &event.lambda_vtx_y);
  tree->Branch("lambda_vtx_z", &event.lambda_vtx_z);
  tree->Branch("lambda_mom_x", &event.lambda_mom_x);
  tree->Branch("lambda_mom_y", &event.lambda_mom_y);
  tree->Branch("lambda_mom_z", &event.lambda_mom_z);
  tree->Branch("lambda_proton_mom_x", &event.lambda_proton_mom_x);
  tree->Branch("lambda_proton_mom_y", &event.lambda_proton_mom_y);
  tree->Branch("lambda_proton_mom_z", &event.lambda_proton_mom_z);
  tree->Branch("lambda_pion_mom_x", &event.lambda_pion_mom_x);
  tree->Branch("lambda_pion_mom_y", &event.lambda_pion_mom_y);
  tree->Branch("lambda_pion_mom_z", &event.lambda_pion_mom_z);

  tree->Branch("lambdap_lambda_id", &event.lambdap_lambda_id);
  tree->Branch("lambdap_proton_id", &event.lambdap_proton_id);
  tree->Branch("lambdap_mass", &event.lambdap_mass);
  tree->Branch("lambdap_close_dist", &event.lambdap_close_dist);
  tree->Branch("lambdap_vtx_x", &event.lambdap_vtx_x);
  tree->Branch("lambdap_vtx_y", &event.lambdap_vtx_y);
  tree->Branch("lambdap_vtx_z", &event.lambdap_vtx_z);
  tree->Branch("lambdap_proton_vtx_x", &event.lambdap_proton_vtx_x);
  tree->Branch("lambdap_proton_vtx_y", &event.lambdap_proton_vtx_y);
  tree->Branch("lambdap_proton_vtx_z", &event.lambdap_proton_vtx_z);
  tree->Branch("lambdap_lambda_mom_x", &event.lambdap_lambda_mom_x);
  tree->Branch("lambdap_lambda_mom_y", &event.lambdap_lambda_mom_y);
  tree->Branch("lambdap_lambda_mom_z", &event.lambdap_lambda_mom_z);
  tree->Branch("lambdap_proton_mom_x", &event.lambdap_proton_mom_x);
  tree->Branch("lambdap_proton_mom_y", &event.lambdap_proton_mom_y);
  tree->Branch("lambdap_proton_mom_z", &event.lambdap_proton_mom_z);


  tree->Branch("X_mass", &event.X_mass);
  tree->Branch("X_mass2", &event.X_mass2);
  tree->Branch("X_vtx_x", &event.X_vtx_x);
  tree->Branch("X_vtx_y", &event.X_vtx_y);
  tree->Branch("X_vtx_z", &event.X_vtx_z);
  tree->Branch("X_mom_x", &event.X_mom_x);
  tree->Branch("X_mom_y", &event.X_mom_y);
  tree->Branch("X_mom_z", &event.X_mom_z);
  tree->Branch("X_beam_px", &event.X_beam_px);
  tree->Branch("X_beam_py", &event.X_beam_py);
  tree->Branch("X_beam_pz", &event.X_beam_pz);
  tree->Branch("X_beam_vp_z", &event.X_beam_vp_z);
  tree->Branch("X_beam_vp_id", &event.X_beam_vp_id);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
