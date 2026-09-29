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
  //constexpr Double_t kBeamMomentumLoss = 0.00579; // GeV/c
  constexpr Double_t kBeamMomentumLoss = 0.0; // GeV/c
  // Require a valid BLC/K18 RK trajectory, but do not require that its
  // VP-constrained TPC helix was reconstructed/tagged as a beam track.
  // The latter is a TPC-reconstruction efficiency cut and is unrelated to
  // the RK trajectory used for the production-vertex calculation below.
  constexpr Bool_t kRequireBeamTrack = false;
  constexpr Bool_t kRequireEffectiveNTracks = true;
  constexpr Int_t kRequiredEffectiveNTracks = 2;
  const std::vector<TString> kUserParamKeys = {};

  constexpr Double_t kTightCloseDist = 5.;
  constexpr Double_t kTightChiSqr = 2.;
  constexpr Double_t kPiMinusMomCut = 0.3;
  constexpr Int_t kSigmaPion = 0x1;
  constexpr Int_t kSigmaProton = 0x4;
  constexpr Double_t kMassPi = 0.13957039;
  constexpr Double_t kMassP = 0.93827208;
  constexpr Double_t kMassLambda = 1.115683;

  // E72 LH2 target: cylinder axis is y, centered at (0, 0, tpc::Z_TARGET).
  // Use guard bands so the two vertex samples cannot contain target-frame
  // boundary vertices.  The physical LH2/frame radii are 40/50 mm.
  constexpr Double_t kLH2InnerRadius = 38.;
  constexpr Double_t kLH2InnerHalfLengthY = 48.;
  constexpr Double_t kFrameOuterRadius = 52.;
  constexpr Double_t kFrameOuterHalfLengthY = 52.;

  Bool_t IsStrictlyInsideLH2(Double_t x, Double_t y, Double_t z)
  {
    const Double_t r_xz = TMath::Hypot(x, z - tpc::Z_TARGET);
    return r_xz < kLH2InnerRadius && TMath::Abs(y) < kLH2InnerHalfLengthY;
  }

  Bool_t IsCompletelyOutsideTargetFrame(Double_t x, Double_t y, Double_t z)
  {
    const Double_t r_xz = TMath::Hypot(x, z - tpc::Z_TARGET);
    return r_xz > kFrameOuterRadius || TMath::Abs(y) > kFrameOuterHalfLengthY;
  }

  Bool_t HasPIDBit(Int_t pid, Int_t bit) { return (pid & bit) != 0; }

  // Local-x displacement from the fitted helix over its stored theta range.
  Bool_t CalcHelixDx(Double_t radius, Double_t theta_min, Double_t theta_max,
                     Double_t& dx)
  {
    if (!std::isfinite(radius) || !std::isfinite(theta_min) ||
        !std::isfinite(theta_max)) return false;
    dx = radius * (std::cos(theta_max) - std::cos(theta_min));
    return std::isfinite(dx);
  }

  Bool_t IsProtonCandidate(Int_t charge, Double_t mom0, Double_t dedx)
  {
    const Double_t nsigma_p = Kinematics::HypTPCdEdxNsigmaProton(dedx, charge * mom0);
    return charge > 0 && dedx >= 25. && nsigma_p > -1.0;
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
  struct BeamLambdaClosestApproach {
    Bool_t found = false;
    TVector3 beam_pos;
    TVector3 lambda_pos;
    TVector3 beam_mom;
    Double_t dca = TMath::QuietNaN();
    Double_t beam_segment_fraction = TMath::QuietNaN();
    Int_t beam_segment_id = -1;
  };

  // Find the production vertex from the closest approach of the K18 RK
  // trajectory (piecewise-linear between consecutive virtual planes) and the
  // Lambda flight line.  The Lambda line is constrained upstream of its decay
  // vertex: decay_vertex + s*pLambda_unit, with s <= 0.
  BeamLambdaClosestApproach FindBeamLambdaProductionVertex(
    const std::vector<std::vector<Double_t>>& xvp,
    const std::vector<std::vector<Double_t>>& yvp,
    const std::vector<std::vector<Double_t>>& zvp,
    const std::vector<std::vector<Double_t>>& pxvp,
    const std::vector<std::vector<Double_t>>& pyvp,
    const std::vector<std::vector<Double_t>>& pzvp,
    const TVector3& lambda_decay_vertex, const TVector3& lambda_mom)
  {
    BeamLambdaClosestApproach best;
    if (xvp.empty() || yvp.empty() || zvp.empty() || pxvp.empty() || pyvp.empty() || pzvp.empty()) return best;
    const auto& xs = xvp.front(); const auto& ys = yvp.front(); const auto& zs = zvp.front();
    const auto& pxs = pxvp.front(); const auto& pys = pyvp.front(); const auto& pzs = pzvp.front();
    const auto n = std::min({xs.size(), ys.size(), zs.size(), pxs.size(), pys.size(), pzs.size()});
    if (n < 2 || !std::isfinite(lambda_mom.Mag2()) || lambda_mom.Mag2() <= 0.) return best;
    const TVector3 lambda_dir = lambda_mom.Unit();
    Double_t best_d2 = std::numeric_limits<Double_t>::max();
    for (std::size_t i = 0; i+1 < n; ++i) {
      const TVector3 b0(xs[i], ys[i], zs[i]);
      const TVector3 b1(xs[i+1], ys[i+1], zs[i+1]);
      if (!std::isfinite(b0.X()) || !std::isfinite(b0.Y()) || !std::isfinite(b0.Z()) ||
          !std::isfinite(b1.X()) || !std::isfinite(b1.Y()) || !std::isfinite(b1.Z())) continue;
      const TVector3 beam_step = b1-b0;
      const Double_t a = beam_step.Mag2();
      if (a <= 0.) continue;
      const TVector3 w = b0-lambda_decay_vertex;
      const Double_t b = beam_step.Dot(lambda_dir);
      const Double_t c = 1.;
      const Double_t d = beam_step.Dot(w);
      const Double_t e = lambda_dir.Dot(w);
      auto consider = [&](Double_t beam_t, Double_t lambda_s) {
        if (beam_t < 0. || beam_t > 1. || lambda_s > 0.) return;
        const TVector3 beam_pos = b0 + beam_t*beam_step;
        const TVector3 lambda_pos = lambda_decay_vertex + lambda_s*lambda_dir;
        const Double_t d2 = (beam_pos-lambda_pos).Mag2();
        if (!std::isfinite(d2) || d2 >= best_d2) return;
        best_d2 = d2;
        best.found = true;
        best.beam_pos = beam_pos;
        best.lambda_pos = lambda_pos;
        best.beam_mom.SetXYZ(pxs[i] + beam_t*(pxs[i+1]-pxs[i]),
                             pys[i] + beam_t*(pys[i+1]-pys[i]),
                             pzs[i] + beam_t*(pzs[i+1]-pzs[i]));
        best.dca = TMath::Sqrt(d2);
        best.beam_segment_id = static_cast<Int_t>(i);
        best.beam_segment_fraction = beam_t;
      };
      // Interior stationary point, when the beam segment and Lambda line are not parallel.
      const Double_t denom = a*c-b*b;
      if (denom > 1.e-12) consider((b*e-c*d)/denom, (a*e-b*d)/denom);
      // Boundary candidates cover the constrained segment endpoints and the
      // Lambda decay endpoint (s=0).
      consider(0., std::min(0., e));
      consider(1., std::min(0., e+b));
      consider(std::clamp(-d/a, 0., 1.), 0.);
    }
    return best;
  }

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
  // Defined only when the beam--Lambda production closest approach is found.
  std::vector<Double_t> lambda_beam_opening_angle_lab;
  std::vector<Double_t> lambda_production_angle_cm;
  std::vector<Double_t> lambda_production_costheta_cm;
  std::vector<Int_t> lambda_proton_track_id;
  std::vector<Int_t> lambda_pion_track_id;
  std::vector<Double_t> lambda_proton_mom_x;
  std::vector<Double_t> lambda_proton_mom_y;
  std::vector<Double_t> lambda_proton_mom_z;
  std::vector<Double_t> lambda_pion_mom_x;
  std::vector<Double_t> lambda_pion_mom_y;
  std::vector<Double_t> lambda_pion_mom_z;
  
  // X_mass uses the beam momentum at the reconstructed Lambda production vertex.
  // The *_decay_z branches retain the previous decay-z prescription for comparison.
  std::vector<Double_t> X_mass;
  std::vector<Double_t> X_mass2;
  std::vector<Double_t> X_mass_decay_z;
  std::vector<Double_t> X_mass2_decay_z;
  std::vector<Int_t> X_prod_found;
  std::vector<Double_t> X_prod_vtx_x;
  std::vector<Double_t> X_prod_vtx_y;
  std::vector<Double_t> X_prod_vtx_z;
  // Beam momentum at the reconstructed Lambda production vertex, after the
  // same material-loss correction used for X_mass and the CM angles.
  std::vector<Double_t> X_prod_beam_px;
  std::vector<Double_t> X_prod_beam_py;
  std::vector<Double_t> X_prod_beam_pz;
  std::vector<Double_t> X_prod_beam_mom;
  std::vector<Double_t> X_beam_lambda_dca;
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

  void clearSelectedLambdaEta() {
    dst::clear_all(
      lambda_mass, lambda_close_dist,
      lambda_vtx_x, lambda_vtx_y, lambda_vtx_z,
      lambda_mom_x, lambda_mom_y, lambda_mom_z,
      lambda_beam_opening_angle_lab, lambda_production_angle_cm, lambda_production_costheta_cm,
      lambda_proton_track_id, lambda_pion_track_id,
      lambda_proton_mom_x, lambda_proton_mom_y, lambda_proton_mom_z,
      lambda_pion_mom_x, lambda_pion_mom_y, lambda_pion_mom_z,
      X_mass, X_mass2, X_mass_decay_z, X_mass2_decay_z,
      X_prod_found, X_prod_vtx_x, X_prod_vtx_y, X_prod_vtx_z,
      X_prod_beam_px, X_prod_beam_py, X_prod_beam_pz, X_prod_beam_mom,
      X_beam_lambda_dca,
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
    clearSelectedLambdaEta();
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_xvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* blc_yvpHS;
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
  Bool_t has_beam_track = false;
  if (blc_entry >= 0) {
    TTreeReaderCont[kBlc]->SetEntry(blc_entry);
    if (**src.blc_runnum == event.runnum && **src.blc_evnum == event.evnum) {
      has_beam_track = CopyBestK18BeamMomentum(**src.blc_px, **src.blc_py, **src.blc_pz,
                                                event.beam_px, event.beam_py, event.beam_pz);
    }
  }
  SanitizeBeamMomentum(event.beam_px, event.beam_py, event.beam_pz);
  has_beam_track = has_beam_track &&
    TMath::Sqrt(event.beam_px*event.beam_px + event.beam_py*event.beam_py +
                event.beam_pz*event.beam_pz) > 0.;

  event.ntTpc = **src.ntTpc;
  event.effective_ntTpc = **src.effective_ntTpc;
  if (kRequireEffectiveNTracks &&
      event.effective_ntTpc != kRequiredEffectiveNTracks) return false;
  event.nhtrack = **src.nhtrack;
  event.is_beam = **src.is_beam;
  // Keep the BLC/RK-beam requirement.  Do not reject events only because
  // the corresponding K18 beam track was not tagged in the TPC.
  if (kRequireBeamTrack && !has_beam_track) return false;
  HF1("Status", event.status++);
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
  for (Int_t it = 0; it < event.ntTpc; ++it) {
    const Double_t track_dz = event.helix_dz.at(it);
    if (std::isfinite(track_dz)) {
      HF1("LambdaEta_Track_HelixDz_All", track_dz);
      if (event.is_beam.at(it)) HF1("LambdaEta_Track_HelixDz_Beam", track_dz);
      if (event.is_k18.at(it)) HF1("LambdaEta_Track_HelixDz_K18", track_dz);
      if (event.is_accidental.at(it)) HF1("LambdaEta_Track_HelixDz_Accidental", track_dz);
      if (!event.is_beam.at(it) && !event.is_k18.at(it) && !event.is_accidental.at(it)) HF1("LambdaEta_Track_HelixDz_Other", track_dz);
    }
    Double_t track_dx = 0.;
    if (CalcHelixDx(event.helix_r.at(it), event.helix_theta_min.at(it),
                    event.helix_theta_max.at(it), track_dx)) {
      auto fill_dx = [&](const char* category) {
        HF1(Form("LambdaEta_Track_HelixDx_%s", category), track_dx);
      };
      fill_dx("All");
      if (event.is_beam.at(it)) fill_dx("Beam");
      if (event.is_k18.at(it)) fill_dx("K18");
      if (event.is_accidental.at(it)) fill_dx("Accidental");
      if (!event.is_beam.at(it) && !event.is_k18.at(it) && !event.is_accidental.at(it)) fill_dx("Other");
    }
    const auto nhit = std::min(event.hitlayer.at(it).size(), event.track_cluster_row_center.at(it).size());
    for (std::size_t ih = 0; ih < nhit; ++ih) {
      if (!std::isfinite(event.hitlayer.at(it).at(ih)) || !std::isfinite(event.track_cluster_row_center.at(it).at(ih))) continue;
      const Int_t layer = static_cast<Int_t>(std::lround(event.hitlayer.at(it).at(ih)));
      const Int_t row = static_cast<Int_t>(std::lround(event.track_cluster_row_center.at(it).at(ih)));
      if (layer < 0 || layer >= NumOfLayersTPC || row < 0) continue;
      const Int_t pad_id = tpc::GetPadId(layer, row);
      if (pad_id < 0) continue;
      auto fill_pad = [&](const char* name) { HF2Poly(name, pad_id + 1, HG2Poly(name, pad_id + 1) + 1.); };
      if (event.is_k18.at(it)) fill_pad("LambdaEta_TrackPad_K18");
      if (event.is_beam.at(it)) fill_pad("LambdaEta_TrackPad_Beam");
      if (event.is_accidental.at(it)) fill_pad("LambdaEta_TrackPad_Accidental");
      if (!event.is_beam.at(it) && !event.is_k18.at(it) && !event.is_accidental.at(it)) fill_pad("LambdaEta_TrackPad_Other");
    }
  }
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

      // Use the beam--Lambda closest approach as the Lambda production vertex.
      // Retain the old decay-z result above solely as a diagnostic comparison.
      const auto production = FindBeamLambdaProductionVertex(
        **src.blc_xvpHS, **src.blc_yvpHS, **src.blc_zvpHS,
        **src.blc_pxvpHS, **src.blc_pyvpHS, **src.blc_pzvpHS,
        TVector3(vtx_x, vtx_y, vtx_z), lambda_mom);
      const Bool_t valid_production = production.found;
      const Double_t qnan = TMath::QuietNaN();
      TLorentzVector missing_production;
      missing_production.SetPxPyPzE(qnan, qnan, qnan, qnan);
      Double_t lambda_beam_opening_angle_lab = TMath::QuietNaN();
      Double_t lambda_production_angle_cm = TMath::QuietNaN();
      Double_t lambda_production_costheta_cm = TMath::QuietNaN();
      // Keep the post-loss vector: this is precisely the beam momentum used
      // below to calculate X_mass and the production angles.
      Double_t prod_beam_px = TMath::QuietNaN();
      Double_t prod_beam_py = TMath::QuietNaN();
      Double_t prod_beam_pz = TMath::QuietNaN();
      if (valid_production) {
        prod_beam_px = production.beam_mom.X();
        prod_beam_py = production.beam_mom.Y();
        prod_beam_pz = production.beam_mom.Z();
        ApplyBeamMaterialLoss(prod_beam_px, prod_beam_py, prod_beam_pz);
        TLorentzVector production_beam;
        production_beam.SetXYZM(prod_beam_px, prod_beam_py, prod_beam_pz, pdg::KaonMass());
        missing_production = MissingFromLambda(lambda_mom, lambda.M(),
                                                prod_beam_px, prod_beam_py, prod_beam_pz);
        if (production_beam.Vect().Mag2() > 0. && lambda_mom.Mag2() > 0.) {
          lambda_beam_opening_angle_lab = production_beam.Vect().Angle(lambda_mom)*TMath::RadToDeg();
          TLorentzVector target; target.SetXYZM(0., 0., 0., pdg::ProtonMass());
          TLorentzVector beam_cm = production_beam;
          TLorentzVector lambda_cm = lambda;
          const TVector3 cm_boost = -(production_beam + target).BoostVector();
          beam_cm.Boost(cm_boost);
          lambda_cm.Boost(cm_boost);
          if (beam_cm.Vect().Mag2() > 0. && lambda_cm.Vect().Mag2() > 0.) {
            lambda_production_angle_cm = beam_cm.Vect().Angle(lambda_cm.Vect())*TMath::RadToDeg();
            lambda_production_costheta_cm = TMath::Cos(lambda_production_angle_cm*TMath::DegToRad());
          }
        }
        HF1("LambdaEta_MissingMass_ProductionVtx", missing_production.M());
        HF1("LambdaEta_BeamLambdaProductionDCA", production.dca);
        HF2("LambdaEta_ProductionVtxZX", production.lambda_pos.Z(), production.lambda_pos.X());
        const Bool_t production_below_220 = production.lambda_pos.Z() < -220.;
        const char* decay_zx_name = production_below_220 ? "LambdaEta_DecayVtxZX_ProductionBelow220" : "LambdaEta_DecayVtxZX_ProductionAbove220";
        const char* decay_zy_name = production_below_220 ? "LambdaEta_DecayVtxZY_ProductionBelow220" : "LambdaEta_DecayVtxZY_ProductionAbove220";
        const char* decay_xy_name = production_below_220 ? "LambdaEta_DecayVtxXY_ProductionBelow220" : "LambdaEta_DecayVtxXY_ProductionAbove220";
        HF2(decay_zx_name, vtx_z, vtx_x);
        HF2(decay_zy_name, vtx_z, vtx_y);
        HF2(decay_xy_name, vtx_x, vtx_y);
        const char* proton_pid_name = production_below_220 ? "LambdaEta_Proton_dEdx_vs_ChargeMom_ProductionBelow220" : "LambdaEta_Proton_dEdx_vs_ChargeMom_ProductionAbove220";
        const char* pion_pid_name = production_below_220 ? "LambdaEta_Pion_dEdx_vs_ChargeMom_ProductionBelow220" : "LambdaEta_Pion_dEdx_vs_ChargeMom_ProductionAbove220";
        HF2(proton_pid_name, event.charge.at(ip) * event.mom0.at(ip), event.dEdx.at(ip));
        HF2(pion_pid_name, event.charge.at(ipi) * event.mom0.at(ipi), event.dEdx.at(ipi));
        const char* proton_pad_name = production_below_220 ? "LambdaEta_ProtonTrackPad_ProductionBelow220" : "LambdaEta_ProtonTrackPad_ProductionAbove220";
        const char* pion_pad_name = production_below_220 ? "LambdaEta_PionTrackPad_ProductionBelow220" : "LambdaEta_PionTrackPad_ProductionAbove220";
        for (const auto& track_info : std::initializer_list<std::pair<Int_t, const char*>>{{ip, proton_pad_name}, {ipi, pion_pad_name}}) {
          const Int_t track_id = track_info.first;
          const auto nhit = std::min(event.hitlayer.at(track_id).size(), event.track_cluster_row_center.at(track_id).size());
          for (std::size_t ih = 0; ih < nhit; ++ih) {
            if (!std::isfinite(event.hitlayer.at(track_id).at(ih)) || !std::isfinite(event.track_cluster_row_center.at(track_id).at(ih))) continue;
            const Int_t layer = static_cast<Int_t>(std::lround(event.hitlayer.at(track_id).at(ih)));
            const Int_t row = static_cast<Int_t>(std::lround(event.track_cluster_row_center.at(track_id).at(ih)));
            if (layer < 0 || layer >= NumOfLayersTPC || row < 0) continue;
            const Int_t pad_id = tpc::GetPadId(layer, row);
            if (pad_id >= 0) HF2Poly(track_info.second, pad_id + 1, HG2Poly(track_info.second, pad_id + 1) + 1.);
          }
        }
        if (production_below_220 && std::isfinite(event.helix_dz.at(ipi))) {
          const Double_t pion_dz = event.helix_dz.at(ipi);
          HF1("LambdaEta_Pion_HelixDz_ProductionBelow220", pion_dz);
          HF1("LambdaEta_Pion_Mom0X_ProductionBelow220", event.mom0_x.at(ipi));
          HF1("LambdaEta_Pion_Mom0Y_ProductionBelow220", event.mom0_y.at(ipi));
          HF1("LambdaEta_Pion_Mom0Z_ProductionBelow220", event.mom0_z.at(ipi));
          HF1("LambdaEta_Pion_Mom0_ProductionBelow220", event.mom0.at(ipi));
          Double_t pion_dx = 0.;
          if (CalcHelixDx(event.helix_r.at(ipi),
                          event.helix_theta_min.at(ipi),
                          event.helix_theta_max.at(ipi), pion_dx))
            HF1("LambdaEta_Pion_HelixDx_ProductionBelow220", pion_dx);
          const Bool_t beamlike = event.charge.at(ipi) < 0 && TMath::Abs(pion_dz) < 0.05;
          HF1(beamlike ? "LambdaEta_Pion_HelixDz_BeamLikePass_ProductionBelow220" : "LambdaEta_Pion_HelixDz_BeamLikeFail_ProductionBelow220", pion_dz);
        }
        if (std::isfinite(lambda_beam_opening_angle_lab)) HF1("LambdaEta_LambdaBeamOpeningAngleLab", lambda_beam_opening_angle_lab);
        if (std::isfinite(lambda_production_angle_cm)) HF1("LambdaEta_LambdaProductionAngleCM", lambda_production_angle_cm);
        if (std::isfinite(lambda_production_costheta_cm)) HF1("LambdaEta_LambdaProductionCosThetaCM", lambda_production_costheta_cm);
      }
      HF1("LambdaEta_MissingMass_DecayZ", missing.M());

      HF1("LambdaEta_Mass", lambda.M());
      if (IsStrictlyInsideLH2(vtx_x, vtx_y, vtx_z)) {
        HF1("LambdaEta_LambdaMass_LH2Inside", lambda.M());
        HF2("LambdaEta_VtxZX_LH2Inside", vtx_z, vtx_x);
      }
      if (IsCompletelyOutsideTargetFrame(vtx_x, vtx_y, vtx_z)) {
        HF1("LambdaEta_LambdaMass_TargetFrameOutside", lambda.M());
        HF2("LambdaEta_VtxZX_TargetFrameOutside", vtx_z, vtx_x);
      }
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
      event.lambda_beam_opening_angle_lab.push_back(lambda_beam_opening_angle_lab);
      event.lambda_production_angle_cm.push_back(lambda_production_angle_cm);
      event.lambda_production_costheta_cm.push_back(lambda_production_costheta_cm);
      event.lambda_proton_track_id.push_back(ip);
      event.lambda_pion_track_id.push_back(ipi);
      event.lambda_proton_mom_x.push_back(p_mom.X());
      event.lambda_proton_mom_y.push_back(p_mom.Y());
      event.lambda_proton_mom_z.push_back(p_mom.Z());
      event.lambda_pion_mom_x.push_back(pi_mom.X());
      event.lambda_pion_mom_y.push_back(pi_mom.Y());
      event.lambda_pion_mom_z.push_back(pi_mom.Z());
      event.X_mass.push_back(missing_production.M());
      event.X_mass2.push_back(missing_production.M2());
      event.X_mass_decay_z.push_back(missing.M());
      event.X_mass2_decay_z.push_back(missing.M2());
      event.X_prod_found.push_back(valid_production);
      event.X_prod_vtx_x.push_back(valid_production ? production.lambda_pos.X() : qnan);
      event.X_prod_vtx_y.push_back(valid_production ? production.lambda_pos.Y() : qnan);
      event.X_prod_vtx_z.push_back(valid_production ? production.lambda_pos.Z() : qnan);
      event.X_prod_beam_px.push_back(prod_beam_px);
      event.X_prod_beam_py.push_back(prod_beam_py);
      event.X_prod_beam_pz.push_back(prod_beam_pz);
      event.X_prod_beam_mom.push_back(std::sqrt(prod_beam_px*prod_beam_px +
                                                 prod_beam_py*prod_beam_py +
                                                 prod_beam_pz*prod_beam_pz));
      event.X_beam_lambda_dca.push_back(valid_production ? production.dca : qnan);
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
  dst::SetBranch(blc_reader, "xvpHS", src.blc_xvpHS);
  dst::SetBranch(blc_reader, "yvpHS", src.blc_yvpHS);
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
  hist::BuildTPCLambdaEta();

  tree = new TTree("tpc", "tree of DstTPCLambdaEta");
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
  tree->Branch("lambda_beam_opening_angle_lab", &event.lambda_beam_opening_angle_lab);
  tree->Branch("lambda_production_angle_cm", &event.lambda_production_angle_cm);
  tree->Branch("lambda_production_costheta_cm", &event.lambda_production_costheta_cm);
  tree->Branch("lambda_proton_track_id", &event.lambda_proton_track_id);
  tree->Branch("lambda_pion_track_id", &event.lambda_pion_track_id);
  tree->Branch("lambda_proton_mom_x", &event.lambda_proton_mom_x);
  tree->Branch("lambda_proton_mom_y", &event.lambda_proton_mom_y);
  tree->Branch("lambda_proton_mom_z", &event.lambda_proton_mom_z);
  tree->Branch("lambda_pion_mom_x", &event.lambda_pion_mom_x);
  tree->Branch("lambda_pion_mom_y", &event.lambda_pion_mom_y);
  tree->Branch("lambda_pion_mom_z", &event.lambda_pion_mom_z);


  tree->Branch("X_mass", &event.X_mass);
  tree->Branch("X_mass2", &event.X_mass2);
  tree->Branch("X_mass_decay_z", &event.X_mass_decay_z);
  tree->Branch("X_mass2_decay_z", &event.X_mass2_decay_z);
  tree->Branch("X_prod_found", &event.X_prod_found);
  tree->Branch("X_prod_vtx_x", &event.X_prod_vtx_x);
  tree->Branch("X_prod_vtx_y", &event.X_prod_vtx_y);
  tree->Branch("X_prod_vtx_z", &event.X_prod_vtx_z);
  tree->Branch("X_prod_beam_px", &event.X_prod_beam_px);
  tree->Branch("X_prod_beam_py", &event.X_prod_beam_py);
  tree->Branch("X_prod_beam_pz", &event.X_prod_beam_pz);
  tree->Branch("X_prod_beam_mom", &event.X_prod_beam_mom);
  tree->Branch("X_beam_lambda_dca", &event.X_beam_lambda_dca);
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
