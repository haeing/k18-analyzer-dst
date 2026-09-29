// -*- C++ -*-

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <TPDGCode.h>
#include <TVector3.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
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

#define RawHit 1
#define RawCluster 1
#define TrackSearch 1
#define TrackCluster 1
#define TruncatedMean 0
#define CalibHist 0 // enable/disable per-pad calibration histograms

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

  const std::vector<TString> kUserParamKeys = {
    // Cluster building (ReCalcTPCHits / MakeUpTPCClusters)
    "MinCDeTPC", "MaxYDifClusterTPC",
    "MinClusterDeTPC", "MinClusterSizeTPC",
    "MinClusterYPosTPC", "MaxClusterYPosTPC",

    // Helix tracking
    "MinLayerTPC", "MaxHoughWindow", "MaxHoughWindowY", "BeamThroughTPC",

    // Helix-track error scaling (TPCLocalTrackHelix ctor)
    "MomResScale", "dZResScale", "PhiResScale",
    
    // Vertex finding (TPCVertex)
    "VertexScanRange",
    
    // Optional parameter (default value is provided in the code)
    // "MaxCenterRowDiffTPC",
    // "BeamLikeMaxAbsDzTPC",  // |dz| threshold for post-fit is_beam tagging (TPCTrackSearch; default 0.05)
  };
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpcHit, kK18, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCHit]", "[K18Tracking]", "[OutFile]" };
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
  std::vector<Double_t> clkTpc;

#if RawHit
  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;
#endif

#if RawCluster
  Int_t nclTpc;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_x_center;
  std::vector<Double_t> cluster_y_center;
  std::vector<Double_t> cluster_z_center;
  std::vector<Int_t> cluster_row_center;
  std::vector<Int_t> cluster_houghflag;
#endif

  // K18 RK VP points and the helix fitted from those points (for event display).
  Int_t ntK18;
  std::vector<Int_t> nh_k18;
  std::vector<std::vector<Double_t>> x_hit_k18;
  std::vector<std::vector<Double_t>> y_hit_k18;
  std::vector<std::vector<Double_t>> z_hit_k18;
  std::vector<std::vector<Double_t>> x_cal_k18;
  std::vector<std::vector<Double_t>> y_cal_k18;
  std::vector<std::vector<Double_t>> z_cal_k18;
  std::vector<Double_t> helix_cx_k18;
  std::vector<Double_t> helix_cy_k18;
  std::vector<Double_t> helix_z0_k18;
  std::vector<Double_t> helix_r_k18;
  std::vector<Double_t> helix_dz_k18;
  std::vector<Double_t> helix_theta_min_k18;
  std::vector<Double_t> helix_theta_max_k18;
  // Momentum of the VP-only helix fit at y = 0, matching TPC mom0_* convention.
  std::vector<Double_t> mom0_x_k18;
  std::vector<Double_t> mom0_y_k18;
  std::vector<Double_t> mom0_z_k18;
  std::vector<Double_t> mom0_k18;
  // Original K18 transfer-matrix state at VO.
  std::vector<Double_t> xout_k18, yout_k18, uout_k18, vout_k18, p_k18;
  Int_t ntTpc; // Number of Tracks
  Int_t effective_ntTpc; // Number of tracks with no beam/accidental/k18 flag
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> is_beam; // 1 = beam-tagged, 0 = not
  std::vector<Int_t> is_k18; // 1 = K18 VP-constrained beam track
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

  std::vector<Double_t> dEdx_0;
  std::vector<Double_t> dEdx_10;
  std::vector<Double_t> dEdx_20;
  std::vector<Double_t> dEdx_30;
  std::vector<Double_t> dEdx_40;
  std::vector<Double_t> dEdx_50;
  std::vector<Double_t> dEdx_60;
  std::vector<Double_t> dEdx_cor_0;
  std::vector<Double_t> dEdx_cor_10;
  std::vector<Double_t> dEdx_cor_20;
  std::vector<Double_t> dEdx_cor_30;
  std::vector<Double_t> dEdx_cor_40;
  std::vector<Double_t> dEdx_cor_50;
  std::vector<Double_t> dEdx_cor_60;

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
  // Outer index: reconstructed TPC helix track. Inner index: every good pad
  // hit in the clusters used by that track (not merely one entry per cluster).
  // transverse_diffusion_track_hit_index maps an inner entry to hitpos_*/
  // track_cluster_* so the existing track-hit conditions can be applied.
  std::vector<std::vector<Double_t>> transverse_diffusion_dist;
  std::vector<std::vector<Double_t>> transverse_diffusion_ratio;
  std::vector<std::vector<Double_t>> transverse_diffusion_drift;
  std::vector<std::vector<Double_t>> transverse_diffusion_track_hit_index;


  
  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc);
  }
  void clearK18HelixTracks() {
    ntK18 = 0;
    dst::clear_all(
      nh_k18, x_hit_k18, y_hit_k18, z_hit_k18, x_cal_k18, y_cal_k18, z_cal_k18,
      helix_cx_k18, helix_cy_k18, helix_z0_k18, helix_r_k18, helix_dz_k18,
      helix_theta_min_k18, helix_theta_max_k18,
      mom0_x_k18, mom0_y_k18, mom0_z_k18, mom0_k18
      , xout_k18, yout_k18, uout_k18, vout_k18, p_k18
    );
  }


#if RawHit
  void clearRawHits() {
    nhTpc = 0;
    dst::clear_all(
      raw_hitpos_x, raw_hitpos_y, raw_hitpos_z,
      raw_de, raw_padid, raw_layer, raw_row
    );
  }
#endif

#if RawCluster
  void clearClusters() {
    nclTpc = 0;
    dst::clear_all(
      cluster_x, cluster_y, cluster_z, cluster_de,
      cluster_size, cluster_layer, cluster_mrow,
      cluster_de_center, cluster_x_center, cluster_y_center,
      cluster_z_center, cluster_row_center, cluster_houghflag
    );
  }
#endif

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

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center,
      transverse_diffusion_dist, transverse_diffusion_ratio, transverse_diffusion_drift, transverse_diffusion_track_hit_index
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
    clearK18HelixTracks();
    clearHelixTracks();
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

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center,
      transverse_diffusion_dist, transverse_diffusion_ratio, transverse_diffusion_drift, transverse_diffusion_track_hit_index
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
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;    // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;     // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time

  TTreeReaderValue<UInt_t>* k18_runnum;
  TTreeReaderValue<UInt_t>* k18_evnum;
  TTreeReaderValue<std::vector<Int_t>>* rkStatusHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* zvpHS;
  TTreeReaderValue<std::vector<Double_t>>* xoutK18;
  TTreeReaderValue<std::vector<Double_t>>* youtK18;
  TTreeReaderValue<std::vector<Double_t>>* uoutK18;
  TTreeReaderValue<std::vector<Double_t>>* voutK18;
  TTreeReaderValue<std::vector<Double_t>>* pK18;
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

  //_____________________________________________________________________________
  std::vector<std::vector<TVector3>>
  MakeK18VPs()
  {
    std::vector<std::vector<TVector3>> vps;
    const Long64_t k18_entry = TTreeCont[kK18]->GetEntryNumberWithIndex(event.runnum, event.evnum);
    if (k18_entry < 0) return vps;

    TTreeReaderCont[kK18]->SetEntry(k18_entry);
    if (**src.k18_runnum != event.runnum || **src.k18_evnum != event.evnum)
      return vps;

    const auto& status = **src.rkStatusHS;
    const auto& xs = **src.xvpHS;
    const auto& ys = **src.yvpHS;
    const auto& zs = **src.zvpHS;
    const std::size_t ntrack = std::min({status.size(), xs.size(), ys.size(), zs.size()});
    vps.reserve(ntrack);

    for (std::size_t it = 0; it < ntrack; ++it) {
      if (status[it] != 1) continue;
      const std::size_t nhit = std::min({xs[it].size(), ys[it].size(), zs[it].size()});
      std::vector<TVector3> vp_track;
      vp_track.reserve(nhit);
      for (std::size_t ih = 0; ih < nhit; ++ih) {
        const Double_t x = xs[it][ih];
        const Double_t y = ys[it][ih];
        const Double_t z = zs[it][ih];
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
          continue;
        vp_track.emplace_back(x, y, z);
      }
      if (vp_track.size() >= 3)
        vps.push_back(vp_track);
    }
    return vps;
  }
  
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
  Bool_t PositionOnHelixAtZ(const TPCLocalTrackHelix* track,
                            Double_t z, TVector3& position)
  {
    if (!track) return false;
    Double_t par[5];
    track->GetParam(par);
    const Double_t radius = par[3];
    if (!std::isfinite(radius) || TMath::Abs(radius) < 1.e-9) return false;

    Double_t sin_theta = (z - tpc::Z_TARGET - par[1])/radius;
    if (sin_theta < -1. || sin_theta > 1.) return false;
    sin_theta = std::max(-1., std::min(1., sin_theta));

    const Double_t min_theta = track->GetMint();
    const Double_t max_theta = track->GetMaxt();
    if (!std::isfinite(min_theta) || !std::isfinite(max_theta) || min_theta > max_theta)
      return false;

    const Double_t theta_roots[2] = {TMath::ASin(sin_theta), TMath::Pi() - TMath::ASin(sin_theta)};
    const Double_t theta_center = 0.5*(min_theta + max_theta);
    Double_t best_theta = std::numeric_limits<Double_t>::quiet_NaN();
    Double_t best_distance = std::numeric_limits<Double_t>::infinity();
    const Double_t two_pi = 2.*TMath::Pi();
    for (const Double_t theta_root : theta_roots) {
      const Int_t n_min = static_cast<Int_t>(std::floor((min_theta - theta_root)/two_pi)) - 1;
      const Int_t n_max = static_cast<Int_t>(std::ceil((max_theta - theta_root)/two_pi)) + 1;
      for (Int_t n = n_min; n <= n_max; ++n) {
        const Double_t theta = theta_root + two_pi*n;
        if (theta < min_theta || theta > max_theta) continue;
        const Double_t distance = TMath::Abs(theta - theta_center);
        if (distance < best_distance) {
          best_distance = distance;
          best_theta = theta;
        }
      }
    }

    if (!std::isfinite(best_theta)) return false;
    position = track->GetPosition(par, best_theta);
    return true;
  }

  //_____________________________________________________________________________
  Double_t MeanUpstreamPlaneYResidual(const TPCLocalTrackHelix* ref_track,
                                     const TPCLocalTrackHelix* cand_track,
                                     Int_t& n_planes)
  {
    n_planes = 0;
    if (!ref_track || !cand_track) return std::numeric_limits<Double_t>::infinity();

    Double_t sum_residual = 0.;
    for (Int_t ivp = 0; ivp < ref_track->GetVPNHit(); ++ivp) {
      const Double_t z = ref_track->GetVPPos(ivp).Z();
      if (!std::isfinite(z) || z >= tpc::Z_TARGET) continue;

      TVector3 ref_pos, cand_pos;
      if (!PositionOnHelixAtZ(ref_track, z, ref_pos) ||
          !PositionOnHelixAtZ(cand_track, z, cand_pos)) continue;
      const Double_t residual = TMath::Abs(cand_pos.Y() - ref_pos.Y());
      if (!std::isfinite(residual)) continue;
      sum_residual += residual;
      ++n_planes;
    }

    if (n_planes == 0) return std::numeric_limits<Double_t>::infinity();
    return sum_residual/static_cast<Double_t>(n_planes);
  }
  //_____________________________________________________________________________
  void FillK18HelixTracks(TPCAnalyzer& tpc_ana)
  {
    event.ntK18 = tpc_ana.GetNTracksTPCHelixVP();
    dst::resize_all(event.ntK18,
      event.nh_k18, event.x_hit_k18, event.y_hit_k18, event.z_hit_k18,
      event.x_cal_k18, event.y_cal_k18, event.z_cal_k18,
      event.helix_cx_k18, event.helix_cy_k18, event.helix_z0_k18,
      event.helix_r_k18, event.helix_dz_k18,
      event.helix_theta_min_k18, event.helix_theta_max_k18,
      event.mom0_x_k18, event.mom0_y_k18, event.mom0_z_k18, event.mom0_k18
      , event.xout_k18, event.yout_k18, event.uout_k18, event.vout_k18, event.p_k18
    );

    const Double_t qnan = TMath::QuietNaN();
    for (Int_t ihs = 0; ihs < event.ntK18; ++ihs) {
      TPCLocalTrackHelix* track = tpc_ana.GetTrackTPCHelixVP(ihs);
      if (!track) continue;

      event.nh_k18[ihs] = track->GetVPNHit();
      event.helix_cx_k18[ihs] = track->Getcx();
      event.helix_cy_k18[ihs] = track->Getcy();
      event.helix_z0_k18[ihs] = track->Getz0();
      event.helix_r_k18[ihs] = track->Getr();
      event.helix_dz_k18[ihs] = track->Getdz();
      event.helix_theta_min_k18[ihs] = track->GetMint();
      event.helix_theta_max_k18[ihs] = track->GetMaxt();
      const TVector3 mom0_k18 = track->GetMom0();
      event.mom0_x_k18[ihs] = mom0_k18.X();
      event.mom0_y_k18[ihs] = mom0_k18.Y();
      event.mom0_z_k18[ihs] = mom0_k18.Z();
      event.mom0_k18[ihs] = mom0_k18.Mag();
      const auto& status = **src.rkStatusHS;
      const auto& xout = **src.xoutK18; const auto& yout = **src.youtK18;
      const auto& uout = **src.uoutK18; const auto& vout = **src.voutK18; const auto& p = **src.pK18;
      std::size_t ik18 = 0, matched = 0;
      for (; ik18 < status.size(); ++ik18) {
        if (status[ik18] != 1) continue;
        if (matched++ == static_cast<std::size_t>(ihs)) break;
      }
      if (ik18 < xout.size() && ik18 < yout.size() && ik18 < uout.size() && ik18 < vout.size() && ik18 < p.size()) {
        event.xout_k18[ihs] = xout[ik18]; event.yout_k18[ihs] = yout[ik18];
        event.uout_k18[ihs] = uout[ik18]; event.vout_k18[ihs] = vout[ik18]; event.p_k18[ihs] = p[ik18];
      }

      dst::resize_all(event.nh_k18[ihs],
        event.x_hit_k18[ihs], event.y_hit_k18[ihs], event.z_hit_k18[ihs],
        event.x_cal_k18[ihs], event.y_cal_k18[ihs], event.z_cal_k18[ihs]
      );
      for (Int_t ivp = 0; ivp < event.nh_k18[ihs]; ++ivp) {
        const TVector3 hit = track->GetVPPos(ivp);
        event.x_hit_k18[ihs][ivp] = hit.X();
        event.y_hit_k18[ihs][ivp] = hit.Y();
        event.z_hit_k18[ihs][ivp] = hit.Z();

        TVector3 cal(qnan, qnan, qnan);
        PositionOnHelixAtZ(track, hit.Z(), cal);
        event.x_cal_k18[ihs][ivp] = cal.X();
        event.y_cal_k18[ihs][ivp] = cal.Y();
        event.z_cal_k18[ihs][ivp] = cal.Z();
      }
    }
  }


  //_____________________________________________________________________________
  void FillResidualsAtVPZ(const TPCLocalTrackHelix* ref_track,
                          const TPCLocalTrackHelix* cand_track)
  {
    for (Int_t ivp = 0; ivp < ref_track->GetVPNHit(); ++ivp) {
      const Double_t z = ref_track->GetVPPos(ivp).Z();
      if (!std::isfinite(z) || z >= tpc::Z_TARGET) continue;

      TVector3 ref_pos, cand_pos;
      if (!PositionOnHelixAtZ(ref_track, z, ref_pos) ||
          !PositionOnHelixAtZ(cand_track, z, cand_pos)) continue;
      const Double_t delta_x = cand_pos.X() - ref_pos.X();
      const Double_t delta_y = cand_pos.Y() - ref_pos.Y();
      const Double_t residual = TMath::Hypot(delta_x, delta_y);
      if (!std::isfinite(delta_x) || !std::isfinite(delta_y) || !std::isfinite(residual))
        continue;

      HF1("TPCHS_DeltaX", delta_x);
      HF1("TPCHS_DeltaY", delta_y);
      HF1("TPCHS_Residual", residual);
      HF2("TPCHS_DeltaX_vs_Z", z, delta_x);
    }
  }

  //_____________________________________________________________________________
  void FillHSHelixResidualMatch(TPCAnalyzer& tpc_ana)
  {
    const Int_t nt_tpc = tpc_ana.GetNTracksTPCHelix();
    const Int_t nt_vp = tpc_ana.GetNTracksTPCHelixVP();

    for (Int_t ik18 = 0; ik18 < nt_vp; ++ik18) {
      TPCLocalTrackHelix* ref_track = tpc_ana.GetTrackTPCHelixVP(ik18);
      if (!ref_track) continue;

      Double_t best_mean_y_residual = std::numeric_limits<Double_t>::infinity();
      TPCLocalTrackHelix* best_track = nullptr;
      for (Int_t it = 0; it < nt_tpc; ++it) {
        TPCLocalTrackHelix* cand_track = tpc_ana.GetTrackTPCHelix(it);
        if (!cand_track || cand_track->GetIsK18() != 0) continue;

        Int_t n_planes = 0;
        const Double_t mean_residual =
          MeanUpstreamPlaneYResidual(ref_track, cand_track, n_planes);
        if (n_planes == 0 || !std::isfinite(mean_residual)) continue;
        if (mean_residual < best_mean_y_residual) {
          best_mean_y_residual = mean_residual;
          best_track = cand_track;
        }
      }

      if (best_track)
        FillResidualsAtVPZ(ref_track, best_track);
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
    const Int_t is_beam = helix_track->GetIsBeam();
    const Int_t is_accidental = helix_track->GetIsAccidental();
    const Int_t is_k18        = helix_track->GetIsK18();
    event.is_k18[it] = is_k18;
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

      // Use the charge-weighted cluster position as the pad-by-pad reference.
      // These flattened inner vectors contain all good pads from this track.
      if (cl_de > 0.) {
        const TPCHitContainer& pad_hits = cl->GetHitContainer();
        for (const auto& pad_hit : pad_hits) {
          if (!pad_hit || !pad_hit->IsGood()) continue;
          const TVector3& pad_pos = pad_hit->GetPosition();
          Double_t trans_dist = std::hypot(hit_pos.x()-pad_pos.x(),
                                           hit_pos.z()-pad_pos.z());
          if (hit_pos.x() < pad_pos.x()) trans_dist = -trans_dist;
          event.transverse_diffusion_dist[it].push_back(trans_dist);
          event.transverse_diffusion_ratio[it].push_back(pad_hit->GetCDe()/cl_de);
          event.transverse_diffusion_drift[it].push_back(std::abs(hit_pos.y()-304.2)); // GEM distance [mm]
          event.transverse_diffusion_track_hit_index[it].push_back(static_cast<Double_t>(ih));
        }
      }

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
  if (TTreeCont[kK18]->BuildIndex("run_number", "event_number") < 0) {
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
  HF1("Status", event.status++);

  
  if (**src.nhTpc == 0) return true;
  HF1("Status", event.status++);
  
  if (!TPCEventAnalyzer::ValidateCoboClocks(event.clkTpc)) {
    return true;
  }
  HF1("Status", event.status++);
  
  TPCAnalyzer tpc_ana;
  tpc_ana.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  HF1("Status", event.status++);

  static TPCEventAnalyzer event_ana;
  event_ana.SetClock(event.clkTpc);

#if RawHit
  FillRawHits(tpc_ana);
  HF1("Status", event.status++);
#endif

#if RawCluster
  FillClusters(tpc_ana);
  HF1("Status", event.status++);
#endif

#if TrackSearch
  const auto k18_vps = MakeK18VPs();
  tpc_ana.TrackSearchTPCHelix(k18_vps);
  FillK18HelixTracks(tpc_ana);
  HF1("Status", event.status++);
#endif

  event.ntTpc = tpc_ana.GetNTracksTPCHelix();
  HF1("NTracks_TPC", event.ntTpc);
  FillHSHelixResidualMatch(tpc_ana);
  if (event.ntTpc == 0)
    return true;
  event.effective_ntTpc = 0;
  event.resizeTracks(event.ntTpc);

  for (Int_t it = 0; it < event.ntTpc; ++it) {
    TPCLocalTrackHelix* helix_track = tpc_ana.GetTrackTPCHelix(it);
    ProcessOneHelixTrack(it, event.ntTpc, helix_track, tpc_ana, event_ana);
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
  if (!dst::SetupReader(kK18, "kK18")) return false;

  const auto& reader = TTreeReaderCont[kTpcHit];

  dst::SetBranch(reader, "run_number",   src.runnum);
  dst::SetBranch(reader, "event_number", src.evnum);
  dst::SetBranch(reader, "trig_pat",     src.trigpat);
  dst::SetBranch(reader, "trig_flag",    src.trigflag);
  dst::SetBranch(reader, "beam_flag",    src.beamflag);
  dst::SetBranch(reader, "npadTpc",      src.npadTpc);
  dst::SetBranch(reader, "nhTpc",        src.nhTpc);
  dst::SetBranch(reader, "layerTpc",     src.layerTpc);
  dst::SetBranch(reader, "rowTpc",       src.rowTpc);
  dst::SetBranch(reader, "padTpc",       src.padTpc);
  dst::SetBranch(reader, "pedTpc",       src.pedTpc);
  dst::SetBranch(reader, "rmsTpc",       src.rmsTpc);
  dst::SetBranch(reader, "deTpc",        src.deTpc);
  dst::SetBranch(reader, "tTpc",         src.tTpc);

  dst::SetBranch(reader, "chisqrTpc",    src.chisqrTpc);
  dst::SetBranch(reader, "clkTpc",       src.clkTpc);

  const auto& k18_reader = TTreeReaderCont[kK18];
  dst::SetBranch(k18_reader, "run_number",   src.k18_runnum);
  dst::SetBranch(k18_reader, "event_number", src.k18_evnum);
  dst::SetBranch(k18_reader, "rk_statusHS",  src.rkStatusHS);
  dst::SetBranch(k18_reader, "xvpHS",        src.xvpHS);
  dst::SetBranch(k18_reader, "yvpHS",        src.yvpHS);
  dst::SetBranch(k18_reader, "zvpHS",        src.zvpHS);
  dst::SetBranch(k18_reader, "xout", src.xoutK18);
  dst::SetBranch(k18_reader, "yout", src.youtK18);
  dst::SetBranch(k18_reader, "uout", src.uoutK18);
  dst::SetBranch(k18_reader, "vout", src.voutK18);
  dst::SetBranch(k18_reader, "pk18", src.pK18);

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
  hist::BuildTPCHSHelixResidual();

  tree = new TTree("tpc", "tree of DstTPCHelixTracking");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("clkTpc", &event.clkTpc);

#if RawHit
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

#if RawCluster
  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x );
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row_center", &event.cluster_row_center );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_x_center", &event.cluster_x_center );
  tree->Branch( "cluster_y_center", &event.cluster_y_center );
  tree->Branch( "cluster_z_center", &event.cluster_z_center );
#endif

  tree->Branch( "ntK18", &event.ntK18 );
  tree->Branch( "nh_k18", &event.nh_k18 );
  tree->Branch( "x_hit_k18", &event.x_hit_k18 );
  tree->Branch( "y_hit_k18", &event.y_hit_k18 );
  tree->Branch( "z_hit_k18", &event.z_hit_k18 );
  tree->Branch( "x_cal_k18", &event.x_cal_k18 );
  tree->Branch( "y_cal_k18", &event.y_cal_k18 );
  tree->Branch( "z_cal_k18", &event.z_cal_k18 );
  tree->Branch( "helix_cx_k18", &event.helix_cx_k18 );
  tree->Branch( "helix_cy_k18", &event.helix_cy_k18 );
  tree->Branch( "helix_z0_k18", &event.helix_z0_k18 );
  tree->Branch( "helix_r_k18", &event.helix_r_k18 );
  tree->Branch( "helix_dz_k18", &event.helix_dz_k18 );
  tree->Branch( "helix_theta_min_k18", &event.helix_theta_min_k18 );
  tree->Branch( "helix_theta_max_k18", &event.helix_theta_max_k18 );
  tree->Branch( "mom0_x_k18", &event.mom0_x_k18 );
  tree->Branch( "mom0_y_k18", &event.mom0_y_k18 );
  tree->Branch( "mom0_z_k18", &event.mom0_z_k18 );
  tree->Branch( "mom0_k18", &event.mom0_k18 );
  tree->Branch( "xout_k18", &event.xout_k18 );
  tree->Branch( "yout_k18", &event.yout_k18 );
  tree->Branch( "uout_k18", &event.uout_k18 );
  tree->Branch( "vout_k18", &event.vout_k18 );
  tree->Branch( "p_k18", &event.p_k18 );

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

#if TruncatedMean
  tree->Branch( "dEdx_0", &event.dEdx_0 );
  tree->Branch( "dEdx_10", &event.dEdx_10 );
  tree->Branch( "dEdx_20", &event.dEdx_20 );
  tree->Branch( "dEdx_30", &event.dEdx_30 );
  tree->Branch( "dEdx_40", &event.dEdx_40 );
  tree->Branch( "dEdx_50", &event.dEdx_50 );
  tree->Branch( "dEdx_60", &event.dEdx_60 );
  tree->Branch( "dEdx_cor_0", &event.dEdx_cor_0 );
  tree->Branch( "dEdx_cor_10", &event.dEdx_cor_10 );
  tree->Branch( "dEdx_cor_20", &event.dEdx_cor_20 );
  tree->Branch( "dEdx_cor_30", &event.dEdx_cor_30 );
  tree->Branch( "dEdx_cor_40", &event.dEdx_cor_40 );
  tree->Branch( "dEdx_cor_50", &event.dEdx_cor_50 );
  tree->Branch( "dEdx_cor_60", &event.dEdx_cor_60 );
#endif
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );
#if TrackSearch
  tree->Branch( "combi_id", &event.combi_id );
  tree->Branch( "closeDistTpc", &event.closeDistTpc );
  tree->Branch( "vtxTpc", &event.vtxTpc );
  tree->Branch( "vtyTpc", &event.vtyTpc );
  tree->Branch( "vtzTpc", &event.vtzTpc );
  tree->Branch( "mom_vtx", &event.mom_vtx );
  tree->Branch( "mom_vty", &event.mom_vty );
  tree->Branch( "mom_vtz", &event.mom_vtz );
#endif

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
#if TrackCluster
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
  tree->Branch( "transverse_diffusion_dist", &event.transverse_diffusion_dist );
  tree->Branch( "transverse_diffusion_ratio", &event.transverse_diffusion_ratio );
  tree->Branch( "transverse_diffusion_drift", &event.transverse_diffusion_drift );
  tree->Branch( "transverse_diffusion_track_hit_index", &event.transverse_diffusion_track_hit_index );
#endif



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

