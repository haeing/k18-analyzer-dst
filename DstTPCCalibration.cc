// -*- C++ -*-

#include <iomanip>
#include <iostream>

#include <TPDGCode.h>

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
  };
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpcHit,  kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCHit]",  "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "" };
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

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> nhrawtrack; //Number of Raw Hits (in 1 tracks)
  std::vector<Int_t> is_beam; // 1 = Beam, 0 = Scat
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
  
  std::vector<std::vector<Double_t>> track_hit_pad;
  std::vector<std::vector<Double_t>> track_hit_de;
  std::vector<std::vector<Double_t>> track_hit_ctTpc;
  
  
  
  
  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc);
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
    dst::clear_all(
		   nhtrack, nhrawtrack, is_beam, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz,
      helix_theta_min, helix_theta_max, dE, dEdx,
      
#if TruncatedMean
      dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, 
      dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60,
#endif
      dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, pid, path,
      combi_id, closeDistTpc, vtxTpc, vtyTpc, vtzTpc, mom_vtx, mom_vty, mom_vtz,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center, track_hit_pad, track_hit_de, track_hit_ctTpc
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
  }

  void resizeTracks(Int_t nTracks) {
    dst::resize_all(nTracks,
		    nhtrack, nhrawtrack, is_beam, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz,
      helix_theta_min, helix_theta_max, dE, dEdx,

#if TruncatedMean
      dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, 
      dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60,
#endif
      dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, pid, path,
      combi_id, closeDistTpc, vtxTpc, vtyTpc, vtzTpc, mom_vtx, mom_vty, mom_vtz,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center, track_hit_pad, track_hit_de, track_hit_ctTpc
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
    Int_t n_rawhits = 0;
    // Per-track basic parameters and containers.
    Int_t n_hits = helix_track->GetNHit();
    Double_t chi_sqr = helix_track->GetChiSquare();
    Double_t helix_cx = helix_track->Getcx(), helix_cy = helix_track->Getcy();
    Double_t helix_z0 = helix_track->Getz0(), helix_r = helix_track->Getr();
    Double_t helix_dz = helix_track->Getdz();
    TVector3 mom0_vec = helix_track->GetMom0();
    Int_t is_beam = helix_track->GetIsBeam();

    event.nhtrack[it] = n_hits;
    event.is_beam[it] = is_beam;
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

      const auto &raw_hit_cont = cl->GetHitContainer();
      for(int j=0;j<raw_hit_cont.size();j++){
	TPCHit *raw_hit = raw_hit_cont[j];
	if(!raw_hit)continue;
	
	event.track_hit_pad[it].push_back(raw_hit->GetPad());
	event.track_hit_de[it].push_back(raw_hit->GetCDe());
	event.track_hit_ctTpc[it].push_back(raw_hit->GetCTime());
	n_rawhits++;
	
      }
      
      
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
    event.nhrawtrack[it] = n_rawhits;
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

  if(ievent%1000==0){
    std::cout << "#D Event Number: " << std::setw(6)
	      << ievent << std::endl;
  }

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
  if (!CheckEntries(TTreeCont)) return false;

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
  tpc_ana.TrackSearchTPCHelix();
  HF1("Status", event.status++);
#endif

  event.ntTpc = tpc_ana.GetNTracksTPCHelix();
  HF1("NTracks_TPC", event.ntTpc);
  if (event.ntTpc == 0)
    return true;
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

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "is_beam", &event.is_beam );
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

  tree->Branch( "nhrawtrack", &event.nhrawtrack);
  tree->Branch( "track_hit_pad", &event.track_hit_pad);
  tree->Branch( "track_hit_de", &event.track_hit_de);
  tree->Branch( "track_hit_ctTpc", &event.track_hit_ctTpc);
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

