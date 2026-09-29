// -*- C++ -*-
// To Do: consider BcOut nhit (not important)

#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TVector3.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "RootHelper.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>

#define Exclusive 1
#define CalibHist 0 // enable/disable per-pad calibration histograms

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;
  const auto& gUnpacker = GUnpacker::get_instance();
  auto&       gConf     = ConfMan::GetInstance();
  const auto& gGeom     = DCGeomMan::GetInstance();
  const auto& gUser     = UserParamMan::GetInstance();
  const auto& gTpcParam = TPCParamMan::GetInstance();
  const auto& gCounter  = debug::ObjectCounter::GetInstance();
  const Double_t MAX_CHISQR_TPC = 50.0;
  const Double_t MAX_CHISQR_BCOUT = 5.0;

  // No TPC-specific UserParam keys are required for DstTPCBcOutTracking,
  // but keep an explicit (empty) manifest for consistency with other Dst modules.
  const std::vector<TString> kUserParamKeys = {};
}

namespace dst
{
  enum kArgc
  {
    kProcess, kConfFile,
    kTpcTracking, kBcOut, kHodo, kOutFile, nArgc
  };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCTracking]", "[BcOut]", "[Hodo]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "bcout", "hodo", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
  std::vector<UInt_t> evnumPerFile;
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
  std::vector<Double_t> cobo_id;

  // TPC Tracking Info
  Int_t ntTpc;
  std::vector<Int_t> nhtrack;
  std::vector<Double_t> chisqrTpc, x0Tpc, y0Tpc, u0Tpc, v0Tpc, theta;
  std::vector<Double_t> dE, dEdx;
    
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x, hitpos_y, hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x, calpos_y, calpos_z;
  std::vector<std::vector<Double_t>> residual, residual_x, residual_y, residual_z;
  std::vector<std::vector<Double_t>> residual_horizontal, residual_vertical;
  std::vector<std::vector<Double_t>> resolution_x, resolution_y, resolution_z;
  std::vector<std::vector<Double_t>> resolution_horizontal, resolution_vertical;
  std::vector<std::vector<Double_t>> pathhit, theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de, track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  std::vector<std::vector<Double_t>> exresidual, exresidual_x, exresidual_y, exresidual_z;
  std::vector<std::vector<Double_t>> exresidual_horizontal, exresidual_vertical;

  Int_t ntTpc_inside;
  Double_t prodvtx_x, prodvtx_y, prodvtx_z;

  // BcOut Info
  Int_t ntBcOut;
  // std::vector<Int_t> nhBcOut;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut, y0BcOut, u0BcOut, v0BcOut;
  std::vector<Double_t> xtgtBcOut, ytgtBcOut, utgtBcOut, vtgtBcOut;

  // Hodo Info
  Double_t btof;
  Double_t ftof;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc, cobo_id);
  }
  
  void clearTPCTracks() {
    ntTpc = 0;
    dst::clear_all(
      nhtrack, chisqrTpc, 
      x0Tpc, y0Tpc, u0Tpc, v0Tpc, 
      theta, dE, dEdx,
      
      hitlayer, 
      hitpos_x, hitpos_y, hitpos_z, 
      calpos_x, calpos_y, calpos_z,
      
      residual, 
      residual_x, residual_y, residual_z, 
      residual_horizontal, residual_vertical,
      
      resolution_x, resolution_y, resolution_z, 
      resolution_horizontal, resolution_vertical,
      
      pathhit, theta_diff, track_cluster_de, track_cluster_size,
      track_cluster_mrow, track_cluster_de_center,
      track_cluster_x_center,       track_cluster_y_center, track_cluster_z_center,
      track_cluster_row_center,

      exresidual,
      exresidual_x, exresidual_y, exresidual_z, 
      exresidual_horizontal, exresidual_vertical
    );
  }

  void clearTPCVertex() {
    ntTpc_inside = 0;
    prodvtx_x = TMath::QuietNaN(); 
    prodvtx_y = TMath::QuietNaN(); 
    prodvtx_z = TMath::QuietNaN();
  }

  void clearBcOut() {
    ntBcOut = 0;
    dst::clear_all(
      // nhBcOut, 
      chisqrBcOut, 
      x0BcOut, y0BcOut, u0BcOut, v0BcOut, 
      xtgtBcOut, ytgtBcOut, utgtBcOut, vtgtBcOut
    );
  }

  void clear() {
    clearBasicInfo();
    clearTPCTracks();
    clearTPCVertex();
    clearBcOut();
    btof = TMath::QuietNaN();
    ftof = TMath::QuietNaN();
  }

};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum_tpc;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;
  TTreeReaderValue<std::vector<Double_t>>* cobo_id;
  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc;
  TTreeReaderValue<std::vector<Double_t>>* x0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* y0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* u0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* v0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* theta;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;

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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* theta_diff;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_vertical;

  TTreeReaderValue<Int_t>* ntTpc_inside;
  TTreeReaderValue<Double_t>* prodvtx_x;
  TTreeReaderValue<Double_t>* prodvtx_y;
  TTreeReaderValue<Double_t>* prodvtx_z;

  TTreeReaderValue<Double_t>* btof;
  TTreeReaderValue<Double_t>* ftof;

  TTreeReaderValue<UInt_t>* evnum_bcout;
  TTreeReaderValue<UInt_t>* evnum_hodo;

  TTreeReaderValue<Int_t>* ntBcOut;
  // TTreeReaderValue<std::vector<Int_t>>* nhBcOut;
  TTreeReaderValue<std::vector<Double_t>>* chisqrBcOut;
  TTreeReaderValue<std::vector<Double_t>>* x0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* y0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* u0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* v0BcOut;
};

//_____________________________________________________________________________
namespace root
{
  Event event;
  Src   src;
  TTree *tree;
}

//_____________________________________________________________________________
namespace
{
  using namespace root;
  using namespace dst;

  //_____________________________________________________________________________
  void CopyTPCTrackingData(Event& event, const Src& src)
  {
    event.ntTpc     = **src.ntTpc;
    event.nhtrack   = **src.nhtrack;
    event.chisqrTpc = **src.chisqrTpc;
    event.x0Tpc     = **src.x0Tpc;
    event.y0Tpc     = **src.y0Tpc;
    event.u0Tpc     = **src.u0Tpc;
    event.v0Tpc     = **src.v0Tpc;
    event.theta     = **src.theta;
    event.dE        = **src.dE;
    event.dEdx      = **src.dEdx;

    event.hitlayer   = **src.hitlayer;
    event.hitpos_x   = **src.hitpos_x;
    event.hitpos_y   = **src.hitpos_y;
    event.hitpos_z   = **src.hitpos_z;
    event.calpos_x   = **src.calpos_x;
    event.calpos_y   = **src.calpos_y;
    event.calpos_z   = **src.calpos_z;
    event.residual   = **src.residual;
    event.residual_x = **src.residual_x;
    event.residual_y = **src.residual_y;
    event.residual_z = **src.residual_z;
    event.residual_horizontal = **src.residual_horizontal;
    event.residual_vertical   = **src.residual_vertical;
    event.resolution_x        = **src.resolution_x;
    event.resolution_y        = **src.resolution_y;
    event.resolution_z        = **src.resolution_z;
    event.resolution_horizontal = **src.resolution_horizontal;
    event.resolution_vertical   = **src.resolution_vertical;
    event.pathhit             = **src.pathhit;
    event.theta_diff          = **src.theta_diff;
    event.track_cluster_de    = **src.track_cluster_de;
    event.track_cluster_size  = **src.track_cluster_size;
    event.track_cluster_mrow       = **src.track_cluster_mrow;
    event.track_cluster_de_center  = **src.track_cluster_de_center;
    event.track_cluster_x_center   = **src.track_cluster_x_center;
    event.track_cluster_y_center   = **src.track_cluster_y_center;
    event.track_cluster_z_center   = **src.track_cluster_z_center;
    event.track_cluster_row_center = **src.track_cluster_row_center;

  #if Exclusive
    event.exresidual   = **src.exresidual;
    event.exresidual_x = **src.exresidual_x;
    event.exresidual_y = **src.exresidual_y;
    event.exresidual_z = **src.exresidual_z;
    event.exresidual_horizontal = **src.exresidual_horizontal;
    event.exresidual_vertical   = **src.exresidual_vertical;
  #endif

    event.ntTpc_inside = **src.ntTpc_inside;
    event.prodvtx_x    = **src.prodvtx_x;
    event.prodvtx_y    = **src.prodvtx_y;
    event.prodvtx_z    = **src.prodvtx_z;
  }

  //_____________________________________________________________________________
  void CopyBcOutData(Event& event, const Src& src)
  {
    static const Double_t ztgt_global = gGeom.GetGlobalPosition("Target").z();

    const Int_t nt_bcout = **src.ntBcOut;
    event.ntBcOut = nt_bcout;
    event.chisqrBcOut = **src.chisqrBcOut;
    event.x0BcOut     = **src.x0BcOut;
    event.y0BcOut     = **src.y0BcOut;
    event.u0BcOut     = **src.u0BcOut;
    event.v0BcOut     = **src.v0BcOut;
    event.utgtBcOut   = **src.u0BcOut;
    event.vtgtBcOut   = **src.v0BcOut;

    event.xtgtBcOut.reserve(nt_bcout);
    event.ytgtBcOut.reserve(nt_bcout);
    for (Int_t it=0; it<nt_bcout; ++it) {
      event.xtgtBcOut.push_back((**src.u0BcOut)[it]*ztgt_global + (**src.x0BcOut)[it]);
      event.ytgtBcOut.push_back((**src.v0BcOut)[it]*ztgt_global + (**src.y0BcOut)[it]);
    }
  }

  //_____________________________________________________________________________
  void SelectBestTrackIndices(Int_t& best_tpc_idx, Int_t& best_bcout_idx)
  {
    best_tpc_idx   = 0;
    best_bcout_idx = 0;
    // NOTE: chisqrTpc is already reduced-chi2-like (chi2/(dof-4)) in TPCLocalTrack.
    if (event.ntTpc > 1) {
      Double_t min_chisqr_tpc = event.chisqrTpc[0];
      for (Int_t it = 1; it < event.ntTpc; ++it) {
        if (event.chisqrTpc[it] < min_chisqr_tpc) {
          min_chisqr_tpc = event.chisqrTpc[it];
          best_tpc_idx = it;
        }
      }
    }
    // NOTE: chisqr in BcOut tree is reduced-chi2-like value, so select closest to 1.0.
    if (event.ntBcOut > 1) {
      Double_t best_score = 1.0e+10;
      Double_t best_chisqr_bcout = 1.0e+10;
      best_bcout_idx = -1;
      for (Int_t it = 0; it < event.ntBcOut; ++it) {
        const Double_t chisqr = event.chisqrBcOut[it];
        if (!std::isfinite(chisqr) || chisqr <= 0.) continue;
        const Double_t score = TMath::Abs(chisqr - 1.0);
        if (score < best_score || (score == best_score && chisqr < best_chisqr_bcout)) {
          best_score = score;
          best_chisqr_bcout = chisqr;
          best_bcout_idx = it;
        }
      }
      if (best_bcout_idx < 0) best_bcout_idx = 0;
    }
  }

  //_____________________________________________________________________________
  void FillTPCTrkSummaryHist(TPCEventAnalyzer& event_ana)
  {
    HF1("TPCTrk_Num_Track", event.ntTpc);
    for (Int_t it = 0; it < event.ntTpc; ++it) {
      event_ana.FillTrkParamHist(
        event.nhtrack[it], event.chisqrTpc[it],
        event.x0Tpc[it], event.y0Tpc[it], event.u0Tpc[it], event.v0Tpc[it]
      );
    }
  }

  //_____________________________________________________________________________
  void FillBcOutSummaryHist(TPCEventAnalyzer& event_ana)
  {
    HF1("BcOut_Num_Track", event.ntBcOut);
    for (Int_t it = 0; it < event.ntBcOut; ++it) {
      event_ana.FillBcOutTrackHist(
        event.chisqrBcOut[it],
        event.x0BcOut[it], event.y0BcOut[it], event.u0BcOut[it], event.v0BcOut[it],
        event.xtgtBcOut[it], event.ytgtBcOut[it], event.utgtBcOut[it], event.vtgtBcOut[it]
      );
    }
  }

  //_____________________________________________________________________________
  std::array<Double_t, 4> CalcBcOutAtTpc(Int_t best_bcout_idx)
  {
    std::array<Double_t, 4> xyuv_bcout_tpc_coor = {};
    static const Double_t ztgt_global = gGeom.GetGlobalPosition("Target").z();
    static const Double_t ra1         = gGeom.GetRotAngle1("HypTPC");
    static const Double_t ra2         = gGeom.GetRotAngle2("HypTPC");
    static const Double_t tan_ra1     = TMath::Tan(ra1*TMath::DegToRad());
    static const Double_t tan_ra2     = TMath::Tan(ra2*TMath::DegToRad());

    ThreeVector global_pos(
      (**src.u0BcOut)[best_bcout_idx]*ztgt_global + (**src.x0BcOut)[best_bcout_idx],
      (**src.v0BcOut)[best_bcout_idx]*ztgt_global + (**src.y0BcOut)[best_bcout_idx],
      ztgt_global
    );

    ThreeVector local_pos = gGeom.Global2LocalPos("HypTPC", global_pos);
    xyuv_bcout_tpc_coor[0] = local_pos.x();
    xyuv_bcout_tpc_coor[1] = local_pos.y();
    xyuv_bcout_tpc_coor[2] = ((**src.u0BcOut)[best_bcout_idx] - tan_ra2)
                            / (1. + (**src.u0BcOut)[best_bcout_idx]*tan_ra2);
    xyuv_bcout_tpc_coor[3] = ((**src.v0BcOut)[best_bcout_idx] + tan_ra1)
                            / (1. - (**src.v0BcOut)[best_bcout_idx]*tan_ra1);
    return xyuv_bcout_tpc_coor;
  }

  //_____________________________________________________________________________
  Int_t GetTrackCenterRow(Int_t it, Int_t ih)
  {
    if (it < static_cast<Int_t>(event.track_cluster_row_center.size())
        && ih < static_cast<Int_t>(event.track_cluster_row_center[it].size())) {
      return static_cast<Int_t>(event.track_cluster_row_center[it][ih]);
    }
    return -1;
  }

  //_____________________________________________________________________________
  void FillTPCBcOutResidualPull(Int_t best_bcout_idx, TPCEventAnalyzer& event_ana)
  {
    for (Int_t it = 0; it < event.ntTpc; ++it) {
      for (Int_t ih = 0; ih < event.nhtrack[it]; ++ih) {
        const Int_t layer = static_cast<Int_t>(event.hitlayer[it][ih]);
        const TVector3 resolution_vec(
          event.resolution_x[it][ih],
          event.resolution_y[it][ih],
          event.resolution_z[it][ih]
        );
        const Bool_t valid_tpc_resolution = !tpc::IsDummyResolutionVec(resolution_vec);

        ThreeVector local_pos_tpc(event.hitpos_x[it][ih], event.hitpos_y[it][ih], event.hitpos_z[it][ih]);
        ThreeVector global_pos_tpc = gGeom.Local2GlobalPos("HypTPC", local_pos_tpc);
        const Double_t z_bcout = global_pos_tpc.z();
        const Double_t x_bcout = (**src.u0BcOut)[best_bcout_idx] * z_bcout + (**src.x0BcOut)[best_bcout_idx];
        const Double_t y_bcout = (**src.v0BcOut)[best_bcout_idx] * z_bcout + (**src.y0BcOut)[best_bcout_idx];
        const Double_t res_x_bcout = global_pos_tpc.x() - x_bcout;
        const Double_t res_y_bcout = global_pos_tpc.y() - y_bcout;

        const Int_t center_row = GetTrackCenterRow(it, ih);
        const TVector3 trk_res_global(event.residual_x[it][ih], event.residual_y[it][ih], event.residual_z[it][ih]);
        const TVector3 trk_res_local(event.residual_horizontal[it][ih], event.residual_vertical[it][ih], 0.0);
        const TVector3 trk_pull_global = valid_tpc_resolution
          ? TVector3(
            event.residual_x[it][ih] / event.resolution_x[it][ih],
            event.residual_y[it][ih] / event.resolution_y[it][ih],
            event.residual_z[it][ih] / event.resolution_z[it][ih]
          )
          : TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
        const TVector3 trk_pull_local = valid_tpc_resolution
          ? TVector3(
            event.residual_horizontal[it][ih] / event.resolution_horizontal[it][ih],
            event.residual_vertical[it][ih] / event.resolution_vertical[it][ih],
            0.0
          )
          : TVector3(TMath::QuietNaN(), TMath::QuietNaN(), 0.0);
        const TVector3 cl_res(res_x_bcout, res_y_bcout, 0.0);

        event_ana.FillTPCBcOutTrackingResidualPullHist(
          layer, center_row, valid_tpc_resolution,
          trk_res_global, trk_res_local, trk_pull_global, trk_pull_local,
          global_pos_tpc.x(), global_pos_tpc.y(), y_bcout, cl_res
        );

        if (center_row < 0) {
          continue;
        }
        const Double_t ctime_trk = 0.0;
        ThreeVector local_trk(event.calpos_x[it][ih], event.calpos_y[it][ih], event.calpos_z[it][ih]);
        ThreeVector global_trk = gGeom.Local2GlobalPos("HypTPC", local_trk);
        const TVector3 local_pos_hit(event.hitpos_x[it][ih], event.hitpos_y[it][ih], event.hitpos_z[it][ih]);
        event_ana.FillCoBoClockTime("TPCTrk", layer, center_row, ctime_trk, local_pos_hit, global_trk.y());
      }
    }
  }

  //_____________________________________________________________________________
  void FillTPCBcOutCorrelation(
    Int_t best_tpc_idx,
    Int_t best_bcout_idx,
    TPCEventAnalyzer& event_ana)
  {
    FillTPCTrkSummaryHist(event_ana);
    FillBcOutSummaryHist(event_ana);
    const auto xyuv_bcout_tpc_coor = CalcBcOutAtTpc(best_bcout_idx);

    event_ana.FillTPCBcOutTgtResidualHist(
      event.x0Tpc[best_tpc_idx], event.y0Tpc[best_tpc_idx],
      event.u0Tpc[best_tpc_idx], event.v0Tpc[best_tpc_idx],
      xyuv_bcout_tpc_coor[0], xyuv_bcout_tpc_coor[1],
      xyuv_bcout_tpc_coor[2], xyuv_bcout_tpc_coor[3]
    );

    FillTPCBcOutResidualPull(best_bcout_idx, event_ana);
  }

} // namespace

//_____________________________________________________________________________
int
main(int argc, char **argv)
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

  std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  DstClose();
  return EXIT_SUCCESS;
}  // main

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
  if (ievent%1000==0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  evnumPerFile = { **src.evnum_tpc, **src.evnum_bcout, **src.evnum_hodo };
  if (!dst::CheckEventNumbers(evnumPerFile, ievent, 
                            {TreeName[kTpcTracking], TreeName[kBcOut], TreeName[kHodo]})) {
    return false;
  }

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum_tpc;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;
  event.cobo_id  = **src.cobo_id;
  HF1("Status", event.status++);

  if (**src.ntTpc == 0) return true;
  HF1("Status", event.status++);

  if (!TPCEventAnalyzer::ValidateCoboClocks(event.clkTpc))
    return true;
  HF1("Status", event.status++);

  CopyTPCTrackingData(event, src);
  static TPCEventAnalyzer event_ana;
  event_ana.SetClock(event.clkTpc);
  HF1("Status", event.status++);

  CopyBcOutData(event, src);
  HF1("Status", event.status++);

  event.btof = **src.btof;
  event.ftof = **src.ftof;
  HF1("Status", event.status++);

  Int_t best_tpc_idx = 0;
  Int_t best_bcout_idx = 0;
  SelectBestTrackIndices(best_tpc_idx, best_bcout_idx);

  if (event.ntTpc == 0 || event.ntBcOut == 0) return true;
  HF1("Status", event.status++);

  if (event.chisqrTpc[best_tpc_idx] > MAX_CHISQR_TPC
     || event.chisqrBcOut[best_bcout_idx] > MAX_CHISQR_BCOUT) return true;
  HF1("Status", event.status++);

  FillTPCBcOutCorrelation(best_tpc_idx, best_bcout_idx, event_ana);
  HF1("Status", event.status++);

  return true;
}  // dst::DstRead


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
  // -------------------------------------------------------
  // TPCTracking
  // -------------------------------------------------------
  if (!dst::SetupReader(kTpcTracking, "kTpcTracking")) return false;

  // Explicitly using TTreeReaderCont[kTpcTracking]
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "run_number",   src.runnum);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "event_number", src.evnum_tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "trig_pat",     src.trigpat);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "trig_flag",    src.trigflag);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "beam_flag",    src.beamflag);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "clkTpc",       src.clkTpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "cobo_id",      src.cobo_id);

  dst::SetBranch(TTreeReaderCont[kTpcTracking], "ntTpc",        src.ntTpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "nhtrack",      src.nhtrack);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "chisqrTpc",    src.chisqrTpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "x0Tpc",        src.x0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "y0Tpc",        src.y0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "u0Tpc",        src.u0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "v0Tpc",        src.v0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "theta",        src.theta);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "dE",           src.dE);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "dEdx",         src.dEdx);

  // Hit & Cal Pos
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitlayer",     src.hitlayer);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_x",     src.hitpos_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_y",     src.hitpos_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_z",     src.hitpos_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_x",     src.calpos_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_y",     src.calpos_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_z",     src.calpos_z);

  // Residuals
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual",            src.residual);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_x",          src.residual_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_y",          src.residual_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_z",          src.residual_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_horizontal", src.residual_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_vertical",   src.residual_vertical);

  // Resolution
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_x",          src.resolution_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_y",          src.resolution_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_z",          src.resolution_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_horizontal", src.resolution_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_vertical",   src.resolution_vertical);

  // Track Clusters
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "pathhit",                  src.pathhit);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "theta_diff",               src.theta_diff);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_de",         src.track_cluster_de);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_size",       src.track_cluster_size);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_mrow",       src.track_cluster_mrow);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_de_center",  src.track_cluster_de_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_x_center",   src.track_cluster_x_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_y_center",   src.track_cluster_y_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_z_center",   src.track_cluster_z_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_row_center", src.track_cluster_row_center);

  // ExResiduals
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual",            src.exresidual);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_x",          src.exresidual_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_y",          src.exresidual_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_z",          src.exresidual_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_horizontal", src.exresidual_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_vertical",   src.exresidual_vertical);

  dst::SetBranch(TTreeReaderCont[kTpcTracking], "ntTpc_target", src.ntTpc_inside);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_x",    src.prodvtx_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_y",    src.prodvtx_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_z",    src.prodvtx_z);

  // -------------------------------------------------------
  // Hodoscope
  // -------------------------------------------------------
  if (!dst::SetupReader(kHodo, "kHodo")) return false;

  dst::SetBranch(TTreeReaderCont[kHodo], "event_number", src.evnum_hodo);
  dst::SetBranch(TTreeReaderCont[kHodo], "btof0", src.btof);
  dst::SetBranch(TTreeReaderCont[kHodo], "ftof0", src.ftof);

  // -------------------------------------------------------
  // BcOut
  // -------------------------------------------------------
  if (!dst::SetupReader(kBcOut, "kBcOut")) return false;

  dst::SetBranch(TTreeReaderCont[kBcOut], "event_number", src.evnum_bcout);
  dst::SetBranch(TTreeReaderCont[kBcOut], "ntrack", src.ntBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "chisqr", src.chisqrBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "x0",     src.x0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "y0",     src.y0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "u0",     src.u0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "v0",     src.v0BcOut);

  return true;
}  // dst::SetupReaders

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
  hist::BuildTPCTrackingCommon();
  hist::BuildTPCLineTrackParam();
  hist::BuildTPCBcOutTracking(TPCEventAnalyzer::GetDstCalibFlag());

  tree = new TTree("tpc", "tree of DstTPCBcOutTracking");
  tree->Branch("status",       &event.status);
  tree->Branch("run_number",   &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat",     &event.trigpat);
  tree->Branch("trig_flag",    &event.trigflag);
  tree->Branch("beam_flag",    &event.beamflag);
  tree->Branch("clkTpc",       &event.clkTpc);
  tree->Branch("cobo_id",      &event.cobo_id);

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("x0Tpc", &event.x0Tpc);
  tree->Branch("y0Tpc", &event.y0Tpc);
  tree->Branch("u0Tpc", &event.u0Tpc);
  tree->Branch("v0Tpc", &event.v0Tpc);
  tree->Branch("theta", &event.theta);
  tree->Branch("hitlayer", &event.hitlayer);
  tree->Branch("hitpos_x", &event.hitpos_x);
  tree->Branch("hitpos_y", &event.hitpos_y);
  tree->Branch("hitpos_z", &event.hitpos_z);
  tree->Branch("calpos_x", &event.calpos_x);
  tree->Branch("calpos_y", &event.calpos_y);
  tree->Branch("calpos_z", &event.calpos_z);
  tree->Branch("residual", &event.residual);
  tree->Branch("residual_x", &event.residual_x);
  tree->Branch("residual_y", &event.residual_y);
  tree->Branch("residual_z", &event.residual_z);
  tree->Branch("residual_horizontal", &event.residual_horizontal);
  tree->Branch("residual_vertical", &event.residual_vertical);
  tree->Branch("resolution_x", &event.resolution_x);
  tree->Branch("resolution_y", &event.resolution_y);
  tree->Branch("resolution_z", &event.resolution_z);
  tree->Branch("resolution_horizontal", &event.resolution_horizontal);
  tree->Branch("resolution_vertical", &event.resolution_vertical);
  tree->Branch("dE", &event.dE );
  tree->Branch("dEdx", &event.dEdx );
  tree->Branch("pathhit", &event.pathhit);
  tree->Branch("theta_diff", &event.theta_diff);
  tree->Branch("track_cluster_de", &event.track_cluster_de);
  tree->Branch("track_cluster_size", &event.track_cluster_size);
  tree->Branch("track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch("track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch("track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch("track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch("track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch("track_cluster_row_center", &event.track_cluster_row_center);

#if Exclusive
  tree->Branch("exresidual", &event.exresidual);
  tree->Branch("exresidual_x", &event.exresidual_x);
  tree->Branch("exresidual_y", &event.exresidual_y);
  tree->Branch("exresidual_z", &event.exresidual_z);
  tree->Branch("exresidual_horizontal", &event.exresidual_horizontal);
  tree->Branch("exresidual_vertical", &event.exresidual_vertical);
#endif

  tree->Branch("btof", &event.btof);
  tree->Branch("ftof", &event.ftof);

  tree->Branch("ntBcOut", &event.ntBcOut);
  // tree->Branch("nhBcOut", &event.nhBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);
  tree->Branch("xtgtBcOut", &event.xtgtBcOut);
  tree->Branch("ytgtBcOut", &event.ytgtBcOut);
  tree->Branch("utgtBcOut", &event.utgtBcOut);
  tree->Branch("vtgtBcOut", &event.vtgtBcOut);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
