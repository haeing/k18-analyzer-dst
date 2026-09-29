// -*- C++ -*-

#include <iomanip>
#include <iostream>

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
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>

#define RawHit 0
#define RawCluster 1
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
  const auto& gTpcParam = TPCParamMan::GetInstance();
  const auto& gCounter  = debug::ObjectCounter::GetInstance();
  const Double_t MAX_RESIDUAL = 20.0; // Residual gate [mm] to accept hits/clusters close to the reference track.

  const std::vector<TString> kUserParamKeys = {
    // Cluster building (ReCalcTPCHits / MakeUpTPCClusters)
    "MinCDeTPC", "MaxYDifClusterTPC",
    "MinClusterDeTPC", "MinClusterSizeTPC",
    "MinClusterYPosTPC", "MaxClusterYPosTPC",

    // Optional parameter (default value is provided in the code)
    // "MaxCenterRowDiffTPC",
  };
}

namespace dst
{
  enum kArgc
  {
    kProcess, kConfFile,
    kTpcHit, kBcOut, kOutFile, nArgc
  };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCHit]", "[BcOut]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "bcout", "" };
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

  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

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

  // BcOut Info
  Int_t ntBcOut;
  Int_t best_bcout_id;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut;
  std::vector<Double_t> y0BcOut;
  std::vector<Double_t> u0BcOut;
  std::vector<Double_t> v0BcOut;

  // Residual Info (BcOut reference)
  std::vector<Double_t> residual_x_cluster, residual_y_cluster;
  std::vector<Double_t> residual_x_hit, residual_y_hit;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc, cobo_id);
  }

  void clearRawHits() {
    nhTpc = 0;
    dst::clear_all(
      raw_hitpos_x, raw_hitpos_y, raw_hitpos_z,
      raw_de, raw_padid, raw_layer, raw_row
    );
  }

  void clearClusters() {
    nclTpc = 0;
    dst::clear_all(
      cluster_x, cluster_y, cluster_z, cluster_de,
      cluster_size, cluster_layer, cluster_mrow,
      cluster_de_center, cluster_x_center, cluster_y_center,
      cluster_z_center, cluster_row_center
    );
  }

  void clearBcOut() {
    ntBcOut = 0;
    best_bcout_id = -1;
    dst::clear_all(
      chisqrBcOut, x0BcOut, y0BcOut, u0BcOut, v0BcOut,
      residual_x_cluster, residual_y_cluster,
      residual_x_hit, residual_y_hit
    );
  }

  void clear()
  {
    clearBasicInfo();
    clearRawHits();
    clearClusters();
    clearBcOut();
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
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time
  TTreeReaderValue<std::vector<Double_t>>* cobo_id;   // CoBo ID

  TTreeReaderValue<UInt_t>* evnum_bcout;
  TTreeReaderValue<Int_t>* ntBcOut;
  TTreeReaderValue<std::vector<Double_t>>* chisqrBcOut;
  TTreeReaderValue<std::vector<Double_t>>* x0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* y0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* u0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* v0BcOut;
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
  using namespace root;
  using namespace dst;


  //_____________________________________________________________________________
  // Pick the best BcOut and extract (u0,v0,x0,y0) tracking parameters.
  // NOTE: "chisqr" in this tree is actually reduced-chi2-like value (chi2/ndf),
  // so we select the candidate whose value is closest to 1.0.
  Bool_t SelectBcOutTracking(
    Int_t& best_bcout_idx,
    Double_t& u0, Double_t& v0,
    Double_t& x0, Double_t& y0)
  {
    if(**src.ntBcOut <= 0) return false;

    Double_t best_score = 1.0e+10;
    Double_t best_chi2 = 1.0e+10;
    best_bcout_idx = -1;
    for(Int_t i_bcout = 0; i_bcout < **src.ntBcOut; ++i_bcout){
      const Double_t chi2 = (**src.chisqrBcOut)[i_bcout];
      if(!std::isfinite(chi2) || chi2 <= 0.) continue;
      const Double_t score = TMath::Abs(chi2 - 1.0);
      if(score < best_score || (score == best_score && chi2 < best_chi2)){
        best_score = score;
        best_chi2 = chi2;
        best_bcout_idx = i_bcout;
      }
    }

    event.best_bcout_id = best_bcout_idx;
    event.ntBcOut       = **src.ntBcOut;
    event.chisqrBcOut   = **src.chisqrBcOut;
    event.x0BcOut       = **src.x0BcOut;
    event.y0BcOut       = **src.y0BcOut;
    event.u0BcOut       = **src.u0BcOut;
    event.v0BcOut       = **src.v0BcOut;

    u0 = TMath::QuietNaN();
    v0 = TMath::QuietNaN();
    x0 = TMath::QuietNaN();
    y0 = TMath::QuietNaN();
    if(best_bcout_idx >= 0){
      u0 = (**src.u0BcOut)[best_bcout_idx];
      v0 = (**src.v0BcOut)[best_bcout_idx];
      x0 = (**src.x0BcOut)[best_bcout_idx];
      y0 = (**src.y0BcOut)[best_bcout_idx];
    }

    return true;
  }

  //_____________________________________________________________________________
  void FillTpcHitsAndResiduals(
    TPCAnalyzer& tpc_ana,
    Int_t best_bcout_idx,
    Double_t u0, Double_t v0,
    Double_t x0, Double_t y0,
    TPCEventAnalyzer& event_ana)
  {
    Int_t nh_tpc = 0;
    for(Int_t layer = 0; layer < NumOfLayersTPC; ++layer){
      const auto& hc_hit = tpc_ana.GetTPCHC(layer);
      for(const auto& hit : hc_hit){
        if(!hit || !hit->IsGood()) continue;

        const Int_t row = hit->GetRow();
        const auto& pos = hit->GetPosition();
        ThreeVector local_pos(pos.X(), pos.Y(), pos.Z());
        const ThreeVector global_pos = gGeom.Local2GlobalPos("HypTPC", local_pos);

        // Raw Hit Info
#if RawHit
        event.raw_hitpos_x.push_back(pos.X());
        event.raw_hitpos_y.push_back(pos.Y());
        event.raw_hitpos_z.push_back(pos.Z());
        event.raw_de.push_back(hit->GetCDe());
        event.raw_padid.push_back(hit->GetPad());
        event.raw_layer.push_back(layer);
        event.raw_row.push_back(row);
#endif

        // Residual calculation
        Double_t res_x = TMath::QuietNaN();
        Double_t res_y = TMath::QuietNaN();
        if(best_bcout_idx >= 0){
          const Double_t z_g = global_pos.z();
          const Double_t x_bc = u0 * z_g + x0;
          const Double_t y_bc = v0 * z_g + y0;
          res_x = global_pos.x() - x_bc;
          res_y = global_pos.y() - y_bc;

          const Bool_t in_window = (TMath::Abs(res_x) <= MAX_RESIDUAL &&
                                    TMath::Abs(res_y) <= MAX_RESIDUAL);
          event_ana.FillResidualHist(
            "TPCHit", layer, row, hit->GetPad(),
            res_x, res_y, x_bc, y_bc, in_window,
            hit->GetCTime(),
            TVector3(local_pos.x(), local_pos.y(), local_pos.z()));
        }
        event.residual_x_hit.push_back(res_x);
        event.residual_y_hit.push_back(res_y);
        
        ++nh_tpc;
      } // for(hit)
    } // for(layer)

    event.nhTpc = nh_tpc;
  }

  //_____________________________________________________________________________
  void FillTpcClustersAndResiduals(
    TPCAnalyzer& tpc_ana,
    Int_t best_bcout_idx,
    Double_t u0, Double_t v0,
    Double_t x0, Double_t y0,
    TPCEventAnalyzer& event_ana)
  {
    Int_t ncl_tpc = 0;
    for(Int_t layer = 0; layer < NumOfLayersTPC; ++layer){
      const auto& hc_cl = tpc_ana.GetTPCClCont(layer);
      for(const auto& cl : hc_cl){
        if(!cl || !cl->IsGood()) continue;

        TPCHit* center_hit = cl->GetCenterHit();
        const Int_t center_row = center_hit ? center_hit->GetRow() : -1;
        const ThreeVector local_pos(cl->GetX(), cl->GetY(), cl->GetZ());
        const ThreeVector global_pos = gGeom.Local2GlobalPos("HypTPC", local_pos);

        // Raw Cluster Info
#if RawCluster
        event.cluster_x.push_back(local_pos.X());
        event.cluster_y.push_back(local_pos.Y());
        event.cluster_z.push_back(local_pos.Z());
        event.cluster_de.push_back(cl->GetDe());
        event.cluster_size.push_back(cl->GetClusterSize());
        event.cluster_layer.push_back(layer);
        event.cluster_mrow.push_back(cl->MeanRow());
        event.cluster_row_center.push_back(center_row);

        if (center_hit) {
          const TVector3& center_pos = center_hit->GetPosition();
          event.cluster_de_center.push_back(center_hit->GetCDe());
          event.cluster_x_center.push_back(center_pos.X());
          event.cluster_y_center.push_back(center_pos.Y());
          event.cluster_z_center.push_back(center_pos.Z());
        } else {
          event.cluster_de_center.push_back(TMath::QuietNaN());
          event.cluster_x_center.push_back(TMath::QuietNaN());
          event.cluster_y_center.push_back(TMath::QuietNaN());
          event.cluster_z_center.push_back(TMath::QuietNaN());
        }
#endif

        // Residual calculation
        Double_t res_x = TMath::QuietNaN();
        Double_t res_y = TMath::QuietNaN();
        if(best_bcout_idx >= 0){
          const Double_t z_g = global_pos.z();
          const Double_t x_bc = u0 * z_g + x0;
          const Double_t y_bc = v0 * z_g + y0;
          res_x = global_pos.x() - x_bc;
          res_y = global_pos.y() - y_bc;

          TPCHit* center_hit = cl->GetCenterHit();
          if(center_hit){
            const Bool_t in_window = (TMath::Abs(res_x) <= MAX_RESIDUAL &&
                                      TMath::Abs(res_y) <= MAX_RESIDUAL);
            event_ana.FillResidualHist(
              "TPCCl", layer, center_row, center_hit->GetPad(),
              res_x, res_y, global_pos.x(), y_bc, in_window,
              center_hit->GetCTime(),
              TVector3(local_pos.x(), local_pos.y(), local_pos.z()));
          }
        }
        event.residual_x_cluster.push_back(res_x);
        event.residual_y_cluster.push_back(res_y);

        ++ncl_tpc;
      } // for(cl)
    } // for(layer)

    event.nclTpc = ncl_tpc;
  }
}

//_____________________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!dst::ValidateUserParams(gUser, kUserParamKeys))
    return EXIT_FAILURE;
  if(!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if(!dst::SetupReaders())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
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
  for(const auto& name : TreeName) if(name != "") n_input_files++;

  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
    if(TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=n_input_files || open_tree!=n_input_files){
    spdlog::error("DstOpen Failed: opened files/trees mismatch based on TreeName definitions."
                  " expected: {}, open_file: {}, open_tree: {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if(ievent%100==0){
    std::cout << "#D Event Number: "
              << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  // Check event numbers
  evnumPerFile = { **src.evnum, **src.evnum_bcout };
  if(!dst::CheckEventNumbers(evnumPerFile, ievent, {TreeName[kTpcHit], TreeName[kBcOut]})){
    return false;
  }

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;
  event.cobo_id  = **src.cobo_id;
  HF1("Status", event.status++);

  if(**src.nhTpc == 0)
    return true;
  HF1("Status", event.status++);

  if(!TPCEventAnalyzer::ValidateCoboClocks(event.clkTpc))
    return true;
  HF1("Status", event.status++);

  // Re-calculate TPC hits
  TPCAnalyzer tpc_ana;
  tpc_ana.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  HF1("Status", event.status++);

  static TPCEventAnalyzer event_ana;
  event_ana.SetClock(**src.clkTpc);

  // check and select best BcOut
  Int_t best_bcout_idx = -1;
  Double_t u0 = TMath::QuietNaN();
  Double_t v0 = TMath::QuietNaN();
  Double_t x0 = TMath::QuietNaN();
  Double_t y0 = TMath::QuietNaN();
  if(!SelectBcOutTracking(best_bcout_idx, u0, v0, x0, y0))
    return true;
  HF1("Status", event.status++);

  FillTpcHitsAndResiduals(tpc_ana, best_bcout_idx, u0, v0, x0, y0, event_ana);
  FillTpcClustersAndResiduals(tpc_ana, best_bcout_idx, u0, v0, x0, y0, event_ana);
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
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i])       delete TTreeCont[i];
    if(TFileCont[i])       delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::SetupReaders()
{
  // -------------------------------------------------------
  // TPC Hit
  // -------------------------------------------------------
  if (!dst::SetupReader(kTpcHit, "kTpcHit")) return false;

  dst::SetBranch(TTreeReaderCont[kTpcHit], "run_number",   src.runnum);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "event_number", src.evnum);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "trig_pat",     src.trigpat);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "trig_flag",    src.trigflag);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "beam_flag",    src.beamflag);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "npadTpc",      src.npadTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "nhTpc",        src.nhTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "layerTpc",     src.layerTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "rowTpc",       src.rowTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "padTpc",       src.padTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "deTpc",        src.deTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "tTpc",         src.tTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "clkTpc",       src.clkTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "cobo_id",      src.cobo_id);

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
  hist::BuildTPCHitBcOutTracking(TPCEventAnalyzer::GetDstCalibFlag());

  tree = new TTree("tpc", "tree of DstTPCTracking");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("clkTpc", &event.clkTpc);
  tree->Branch("cobo_id", &event.cobo_id);

  // BcOut Info
  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("best_bcout_id", &event.best_bcout_id);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);

  // Residual Info
  tree->Branch("residual_x_cluster", &event.residual_x_cluster);
  tree->Branch("residual_y_cluster", &event.residual_y_cluster);
  tree->Branch("residual_x_hit", &event.residual_x_hit);
  tree->Branch("residual_y_hit", &event.residual_y_hit);

#if RawHit
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("raw_hitpos_x", &event.raw_hitpos_x);
  tree->Branch("raw_hitpos_y", &event.raw_hitpos_y);
  tree->Branch("raw_hitpos_z", &event.raw_hitpos_z);
  tree->Branch("raw_de", &event.raw_de);
  tree->Branch("raw_padid", &event.raw_padid);
  tree->Branch("raw_layer", &event.raw_layer);
  tree->Branch("raw_row", &event.raw_row);
#endif

#if RawCluster
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("cluster_x", &event.cluster_x);
  tree->Branch("cluster_y", &event.cluster_y);
  tree->Branch("cluster_z", &event.cluster_z);
  tree->Branch("cluster_de", &event.cluster_de);
  tree->Branch("cluster_size", &event.cluster_size);
  tree->Branch("cluster_layer", &event.cluster_layer);
  tree->Branch("cluster_row_center", &event.cluster_row_center);
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_x_center", &event.cluster_x_center);
  tree->Branch("cluster_y_center", &event.cluster_y_center);
  tree->Branch("cluster_z_center", &event.cluster_z_center);
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
