#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <TTree.h>
#include <TTreeReaderValue.h>
#include <TMath.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "RootHelper.hh"

#include <UnpackerManager.hh>

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto& gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();

constexpr Double_t kMinClusterSize = 3.;
constexpr Double_t kAlphaMaxInner = 0.10;      // rad, |track-pad angle| for inner layers.
constexpr Double_t kAlphaMaxOuter = 0.20;      // rad, |track-pad angle| for outer layers.
constexpr Int_t    kPionPidMask = 0x1;         // TPCLocalTrackHelix pion PID bit.
constexpr Double_t kPionMomentumMin = 0.10;    // GeV/c, inclusive.
constexpr Double_t kPionMomentumMax = 1.00;    // GeV/c, inclusive.
constexpr Double_t kDriftBinWidth = 100.;      // mm = 10 cm.
constexpr Int_t    kNumDriftBins = 7;
}

namespace dst
{
enum kArgc { kProcess, kConfFile, kInput, kOutFile, nArgc };
std::vector<TString> ArgName = {"[Process]", "[ConfFile]", "[DstTPCK18HelixTracking]", "[OutFile]"};
std::vector<TString> TreeName = {"", "", "tpc", ""};
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
Bool_t SetupReaders();
}

struct Event
{
  Int_t status = 0;
  UInt_t runnum = 0;
  UInt_t evnum = 0;
  Int_t n_pad_hits = 0;
  void clear() { status = 0; runnum = 0; evnum = 0; n_pad_hits = 0; }
};

struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum;
  TTreeReaderValue<std::vector<Double_t>>* mom0;
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* theta_diff;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* diffusion_dist;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* diffusion_ratio;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* diffusion_drift;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* diffusion_hit_index;
};

namespace root { Event event; Src src; TTree* tree; }

Int_t main(Int_t argc, char** argv)
{
  std::vector<std::string> arg(argv, argv+argc);
  if (!CheckArg(arg) || !dst::DstOpen(arg) || !gConf.Initialize(arg[dst::kConfFile]) ||
      !gConf.InitializeHistograms() || !gConf.InitializeUnpacker() || !dst::SetupReaders())
    return EXIT_FAILURE;

  Int_t skip = std::max(0, gUnpacker.get_skip());
  Int_t nevent = GetEntries(dst::TTreeCont);
  const Int_t max_loop = gUnpacker.get_max_loop();
  if (max_loop > 0) nevent = skip + max_loop;
  CatchSignal::Set();
  Int_t ievent = skip;
  for (; ievent<nevent && !CatchSignal::Stop(); ++ievent) {
    gCounter.check();
    dst::InitializeEvent();
    if (dst::DstRead(ievent)) root::tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  dst::DstClose();
  return EXIT_SUCCESS;
}

Bool_t dst::InitializeEvent() { event.clear(); return true; }

Bool_t dst::DstOpen(std::vector<std::string> arg)
{
  if (!OpenFile(TFileCont[kInput], arg[kInput]) || !OpenTree(TFileCont[kInput], TTreeCont[kInput], TreeName[kInput]))
    return false;
  const std::vector<TTree*> aligned_trees = {TTreeCont[kInput]};
  if (!CheckEntries(aligned_trees)) return false;
  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  return true;
}

Bool_t dst::SetupReaders()
{
  if (!SetupReader(kInput, "kInput")) return false;
  const auto& reader = TTreeReaderCont[kInput];
  SetBranch(reader, "run_number", src.runnum);
  SetBranch(reader, "event_number", src.evnum);
  // DstTPCTracking does not provide momentum/PID branches.  Keep these
  // inputs optional so diffusion can still be analyzed from that DST.
  if (TTreeCont[kInput]->GetBranch("mom0")) SetBranch(reader, "mom0", src.mom0);
  if (TTreeCont[kInput]->GetBranch("pid"))  SetBranch(reader, "pid", src.pid);
  SetBranch(reader, "theta_diff", src.theta_diff);
  SetBranch(reader, "hitlayer", src.hitlayer);
  SetBranch(reader, "track_cluster_size", src.cluster_size);
  SetBranch(reader, "transverse_diffusion_dist", src.diffusion_dist);
  SetBranch(reader, "transverse_diffusion_ratio", src.diffusion_ratio);
  SetBranch(reader, "transverse_diffusion_drift", src.diffusion_drift);
  SetBranch(reader, "transverse_diffusion_track_hit_index", src.diffusion_hit_index);
  return true;
}

Bool_t dst::DstRead(Int_t ievent)
{
  GetEntry(ievent);
  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  HF1("Status", event.status++);

  const std::vector<Double_t> empty_mom;
  const std::vector<Int_t> empty_pid;
  const auto& mom = src.mom0 ? **src.mom0 : empty_mom;
  const auto& pid = src.pid ? **src.pid : empty_pid;
  const auto& alpha = **src.theta_diff;
  const auto& hitlayer = **src.hitlayer;
  const auto& cluster_size = **src.cluster_size;
  const auto& distance = **src.diffusion_dist;
  const auto& ratio = **src.diffusion_ratio;
  const auto& drift = **src.diffusion_drift;
  const auto& hit_index = **src.diffusion_hit_index;
  const std::size_t ntrack_required = std::min({alpha.size(), hitlayer.size(), cluster_size.size(),
                                                distance.size(), ratio.size(), drift.size(), hit_index.size()});
  const std::size_t ntrack = (src.mom0 && src.pid)
    ? std::min(ntrack_required, std::min(mom.size(), pid.size()))
    : ntrack_required;
  for (std::size_t it=0; it<ntrack; ++it) {
    const Bool_t have_pid_mom = src.pid && src.mom0 && it < pid.size() && it < mom.size();
    const Bool_t pass_pid = !have_pid_mom || ((pid[it] & kPionPidMask) &&
      mom[it] >= kPionMomentumMin && mom[it] <= kPionMomentumMax);
    const std::size_t npad = std::min({distance[it].size(), ratio[it].size(), drift[it].size(), hit_index[it].size()});
    for (std::size_t ipad=0; ipad<npad; ++ipad) {
      const Int_t ihit = static_cast<Int_t>(std::lround(hit_index[it][ipad]));
      if (ihit < 0 || ihit >= static_cast<Int_t>(alpha[it].size()) ||
          ihit >= static_cast<Int_t>(cluster_size[it].size()) ||
          ihit >= static_cast<Int_t>(hitlayer[it].size())) continue;
      const Double_t r = distance[it][ipad];
      const Double_t qfrac = ratio[it][ipad];
      const Double_t drift_mm = drift[it][ipad];
      const Double_t abs_alpha = std::abs(alpha[it][ihit]);
      const Double_t ncluster = cluster_size[it][ihit];
      if (!std::isfinite(r) || !std::isfinite(qfrac) || !std::isfinite(drift_mm) ||
          !std::isfinite(abs_alpha) || !std::isfinite(ncluster) || ncluster < kMinClusterSize) continue;
      const Int_t layer = static_cast<Int_t>(hitlayer[it][ihit]);
      const Double_t alpha_limit = layer <= 9 ? kAlphaMaxInner : kAlphaMaxOuter;
      const Bool_t pass_alpha = abs_alpha < alpha_limit;


      if(ncluster == 1)HF2("TPCDiffusion_Cluster1", r, qfrac);
      else if(ncluster == 2)HF2("TPCDiffusion_Cluster2", r, qfrac);
      else if(ncluster == 3)HF2("TPCDiffusion_Cluster3", r, qfrac);
      else if(ncluster == 4)HF2("TPCDiffusion_Cluster4", r, qfrac);
      
      HF1("TPCDiffusion_Alpha", abs_alpha*TMath::RadToDeg());
      HF1("TPCDiffusion_ClusterSize", ncluster);
      HF2("TPCDiffusion_All", r, qfrac);
      if (pass_alpha) HF2("TPCDiffusion_AlphaCut", r, qfrac);
      if (pass_pid) HF2("TPCDiffusion_PidCut", r, qfrac);
      if (pass_alpha && pass_pid) {
        HF2("TPCDiffusion_AlphaPidCut", r, qfrac);
        HF1("TPCDiffusion_Drift", drift_mm/10.);
        const Int_t drift_bin = static_cast<Int_t>(drift_mm/kDriftBinWidth);
        if (drift_bin >= 0 && drift_bin < kNumDriftBins)
          HF2(Form("TPCDiffusion_Drift%02d_%02dcm", 10*drift_bin, 10*(drift_bin+1)), r, qfrac);
      }
      ++event.n_pad_hits;
    }
  }
  HF1("Status", event.status++);
  HF1("TPCDiffusion_NPadHits", event.n_pad_hits);
  return true;
}

Bool_t dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  TFileCont[kOutFile]->Close();
  for (auto* reader : TTreeReaderCont) delete reader;
  for (auto* input : TFileCont) if (input) delete input;
  return true;
}

Bool_t ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTPCDiffusion();
  root::tree = new TTree("tpc", "tree of DstTPCDiffusion");
  root::tree->Branch("status", &root::event.status);
  root::tree->Branch("run_number", &root::event.runnum);
  root::tree->Branch("event_number", &root::event.evnum);
  root::tree->Branch("n_diffusion_pad_hits", &root::event.n_pad_hits);
  return true;
}

Bool_t ConfMan::InitializeParameterFiles() { return true; }
Bool_t ConfMan::FinalizeProcess() { return true; }
