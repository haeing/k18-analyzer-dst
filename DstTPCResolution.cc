#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include "ConfMan.hh"
#include "TPCHit.hh"
#include "TPCLTrackHit.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

namespace dst {
enum kArgc { kProcess, kConfFile, kInput, kOutFile, nArgc };
std::vector<TString> ArgName = {"[Process]", "[ConfFile]", "[DstTPCK18HelixTracking]", "[OutFile]"};
std::vector<TString> TreeName = {"", "", "tpc", ""};
TFile* input_file = nullptr; TFile* output_file = nullptr;
TTree* input_tree = nullptr; TTree* output_tree = nullptr;
}
namespace dst {
Bool_t DstOpen(std::vector<std::string>);
Bool_t SetupReaders();
Bool_t DstRead(Int_t);
Bool_t DstClose();
}
namespace {
using V = std::vector<Double_t>; using VV = std::vector<V>;
struct Src {
  UInt_t run=0, event=0; Int_t ntrack=0;
  VV *x=nullptr,*y=nullptr,*z=nullptr,*layer=nullptr,*mrow=nullptr;
  V *cx=nullptr,*cy=nullptr,*z0=nullptr,*r=nullptr,*dz=nullptr;
} src;
VV exclusive_x, exclusive_y, exclusive_z, residual_x, residual_y, residual_z, residual_chi, exclusive_layer;
}

Bool_t dst::DstOpen(std::vector<std::string> arg) {
  input_file = TFile::Open(arg[kInput].c_str(), "READ");
  if (!input_file || input_file->IsZombie()) return false;
  input_tree = dynamic_cast<TTree*>(input_file->Get("tpc"));
  if (!input_tree) return false;
  output_file = TFile::Open(arg[kOutFile].c_str(), "RECREATE");
  if (!output_file || output_file->IsZombie()) return false;
  input_tree->SetBranchAddress("run_number", &src.run);
  input_tree->SetBranchAddress("event_number", &src.event);
  input_tree->SetBranchAddress("ntTpc", &src.ntrack);
  input_tree->SetBranchAddress("hitpos_x", &src.x);
  input_tree->SetBranchAddress("hitpos_y", &src.y);
  input_tree->SetBranchAddress("hitpos_z", &src.z);
  input_tree->SetBranchAddress("hitlayer", &src.layer);
  input_tree->SetBranchAddress("track_cluster_mrow", &src.mrow);
  input_tree->SetBranchAddress("helix_cx", &src.cx);
  input_tree->SetBranchAddress("helix_cy", &src.cy);
  input_tree->SetBranchAddress("helix_z0", &src.z0);
  input_tree->SetBranchAddress("helix_r", &src.r);
  input_tree->SetBranchAddress("helix_dz", &src.dz);
  output_tree = input_tree->CloneTree(0);
  output_tree->Branch("exclusive_track_x", &exclusive_x);
  output_tree->Branch("exclusive_track_y", &exclusive_y);
  output_tree->Branch("exclusive_track_z", &exclusive_z);
  output_tree->Branch("exclusive_residual_x", &residual_x);
  output_tree->Branch("exclusive_residual_y", &residual_y);
  output_tree->Branch("exclusive_residual_z", &residual_z);
  output_tree->Branch("exclusive_chisqr", &residual_chi);
  output_tree->Branch("exclusive_layer", &exclusive_layer);
  return true;
}
Bool_t dst::SetupReaders() { return input_tree != nullptr; }

static Bool_t RefitOne(const VV& x, const VV& y, const VV& z, const VV& layers, const VV& mrows,
                       const V& cx, const V& cy, const V& z0, const V& rr, const V& dz,
                       std::size_t it, std::size_t omit, TVector3& pred, Double_t& chi) {
  const std::size_t n=std::min({x[it].size(),y[it].size(),z[it].size(),layers[it].size(),mrows[it].size()});
  if (n < 7 || omit >= n) return false;
  TPCLocalTrackHelix fit;
  if (it<cx.size()&&it<cy.size()&&it<z0.size()&&it<rr.size()&&it<dz.size()) {
    Double_t par[5]={cx[it],cy[it],z0[it],rr[it],dz[it]};
    if (std::all_of(par,par+5,[](Double_t q){return std::isfinite(q);})) fit.SetParam(par);
  }
  // TPCHit owns registered TPCLTrackHit objects; keep only TPCHit owners here.
  std::vector<std::unique_ptr<TPCHit>> hs;
  for(std::size_t ih=0; ih<n; ++ih) if(ih!=omit) {
    auto h=std::make_unique<TPCHit>(static_cast<Int_t>(std::lround(layers[it][ih])),mrows[it][ih]);
    h->AddHit(0.,0.); h->SetPosition(TVector3(x[it][ih],y[it][ih],z[it][ih])); h->SetIsGood(true);
    auto th=std::make_unique<TPCLTrackHit>(h.get()); fit.AddTPCHit(th.release()); hs.push_back(std::move(h));
  }
  if(!fit.DoFit()) return false;
  TVector3 measured(x[it][omit],y[it][omit],z[it][omit]);
  TVector3 res=fit.CalcResidual(measured); if(!std::isfinite(res.X())||!std::isfinite(res.Y())||!std::isfinite(res.Z())) return false;
  pred=measured-res; chi=fit.GetChiSquare(); return true;
}

Bool_t dst::DstRead(Int_t ievent) {
  input_tree->GetEntry(ievent);
  exclusive_x.clear(); exclusive_y.clear(); exclusive_z.clear(); residual_x.clear(); residual_y.clear(); residual_z.clear(); residual_chi.clear(); exclusive_layer.clear();
  const std::size_t nt=std::min<std::size_t>(std::max(0,src.ntrack),std::min({src.x?src.x->size():0,src.y?src.y->size():0,src.z?src.z->size():0}));
  exclusive_x.resize(nt); exclusive_y.resize(nt); exclusive_z.resize(nt); residual_x.resize(nt); residual_y.resize(nt); residual_z.resize(nt); residual_chi.resize(nt); exclusive_layer.resize(nt);
  for(std::size_t it=0;it<nt;++it) {
    const std::size_t n=(*src.x)[it].size(); exclusive_x[it].assign(n,NAN); exclusive_y[it].assign(n,NAN); exclusive_z[it].assign(n,NAN); residual_x[it].assign(n,NAN); residual_y[it].assign(n,NAN); residual_z[it].assign(n,NAN); residual_chi[it].assign(n,NAN); exclusive_layer[it]=(*src.layer)[it];
    for(std::size_t ih=0;ih<n;++ih) { TVector3 p; Double_t chi=NAN; if(!RefitOne(*src.x,*src.y,*src.z,*src.layer,*src.mrow,*src.cx,*src.cy,*src.z0,*src.r,*src.dz,it,ih,p,chi)) continue; TVector3 m((*src.x)[it][ih],(*src.y)[it][ih],(*src.z)[it][ih]); TVector3 q=p-m; exclusive_x[it][ih]=p.X(); exclusive_y[it][ih]=p.Y(); exclusive_z[it][ih]=p.Z(); residual_x[it][ih]=q.X(); residual_y[it][ih]=q.Y(); residual_z[it][ih]=q.Z(); residual_chi[it][ih]=chi; }
  }
  output_tree->Fill(); return true;
}
Bool_t dst::DstClose() { output_file->cd(); output_tree->Write(); output_file->Close(); input_file->Close(); delete output_file; delete input_file; return true; }
Bool_t ConfMan::InitializeParameterFiles() {
  return InitializeParameter<TPCParamMan>("TPCPRM") &&
         InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
         InitializeParameter<UserParamMan>("USER");
}
Bool_t ConfMan::FinalizeProcess() { return true; }

int main(int argc,char** argv) {
  if(argc!=4) { std::cerr<<"Usage: DstTPCResolution [ConfFile] [DstTPCK18HelixTracking] [OutFile]\n"; return 1; }
  if(!ConfMan::GetInstance().Initialize(argv[1])) return 1;
  std::vector<std::string> a(argv,argv+argc); if(!dst::DstOpen(a)||!dst::SetupReaders()) return 1;
  for(Long64_t i=0;i<dst::input_tree->GetEntries();++i) if(!dst::DstRead(i)) break;
  dst::DstClose(); return 0;
}
