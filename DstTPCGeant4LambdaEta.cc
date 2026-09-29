// -*- C++ -*-

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "Kinematics.hh"
#include "RootHelper.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>
#include <UnpackerManager.hh>

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;

  const auto& gUnpacker = GUnpacker::get_instance();
  auto&       gConf = ConfMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();

  constexpr Double_t kMassPi = 0.13957039;
  constexpr Double_t kMassP = 0.93827208;
  constexpr Double_t kMassLambda = 1.115683;
  constexpr Double_t kTightCloseDist = 5.;
  constexpr Double_t kTightChiSqr = 2.;
  constexpr Double_t kPiMinusMomCut = .3;
  constexpr Int_t kSigmaPion = 0x1;
  constexpr Int_t kSigmaProton = 0x4;
  constexpr Double_t kMeVToGeV = 1.e-3;

  // Beam definitions used for the source-comparison histograms.
  // The reaction-generator value is BEAM[0] (MeV/c); TPCReco is the beam-only
  // helix-fit result (GeV/c).  Do not use TPC hit truth for a target vertex.
  enum class BeamMomentumSource : Int_t {
    kReactionGeneratorTruth = 0,
    kTPCHelixReconstruction = 1
  };
  constexpr Int_t kNBeamMomentumSources = 2;

  const char* BeamMomentumSourceName(BeamMomentumSource source)
  {
    switch (source) {
    case BeamMomentumSource::kReactionGeneratorTruth: return "BeamTruth";
    case BeamMomentumSource::kTPCHelixReconstruction: return "BeamRecon";
    }
    return "Unknown";
  }
  Bool_t CalcHelixDx(Double_t radius, Double_t theta_min, Double_t theta_max, Double_t& dx)
  {
    if (!std::isfinite(radius) || !std::isfinite(theta_min) || !std::isfinite(theta_max)) return false;
    dx = radius * (std::cos(theta_max) - std::cos(theta_min));
    return std::isfinite(dx);
  }

  Bool_t IsProton(Int_t q, Double_t p, Double_t dedx)
  {
    return (Kinematics::HypTPCdEdxPID(dedx, q*p) & kSigmaProton) &&
      q > 0 && dedx >= 25.;
  }

  Bool_t IsPion(Int_t q, Double_t p, Double_t dedx)
  {
    return (Kinematics::HypTPCdEdxPID(dedx, q*p) & kSigmaPion) && q < 0;
  }

  Bool_t IsTightCandidate(Double_t close_dist, Double_t chisqr_p,
                          Double_t chisqr_pi, Double_t mom_pi)
  {
    return std::isfinite(close_dist) && close_dist < kTightCloseDist &&
      chisqr_p <= kTightChiSqr && chisqr_pi <= kTightChiSqr &&
      mom_pi < kPiMinusMomCut;
  }

  Bool_t BeamMomentumAtZ(const std::vector<Double_t>& z,
                         const std::vector<Double_t>& px,
                         const std::vector<Double_t>& py,
                         const std::vector<Double_t>& pz,
                         Double_t vertex_z,
                         Double_t& bx, Double_t& by, Double_t& bz,
                         Double_t& selected_z, Int_t& selected_id)
  {
    selected_z = TMath::QuietNaN();
    selected_id = -1;
    if (z.size() != px.size() || z.size() != py.size() || z.size() != pz.size())
      return false;

    Double_t best_dz = std::numeric_limits<Double_t>::max();
    for (std::size_t i=0; i<z.size(); ++i) {
      if (!std::isfinite(z[i]) || !std::isfinite(px[i]) ||
          !std::isfinite(py[i]) || !std::isfinite(pz[i]))
        continue;
      const Double_t dz = TMath::Abs(z[i] - vertex_z);
      if (dz >= best_dz)
        continue;
      best_dz = dz;
      bx = px[i]; by = py[i]; bz = pz[i];
      selected_z = z[i];
      selected_id = static_cast<Int_t>(i);
    }
    return selected_id >= 0;
  }

  Bool_t SelectBeamMomentum(BeamMomentumSource source, Double_t vertex_z,
                            Double_t initial_px, Double_t initial_py, Double_t initial_pz,
                            const std::vector<Double_t>& reco_z,
                            const std::vector<Double_t>& reco_px,
                            const std::vector<Double_t>& reco_py,
                            const std::vector<Double_t>& reco_pz,
                            Double_t& bx, Double_t& by, Double_t& bz,
                            Double_t& selected_z, Int_t& selected_id)
  {
    if (source == BeamMomentumSource::kReactionGeneratorTruth) {
      if (!std::isfinite(initial_px) || !std::isfinite(initial_py) || !std::isfinite(initial_pz))
        return false;
      bx = initial_px*kMeVToGeV;
      by = initial_py*kMeVToGeV;
      bz = initial_pz*kMeVToGeV;
      selected_z = TMath::QuietNaN();
      selected_id = -1;
      return true;
    }
    return BeamMomentumAtZ(reco_z, reco_px, reco_py, reco_pz, vertex_z,
                           bx, by, bz, selected_z, selected_id);
  }

  enum class TruthTargetRegion { kOther, kLH2Inside, kLH2Outside };
  TruthTargetRegion ClassifyTruthTargetVertex(Double_t x, Double_t y, Double_t z)
  {
    constexpr Double_t kTargetZ = -143., kInnerR = 38., kOuterR = 42., kHalfY = 48.;
    if (TMath::Abs(y) >= kHalfY) return TruthTargetRegion::kOther;
    const Double_t r = TMath::Hypot(x, z-kTargetZ);
    if (r < kInnerR) return TruthTargetRegion::kLH2Inside;
    if (r > kOuterR) return TruthTargetRegion::kLH2Outside;
    return TruthTargetRegion::kOther;
  }

  Bool_t FindTruthLambdaDecay(const std::vector<Int_t>& vtx_type,
                              const std::vector<Int_t>& vtx_motherpid,
                              const std::vector<std::vector<Double_t>>& vtx_pid,
                              const std::vector<std::vector<Double_t>>& vtx_px,
                              const std::vector<std::vector<Double_t>>& vtx_py,
                              const std::vector<std::vector<Double_t>>& vtx_pz,
                              TVector3& proton_mom, TVector3& pion_mom, Int_t& vtx_index)
  {
    for (std::size_t iv=0; iv<vtx_type.size(); ++iv) {
      if (vtx_type[iv] != 1 || iv >= vtx_motherpid.size() ||
          vtx_motherpid[iv] != 3122 || iv >= vtx_pid.size() || iv >= vtx_px.size() ||
          iv >= vtx_py.size() || iv >= vtx_pz.size()) continue;
      const auto n = std::min({vtx_pid[iv].size(), vtx_px[iv].size(),
                               vtx_py[iv].size(), vtx_pz[iv].size()});
      Int_t proton = -1, pion = -1;
      for (std::size_t ip=0; ip<n; ++ip) {
        if (vtx_pid[iv][ip] == 2212) proton = static_cast<Int_t>(ip);
        if (vtx_pid[iv][ip] == -211) pion = static_cast<Int_t>(ip);
      }
      if (proton < 0 || pion < 0) continue;
      proton_mom.SetXYZ(vtx_px[iv][proton], vtx_py[iv][proton], vtx_pz[iv][proton]);
      pion_mom.SetXYZ(vtx_px[iv][pion], vtx_py[iv][pion], vtx_pz[iv][pion]);
      vtx_index = static_cast<Int_t>(iv);
      return std::isfinite(proton_mom.X()) && std::isfinite(pion_mom.X());
    }
    return false;
  }

  Bool_t FindReactionVertexZ(const std::vector<Int_t>& vtx_type,
                             const std::vector<Double_t>& vtx_z, Double_t& z)
  {
    for (std::size_t iv=0; iv<vtx_type.size() && iv<vtx_z.size(); ++iv)
      if (vtx_type[iv] == 0 && std::isfinite(vtx_z[iv])) { z = vtx_z[iv]; return true; }
    return false;
  }
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kGeant4Tpc, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[DstTPCGeant4HelixTracking]", "[OutFile]" };
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
  Int_t trigflag;
  Int_t beam_reco_found;
  std::vector<Double_t> lambda_mass, lambda_close_dist;
  std::vector<Double_t> lambda_vtx_x, lambda_vtx_y, lambda_vtx_z;
  std::vector<Double_t> lambda_mom_x, lambda_mom_y, lambda_mom_z;
  std::vector<Double_t> X_mass, X_mass2;
  std::vector<Double_t> X_vtx_x, X_vtx_y, X_vtx_z;
  std::vector<Double_t> X_mom_x, X_mom_y, X_mom_z;
  std::vector<Double_t> X_beam_px, X_beam_py, X_beam_pz, X_beam_vp_z;
  std::vector<Int_t> X_beam_vp_id, X_beam_source;

  void clear()
  {
    status = 0; runnum = 0; evnum = 0; trigflag = 0; beam_reco_found = 0;
    dst::clear_all(lambda_mass, lambda_close_dist,
                   lambda_vtx_x, lambda_vtx_y, lambda_vtx_z,
                   lambda_mom_x, lambda_mom_y, lambda_mom_z,
                   X_mass, X_mass2, X_vtx_x, X_vtx_y, X_vtx_z,
                   X_mom_x, X_mom_y, X_mom_z,
                   X_beam_px, X_beam_py, X_beam_pz, X_beam_vp_z,
                   X_beam_vp_id, X_beam_source);
  }
};

//_____________________________________________________________________________
struct Src
{
    TTreeReaderValue<Int_t>* event_number;
    TTreeReaderValue<Int_t>* trig_flag;
    TTreeReaderValue<Int_t>* nt_tpc;
    TTreeReaderValue<Int_t>* beam_found;
    TTreeReaderValue<std::vector<Int_t>>* is_beam;
    TTreeReaderValue<std::vector<Int_t>>* is_accidental;
    TTreeReaderValue<std::vector<Int_t>>* charge;
    TTreeReaderValue<std::vector<Double_t>>* chisqr;
    TTreeReaderValue<std::vector<Double_t>>* dedx;
    TTreeReaderValue<std::vector<Double_t>>* mom0;
    TTreeReaderValue<std::vector<Double_t>>* dE;
    TTreeReaderValue<std::vector<Double_t>>* mom0_x;
    TTreeReaderValue<std::vector<Double_t>>* mom0_y;
    TTreeReaderValue<std::vector<Double_t>>* mom0_z;
    TTreeReaderValue<std::vector<Int_t>>* pid;
    TTreeReaderValue<std::vector<Double_t>>* helix_r;
    TTreeReaderValue<std::vector<Double_t>>* helix_dz;
    TTreeReaderValue<std::vector<Double_t>>* helix_theta_min;
    TTreeReaderValue<std::vector<Double_t>>* helix_theta_max;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* close;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtx_x;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtx_y;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtx_z;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_x;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_y;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_z;
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
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;
    TTreeReaderValue<Double_t>* beam_initial_px_g4;
    TTreeReaderValue<Double_t>* beam_initial_py_g4;
    TTreeReaderValue<Double_t>* beam_initial_pz_g4;
    TTreeReaderValue<Double_t>* reaction_beam_px_g4;
    TTreeReaderValue<Double_t>* reaction_beam_py_g4;
    TTreeReaderValue<Double_t>* reaction_beam_pz_g4;
    TTreeReaderValue<std::vector<Double_t>>* beam_reco_hit_z;
    TTreeReaderValue<std::vector<Double_t>>* beam_reco_px;
    TTreeReaderValue<std::vector<Double_t>>* beam_reco_py;
    TTreeReaderValue<std::vector<Double_t>>* beam_reco_pz;
    TTreeReaderValue<std::vector<Int_t>>* reaction_vtx_type;
    TTreeReaderValue<std::vector<Int_t>>* reaction_vtx_motherpid;
    TTreeReaderValue<std::vector<Double_t>>* reaction_vtx_x;
    TTreeReaderValue<std::vector<Double_t>>* reaction_vtx_y;
    TTreeReaderValue<std::vector<Double_t>>* reaction_vtx_z;
    TTreeReaderValue<std::vector<Int_t>>* truth_found;
    TTreeReaderValue<std::vector<Int_t>>* truth_source;
    TTreeReaderValue<std::vector<Int_t>>* truth_pid;
    TTreeReaderValue<std::vector<Int_t>>* truth_vtx_index;
    TTreeReaderValue<std::vector<Int_t>>* truth_vtx_type;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* reaction_vtx_trackpid;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* reaction_vtx_px;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* reaction_vtx_py;
    TTreeReaderValue<std::vector<std::vector<Double_t>>>* reaction_vtx_pz;
};

namespace root
{
  Event event;
  Src src;
  TTree *tree;
}

using namespace root;

namespace
{
  void FillCandidate(Double_t close_dist, Double_t x, Double_t y, Double_t z,
                     const TVector3& proton_mom, const TVector3& pion_mom,
                     Double_t bx, Double_t by, Double_t bz,
                     Double_t beam_z, Int_t beam_id, Int_t beam_source,
                     const char* source_name)
  {
    TLorentzVector proton, pion, beam, target, lambda_for_missing;
    proton.SetXYZM(proton_mom.X(), proton_mom.Y(), proton_mom.Z(), kMassP);
    pion.SetXYZM(pion_mom.X(), pion_mom.Y(), pion_mom.Z(), kMassPi);
    beam.SetXYZM(bx, by, bz, pdg::KaonMass());
    target.SetXYZM(0., 0., 0., pdg::ProtonMass());
    const TLorentzVector lambda = proton + pion;
    lambda_for_missing.SetXYZM(lambda.Px(), lambda.Py(), lambda.Pz(), kMassLambda);
    const TLorentzVector missing = beam + target - lambda_for_missing;

    event.lambda_mass.push_back(lambda.M());
    event.lambda_close_dist.push_back(close_dist);
    event.lambda_vtx_x.push_back(x); event.lambda_vtx_y.push_back(y); event.lambda_vtx_z.push_back(z);
    event.lambda_mom_x.push_back(lambda.Px()); event.lambda_mom_y.push_back(lambda.Py()); event.lambda_mom_z.push_back(lambda.Pz());
    event.X_mass.push_back(missing.M()); event.X_mass2.push_back(missing.M2());
    event.X_vtx_x.push_back(x); event.X_vtx_y.push_back(y); event.X_vtx_z.push_back(z);
    event.X_mom_x.push_back(missing.Px()); event.X_mom_y.push_back(missing.Py()); event.X_mom_z.push_back(missing.Pz());
    event.X_beam_px.push_back(bx); event.X_beam_py.push_back(by); event.X_beam_pz.push_back(bz);
    event.X_beam_vp_z.push_back(beam_z); event.X_beam_vp_id.push_back(beam_id);
    event.X_beam_source.push_back(beam_source);

    const TString prefix = Form("Geant4LambdaEta_%s_", source_name);
    HF1(prefix+"ReactionReconProtonMomX", proton.Px());
    HF1(prefix+"ReactionReconProtonMomY", proton.Py());
    HF1(prefix+"ReactionReconProtonMomZ", proton.Pz());
    HF1(prefix+"ReactionReconProtonMom", proton.P());
    HF1(prefix+"ReactionReconPionMomX", pion.Px());
    HF1(prefix+"ReactionReconPionMomY", pion.Py());
    HF1(prefix+"ReactionReconPionMomZ", pion.Pz());
    HF1(prefix+"ReactionReconPionMom", pion.P());
    HF1(prefix+"ReactionReconLambdaMass", lambda.M());
    HF1(prefix+"ReactionReconLambdaMomX", lambda.Px());
    HF1(prefix+"ReactionReconLambdaMomY", lambda.Py());
    HF1(prefix+"ReactionReconLambdaMomZ", lambda.Pz());
    HF1(prefix+"ReactionReconLambdaMom", lambda.P());
    HF1(prefix+"ReactionReconMissingMass", missing.M());
    HF1(prefix+"ReactionReconMissingMass2", missing.M2());
  }

  void FillTruthCandidate(const TVector3& proton_mom, const TVector3& pion_mom,
                          Double_t bx, Double_t by, Double_t bz,
                          const char* source_name)
  {
    TLorentzVector proton, pion, beam, target, lambda_for_missing;
    proton.SetXYZM(proton_mom.X(), proton_mom.Y(), proton_mom.Z(), kMassP);
    pion.SetXYZM(pion_mom.X(), pion_mom.Y(), pion_mom.Z(), kMassPi);
    beam.SetXYZM(bx, by, bz, pdg::KaonMass());
    target.SetXYZM(0., 0., 0., pdg::ProtonMass());
    const TLorentzVector lambda = proton + pion;
    lambda_for_missing.SetXYZM(lambda.Px(), lambda.Py(), lambda.Pz(), kMassLambda);
    const TLorentzVector missing = beam + target - lambda_for_missing;

    const TString prefix = Form("Geant4LambdaEta_%s_", source_name);
    HF1(prefix+"ReactionTruthProtonMomX", proton.Px());
    HF1(prefix+"ReactionTruthProtonMomY", proton.Py());
    HF1(prefix+"ReactionTruthProtonMomZ", proton.Pz());
    HF1(prefix+"ReactionTruthProtonMom", proton.P());
    HF1(prefix+"ReactionTruthPionMomX", pion.Px());
    HF1(prefix+"ReactionTruthPionMomY", pion.Py());
    HF1(prefix+"ReactionTruthPionMomZ", pion.Pz());
    HF1(prefix+"ReactionTruthPionMom", pion.P());
    HF1(prefix+"ReactionTruthLambdaMass", lambda.M());
    HF1(prefix+"ReactionTruthLambdaMomX", lambda.Px());
    HF1(prefix+"ReactionTruthLambdaMomY", lambda.Py());
    HF1(prefix+"ReactionTruthLambdaMomZ", lambda.Pz());
    HF1(prefix+"ReactionTruthLambdaMom", lambda.P());
    HF1(prefix+"ReactionTruthMissingMass", missing.M());
    HF1(prefix+"ReactionTruthMissingMass2", missing.M2());
  }
}

//_____________________________________________________________________________
Int_t
main(Int_t argc, char** argv)
{
  std::vector<std::string> arg(argv, argv+argc);
  if (!CheckArg(arg)) return EXIT_FAILURE;
  if (!DstOpen(arg)) return EXIT_FAILURE;
  if (!gConf.Initialize(arg[kConfFile])) return EXIT_FAILURE;
  if (!gConf.InitializeHistograms()) return EXIT_FAILURE;
  if (!gConf.InitializeUnpacker()) return EXIT_FAILURE;
  if (!dst::SetupReaders()) return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t nevent = GetEntries(TTreeCont);
  const Int_t max_loop = gUnpacker.get_max_loop();
  // Never read past the end of the input tree.  TTreeReader keeps its last
  // valid value after an out-of-range SetEntry(), which would otherwise
  // duplicate the final event when MaxLoop exceeds the available entries.
  if (max_loop > 0) nevent = std::min(nevent, skip + max_loop);

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
  if (!OpenFile(TFileCont[kGeant4Tpc], arg[kGeant4Tpc])) return false;
  if (!OpenTree(TFileCont[kGeant4Tpc], TTreeCont[kGeant4Tpc], TreeName[kGeant4Tpc])) return false;
  if (!CheckEntries(TTreeCont)) return false;
  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  return TFileCont[kOutFile] && !TFileCont[kOutFile]->IsZombie();
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if (ievent%1 == 0) std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  GetEntry(ievent);

  event.evnum = **src.event_number;
  event.trigflag = **src.trig_flag;
  event.beam_reco_found = **src.beam_found;
  HF1("Status", event.status++);

  // The Geant4 beam truth is stored hit-by-hit for the selected 7201 track.
  // Select the momentum at the point nearest to the Lambda vertex z.
  const auto& zbeam_reco = **src.beam_reco_hit_z;
  const auto& pxbeam_reco = **src.beam_reco_px;
  const auto& pybeam_reco = **src.beam_reco_py;
  const auto& pzbeam_reco = **src.beam_reco_pz;
  // 7202 BEAM[0]: beam momentum actually supplied to the reaction generator.
  const Double_t initial_px = **src.reaction_beam_px_g4;
  const Double_t initial_py = **src.reaction_beam_py_g4;
  const Double_t initial_pz = **src.reaction_beam_pz_g4;
  const Double_t initial_mom = TMath::Sqrt(initial_px*initial_px +
                                           initial_py*initial_py +
                                           initial_pz*initial_pz);
  const auto& reaction_vtx_type = **src.reaction_vtx_type;
  const auto& reaction_vtx_motherpid = **src.reaction_vtx_motherpid;
  const auto& reaction_vtx_x = **src.reaction_vtx_x;
  const auto& reaction_vtx_y = **src.reaction_vtx_y;
  const auto& reaction_vtx_z = **src.reaction_vtx_z;
  const auto& reaction_vtx_trackpid = **src.reaction_vtx_trackpid;
  const auto& reaction_vtx_px = **src.reaction_vtx_px;
  const auto& reaction_vtx_py = **src.reaction_vtx_py;
  const auto& reaction_vtx_pz = **src.reaction_vtx_pz;
  const auto& truth_found = **src.truth_found;
  const auto& truth_source = **src.truth_source;
  const auto& truth_pid = **src.truth_pid;
  const auto& truth_vtx_index = **src.truth_vtx_index;
  const auto& truth_vtx_type = **src.truth_vtx_type;

  // Fill truth spectra from the original Geant4 decay-vertex record, not from
  // a TPC helix matched to truth hits.
  TVector3 truth_proton_mom, truth_pion_mom;
  Int_t truth_lambda_vtx_index = -1;
  Double_t reaction_vertex_z = TMath::QuietNaN();
  const Bool_t has_truth_lambda = FindTruthLambdaDecay(reaction_vtx_type, reaction_vtx_motherpid,
                                                        reaction_vtx_trackpid,
                                                        reaction_vtx_px, reaction_vtx_py, reaction_vtx_pz,
                                                        truth_proton_mom, truth_pion_mom, truth_lambda_vtx_index);
  const Bool_t has_reaction_vertex = FindReactionVertexZ(reaction_vtx_type, reaction_vtx_z,
                                                          reaction_vertex_z);
  if (has_truth_lambda) {
    for (Int_t isource=0; isource<kNBeamMomentumSources; ++isource) {
      const auto source = static_cast<BeamMomentumSource>(isource);
      Double_t bx = 0., by = 0., bz = 0., selected_z;
      Int_t selected_id;
      if (!SelectBeamMomentum(source, reaction_vertex_z, initial_px, initial_py, initial_pz,
                              zbeam_reco, pxbeam_reco, pybeam_reco, pzbeam_reco,
                              bx, by, bz, selected_z, selected_id)) continue;
      if (source == BeamMomentumSource::kTPCHelixReconstruction && !has_reaction_vertex) continue;
      FillTruthCandidate(truth_proton_mom, truth_pion_mom, bx, by, bz,
                         BeamMomentumSourceName(source));
    }
    // These target-effect histograms do not depend on any beam definition.
    if (truth_lambda_vtx_index >= 0 && static_cast<std::size_t>(truth_lambda_vtx_index) < reaction_vtx_x.size() && static_cast<std::size_t>(truth_lambda_vtx_index) < reaction_vtx_y.size() && static_cast<std::size_t>(truth_lambda_vtx_index) < reaction_vtx_z.size()) {
      const auto region = ClassifyTruthTargetVertex(reaction_vtx_x[truth_lambda_vtx_index], reaction_vtx_y[truth_lambda_vtx_index], reaction_vtx_z[truth_lambda_vtx_index]);
      TLorentzVector p, pi;
      p.SetXYZM(truth_proton_mom.X(), truth_proton_mom.Y(), truth_proton_mom.Z(), kMassP);
      pi.SetXYZM(truth_pion_mom.X(), truth_pion_mom.Y(), truth_pion_mom.Z(), kMassPi);
      if (region == TruthTargetRegion::kLH2Inside) HF1("Geant4LambdaEta_ReactionTruthLambdaMass_TruthVtxLH2Inside", (p+pi).M());
      if (region == TruthTargetRegion::kLH2Outside) HF1("Geant4LambdaEta_ReactionTruthLambdaMass_TruthVtxLH2Outside", (p+pi).M());
    }
  }

  const Int_t nt_tpc = **src.nt_tpc;
  for (Int_t it = 0; it < nt_tpc; ++it) {
    if (it >= static_cast<Int_t>((**src.helix_dz).size())) continue;
    const Double_t dz = (**src.helix_dz).at(it);
    if (std::isfinite(dz)) {
      HF1("Geant4LambdaEta_Track_HelixDz_All", dz);
      if ((**src.is_beam).at(it)) HF1("Geant4LambdaEta_Track_HelixDz_Beam", dz);
      if ((**src.is_accidental).at(it)) HF1("Geant4LambdaEta_Track_HelixDz_Accidental", dz);
      if (!(**src.is_beam).at(it) && !(**src.is_accidental).at(it)) HF1("Geant4LambdaEta_Track_HelixDz_Other", dz);
    }
    if (it >= static_cast<Int_t>((**src.helix_r).size()) ||
        it >= static_cast<Int_t>((**src.helix_theta_min).size()) ||
        it >= static_cast<Int_t>((**src.helix_theta_max).size())) continue;
    Double_t dx = 0.;
    if (!CalcHelixDx((**src.helix_r).at(it), (**src.helix_theta_min).at(it),
                     (**src.helix_theta_max).at(it), dx)) continue;
    HF1("Geant4LambdaEta_Track_HelixDx_All", dx);
    if ((**src.is_beam).at(it)) HF1("Geant4LambdaEta_Track_HelixDx_Beam", dx);
    if ((**src.is_accidental).at(it)) HF1("Geant4LambdaEta_Track_HelixDx_Accidental", dx);
    if (!(**src.is_beam).at(it) && !(**src.is_accidental).at(it)) HF1("Geant4LambdaEta_Track_HelixDx_Other", dx);
  }
  for (Int_t ip=0; ip<nt_tpc; ++ip) {
    if ((**src.is_beam).at(ip) || (**src.is_accidental).at(ip)) continue;
    if (!IsProton((**src.charge).at(ip), (**src.mom0).at(ip), (**src.dedx).at(ip))) continue;
    for (Int_t ipi=0; ipi<nt_tpc; ++ipi) {
      if (ip == ipi || (**src.is_beam).at(ipi) || (**src.is_accidental).at(ipi)) continue;
      if (!IsPion((**src.charge).at(ipi), (**src.mom0).at(ipi), (**src.dedx).at(ipi))) continue;

      const Double_t close_dist = (**src.close).at(ip).at(ipi);
      const Double_t z = (**src.vtx_z).at(ip).at(ipi);
      if (!IsTightCandidate(close_dist, (**src.chisqr).at(ip), (**src.chisqr).at(ipi), (**src.mom0).at(ipi))) continue;
      if (!std::isfinite(z)) continue;

      const TVector3 proton_mom((**src.mom_x).at(ip).at(ipi), (**src.mom_y).at(ip).at(ipi), (**src.mom_z).at(ip).at(ipi));
      const TVector3 pion_mom((**src.mom_x).at(ipi).at(ip), (**src.mom_y).at(ipi).at(ip), (**src.mom_z).at(ipi).at(ip));
      if (!std::isfinite(proton_mom.X()) || !std::isfinite(pion_mom.X())) continue;

      // Both reconstructed daughters must match the same direct 7202 Lambda
      // decay vertex before that truth vertex can classify this candidate.
      const Bool_t truth_matched_lambda =
        ip < static_cast<Int_t>(truth_found.size()) && ipi < static_cast<Int_t>(truth_found.size()) &&
        ip < static_cast<Int_t>(truth_source.size()) && ipi < static_cast<Int_t>(truth_source.size()) &&
        ip < static_cast<Int_t>(truth_pid.size()) && ipi < static_cast<Int_t>(truth_pid.size()) &&
        ip < static_cast<Int_t>(truth_vtx_index.size()) && ipi < static_cast<Int_t>(truth_vtx_index.size()) &&
        ip < static_cast<Int_t>(truth_vtx_type.size()) && ipi < static_cast<Int_t>(truth_vtx_type.size()) &&
        truth_found[ip] && truth_found[ipi] &&
        truth_source[ip] == 7202 && truth_source[ipi] == 7202 &&
        truth_pid[ip] == 2212 && truth_pid[ipi] == -211 &&
        truth_vtx_type[ip] == 1 && truth_vtx_type[ipi] == 1 &&
        truth_vtx_index[ip] >= 0 && truth_vtx_index[ip] == truth_vtx_index[ipi];

      // Compare the generator beam to itself: this is the beam definition used
      // by the reaction generator, not a nearest-hit approximation.
      const Double_t generator_mom = initial_mom;
      if (std::isfinite(initial_mom) && std::isfinite(generator_mom)) {
        HF1("Geant4LambdaEta_BeamMomInitial", initial_mom);
        HF1("Geant4LambdaEta_BeamMomAtVertex", generator_mom);
        HF2("Geant4LambdaEta_BeamMomAtVertex_vs_Initial", initial_mom, generator_mom);
        HF1("Geant4LambdaEta_BeamMomDiff", generator_mom - initial_mom);
      }

      for (Int_t isource=0; isource<kNBeamMomentumSources; ++isource) {
        const auto source = static_cast<BeamMomentumSource>(isource);
        Double_t bx = 0., by = 0., bz = 0., selected_z;
        Int_t selected_id;
        if (!SelectBeamMomentum(source, z, initial_px, initial_py, initial_pz,
                                zbeam_reco, pxbeam_reco, pybeam_reco, pzbeam_reco,
                                bx, by, bz, selected_z, selected_id)) continue;
        FillCandidate(close_dist, (**src.vtx_x).at(ip).at(ipi), (**src.vtx_y).at(ip).at(ipi), z,
                      proton_mom, pion_mom, bx, by, bz, selected_z, selected_id,
                      isource, BeamMomentumSourceName(source));
      }

      // Truth-beam Lambda production vertex and related diagnostics.
      if (truth_matched_lambda) {
        const Int_t ivtx = truth_vtx_index[ip];
        if (ivtx >= 0 && ivtx < static_cast<Int_t>(reaction_vtx_x.size()) &&
            ivtx < static_cast<Int_t>(reaction_vtx_y.size()) &&
            ivtx < static_cast<Int_t>(reaction_vtx_z.size())) {
          const TVector3 truth_prod(reaction_vtx_x[ivtx], reaction_vtx_y[ivtx], reaction_vtx_z[ivtx]);
          const TVector3 decay((**src.vtx_x).at(ip).at(ipi), (**src.vtx_y).at(ip).at(ipi), z);
          const TruthProduction prod = FindTruthProduction(truth_prod, TVector3(initial_px, initial_py, initial_pz), decay, proton_mom+pion_mom);
          if (prod.found) {
            HF1("Geant4LambdaEta_ProductionVtxDCA", prod.dca);
            HF2("Geant4LambdaEta_ProductionVtxZX", prod.pos.Z(), prod.pos.X());
            HF2("Geant4LambdaEta_ProductionVtxZY", prod.pos.Z(), prod.pos.Y());
            HF2("Geant4LambdaEta_ProductionVtxXY", prod.pos.X(), prod.pos.Y());
            TLorentzVector beam_truth, target, lam;
            beam_truth.SetXYZM(initial_px, initial_py, initial_pz, pdg::KaonMass());
            target.SetXYZM(0., 0., 0., pdg::ProtonMass());
            lam.SetXYZM((proton_mom+pion_mom).X(), (proton_mom+pion_mom).Y(), (proton_mom+pion_mom).Z(), kMassLambda);
            HF1("Geant4LambdaEta_MissingMass_ProductionVtx", (beam_truth+target-lam).M());
          }
        }
      }

      // Unlike the regular reconstruction spectra above, this has no beam
      // source dimension: Lambda mass and its truth target vertex are beam-independent.
      if (truth_matched_lambda) {
        const Int_t ivtx = truth_vtx_index[ip];
        if (ivtx < static_cast<Int_t>(reaction_vtx_x.size()) &&
            ivtx < static_cast<Int_t>(reaction_vtx_y.size()) &&
            ivtx < static_cast<Int_t>(reaction_vtx_z.size())) {
          const auto region = ClassifyTruthTargetVertex(reaction_vtx_x[ivtx], reaction_vtx_y[ivtx], reaction_vtx_z[ivtx]);
          TLorentzVector p, pi;
          p.SetXYZM(proton_mom.X(), proton_mom.Y(), proton_mom.Z(), kMassP);
          pi.SetXYZM(pion_mom.X(), pion_mom.Y(), pion_mom.Z(), kMassPi);
          if (region == TruthTargetRegion::kLH2Inside) HF1("Geant4LambdaEta_ReactionReconLambdaMass_TruthVtxLH2Inside", (p+pi).M());
          if (region == TruthTargetRegion::kLH2Outside) HF1("Geant4LambdaEta_ReactionReconLambdaMass_TruthVtxLH2Outside", (p+pi).M());
        }
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
  for (Int_t i=0; i<TFileCont.size(); ++i) {
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
  if (!dst::SetupReader(kGeant4Tpc, "kGeant4Tpc")) return false;
  const auto& reader = TTreeReaderCont[kGeant4Tpc];
  dst::SetBranch(reader, "event_number", src.event_number);
  dst::SetBranch(reader, "trig_flag_g4", src.trig_flag);
  dst::SetBranch(reader, "ntTpc", src.nt_tpc);
  dst::SetBranch(reader, "beam_reco_found", src.beam_found);
  dst::SetBranch(reader, "is_beam", src.is_beam);
  dst::SetBranch(reader, "is_accidental", src.is_accidental);
  dst::SetBranch(reader, "charge", src.charge);
  dst::SetBranch(reader, "chisqr", src.chisqr);
  dst::SetBranch(reader, "dEdx", src.dedx);
  dst::SetBranch(reader, "mom0", src.mom0);
  dst::SetBranch(reader, "dE", src.dE);
  dst::SetBranch(reader, "mom0_x", src.mom0_x);
  dst::SetBranch(reader, "mom0_y", src.mom0_y);
  dst::SetBranch(reader, "mom0_z", src.mom0_z);
  dst::SetBranch(reader, "pid", src.pid);
  dst::SetBranch(reader, "helix_r", src.helix_r);
  dst::SetBranch(reader, "helix_dz", src.helix_dz);
  dst::SetBranch(reader, "helix_theta_min", src.helix_theta_min);
  dst::SetBranch(reader, "helix_theta_max", src.helix_theta_max);
  dst::SetBranch(reader, "closeDistTpc", src.close);
  dst::SetBranch(reader, "vtxTpc", src.vtx_x);
  dst::SetBranch(reader, "vtyTpc", src.vtx_y);
  dst::SetBranch(reader, "vtzTpc", src.vtx_z);
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
  dst::SetBranch(reader, "track_cluster_de", src.track_cluster_de);
  dst::SetBranch(reader, "track_cluster_size", src.track_cluster_size);
  dst::SetBranch(reader, "track_cluster_mrow", src.track_cluster_mrow);
  dst::SetBranch(reader, "track_cluster_de_center", src.track_cluster_de_center);
  dst::SetBranch(reader, "track_cluster_x_center", src.track_cluster_x_center);
  dst::SetBranch(reader, "track_cluster_y_center", src.track_cluster_y_center);
  dst::SetBranch(reader, "track_cluster_z_center", src.track_cluster_z_center);
  dst::SetBranch(reader, "track_cluster_row_center", src.track_cluster_row_center);
  dst::SetBranch(reader, "mom_vtx", src.mom_x);
  dst::SetBranch(reader, "mom_vty", src.mom_y);
  dst::SetBranch(reader, "mom_vtz", src.mom_z);
  dst::SetBranch(reader, "beam_initial_px_g4", src.beam_initial_px_g4);
  dst::SetBranch(reader, "beam_initial_py_g4", src.beam_initial_py_g4);
  dst::SetBranch(reader, "beam_initial_pz_g4", src.beam_initial_pz_g4);
  dst::SetBranch(reader, "reaction_beam_px_g4", src.reaction_beam_px_g4);
  dst::SetBranch(reader, "reaction_beam_py_g4", src.reaction_beam_py_g4);
  dst::SetBranch(reader, "reaction_beam_pz_g4", src.reaction_beam_pz_g4);
  dst::SetBranch(reader, "beam_reco_hit_z", src.beam_reco_hit_z);
  dst::SetBranch(reader, "beam_reco_px", src.beam_reco_px);
  dst::SetBranch(reader, "beam_reco_py", src.beam_reco_py);
  dst::SetBranch(reader, "beam_reco_pz", src.beam_reco_pz);
  dst::SetBranch(reader, "reaction_vtx_type_g4", src.reaction_vtx_type);
  dst::SetBranch(reader, "reaction_vtx_motherpid_g4", src.reaction_vtx_motherpid);
  dst::SetBranch(reader, "reaction_vtx_x_g4", src.reaction_vtx_x);
  dst::SetBranch(reader, "reaction_vtx_y_g4", src.reaction_vtx_y);
  dst::SetBranch(reader, "reaction_vtx_z_g4", src.reaction_vtx_z);
  dst::SetBranch(reader, "truth_found", src.truth_found);
  dst::SetBranch(reader, "truth_source_g4", src.truth_source);
  dst::SetBranch(reader, "truth_pid_g4", src.truth_pid);
  dst::SetBranch(reader, "truth_vtx_index_g4", src.truth_vtx_index);
  dst::SetBranch(reader, "truth_vtx_type_g4", src.truth_vtx_type);
  dst::SetBranch(reader, "reaction_vtx_trackpid_g4", src.reaction_vtx_trackpid);
  dst::SetBranch(reader, "reaction_vtx_px_g4", src.reaction_vtx_px);
  dst::SetBranch(reader, "reaction_vtx_py_g4", src.reaction_vtx_py);
  dst::SetBranch(reader, "reaction_vtx_pz_g4", src.reaction_vtx_pz);
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return InitializeParameter<TPCParamMan>("TPCPRM") &&
    InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
    InitializeParameter<UserParamMan>("USER");
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  TPCEventAnalyzer::SetDstCalibFlag(false);
  hist::BuildStatus();
  hist::BuildTPCLambdaEta();
  hist::BuildTPCGeant4LambdaEta();

  tree = new TTree("tpc", "tree of DstTPCGeant4LambdaEta");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_flag_g4", &event.trigflag);
  tree->Branch("beam_reco_found", &event.beam_reco_found);
  tree->Branch("lambda_mass", &event.lambda_mass);
  tree->Branch("lambda_close_dist", &event.lambda_close_dist);
  tree->Branch("lambda_vtx_x", &event.lambda_vtx_x);
  tree->Branch("lambda_vtx_y", &event.lambda_vtx_y);
  tree->Branch("lambda_vtx_z", &event.lambda_vtx_z);
  tree->Branch("lambda_mom_x", &event.lambda_mom_x);
  tree->Branch("lambda_mom_y", &event.lambda_mom_y);
  tree->Branch("lambda_mom_z", &event.lambda_mom_z);
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
  tree->Branch("X_beam_source", &event.X_beam_source);
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
