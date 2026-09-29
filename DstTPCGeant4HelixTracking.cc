// -*- C++ -*-

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCVertex.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>

namespace
{
  auto& gConf = ConfMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();

  constexpr Int_t kBeamGenerator = 7201;

  // Change only this line for analyzer tests with a single generator entry.
  enum class Geant4InputMode { kBeamOnly, kReactionOnly, kBeamAndReaction };
  constexpr Geant4InputMode kGeant4InputMode = Geant4InputMode::kBeamAndReaction;

  struct InputBranches
  {
    Int_t evnum = -1;
    Int_t effective_evnum = -1;
    Int_t generator = -1;
    Int_t trig_flag = 0;
    std::vector<Int_t>* pad = nullptr;
    std::vector<Int_t>* layer = nullptr;
    std::vector<Int_t>* row = nullptr;
    std::vector<Int_t>* trackid = nullptr;
    std::vector<Int_t>* pid = nullptr;
    std::vector<Int_t>* parentid = nullptr;
    std::vector<Int_t>* parentpid = nullptr;
    std::vector<Double_t>* x = nullptr;
    std::vector<Double_t>* y = nullptr;
    std::vector<Double_t>* z = nullptr;
    std::vector<Double_t>* de = nullptr;
    std::vector<Double_t>* px = nullptr;
    std::vector<Double_t>* py = nullptr;
    std::vector<Double_t>* pz = nullptr;
    std::vector<TParticle>* beam_particles = nullptr;
    std::vector<Int_t>* vtx_type = nullptr;
    std::vector<Int_t>* vtx_motherpid = nullptr;
    std::vector<Double_t>* vtx_x = nullptr;
    std::vector<Double_t>* vtx_y = nullptr;
    std::vector<Double_t>* vtx_z = nullptr;
    std::vector<std::vector<Int_t>>* vtx_trackid = nullptr;
    std::vector<std::vector<Int_t>>* vtx_trackpid = nullptr;
    std::vector<std::vector<Double_t>>* vtx_px = nullptr;
    std::vector<std::vector<Double_t>>* vtx_py = nullptr;
    std::vector<std::vector<Double_t>>* vtx_pz = nullptr;
  } src;

  struct GeantClusters
  {
    Int_t evnum = -1;
    Int_t effective_evnum = -1;
    Int_t trig_flag = 0;
    Int_t generator = -1;
    Int_t initial_beam_pid = 0;
    Double_t initial_beam_px = TMath::QuietNaN();
    Double_t initial_beam_py = TMath::QuietNaN();
    Double_t initial_beam_pz = TMath::QuietNaN();
    std::vector<Int_t> pad, layer, row, trackid, pid, parentid, parentpid;
    std::vector<Double_t> x, y, z, de, px, py, pz;
    std::vector<Int_t> vtx_type, vtx_motherpid;
    std::vector<Double_t> vtx_x, vtx_y, vtx_z;
    std::vector<std::vector<Int_t>> vtx_trackid, vtx_trackpid;
    std::vector<std::vector<Double_t>> vtx_px, vtx_py, vtx_pz;

    void clear()
    {
      evnum = effective_evnum = -1; trig_flag = 0; generator = -1;
      initial_beam_pid = 0;
      initial_beam_px = initial_beam_py = initial_beam_pz = TMath::QuietNaN();
      pad.clear(); layer.clear(); row.clear(); trackid.clear(); pid.clear();
      parentid.clear(); parentpid.clear(); x.clear(); y.clear(); z.clear(); de.clear();
      px.clear(); py.clear(); pz.clear();
      vtx_type.clear(); vtx_motherpid.clear(); vtx_x.clear(); vtx_y.clear(); vtx_z.clear();
      vtx_trackid.clear(); vtx_trackpid.clear();
      vtx_px.clear(); vtx_py.clear(); vtx_pz.clear();

    }

    void copyInitialBeam(const InputBranches& in)
    {
      if (!in.beam_particles || in.beam_particles->empty()) return;
      const auto& beam = in.beam_particles->front();
      initial_beam_pid = beam.GetPdgCode();
      initial_beam_px = beam.Px();
      initial_beam_py = beam.Py();
      initial_beam_pz = beam.Pz();
    }

    void append(const InputBranches& in, Bool_t upstream_only)
    {
      const auto nhits = in.pad ? in.pad->size() : 0;
      for (std::size_t ih = 0; ih < nhits; ++ih) {
        if (upstream_only && in.z->at(ih) >= tpc::Z_TARGET) continue;
        pad.push_back(in.pad->at(ih));
        layer.push_back(in.layer->at(ih));
        row.push_back(in.row->at(ih));
        trackid.push_back(in.trackid->at(ih));
        pid.push_back(in.pid->at(ih));
        parentid.push_back(in.parentid->at(ih));
        parentpid.push_back(in.parentpid->at(ih));
        x.push_back(in.x->at(ih));
        y.push_back(in.y->at(ih));
        z.push_back(in.z->at(ih));
        de.push_back(in.de->at(ih));
        px.push_back(in.px->at(ih));
        py.push_back(in.py->at(ih));
        pz.push_back(in.pz->at(ih));
      }

    }
  };

  void SetVertices(GeantClusters& dst, const InputBranches& in)
  {
    dst.vtx_type = *in.vtx_type;
    dst.vtx_motherpid = *in.vtx_motherpid;
    dst.vtx_x = *in.vtx_x; dst.vtx_y = *in.vtx_y; dst.vtx_z = *in.vtx_z;
    dst.vtx_trackid = *in.vtx_trackid; dst.vtx_trackpid = *in.vtx_trackpid;
    dst.vtx_px = *in.vtx_px; dst.vtx_py = *in.vtx_py; dst.vtx_pz = *in.vtx_pz;
  }

  struct Event
  {
    Int_t status = 0;
    Int_t event_number = -1; // Geant4 effective_evnum: one entry per beam+reaction pair
    Int_t beam_event_number = -1;
    Int_t reaction_event_number = -1;
    Int_t trig_flag = 0;
    Int_t nhTpc = 0;
    std::vector<Int_t> raw_source, raw_padid, raw_layer, raw_row;
    std::vector<Int_t> raw_trackid, raw_pid, raw_parentid, raw_parentpid;
    std::vector<Double_t> raw_hitpos_x, raw_hitpos_y, raw_hitpos_z, raw_de;
    std::vector<Double_t> raw_truth_px, raw_truth_py, raw_truth_pz;

    // Initial 7201 BEAM[0] truth before transport through the apparatus.
    Int_t beam_initial_pid_g4 = 0;
    Double_t beam_initial_px_g4 = TMath::QuietNaN();
    Double_t beam_initial_py_g4 = TMath::QuietNaN();
    Double_t beam_initial_pz_g4 = TMath::QuietNaN();
    // Beam four-vector used by the 7202 reaction generator (BEAM[0]).
    Int_t reaction_beam_pid_g4 = 0;
    Double_t reaction_beam_px_g4 = TMath::QuietNaN();
    Double_t reaction_beam_py_g4 = TMath::QuietNaN();
    Double_t reaction_beam_pz_g4 = TMath::QuietNaN();

    // Original 7201 beam TPC-hit truth, retained independently from reaction tracks.
    Int_t beam_truth_trackid = -1, beam_truth_pid = 0;
    std::vector<Double_t> beam_truth_x, beam_truth_y, beam_truth_z;
    std::vector<Double_t> beam_truth_px, beam_truth_py, beam_truth_pz;

    // A second common helix tracking pass using only all-z 7201 beam hits.
    Int_t beam_reco_found = 0, beam_reco_nhit = 0, beam_reco_charge = 0;
    Double_t beam_reco_chisqr = TMath::QuietNaN();
    Double_t beam_reco_cx = TMath::QuietNaN(), beam_reco_cy = TMath::QuietNaN();
    Double_t beam_reco_z0 = TMath::QuietNaN(), beam_reco_r = TMath::QuietNaN(), beam_reco_dz = TMath::QuietNaN();
    std::vector<Double_t> beam_reco_hit_x, beam_reco_hit_y, beam_reco_hit_z;
    std::vector<Double_t> beam_reco_cal_x, beam_reco_cal_y, beam_reco_cal_z;
    std::vector<Double_t> beam_reco_px, beam_reco_py, beam_reco_pz;

    // Original Geant4 reaction-vertex record, independent of any TPC match.
    std::vector<Int_t> reaction_vtx_type_g4, reaction_vtx_motherpid_g4;
    std::vector<Double_t> reaction_vtx_x_g4, reaction_vtx_y_g4, reaction_vtx_z_g4;
    std::vector<std::vector<Double_t>> reaction_vtx_trackpid_g4;
    std::vector<std::vector<Double_t>> reaction_vtx_px_g4, reaction_vtx_py_g4, reaction_vtx_pz_g4;

    Int_t ntTpc = 0;
    std::vector<Int_t> nhtrack, is_beam, is_accidental, charge, pid;
    std::vector<Double_t> chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz;
    std::vector<Double_t> helix_theta_min, helix_theta_max;
    std::vector<Double_t> mom0_x, mom0_y, mom0_z, mom0, dE, dEdx;
    std::vector<std::vector<Double_t>> closeDistTpc, vtxTpc, vtyTpc, vtzTpc;
    std::vector<std::vector<Double_t>> mom_vtx, mom_vty, mom_vtz;
    // Geant4 vertex truth matched to each reconstructed helix track.
    std::vector<Int_t> truth_found, truth_source, truth_trackid, truth_pid, truth_vtx_index, truth_vtx_type, truth_match_nhit;
    std::vector<Double_t> truth_match_purity;
    std::vector<Double_t> truth_vtx_x, truth_vtx_y, truth_vtx_z;
    std::vector<Double_t> truth_vtx_px, truth_vtx_py, truth_vtx_pz;
    std::vector<Double_t> reco_vtx_px, reco_vtx_py, reco_vtx_pz;
    std::vector<Double_t> delta_vtx_px, delta_vtx_py, delta_vtx_pz, delta_vtx_p, delta_vtx_pt;
    std::vector<std::vector<Double_t>> hitlayer, hitpos_x, hitpos_y, hitpos_z;
    std::vector<std::vector<Double_t>> calpos_x, calpos_y, calpos_z;
    std::vector<std::vector<Double_t>> residual, residual_x, residual_y, residual_z;
    std::vector<std::vector<Double_t>> track_cluster_de, track_cluster_size, track_cluster_mrow;
    std::vector<std::vector<Double_t>> track_cluster_de_center;
    std::vector<std::vector<Double_t>> track_cluster_x_center, track_cluster_y_center;
    std::vector<std::vector<Double_t>> track_cluster_z_center, track_cluster_row_center;

    void clear()
    {
      status = 0; event_number = beam_event_number = reaction_event_number = -1; trig_flag = 0;
      nhTpc = 0; ntTpc = 0;
      raw_source.clear(); raw_padid.clear(); raw_layer.clear(); raw_row.clear();
      raw_trackid.clear(); raw_pid.clear(); raw_parentid.clear(); raw_parentpid.clear();
      raw_hitpos_x.clear(); raw_hitpos_y.clear(); raw_hitpos_z.clear(); raw_de.clear();
      raw_truth_px.clear(); raw_truth_py.clear(); raw_truth_pz.clear();
      beam_initial_pid_g4 = 0;
      beam_initial_px_g4 = beam_initial_py_g4 = beam_initial_pz_g4 = TMath::QuietNaN();
      reaction_beam_pid_g4 = 0;
      reaction_beam_px_g4 = reaction_beam_py_g4 = reaction_beam_pz_g4 = TMath::QuietNaN();
      beam_truth_trackid = -1; beam_truth_pid = 0;
      beam_truth_x.clear(); beam_truth_y.clear(); beam_truth_z.clear();
      beam_truth_px.clear(); beam_truth_py.clear(); beam_truth_pz.clear();
      beam_reco_found = beam_reco_nhit = beam_reco_charge = 0;
      beam_reco_chisqr = TMath::QuietNaN();
      beam_reco_cx = beam_reco_cy = beam_reco_z0 = beam_reco_r = beam_reco_dz = TMath::QuietNaN();
      beam_reco_hit_x.clear(); beam_reco_hit_y.clear(); beam_reco_hit_z.clear();
      beam_reco_cal_x.clear(); beam_reco_cal_y.clear(); beam_reco_cal_z.clear();
      beam_reco_px.clear(); beam_reco_py.clear(); beam_reco_pz.clear();
      reaction_vtx_type_g4.clear(); reaction_vtx_motherpid_g4.clear();
      reaction_vtx_x_g4.clear(); reaction_vtx_y_g4.clear(); reaction_vtx_z_g4.clear();
      reaction_vtx_trackpid_g4.clear();
      reaction_vtx_px_g4.clear(); reaction_vtx_py_g4.clear(); reaction_vtx_pz_g4.clear();
      nhtrack.clear(); is_beam.clear(); is_accidental.clear(); charge.clear(); pid.clear();
      chisqr.clear(); helix_cx.clear(); helix_cy.clear(); helix_z0.clear(); helix_r.clear(); helix_dz.clear();
      helix_theta_min.clear(); helix_theta_max.clear();
      mom0_x.clear(); mom0_y.clear(); mom0_z.clear(); mom0.clear(); dE.clear(); dEdx.clear();
      closeDistTpc.clear(); vtxTpc.clear(); vtyTpc.clear(); vtzTpc.clear();
      mom_vtx.clear(); mom_vty.clear(); mom_vtz.clear();
      truth_found.clear(); truth_source.clear(); truth_trackid.clear(); truth_pid.clear(); truth_vtx_index.clear(); truth_vtx_type.clear(); truth_match_nhit.clear();
      truth_match_purity.clear(); truth_vtx_x.clear(); truth_vtx_y.clear(); truth_vtx_z.clear();
      truth_vtx_px.clear(); truth_vtx_py.clear(); truth_vtx_pz.clear();
      reco_vtx_px.clear(); reco_vtx_py.clear(); reco_vtx_pz.clear();
      delta_vtx_px.clear(); delta_vtx_py.clear(); delta_vtx_pz.clear(); delta_vtx_p.clear(); delta_vtx_pt.clear();
      hitlayer.clear(); hitpos_x.clear(); hitpos_y.clear(); hitpos_z.clear();
      calpos_x.clear(); calpos_y.clear(); calpos_z.clear();
      residual.clear(); residual_x.clear(); residual_y.clear(); residual_z.clear();
      track_cluster_de.clear(); track_cluster_size.clear(); track_cluster_mrow.clear();
      track_cluster_de_center.clear(); track_cluster_x_center.clear(); track_cluster_y_center.clear();
      track_cluster_z_center.clear(); track_cluster_row_center.clear();

    }
  } event;

  template <class T>
  Bool_t SetBranch(TChain& tree, const char* name, T* address)
  {
    if (!tree.GetBranch(name)) {
      std::cerr << "missing Geant4 branch: " << name << std::endl;
      return false;

    }
    tree.SetBranchAddress(name, address);
    return true;
  }

  Bool_t SetupInput(TChain& tree)
  {
    return SetBranch(tree, "evnum", &src.evnum)
      && SetBranch(tree, "effective_evnum", &src.effective_evnum)
      && SetBranch(tree, "generator", &src.generator)
      && SetBranch(tree, "trig_flag", &src.trig_flag)
      && SetBranch(tree, "padtpc", &src.pad)
      && SetBranch(tree, "layertpc", &src.layer)
      && SetBranch(tree, "rowtpc", &src.row)
      && SetBranch(tree, "trackidtpc", &src.trackid)
      && SetBranch(tree, "pidtpc", &src.pid)
      && SetBranch(tree, "parentidtpc", &src.parentid)
      && SetBranch(tree, "parentpidtpc", &src.parentpid)
      && SetBranch(tree, "xtpc_pad", &src.x)
      && SetBranch(tree, "ytpc_pad", &src.y)
      && SetBranch(tree, "ztpc_pad", &src.z)
      && SetBranch(tree, "edeptpc", &src.de)
      && SetBranch(tree, "pxtpc", &src.px)
      && SetBranch(tree, "pytpc", &src.py)
      && SetBranch(tree, "pztpc", &src.pz)
      && SetBranch(tree, "BEAM", &src.beam_particles)
      && SetBranch(tree, "vtx_type", &src.vtx_type)
      && SetBranch(tree, "vtx_motherpid", &src.vtx_motherpid)
      && SetBranch(tree, "vtx_x", &src.vtx_x)
      && SetBranch(tree, "vtx_y", &src.vtx_y)
      && SetBranch(tree, "vtx_z", &src.vtx_z)
      && SetBranch(tree, "vtx_trackid", &src.vtx_trackid)
      && SetBranch(tree, "vtx_trackpid", &src.vtx_trackpid)
      && SetBranch(tree, "vtx_px", &src.vtx_px)
      && SetBranch(tree, "vtx_py", &src.vtx_py)
      && SetBranch(tree, "vtx_pz", &src.vtx_pz);
  }

  void BookTree(TTree& tree)
  {
    tree.Branch("status", &event.status);
    tree.Branch("event_number", &event.event_number);
    tree.Branch("beam_event_number", &event.beam_event_number);
    tree.Branch("reaction_event_number", &event.reaction_event_number);
    tree.Branch("trig_flag_g4", &event.trig_flag);
    tree.Branch("nhTpc", &event.nhTpc);
    tree.Branch("raw_source", &event.raw_source); // original Geant4 generator number
    tree.Branch("raw_padid", &event.raw_padid);
    tree.Branch("raw_layer", &event.raw_layer);
    tree.Branch("raw_row", &event.raw_row);
    tree.Branch("raw_trackid", &event.raw_trackid);
    tree.Branch("raw_pid", &event.raw_pid);
    tree.Branch("raw_parentid", &event.raw_parentid);
    tree.Branch("raw_parentpid", &event.raw_parentpid);
    tree.Branch("raw_hitpos_x", &event.raw_hitpos_x);
    tree.Branch("raw_hitpos_y", &event.raw_hitpos_y);
    tree.Branch("raw_hitpos_z", &event.raw_hitpos_z);
    tree.Branch("raw_de", &event.raw_de);
    tree.Branch("raw_px_g4", &event.raw_truth_px);
    tree.Branch("raw_py_g4", &event.raw_truth_py);
    tree.Branch("raw_pz_g4", &event.raw_truth_pz);
    tree.Branch("beam_initial_pid_g4", &event.beam_initial_pid_g4);
    tree.Branch("beam_initial_px_g4", &event.beam_initial_px_g4);
    tree.Branch("beam_initial_py_g4", &event.beam_initial_py_g4);
    tree.Branch("beam_initial_pz_g4", &event.beam_initial_pz_g4);
    tree.Branch("reaction_beam_pid_g4", &event.reaction_beam_pid_g4);
    tree.Branch("reaction_beam_px_g4", &event.reaction_beam_px_g4);
    tree.Branch("reaction_beam_py_g4", &event.reaction_beam_py_g4);
    tree.Branch("reaction_beam_pz_g4", &event.reaction_beam_pz_g4);
    tree.Branch("beam_trackid_g4", &event.beam_truth_trackid);
    tree.Branch("beam_pid_g4", &event.beam_truth_pid);
    tree.Branch("beam_x_g4", &event.beam_truth_x);
    tree.Branch("beam_y_g4", &event.beam_truth_y);
    tree.Branch("beam_z_g4", &event.beam_truth_z);
    tree.Branch("beam_px_g4", &event.beam_truth_px);
    tree.Branch("beam_py_g4", &event.beam_truth_py);
    tree.Branch("beam_pz_g4", &event.beam_truth_pz);
    tree.Branch("beam_reco_found", &event.beam_reco_found);
    tree.Branch("beam_reco_nhit", &event.beam_reco_nhit);
    tree.Branch("beam_reco_charge", &event.beam_reco_charge);
    tree.Branch("beam_reco_chisqr", &event.beam_reco_chisqr);
    tree.Branch("beam_reco_cx", &event.beam_reco_cx);
    tree.Branch("beam_reco_cy", &event.beam_reco_cy);
    tree.Branch("beam_reco_z0", &event.beam_reco_z0);
    tree.Branch("beam_reco_r", &event.beam_reco_r);
    tree.Branch("beam_reco_dz", &event.beam_reco_dz);
    tree.Branch("beam_reco_hit_x", &event.beam_reco_hit_x);
    tree.Branch("beam_reco_hit_y", &event.beam_reco_hit_y);
    tree.Branch("beam_reco_hit_z", &event.beam_reco_hit_z);
    tree.Branch("beam_reco_cal_x", &event.beam_reco_cal_x);
    tree.Branch("beam_reco_cal_y", &event.beam_reco_cal_y);
    tree.Branch("beam_reco_cal_z", &event.beam_reco_cal_z);
    tree.Branch("beam_reco_px", &event.beam_reco_px);
    tree.Branch("beam_reco_py", &event.beam_reco_py);
    tree.Branch("beam_reco_pz", &event.beam_reco_pz);
    tree.Branch("reaction_vtx_type_g4", &event.reaction_vtx_type_g4);
    tree.Branch("reaction_vtx_motherpid_g4", &event.reaction_vtx_motherpid_g4);
    tree.Branch("reaction_vtx_x_g4", &event.reaction_vtx_x_g4);
    tree.Branch("reaction_vtx_y_g4", &event.reaction_vtx_y_g4);
    tree.Branch("reaction_vtx_z_g4", &event.reaction_vtx_z_g4);
    tree.Branch("reaction_vtx_trackpid_g4", &event.reaction_vtx_trackpid_g4);
    tree.Branch("reaction_vtx_px_g4", &event.reaction_vtx_px_g4);
    tree.Branch("reaction_vtx_py_g4", &event.reaction_vtx_py_g4);
    tree.Branch("reaction_vtx_pz_g4", &event.reaction_vtx_pz_g4);
    tree.Branch("ntTpc", &event.ntTpc);
    tree.Branch("nhtrack", &event.nhtrack);
    tree.Branch("is_beam", &event.is_beam);
    tree.Branch("is_accidental", &event.is_accidental);
    tree.Branch("charge", &event.charge);
    tree.Branch("pid", &event.pid);
    tree.Branch("chisqr", &event.chisqr);
    tree.Branch("helix_cx", &event.helix_cx);
    tree.Branch("helix_cy", &event.helix_cy);
    tree.Branch("helix_z0", &event.helix_z0);
    tree.Branch("helix_r", &event.helix_r);
    tree.Branch("helix_dz", &event.helix_dz);
    tree.Branch("helix_theta_min", &event.helix_theta_min);
    tree.Branch("helix_theta_max", &event.helix_theta_max);
    tree.Branch("mom0_x", &event.mom0_x);
    tree.Branch("mom0_y", &event.mom0_y);
    tree.Branch("mom0_z", &event.mom0_z);
    tree.Branch("mom0", &event.mom0);
    tree.Branch("dE", &event.dE);
    tree.Branch("dEdx", &event.dEdx);
    tree.Branch("closeDistTpc", &event.closeDistTpc);
    tree.Branch("vtxTpc", &event.vtxTpc);
    tree.Branch("vtyTpc", &event.vtyTpc);
    tree.Branch("vtzTpc", &event.vtzTpc);
    tree.Branch("mom_vtx", &event.mom_vtx);
    tree.Branch("mom_vty", &event.mom_vty);
    tree.Branch("mom_vtz", &event.mom_vtz);
    tree.Branch("truth_found", &event.truth_found);
    tree.Branch("truth_source_g4", &event.truth_source);
    tree.Branch("truth_trackid_g4", &event.truth_trackid);
    tree.Branch("truth_pid_g4", &event.truth_pid);
    tree.Branch("truth_vtx_index_g4", &event.truth_vtx_index);
    tree.Branch("truth_vtx_type_g4", &event.truth_vtx_type);
    tree.Branch("truth_match_nhit", &event.truth_match_nhit);
    tree.Branch("truth_match_purity", &event.truth_match_purity);
    tree.Branch("truth_vtx_x_g4", &event.truth_vtx_x);
    tree.Branch("truth_vtx_y_g4", &event.truth_vtx_y);
    tree.Branch("truth_vtx_z_g4", &event.truth_vtx_z);
    tree.Branch("truth_vtx_px_g4", &event.truth_vtx_px);
    tree.Branch("truth_vtx_py_g4", &event.truth_vtx_py);
    tree.Branch("truth_vtx_pz_g4", &event.truth_vtx_pz);
    tree.Branch("reco_vtx_px", &event.reco_vtx_px);
    tree.Branch("reco_vtx_py", &event.reco_vtx_py);
    tree.Branch("reco_vtx_pz", &event.reco_vtx_pz);
    tree.Branch("delta_vtx_px", &event.delta_vtx_px);
    tree.Branch("delta_vtx_py", &event.delta_vtx_py);
    tree.Branch("delta_vtx_pz", &event.delta_vtx_pz);
    tree.Branch("delta_vtx_p", &event.delta_vtx_p);
    tree.Branch("delta_vtx_pt", &event.delta_vtx_pt);
    tree.Branch("hitlayer", &event.hitlayer);
    tree.Branch("hitpos_x", &event.hitpos_x);
    tree.Branch("hitpos_y", &event.hitpos_y);
    tree.Branch("hitpos_z", &event.hitpos_z);
    tree.Branch("calpos_x", &event.calpos_x);
    tree.Branch("calpos_y", &event.calpos_y);
    tree.Branch("calpos_z", &event.calpos_z);
    tree.Branch("residual", &event.residual);
    tree.Branch("residual_x", &event.residual_x);
    tree.Branch("residual_y", &event.residual_y);
    tree.Branch("residual_z", &event.residual_z);
    tree.Branch("track_cluster_de", &event.track_cluster_de);
    tree.Branch("track_cluster_size", &event.track_cluster_size);
    tree.Branch("track_cluster_mrow", &event.track_cluster_mrow);
    tree.Branch("track_cluster_de_center", &event.track_cluster_de_center);
    tree.Branch("track_cluster_x_center", &event.track_cluster_x_center);
    tree.Branch("track_cluster_y_center", &event.track_cluster_y_center);
    tree.Branch("track_cluster_z_center", &event.track_cluster_z_center);
    tree.Branch("track_cluster_row_center", &event.track_cluster_row_center);
  }

  void FillRaw(const GeantClusters& beam, const GeantClusters& reaction)
  {
    const auto append = [](const GeantClusters& in, Int_t source, Bool_t upstream_only) {
      for (std::size_t ih = 0; ih < in.pad.size(); ++ih) {
        if (upstream_only && in.z[ih] >= tpc::Z_TARGET) continue;
        event.raw_source.push_back(source);
        event.raw_padid.push_back(in.pad[ih]); event.raw_layer.push_back(in.layer[ih]); event.raw_row.push_back(in.row[ih]);
        event.raw_trackid.push_back(in.trackid[ih]); event.raw_pid.push_back(in.pid[ih]);
        event.raw_parentid.push_back(in.parentid[ih]); event.raw_parentpid.push_back(in.parentpid[ih]);
        event.raw_hitpos_x.push_back(in.x[ih]); event.raw_hitpos_y.push_back(in.y[ih]); event.raw_hitpos_z.push_back(in.z[ih]);
        event.raw_de.push_back(in.de[ih]);
        event.raw_truth_px.push_back(in.px[ih]); event.raw_truth_py.push_back(in.py[ih]); event.raw_truth_pz.push_back(in.pz[ih]);
      }
    };
    append(beam, beam.generator, kGeant4InputMode == Geant4InputMode::kBeamAndReaction);
    append(reaction, reaction.generator, false);
    event.nhTpc = event.raw_padid.size();
  }

  void FillInitialBeamTruth(const GeantClusters& beam)
  {
    event.beam_initial_pid_g4 = beam.initial_beam_pid;
    event.beam_initial_px_g4 = beam.initial_beam_px;
    event.beam_initial_py_g4 = beam.initial_beam_py;
    event.beam_initial_pz_g4 = beam.initial_beam_pz;
  }

  void FillReactionBeamTruth(const GeantClusters& reaction)
  {
    // For a 7202 event, BEAM[0] is the momentum passed to the reaction
    // generator at its production vertex; unlike beam_initial_*, it is not
    // the upstream 7201 beam event.
    event.reaction_beam_pid_g4 = reaction.initial_beam_pid;
    event.reaction_beam_px_g4 = reaction.initial_beam_px;
    event.reaction_beam_py_g4 = reaction.initial_beam_py;
    event.reaction_beam_pz_g4 = reaction.initial_beam_pz;
  }

  void FillBeamTruth(const GeantClusters& beam)
  {
    // Prefer the primary K-; fall back to the most populated truth track.
    std::map<Int_t, Int_t> count;
    for (std::size_t ih = 0; ih < beam.trackid.size(); ++ih)
      if (beam.pid[ih] == -321) ++count[beam.trackid[ih]];
    if (count.empty())
      for (const Int_t trackid : beam.trackid) ++count[trackid];
    if (count.empty()) return;

    const auto best = std::max_element(
      count.begin(), count.end(),
      [](const auto& a, const auto& b) { return a.second < b.second; });
    event.beam_truth_trackid = best->first;
    for (std::size_t ih = 0; ih < beam.trackid.size(); ++ih) {
      if (beam.trackid[ih] != event.beam_truth_trackid) continue;
      event.beam_truth_pid = beam.pid[ih];
      event.beam_truth_x.push_back(beam.x[ih]);
      event.beam_truth_y.push_back(beam.y[ih]);
      event.beam_truth_z.push_back(beam.z[ih]);
      event.beam_truth_px.push_back(beam.px[ih]);
      event.beam_truth_py.push_back(beam.py[ih]);
      event.beam_truth_pz.push_back(beam.pz[ih]);

    }
  }

  void FillReactionVertexTruth(const GeantClusters& reaction)
  {
    // Preserve the generator vertex record verbatim.  Unlike truth_vtx_* this
    // is not selected through a reconstructed TPC helix or hit-vote match.
    event.reaction_vtx_type_g4 = reaction.vtx_type;
    event.reaction_vtx_motherpid_g4 = reaction.vtx_motherpid;
    event.reaction_vtx_x_g4 = reaction.vtx_x;
    event.reaction_vtx_y_g4 = reaction.vtx_y;
    event.reaction_vtx_z_g4 = reaction.vtx_z;
    event.reaction_vtx_trackpid_g4.clear();
    event.reaction_vtx_trackpid_g4.reserve(reaction.vtx_trackpid.size());
    for (const auto& pid : reaction.vtx_trackpid)
      event.reaction_vtx_trackpid_g4.emplace_back(pid.begin(), pid.end());
    event.reaction_vtx_px_g4 = reaction.vtx_px;
    event.reaction_vtx_py_g4 = reaction.vtx_py;
    event.reaction_vtx_pz_g4 = reaction.vtx_pz;
  }

  void FillBeamReconstruction(const GeantClusters& beam)
  {
    if (beam.pad.empty()) return;
    TPCAnalyzer analyzer;
    analyzer.ReCalcTPCHitsGeant4(beam.pad, beam.de, beam.x, beam.y, beam.z);
    analyzer.TrackSearchTPCHelix(); // unchanged common helix tracking, beam hits only

    TPCLocalTrackHelix* best = nullptr;
    for (Int_t it = 0; it < analyzer.GetNTracksTPCHelix(); ++it) {
      auto* track = analyzer.GetTrackTPCHelix(it);
      if (track && (!best || track->GetNHit() > best->GetNHit())) best = track;

    }
    if (!best) return;

    event.beam_reco_found = 1;
    event.beam_reco_nhit = best->GetNHit();
    event.beam_reco_charge = best->GetCharge();
    event.beam_reco_chisqr = best->GetChiSquare();
    event.beam_reco_cx = best->Getcx(); event.beam_reco_cy = best->Getcy();
    event.beam_reco_z0 = best->Getz0(); event.beam_reco_r = best->Getr(); event.beam_reco_dz = best->Getdz();
    Double_t par[5]; best->GetParam(par);
    for (Int_t ih = 0; ih < best->GetNHit(); ++ih) {
      auto* hit = best->GetHitInOrder(ih);
      if (!hit) continue;
      const TVector3& pos = hit->GetLocalHitPos();
      const Double_t theta = hit->GetTheta();
      const TVector3 cal = best->GetPosition(par, theta);
      const TVector3 mom = best->CalcHelixMom(par, theta);
      event.beam_reco_hit_x.push_back(pos.X()); event.beam_reco_hit_y.push_back(pos.Y()); event.beam_reco_hit_z.push_back(pos.Z());
      event.beam_reco_cal_x.push_back(cal.X()); event.beam_reco_cal_y.push_back(cal.Y()); event.beam_reco_cal_z.push_back(cal.Z());
      event.beam_reco_px.push_back(mom.X()); event.beam_reco_py.push_back(mom.Y()); event.beam_reco_pz.push_back(mom.Z());

    }
  }

  void FillTracks(TPCAnalyzer& analyzer)
  {
    event.ntTpc = analyzer.GetNTracksTPCHelix();
    for (Int_t it = 0; it < event.ntTpc; ++it) {
      auto* track = analyzer.GetTrackTPCHelix(it);
      if (!track) continue;
      event.nhtrack.push_back(track->GetNHit());
      event.is_beam.push_back(track->GetIsBeam());
      event.is_accidental.push_back(track->GetIsAccidental());
      event.charge.push_back(track->GetCharge());
      event.pid.push_back(track->GetPid());
      event.chisqr.push_back(track->GetChiSquare());
      event.helix_cx.push_back(track->Getcx()); event.helix_cy.push_back(track->Getcy());
      event.helix_z0.push_back(track->Getz0()); event.helix_r.push_back(track->Getr()); event.helix_dz.push_back(track->Getdz());
      if (track->IsThetaCalculated()) {
        event.helix_theta_min.push_back(track->GetMint());
        event.helix_theta_max.push_back(track->GetMaxt());
      } else {
        event.helix_theta_min.push_back(TMath::QuietNaN());
        event.helix_theta_max.push_back(TMath::QuietNaN());
      }
      const TVector3 mom = track->GetMom0();
      event.mom0_x.push_back(mom.X()); event.mom0_y.push_back(mom.Y()); event.mom0_z.push_back(mom.Z()); event.mom0.push_back(mom.Mag());
      event.dE.push_back(track->GetTrackdE()); event.dEdx.push_back(track->GetdEdx(0.8));
      event.hitlayer.emplace_back(); event.hitpos_x.emplace_back(); event.hitpos_y.emplace_back(); event.hitpos_z.emplace_back();
      event.calpos_x.emplace_back(); event.calpos_y.emplace_back(); event.calpos_z.emplace_back();
      event.residual.emplace_back(); event.residual_x.emplace_back(); event.residual_y.emplace_back(); event.residual_z.emplace_back();
      event.track_cluster_de.emplace_back(); event.track_cluster_size.emplace_back(); event.track_cluster_mrow.emplace_back();
      event.track_cluster_de_center.emplace_back(); event.track_cluster_x_center.emplace_back();
      event.track_cluster_y_center.emplace_back(); event.track_cluster_z_center.emplace_back(); event.track_cluster_row_center.emplace_back();
      for (Int_t ih = 0; ih < track->GetNHit(); ++ih) {
        auto* hit = track->GetHit(ih);
        if (!hit) continue;
        const auto& pos = hit->GetLocalHitPos();
        const auto& cal = hit->GetLocalCalPosHelix();
        const auto& res = hit->GetResidualVect();
        event.hitlayer.back().push_back(hit->GetLayer());
        event.hitpos_x.back().push_back(pos.X()); event.hitpos_y.back().push_back(pos.Y()); event.hitpos_z.back().push_back(pos.Z());
        event.calpos_x.back().push_back(cal.X()); event.calpos_y.back().push_back(cal.Y()); event.calpos_z.back().push_back(cal.Z());
        event.residual.back().push_back(hit->GetResidual());
        event.residual_x.back().push_back(res.X()); event.residual_y.back().push_back(res.Y()); event.residual_z.back().push_back(res.Z());
        auto* cluster_hit = hit->GetHit();
        auto* cluster = cluster_hit ? cluster_hit->GetParentCluster() : nullptr;
        auto* center_hit = cluster ? cluster->GetCenterHit() : nullptr;
        event.track_cluster_de.back().push_back(cluster ? cluster->GetDe() : TMath::QuietNaN());
        event.track_cluster_size.back().push_back(cluster ? cluster->GetClusterSize() : TMath::QuietNaN());
        event.track_cluster_mrow.back().push_back(cluster ? cluster->MeanRow() : TMath::QuietNaN());
        event.track_cluster_de_center.back().push_back(center_hit ? center_hit->GetCDe() : TMath::QuietNaN());
        event.track_cluster_x_center.back().push_back(center_hit ? center_hit->GetPosition().X() : TMath::QuietNaN());
        event.track_cluster_y_center.back().push_back(center_hit ? center_hit->GetPosition().Y() : TMath::QuietNaN());
        event.track_cluster_z_center.back().push_back(center_hit ? center_hit->GetPosition().Z() : TMath::QuietNaN());
        event.track_cluster_row_center.back().push_back(center_hit ? center_hit->GetRow() : TMath::QuietNaN());
      }

    }
  }

  void FillHelixPairKinematics(TPCAnalyzer& analyzer)
  {
    const Int_t ntracks = analyzer.GetNTracksTPCHelix();
    const Double_t qnan = TMath::QuietNaN();
    event.closeDistTpc.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.vtxTpc.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.vtyTpc.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.vtzTpc.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.mom_vtx.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.mom_vty.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    event.mom_vtz.assign(ntracks, std::vector<Double_t>(ntracks, qnan));
    for (Int_t i = 0; i < ntracks; ++i) {
      for (Int_t j = 0; j < ntracks; ++j) {
        if (i == j) continue;
        TPCVertex* vertex = analyzer.FindVertexTPC(i, j);
        if (!vertex || !vertex->IsCalculated()) continue;
        Int_t vertex_track_index = -1;
        for (Int_t iv = 0; iv < vertex->GetNTracks(); ++iv)
          if (vertex->GetTrackId(iv) == i) { vertex_track_index = iv; break; }
        if (vertex_track_index < 0) continue;
        const TVector3& pos = vertex->GetVertex();
        const TVector3 mom = vertex->GetTrackMom(vertex_track_index);
        event.closeDistTpc[i][j] = vertex->GetClosestDist();
        event.vtxTpc[i][j] = pos.X(); event.vtyTpc[i][j] = pos.Y(); event.vtzTpc[i][j] = pos.Z();
        event.mom_vtx[i][j] = mom.X(); event.mom_vty[i][j] = mom.Y(); event.mom_vtz[i][j] = mom.Z();
      }

    }
  }

  struct TruthVertexParticle
  {
    Int_t source, vtx_index, vtx_type, pid;
    Double_t x, y, z, px, py, pz;
  };

  void FillTruthVertexComparison(TPCAnalyzer& analyzer, const GeantClusters& beam,
                                 const GeantClusters& reaction)
  {
    const Int_t ntracks = analyzer.GetNTracksTPCHelix();
    const Double_t qnan = TMath::QuietNaN();
    event.truth_found.assign(ntracks, 0); event.truth_source.assign(ntracks, -1);
    event.truth_trackid.assign(ntracks, -1); event.truth_pid.assign(ntracks, 0);
    event.truth_vtx_index.assign(ntracks, -1); event.truth_vtx_type.assign(ntracks, -1);
    event.truth_match_nhit.assign(ntracks, 0); event.truth_match_purity.assign(ntracks, qnan);
    event.truth_vtx_x.assign(ntracks, qnan); event.truth_vtx_y.assign(ntracks, qnan); event.truth_vtx_z.assign(ntracks, qnan);
    event.truth_vtx_px.assign(ntracks, qnan); event.truth_vtx_py.assign(ntracks, qnan); event.truth_vtx_pz.assign(ntracks, qnan);
    event.reco_vtx_px.assign(ntracks, qnan); event.reco_vtx_py.assign(ntracks, qnan); event.reco_vtx_pz.assign(ntracks, qnan);
    event.delta_vtx_px.assign(ntracks, qnan); event.delta_vtx_py.assign(ntracks, qnan); event.delta_vtx_pz.assign(ntracks, qnan);
    event.delta_vtx_p.assign(ntracks, qnan); event.delta_vtx_pt.assign(ntracks, qnan);

    std::map<std::pair<Int_t, Int_t>, TruthVertexParticle> truth;
    const auto add_truth = [&](const GeantClusters& in) {
      for (std::size_t iv = 0; iv < in.vtx_trackid.size(); ++iv) {
        if (iv >= in.vtx_type.size() || iv >= in.vtx_x.size() || iv >= in.vtx_y.size() || iv >= in.vtx_z.size() ||
            iv >= in.vtx_trackpid.size() || iv >= in.vtx_px.size() || iv >= in.vtx_py.size() || iv >= in.vtx_pz.size()) continue;
        const auto nmember = std::min({in.vtx_trackid[iv].size(), in.vtx_trackpid[iv].size(), in.vtx_px[iv].size(), in.vtx_py[iv].size(), in.vtx_pz[iv].size()});
        for (std::size_t ip = 0; ip < nmember; ++ip) {
          const Int_t trackid = in.vtx_trackid[iv][ip];
          if (trackid <= 0) continue;
          truth[{in.generator, trackid}] = {in.generator, static_cast<Int_t>(iv), in.vtx_type[iv], in.vtx_trackpid[iv][ip],
                                             in.vtx_x[iv], in.vtx_y[iv], in.vtx_z[iv],
                                             in.vtx_px[iv][ip], in.vtx_py[iv][ip], in.vtx_pz[iv][ip]};
        }
      }
    };
    add_truth(beam); add_truth(reaction);

    for (Int_t it = 0; it < ntracks; ++it) {
      auto* track = analyzer.GetTrackTPCHelix(it);
      if (!track) continue;
      std::map<std::pair<Int_t, Int_t>, Int_t> votes;
      for (Int_t ih = 0; ih < track->GetNHit(); ++ih) {
        const auto* hit = track->GetHit(ih);
        if (!hit) continue;
        const TVector3& pos = hit->GetLocalHitPos();
        Double_t best_d2 = 1.e99; Int_t best = -1;
        for (std::size_t iraw = 0; iraw < event.raw_trackid.size(); ++iraw) {
          const Double_t dx = pos.X()-event.raw_hitpos_x[iraw];
          const Double_t dy = pos.Y()-event.raw_hitpos_y[iraw];
          const Double_t dz = pos.Z()-event.raw_hitpos_z[iraw];
          const Double_t d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best_d2) { best_d2 = d2; best = static_cast<Int_t>(iraw); }
        }
        if (best >= 0 && best_d2 < 1.e-8) ++votes[{event.raw_source[best], event.raw_trackid[best]}];
      }
      if (votes.empty()) continue;
      const auto best_vote = std::max_element(votes.begin(), votes.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });
      event.truth_match_nhit[it] = best_vote->second;
      event.truth_match_purity[it] = static_cast<Double_t>(best_vote->second)/track->GetNHit();
      const auto truth_it = truth.find(best_vote->first);
      if (truth_it == truth.end()) continue;
      const auto& t = truth_it->second;
      event.truth_found[it] = 1; event.truth_source[it] = t.source; event.truth_trackid[it] = best_vote->first.second;
      event.truth_pid[it] = t.pid; event.truth_vtx_index[it] = t.vtx_index; event.truth_vtx_type[it] = t.vtx_type;
      event.truth_vtx_x[it] = t.x; event.truth_vtx_y[it] = t.y; event.truth_vtx_z[it] = t.z;
      event.truth_vtx_px[it] = t.px; event.truth_vtx_py[it] = t.py; event.truth_vtx_pz[it] = t.pz;
    }

    for (Int_t it = 0; it < ntracks; ++it) {
      if (!event.truth_found[it]) continue;
      Int_t partner = -1;
      for (Int_t jt = 0; jt < ntracks; ++jt) {
        if (jt != it && event.truth_found[jt] && event.truth_source[jt] == event.truth_source[it] &&
            event.truth_vtx_index[jt] == event.truth_vtx_index[it]) { partner = jt; break; }
      }
      if (partner < 0 || !std::isfinite(event.mom_vtx[it][partner])) continue;
      const Double_t px = event.mom_vtx[it][partner], py = event.mom_vty[it][partner], pz = event.mom_vtz[it][partner];
      event.reco_vtx_px[it] = px; event.reco_vtx_py[it] = py; event.reco_vtx_pz[it] = pz;
      const Double_t dx = px-event.truth_vtx_px[it], dy = py-event.truth_vtx_py[it], dz = pz-event.truth_vtx_pz[it];
      event.delta_vtx_px[it] = dx; event.delta_vtx_py[it] = dy; event.delta_vtx_pz[it] = dz;
      event.delta_vtx_p[it] = std::sqrt(dx*dx + dy*dy + dz*dz);
      event.delta_vtx_pt[it] = std::sqrt(dx*dx + dz*dz);
    }
  }

  void ProcessEvent(const GeantClusters* beam, const GeantClusters* reaction,
                    TTree& tree)
  {
    static const GeantClusters empty;
    const auto& beam_data = beam ? *beam : empty;
    const auto& reaction_data = reaction ? *reaction : empty;
    const Int_t effective_evnum = reaction ? reaction->effective_evnum : beam->effective_evnum;

    event.clear();
    event.event_number = effective_evnum;
    event.beam_event_number = beam ? beam->evnum : -1;
    event.reaction_event_number = reaction ? reaction->evnum : -1;
    event.trig_flag = reaction ? reaction->trig_flag : beam_data.trig_flag;
    FillRaw(beam_data, reaction_data);
    FillReactionVertexTruth(reaction_data);
    FillInitialBeamTruth(beam_data);
    FillReactionBeamTruth(reaction_data);
    FillBeamTruth(beam_data);
    FillBeamReconstruction(beam_data);
    event.status = 1;
    if (event.nhTpc > 0) {
      TPCAnalyzer analyzer;
      analyzer.ReCalcTPCHitsGeant4(event.raw_padid, event.raw_de,
                                   event.raw_hitpos_x, event.raw_hitpos_y, event.raw_hitpos_z);
      event.status = 2;
      analyzer.TrackSearchTPCHelix(); // unchanged common helix tracking algorithm
      event.status = 3;
      FillTracks(analyzer);
      FillHelixPairKinematics(analyzer);
      FillTruthVertexComparison(analyzer, beam_data, reaction_data);

    }
    tree.Fill();
  }
}

namespace dst
{
  enum kArgc
  {
    kProcess, kConfFile, kGeant4, kOutFile, nArgc
  };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[Geant4]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "g4hyptpc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
}

Int_t main(Int_t argc, char** argv)
{
  std::vector<std::string> arg(argv, argv+argc);
  if (!dst::CheckArg(arg)) return EXIT_FAILURE;
  if (!gConf.Initialize(arg[dst::kConfFile]) || !gConf.InitializeHistograms()
      || !gConf.InitializeUnpacker()) return EXIT_FAILURE;

  TChain input(dst::TreeName[dst::kGeant4]);
  if (input.Add(arg[dst::kGeant4].c_str()) == 0 || !SetupInput(input)) return EXIT_FAILURE;
  TFile output(arg[dst::kOutFile].c_str(), "recreate");
  TTree tree("tpc", "Geant4 beam+reaction TPC helix tracking");
  BookTree(tree);

  CatchSignal::Set();
  GeantClusters pending_beam;
  const Long64_t nentries = input.GetEntries();
  Long64_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Long64_t stop = nentries;
  const Long64_t max_loop = gUnpacker.get_max_loop();
  if (max_loop > 0) stop = std::min(nentries, skip + max_loop);
  // A job may begin at a 7202 entry. Read its immediately preceding 7201
  // only as pairing context; output remains restricted to [skip, stop).
  const Long64_t first = skip > 0 ? skip - 1 : 0;
  for (Long64_t ie = first; ie < stop && !CatchSignal::Stop(); ++ie) {
    gCounter.check();
    input.GetEntry(ie);
    if (src.generator == kBeamGenerator) {
      GeantClusters beam;
      beam.evnum = src.evnum;
      beam.effective_evnum = src.effective_evnum;
      beam.trig_flag = src.trig_flag;
      beam.generator = src.generator;
      beam.copyInitialBeam(src);
      beam.append(src, false); // retain the full beam for beam momentum reconstruction
      SetVertices(beam, src);
      if (kGeant4InputMode == Geant4InputMode::kBeamOnly) {
        if (ie >= skip) ProcessEvent(&beam, nullptr, tree);
      } else {
        pending_beam = std::move(beam);
      }
      continue;

    }
    // The second generator depends on the reaction channel (e.g. 7202, 7205).
    // Every non-7201 entry paired with the preceding beam is a reaction entry.
    if (ie < skip) continue;

    GeantClusters reaction;
    reaction.evnum = src.evnum;
    reaction.effective_evnum = src.effective_evnum;
    reaction.trig_flag = src.trig_flag;
    reaction.generator = src.generator;
    reaction.copyInitialBeam(src);
    reaction.append(src, false); // preserve every 7202 cluster
    SetVertices(reaction, src);
    if (kGeant4InputMode == Geant4InputMode::kReactionOnly) {
      ProcessEvent(nullptr, &reaction, tree);
    } else if (kGeant4InputMode == Geant4InputMode::kBeamAndReaction
               && pending_beam.effective_evnum == src.effective_evnum) {
      ProcessEvent(&pending_beam, &reaction, tree);
      pending_beam.clear();

    }
  }
  output.Write();
  output.Close();
  return EXIT_SUCCESS;
}

Bool_t ConfMan::InitializeParameterFiles()
{
  return InitializeParameter<TPCParamMan>("TPCPRM")
    && InitializeParameter<TPCPositionCorrector>("TPCPOS")
    && InitializeParameter<UserParamMan>("USER");
}

Bool_t ConfMan::InitializeHistograms()
{
  TPCEventAnalyzer::SetDstCalibFlag(false);
  hist::BuildStatus();
  hist::BuildTPCHelixTracking(false);
  return true;
}

Bool_t ConfMan::FinalizeProcess()
{
  return true;
}
