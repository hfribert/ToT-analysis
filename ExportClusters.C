///   clusters  -- one entry per cluster
///     det             0 = TRK, 1 = FT3
///     layer           Trk: layer 0..10 (VD 0-2, ML/OT 3-10); FT3: disk 0..11
///     chipID          O2 sensor index
///     row, col        pixel address of the seed (highest-charge) pixel on the chip
///     qSeed           collected charge of the seed pixel [e-]
///     x, y, z         global position of the cluster centre [cm]
///     event, trackID  MC label (-1 for noise clusters)
///     pdg             PDG code of the MC particle (0 for noise clusters)
///
///   particles -- every MC particle of the job
///     event, trackID, pdg, primary, mother, vx, vy, vz [cm], px, py, pz [GeV/c]
///
/// Usage (in the O2 environment):
///   root -l -b -q 'ExportClusters.C+("<jobdir>", "<out.root>")'
///   root -l -b -q 'ExportClusters.C+("<jobdir>", "<out.root>", true)'   // also noise clusters

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsTRKFT3/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "FT3Base/GeometryTGeo.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKReconstruction/Clusterer.h"
#endif

R__LOAD_LIBRARY(libO2TRKReconstruction)

namespace
{
struct Row {
  int det, layer, chipID, row, col;
  float qSeed, x, y, z;
  int event, trackID, pdg;
};

template <int DetID, typename Geo>
void exportDetector(const std::string& clusFile, Geo* geo, const std::vector<std::vector<int>>& pdgOf,
                    TTree* out, Row& r, bool writeNoise, long& nSignal, long& nNoise)
{
  constexpr bool isTRK = DetID == o2::detectors::DetID::TRK;
  const std::string det = isTRK ? "TRK" : "FT3";
  std::unique_ptr<TFile> f(TFile::Open(clusFile.c_str()));
  if (!f || f->IsZombie()) {
    LOGP(warning, "Cannot open {}, skipping {}", clusFile, det);
    return;
  }
  auto* t = (TTree*)f->Get("o2sim");

  int nLayers = 0;
  while (t->GetBranch(Form("%sClusterComp_%d", det.c_str(), nLayers))) {
    ++nLayers;
  }
  std::vector<std::vector<o2::trkft3::Cluster<DetID>>*> clus(nLayers, nullptr);
  std::vector<std::vector<unsigned char>*> patt(nLayers, nullptr);
  std::vector<std::vector<uint32_t>*> pixQ(nLayers, nullptr);
  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> labels(nLayers, nullptr);
  for (int l = 0; l < nLayers; ++l) {
    for (const auto* b : {"ClusterPatt", "ClusterPixCharge", "ClusterMCTruth"}) {
      if (!t->GetBranch(Form("%s%s_%d", det.c_str(), b, l))) {
        LOGP(fatal, "Missing branch {}{}_{} in {}", det, b, l, clusFile);
      }
    }
    t->SetBranchAddress(Form("%sClusterComp_%d", det.c_str(), l), &clus[l]);
    t->SetBranchAddress(Form("%sClusterPatt_%d", det.c_str(), l), &patt[l]);
    t->SetBranchAddress(Form("%sClusterPixCharge_%d", det.c_str(), l), &pixQ[l]);
    t->SetBranchAddress(Form("%sClusterMCTruth_%d", det.c_str(), l), &labels[l]);
  }

  for (Long64_t entry = 0; entry < t->GetEntries(); ++entry) {
    t->GetEntry(entry);
    for (int l = 0; l < nLayers; ++l) {
      const auto& cls = *clus[l];
      const unsigned char* p = patt[l]->data();
      size_t pixOff = 0;
      for (size_t ic = 0; ic < cls.size(); ++ic) {
        const auto& cl = cls[ic];
        const unsigned char* pattern = p;
        const int rowSpan = p[0], colSpan = p[1];
        p += 2 + (rowSpan * colSpan + 7) / 8;
        const size_t pixBeg = pixOff;
        pixOff += cl.size;

        const auto& labs = labels[l]->getLabels(ic);
        const bool noise = labs.empty() || !labs[0].isValid() || labs[0].isNoise();
        if (noise && !writeNoise) {
          continue;
        }

        int kSeed = 0;
        uint32_t qSeed = 0;
        for (int k = 0; k < cl.size && pixBeg + k < pixQ[l]->size(); ++k) {
          if ((*pixQ[l])[pixBeg + k] > qSeed) {
            qSeed = (*pixQ[l])[pixBeg + k];
            kSeed = k;
          }
        }
        int seedBit = 0;
        for (int bit = 0, nSet = 0; bit < rowSpan * colSpan; ++bit) {
          if (pattern[2 + bit / 8] & (0x80 >> (bit % 8))) {
            if (nSet++ == kSeed) {
              seedBit = bit;
              break;
            }
          }
        }

        const auto loc = o2::trk::Clusterer<DetID>::getClusterLocalCoordinates(cl, pattern, 0.f);
        const auto glo = geo->getMatrixL2G(cl.chipID) * loc;

        r.det = isTRK ? 0 : 1;
        if constexpr (isTRK) {
          r.layer = geo->getLayerTRK(cl.chipID);
        } else {
          r.layer = geo->getLayer(cl.chipID);
        }
        r.chipID = cl.chipID;
        r.row = cl.row + seedBit / colSpan;
        r.col = cl.col + seedBit % colSpan;
        r.qSeed = qSeed;
        r.x = glo.X();
        r.y = glo.Y();
        r.z = glo.Z();
        r.event = r.trackID = -1;
        r.pdg = 0;
        if (!noise) {
          r.event = labs[0].getEventID();
          r.trackID = labs[0].getTrackID();
          if (r.event < (int)pdgOf.size() && r.trackID >= 0 && r.trackID < (int)pdgOf[r.event].size()) {
            r.pdg = pdgOf[r.event][r.trackID];
          }
          ++nSignal;
        } else {
          ++nNoise;
        }
        out->Fill();
      }
    }
  }
  t->ResetBranchAddresses();
}
} // namespace

void ExportClusters(const std::string& jobDir, const std::string& outFile, bool writeNoise = false)
{
  // same geometry configuration as the simulation of this job
  const std::string ini = jobDir + "/o2sim_configuration.ini";
  if (!gSystem->AccessPathName(ini.c_str())) {
    o2::conf::ConfigurableParam::updateFromFile(ini, "TRKBase,FT3Base");
  }
  o2::base::GeometryManager::loadGeometry(jobDir + "/o2sim");
  const auto l2g = o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G);
  auto* geoTRK = o2::trk::GeometryTGeo::Instance();
  geoTRK->fillMatrixCache(l2g);
  const bool haveFT3 = !gSystem->AccessPathName((jobDir + "/o2clus_ft3.root").c_str());
  o2::ft3::GeometryTGeo* geoFT3 = nullptr;
  if (haveFT3) {
    geoFT3 = o2::ft3::GeometryTGeo::Instance();
    geoFT3->fillMatrixCache(l2g);
  }

  std::unique_ptr<TFile> fOut(TFile::Open(outFile.c_str(), "RECREATE"));
  auto* tPart = new TTree("particles", "MC particles");
  int pEvent, pTrack, pPdg, pMother;
  bool pPrimary;
  float pVx, pVy, pVz, pPx, pPy, pPz;
  tPart->Branch("event", &pEvent, "event/I");
  tPart->Branch("trackID", &pTrack, "trackID/I");
  tPart->Branch("pdg", &pPdg, "pdg/I");
  tPart->Branch("primary", &pPrimary, "primary/O");
  tPart->Branch("mother", &pMother, "mother/I");
  tPart->Branch("vx", &pVx, "vx/F");
  tPart->Branch("vy", &pVy, "vy/F");
  tPart->Branch("vz", &pVz, "vz/F");
  tPart->Branch("px", &pPx, "px/F");
  tPart->Branch("py", &pPy, "py/F");
  tPart->Branch("pz", &pPz, "pz/F");

  std::vector<std::vector<int>> pdgOf;
  {
    std::unique_ptr<TFile> fKine(TFile::Open((jobDir + "/o2sim_Kine.root").c_str()));
    if (!fKine || fKine->IsZombie()) {
      LOGP(fatal, "Cannot open {}/o2sim_Kine.root", jobDir);
    }
    auto* tKine = (TTree*)fKine->Get("o2sim");
    std::vector<o2::MCTrack>* tracks = nullptr;
    tKine->SetBranchAddress("MCTrack", &tracks);
    pdgOf.resize(tKine->GetEntries());
    for (Long64_t ev = 0; ev < tKine->GetEntries(); ++ev) {
      tKine->GetEntry(ev);
      for (size_t i = 0; i < tracks->size(); ++i) {
        const auto& trk = (*tracks)[i];
        pdgOf[ev].push_back(trk.GetPdgCode());
        pEvent = ev;
        pTrack = i;
        pPdg = trk.GetPdgCode();
        pPrimary = trk.isPrimary();
        pMother = trk.getMotherTrackId();
        pVx = trk.Vx();
        pVy = trk.Vy();
        pVz = trk.Vz();
        pPx = trk.Px();
        pPy = trk.Py();
        pPz = trk.Pz();
        tPart->Fill();
      }
    }
    tKine->ResetBranchAddresses();
    delete tracks;
  }

  fOut->cd();
  auto* tClus = new TTree("clusters", "TRK and FT3 clusters");
  Row r{};
  tClus->Branch("det", &r.det, "det/I");
  tClus->Branch("layer", &r.layer, "layer/I");
  tClus->Branch("chipID", &r.chipID, "chipID/I");
  tClus->Branch("row", &r.row, "row/I");
  tClus->Branch("col", &r.col, "col/I");
  tClus->Branch("qSeed", &r.qSeed, "qSeed/F");
  tClus->Branch("x", &r.x, "x/F");
  tClus->Branch("y", &r.y, "y/F");
  tClus->Branch("z", &r.z, "z/F");
  tClus->Branch("event", &r.event, "event/I");
  tClus->Branch("trackID", &r.trackID, "trackID/I");
  tClus->Branch("pdg", &r.pdg, "pdg/I");

  long nSignal = 0, nNoise = 0;
  exportDetector<o2::detectors::DetID::TRK>(jobDir + "/o2clus_trk.root", geoTRK, pdgOf, tClus, r, writeNoise, nSignal, nNoise);
  if (haveFT3) {
    exportDetector<o2::detectors::DetID::FT3>(jobDir + "/o2clus_ft3.root", geoFT3, pdgOf, tClus, r, writeNoise, nSignal, nNoise);
  }

  // simulation settings of the job for reference
  fOut->cd();
  for (const auto* cfg : {"sim_config.json", "sim+digi_config.json"}) {
    std::ifstream in(jobDir + "/" + cfg);
    if (in) {
      std::stringstream ss;
      ss << in.rdbuf();
      TNamed(cfg, ss.str().c_str()).Write();
    }
  }
  tPart->Write();
  tClus->Write();
  LOGP(info, "{}: {} particles, {} signal and {} noise clusters written to {}",
       jobDir, tPart->GetEntries(), nSignal, nNoise, outFile);
}
