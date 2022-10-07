// system include files
#include <memory>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class TriggerObjectTableProducer : public edm::stream::EDProducer<> {
public:
  explicit TriggerObjectTableProducer(const edm::ParameterSet &iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        src_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("src"))),
        l1EG_(consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("l1EG"))),
        l1Sum_(consumes<l1t::EtSumBxCollection>(iConfig.getParameter<edm::InputTag>("l1Sum"))),
        l1Jet_(consumes<l1t::JetBxCollection>(iConfig.getParameter<edm::InputTag>("l1Jet"))),
        l1Muon_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1Muon"))),
        l1Tau_(consumes<l1t::TauBxCollection>(iConfig.getParameter<edm::InputTag>("l1Tau"))) {
    std::vector<edm::ParameterSet> selPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("selections");
    sels_.reserve(selPSets.size());
    std::stringstream idstr, qualitystr;
    idstr << "ID of the object: ";
    for (auto &pset : selPSets) {
      sels_.emplace_back(pset);
      idstr << sels_.back().id << " = " << sels_.back().name;
      if (sels_.size() < selPSets.size())
        idstr << ", ";
      if (!sels_.back().qualityBitsDoc.empty()) {
        qualitystr << sels_.back().qualityBitsDoc << " for " << sels_.back().name << "; ";
      }
    }
    idDoc_ = idstr.str();
    bitsDoc_ = qualitystr.str();

    produces<nanoaod::FlatTable>();
  }

  ~TriggerObjectTableProducer() override {}

private:
  void produce(edm::Event &, edm::EventSetup const &) override;

  std::string name_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> src_;
  std::string idDoc_, bitsDoc_;

  edm::EDGetTokenT<l1t::EGammaBxCollection> l1EG_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> l1Sum_;
  edm::EDGetTokenT<l1t::JetBxCollection> l1Jet_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1Muon_;
  edm::EDGetTokenT<l1t::TauBxCollection> l1Tau_;

  struct SelectedObject {
    std::string name;
    int id;
    StringCutObjectSelector<pat::TriggerObjectStandAlone> cut;
    StringCutObjectSelector<pat::TriggerObjectStandAlone> l1cut, l1cut_2, l2cut;
    float l1DR2, l1DR2_2, l2DR2;
    bool skipObjectsNotPassingQualityBits;
    StringObjectFunction<pat::TriggerObjectStandAlone> qualityBits;
    std::string qualityBitsDoc;

    SelectedObject(const edm::ParameterSet &pset)
        : name(pset.getParameter<std::string>("name")),
          id(pset.getParameter<int>("id")),
          cut(pset.getParameter<std::string>("sel")),
          l1cut(""),
          l1cut_2(""),
          l2cut(""),
          l1DR2(-1),
          l1DR2_2(-1),
          l2DR2(-1),
          skipObjectsNotPassingQualityBits(pset.getParameter<bool>("skipObjectsNotPassingQualityBits")),
          qualityBits(pset.getParameter<std::string>("qualityBits")),
          qualityBitsDoc(pset.getParameter<std::string>("qualityBitsDoc")) {
      if (pset.existsAs<std::string>("l1seed")) {
        l1cut = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l1seed"));
        l1DR2 = std::pow(pset.getParameter<double>("l1deltaR"), 2);
      }
      if (pset.existsAs<std::string>("l1seed_2")) {
        l1cut_2 = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l1seed_2"));
        l1DR2_2 = std::pow(pset.getParameter<double>("l1deltaR_2"), 2);
      }
      if (pset.existsAs<std::string>("l2seed")) {
        l2cut = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l2seed"));
        l2DR2 = std::pow(pset.getParameter<double>("l2deltaR"), 2);
      }
    }

    bool match(const pat::TriggerObjectStandAlone &obj) const { return cut(obj); }
  };

  std::vector<SelectedObject> sels_;
};

// ------------ method called to produce the data  ------------
void TriggerObjectTableProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  const auto& trigObjs = iEvent.get(src_);

  std::vector<std::pair<const pat::TriggerObjectStandAlone *, const SelectedObject *>> selected;
  for (const auto &obj : trigObjs) {
    for (const auto &sel : sels_) {
      if (sel.match(obj) && (sel.skipObjectsNotPassingQualityBits ? (int(sel.qualityBits(obj)) > 0) : true)) {
        selected.emplace_back(&obj, &sel);
        break;
      }
    }
  }

  // Self-cleaning
  std::map<const pat::TriggerObjectStandAlone *, int> selected_bits;
  for (unsigned int i = 0; i < selected.size(); ++i) {
    const auto &obj = *selected[i].first;
    const auto &sel = *selected[i].second;
    selected_bits[&obj] = int(sel.qualityBits(obj));

    for (unsigned int j = 0; j < i; ++j) {
      const auto &obj2 = *selected[j].first;
      const auto &sel2 = *selected[j].second;
      if (sel.id == sel2.id && abs(obj.pt() - obj2.pt()) < 1e-6 && deltaR2(obj, obj2) < 1e-6) {
        selected_bits[&obj2] |= selected_bits[&obj];  //Keep filters from all the objects
        selected.erase(selected.begin() + i);
        i--;
      }
    }
  }

  const auto& l1EG = iEvent.get(l1EG_);
  const auto& l1Sum = iEvent.get(l1Sum_);
  const auto& l1Jet = iEvent.get(l1Jet_);
  const auto& l1Muon = iEvent.get(l1Muon_);
  const auto& l1Tau = iEvent.get(l1Tau_);

  std::vector<pair<pat::TriggerObjectStandAlone, int>> l1Objects;
  l1Objects.reserve(l1EG.size() + l1Sum.size() + l1Jet.size() + l1Muon.size() + l1Tau.size());

  for (const auto& obj: l1EG) {
    pat::TriggerObjectStandAlone l1obj(obj.p4());
    l1obj.setCollection("L1EG");
    l1obj.addTriggerObjectType(trigger::TriggerL1EG);
    l1Objects.emplace_back(l1obj, obj.hwIso());
  }

  for (const auto& obj: l1Sum) {
    pat::TriggerObjectStandAlone l1obj(obj.p4());

    switch (obj.getType()) {
      case l1t::EtSum::EtSumType::kMissingEt:
        l1obj.addTriggerObjectType(trigger::TriggerL1ETM);
        l1obj.setCollection("L1ETM");
        break;

      case l1t::EtSum::EtSumType::kMissingEtHF:
        l1obj.addTriggerObjectType(trigger::TriggerL1ETM);
        l1obj.setCollection("L1ETMHF");
        break;

      case l1t::EtSum::EtSumType::kTotalEt:
        l1obj.addTriggerObjectType(trigger::TriggerL1ETT);
        l1obj.setCollection("L1ETT");
        break;

      case l1t::EtSum::EtSumType::kTotalEtEm:
        l1obj.addTriggerObjectType(trigger::TriggerL1ETT);
        l1obj.setCollection("L1ETEm");
        break;

      case l1t::EtSum::EtSumType::kTotalHt:
        l1obj.addTriggerObjectType(trigger::TriggerL1HTT);
        l1obj.setCollection("L1HTT");
        break;

      case l1t::EtSum::EtSumType::kTotalHtHF:
        l1obj.addTriggerObjectType(trigger::TriggerL1HTT);
        l1obj.setCollection("L1HTTHF");
        break;

      case l1t::EtSum::EtSumType::kMissingHt:
        l1obj.addTriggerObjectType(trigger::TriggerL1HTM);
        l1obj.setCollection("L1HTM");
        break;

      case l1t::EtSum::EtSumType::kMissingHtHF:
        l1obj.addTriggerObjectType(trigger::TriggerL1HTM);
        l1obj.setCollection("L1HTMHF");
        break;

      default:
        continue;
    }

    l1Objects.emplace_back(l1obj, obj.hwIso());
  }

  for (const auto& obj: l1Jet) {
    pat::TriggerObjectStandAlone l1obj(obj.p4());
    l1obj.setCollection("L1Jet");
    l1obj.addTriggerObjectType(trigger::TriggerL1Jet);
    l1Objects.emplace_back(l1obj, obj.hwIso());
  }

  for (const auto& obj: l1Muon) {
    pat::TriggerObjectStandAlone l1obj(obj.p4());
    l1obj.setCollection("L1Mu");
    l1obj.addTriggerObjectType(trigger::TriggerL1Mu);
    l1obj.setCharge(obj.charge());
    l1Objects.emplace_back(l1obj, obj.hwIso());
  }

  for (const auto& obj: l1Tau) {
    pat::TriggerObjectStandAlone l1obj(obj.p4());
    l1obj.setCollection("L1Tau");
    l1obj.addTriggerObjectType(trigger::TriggerL1Tau);
    l1Objects.emplace_back(l1obj, obj.hwIso());
  }

  unsigned int nobj = selected.size();
  std::vector<float> pt(nobj, 0), eta(nobj, 0), phi(nobj, 0), l1pt(nobj, 0), l1pt_2(nobj, 0), l2pt(nobj, 0);
  std::vector<int> id(nobj, 0), bits(nobj, 0), l1iso(nobj, 0), l1charge(nobj, 0);
  for (unsigned int i = 0; i < nobj; ++i) {
    const auto &obj = *selected[i].first;
    const auto &sel = *selected[i].second;
    pt[i] = obj.pt();
    eta[i] = obj.eta();
    phi[i] = obj.phi();
    id[i] = sel.id;
    bits[i] = selected_bits[&obj];
    if (sel.l1DR2 > 0) {
      float best = sel.l1DR2;
      for (const auto &l1obj : l1Objects) {
        const auto &seed = l1obj.first;
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l1cut(seed)) {
          l1pt[i] = seed.pt();
          l1iso[i] = l1obj.second;
          l1charge[i] = seed.charge();
        }
      }
    }
    if (sel.l1DR2_2 > 0) {
      float best = sel.l1DR2_2;
      for (const auto &l1obj : l1Objects) {
        const auto &seed = l1obj.first;
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l1cut_2(seed)) {
          l1pt_2[i] = seed.pt();
        }
      }
    }
    if (sel.l2DR2 > 0) {
      float best = sel.l2DR2;
      for (const auto &seed : trigObjs) {
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l2cut(seed)) {
          l2pt[i] = seed.pt();
        }
      }
    }
  }

  auto tab = std::make_unique<nanoaod::FlatTable>(nobj, name_, false, false);
  tab->addColumn<int>("id", id, idDoc_);
  tab->addColumn<float>("pt", pt, "pt", 12);
  tab->addColumn<float>("eta", eta, "eta", 12);
  tab->addColumn<float>("phi", phi, "phi", 12);
  tab->addColumn<float>("l1pt", l1pt, "pt of associated L1 seed", 8);
  tab->addColumn<int>("l1iso", l1iso, "iso of associated L1 seed");
  tab->addColumn<int>("l1charge", l1charge, "charge of associated L1 seed");
  tab->addColumn<float>("l1pt_2", l1pt_2, "pt of associated secondary L1 seed", 8);
  tab->addColumn<float>("l2pt", l2pt, "pt of associated 'L2' seed (i.e. HLT before tracking/PF)", 10);
  tab->addColumn<int>("filterBits", bits, "extra bits of associated information: " + bitsDoc_);
  iEvent.put(std::move(tab));
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerObjectTableProducer);
