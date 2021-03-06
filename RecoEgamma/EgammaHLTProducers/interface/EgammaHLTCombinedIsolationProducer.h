// -*- C++ -*-
//
// Package:    EgammaHLTProducers
// Class:      EgammaHLTCombinedIsolationProducer
//
/**\class EgammaHLTCombinedIsolationProducer EgammaHLTCombinedIsolationProducer.cc RecoEgamma/EgammaHLTProducers/interface/EgammaHLTCombinedIsolationProducer.h
*/
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

namespace edm {
  class ConfigurationDescriptions;
}

class EgammaHLTCombinedIsolationProducer : public edm::global::EDProducer<> {
public:
  explicit EgammaHLTCombinedIsolationProducer(const edm::ParameterSet&);
  ~EgammaHLTCombinedIsolationProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  // ----------member data ---------------------------

  edm::EDGetTokenT<reco::RecoEcalCandidateCollection> recoEcalCandidateProducer_;
  std::vector<edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> > IsolTag_;
  std::vector<double> IsolWeight_;
  const edm::ParameterSet conf_;
};
