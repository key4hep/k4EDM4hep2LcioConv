#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoCaloParticleAssociationCollection.h"
#include "edm4hep/MCRecoClusterParticleAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociationCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/RecoParticleVertexAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TPCHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/VertexCollection.h"

#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ParticleID.h>
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TPCHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/Vertex.h>
#include <lcio.h>

// bool compare(const EVENT::CaloHitContribution *lcio,
//             const edm4hep::CaloHitContribution &edm4hep);
// bool compare(const lcio::LCCollection *lcioCollection,
//             const edm4hep::CaloHitContributionCollection &edm4hepCollection);

bool compare(const EVENT::CalorimeterHit* lcio, const edm4hep::CalorimeterHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::CalorimeterHitCollection& edm4hepCollection);

bool compare(const EVENT::Cluster* lcio, const edm4hep::Cluster& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::ClusterCollection& edm4hepCollection);

// bool compare(const EVENT::EventHeader *lcio,
//              const edm4hep::EventHeader &edm4hep);
// bool compare(const lcio::LCCollection *lcioCollection,
//              const edm4hep::EventHeaderCollection &edm4hepCollection);

bool compare(const EVENT::MCParticle* lcio, const edm4hep::MCParticle& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::MCParticleCollection& edm4hepCollection);

bool compare(const EVENT::RawCalorimeterHit* lcio, const edm4hep::RawCalorimeterHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::RawCalorimeterHitCollection& edm4hepCollection);

bool compare(const EVENT::ReconstructedParticle* lcio, const edm4hep::ReconstructedParticle& edm4hep);
bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::ReconstructedParticleCollection& edm4hepCollection);

bool compare(const EVENT::SimCalorimeterHit* lcio, const edm4hep::SimCalorimeterHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimCalorimeterHitCollection& edm4hepCollection);

bool compare(const EVENT::SimTrackerHit* lcio, const edm4hep::SimTrackerHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimTrackerHitCollection& edm4hepCollection);

bool compare(const EVENT::TPCHit* lcio, const edm4hep::TPCHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TPCHitCollection& edm4hepCollection);

bool compare(const EVENT::TrackerHit* lcio, const edm4hep::TrackerHit& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHitCollection& edm4hepCollection);

bool compare(const EVENT::TrackerHitPlane* lcio, const edm4hep::TrackerHitPlane& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHitPlaneCollection& edm4hepCollection);

bool compare(const EVENT::Track* lcio, const edm4hep::Track& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackCollection& edm4hepCollection);

bool compare(const EVENT::Vertex* lcio, const edm4hep::Vertex& edm4hep);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::VertexCollection& edm4hepCollection);

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
