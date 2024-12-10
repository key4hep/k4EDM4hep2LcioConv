#ifndef K4EDM4HEP2LCIOCONV_TEST_EDM4HEP2LCIOUTILITIES_H
#define K4EDM4HEP2LCIOCONV_TEST_EDM4HEP2LCIOUTILITIES_H

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include <edm4hep/CaloHitMCParticleLinkCollection.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/ClusterMCParticleLinkCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/ParticleIDCollection.h>
#include <edm4hep/RawTimeSeriesCollection.h>
#include <edm4hep/RecDqdxCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/VertexCollection.h>
#include <edm4hep/VertexRecoParticleLinkCollection.h>

#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

namespace podio {
class Frame;
}

namespace test_config {
constexpr static int nMCParticles = 5; ///< The number of MCParticles to generate

using IdxPair = std::pair<int, int>;
/// How to create the MC particle hierarchy, e.g. {4, 0} means that mc[4] will
/// have mc[0] as a parent, and mc[0] will get mc[4] as a daughter
const static std::vector<IdxPair> mcpParentIdcs = {{4, 0}, {4, 1}, {3, 2}, {3, 0}, {3, 1}, {2, 1}};

constexpr static int nCaloHits = 2;    ///< The number of CalorimeterHits to generate
constexpr static int nRawCaloHits = 2; ///< The number of RawCalorimeterHits to generate

constexpr static int nTPCHits = 4;     ///< The number of TPCHits (RawTimeSeries) to generate
constexpr static int nTPCRawWords = 6; ///< The number of raw words to put into each TPCHit

constexpr static int nTrackerHits = 5; ///< The number of TrackerHits to generate

constexpr static int nTracks = 4;                ///< The number of tracks to generate
constexpr static int nSubdetectorHitNumbers = 4; ///< The number of subdetector hits for each track
/// The tracker hits that should be added to each track
const static std::vector<std::size_t> trackTrackerHitIdcs = {0, 2, 4};
constexpr static int nTrackStates = 5; ///< The number of track states for each track
/// The tracks that should be linked, first index is the track to which the
/// second index will be added
const static std::vector<IdxPair> trackTrackIdcs = {{0, 2}, {1, 3}, {2, 3}, {3, 2}, {3, 0}};

constexpr static int nSimCaloHits = 3;          ///< The number of SimCalorimeterHits
constexpr static int nCaloHitContributions = 4; ///< The number of CalorimeterHit Contributions
/// idcs to setup relations between SimCalorimeterHits, CaloHitContributions
/// and MCParticles (in this order)
using CaloContIdx = std::tuple<int, int, int>;
/// The index values to use for setting up the relations
const static std::vector<CaloContIdx> simCaloHitMCIdcs = {{0, 0, 0}, {0, 1, 2}, {0, 2, 1}, {0, 3, 4},
                                                          {1, 0, 1}, {1, 1, 3}, {1, 2, 4}, {1, 3, 4},
                                                          {2, 0, 0}, {2, 1, 3}, {2, 2, 2}, {2, 3, 0}};

/// The number of clusters to create
constexpr static int nClusters = 5;
/// The number of subdetector energy entries to create
constexpr static int nSubdetectorEnergies = 6;
/// The calorimeter hits that should be associated with each cluster. First
/// index is the cluster, second is the calorimeter hit
const static std::vector<IdxPair> clusterHitIdcs = {{0, 0}, {0, 1}, {1, 0}, {2, 1}, {2, 0}, {3, 0}, {3, 0}};
/// The clustes (from inside the same collection) that should be added to each
/// cluster. First index is the cluster to which the second index cluster will
/// be added
const static std::vector<IdxPair> clusterClusterIdcs = {{0, 4}, {0, 3}, {0, 1}, {4, 3}, {4, 2}, {2, 3}, {1, 1}};

/// The number of reco particles to create
constexpr static int nRecoParticles = 6;
/// The related clusters that should be associated to the reco particles.
/// First index is the reco particle, second is the index of the cluster in
/// its collection
const static std::vector<IdxPair> recoClusterIdcs = {{0, 4}, {0, 3}, {2, 2}, {5, 1}, {5, 0}};

/// The related tracks that should be associated to the reco particles. First
/// index is the reco particle, second is the index of the track in its
/// collection
const static std::vector<IdxPair> recoTrackIdcs = {{0, 3}, {2, 2}, {5, 1}, {5, 0}};

/// The related reco partiles that should be associated to the reco particles.
/// First index is the reco particle, second is the index of the related reco
/// particle
const static std::vector<IdxPair> recoRecoIdcs = {{0, 3}, {2, 2}, {5, 1}, {5, 0}, {1, 2}, {2, 4}};

/// The number of entries for the generated ParticleID collections
const static std::vector<std::vector<int>> pidRecoIdcs = {{1, 3, 4}, {2, 3}, {0, 1, 2, 3, 4, 5}};

/// The number of vertices to create
constexpr static int nVertices = 3;
/// The particles that should be associated to each vertex. First index is the
/// vertex, second index is the reconstructed particle.
const static std::vector<IdxPair> vtxParticleIdcs = {{0, 1}, {0, 2}, {1, 1}, {1, 3}, {1, 4}, {2, 1}, {2, 3}, {2, 5}};
/// The (high level) reconstructed particles that will have a vertex assigned as
/// their decay vertex. First index is the reconstruced particle, second one is
/// the assigned vertex
const static std::vector<IdxPair> recoVtxIdcs = {{1, 0}, {0, 1}, {2, 2}};
} // namespace test_config

/**
 * Create an MCParticle collection with num_elements and create a parent
 * hierarchy using the passed mcp_parents_idx indices. The .first index will be
 * used to determine the daughter element, while the .second index will be used
 * to determine the parent.
 */
edm4hep::MCParticleCollection createMCParticles(const int num_elements,
                                                const std::vector<test_config::IdxPair>& mcp_parents_idx);

/**
 * Create a CalorimeterHit collection
 */
edm4hep::CalorimeterHitCollection createCalorimeterHits(const int num_elements);

/**
 * Create a RawCalorimeterHit collection
 */
edm4hep::RawCalorimeterHitCollection createRawCalorimeterHits(const int num_elements);

/**
 * Create a TPCHit collection
 */
edm4hep::RawTimeSeriesCollection createTPCHits(const int num_elements, const int num_rawwords);

/**
 * Create a TrackerHit collection
 */
edm4hep::TrackerHit3DCollection createTrackerHits(const int num_elements);

/**
 * Create a track collection with tracks that have links to other tracks (in the
 * same collection) and tracker hits
 */
edm4hep::TrackCollection createTracks(const int num_elements, const int subdetectorhitnumbers,
                                      const int num_track_states, const edm4hep::TrackerHit3DCollection& trackerHits,
                                      const std::vector<std::size_t>& link_trackerhit_idcs,
                                      const std::vector<test_config::IdxPair>& track_link_tracks_idcs);

/**
 * Create a SimCalorimeterHit (and an accompanying CaloHitContribution)
 * collection with links to MCParticles
 */
std::pair<edm4hep::SimCalorimeterHitCollection, edm4hep::CaloHitContributionCollection>
createSimCalorimeterHits(const int num_elements, const int num_contributions,
                         const edm4hep::MCParticleCollection& mcParticles,
                         const std::vector<test_config::CaloContIdx>& link_mcparticles_idcs);

edm4hep::EventHeaderCollection createEventHeader();

edm4hep::ClusterCollection createClusters(const int num_elements, const edm4hep::CalorimeterHitCollection& caloHits,
                                          const int num_subdet_energies,
                                          const std::vector<test_config::IdxPair>& clusterHitIdcs,
                                          const std::vector<test_config::IdxPair>& clusterClusterIdcs);

edm4hep::ReconstructedParticleCollection createRecoParticles(const int nRecos, const edm4hep::TrackCollection& tracks,
                                                             const std::vector<test_config::IdxPair>& trackIdcs,
                                                             const edm4hep::ClusterCollection& clusters,
                                                             const std::vector<test_config::IdxPair>& clusterIdcs,
                                                             const std::vector<test_config::IdxPair>& recIdcs);

std::vector<edm4hep::ParticleIDCollection>
createParticleIDs(const std::vector<std::vector<int>>& recoIdcs,
                  const edm4hep::ReconstructedParticleCollection& recoParticles);

edm4hep::RecoMCParticleLinkCollection
createMCRecoParticleLinks(const edm4hep::MCParticleCollection& mcParticles,
                          const edm4hep::ReconstructedParticleCollection& recoParticles);

edm4hep::CaloHitSimCaloHitLinkCollection createMCCaloAssocs(const edm4hep::SimCalorimeterHitCollection& simHits,
                                                            const edm4hep::CalorimeterHitCollection& caloHits);

/**
 * Create a Vertex collection and an accompanyin ReconstructedParticle
 * collection that uses them as decay vertices.
 */
std::tuple<edm4hep::VertexCollection, edm4hep::ReconstructedParticleCollection,
           edm4hep::VertexRecoParticleLinkCollection>
createVertices(const int nVertices, const edm4hep::ReconstructedParticleCollection& particles,
               const std::vector<test_config::IdxPair>& recoIdcs, const std::vector<test_config::IdxPair>& vtxRecoIdcs);

/**
 * Create an example event that can be used to test the converter. Also populate
 * a metadata frame that can be used in the conversion.
 *
 * Content:
 *
 * The following table gives an overview of the contents. The arguments for
 * calling the individual creation functions are taken from the test_config
 * namespace.
 *
 * | Name                 | Data type                       | comment                  |
 * |----------------------|---------------------------------|--------------------------|
 * | mcParticles          | MCParticle                      | createMCParticles        |
 * | caloHits             | CalorimeterHit                  | createCalorimeterHits    |
 * | rawCaloHits          | RawCalorimeterHit               | createRawCalorimeterHits |
 * | tpcHits              | RawTimeSeries                   | createTPCHits            |
 * | trackerHits          | TrackerHit                      | createTrackerHits        |
 * | simCaloHits          | SimCalorimeterHit               | createSimCalorimeterHits |
 * | caloHitContributions | CaloHitContribution             | createSimCalorimeterHits |
 * | clusters             | ClusterCollection               | createClusters           |
 * | recos                | ReconstructedParticleCollection | createRecoParticles      |
 * | recos_particleIDs    | ParticleIDCollection            | createRecoParticles      |
 */
std::tuple<podio::Frame, podio::Frame> createExampleEvent();

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
