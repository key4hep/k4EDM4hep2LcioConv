#ifndef K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H
#define K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H

#include "k4EDM4hep2LcioConv/MappingUtils.h"

// EDM4hep
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ClusterMCParticleLinkCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/VertexRecoParticleLinkCollection.h"
#include "edm4hep/utils/ParticleIDUtils.h"

// LCIO
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCIntVec.h>
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
#include <UTIL/LCIterator.h>
#include <lcio.h>

#include "podio/Frame.h"
#include "podio/UserDataCollection.h"

#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace LCIO2EDM4hepConv {
template <typename LcioT, typename EdmT>
using ObjectMapT = k4EDM4hep2LcioConv::MapT<LcioT, EdmT>;

/**
 * Maping holding all the original and converted objects in a 1:1 mapping in a
 * way that makes the lookup from LCIO to EDM4hep easy.
 */
struct LcioEdmTypeMapping {
  ObjectMapT<lcio::Track*, edm4hep::MutableTrack> tracks{};
  ObjectMapT<lcio::TrackerHit*, edm4hep::MutableTrackerHit3D> trackerHits{};
  ObjectMapT<lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit> simTrackerHits{};
  ObjectMapT<lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit> caloHits{};
  ObjectMapT<lcio::RawCalorimeterHit*, edm4hep::MutableRawCalorimeterHit> rawCaloHits{};
  ObjectMapT<lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit> simCaloHits{};
  ObjectMapT<lcio::TPCHit*, edm4hep::MutableRawTimeSeries> tpcHits{};
  ObjectMapT<lcio::Cluster*, edm4hep::MutableCluster> clusters{};
  ObjectMapT<lcio::Vertex*, edm4hep::MutableVertex> vertices{};
  ObjectMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle> recoParticles{};
  ObjectMapT<lcio::MCParticle*, edm4hep::MutableMCParticle> mcParticles{};
  ObjectMapT<lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane> trackerHitPlanes{};
};

using CollNamePair = std::tuple<std::string, std::unique_ptr<podio::CollectionBase>>;

/*
 * Convert a LCRunHeader to EDM4hep as a frame.
 */
podio::Frame convertRunHeader(EVENT::LCRunHeader* rheader);

/**
 * Convert a complete LCEvent from LCIO to EDM4hep.
 *
 * A second, optional argument can be passed to limit the collections to convert
 * to the subset that is passed. Additionally, it allows to rename collections
 * on the fly where the first element of each pair is the (original) LCIO name
 * and the second one is the one that is used for the EDM4hep collection.
 *
 * NOTE: There is an implicit assumption here that collsToConvert only contains
 * collection names that are present in the passed evt. There is no exception
 * handling internally to guard against collections that are missing.
 */
podio::Frame convertEvent(EVENT::LCEvent* evt, const std::vector<std::pair<std::string, std::string>>& = {});

/**
 * Convert an LCIOCollection by dispatching to the specific conversion
 * function for the corresponding type (after querying the input collection).
 * Populates the correct object mapping along the way.
 *
 * Returns a vector of names and collections (since some LCIO collections will
 * result in more than one EDM4hep collection)
 */
template <typename ObjectMappingT>
std::vector<CollNamePair> convertCollection(const std::string& name, EVENT::LCCollection* LCCollection,
                                            ObjectMappingT& typeMapping);

/**
 * Resolve all relations in all converted objects that are held in the map.
 * Dispatch to the correpsonding implementation for all the types that have
 * relations
 */
template <typename ObjectMappingT>
void resolveRelations(ObjectMappingT& typeMapping);

template <typename ObjectMappingT, typename ObjectMappingU>
void resolveRelations(ObjectMappingT& updateMaps, const ObjectMappingU& lookupMaps);

/**
 * Convert LCRelation collections into the corresponding Link collections in
 * EDM4hep
 */
template <typename ObjectMappingT>
std::vector<CollNamePair> createLinks(const ObjectMappingT& typeMapping,
                                      const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation);

/**
 * Convert a subset collection, dispatching to the correct function for the
 * type of the input collection
 */
template <typename ObjectMappingT>
std::unique_ptr<podio::CollectionBase> fillSubset(EVENT::LCCollection* LCCollection, const ObjectMappingT& typeMapping,
                                                  const std::string& type);

/*
 * Converts a LCIntVec or LCFloatVec Collection into a podio::UserDataCollection
 * of the appropriate type.
 *
 * NOTE: LC[Int|Float]Vec are nested, but podio::UserDataCollection are flat.
 * Hence, this will put all contents into one collection, and the [begin, end)
 * indices in this collection into a second (flat) collection (with the suffix
 * "_VecLengths" added to its name), such that the elements at position i, resp.
 * (i + 1) form the [begin, end) indices for each of the original vector
 * collections.
 */
template <typename LCVecType>
std::vector<CollNamePair> convertLCVec(const std::string& name, EVENT::LCCollection* LCCollection);

/**
 * Helper struct to wrap the functionality of putting parameters into a Frame in
 * a semi type-erases way, such that it can easily be used as a template
 * parameter
 */
struct ParamFramePutter {
  ParamFramePutter() = delete;
  ParamFramePutter(podio::Frame& f) : m_frame(f) {}

  template <typename T>
  void operator()(const std::string& key, const T& value) {
    m_frame.putParameter(key, value);
  }

private:
  podio::Frame& m_frame;
};

/**
 * Converting all parameters of an LCIO Object and passing them to the PutParamF
 * function that takes care of storing them appropriately.
 *
 * The indirection is necessary for better integration with k4FWCore where
 * direct access to a Frame is not possible, but a putParameter method is
 * available instead.
 */
template <typename LCIOType, typename PutParamF = ParamFramePutter>
void convertObjectParameters(LCIOType* lcioobj, PutParamF putParamFun);

/**
 * Converting all parameters of an LCIO Object and attaching them to the
 * passed podio::Frame.
 */
template <typename LCIOType>
void convertObjectParameters(LCIOType* lcioobj, podio::Frame& event) {
  convertObjectParameters(lcioobj, ParamFramePutter{event});
}

inline edm4hep::Vector3f Vector3fFrom(const double* v) { return edm4hep::Vector3f(v[0], v[1], v[2]); }

inline edm4hep::Vector3f Vector3fFrom(const EVENT::FloatVec& v) { return edm4hep::Vector3f(v[0], v[1], v[2]); }

/**
 * Get the name of a ParticleID collection from the name of the reco
 * collection (from which it is created) and the PID algorithm name.
 */
inline std::string getPIDCollName(const std::string& recoCollName, const std::string& algoName) {
  return recoCollName + "_PID_" + algoName;
}

/**
 * Get the meta information for all particle id collections that are available
 * from the PIDHandler
 */
std::vector<edm4hep::utils::ParticleIDMeta> getPIDMetaInfo(const EVENT::LCCollection* recoColl);

/**
 * Convert a TrackState
 */
edm4hep::TrackState convertTrackState(const EVENT::TrackState* trackState);

/**
 * Convert a ParticleID object.
 *
 * In LCIO ParticleIDs are persisted as part of the ReconstructedParticle they
 * are attached to. There are no ParticleID collections in the event. Hence,
 * this function converts single particle ID objects and the management of
 * putting them into a collection and of creating the LCIO to EDM4hep mapping
 * is done in the conversion of the ReconstructedParticles.
 */
edm4hep::MutableParticleID convertParticleID(const EVENT::ParticleID* pid);

/**
 * Convert an MCParticle collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename MCParticleMapT>
std::unique_ptr<edm4hep::MCParticleCollection>
convertMCParticles(const std::string& name, EVENT::LCCollection* LCCollection, MCParticleMapT& mcparticlesMap);

/**
 * Convert a ReconstructedParticle collection and return the resulting
 * collection. Simultaneously populates the mapping from LCIO to EDM4hep
 * objects.
 *
 * @note: Also populates ParticleID collections, as those are persisted as
 * part of the ReconstructedParticles in LCIO. The name of this collection is
 * <name>_PID_<pid_algo_name> (see getPIDCollName)
 *
 * @note: Also populates one (partially filled) VertexRecoParticleLink
 * collection for keeping the startVertex information
 */
template <typename RecoMapT>
std::vector<CollNamePair> convertReconstructedParticles(const std::string& name, EVENT::LCCollection* LCCollection,
                                                        RecoMapT& recoparticlesMap);

/**
 * Convert a Vertex collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 *
 * @note: Also creates a (partially filled) VertexRecoParticleLink
 * collection for keeping the associatedParticle information
 */
template <typename VertexMapT>
std::vector<CollNamePair> convertVertices(const std::string& name, EVENT::LCCollection* LCCollection,
                                          VertexMapT& vertexMap);

/**
 * Convert a SimTrackerHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename SimTrHitMapT>
std::unique_ptr<edm4hep::SimTrackerHitCollection>
convertSimTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, SimTrHitMapT& SimTrHitMap);

/**
 * Convert a TPCHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::RawTimeSeriesCollection> convertTPCHits(const std::string& name,
                                                                 EVENT::LCCollection* LCCollection, HitMapT& TPCHitMap);

/**
 * Convert a TrackerHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::TrackerHit3DCollection>
convertTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitMap);

/**
 * Convert a TrackerHitPlane collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::TrackerHitPlaneCollection>
convertTrackerHitPlanes(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitPlaneMap);

/**
 * Convert a Track collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename TrackMapT>
std::vector<CollNamePair> convertTracks(const std::string& name, EVENT::LCCollection* LCCollection,
                                        TrackMapT& TrackMap);

/**
 * Convert a SimCalorimeterHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::SimCalorimeterHitCollection>
convertSimCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& SimCaloHitMap);

/**
 * Convert a RawCalorimeterHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::RawCalorimeterHitCollection>
convertRawCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& rawCaloHitMap);

/**
 * Convert a CalorimeterHit collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename HitMapT>
std::unique_ptr<edm4hep::CalorimeterHitCollection>
convertCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& caloHitMap);

/**
 * Convert a Cluster collection and return the resulting collection.
 * Simultaneously populates the mapping from LCIO to EDM4hep objects.
 */
template <typename ClusterMapT>
std::unique_ptr<edm4hep::ClusterCollection> convertClusters(const std::string& name, EVENT::LCCollection* LCCollection,
                                                            ClusterMapT& clusterMap);

/**
 * Create an EventHeaderCollection and fills it with the Metadata.
 */

std::unique_ptr<edm4hep::EventHeaderCollection> createEventHeader(const EVENT::LCEvent* evt);

/**
 * Helper function to create a subset collection from an existing (LCIO)
 * subset collection. Needs the object mapping as input to resolve to the
 * correct EDM4hep object for a given LCIO object.
 *
 * NOTE: Users responsibility to call this with the right inputs (i.e.
 * matching types)
 */
template <typename CollT, typename ObjectMapT,
          typename LcioT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<ObjectMapT>>,
          typename Edm4hepT = k4EDM4hep2LcioConv::detail::mapped_t<ObjectMapT>>
auto handleSubsetColl(EVENT::LCCollection* lcioColl, const ObjectMapT& elemMap);

/**
 * Create an Link collection from an LCRelations collection. Templated
 * on the From and To types as well as the direction of the relations in the
 * input LCRelations collection with respect to the order in which they are
 * mentioned in the Link collection of EDM4hep (since those are not
 * directed).
 *
 * Necessary inputs apart from the LCRelations collection are the correct LCIO
 * to EDM4hep object mappings to actually resolve the necessary relations.
 */
template <typename CollT, bool Reverse, typename FromMapT, typename ToMapT,
          typename FromLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<FromMapT>>,
          typename ToLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<ToMapT>>,
          typename FromEDM4hepT = k4EDM4hep2LcioConv::detail::mapped_t<FromMapT>,
          typename ToEDM4hepT = k4EDM4hep2LcioConv::detail::mapped_t<ToMapT>>
std::unique_ptr<CollT> createLinkCollection(EVENT::LCCollection* relations, const FromMapT& fromMap,
                                            const ToMapT& toMap);

/**
 * Creates the CaloHitContributions for all SimCaloHits.
 * has to be done this way, since the converted McParticles are needeed.
 * The contributions are also attached to their corresponding
 * SimCalorimeterHits.
 */
template <typename HitMapT, typename MCParticleMapT>
std::unique_ptr<edm4hep::CaloHitContributionCollection>
createCaloHitContributions(HitMapT& SimCaloHitMap, const MCParticleMapT& mcparticlesMap);

/**
 * Resolve the relations for the MCParticles.
 */
template <typename MCParticleMapT, typename MCParticleLookupMapT>
void resolveRelationsMCParticles(MCParticleMapT& mcparticlesMap, const MCParticleLookupMapT& lookupMap);

/**
 * Resolve the relations for SimTrackerHits
 */
template <typename HitMapT, typename MCParticleMapT>
void resolveRelationsSimTrackerHits(HitMapT& SimTrHitMap, const MCParticleMapT& mcparticlesMap);

/**
 * Resolve the relations for ReconstructedParticles
 */
template <typename RecoParticleMapT, typename RecoParticleLookupMapT, typename ClusterMapT, typename TrackMapT>
void resolveRelationsRecoParticles(RecoParticleMapT& recoparticlesMap, const RecoParticleLookupMapT& recoLookupMap,
                                   const ClusterMapT& clusterMap, const TrackMapT& tracksMap);

/**
 * Resolve the relations for Clusters
 */
template <typename ClusterMapT, typename CaloHitMapT>
void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap);

/**
 * Resolve the relations for Tracks
 */
template <typename TrackMapT, typename TrackHitMapT, typename THPlaneHitMapT, typename TPCHitMapT>
void resolveRelationsTracks(TrackMapT& tracksMap, const TrackHitMapT& trackerHitMap, const THPlaneHitMapT&,
                            const TPCHitMapT&);

/**
 * Resolve the relations for Vertices. Vertex related information in
 * reconstructed particles will only be mutated in the updateRPMap.
 */
template <typename VertexMapT, typename URecoParticleMapT, typename LURecoParticleMapT>
void resolveRelationsVertices(VertexMapT& vertexMap, URecoParticleMapT& updateRPMap,
                              const LURecoParticleMapT& lookupRPMap);

template <typename VertexMapT, typename RecoParticleMapT>
void finalizeVertexRecoParticleLinks(edm4hep::VertexRecoParticleLinkCollection& links, const VertexMapT& vertexMap,
                                     const RecoParticleMapT& recoParticleMap);

/**
 * Go from chi^2 and probability (1 - CDF(chi^2, ndf)) to ndf by a binary search
 */
int find_ndf(double chi2, double prob);

} // namespace LCIO2EDM4hepConv

#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.ipp"

#endif // K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H
