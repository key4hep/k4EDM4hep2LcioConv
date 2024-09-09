#include "k4EDM4hep2LcioConv/MappingUtils.h"

#include <edm4hep/CaloHitMCParticleLinkCollection.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/ClusterMCParticleLinkCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/VertexRecoParticleLinkCollection.h>

#include <UTIL/LCRelationNavigator.h>

#include "TMath.h"

namespace EDM4hep2LCIOConv {

template <typename TrackMapT>
std::unique_ptr<lcio::LCCollectionVec> convertTracks(const edm4hep::TrackCollection* const edmCollection,
                                                     TrackMapT& trackMap) {
  auto tracks = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::TRACK);
  // Loop over EDM4hep tracks converting them to lcio tracks.
  for (const auto& edm_tr : (*edmCollection)) {
    if (edm_tr.isAvailable()) {
      auto* lcio_tr = new lcio::TrackImpl();
      // The Type of the Tracks need to be set bitwise in LCIO since the
      // setType(int) function is private for the LCIO TrackImpl and only a
      // setTypeBit(bitnumber) function can be used to set the Type bit by bit.
      int type = edm_tr.getType();
      for (auto i = 0u; i < sizeof(int) * 8; i++) {
        lcio_tr->setTypeBit(i, type & (1 << i));
      }
      lcio_tr->setChi2(edm_tr.getChi2());
      lcio_tr->setNdf(edm_tr.getNdf());
      lcio_tr->setRadiusOfInnermostHit(getRadiusOfStateAtFirstHit(edm_tr).value_or(-1.0));

      // Loop over the hit Numbers in the track
      lcio_tr->subdetectorHitNumbers().resize(edm_tr.subdetectorHitNumbers_size());
      for (auto i = 0u; i < edm_tr.subdetectorHitNumbers_size(); ++i) {
        lcio_tr->subdetectorHitNumbers()[i] = edm_tr.getSubdetectorHitNumbers(i);
      }

      // Pad until 50 hitnumbers are resized
      const int hit_number_limit = 50;
      if (edm_tr.subdetectorHitNumbers_size() < hit_number_limit) {
        lcio_tr->subdetectorHitNumbers().resize(hit_number_limit);
        for (auto i = edm_tr.subdetectorHitNumbers_size(); i < hit_number_limit; ++i) {
          lcio_tr->subdetectorHitNumbers()[i] = 0;
        }
      }

      lcio_tr->setNholes(edm_tr.getNholes());
      const auto edmHoleNumbers = edm_tr.getSubdetectorHoleNumbers();
      lcio_tr->subdetectorHoleNumbers().resize(edmHoleNumbers.size());
      for (auto i = 0u; i < edmHoleNumbers.size(); ++i) {
        lcio_tr->subdetectorHoleNumbers()[i] = edmHoleNumbers[i];
      }

      // Loop over the track states in the track
      const auto edm_track_states = edm_tr.getTrackStates();
      for (const auto& tr_state : edm_track_states) {
        const auto& cov = tr_state.covMatrix;
        std::array<float, 3> refP = {tr_state.referencePoint.x, tr_state.referencePoint.y, tr_state.referencePoint.z};

        auto* lcio_tr_state = new lcio::TrackStateImpl(tr_state.location, tr_state.D0, tr_state.phi, tr_state.omega,
                                                       tr_state.Z0, tr_state.tanLambda, cov.data(), refP.data());

        lcio_tr->addTrackState(lcio_tr_state);
      }

      // Save intermediate tracks ref
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_tr, edm_tr, trackMap);

      // Add to lcio tracks collection
      tracks->addElement(lcio_tr);
    }
  }

  return tracks;
}

template <typename TrackerHitMapT>
std::unique_ptr<lcio::LCCollectionVec> convertTrackerHits(const edm4hep::TrackerHit3DCollection* const edmColection,
                                                          const CellIDStrType& cellIDstr,
                                                          TrackerHitMapT& trackerHitMap) {
  auto trackerhits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::TRACKERHIT);

  if (cellIDstr.has_value()) {
    lcio::CellIDEncoder<lcio::SimCalorimeterHitImpl> idEnc(cellIDstr.value(), trackerhits.get());
  }

  // Loop over EDM4hep trackerhits converting them to lcio trackerhits
  for (const auto& edm_trh : (*edmColection)) {
    if (edm_trh.isAvailable()) {
      auto* lcio_trh = new lcio::TrackerHitImpl();

      uint64_t combined_value = edm_trh.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_trh->setCellID0(combined_value_ptr[0]);
      lcio_trh->setCellID1(combined_value_ptr[1]);
      lcio_trh->setType(edm_trh.getType());
      std::array<double, 3> positions{edm_trh.getPosition()[0], edm_trh.getPosition()[1], edm_trh.getPosition()[2]};
      lcio_trh->setPosition(positions.data());
      lcio_trh->setCovMatrix(edm_trh.getCovMatrix().data());
      lcio_trh->setEDep(edm_trh.getEDep());
      lcio_trh->setEDepError(edm_trh.getEDepError());
      lcio_trh->setTime(edm_trh.getTime());
      lcio_trh->setQuality(edm_trh.getQuality());

      // Save intermediate trackerhits ref
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_trh, edm_trh, trackerHitMap);

      // Add to lcio trackerhits collection
      trackerhits->addElement(lcio_trh);
    }
  }

  return trackerhits;
}

template <typename TrackerHitPlaneMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertTrackerHitPlanes(const edm4hep::TrackerHitPlaneCollection* const edmCollection, const CellIDStrType& cellIDstr,
                        TrackerHitPlaneMapT& trackerHitsMap) {
  auto trackerHitPlanes = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::TRACKERHITPLANE);

  if (cellIDstr.has_value()) {
    lcio::CellIDEncoder<lcio::SimCalorimeterHitImpl> idEnc(cellIDstr.value(), trackerHitPlanes.get());
  }

  for (const auto& edm_trh : (*edmCollection)) {
    if (edm_trh.isAvailable()) {
      auto* lcio_trh = new lcio::TrackerHitPlaneImpl();

      uint64_t combined_value = edm_trh.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_trh->setCellID0(combined_value_ptr[0]);
      lcio_trh->setCellID1(combined_value_ptr[1]);
      lcio_trh->setType(edm_trh.getType());
      const std::array positions{edm_trh.getPosition()[0], edm_trh.getPosition()[1], edm_trh.getPosition()[2]};
      lcio_trh->setPosition(positions.data());
      // No public setter in LCIO
      // lcio_trh->setCovMatrix(edm_trh.getCovMatrix().data());
      lcio_trh->setEDep(edm_trh.getEDep());
      lcio_trh->setEDepError(edm_trh.getEDepError());
      lcio_trh->setTime(edm_trh.getTime());
      lcio_trh->setQuality(edm_trh.getQuality());

      const std::array posU{edm_trh.getU()[0], edm_trh.getU()[1]};
      lcio_trh->setU(posU.data());
      lcio_trh->setdU(edm_trh.getDu());
      const std::array posV{edm_trh.getV()[0], edm_trh.getV()[1]};
      lcio_trh->setV(posV.data());
      lcio_trh->setdV(edm_trh.getDv());

      k4EDM4hep2LcioConv::detail::mapInsert(lcio_trh, edm_trh, trackerHitsMap);

      trackerHitPlanes->addElement(lcio_trh);
    }
  }

  return trackerHitPlanes;
}

template <typename SimTrHitMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertSimTrackerHits(const edm4hep::SimTrackerHitCollection* const edmCollection, const CellIDStrType& cellIDstr,
                      SimTrHitMapT& simTrHitMap) {
  auto simtrackerhits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::SIMTRACKERHIT);

  if (cellIDstr.has_value()) {
    lcio::CellIDEncoder<lcio::SimCalorimeterHitImpl> idEnc(cellIDstr.value(), simtrackerhits.get());
  }

  // Loop over EDM4hep simtrackerhits converting them to LCIO simtrackerhits
  for (const auto& edm_strh : (*edmCollection)) {
    if (edm_strh.isAvailable()) {
      auto* lcio_strh = new lcio::SimTrackerHitImpl();

      uint64_t combined_value = edm_strh.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_strh->setCellID0(combined_value_ptr[0]);
      lcio_strh->setCellID1(combined_value_ptr[1]);
      std::array<double, 3> positions{edm_strh.getPosition()[0], edm_strh.getPosition()[1], edm_strh.getPosition()[2]};
      lcio_strh->setPosition(positions.data());
      lcio_strh->setEDep(edm_strh.getEDep());
      lcio_strh->setTime(edm_strh.getTime());
      lcio_strh->setMomentum(edm_strh.getMomentum()[0], edm_strh.getMomentum()[1], edm_strh.getMomentum()[2]);
      lcio_strh->setPathLength(edm_strh.getPathLength());
      lcio_strh->setQuality(edm_strh.getQuality());
      // lcio_strh->setQualityBit( int bit , bool val=true ) ;
      lcio_strh->setOverlay(edm_strh.isOverlay());
      lcio_strh->setProducedBySecondary(edm_strh.isProducedBySecondary());

      // Save intermediate simtrackerhits ref
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_strh, edm_strh, simTrHitMap);

      // Add to lcio simtrackerhits collection
      simtrackerhits->addElement(lcio_strh);
    }
  }

  return simtrackerhits;
}

// Convert EDM4hep Calorimeter Hits to LCIO
// Add converted LCIO ptr and original EDM4hep collection to vector of pairs
// Add converted LCIO Collection Vector to LCIO event
template <typename CaloHitMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertCalorimeterHits(const edm4hep::CalorimeterHitCollection* const edmCollection, const CellIDStrType& cellIDstr,
                       CaloHitMapT& caloHitMap) {
  auto calohits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::CALORIMETERHIT);

  if (cellIDstr.has_value()) {
    lcio::CellIDEncoder<lcio::SimCalorimeterHitImpl> idEnc(cellIDstr.value(), calohits.get());
  }

  for (const auto& edm_calohit : (*edmCollection)) {
    if (edm_calohit.isAvailable()) {
      auto* lcio_calohit = new lcio::CalorimeterHitImpl();

      uint64_t combined_value = edm_calohit.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_calohit->setCellID0(combined_value_ptr[0]);
      lcio_calohit->setCellID1(combined_value_ptr[1]);
      lcio_calohit->setEnergy(edm_calohit.getEnergy());
      lcio_calohit->setEnergyError(edm_calohit.getEnergyError());
      lcio_calohit->setTime(edm_calohit.getTime());
      std::array<float, 3> positions{edm_calohit.getPosition()[0], edm_calohit.getPosition()[1],
                                     edm_calohit.getPosition()[2]};
      lcio_calohit->setPosition(positions.data());
      lcio_calohit->setType(edm_calohit.getType());

      // TODO
      // lcio_calohit->setRawHit(EVENT::LCObject* rawHit );

      // Save Calorimeter Hits LCIO and EDM4hep collections
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_calohit, edm_calohit, caloHitMap);

      // Add to lcio tracks collection
      calohits->addElement(lcio_calohit);
    }
  }

  return calohits;
}

template <typename RawCaloHitMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertRawCalorimeterHits(const edm4hep::RawCalorimeterHitCollection* const edmCollection,
                          RawCaloHitMapT& rawCaloHitMap) {
  auto rawcalohits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::RAWCALORIMETERHIT);

  for (const auto& edm_raw_calohit : (*edmCollection)) {
    if (edm_raw_calohit.isAvailable()) {
      auto* lcio_rawcalohit = new lcio::RawCalorimeterHitImpl();

      uint64_t combined_value = edm_raw_calohit.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_rawcalohit->setCellID0(combined_value_ptr[0]);
      lcio_rawcalohit->setCellID1(combined_value_ptr[1]);
      lcio_rawcalohit->setAmplitude(edm_raw_calohit.getAmplitude());
      lcio_rawcalohit->setTimeStamp(edm_raw_calohit.getTimeStamp());

      // Save Raw Calorimeter Hits LCIO and EDM4hep collections
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_rawcalohit, edm_raw_calohit, rawCaloHitMap);

      // Add to lcio tracks collection
      rawcalohits->addElement(lcio_rawcalohit);
    }
  }

  return rawcalohits;
}

template <typename SimCaloHitMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertSimCalorimeterHits(const edm4hep::SimCalorimeterHitCollection* const edmCollection,
                          const CellIDStrType& cellIDstr, SimCaloHitMapT& simCaloHitMap) {
  auto simcalohits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::SIMCALORIMETERHIT);

  if (cellIDstr.has_value()) {
    lcio::CellIDEncoder<lcio::SimCalorimeterHitImpl> idEnc(cellIDstr.value(), simcalohits.get());
  }

  for (const auto& edm_sim_calohit : (*edmCollection)) {
    if (edm_sim_calohit.isAvailable()) {
      auto* lcio_simcalohit = new lcio::SimCalorimeterHitImpl();

      uint64_t combined_value = edm_sim_calohit.getCellID();
      uint32_t* combined_value_ptr = reinterpret_cast<uint32_t*>(&combined_value);
      lcio_simcalohit->setCellID0(combined_value_ptr[0]);
      lcio_simcalohit->setCellID1(combined_value_ptr[1]);
      lcio_simcalohit->setEnergy(edm_sim_calohit.getEnergy());
      std::array<float, 3> positions{edm_sim_calohit.getPosition()[0], edm_sim_calohit.getPosition()[1],
                                     edm_sim_calohit.getPosition()[2]};
      lcio_simcalohit->setPosition(positions.data());

      // Contributions are converted in resolveRelations to make it a higher
      // probability that we have the MCParticles converted

      // Save Sim Calorimeter Hits LCIO and EDM4hep collections
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_simcalohit, edm_sim_calohit, simCaloHitMap);

      // Add to sim calo hits collection
      simcalohits->addElement(lcio_simcalohit);
    }
  }

  return simcalohits;
}

template <typename TPCHitMapT>
std::unique_ptr<lcio::LCCollectionVec> convertTPCHits(const edm4hep::RawTimeSeriesCollection* const edmCollection,
                                                      TPCHitMapT& tpcHitMap) {
  auto tpchits = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::TPCHIT);

  for (const auto& edm_tpchit : (*edmCollection)) {
    if (edm_tpchit.isAvailable()) {
      auto* lcio_tpchit = new lcio::TPCHitImpl();

#pragma message "unsigned long long conversion to int"

      lcio_tpchit->setCellID(edm_tpchit.getCellID());
      lcio_tpchit->setTime(edm_tpchit.getTime());
      lcio_tpchit->setCharge(edm_tpchit.getCharge());
      lcio_tpchit->setQuality(edm_tpchit.getQuality());

      std::vector<int> rawdata;
      for (auto i = 0u; i < edm_tpchit.adcCounts_size(); ++i) {
        rawdata.push_back(edm_tpchit.getAdcCounts(i));
      }

      lcio_tpchit->setRawData(rawdata.data(), edm_tpchit.adcCounts_size());

      // Save TPC Hits LCIO and EDM4hep collections
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_tpchit, edm_tpchit, tpcHitMap);

      // Add to lcio tracks collection
      tpchits->addElement(lcio_tpchit);
    }
  }

  return tpchits;
}

template <typename ClusterMapT>
std::unique_ptr<lcio::LCCollectionVec> convertClusters(const edm4hep::ClusterCollection* const edmCollection,
                                                       ClusterMapT& clusterMap) {
  auto clusters = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::CLUSTER);

  // Loop over EDM4hep clusters converting them to lcio clusters
  for (const auto& edm_cluster : (*edmCollection)) {
    if (edm_cluster.isAvailable()) {
      auto* lcio_cluster = new lcio::ClusterImpl();

      std::bitset<sizeof(uint32_t)> type_bits = edm_cluster.getType();
      for (auto j = 0u; j < sizeof(uint32_t); j++) {
        lcio_cluster->setTypeBit(j, (type_bits[j] == 0) ? false : true);
      }
      lcio_cluster->setEnergy(edm_cluster.getEnergy());
      lcio_cluster->setEnergyError(edm_cluster.getEnergyError());

      std::array<float, 3> edm_cluster_pos = {edm_cluster.getPosition().x, edm_cluster.getPosition().y,
                                              edm_cluster.getPosition().z};
      lcio_cluster->setPosition(edm_cluster_pos.data());

      lcio_cluster->setPositionError(edm_cluster.getPositionError().data());
      lcio_cluster->setITheta(edm_cluster.getITheta());
      lcio_cluster->setIPhi(edm_cluster.getPhi());
      std::array<float, 3> edm_cluster_dir_err = {edm_cluster.getPosition().x, edm_cluster.getPosition().y,
                                                  edm_cluster.getPosition().z};
      lcio_cluster->setDirectionError(edm_cluster_dir_err.data());

      EVENT::FloatVec shape_vec{};
      for (auto& param : edm_cluster.getShapeParameters()) {
        shape_vec.push_back(param);
      }
      lcio_cluster->setShape(shape_vec);

      auto& subdetEnergies = lcio_cluster->subdetectorEnergies();
      for (const auto edmEnergy : edm_cluster.getSubdetectorEnergies()) {
        subdetEnergies.push_back(edmEnergy);
      }

      // Add LCIO and EDM4hep pair collections to vec
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_cluster, edm_cluster, clusterMap);

      // Add to lcio tracks collection
      clusters->addElement(lcio_cluster);
    }
  }

  return clusters;
}

template <typename VertexMapT>
std::unique_ptr<lcio::LCCollectionVec> convertVertices(const edm4hep::VertexCollection* const edmCollection,
                                                       VertexMapT& vertexMap) {
  auto vertices = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::VERTEX);

  // Loop over EDM4hep vertex converting them to lcio vertex
  for (const auto& edm_vertex : (*edmCollection)) {
    if (edm_vertex.isAvailable()) {
      auto* lcio_vertex = new lcio::VertexImpl();
      lcio_vertex->setPrimary(edm_vertex.isPrimary());
      lcio_vertex->setAlgorithmType(std::to_string(edm_vertex.getAlgorithmType()));
      lcio_vertex->setChi2(edm_vertex.getChi2());
      lcio_vertex->setProbability(TMath::Prob(edm_vertex.getChi2(), edm_vertex.getNdf()));
      lcio_vertex->setPosition(edm_vertex.getPosition()[0], edm_vertex.getPosition()[1], edm_vertex.getPosition()[2]);
      lcio_vertex->setCovMatrix(edm_vertex.getCovMatrix().data());

      for (auto& param : edm_vertex.getParameters()) {
        lcio_vertex->addParameter(param);
      }

      // Add LCIO and EDM4hep pair collections to vec
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_vertex, edm_vertex, vertexMap);

      // Add to lcio tracks collection
      vertices->addElement(lcio_vertex);
    }
  }

  return vertices;
}

template <typename RecoPartMapT>
std::unique_ptr<lcio::LCCollectionVec>
convertReconstructedParticles(const edm4hep::ReconstructedParticleCollection* const recos_coll,
                              RecoPartMapT& recoparticles_vec) {
  auto recops = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::RECONSTRUCTEDPARTICLE);

  for (const auto& edm_rp : (*recos_coll)) {
    auto* lcio_recp = new lcio::ReconstructedParticleImpl;
    if (edm_rp.isAvailable()) {
      lcio_recp->setType(edm_rp.getPDG());
      float m[3] = {edm_rp.getMomentum()[0], edm_rp.getMomentum()[1], edm_rp.getMomentum()[2]};
      lcio_recp->setMomentum(m);
      lcio_recp->setEnergy(edm_rp.getEnergy());
      lcio_recp->setCovMatrix(edm_rp.getCovMatrix().data()); // TODO Check lower or upper
      lcio_recp->setMass(edm_rp.getMass());
      lcio_recp->setCharge(edm_rp.getCharge());
      float rp[3] = {edm_rp.getReferencePoint()[0], edm_rp.getReferencePoint()[1], edm_rp.getReferencePoint()[2]};
      lcio_recp->setReferencePoint(rp);
      lcio_recp->setGoodnessOfPID(edm_rp.getGoodnessOfPID());
      // Add LCIO and EDM4hep pair collections to vec
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_recp, edm_rp, recoparticles_vec);

      // Add to reconstructed particles collection
      recops->addElement(lcio_recp);
    }
  }

  return recops;
}

template <typename MCPartMapT>
std::unique_ptr<lcio::LCCollectionVec> convertMCParticles(const edm4hep::MCParticleCollection* const edmCollection,
                                                          MCPartMapT& mcParticleMap) {
  auto mcparticles = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::MCPARTICLE);

  for (const auto& edm_mcp : (*edmCollection)) {
    auto* lcio_mcp = new lcio::MCParticleImpl;
    if (edm_mcp.isAvailable()) {
      lcio_mcp->setPDG(edm_mcp.getPDG());
      lcio_mcp->setGeneratorStatus(edm_mcp.getGeneratorStatus());
      // Note LCIO sets some Bits during writing which makes a trivial integer
      // conversion afterwards not work
      int status = edm_mcp.getSimulatorStatus();
      lcio_mcp->setSimulatorStatus(status);

      double vertex[3] = {edm_mcp.getVertex()[0], edm_mcp.getVertex()[1], edm_mcp.getVertex()[2]};
      lcio_mcp->setVertex(vertex);
      lcio_mcp->setTime(edm_mcp.getTime());
      double endpoint[3] = {edm_mcp.getEndpoint()[0], edm_mcp.getEndpoint()[1], edm_mcp.getEndpoint()[2]};
      lcio_mcp->setEndpoint(endpoint);
      double momentum[3] = {edm_mcp.getMomentum()[0], edm_mcp.getMomentum()[1], edm_mcp.getMomentum()[2]};
      lcio_mcp->setMomentum(momentum);
      double momentumEndpoint[3] = {edm_mcp.getMomentumAtEndpoint()[0], edm_mcp.getMomentumAtEndpoint()[1],
                                    edm_mcp.getMomentumAtEndpoint()[2]};
      lcio_mcp->setMomentumAtEndpoint(momentumEndpoint);
      lcio_mcp->setMass(edm_mcp.getMass());
      lcio_mcp->setCharge(edm_mcp.getCharge());
      float spin[3] = {edm_mcp.getSpin()[0], edm_mcp.getSpin()[1], edm_mcp.getSpin()[2]};
      lcio_mcp->setSpin(spin);
      int colorflow[2] = {edm_mcp.getColorFlow()[0], edm_mcp.getColorFlow()[1]};
      lcio_mcp->setColorFlow(colorflow);

      lcio_mcp->setCreatedInSimulation(edm_mcp.isCreatedInSimulation());
      lcio_mcp->setBackscatter(edm_mcp.isBackscatter());
      lcio_mcp->setVertexIsNotEndpointOfParent(edm_mcp.vertexIsNotEndpointOfParent());
      lcio_mcp->setDecayedInTracker(edm_mcp.isDecayedInTracker());
      lcio_mcp->setDecayedInCalorimeter(edm_mcp.isDecayedInCalorimeter());
      lcio_mcp->setHasLeftDetector(edm_mcp.hasLeftDetector());
      lcio_mcp->setStopped(edm_mcp.isStopped());
      lcio_mcp->setOverlay(edm_mcp.isOverlay());

      // Add LCIO and EDM4hep pair collections to vec
      k4EDM4hep2LcioConv::detail::mapInsert(lcio_mcp, edm_mcp, mcParticleMap);

      // Add to reconstructed particles collection
      mcparticles->addElement(lcio_mcp);
    }
  }

  return mcparticles;
}

template <typename PidMapT>
void convertParticleIDs(const edm4hep::ParticleIDCollection* const edmCollection, PidMapT& pidMap, const int algoId) {
  for (const auto& edmPid : (*edmCollection)) {
    auto [lcioPid, _] = k4EDM4hep2LcioConv::detail::mapInsert(new lcio::ParticleIDImpl(), edmPid, pidMap).first;

    lcioPid->setType(edmPid.getType());
    lcioPid->setPDG(edmPid.getPDG());
    lcioPid->setLikelihood(edmPid.getLikelihood());
    lcioPid->setAlgorithmType(algoId);
    for (const auto& param : edmPid.getParameters()) {
      lcioPid->addParameter(param);
    }
  }
}

template <typename MCParticleMapT, typename MCParticleLookupMapT>
void resolveRelationsMCParticles(MCParticleMapT& mcParticlesMap, const MCParticleLookupMapT& lookupMap) {
  // Add parent MCParticles after converting all MCparticles
  for (auto& [lcio_mcp, edm_mcp] : mcParticlesMap) {
    for (const auto& emd_parent_mcp : edm_mcp.getParents()) {
      if (emd_parent_mcp.isAvailable()) {
        // Search for the parent mcparticle in the converted vector
        if (const auto lcio_mcp_linked = k4EDM4hep2LcioConv::detail::mapLookupFrom(emd_parent_mcp, lookupMap)) {
          lcio_mcp->addParent(lcio_mcp_linked.value());
        }
      }
    }
  }
}

template <typename TrackMapT, typename TrackHitMapT, typename THPlaneHitMapT, typename TPCHitMapT>
void resolveRelationsTracks(TrackMapT& tracksMap, const TrackHitMapT& trackerHitMap,
                            const THPlaneHitMapT& trackerHitPlaneMap, const TPCHitMapT&) {
  for (auto& [lcio_tr, edm_tr] : tracksMap) {
    auto tracks = edm_tr.getTracks();
    for (auto t : tracks) {
      if (!t.isAvailable()) {
        continue;
      }
      if (const auto track = k4EDM4hep2LcioConv::detail::mapLookupFrom(t, tracksMap)) {
        lcio_tr->addTrack(track.value());
      }
    }

    auto trackerHits = edm_tr.getTrackerHits();
    for (auto th : trackerHits) {
      if (!th.isAvailable()) {
        continue;
      }
      if (const auto hit = k4EDM4hep2LcioConv::detail::mapLookupFrom(th, trackerHitMap)) {
        lcio_tr->addHit(hit.value());
      } else if (const auto hit = k4EDM4hep2LcioConv::detail::mapLookupFrom(th, trackerHitPlaneMap)) {
        lcio_tr->addHit(hit.value());
      }
    }
  }
}

template <typename SimTrHitMapT, typename MCParticleMapT>
void resolveRelationsSimTrackerHits(SimTrHitMapT& simTrHitMap, const MCParticleMapT& mcParticleMap) {
  for (auto& [lcio_strh, edm_strh] : simTrHitMap) {
    auto edm_strh_mcp = edm_strh.getParticle();
    if (edm_strh_mcp.isAvailable()) {
      if (const auto lcio_mcp = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm_strh_mcp, mcParticleMap)) {
        lcio_strh->setMCParticle(lcio_mcp.value());
      }
    }
  }
}

template <typename VertexMapT, typename URecoParticleMapT, typename LURecoParticleMapT>
void resolveRelationsVertices(VertexMapT& vertexMap, URecoParticleMapT& updateRPMap,
                              const LURecoParticleMapT& lookupRPMap) {
  // "Invert" the relation to accomodate the different conventions
  for (const auto& [lcio_reco, edm_reco] : updateRPMap) {
    const auto decayVtx = edm_reco.getDecayVertex();
    if (!decayVtx.isAvailable()) {
      continue;
    }
    if (const auto lcio_vtx = k4EDM4hep2LcioConv::detail::mapLookupFrom(decayVtx, vertexMap)) {
      lcio_vtx.value()->setAssociatedParticle(lcio_reco);
    }
    for (const auto& edm_p : decayVtx.getParticles()) {
      if (const auto lcio_p = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm_p, lookupRPMap)) {
        lcio_reco->addParticle(lcio_p.value());
      }
    }
  }

  for (const auto& [lcio_vertex, edm_vertex] : vertexMap) {
    for (const auto& particle : edm_vertex.getParticles()) {
      if (const auto lcio_part = k4EDM4hep2LcioConv::detail::mapLookupFrom(particle, updateRPMap)) {
        lcio_part.value()->setStartVertex(lcio_vertex);
      }
    }
  }
}

template <typename RecoParticleMapT, typename RecoParticleLookupMapT, typename ClusterMapT, typename TrackMapT>
void resolveRelationsRecoParticles(RecoParticleMapT& recoParticleMap, const RecoParticleLookupMapT& recoLookupMap,
                                   const ClusterMapT& clusterMap, const TrackMapT& trackMap) {
  for (auto& [lcio_rp, edm_rp] : recoParticleMap) {
    // Link Tracks
    const auto edmTracks = edm_rp.getTracks();
    for (const auto& t : edmTracks) {
      if (t.isAvailable()) {
        if (const auto lcio_tr = k4EDM4hep2LcioConv::detail::mapLookupFrom(t, trackMap)) {
          lcio_rp->addTrack(lcio_tr.value());
        }
      }
    }

    // Link Clusters
    const auto edmClusters = edm_rp.getClusters();
    for (const auto& c : edmClusters) {
      if (c.isAvailable()) {
        if (const auto lcio_cluster = k4EDM4hep2LcioConv::detail::mapLookupFrom(c, clusterMap)) {
          lcio_rp->addCluster(lcio_cluster.value());
        }
      }
    }

    // Link particles
    const auto edmParticles = edm_rp.getParticles();
    for (const auto& p : edmParticles) {
      if (p.isAvailable()) {
        if (const auto lcio_p = k4EDM4hep2LcioConv::detail::mapLookupFrom(p, recoLookupMap)) {
          lcio_rp->addParticle(lcio_p.value());
        }
      }
    }
  }
}

template <typename ClusterMapT, typename CaloHitMapT>
void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap) {
  // Resolve relations for clusters
  for (auto& [lcio_cluster, edm_cluster] : clustersMap) {
    for (const auto& edm_linked_cluster : edm_cluster.getClusters()) {
      if (edm_linked_cluster.isAvailable()) {
        if (const auto lcio_cluster_linked =
                k4EDM4hep2LcioConv::detail::mapLookupFrom(edm_linked_cluster, clustersMap)) {
          lcio_cluster->addCluster(lcio_cluster_linked.value());
        }
      }
    }

    for (const auto& edm_calohit : edm_cluster.getHits()) {
      if (edm_calohit.isAvailable()) {
        if (const auto lcio_calohit = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm_calohit, caloHitMap)) {
          lcio_cluster->addHit(lcio_calohit.value(), 1.0);
        }
      }
    }
  }
}

template <typename SimCaloHitMapT, typename MCParticleMapT>
void resolveRelationsSimCaloHit(SimCaloHitMapT& simCaloHitMap, const MCParticleMapT& mcParticleMap) {
  // Fill SimCaloHit collections with contributions
  //
  // We loop over all pairs of lcio and edm4hep simcalo hits and add the
  // contributions, by now MCParticle collection(s) should be converted!
  for (auto& [lcio_sch, edm_sch] : simCaloHitMap) {
    const auto contribs = edm_sch.getContributions();
    for (const auto& contrib : contribs) {
      if (not contrib.isAvailable()) {
        // We need a logging library independent of Gaudi for this!
        // std::cout << "WARNING: CaloHit Contribution is not available!" <<
        // std::endl;
        continue;
      }
      auto edm_contrib_mcp = contrib.getParticle();
      std::array<float, 3> step_position{contrib.getStepPosition()[0], contrib.getStepPosition()[1],
                                         contrib.getStepPosition()[2]};
      EVENT::MCParticle* lcio_mcp = nullptr;
      if (edm_contrib_mcp.isAvailable()) {
        // if we have the MCParticle we look for its partner
        lcio_mcp = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm_contrib_mcp, mcParticleMap).value_or(nullptr);
      }
      // add associated Contributions (MCParticles)
      // we add contribution with whatever lcio mc particle we found
      lcio_sch->addMCParticleContribution(lcio_mcp, contrib.getEnergy(), contrib.getTime(), contrib.getPDG(),
                                          step_position.data());
    }
    // We need to reset the energy to the original one, because adding
    // contributions alters the energy in LCIO
    lcio_sch->setEnergy(edm_sch.getEnergy());
  }
}

template <typename PidMapT, typename RecoParticleMapT>
void resolveRelationsParticleIDs(PidMapT& pidMap, const RecoParticleMapT& recoMap) {
  for (auto& [lcioPid, edmPid] : pidMap) {
    const auto edmReco = edmPid.getParticle();
    const auto lcioReco = k4EDM4hep2LcioConv::detail::mapLookupFrom(edmReco, recoMap);
    if (lcioReco) {
      lcioReco.value()->addParticleID(lcioPid);
    } else {
      std::cerr << "Cannot find a reconstructed particle to attach a ParticleID to" << std::endl;
    }
  }
}

template <typename TrackMapT>
void attachDqdxInfo(TrackMapT& trackMap, const std::vector<TrackDqdxConvData>& dQdxCollections) {
  for (const auto& coll : dQdxCollections) {
    attachDqdxInfo(trackMap, coll);
  }
}

template <typename TrackMapT>
void attachDqdxInfo(TrackMapT& trackMap, const TrackDqdxConvData& dQdxCollection) {
  const auto& [name, coll] = dQdxCollection;

  for (const auto& elem : *coll) {
    const auto dQdxTrack = elem.getTrack();
    const auto lcioTrack = k4EDM4hep2LcioConv::detail::mapLookupFrom(dQdxTrack, trackMap);
    if (lcioTrack) {
      lcioTrack.value()->setdEdx(elem.getDQdx().value);
      lcioTrack.value()->setdEdxError(elem.getDQdx().error);
    } else {
      std::cerr << "Cannot find a track to attach dQ/dx information to for collection: " << name << std::endl;
    }
  }
}

template <typename ObjectMappingT>
void resolveRelations(ObjectMappingT& collection_pairs) {
  resolveRelations(collection_pairs, collection_pairs);
}

// Depending on the order of the collections in the parameters,
// and for the mutual dependencies between some collections,
// go over the possible missing associated collections and fill them.
template <typename ObjectMappingT, typename ObjectMappingU>
void resolveRelations(ObjectMappingT& update_pairs, const ObjectMappingU& lookup_pairs) {
  resolveRelationsMCParticles(update_pairs.mcParticles, lookup_pairs.mcParticles);
  resolveRelationsTracks(update_pairs.tracks, lookup_pairs.trackerHits, lookup_pairs.trackerHitPlanes,
                         lookup_pairs.tpcHits);
  resolveRelationsRecoParticles(update_pairs.recoParticles, lookup_pairs.recoParticles, lookup_pairs.clusters,
                                lookup_pairs.tracks);
  resolveRelationsParticleIDs(lookup_pairs.particleIDs, update_pairs.recoParticles);
  resolveRelationsVertices(update_pairs.vertices, update_pairs.recoParticles, lookup_pairs.recoParticles);
  resolveRelationsSimCaloHit(update_pairs.simCaloHits, lookup_pairs.mcParticles);
  resolveRelationsSimTrackerHits(update_pairs.simTrackerHits, lookup_pairs.mcParticles);
  resolveRelationsClusters(update_pairs.clusters, lookup_pairs.caloHits);
}

template <typename ObjectMappingT>
std::vector<std::tuple<std::string, std::unique_ptr<lcio::LCCollection>>>
createLCRelationCollections(const std::vector<std::tuple<std::string, const podio::CollectionBase*>>& linkCollections,
                            const ObjectMappingT& objectMaps) {
  std::vector<std::tuple<std::string, std::unique_ptr<lcio::LCCollection>>> relationColls{};
  relationColls.reserve(linkCollections.size());

  for (const auto& [name, coll] : linkCollections) {
    if (const auto links = dynamic_cast<const edm4hep::RecoMCParticleLinkCollection*>(coll)) {
      relationColls.emplace_back(name,
                                 createLCRelationCollection(*links, objectMaps.recoParticles, objectMaps.mcParticles));
    } else if (const auto links = dynamic_cast<const edm4hep::CaloHitSimCaloHitLinkCollection*>(coll)) {
      relationColls.emplace_back(name, createLCRelationCollection(*links, objectMaps.caloHits, objectMaps.simCaloHits));
    } else if (const auto links = dynamic_cast<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>(coll)) {
      relationColls.emplace_back(name,
                                 createLCRelationCollection(*links, objectMaps.trackerHits, objectMaps.simTrackerHits));
    } else if (const auto links = dynamic_cast<const edm4hep::CaloHitMCParticleLinkCollection*>(coll)) {
      relationColls.emplace_back(name, createLCRelationCollection(*links, objectMaps.caloHits, objectMaps.mcParticles));
    } else if (const auto links = dynamic_cast<const edm4hep::ClusterMCParticleLinkCollection*>(coll)) {
      relationColls.emplace_back(name, createLCRelationCollection(*links, objectMaps.clusters, objectMaps.mcParticles));
    } else if (const auto links = dynamic_cast<const edm4hep::TrackMCParticleLinkCollection*>(coll)) {
      relationColls.emplace_back(name, createLCRelationCollection(*links, objectMaps.tracks, objectMaps.mcParticles));
    } else if (const auto links = dynamic_cast<const edm4hep::VertexRecoParticleLinkCollection*>(coll)) {
      relationColls.emplace_back(name,
                                 createLCRelationCollection(*links, objectMaps.vertices, objectMaps.recoParticles));
    } else {
      std::cerr << "Trying to create an LCRelation collection from a " << coll->getTypeName()
                << " which is not supported" << std::endl;
    }
  }

  return relationColls;
}

namespace detail {
  template <typename T>
  consteval const char* getTypeName();

#define DEFINE_TYPE_NAME(type)                                                                                         \
  template <>                                                                                                          \
  consteval const char* getTypeName<IMPL::type##Impl>() {                                                              \
    return #type;                                                                                                      \
  }                                                                                                                    \
  template <>                                                                                                          \
  consteval const char* getTypeName<EVENT::type>() {                                                                   \
    return #type;                                                                                                      \
  }

  DEFINE_TYPE_NAME(MCParticle);
  DEFINE_TYPE_NAME(SimTrackerHit);
  DEFINE_TYPE_NAME(SimCalorimeterHit);
  DEFINE_TYPE_NAME(Track);
  DEFINE_TYPE_NAME(TrackerHit);
  DEFINE_TYPE_NAME(Vertex);
  DEFINE_TYPE_NAME(ReconstructedParticle);
  DEFINE_TYPE_NAME(Cluster);
  DEFINE_TYPE_NAME(CalorimeterHit);

#undef DEFINE_TYPE_NAME
} // namespace detail

template <typename LinkCollT, typename FromMapT, typename ToMapT>
std::unique_ptr<lcio::LCCollection> createLCRelationCollection(const LinkCollT& links, const FromMapT& fromMap,
                                                               const ToMapT& toMap) {
  using FromLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<FromMapT>>;
  using ToLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<ToMapT>>;

  auto lcioColl = std::make_unique<lcio::LCCollectionVec>(lcio::LCIO::LCRELATION);
  lcioColl->parameters().setValue("FromType", detail::getTypeName<FromLCIOT>());
  lcioColl->parameters().setValue("ToType", detail::getTypeName<ToLCIOT>());

  for (const auto link : links) {
    auto lcioRel = new lcio::LCRelationImpl{};
    lcioRel->setWeight(link.getWeight());

    const auto edm4hepFrom = link.getFrom();
    const auto lcioFrom = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm4hepFrom, fromMap);
    if (lcioFrom) {
      lcioRel->setFrom(lcioFrom.value());
    } else {
      std::cerr << "Cannot find an object for building an LCRelation of type " << detail::getTypeName<FromLCIOT>()
                << std::endl;
    }

    const auto edm4hepTo = link.getTo();
    const auto lcioTo = k4EDM4hep2LcioConv::detail::mapLookupFrom(edm4hepTo, toMap);
    if (lcioTo) {
      lcioRel->setTo(lcioTo.value());
    } else {
      std::cerr << "Cannot find an objects for building an LCRelation of type " << detail::getTypeName<ToLCIOT>()
                << std::endl;
    }

    lcioColl->addElement(lcioRel);
  }

  return lcioColl;
}

} // namespace EDM4hep2LCIOConv
