# EDM4hep to LCIO in-memory converter

<p align="center">
  <img src="doc/k4EDM4hep2LcioConv_logo.svg"/>
</p>

Sample usage:

```cpp
// Declare struct to hold the converted collections
CollectionsPairVectors collection_pairs{};

// Pointer to edm4hep collection
edm4hep::ReconstructedParticleCollection* rp_collection;

// Converting a collection of Reconstructed Particles to LCIO
auto* lcio_converted_reconstructed_particle_ptr = convReconstructedParticles(
  rp_collection, // input collection to be converted
  collection_pairs.recoparticles, // vector holding converted and original collections
  collection_pairs.tracks, // Tracks related to Reconstructed Particles to link them
  collection_pairs.vertices, // Vertices related to Reconstructed Particles to link them
  collection_pairs.clusters); // Clusters related to Reconstructed Particles to link them

// Some collections need the cellID string
edm4hep::SimTrackerHitCollection* sth_collection;
// Example from k4FWCore Data Handle
const auto collID    = sth_collection->getID();
const auto cellIDstr = simtrackerhits_handle.getCollMetadataCellID(collID);
auto* lcio_converted_sim_tracker_hit_ptr = convSimTrackerHits(
  sth_collection,
  cellIDstr,
  collection_pairs.simtrackerhits,
  collection_pairs.mcparticles);

// Run function to fix missing links between collections.
// Some collections that need to be linked to other collections may be converted
// after these are linked. Running this function after all conversions guarantees correct links
// between collections.
resolveRelations(collection_pairs);
```
