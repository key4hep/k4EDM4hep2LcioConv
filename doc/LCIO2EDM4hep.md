There exists a convert function for every data collection type (e.g. convertReconstructedParticleconvert). These need to be called before the relations can be handled, since they fill the maps linking the particle in LCIO to their edm4HEP equivalents. Every type has a seperate map. The maps are grouped in a struct for ease of use. They can be defined and stored seperatly.
The order in which the data is converted does not matter because converting data and resolving relations are two differet steps that are carried out in sequence. 
Subset collections are also handled at the same step as relations using the fuction "fillSubset()". Alternatively "handleSubsetColl" can also be called to convert a Subset Collection. This way the unique pointers can be obtained directly.
The OneToMany and OnToOne Relations can be resolved using "resolveRelations()".  There is a resolveRelations function for each type.
The AssociationCollections in EDM4hep are then created using "createAssociations()", the exception here are the CaloHitContributions. Since they are part of the SimCalorimeterHits in LCIO while being a seperate Association in EDM4hep. They are created by "createCaloHitContributions()".
The EventHeader Colletion can be created using "EventHeaderCollection()".

Particle IDs are converted during the conversion of the the reconstructed Particle collection.


Converting an entire event can be done calling the "convertEvent()". 