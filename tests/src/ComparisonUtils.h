#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H

#include "ObjectMapping.h"

#include "edm4hep/Vector2f.h"
#include "edm4hep/Vector2i.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector3f.h"
#include <edm4hep/CovMatrix2f.h>
#include <edm4hep/CovMatrix3f.h>
#include <edm4hep/CovMatrix4f.h>
#include <edm4hep/CovMatrix6f.h>

#include "EVENT/LCCollection.h"
#include "UTIL/LCIterator.h"

#include "podio/RelationRange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

template <typename T>
inline std::ostream& printContainer(std::ostream& os, const T& cont) {
  os << "(";
  if (!cont.empty()) {
    os << cont[0];
    for (size_t i = 1; i < cont.size(); ++i) {
      os << ", " << cont[i];
    }
  }

  return os << ")";
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  return printContainer(os, vec);
}

template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
  return printContainer(os, arr);
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const podio::RelationRange<T>& range) {
  return printContainer(os, range);
}

template <typename T, typename U>
inline bool nanSafeComp(T x, U y) {
  return (x == y) || (std::isnan(x) && std::isnan(y));
}

// Only enable for vectors
template <typename T, typename = void>
struct has_size_method : std::false_type {};

template <typename T>
struct has_size_method<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};

template <typename T, typename... Ts>
constexpr bool isAnyOf = (std::is_same_v<T, Ts> || ...);

// Helper function for comparing values and vectors of values element by
// element, ignoring cases where both values are nan
template <typename LCIO, typename EDM4hepT>
bool compareValuesNanSafe(LCIO lcioV, EDM4hepT edm4hepV, const std::string& msg) {
  constexpr auto isVectorLike =
      has_size_method<EDM4hepT>::value ||
      isAnyOf<EDM4hepT, edm4hep::Vector3f, edm4hep::Vector3d, edm4hep::Vector2f, edm4hep::Vector2i,
              // These also effectively behave like vectors for
              // the purposes of this function
              edm4hep::CovMatrix2f, edm4hep::CovMatrix3f, edm4hep::CovMatrix4f, edm4hep::CovMatrix6f>;

  if constexpr (isVectorLike) {
    const auto vecSize = [](EDM4hepT& edm4hepVal) -> std::size_t {
      if constexpr (has_size_method<EDM4hepT>::value) {
        return edm4hepVal.size();
      } else if constexpr (std::is_same_v<EDM4hepT, edm4hep::Vector3f> || std::is_same_v<EDM4hepT, edm4hep::Vector3d>) {
        return 3u;
      } else if constexpr (std::is_same_v<EDM4hepT, edm4hep::Vector2i> || std::is_same_v<EDM4hepT, edm4hep::Vector2f>) {
        return 2;
      }
      return 0;
    }(edm4hepV);
    for (size_t i = 0; i < vecSize; ++i) {
      if (!nanSafeComp(lcioV[i], edm4hepV[i])) {
        std::cerr << msg << " at index " << i << ": (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")"
                  << std::endl;
        return false;
      }
    }
  } else {
    if (!nanSafeComp(lcioV, edm4hepV)) {
      std::cerr << msg << " (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;
      return false;
    }
  }

  return true;
}

// Macro for comparing the return types of the different functions and return
// false if they are not equal while also emitting a message
#define ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg)                                                                      \
  if (!compareValuesNanSafe(lcioV, edm4hepV, msg)) {                                                                   \
    return false;                                                                                                      \
  }

#define ASSERT_COMPARE_VALS_FLOAT(lcioV, edm4hepV, tol, msg)                                                           \
  if (std::abs(lcioV - edm4hepV) > tol) {                                                                              \
    std::cerr << msg << " (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;                      \
    return false;                                                                                                      \
  }

#define ASSERT_COMPARE(lcioE, edm4hepE, func, msg)                                                                     \
  {                                                                                                                    \
    const auto lcioV = lcioE->func();                                                                                  \
    const auto edm4hepV = edm4hepE.func();                                                                             \
    ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg)                                                                          \
  }

/**
 * Compare a single relation by checking whether the LCIO object points to the
 * correct EDM4hep element (using the ObjectIDs)
 */
template <typename LcioT, typename EDM4hepT, typename MapT>
inline bool compareRelation(const LcioT* lcioElem, const EDM4hepT& edm4hepElem, const MapT& objectMap,
                            const std::string& msg) {
  if (lcioElem == nullptr && !edm4hepElem.isAvailable()) {
    // Both elements are "empty". Nothing more to do here
    return true;
  }
  if ((lcioElem == nullptr && edm4hepElem.isAvailable()) || (lcioElem != nullptr && !edm4hepElem.isAvailable())) {
    const auto emptyOrNot = [](const bool b) { return b ? "not empty" : "empty"; };
    std::cerr << msg << " LCIO element is " << emptyOrNot(lcioElem) << " but edm4hep element is "
              << emptyOrNot(edm4hepElem.isAvailable()) << std::endl;
    return false;
  }

  // Now we know for sure that
  if (const auto it = objectMap.find(lcioElem); it != objectMap.end()) {
    if (!(it->second == edm4hepElem.getObjectID())) {
      std::cerr << msg << " LCIO element " << lcioElem << " points to " << it->second << " but should point to "
                << edm4hepElem.getObjectID() << std::endl;
      return false;
    }
  } else {
    std::cerr << msg << " cannot find LCIO object " << lcioElem << " in object map for relation checking" << std::endl;
    return false;
  }

  return true;
}

/**
 * Compare the relations in the form of a range of LCIO objects to the
 * corresponding range of EDM4hep objects. This uses the LCIO object for lookup
 * in the object map and then compares the object ID of the EDM4hep object.
 *
 * Naming is using the singular form to have an overload set in the macro below
 * that dispatches this.
 */
template <typename LcioT, typename EDM4hepT, typename MapT>
inline bool compareRelation(const std::vector<LcioT*>& lcioRange, const podio::RelationRange<EDM4hepT>& edm4hepRange,
                            const MapT& objectMap, const std::string& msg) {
  if (lcioRange.size() != edm4hepRange.size()) {
    // Make sure to take into account that the LCIO -> EDM4hep conversion does
    // not fill relations if the original relations in LCIO were empty
    const auto nonNullLcio =
        std::count_if(lcioRange.begin(), lcioRange.end(), [](const auto e) { return e != nullptr; });
    if ((unsigned)nonNullLcio != edm4hepRange.size()) {
      std::cerr << msg
                << " different sizes (even after taking null values into account): (expected: " << edm4hepRange.size()
                << " actual: " << lcioRange.size() << " | " << nonNullLcio << " cleaned)" << std::endl;
      return false;
    }
  }

  for (size_t i = 0; i < edm4hepRange.size(); ++i) {
    const auto lcioElem = lcioRange[i];
    const auto edm4hepElem = edm4hepRange[i];
    if (!compareRelation(lcioElem, edm4hepElem, objectMap, msg)) {
      return false;
    }
  }

  return true;
}

#define ASSERT_COMPARE_RELATION(lcioE, edm4hepE, func, map, msg)                                                       \
  {                                                                                                                    \
    const auto& lcioRel = lcioE->func();                                                                               \
    const auto edm4hepRel = edm4hepE.func();                                                                           \
    if (!compareRelation(lcioRel, edm4hepRel, map, msg)) {                                                             \
      return false;                                                                                                    \
    }                                                                                                                  \
  }

// Compare an LCIO collection and an EDM4hep collection. Assumes that a compare
// function working with the element types is available
template <typename LCIOT, typename EDM4hepCollT>
inline bool compareCollection(const lcio::LCCollection* lcioCollection, const EDM4hepCollT& edm4hepCollection,
                              const ObjectMappings& objectMaps) {
  UTIL::LCIterator<LCIOT> lcioIt(lcioCollection);
  int counter = 0;
  for (const auto edm4hepElem : edm4hepCollection) {
    if (!compare(lcioIt.next(), edm4hepElem, objectMaps)) {
      std::cerr << "in Element " << counter << std::endl;
      return false;
    }
    counter++;
  }

  return true;
}

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
