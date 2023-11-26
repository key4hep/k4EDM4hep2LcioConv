#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H

#include "ObjectMapping.h"

#include "edm4hep/Vector2f.h"
#include "edm4hep/Vector2i.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector3f.h"

#include "UTIL/LCIterator.h"
#include "EVENT/LCCollection.h"

#include "podio/RelationRange.h"

#include <cmath>
#include <array>
#include <iostream>
#include <vector>
#include <algorithm>

template<typename T>
inline std::ostream& printContainer(std::ostream& os, const T& cont)
{
  os << "(";
  if (!cont.empty()) {
    os << cont[0];
    for (size_t i = 1; i < cont.size(); ++i) {
      os << ", " << cont[i];
    }
  }

  return os << ")";
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  return printContainer(os, vec);
}

template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
  return printContainer(os, arr);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const podio::RelationRange<T>& range)
{
  return printContainer(os, range);
}

template<typename T, size_t N>
inline bool operator==(const std::vector<T>& vec, const std::array<T, N>& arr)
{
  if (vec.size() != N) {
    return false;
  }
  for (size_t i = 0; i < N; ++i) {
    if (vec[i] != arr[i]) {
      return false;
    }
  }
  return true;
}

template<typename T, size_t N>
inline bool operator!=(const std::vector<T>& vec, const std::array<T, N>& arr)
{
  return !(vec == arr);
}

template<typename T>
inline bool operator==(const std::vector<T>& vec, const podio::RelationRange<T>& range)
{
  if (vec.size() != range.size()) {
    return false;
  }
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i] != range[i]) {
      return false;
    }
  }
  return true;
}

template<typename T>
inline bool operator!=(const std::vector<T>& vec, const podio::RelationRange<T>& range)
{
  return !(vec == range);
}

template<typename T, typename U>
inline bool nanSafeComp(T x, U y)
{
  return (x == y) || (std::isnan(x) && std::isnan(y));
}

// Macro for defining the comparison operators for edm4hep::Vector3X and
// different return types (X* or vector<X> from LCIO)
#define VECTOR3_COMPARE(FT, VT)                                                                          \
  inline bool operator==(const FT* vals, const VT& vec)                                                  \
  {                                                                                                      \
    return nanSafeComp(vals[0], vec[0]) && nanSafeComp(vals[1], vec[1]) && nanSafeComp(vals[2], vec[2]); \
  }                                                                                                      \
  inline bool operator!=(const FT* vals, const VT& vec) { return !(vals == vec); }                       \
  inline bool operator==(const std::vector<FT>& vals, const VT& vec)                                     \
  {                                                                                                      \
    if (vals.size() != 3) {                                                                              \
      return false;                                                                                      \
    }                                                                                                    \
    return vals.data() == vec;                                                                           \
  }                                                                                                      \
  inline bool operator!=(const std::vector<FT>& vals, const VT& vec) { return !(vals == vec); }

VECTOR3_COMPARE(float, edm4hep::Vector3f)
VECTOR3_COMPARE(double, edm4hep::Vector3d)
// Necessary in some MCParticle return types
VECTOR3_COMPARE(double, edm4hep::Vector3f)
#undef VECTOR3_COMPARE

// Macro for defining the comparison operators for edm4hep::Vector3X and
// different return types (X* or vector<X> from LCIO)
#define VECTOR2_COMPARE(FT, VT)                                                                            \
  inline bool operator==(const FT* vals, const VT& vec) { return vals[0] == vec[0] && vals[1] == vec[1]; } \
  inline bool operator!=(const FT* vals, const VT& vec) { return !(vals == vec); }                         \
  inline bool operator==(const std::vector<FT>& vals, const VT& vec)                                       \
  {                                                                                                        \
    if (vals.size() != 2) {                                                                                \
      return false;                                                                                        \
    }                                                                                                      \
    return vals.data() == vec;                                                                             \
  }                                                                                                        \
  inline bool operator!=(const std::vector<FT>& vals, const VT& vec) { return !(vals == vec); }

VECTOR2_COMPARE(int, edm4hep::Vector2i)
VECTOR2_COMPARE(float, edm4hep::Vector2f)
#undef VECTOR2_COMPARE

// Macro for comparing the return types of the different functions and return
// false if they are not equal while also emitting a message

// Only enable for vectors
template<typename LCIO, typename EDM4hep, typename = void>
struct has_size_method : std::false_type {};

template<typename LCIO, typename EDM4hep>
struct has_size_method<LCIO, EDM4hep, std::void_t<decltype(std::declval<EDM4hep>().size())>> : std::true_type {};

template<typename LCIO, typename EDM4hep>
bool ASSERT_COMPARE_VALS_TEMPLATE(LCIO lcioV, EDM4hep edm4hepV, const std::string& msg)
{
  if constexpr (has_size_method<LCIO, EDM4hep>::value) {
    for (size_t i = 0; i < (edm4hepV).size(); ++i) {
      if ((std::isnan((lcioV)[i]) != std::isnan((edm4hepV)[i])) || ((lcioV)[i] != (edm4hepV)[i])) {
        std::cerr << msg << " at index i: (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;
        return false;
      }
    }
  }
  else if constexpr ((std::is_same_v<EDM4hep, edm4hep::Vector3f>) || std::is_same_v<EDM4hep, edm4hep::Vector3d>) {
    for (size_t i = 0; i < 3; ++i) {
      if ((std::isnan((lcioV)[i]) != std::isnan((edm4hepV)[i])) || ((lcioV)[i] != (edm4hepV)[i])) {
        std::cerr << msg << " at index i: (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;
        return false;
      }
    }
  }
  else if constexpr ((std::is_same_v<EDM4hep, edm4hep::Vector2i>) || std::is_same_v<EDM4hep, edm4hep::Vector2f>) {
    for (size_t i = 0; i < 2; ++i) {
      if ((std::isnan((lcioV)[i]) != std::isnan((edm4hepV)[i])) || ((lcioV)[i] != (edm4hepV)[i])) {
        std::cerr << msg << " at index i: (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;
        return false;
      }
    }
  }
  else {
    if ((std::isnan((lcioV)) != std::isnan((edm4hepV))) || ((lcioV) != (edm4hepV))) {
      std::cerr << msg << " (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl;
      return false;
    }
  }
  return true;
}

#define ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg) return ASSERT_COMPARE_VALS_TEMPLATE(lcioV, edm4hepV, msg);

#define ASSERT_COMPARE(lcioE, edm4hepE, func, msg) \
  {                                                \
    const auto lcioV = lcioE->func();              \
    const auto edm4hepV = edm4hepE.func();         \
    ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg)      \
  }

/**
 * Compare a single relation by checking whether the LCIO object points to the
 * correct EDM4hep element (using the ObjectIDs)
 */
template<typename LcioT, typename EDM4hepT, typename MapT>
inline bool
compareRelation(const LcioT* lcioElem, const EDM4hepT& edm4hepElem, const MapT& objectMap, const std::string& msg)
{
  if (lcioElem == nullptr && edm4hepElem.isAvailable()) {
    std::cerr << msg << " LCIO element is empty but edm4hep element is not" << std::endl;
  }
  if (const auto it = objectMap.find(lcioElem); it != objectMap.end()) {
    if (!(it->second == edm4hepElem.getObjectID())) {
      std::cerr << msg << " LCIO element " << lcioElem << " points to " << it->second << " but should point to "
                << edm4hepElem.getObjectID() << std::endl;
      return false;
    }
  }
  else {
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
template<typename LcioT, typename EDM4hepT, typename MapT>
inline bool compareRelation(
  const std::vector<LcioT*>& lcioRange,
  const podio::RelationRange<EDM4hepT>& edm4hepRange,
  const MapT& objectMap,
  const std::string& msg)
{
  if (lcioRange.size() != edm4hepRange.size()) {
    // Make sure to take into account that the LCIO -> EDM4hep conversion does
    // not fill relations if the original relations in LCIO were empty
    const auto nonNullLcio =
      std::count_if(lcioRange.begin(), lcioRange.end(), [](const auto e) { return e != nullptr; });
    if (nonNullLcio != edm4hepRange.size()) {
      std::cerr << msg << " different sizes (even after taking null values into account)" << std::endl;
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

#define ASSERT_COMPARE_RELATION(lcioE, edm4hepE, func, map, msg) \
  {                                                              \
    const auto& lcioRel = lcioE->func();                         \
    const auto edm4hepRel = edm4hepE.func();                     \
    return compareRelation(lcioRel, edm4hepRel, map, msg);       \
  }

// Compare an LCIO collection and an EDM4hep collection. Assumes that a compare
// function working with the element types is available
template<typename LCIOT, typename EDM4hepCollT>
inline bool compareCollection(
  const lcio::LCCollection* lcioCollection,
  const EDM4hepCollT& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
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

#define ASSERT_COMPARE_OR_EXIT(collType)                   \
  if (type == #collType) {                                 \
    auto& edmcoll = edmEvent.get<collType>(name);          \
    if (!compare(lcioColl, edmcoll, objectMapping)) {      \
      std::cerr << "in collection: " << name << std::endl; \
      return 1;                                            \
    }                                                      \
  }

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
