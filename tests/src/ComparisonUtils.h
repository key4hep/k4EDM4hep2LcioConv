#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H

#include "edm4hep/Vector2f.h"
#include "edm4hep/Vector2i.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector3f.h"

#include "UTIL/LCIterator.h"
#include "EVENT/LCCollection.h"

#include <cmath>
#include <array>
#include <iostream>
#include <vector>

template<typename T>
std::ostream& printContainer(std::ostream& os, const T& cont)
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
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  return printContainer(os, vec);
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
  return printContainer(os, arr);
}

template<typename T, size_t N>
bool operator==(const std::vector<T>& vec, const std::array<T, N>& arr)
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
bool operator!=(const std::vector<T>& vec, const std::array<T, N>& arr)
{
  return !(vec == arr);
}

template<typename T, typename U>
bool nanSafeComp(T x,U y){
  return (x == y) || (std::isnan(x) && std::isnan(y));
}

// Macro for defining the comparison operators for edm4hep::Vector3X and
// different return types (X* or vector<X> from LCIO)
#define VECTOR3_COMPARE(FT, VT)                                             \
  bool operator==(const FT* vals, const VT& vec)                            \
  {                                                                         \
    return nanSafeComp(vals[0],vec[0]) && nanSafeComp(vals[1],vec[1]) && nanSafeComp(vals[2],vec[2]); \
  }                                                                         \
  bool operator!=(const FT* vals, const VT& vec) { return !(vals == vec); } \
  bool operator==(const std::vector<FT>& vals, const VT& vec)               \
  {                                                                         \
    if (vals.size() != 3) {                                                 \
      return false;                                                         \
    }                                                                       \
    return vals.data() == vec;                                              \
  }                                                                         \
  bool operator!=(const std::vector<FT>& vals, const VT& vec) { return !(vals == vec); }

VECTOR3_COMPARE(float, edm4hep::Vector3f)
VECTOR3_COMPARE(double, edm4hep::Vector3d)
// Necessary in some MCParticle return types
VECTOR3_COMPARE(double, edm4hep::Vector3f)
#undef VECTOR3_COMPARE

// Macro for defining the comparison operators for edm4hep::Vector3X and
// different return types (X* or vector<X> from LCIO)
#define VECTOR2_COMPARE(FT, VT)                                                                     \
  bool operator==(const FT* vals, const VT& vec) { return vals[0] == vec[0] && vals[1] == vec[1]; } \
  bool operator!=(const FT* vals, const VT& vec) { return !(vals == vec); }                         \
  bool operator==(const std::vector<FT>& vals, const VT& vec)                                       \
  {                                                                                                 \
    if (vals.size() != 2) {                                                                         \
      return false;                                                                                 \
    }                                                                                               \
    return vals.data() == vec;                                                                      \
  }                                                                                                 \
  bool operator!=(const std::vector<FT>& vals, const VT& vec) { return !(vals == vec); }

VECTOR2_COMPARE(int, edm4hep::Vector2i)
VECTOR2_COMPARE(float, edm4hep::Vector2f)
#undef VECTOR2_COMPARE

// Macro for comparing the return types of the different functions and return
// false if they are not equal while also emitting a message
#define ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg)                                                 \
  if ((lcioV) != (edm4hepV)) {                                                                    \
    std::cerr << msg << " (LCIO: " << (lcioV) << ", EDM4hep: " << (edm4hepV) << ")" << std::endl; \
    return false;                                                                                 \
  }

#define ASSERT_COMPARE(lcioE, edm4hepE, func, msg) \
  {                                                \
    const auto lcioV = lcioE->func();              \
    const auto edm4hepV = edm4hepE.func();         \
    ASSERT_COMPARE_VALS(lcioV, edm4hepV, msg)      \
  }

// Compare an LCIO collection and an EDM4hep collection. Assumes that a compare
// function working with the element types is available
template<typename LCIOT, typename EDM4hepCollT>
bool compareCollection(const lcio::LCCollection* lcioCollection, const EDM4hepCollT& edm4hepCollection)
{
  if (lcioCollection->getNumberOfElements() != edm4hepCollection.size()) {
    return false;
  }

  UTIL::LCIterator<LCIOT> lcioIt(lcioCollection);
  int counter = 0;
  for (const auto edm4hepElem : edm4hepCollection) {
    if (!compare(lcioIt.next(), edm4hepElem)) {
      std::cerr << "in Element " << counter << std::endl;
      return false;
    }
    counter++;
  }

  return true;
}

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPARISONUTILS_H
