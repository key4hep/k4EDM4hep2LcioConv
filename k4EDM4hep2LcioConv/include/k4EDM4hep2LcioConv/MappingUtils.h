#ifndef K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H
#define K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H

#include <optional>
#include <algorithm>

#if __has_include("experimental/type_traits.h")
#include <experimental/type_traits>
namespace det {
  using namespace std::experimental;
}
#else
// Implement the minimal feature set we need
#include <type_traits>

namespace det {
  namespace detail {
    template<typename DefT, typename AlwaysVoidT, template<typename...> typename Op, typename... Args>
    struct detector {
      using value_t = std::false_type;
      using type = DefT;
    };

    template<typename DefT, template<typename...> typename Op, typename... Args>
    struct detector<DefT, std::void_t<Op<Args...>>, Op, Args...> {
      using value_t = std::true_type;
      using type = Op<Args...>;
    };
  } // namespace detail

  struct nonesuch {
    ~nonesuch() = delete;
    nonesuch(const nonesuch&) = delete;
    void operator=(const nonesuch&) = delete;
  };

  template<template<typename...> typename Op, typename... Args>
  using is_detected = typename detail::detector<nonesuch, void, Op, Args...>::value_t;

  template<template<typename...> typename Op, typename... Args>
  static constexpr bool is_detected_v = is_detected<Op, Args...>::value;
} // namespace det
#endif

namespace k4EDM4hep2LCIOConv::detail {
  /// Very minimal detector for whether T is a map or not. We simply assume that
  /// if it has a key_type we can also call find on T.
  template<typename T>
  using has_key_t = typename T::key_type;

  template<typename LcioT, typename EDM4hepT, template<typename...> typename MapT>
  std::optional<EDM4hepT> mapLookup(LcioT lcioObj, const MapT<LcioT, EDM4hepT>& map)
  {
    if constexpr (det::is_detected_v<has_key_t, MapT<LcioT, EDM4hepT>>) {
      if (const auto& it = map.find(lcioObj); it != map.end()) {
        return it->second;
      }
    }
    else {
      if (const auto& it = std::find_if(
            map.begin(), map.end(), [&lcioObj](const auto& mapElem) { return std::get<0>(mapElem) == lcioObj; });
          it != map.end()) {
        return it->second;
      }
    }

    return std::nullopt;
  }

} // namespace k4EDM4hep2LCIOConv::detail

#endif // K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H
