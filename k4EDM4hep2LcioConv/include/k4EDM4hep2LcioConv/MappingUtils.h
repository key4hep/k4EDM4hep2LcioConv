#ifndef K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H
#define K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H

#include <optional>
#include <algorithm>
#include <vector>
#include <tuple>
#include <unordered_map>

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

namespace k4EDM4hep2LcioConv {

  namespace detail {
    /// Very minimal detector for whether T is a map or not. We simply assume that
    /// if it has a key_type we can also call find on T.
    template<typename T>
    using has_key_t = typename T::key_type;

    template<typename T>
    constexpr static bool is_map_v = det::is_detected_v<has_key_t, T>;

    /// Helper struct to determine the key and mapped types for map-like types or
    /// maps
    template<typename T, typename IsMap = std::bool_constant<is_map_v<T>>>
    struct map_t_helper {
    };

    template<typename T>
    struct map_t_helper<T, std::bool_constant<true>> {
      using key_type = typename T::key_type;
      using mapped_type = typename T::mapped_type;
    };

    template<typename T>
    struct map_t_helper<T, std::bool_constant<false>> {
      using key_type = typename std::tuple_element<0, typename T::value_type>::type;
      using mapped_type = typename std::tuple_element<1, typename T::value_type>::type;
    };

    template<typename T>
    using key_t = typename map_t_helper<T>::key_type;

    template<typename T>
    using mapped_t = typename map_t_helper<T>::mapped_type;

    /**
     * Find the mapped-to object in a map provided a key object
     *
     * NOTE: This will use a potentially more efficient lookup for actual map
     * types (i.e. MapT::find). In that case it will have the time complexity of
     * that. In case of a "map-like" (e.g. vector<tuple<K, V>>) it will be O(N).
     */
    template<typename FromT, typename MapT, typename = std::enable_if_t<std::is_same_v<FromT, key_t<MapT>>>>
    auto mapLookupTo(FromT keyObj, const MapT& map) -> std::optional<mapped_t<MapT>>
    {
      if constexpr (is_map_v<MapT>) {
        if (const auto& it = map.find(keyObj); it != map.end()) {
          return it->second;
        }
      }
      else {
        if (const auto& it = std::find_if(
              map.begin(), map.end(), [&keyObj](const auto& mapElem) { return std::get<0>(mapElem) == keyObj; });
            it != map.end()) {
          return std::get<1>(*it);
        }
      }

      return std::nullopt;
    }

    /**
     * Find the mapped-from (or key object) in a "map" provided a mapped-to object
     *
     * NOTE: This will always loop over potentially all elements in the provided
     * map, so it is definitely O(N) regardless of the provided map type
     */
    template<typename ToT, typename MapT, typename = std::enable_if_t<std::is_same_v<ToT, mapped_t<MapT>>>>
    auto mapLookupFrom(ToT mappedObj, const MapT& map) -> std::optional<key_t<MapT>>
    {
      // In this case we cannot use a potential find method for an actual map, but
      // looping over the map and doing the actual comparison will work
      if (const auto& it = std::find_if(
            map.begin(), map.end(), [&mappedObj](const auto& mapElem) { return std::get<1>(mapElem) == mappedObj; });
          it != map.end()) {
        return std::get<0>(*it);
      }

      return std::nullopt;
    }

    enum class InsertMode { Unchecked, Checked };

    /**
     * Insert a key-value pair into a "map"
     *
     * The InsertMode argument can be use to check whether the Key already
     * exists in the map before inserting. This is only useful for maps using a
     * vector as backing, since the usual emplace already does this check and
     * does not insert if a key already exists
     */
    template<typename MapT, typename KeyT = key_t<MapT>, typename MappedT = mapped_t<MapT>>
    auto mapInsert(KeyT&& key, MappedT&& mapped, MapT& map, InsertMode insertMode = InsertMode::Unchecked)
    {
      if constexpr (is_map_v<MapT>) {
        return map.emplace(std::forward<KeyT>(key), std::forward<MappedT>(mapped));
      }
      else {
        if (insertMode == InsertMode::Checked) {
          if (auto existing = mapLookupTo(key, map)) {
            return std::make_pair(std::make_tuple(key, existing.value()), false);
          }
        }
        return std::make_pair(map.emplace_back(std::forward<KeyT>(key), std::forward<MappedT>(mapped)), true);
      }
    }

    /// Helper type alias that can be used to detect whether a T can be used
    /// with std::get directly or whether it has to be dereferenced first
    template<typename T>
    using std_get_usable = decltype(std::get<0>(std::declval<T>()));

    /**
     * Helper function to get the Key from an Iterator (e.g. returned by
     * mapInsert). This is necessary because map::emplace returns an iterator,
     * while vector::emplace_back returns a reference to the inserted element.
     * Simply providing two overloads here does the trick.
     */
    template<typename It>
    auto getKey(const It& it)
    {
      if constexpr (det::is_detected_v<std_get_usable, It>) {
        return std::get<0>(it);
      }
      else {
        return std::get<0>(*it);
      }
    }

    /**
     * Helper function to get the Value from an Iterator (e.g. returned by
     * mapInsert). This is necessary because map::emplace returns an iterator,
     * while vector::emplace_back returns a reference to the inserted element.
     * Simply providing two overloads here does the trick.
     */
    template<typename It>
    auto getMapped(const It& it)
    {
      if constexpr (det::is_detected_v<std_get_usable, It>) {
        return std::get<1>(it);
      }
      else {
        return std::get<1>(*it);
      }
    }
  } // namespace detail

  template<typename K, typename V>
  using MapT = std::unordered_map<K, V>;

  template<typename K, typename V>
  using VecMapT = std::vector<std::tuple<K, V>>;

} // namespace k4EDM4hep2LcioConv

#endif // K4EDM4HEP2LCIOCONV_MAPPINGUTILS_H
