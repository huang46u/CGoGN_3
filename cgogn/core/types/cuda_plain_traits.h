#pragma once

#include <type_traits>

namespace cgogn::cuda
{

// Default trait: PlainType is the same as the value type. The default
// assumes the value is trivially copyable so we can memcpy between host buffers
// and device buffers without additional work.
template <typename T, typename Enable = void>
struct CudaPlainValueTraits
{
    using PlainType = T;

    static PlainType to_plain(const T& value)
    {
        static_assert(std::is_trivially_copyable_v<T>,
                      "CudaPlainValueTraits<T>: provide a specialisation for non-trivial types");
        return value;
    }

    static T from_plain(const PlainType& plain)
    {
        return plain;
    }

    static constexpr bool can_write_back = true;
};

} // namespace cgogn::cuda
