#ifndef ADHOC_PASSIVE_ID_HPP
#define ADHOC_PASSIVE_ID_HPP

#include <cstddef>
#include <limits>

namespace adhoc {

template<class CounterType>
static constexpr std::size_t passive_id = std::numeric_limits<CounterType>::max();

} // namespace adhoc

#endif // ADHOC_PASSIVE_ID_HPP
