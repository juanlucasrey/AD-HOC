#ifndef ADHOC3_TAPE_BACKPROPAGATOR_HPP
#define ADHOC3_TAPE_BACKPROPAGATOR_HPP

#include <tape_buffer_size.hpp>
#include <tuple>
#include <tuple_utils.hpp>

#include "../tests3/type_name.hpp"

#include <array>
#include <type_traits>

namespace adhoc3::detail {

class available_t {};

template <class NodesDerivativeRest, class InterfaceTypes, class DersActive,
          class CandidateNode>
constexpr auto conditional_node2(NodesDerivativeRest dernodes,
                                 InterfaceTypes dersin,
                                 CandidateNode candnode) {
    if constexpr (contains(dernodes, candnode) || contains(dersin, candnode)) {
        return std::make_tuple(candnode);
    } else {
        return std::tuple<>{};
    }
}

template <std::size_t... I, class... Ids, class NodeValues>
constexpr auto get_op(std::index_sequence<I...> /* prev */,
                      std::tuple<Ids...> /* in */, NodeValues nodes) {
    return multiply_ordered(nodes, pow<I>(Ids{})...);
}

template <std::size_t Idx, std::size_t... I>
constexpr auto get(std::index_sequence<I...>) {
    return std::array<std::size_t, sizeof...(I)>{I...}[Idx];
}

// template <std::size_t... I, class UnivariateArray>
// constexpr void mult_univariate(std::index_sequence<> /* is */,
//                                UnivariateArray const & /* ua */,
//                                double & /* out */) {}

template <std::size_t I, class PartitionIntegerSequence, class UnivariateArray>
constexpr void mult_univariate(PartitionIntegerSequence is,
                               UnivariateArray const &ua, double &out) {
    if constexpr (I < is.size()) {
        if constexpr (get<I>(is) == 1) {
            std::cout << is.size() - 1 - I << std::endl;
            out *= ua[is.size() - 1 - I];
        } else if constexpr (get<I>(is) > 1) {
            std::cout << is.size() - 1 - I << std::endl;
            out *= std::pow(ua[is.size() - 1 - I], get<I>(is));
        }

        mult_univariate<I + 1>(is, ua, out);
    }
}

template <std::size_t N, class PartitionIntegerSequence, class OpExpansion,
          class NodeValues, class UnivariateArray, class ArrayOut>
void backpropagate_process_univariate(PartitionIntegerSequence const current,
                                      OpExpansion in, NodeValues nodes,
                                      UnivariateArray const &ua,
                                      ArrayOut &out) {

    // std::cout << get<N>(current);
    out[N] = BellCoeff(current);
    mult_univariate<0>(current, ua, out[N]);
    std::cout << out[N] << std::endl;

    constexpr auto this_op = get_op(current, in, nodes);
    std::cout << type_name2<decltype(this_op)>() << std::endl;

    auto constexpr next = NextPartitionIS(current);
    if constexpr (!std::is_same_v<decltype(current), decltype(next)>) {
        return backpropagate_process_univariate<N + 1>(next, in, nodes, ua,
                                                       out);
    }
}

template <std::size_t Power, std::size_t Order,
          template <class> class Univariate, class NodeDerivative,
          class FirstNodesDerivativeRest, class CalcTreeValue,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, class UnivariateType, class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, Univariate<NodeDerivative>>> /* nd */,
    FirstNodesDerivativeRest ndr, CalcTreeValue ct, InterfaceTypes it,
    InterfaceArray ia, BufferTypes bt, BufferArray ba, UnivariateType /* ut*/,
    UnivariateArray &ua) {
    using OutId = Univariate<NodeDerivative>;
    using InId = NodeDerivative;
    constexpr bool is_node_in_interface =
        has_type_tuple_v<InId, InterfaceTypes>;

    // constexpr auto position_interface = get_idx_first2<InId>(it);
    // constexpr auto position_buffer = get_idx_first2<InId>(bt);
    // double &node_val =
    //     is_node_in_interface ? ia[position_interface] : ba[position_buffer];

    constexpr bool is_repeated_univariate =
        std::is_same_v<OutId, UnivariateType>;
    if constexpr (!is_repeated_univariate) {
        // we need to calculate derivative values
        auto univ_ders =
            OutId::template d2<Order>(ct.val(OutId{}), ct.val(InId{}));

        std::copy(univ_ders.begin(), univ_ders.end(), ua.begin());
    }

    using NodesValue = decltype(ct)::ValuesTupleInverse;

    constexpr auto in = detail::expand_univariate_initial<InId, Order>(
        std::make_index_sequence<Order>{});

    std::array<double, partition_function(Order)> univariate_expansion;

    backpropagate_process_univariate<0>(FirstPartitionIS<Order>(), in,
                                        NodesValue{}, ua, univariate_expansion);
}

template <class NodeDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue ct, InterfaceTypes it,
                       InterfaceArray ia, BufferTypes bt, BufferArray ba,
                       UnivariateType ut, UnivariateArray ua) {

    if constexpr (std::tuple_size_v<NodeDerivatives>) {
        constexpr auto current_derivative_node = head(nd);

        constexpr auto current_derivative_subnode =
            head(current_derivative_node);
        constexpr auto current_derivative_subnode_rest =
            tail(current_derivative_node);

        backpropagate_process(current_derivative_subnode,
                              current_derivative_subnode_rest, ct, it, ia, bt,
                              ba, ut, ua);

        backpropagate_aux(tail(nd), ct, it, ia, bt, ba, ut, ua);
    }
}

} // namespace adhoc3::detail

#endif // ADHOC3_TAPE_BACKPROPAGATOR_HPP
