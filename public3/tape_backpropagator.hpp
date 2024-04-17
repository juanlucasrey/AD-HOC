#ifndef ADHOC3_TAPE_BACKPROPAGATOR_HPP
#define ADHOC3_TAPE_BACKPROPAGATOR_HPP

#include <tape_buffer_size.hpp>
#include <tuple>
#include <tuple_utils.hpp>

#include <partition/multinomial_coefficient_index_sequence.hpp>

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

template <std::size_t I = 0, class IntegerSequence, class InputArray>
constexpr void mult_univariate(IntegerSequence is, InputArray const &ua,
                               double &out) {
    if constexpr (I < is.size()) {
        if constexpr (get<I>(is) == 1) {
            out *= ua[is.size() - 1 - I];
        } else if constexpr (get<I>(is) > 1) {
            out *= std::pow(ua[is.size() - 1 - I], get<I>(is));
        }

        mult_univariate<I + 1>(is, ua, out);
    }
}

template <std::size_t I = 0, class IntegerSequence, class InputArray>
constexpr void mult_multinomial(IntegerSequence is, InputArray const &in,
                                double &out) {
    if constexpr (I < is.size()) {
        if constexpr (get<I>(is) == 1) {
            out *= in[I];
        } else if constexpr (get<I>(is) > 1) {
            out *= std::pow(in[I], get<I>(is));
        }

        mult_multinomial<I + 1>(is, in, out);
    }
}

template <std::size_t N, std::size_t Power, class PartitionIntegerSequence,
          class CalculateFlags, class Array, class ArrayOut>
auto backpropagate_process_power(PartitionIntegerSequence const current,
                                 CalculateFlags calc_flags, Array const &vals,
                                 ArrayOut &arrayout) {

    if constexpr (std::get<N>(calc_flags)) {
        constexpr auto factorials = factorialarr<Power>();
        arrayout[N] = MultinomialCoeff<Power>(factorials, current);
        mult_multinomial(current, vals, arrayout[N]);
    }

    auto constexpr next = NextMultinomialIS(current);

    if constexpr (!std::is_same_v<decltype(current), decltype(next)>) {
        return backpropagate_process_power<N + 1, Power>(next, calc_flags, vals,
                                                         arrayout);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences,
          class CalculateFlags, class Array, class ArrayOut>
auto backpropagate_univariate(PartitionIntegerSequences const sequences,
                              CalculateFlags const calc_flags,
                              Array const &vals, ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        if constexpr (std::get<N>(calc_flags)) {
            constexpr auto current_sequence = std::get<N>(sequences);
            arrayout[N] = BellCoeff(current_sequence);
            mult_univariate(current_sequence, vals, arrayout[N]);
        }

        return backpropagate_univariate<N + 1>(sequences, calc_flags, vals,
                                               arrayout);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences, class Array,
          class ArrayOut>
auto backpropagate_multivariate(PartitionIntegerSequences const sequences,
                                Array const &vals, ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        arrayout[N] = MultinomialCoeff2(current_sequence);
        mult_multinomial(current_sequence, vals, arrayout[N]);

        backpropagate_multivariate<N + 1>(sequences, vals, arrayout);
    }
}

template <std::size_t Power, std::size_t Order,
          template <class> class Univariate, class NodeDerivative,
          class FirstNodesDerivativeRest, class NextNodesDerivatives,
          class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class UnivariateType,
          class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, Univariate<NodeDerivative>>> /* nd */,
    FirstNodesDerivativeRest ndr, NextNodesDerivatives next_derivative_nodes,
    CalcTreeValue ct, InterfaceTypes it, InterfaceArray ia, BufferTypes bt,
    BufferArray ba, UnivariateType /* ut*/, UnivariateArray &ua) {

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

    constexpr auto univariate_size = partition_function(Order);
    std::array<double, univariate_size> univariate_expansion_values{0};

    constexpr auto univariate_expansion_types =
        expand_univariate<InId, Order>(NodesValue{});

    constexpr auto multinomial_expansion_types =
        expand_multinomial<Power>(NodesValue{}, univariate_expansion_types);

    constexpr auto multinomial_expansion_types_full =
        multiply_ordered_tuple(NodesValue{}, multinomial_expansion_types, ndr);

    constexpr auto multinomial_calc_flags = std::apply(
        [next_derivative_nodes, it](auto... type) {
            return std::tuple_cat(
                std::conditional_t < contains(next_derivative_nodes, type) ||
                    contains(it, type),
                std::tuple<std::true_type>,
                std::tuple < std::false_type >> {}...);
        },
        multinomial_expansion_types_full);

    constexpr auto multinomial_sequences =
        MultinomialSequences<univariate_size, Power>();

    constexpr auto multinomial_sequences_filtered =
        filter(multinomial_sequences, multinomial_calc_flags);
    constexpr auto multinomial_sequences_filtered_size =
        std::tuple_size_v<decltype(multinomial_sequences_filtered)>;

    constexpr auto univariate_calc_flags =
        precedent_required(multinomial_sequences_filtered);

    constexpr auto partition_sequences = PartitionSequences<Order>();

    backpropagate_univariate(partition_sequences, univariate_calc_flags, ua,
                             univariate_expansion_values);

    constexpr auto multinomial_size = combinations(univariate_size, Power);
    std::array<double, multinomial_size> multinomial_expansion_values{0};
    backpropagate_process_power<0, Power>(
        FirstMultinomialIS<univariate_size, Power>(), multinomial_calc_flags,
        univariate_expansion_values, multinomial_expansion_values);

    std::array<double, multinomial_sequences_filtered_size>
        multinomial_expansion_values2{0};

    backpropagate_multivariate(multinomial_sequences_filtered,
                               univariate_expansion_values,
                               multinomial_expansion_values2);

    constexpr auto multinomial_expansion_types_full_filtered =
        filter(multinomial_expansion_types_full, multinomial_calc_flags);

    std::cout
        << type_name2<decltype(multinomial_expansion_types_full_filtered)>()
        << std::endl;
}

template <class NodeDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue ct, InterfaceTypes it,
                       InterfaceArray ia, BufferTypes bt, BufferArray ba,
                       UnivariateType ut, UnivariateArray ua) {

    if constexpr (std::tuple_size_v<NodeDerivatives>) {
        constexpr auto current_derivative_node = head(nd);
        constexpr auto next_derivative_nodes = tail(nd);

        constexpr auto current_derivative_subnode =
            head(current_derivative_node);
        constexpr auto current_derivative_subnode_rest =
            tail(current_derivative_node);

        backpropagate_process(
            current_derivative_subnode, current_derivative_subnode_rest,
            next_derivative_nodes, ct, it, ia, bt, ba, ut, ua);

        backpropagate_aux(next_derivative_nodes, ct, it, ia, bt, ba, ut, ua);
    }
}

} // namespace adhoc3::detail

#endif // ADHOC3_TAPE_BACKPROPAGATOR_HPP
