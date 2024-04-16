#ifndef ADHOC3_TAPE_BACKPROPAGATOR_HPP
#define ADHOC3_TAPE_BACKPROPAGATOR_HPP

#include <tape_buffer_size.hpp>
#include <tuple>
#include <tuple_utils.hpp>

#include "../tests3/type_name.hpp"

#include <array>

namespace adhoc3 {

namespace detail {
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

// template <std::size_t Power, std::size_t Order,
//           template <class> class Univariate, class NodeDerivative,
//           class FirstNodesDerivativeRest, class NodesDerivativeRest,
//           class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
//           class BufferTypes, class BufferArray, class UnivariateType,
//           class UnivariateArray>
// void backpropagate_node_sep(
//     der2::p<Power, der2::d<Order, Univariate<NodeDerivative>>> /* nd */,
//     FirstNodesDerivativeRest /* fndr */, NodesDerivativeRest /* ndr */,
//     CalcTreeValue ct, InterfaceTypes it, InterfaceArray ia, BufferTypes bt,
//     BufferArray ba, UnivariateType /* ut
//                                     */
//     ,
//     UnivariateArray ua) {
//     using OutId = Univariate<NodeDerivative>;
//     using InId = NodeDerivative;
//     constexpr bool is_node_in_interface =
//         has_type_tuple_v<InId, InterfaceTypes>();

//     // constexpr auto position_interface = get_idx_first2<InId>(it);
//     // constexpr auto position_buffer = get_idx_first2<InId>(bt);
//     // double &node_val =
//     //     is_node_in_interface ? ia[position_interface] :
//     ba[position_buffer];

//     // constexpr bool is_repeated_univariate =
//     //     std::is_same_v<OutId, UnivariateType>;
//     // if constexpr (!is_repeated_univariate) {
//     //     // we need to calculate derivative values
//     //     auto univ_ders = OutId::d2<Order>(ct.val(OutId{}),
//     ct.val(InId{}));

//     //     std::copy(univ_ders.begin(), univ_ders.end(), ua.begin);
//     // }

//     // using NodesValue = decltype(ct)::ValuesTupleInverse;
//     // constexpr auto nodes_value = NodesValue{};
//     // // constexpr auto expansion_result =
//     // //     expand_single(nodes_value, Univariate<NodeDerivative>{});

//     // // constexpr auto result_filtered = std::apply(
//     // //     [nd, ia](auto... single_output) {
//     // //         return std::tuple_cat(conditional_node2(nd, ia,
//     // //         single_output)...);
//     // //     },
//     // //     expansion_result);
//     // {
//     //     constexpr auto in = detail::expand_univariate_initial<InId,
//     Order>(
//     //         std::make_index_sequence<Order>{});
//     //     std::cout << type_name2<decltype(in)>() << std::endl;
//     // }
// }

// template <class FirstNodesDerivative, class... FirstNodesDerivativeRest,
//           class NodesDerivativeRest, class CalcTreeValue, class
//           InterfaceTypes, class InterfaceArray, class BufferTypes, class
//           BufferArray, class UnivariateType, class UnivariateArray>
// void backpropagate_node(
//     std::tuple<FirstNodesDerivative, FirstNodesDerivativeRest...> /* nd */,
//     NodesDerivativeRest ndr, CalcTreeValue ct, InterfaceTypes it,
//     InterfaceArray ia, BufferTypes bt, BufferArray ba, UnivariateType ut,
//     UnivariateArray ua) {
//     backpropagate_node_sep(FirstNodesDerivative{},
//                            std::tuple<FirstNodesDerivativeRest...>{}, ndr,
//                            ct, it, ia, bt, ba, ut, ua);
// }

// template <class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
//           class BufferTypes, class BufferArray, class UnivariateType,
//           class UnivariateArray>
// void backpropagate_aux(std::tuple<> /* nd */, CalcTreeValue /* ct */,
//                        InterfaceTypes /* it */, InterfaceArray /* ia */,
//                        BufferTypes /* bt */, BufferArray /* ba */,
//                        UnivariateType /* ut */, UnivariateArray /* ua */) {}

template <std::size_t Power, std::size_t Order,
          template <class> class Univariate, class NodeDerivative,
          class FirstNodesDerivativeRest, class CalcTreeValue,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, class UnivariateType, class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, Univariate<NodeDerivative>>> /* nd */,
    FirstNodesDerivativeRest /* ndr */, CalcTreeValue ct, InterfaceTypes it,
    InterfaceArray ia, BufferTypes bt, BufferArray ba, UnivariateType /* ut
                                                                       */
    ,
    UnivariateArray ua) {
    using OutId = Univariate<NodeDerivative>;
    using InId = NodeDerivative;
    constexpr bool is_node_in_interface =
        has_type_tuple_v<InId, InterfaceTypes>;

    // constexpr auto position_interface = get_idx_first2<InId>(it);
    // constexpr auto position_buffer = get_idx_first2<InId>(bt);
    // double &node_val =
    //     is_node_in_interface ? ia[position_interface] : ba[position_buffer];

    // constexpr bool is_repeated_univariate =
    //     std::is_same_v<OutId, UnivariateType>;
    // if constexpr (!is_repeated_univariate) {
    //     // we need to calculate derivative values
    //     auto univ_ders = OutId::d2<Order>(ct.val(OutId{}), ct.val(InId{}));

    //     std::copy(univ_ders.begin(), univ_ders.end(), ua.begin);
    // }

    // using NodesValue = decltype(ct)::ValuesTupleInverse;
    // constexpr auto nodes_value = NodesValue{};
    // // constexpr auto expansion_result =
    // //     expand_single(nodes_value, Univariate<NodeDerivative>{});

    // // constexpr auto result_filtered = std::apply(
    // //     [nd, ia](auto... single_output) {
    // //         return std::tuple_cat(conditional_node2(nd, ia,
    // //         single_output)...);
    // //     },
    // //     expansion_result);
    {
        constexpr auto in = detail::expand_univariate_initial<InId, Order>(
            std::make_index_sequence<Order>{});
        std::cout << type_name2<decltype(in)>() << std::endl;
    }
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

} // namespace detail

} // namespace adhoc3

#endif // ADHOC3_TAPE_BACKPROPAGATOR_HPP
