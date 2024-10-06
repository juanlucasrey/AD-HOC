#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>
#include <format>
#include <utils/bivariate.hpp>

#include <utils/tuple.hpp>

#include <sort.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

#include <random>
#include <tuple>

namespace adhoc4 {

namespace {

// namespace detail2 {

// template <class... Ids, std::size_t... Orders, class Node>
// constexpr auto
// create_single_tuple(std::tuple<der::d<Orders, Ids>...> diff_operator, Node
// in) {

//     std::size_t constexpr pos = get_first_type_idx(std::tuple<Ids...>{}, in);
//     if constexpr (pos < sizeof...(Ids)) {
//         constexpr auto v = std::get<pos>(diff_operator);
//         return std::tuple<decltype(v)>{};
//     } else {
//         return std::tuple<>{};
//     }
// };

// } // namespace detail2

// template <class... Ids, std::size_t... Orders, class... CalculationNodes>
// constexpr auto
// order_differential_operator(std::tuple<der::d<Orders, Ids>...> diff_operator,
//                             std::tuple<CalculationNodes...> nodes) {
//     return std::apply(
//         [diff_operator](auto &&...args) {
//             return std::tuple_cat(
//                 detail2::create_single_tuple(diff_operator, args)...);
//         },
//         nodes);
// }

} // namespace

TEST(BackPropagator, IsInputDer) {
    ADHOC(var);
    ADHOC(var2);
    auto calc = exp(var);
    auto res = d(var) * d<2>(var2);

    constexpr auto resb1 = detail::is_derivative_input(res);
    static_assert(resb1);

    auto res2 = d(calc) * d<2>(var2);

    constexpr auto resb2 = detail::is_derivative_input(res2);
    static_assert(!resb2);
}

namespace {

// template <class Result>
// constexpr auto calc_tree_inverse_no_input_aux(Result result,
//                                               std::tuple<> /* nodes */) {
//     return result;
// }

// template <class... ResultNodes, class Node, class... Nodes>
// constexpr auto
// calc_tree_inverse_no_input_aux(std::tuple<ResultNodes...> res,
//                                std::tuple<Node, Nodes...> /* nodes */) {

//     constexpr bool other_nodes_depend_on_this_node =
//         (depends(Nodes{}, Node{}) || ...);

//     constexpr bool is_const = is_constant(Node{});

//     if constexpr (other_nodes_depend_on_this_node || is_const) {
//         // this node will come up again because it is included on the Nodes
//         return calc_tree_inverse_no_input_aux(res, std::tuple<Nodes...>{});
//     } else {
//         constexpr auto deps = detail::get_dependents(Node{});
//         constexpr auto new_nodes = std::tuple_cat(deps,
//         std::tuple<Nodes...>{}); if constexpr (is_input(Node{})) {
//             // we don't add it
//             return calc_tree_inverse_no_input_aux(res, new_nodes);
//         } else {
//             // we add it to the result
//             return calc_tree_inverse_no_input_aux(
//                 std::tuple<ResultNodes..., Node>{}, new_nodes);
//         }
//     }
// }

// template <class... Outputs>
// constexpr auto calc_tree_inverse_no_input(CalcTree<Outputs...> /* o */) {
//     return calc_tree_inverse_no_input_aux(std::tuple<>{},
//                                           std::tuple<Outputs...>{});
// }

} // namespace

TEST(BackPropagator, UnivariateLogExp) {
    ADHOC(x);
    auto res = log(exp(x));

    CalcTree ct(res);

    const std::size_t evals = 100;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(-3., 3.);

    for (std::size_t i = 0; i < evals; i++) {
        ct.set(x) = dis(gen);
        ct.evaluate();

        BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x), d<4>(x), d<5>(x));
        bp.set(d(res)) = 1.;
        bp.backpropagate(ct);

        EXPECT_NEAR(bp.get(d(x)), 1., 1e-14);
        EXPECT_NEAR(bp.get(d<2>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<3>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<4>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<5>(x)), 0., 1e-14);
    }
}

TEST(BackPropagator, UnivariateExpLog) {
    ADHOC(x);
    auto res = exp(log(x));

    CalcTree ct(res);

    const std::size_t evals = 100;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(0.1, 3.);

    for (std::size_t i = 0; i < evals; i++) {
        ct.set(x) = dis(gen);
        ct.evaluate();

        BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x), d<4>(x), d<5>(x));
        bp.set(d(res)) = 1.;
        bp.backpropagate(ct);

        EXPECT_NEAR(bp.get(d(x)), 1., 1e-14);
        EXPECT_NEAR(bp.get(d<2>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<3>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<4>(x)), 0., 1e-13);
        EXPECT_NEAR(bp.get(d<5>(x)), 0., 1e-12);
    }
}

TEST(BackPropagator, UnivariateSingle) {
    ADHOC(x);
    auto res = erfc(x);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = erfc(x)

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -0.267344347003539, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), 0.320813216404247, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -0.167535790788885, 1e-14);
}

TEST(BackPropagator, UnivariateDouble) {
    ADHOC(x);
    auto res = log(erfc(x));

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = log(erfc(x))

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -2.98089202449025, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), -0.865788201446497, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -0.0342861937391300, 1e-14);
}

// // default for sum and subs
// template <std::size_t Order, template <class, class> class Bivariate, class
// Id1,
//           class Id2>
// constexpr auto expand_tree_bivariate(Bivariate<Id1, Id2> /* id */) {
//     return std::tuple<std::tuple<der::d<Order, Id1>>,
//                       std::tuple<der::d<Order, Id2>>>{};
// }

// template <std::size_t Order, class IdFirst, class IdSecond, std::size_t... I>
// constexpr auto generate_operators_dflt(std::index_sequence<I...> /* i */) {
//     return std::make_tuple(std::tuple<der::d<Order, IdFirst>>{},
//                            std::tuple<der::d<Order - I - 1, IdFirst>,
//                                       der::d<I + 1, IdSecond>>{}...,
//                            std::tuple<der::d<Order, IdSecond>>{});
// }

// // default for sum and subs
// template <std::size_t Order, template <class, class> class Bivariate, class
// Id1,
//           class Id2>
// constexpr auto expand_tree_bivariate(Bivariate<Id1, Id2> /* id */) {
//     constexpr auto seq = std::make_index_sequence<Order - 1>{};
//     return generate_operators_dflt<Order, Id1, Id2>(seq);
// }

// template <template <class, class> class Bivariate, class Id1, class Id2,
//           std::size_t Order, class Nodes>
// constexpr auto expand_tree_single(der::d<Order, Bivariate<Id1, Id2>>
//                                   /* id */,
//                                   Nodes nodes) {
//     if constexpr (is_constant(Id1{})) {
//         return std::tuple<std::tuple<der::d<Order, Id2>>>{};
//     } else if constexpr (is_constant(Id2{})) {
//         return std::tuple<std::tuple<der::d<Order, Id1>>>{};
//     } /* else if constexpr (std::is_same_v<Id1, Id2>) {
//         // todo
//     } */
//     else {
//         constexpr auto ordered_bivariate =
//             detail::sort_bivariate(Bivariate<Id1, Id2>{}, nodes);
//         return expand_tree_bivariate<Order>(ordered_bivariate);
//     }
// }

TEST(BackPropagator, ExpandSum) {
    ADHOC(x);
    ADHOC(y);
    auto res = x + y;

    CalcTree ct(res);
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;
    constexpr auto primal_nodes_inverted = PrimalNodes{};

    auto der = d<4>(res);
    auto der_single = head(der);

    constexpr auto expanded_d =
        expand_bivariate(der_single, primal_nodes_inverted);
    std::cout << type_name2<decltype(expanded_d)>() << std::endl;
}

TEST(BackPropagator, Sort) {
    ADHOC(x);
    ADHOC(y);
    ADHOC(z);
    auto res = x + (y * z);
    auto temp = (z * y);
    CalcTree ct(res);
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;
    constexpr auto primal_nodes_inverted = PrimalNodes{};

    constexpr auto new_bi = detail::sort_bivariate(temp, primal_nodes_inverted);
    std::cout << type_name2<decltype(new_bi)>() << std::endl;

    auto d1 = d(res);
    auto d2 = d(res) * d(x);
    auto d3 = d(res) * d(y);
    auto d4 = d(res) * d(z);
    auto d5 = d<2>(res);
    auto d6 = d<2>(res) * d<2>(x);

    constexpr auto newtuple = sort_differential_operators(
        std::make_tuple(d4, d3, d2, d1, d5, d6), primal_nodes_inverted);
    std::cout << type_name2<decltype(newtuple)>() << std::endl;
}

namespace detail {

// template <class In1, class In2, class Nodes>
// constexpr auto less_than_check_first2(In1 in1, In2 in2, Nodes nodes) {
//     return less_than_single(head(in1), head(in2), tail(in1), tail(in2),
//     nodes);
// }

// template <class Op1, class Op2, class Nodes>
// constexpr auto less_than_check_empty2(Op1 in1, Op2 in2, Nodes nodes) {
//     if constexpr (std::is_same_v<const Op2, const std::tuple<>>) {
//         return false;
//     } else if constexpr (std::is_same_v<const Op1, const std::tuple<>>) {
//         return true;
//     } else {
//         // we can now assume that both ids are different and non empty
//         return less_than_check_first2(in1, in2, nodes);
//     }
// }

// template <class Op1, class Op2, class Nodes>
// constexpr auto less_than2(Op1 in1, Op2 in2, Nodes nodes) {
//     if constexpr (std::is_same_v<Op1, Op2>) {
//         return false;
//     } else {
//         return less_than_check_empty2(in1, in2, nodes);
//     }
// }

// template <class Tuple, class Nodes>
// constexpr auto is_sorted2(Tuple in, Nodes nodes) {
//     if constexpr (size(in) <= 1) {
//         return true;
//     } else {
//         if constexpr (!less_than_single(std::get<1>(in), std::get<0>(in),
//                                         nodes)) {
//             return is_sorted2(tail(in), nodes);
//         } else {
//             return false;
//         }
//     }
// }
} // namespace detail

// TEST(BackPropagator, Sorted) {
//     ADHOC(y);
//     ADHOC(x);
//     auto res = log(x + y);
//     CalcTree ct(res);
//     using NodesValue = decltype(ct)::ValuesTupleInverse;

//     constexpr auto d1 = d(x) * d(y);
//     std::cout << "unsorted" << std::endl;
//     std::cout << type_name2<decltype(d1)>() << std::endl;
//     std::cout << "NodesValue" << std::endl;
//     std::cout << type_name2<NodesValue>() << std::endl;
//     auto resbool = detail::is_sorted2(d1, NodesValue{});
//     static_assert(detail::is_sorted2(d1, NodesValue{}));
// }

TEST(BackPropagator, Size) {
    ADHOC(x);
    ADHOC(y);
    auto c = constants::CD<0.12345>();
    auto z = cos(exp(x) * c);
    // std::cout << decltype(z) << std::endl;

    // static_assert(sizeof(z) == 1);
    // using constants::CD;
    // auto c = CD<0.12345>();
    // static_assert(sizeof(c) == 1);

    // ADHOC(S);
    // ADHOC(K);
    // ADHOC(v);
    // ADHOC(T);

    // auto res = call_price(S, K, v, T);
    // std::cout << type_name2<decltype(res)>() << std::endl;
}

TEST(BackPropagator, SumSingle) {
    ADHOC(x);
    ADHOC(y);
    // auto res = x + y;
    auto res = log(x + y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    std::cout << bp.get(d(x)) << std::endl;
    std::cout << bp.get(d(y)) << std::endl;
    std::cout << bp.get(d(x) * d(y)) << std::endl;

    // EXPECT_NEAR(bp.get(d(x)), -2.98089202449025, 1e-14);
    // EXPECT_NEAR(bp.get(d<2>(x)), -0.865788201446497, 1e-14);
    // EXPECT_NEAR(bp.get(d<3>(x)), -0.0342861937391300, 1e-14);
}

// namespace {

// template <class Numerator, class... Denominators>
// /* constexpr */ auto
// convex_hull_orders_2(Numerator num, std::tuple<Denominators...> /* den */) {
//     std::cout << "num" << std::endl;
//     std::cout << type_name2<decltype(num)>() << std::endl;

//     constexpr auto expanded = expand(num);
//     std::cout << "expanded" << std::endl;
//     std::cout << type_name2<decltype(expanded)>() << std::endl;

//     constexpr auto orders = order_vars2_aux(expanded, Denominators{}...);
//     std::cout << "orders" << std::endl;
//     std::cout << type_name2<decltype(orders)>() << std::endl;

//     return filter_smaller(orders);
// }

// template <class Res, class Inputs>
// /* constexpr */ auto expand_aux_2(Res res, std::tuple<> /* seqs */,
//                                   Inputs /* inputs */) {
//     return res;
// }

// template <class Res, class DiffOpOut, class Inputs>
// /* constexpr */ auto expand_aux_2(Res res, DiffOpOut diff_op_out,
//                                   Inputs inputs) {
//     constexpr auto first = head(diff_op_out);
//     std::cout << "first" << std::endl;
//     std::cout << type_name2<decltype(first)>() << std::endl;

//     constexpr auto last = tail(diff_op_out);
//     std::cout << "last" << std::endl;
//     std::cout << type_name2<decltype(last)>() << std::endl;

//     /* constexpr */ auto ch = convex_hull_orders_2(get_id(first), inputs);
//     std::cout << "ch" << std::endl;
//     std::cout << type_name2<decltype(ch)>() << std::endl;

//     constexpr auto new_res = std::tuple_cat(
//         res, expand_aux2(ch, std::make_index_sequence<get_power(first)>()));

//     std::cout << "new_res" << std::endl;
//     std::cout << type_name2<decltype(new_res)>() << std::endl;

//     return expand_aux_2(new_res, last, inputs);
// }

// template <class DiffOptOut, class Inputs>
// /* constexpr */ auto expand_2(DiffOptOut diff_op_out, Inputs inputs) {
//     return expand_aux_2(std::tuple<>{}, diff_op_out, inputs);
// }

// template <class... IndexSequence, class SingleIndexSequence>
// /* constexpr */ auto flatten_aux4_2(std::tuple<IndexSequence...> /* in1 */,
//                                     SingleIndexSequence in2) {
//     return std::make_tuple(std::tuple_cat(std::make_tuple(IndexSequence{}),
//                                           std::make_tuple(in2))...);
// }

// template <class TupleRes, class... IndexSequence>
// /* constexpr */ auto flatten_aux3_2(TupleRes in1,
//                                     std::tuple<IndexSequence...> /* in2 */) {
//     return std::tuple_cat(flatten_aux4_2(in1, IndexSequence{})...);
// }

// template <class... TupleRes, class... IndexSequence>
// /* constexpr */ auto flatten_aux2_2(std::tuple<TupleRes...> /* res */,
//                                     std::tuple<IndexSequence...> newinputs) {
//     return std::tuple_cat(flatten_aux3_2(TupleRes{}, newinputs)...);
// }

// template <class Res>
// /* constexpr */ auto flatten_aux_2(Res res, std::tuple<> /* idx_seqs */) {
//     return res;
// }

// template <class Res, class... IndexSequences, class... TupleIndexSequences>
// /* constexpr */ auto
// flatten_aux_2(Res res,
//               std::tuple<std::tuple<IndexSequences...>,
//               TupleIndexSequences...>
//                   idx_seqs) {
//     std::cout << "idx_seqs" << std::endl;
//     std::cout << type_name2<decltype(idx_seqs)>() << std::endl;
//     if constexpr (std::tuple_size_v<Res> == 0) {
//         constexpr auto new_res =
//             std::make_tuple(std::tuple<IndexSequences...>{});
//         std::cout << "new_res" << std::endl;
//         std::cout << type_name2<decltype(new_res)>() << std::endl;
//         return flatten_aux_2(new_res, std::tuple<TupleIndexSequences...>{});
//     } else {
//         /* constexpr */ auto new_res =
//             flatten_aux2_2(res, std::tuple<IndexSequences...>{});
//         std::cout << "new_res" << std::endl;
//         std::cout << type_name2<decltype(new_res)>() << std::endl;
//         return flatten_aux_2(new_res, std::tuple<TupleIndexSequences...>{});
//     }
// }

// template <class... IndexSequences, class... TupleIndexSequences>
// /* constexpr */ auto
// flatten_2(std::tuple<std::tuple<IndexSequences...>, TupleIndexSequences...>
//               idx_seqs) {
//     return flatten_aux_2(std::tuple<>{}, idx_seqs);
// }

// template <class Single, class Many>
// /* constexpr */ auto monomial_included2(Single diff_op_out, Many
// diff_op_inputs)
//     -> bool {
//     std::cout << "diff_op_out" << std::endl;
//     std::cout << type_name2<decltype(diff_op_out)>() << std::endl;

//     std::cout << "diff_op_inputs" << std::endl;
//     std::cout << type_name2<decltype(diff_op_inputs)>() << std::endl;

//     constexpr auto inputs = get_vars(diff_op_inputs);

//     std::cout << "inputs" << std::endl;
//     std::cout << type_name2<decltype(inputs)>() << std::endl;

//     constexpr auto orders = expand(diff_op_out, inputs);
//     // auto orders = expand_2(diff_op_out, inputs);

//     std::cout << "orders" << std::endl;
//     std::cout << type_name2<decltype(orders)>() << std::endl;

//     /* constexpr */ auto flat_orders_output = flatten_2(orders);

//     std::cout << "flat_orders_output" << std::endl;
//     std::cout << type_name2<decltype(flat_orders_output)>() << std::endl;

//     constexpr auto order_inputs_vals = order_inputs(diff_op_inputs, inputs);

//     std::cout << "order_inputs_vals" << std::endl;
//     std::cout << type_name2<decltype(order_inputs_vals)>() << std::endl;

//     return check_included(flat_orders_output, order_inputs_vals);
// }

// } // namespace

TEST(BackPropagator, MonomialIncluded) {
    ADHOC(x);
    ADHOC(y);
    constexpr auto dnin = std::make_tuple(d(x), d(y), d(x) * d(y));
    constexpr auto node = d(x) * d<2>(y);
    constexpr auto resbool = monomial_included(node, dnin);
    static_assert(!resbool);
}

// TEST(BackPropagator, MonomialIncluded) {
//     ADHOC(x);
//     ADHOC(y);
//     // auto mul = x * y;
//     // auto res = log(mul);

//     // CalcTree ct(res);
//     // decltype(ct)::ValuesTuple tree;
//     // using NodesValue = decltype(ct)::ValuesTupleInverse;

//     constexpr auto dnin = std::make_tuple(d(x), d(y), d(x) * d(y));

//     // std::cout << "dnin" << std::endl;
//     // std::cout << type_name2<decltype(dnin)>() << std::endl;

//     // constexpr auto current_node_der = d<2>(mul);
//     // constexpr auto current_derivative_subnode = head(current_node_der);

//     // std::cout << "current_derivative_subnode" << std::endl;
//     // std::cout << type_name2<decltype(current_derivative_subnode)>()
//     //           << std::endl;

//     // constexpr auto expansion_types =
//     //     expand_bivariate(current_derivative_subnode, NodesValue{});

//     // std::cout << "expansion_types" << std::endl;
//     // std::cout << type_name2<decltype(expansion_types)>() << std::endl;

//     // constexpr auto flags_next_derivatives = std::apply(
//     //     [dnin](auto... next_derivative) {
//     //         return std::tuple_cat(
//     //             std::conditional_t<monomial_included(next_derivative,
//     dnin),
//     //                                std::tuple<std::true_type>,
//     //                                std::tuple<std::false_type>>{}...);
//     //     },
//     //     expansion_types);

//     // std::cout << "flags_next_derivatives" << std::endl;
//     // std::cout << type_name2<decltype(flags_next_derivatives)>() <<
//     std::endl;

//     constexpr auto node = d(x) * d<2>(y);
//     constexpr auto resbool = monomial_included(node, dnin);
//     std::cout << resbool << std::endl;
// }

TEST(BackPropagator, MultSingle) {
    ADHOC(x);
    ADHOC(y);
    // auto res = x + y;
    auto res = log(x * y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    std::cout << bp.get(d(x)) << std::endl;
    std::cout << bp.get(d(y)) << std::endl;
    std::cout << bp.get(d(x) * d(y)) << std::endl;

    // EXPECT_NEAR(bp.get(d(x)), -2.98089202449025, 1e-14);
    // EXPECT_NEAR(bp.get(d<2>(x)), -0.865788201446497, 1e-14);
    // EXPECT_NEAR(bp.get(d<3>(x)), -0.0342861937391300, 1e-14);
}

// template <class Tuple, class Old, class New>
// auto constexpr replace_first(Tuple tuple, Old old_value, New new_value) {
//     return detail::replace_first_aux(tuple, old_value, new_value,
//                                      std::tuple<>{});
// }

// template <class DiffOp, class Node>
// constexpr auto first_type_is2(DiffOp diffop, Node /* node */) -> bool {
//     constexpr auto first = head(diffop);
//     // return first_type_is_aux(head(diffop), node);
// }

TEST(BackPropagator, ReplaceFirst) {

    constexpr std::tuple<double, int, char, int> temp;
    constexpr int v1 = 0;
    constexpr float v2 = 0;

    // constexpr auto first = head(temp);

    // constexpr auto newtuple = first_type_is2(temp, v1);
    constexpr auto newtuple = replace_first(temp, v1, v2);
    std::cout << type_name2<decltype(newtuple)>() << std::endl;

    constexpr auto newtuple2 = replace_first(temp, v2, v1);
    std::cout << type_name2<decltype(newtuple2)>() << std::endl;
}

// template <std::size_t I, class T> using Type = T;

// template <class T, std::size_t... Is>
// constexpr auto make_tuple_aux(std::index_sequence<Is...>) {
//     return std::tuple<Type<Is, T>...>{};
// }

// template <std::size_t N, class T> constexpr auto make_tuple(T /* t */) {
//     return make_tuple_aux<T>(std::make_index_sequence<N>{});
// }

TEST(BackPropagator, SumExpansion) {
    ADHOC(v1);
    ADHOC(v2);

    auto r1 = v1 * v2;
    auto r2 = cos(v1);
    auto r = r1 + r2;
    CalcTree ct(r);
    decltype(ct)::ValuesTuple tree;
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;

    {
        constexpr auto current_node_der = d(r1 + r2);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<2, 1>();
        std::cout << "1" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result = std::make_tuple(d(r1), d(r2));
        static_assert(expansion_types == result);

        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }

    {
        constexpr auto current_node_der = d<2>(r1 + r2);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<2, 2>();
        std::cout << "2" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result =
            std::make_tuple(d(r1) * d(r2), d<2>(r1), d<2>(r2));
        static_assert(expansion_types == result);

        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }

    {
        constexpr auto current_node_der = d<3>(r1 + r2);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<2, 3>();
        std::cout << "3" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result = std::make_tuple(
            d(r1) * d<2>(r2), d<2>(r1) * d(r2), d<3>(r1), d<3>(r2));
        static_assert(expansion_types == result);

        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }
}

TEST(BackPropagator, MultExpansion) {
    ADHOC(v1);
    ADHOC(v2);

    auto r1 = v1 * v2;
    auto r2 = cos(v1);
    auto r = r1 * r2;
    CalcTree ct(r);
    decltype(ct)::ValuesTuple tree;
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;

    {
        constexpr auto current_node_der = d(r1 * r2);
        constexpr auto current_derivative_subnode = head(current_node_der);

        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<3, 1>();
        std::cout << "1" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result = std::make_tuple(d(r1) * d(r2), d(r1), d(r2));
        static_assert(expansion_types == result);

        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }

    {
        constexpr auto current_node_der = d<2>(r1 * r2);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<3, 2>();
        std::cout << "2" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result =
            std::make_tuple(d(r1) * d(r2), d(r1) * d<2>(r2), d<2>(r1) * d(r2),
                            d<2>(r1) * d<2>(r2), d<2>(r1), d<2>(r2));
        static_assert(expansion_types == result);

        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }

    {
        constexpr auto current_node_der = d<3>(r1 * r2);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, PrimalNodes{});

        constexpr auto multinomial_sequences = MultinomialSequences<3, 3>();
        std::cout << "3" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        std::cout << type_name2<decltype(multinomial_sequences)>() << std::endl;

        constexpr auto result = std::make_tuple(
            d(r1) * d<2>(r2), d(r1) * d<3>(r2), d<2>(r1) * d(r2),
            d<2>(r1) * d<2>(r2), d<2>(r1) * d<3>(r2), d<3>(r1) * d(r2),
            d<3>(r1) * d<2>(r2), d<3>(r1) * d<3>(r2), d<3>(r1), d<3>(r2));
        static_assert(expansion_types == result);
        static_assert(detail::is_sorted(expansion_types, PrimalNodes{}));
    }
}

/*TEST(BackPropagator, CallPrice) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);
    ct.set(S) = 100.0;
    ct.set(K) = 102.0;
    ct.set(v) = 0.15;
    ct.set(T) = 0.5;
    ct.evaluate();

    // constexpr decltype(ct)::ValuesTuple tree;
    // std::cout << type_name2<decltype(tree)>() << std::endl;

    BackPropagator bp(d(S), d(K), d<2>(res), d(res), d<2>(S));
    bp.get(d(S));
    // bp.get(d(T));
    bp.backpropagate(ct);

    auto tupleders =
        std::make_tuple(d(S), d(K), d<4>(res), d<2>(S), d<2>(res) * d(S));
    constexpr auto maxorder = detail::max_orders(tupleders);
    std::cout << maxorder << std::endl;

    constexpr decltype(ct)::ValuesTuple tree;
    auto diffop = d<3>(res) * d<2>(S);
    auto newdiffop = detail::order_differential_operator(diffop, tree);
    std::cout << type_name2<decltype(diffop)>() << std::endl;
    std::cout << type_name2<decltype(newdiffop)>() << std::endl;

    // d<2>(S) *
}*/

TEST(BackPropagator, Remove) {
    std::tuple<double, int, char> temp;
    std::cout << type_name2<decltype(temp)>() << std::endl;
    // double a = 0;
    // auto result = remove(temp, a);
    // std::cout << type_name2<decltype(result)>() << std::endl;

    auto inverted_temp = invert(temp);
    std::cout << type_name2<decltype(inverted_temp)>() << std::endl;
}

namespace detail {

// template <class Something> struct rettrue {
//     constexpr auto return_true(Something /* something */) { return true; }
// };

// template <class Something>
// constexpr auto return_true(Something /* something */) {
//     return true;
// }

// template <class Input, template <class> class Func, class FuncArg>
// constexpr auto filter_something(Input input, Func<FuncArg> func) {
//     constexpr auto ordered_inputs = std::apply(
//         [func](auto &&...args) {
//             return std::tuple_cat(std::make_tuple(func(args))...);
//         },
//         input);
// }

} // namespace detail

TEST(BackPropagator, ConstexprFunc) {
    std::tuple<double, int, char> temp;
    std::cout << type_name2<decltype(temp)>() << std::endl;
    // double a = 0;
    // auto result = remove(temp, a);
    // std::cout << type_name2<decltype(result)>() << std::endl;

    auto inverted_temp = invert(temp);
    std::cout << type_name2<decltype(inverted_temp)>() << std::endl;

    // detail::filter_something(temp, detail::return_true);
}

// TEST(BackPropagator, FindLargest) {
//     ADHOC(v1);
//     ADHOC(v2);

//     auto r = v1 * v2;
//     CalcTree ct(r);
//     decltype(ct)::ValuesTuple tree;

//     auto d1 = d<2>(v1);
//     auto d2 = d<2>(r);
//     auto d3 = d(r) * d(v1);
//     auto d3o = detail::order_differential_operator(d3, tree);

//     // auto ders = std::make_tuple(d1, d3o, d2);
//     auto ders = std::make_tuple(d2, d1);

//     auto derso = sort_differential_operators(ders, tree);
//     // std::cout << type_name2<decltype(ders)>() << std::endl;
//     // auto derso = detail::find_largest(ders, tree);
//     std::cout << type_name2<decltype(derso)>() << std::endl;
// }

namespace {

constexpr auto half_square_size(std::size_t max) -> std::size_t {
    return (max * (max + 1)) / 2;
}

constexpr auto half_square_index(std::size_t y, std::size_t x, std::size_t max)
    -> std::size_t {
    return x + max * y - (y * (y - 1)) / 2;
}

template <std::size_t, class Type> using return_type = Type;

template <std::size_t Size, std::size_t... I>
constexpr auto init_flags_aux(std::index_sequence<I...> /* i */) {
    return std::make_tuple(return_type<I, std::false_type>{}...);
}

template <std::size_t Size> constexpr auto init_flags() {
    return init_flags_aux<Size>(std::make_index_sequence<Size>{});
}

template <std::size_t Loc, class TupleBool, std::size_t... First,
          std::size_t... Last>
constexpr auto switch_flag_aux(TupleBool tb,
                               std::index_sequence<First...> /* i */,
                               std::index_sequence<Last...> /* i */) {
    return std::tuple_cat(
        std::make_tuple(std::get<First>(tb))...,
        std::make_tuple(
            std::conditional_t<std::tuple_element<Loc, TupleBool>::type::value,
                               std::false_type, std::true_type>{}),
        std::make_tuple(std::get<Last + Loc + 1>(tb))...);
}

template <std::size_t Loc, class TupleBool>
constexpr auto switch_flag(TupleBool tb) {
    constexpr auto total_size = std::tuple_size_v<TupleBool>;
    static_assert(total_size > Loc);
    return switch_flag_aux<Loc>(
        tb, std::make_index_sequence<Loc>{},
        std::make_index_sequence<total_size - Loc - 1>{});
}

template <std::size_t Row, std::size_t Col, std::size_t Side, std::size_t Size,
          class TupleBool>
auto set(std::array<double, Size> & /* vals */, TupleBool tb) {
    if constexpr (std::tuple_element<half_square_index(Row, Col, Side),
                                     TupleBool>::type::value) {
        return tb;
    } else {
        return tb;
    }
}

template <std::size_t Loc, std::size_t End, std::size_t Size>
auto fill_side1(std::array<double, Size> &vals, double val) {
    if constexpr (Loc < End) {
        vals[half_square_index(0, Loc, End)] =
            vals[half_square_index(0, Loc - 1, End)] * val;
        fill_side1<Loc + 1, End>(vals, val);
    }
}

template <std::size_t Loc, std::size_t End, std::size_t Size>
auto fill_side2(std::array<double, Size> &vals, double val) {
    if constexpr (Loc < End) {
        vals[half_square_index(Loc, 0, End)] =
            vals[half_square_index(Loc - 1, 0, End)] * val;
        fill_side2<Loc + 1, End>(vals, val);
    }
}

template <std::size_t Loc, std::size_t Size>
auto fill_inside(std::array<double, Size> & /* vals */) {}

template <std::size_t Order> auto mult_coeffs(double a, double b) {
    static_assert(Order > 0);

    std::array<double, half_square_size(Order + 1)> vals;
    vals[0] = 1.;

    vals[half_square_index(0, 1, Order + 1)] = a;
    fill_side1<2, Order + 1>(vals, a);

    vals[half_square_index(1, 0, Order + 1)] = b;
    fill_side2<2, Order + 1>(vals, b);

    fill_inside<Order + 1>(vals);

    return vals;
}

// template <std::size_t Loc, std::size_t Size>
// auto powers_fill(std::array<double, Size> &vals, double val) {
//     if constexpr (Loc < Size) {
//         vals[Loc] = vals[Loc - 1] * val;
//         powers_fill<Loc + 1>(vals, val);
//     }
// }

// template <std::size_t Order> auto powers(double a) {
//     static_assert(Order > 0);
//     std::array<double, Order> vals;
//     vals[0] = a;
//     powers_fill<1>(vals, a);
//     return vals;
// }

} // namespace

TEST(BackPropagator, TriangleIndex) {
    static_assert(half_square_size(5) == 15);
    static_assert(half_square_index(0, 0, 5) == 0);
    static_assert(half_square_index(0, 1, 5) == 1);
    static_assert(half_square_index(0, 2, 5) == 2);
    static_assert(half_square_index(0, 3, 5) == 3);
    static_assert(half_square_index(0, 4, 5) == 4);
    static_assert(half_square_index(1, 0, 5) == 5);
    static_assert(half_square_index(1, 1, 5) == 6);
    static_assert(half_square_index(1, 2, 5) == 7);
    static_assert(half_square_index(1, 3, 5) == 8);
    static_assert(half_square_index(2, 0, 5) == 9);
    static_assert(half_square_index(2, 1, 5) == 10);
    static_assert(half_square_index(2, 2, 5) == 11);
    static_assert(half_square_index(3, 0, 5) == 12);
    static_assert(half_square_index(3, 1, 5) == 13);
    static_assert(half_square_index(4, 0, 5) == 14);

    double a = 0.1;
    auto const res = detail::powers<5>(a);
    for (std::size_t i = 0; i < res.size(); i++) {
        std::cout << res[i] << ", ";
    }

    // double a = 0.1;
    // double b = 3.14;
    // auto res = mult_coeffs<4>(a, b);

    // for (std::size_t i = 0; i < res.size(); i++) {
    //     std::cout << res[i] << ", ";
    // }
    // std::cout << std::endl;

    // constexpr auto size = half_square_size(5);
    // constexpr auto flags = init_flags<size>();
    // constexpr auto flags1 = switch_flag<half_square_index(0, 0, 5)>(flags);
    // constexpr auto flags2 = switch_flag<half_square_index(1, 0, 5)>(flags1);
    // constexpr auto flags3 = switch_flag<half_square_index(1, 0, 5)>(flags2);

    // std::array<double, size> vals;
    // vals[half_square_index(0, 0, 5)] = 1.;
    // vals[half_square_index(1, 0, 5)] = a;
    // vals[half_square_index(0, 1, 5)] = b;

    // auto flags4 = set<0, 0, 5>(vals, flags3);

    // auto new_flags = switch_flag<half_square_index(0, 0, size)>(flags);
    // {
    //     auto flags = init_flags<5>();
    //     std::cout << type_name2<decltype(flags)>() << std::endl;

    //     auto new_flags = switch_flag<3>(flags);
    //     std::cout << type_name2<decltype(new_flags)>() << std::endl;

    //     auto new_flags2 = switch_flag<2>(new_flags);
    //     std::cout << type_name2<decltype(new_flags2)>() << std::endl;

    //     auto new_flags3 = switch_flag<0>(new_flags2);
    //     std::cout << type_name2<decltype(new_flags3)>() << std::endl;

    //     auto new_flags4 = switch_flag<4>(new_flags3);
    //     std::cout << type_name2<decltype(new_flags4)>() << std::endl;
    // }
}

} // namespace adhoc4
