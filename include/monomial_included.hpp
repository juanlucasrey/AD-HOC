/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MONOMIAL_INCLUDED_HPP
#define MONOMIAL_INCLUDED_HPP

#include "add_sat.hpp"
#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "differential_operator.hpp"
#include "utils/flatten.hpp"

#include <tuple>

namespace adhoc4 {

template <constants::ArgType N, class Denominator>
constexpr auto order_var(constants::CD<N> /* num */, Denominator /* den */)
    -> std::size_t {
    // why 0? Because a constant dos not depend on any variable
    return 0;
};

template <detail::StringLiteral literal1, class Denominator>
constexpr auto order_var(double_t<literal1> num, Denominator den)
    -> std::size_t {
    return num == den ? 1U : 0U;
};

template <template <class> class Univariate, class Input, class Denominator>
constexpr auto order_var(Univariate<Input> /* num */, Denominator den)
    -> std::size_t {
    if constexpr (depends(Input{}, den)) {
        return std::numeric_limits<std::size_t>::max();
    } else {
        return 0;
    }
}

template <class Input1, class Input2, class Denominator>
constexpr auto order_var(mul_t<Input1, Input2> /* num */, Denominator den)
    -> std::size_t {
    return detail::add_sat(order_var(Input1{}, den), order_var(Input2{}, den));
}

template <class Numerator, class... Denominators>
constexpr auto order_vars(Numerator num, Denominators... den) {
    return std::index_sequence<order_var(num, den)...>{};
}

template <constants::ArgType N> constexpr auto expand(constants::CD<N> num) {
    return num;
};

template <detail::StringLiteral literal1>
constexpr auto expand(double_t<literal1> num) {
    return num;
};

template <template <class> class Univariate, class Input>
constexpr auto expand(Univariate<Input> num) {
    return num;
}

// default for add_t and sub_t
template <template <class, class> class Bivariate, class Input1, class Input2>
constexpr auto expand(Bivariate<Input1, Input2> /* num */) {

    constexpr auto expanded1 = expand(Input1{});
    constexpr auto expanded2 = expand(Input2{});
    constexpr bool same1 =
        std::is_same_v<const decltype(expanded1), const Input1>;
    constexpr bool same2 =
        std::is_same_v<const decltype(expanded2), const Input2>;

    if constexpr (same1 && same2) {
        return std::make_tuple(expanded1, expanded2);
    } else if constexpr (same1 && !same2) {
        return std::tuple_cat(std::make_tuple(expanded1), expanded2);
    } else if constexpr (!same1 && same2) {
        return std::tuple_cat(expanded1, std::make_tuple(expanded2));
    } else {
        return std::tuple_cat(std::make_tuple(expanded1),
                              std::make_tuple(expanded2));
    }
}

template <class Input1, class... Input2>
constexpr auto mult_aux1(Input1 in1, std::tuple<Input2...> /* in2 */) {
    return std::tuple{(in1 * Input2{})...};
}

template <class... Input1, class Input2>
constexpr auto mult_aux2(std::tuple<Input1...> /* in1 */, Input2 in2) {
    return std::tuple{(Input1{} * in2)...};
}

template <class... Input1, class... Input2>
constexpr auto mult_aux3(std::tuple<Input1...> /* in1 */,
                         std::tuple<Input2...> in2) {
    return std::tuple_cat(mult_aux1(Input1{}, in2)...);
}

template <class Input1, class Input2>
constexpr auto expand(mul_t<Input1, Input2> num) {
    constexpr auto expanded1 = expand(Input1{});
    constexpr auto expanded2 = expand(Input2{});
    constexpr bool same1 =
        std::is_same_v<const decltype(expanded1), const Input1>;
    constexpr bool same2 =
        std::is_same_v<const decltype(expanded2), const Input2>;
    if constexpr (same1 && same2) {
        return num;
    } else if constexpr (same1 && !same2) {
        return mult_aux1(Input1{}, expanded2);
    } else if constexpr (!same1 && same2) {
        return mult_aux2(expanded1, Input2{});
    } else {
        return mult_aux3(expanded1, expanded2);
    }
}

template <class... Numerators, class... Denominators>
constexpr auto order_vars2_aux(std::tuple<Numerators...> /* num */,
                               Denominators... den) {
    return std::tuple{order_vars(Numerators{}, den...)...};
}

template <class Numerator, class... Denominators>
constexpr auto order_vars2_aux(Numerator num, Denominators... den) {
    return std::tuple{order_vars(num, den...)};
}

constexpr auto operator<=(std::index_sequence<> /* lhs */,
                          std::index_sequence<> /* rhs */) -> bool {
    return true;
}

template <std::size_t First1, std::size_t... Seq1, std::size_t First2,
          std::size_t... Seq2>
constexpr auto operator<=(std::index_sequence<First1, Seq1...> /* lhs */,
                          std::index_sequence<First2, Seq2...> /* rhs */)
    -> bool {
    return (First1 <= First2) &&
           (std::index_sequence<Seq1...>{} <= std::index_sequence<Seq2...>{});
}

template <class IS>
constexpr auto smaller(IS /* seq */, std::tuple<> /* seqs */) -> bool {
    return false;
}

template <class IS, class... ISs>
constexpr auto smaller(IS seq, std::tuple<ISs...> /* seqs */) -> bool {
    return ((seq <= ISs{}) || ...);
}

template <class Res>
constexpr auto filter_smaller_aux(Res res, std::tuple<> /* seqs */) {
    return res;
}

template <class Res, class IS, class... ISs>
constexpr auto filter_smaller_aux(Res res, std::tuple<IS, ISs...> /* seqs */) {
    constexpr bool smaller1 = smaller(IS{}, std::tuple<ISs...>{});
    constexpr bool smaller2 = smaller(IS{}, res);
    if constexpr (smaller1 || smaller2) {
        return filter_smaller_aux(res, std::tuple<ISs...>{});
    } else {
        constexpr auto new_result = std::tuple_cat(res, std::make_tuple(IS{}));
        return filter_smaller_aux(new_result, std::tuple<ISs...>{});
    }
}

template <class In> constexpr auto filter_smaller(In seqs) {
    return filter_smaller_aux(std::tuple<>{}, seqs);
}

template <class Numerator, class... Denominators>
constexpr auto convex_hull_orders(Numerator num,
                                  std::tuple<Denominators...> /* den */) {
    constexpr auto expanded = expand(num);
    constexpr auto orders = order_vars2_aux(expanded, Denominators{}...);
    return filter_smaller(orders);
}

template <class Tuple> constexpr auto get_vars_single(Tuple diff_ops) {
    return std::apply(
        [](auto... dop) { return std::make_tuple(get_id(dop)...); }, diff_ops);
}

template <class IS>
constexpr auto equal(IS /* seq */, std::tuple<> /* seqs */) -> bool {
    return false;
}

template <class IS, class... ISs>
constexpr auto equal(IS seq, std::tuple<ISs...> /* seqs */) -> bool {
    return ((seq == ISs{}) || ...);
}

template <class Res>
constexpr auto filter_equal_aux(Res res, std::tuple<> /* seqs */) {
    return res;
}

template <class Res, class IS, class... ISs>
constexpr auto filter_equal_aux(Res res, std::tuple<IS, ISs...> /* seqs */) {
    if constexpr (equal(IS{}, res)) {
        return filter_equal_aux(res, std::tuple<ISs...>{});
    } else {
        constexpr auto new_result = std::tuple_cat(res, std::make_tuple(IS{}));
        return filter_equal_aux(new_result, std::tuple<ISs...>{});
    }
}

template <class In> constexpr auto filter_equal(In seqs) {
    return filter_equal_aux(std::tuple<>{}, seqs);
}

template <class... DiffOpInputs, class... Positions>
constexpr auto
get_vars(std::tuple<std::pair<DiffOpInputs, Positions>...> /* diff_ops */) {
    constexpr auto all_vars =
        std::tuple_cat(get_vars_single(DiffOpInputs{})...);
    return filter_equal(all_vars);
}

template <class Tuple>
constexpr auto expand_aux3(Tuple tup, std::size_t /* idx */) {
    return tup;
}

template <class Tuple, size_t... Indices>
constexpr auto expand_aux2(Tuple tup,
                           std::index_sequence<Indices...> /* idx */) {
    return std::make_tuple(expand_aux3(tup, Indices)...);
}

template <class Res, class Inputs>
constexpr auto expand_aux(Res res, std::tuple<> /* seqs */,
                          Inputs /* inputs */) {
    return res;
}

template <class Res, class DiffOpOut, class Inputs>
constexpr auto expand_aux(Res res, DiffOpOut diff_op_out, Inputs inputs) {
    constexpr auto first = head(diff_op_out);
    constexpr auto last = tail(diff_op_out);
    constexpr auto new_res = std::tuple_cat(
        res, expand_aux2(convex_hull_orders(get_id(first), inputs),
                         std::make_index_sequence<get_power(first)>()));

    return expand_aux(new_res, last, inputs);
}

template <class DiffOptOut, class Inputs>
constexpr auto expand(DiffOptOut diff_op_out, Inputs inputs) {
    return expand_aux(std::tuple<>{}, diff_op_out, inputs);
}

template <std::size_t N, std::size_t... I>
constexpr auto AllZeroIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto first = std::array<std::size_t, N>{0};
    return std::index_sequence<first[I]...>{};
}

template <std::size_t N> constexpr auto AllZeroIS() {
    return AllZeroIS_aux<N>(std::make_index_sequence<N>{});
}

template <std::size_t... I>
constexpr auto element_wise_sum_aux(std::index_sequence<I...> ret,
                                    std::tuple<> /* idx_seqs */) {
    return ret;
}

template <std::size_t... I, std::size_t... IS, class... ISs>
constexpr auto element_wise_sum_aux(
    std::index_sequence<I...> /* ret */,
    std::tuple<std::index_sequence<IS...>, ISs...> /* idx_seqs */) {
    return element_wise_sum_aux(
        std::index_sequence<detail::add_sat(I, IS)...>{}, std::tuple<ISs...>{});
}

template <class... IS>
constexpr auto element_wise_sum(std::tuple<> /* idx_seqs */) {
    return std::make_index_sequence<0>{};
}

template <class First, class... IS>
constexpr auto element_wise_sum(std::tuple<First, IS...> idx_seqs) {
    return element_wise_sum_aux(AllZeroIS<First{}.size()>(), idx_seqs);
}

template <class Input>
constexpr auto order_input_aux2(std::tuple<> /* diff_op */, Input /* input */)
    -> std::size_t {
    return 0;
}

template <class DiffOp, class Input>
constexpr auto order_input_aux2(DiffOp diff_op, Input input) -> std::size_t {
    constexpr auto first = head(diff_op);
    constexpr auto last = tail(diff_op);

    if constexpr (get_id(first) == input) {
        return get_power(first);
    } else {
        return order_input_aux2(last, input);
    }
}

template <class DiffOp, class... Inputs>
constexpr auto order_input(DiffOp diff_op, std::tuple<Inputs...> /* input */)
    -> auto {
    return std::index_sequence<order_input_aux2(diff_op, Inputs{})...>{};
}

template <class... DiffOpInputs, class... Positions, class... Inputs>
constexpr auto order_inputs(
    std::tuple<std::pair<DiffOpInputs, Positions>...> /* diff_op_inputs */,
    std::tuple<Inputs...> input) -> auto {
    return std::make_tuple(order_input(DiffOpInputs{}, input)...);
}

template <std::size_t Size, std::size_t Pos>
constexpr auto make_tuple_one_aux_arr() {
    auto result = std::array<std::size_t, Size>{0};
    result[Pos] = 1;
    return result;
}

template <std::size_t Pos, std::size_t... I>
constexpr auto make_tuple_one_aux(std::index_sequence<I...> /* i */) {
    constexpr auto result = make_tuple_one_aux_arr<sizeof...(I), Pos>();
    return std::index_sequence<result[I]...>{};
}

template <std::size_t Size, std::size_t Pos> constexpr auto make_tuple_one() {
    return make_tuple_one_aux<Pos>(std::make_index_sequence<Size>{});
}

template <bool Cond, std::size_t Size, std::size_t Pos>
constexpr auto conditional_IS() {
    if constexpr (Cond) {
        return std::tuple<>{};
    } else {
        return std::make_tuple(make_tuple_one<Size, Pos>());
    }
}

template <std::size_t... I, std::size_t... Pos>
constexpr auto
get_minimum_polys_single(std::index_sequence<I...> /* diff_op_output */,
                         std::index_sequence<Pos...> /* position */) {
    constexpr std::size_t ISsize = sizeof...(I);
    return std::tuple_cat(conditional_IS<I == 0, ISsize, Pos>()...);
}

template <std::size_t... I>
constexpr auto
get_minimum_polys_single(std::index_sequence<I...> diff_op_output) {
    constexpr std::size_t ISsize = sizeof...(I);
    return get_minimum_polys_single(diff_op_output,
                                    std::make_index_sequence<ISsize>{});
}

template <class TupleRes, class... IndexSequence>
constexpr auto flatten2_aux3(TupleRes in1,
                             std::tuple<IndexSequence...> /* in2 */) {
    auto result = std::make_tuple(
        std::tuple_cat(in1, std::make_tuple(IndexSequence{}))...);
    return result;
}

template <class... TupleRes, class... IndexSequence>
constexpr auto flatten2_aux2(std::tuple<TupleRes...> /* res */,
                             std::tuple<IndexSequence...> newinputs) {
    return std::tuple_cat(flatten2_aux3(TupleRes{}, newinputs)...);
}

template <class Res>
constexpr auto flatten2_aux(Res res, std::tuple<> /* idx_seqs */) {
    return res;
}

template <class Res, class... IndexSequences, class... TupleIndexSequences>
constexpr auto flatten2_aux(Res res,
                            std::tuple<std::tuple<IndexSequences...>,
                                       TupleIndexSequences...> /* idx_seqs */) {
    if constexpr (std::tuple_size_v<Res> == 0) {
        constexpr auto new_res =
            std::make_tuple(std::make_tuple(IndexSequences{})...);
        return flatten2_aux(new_res, std::tuple<TupleIndexSequences...>{});
    } else {
        constexpr auto new_res =
            flatten2_aux2(res, std::tuple<IndexSequences...>{});
        return flatten2_aux(new_res, std::tuple<TupleIndexSequences...>{});
    }
}

template <class... IndexSequences, class... TupleIndexSequences>
constexpr auto
flatten2(std::tuple<std::tuple<IndexSequences...>, TupleIndexSequences...>
             idx_seqs) {
    return flatten2_aux(std::tuple<>{}, idx_seqs);
}

template <class... TupleIndexSequences>
constexpr auto element_wise_sum_mult(std::tuple<TupleIndexSequences...>
                                     /* idx_seqs */) {
    return std::make_tuple(element_wise_sum(TupleIndexSequences{})...);
}

template <class... IndexSequences>
constexpr auto get_minimum_polys(std::tuple<IndexSequences...> /* poly */) {
    constexpr auto to_flatten = std::tuple_cat(
        std::make_tuple(get_minimum_polys_single(IndexSequences{}))...);

    constexpr auto flattenned = flatten2(to_flatten);
    return element_wise_sum_mult(flattenned);
}

template <class... IndexSequences, class IndexSequence>
constexpr auto compare_against_input(std::tuple<IndexSequences...> /* poly */,
                                     IndexSequence diff_op_input) {

    return ((IndexSequences{} <= diff_op_input) || ...);
}

template <class... IndexSequences>
constexpr auto any_sequence_is_zero(std::tuple<IndexSequences...> /* poly */)
    -> bool {
    return (all_equal_to<0U>(IndexSequences{}) || ...);
}

template <class TupleIndexSequence, class IndexSequence>
constexpr auto check_included_final(TupleIndexSequence poly,
                                    IndexSequence diff_op_input) -> bool {
    constexpr std::size_t order_monomial = sum(diff_op_input);
    constexpr std::size_t size_poly = std::tuple_size_v<TupleIndexSequence>;

    // if polynomial minimal order is larger than monomial it can't be included
    if constexpr (order_monomial < size_poly) {
        return false;
    } else {
        constexpr auto element_sum = element_wise_sum(poly);
        constexpr bool smallerthan = diff_op_input <= element_sum;
        // if any of the monomial order is larger than all polynomials, it can't
        // be included
        if constexpr (!smallerthan) {
            return false;
        } else {
            if constexpr (any_sequence_is_zero(poly)) {
                // is any sequence has sum 0 it means there is a monomial
                // variable outside the inputs, so it can't be included
                return false;
            } else {
                constexpr auto min_poly = get_minimum_polys(poly);

                return compare_against_input(min_poly, diff_op_input);
                // analyse all polynomials have at least one
            }
        }
    }
}

template <class TupleIndexSequence, class... IndexSequences>
constexpr auto
check_included_single_output(TupleIndexSequence poly,
                             std::tuple<IndexSequences...> /* diff_op_input */)
    -> bool {

    return (check_included_final(poly, IndexSequences{}) || ...);
}

template <class... TupleIndexSequences, class... IndexSequences>
constexpr auto check_included(std::tuple<TupleIndexSequences...> /* poly */,
                              std::tuple<IndexSequences...> diff_op_input)
    -> bool {

    return (
        check_included_single_output(TupleIndexSequences{}, diff_op_input) ||
        ...);
}

template <class Single, class Many>
constexpr auto monomial_included(Single diff_op_out, Many diff_op_inputs)
    -> bool {
    constexpr auto inputs = get_vars(diff_op_inputs);
    constexpr auto orders = expand(diff_op_out, inputs);
    constexpr auto flat_orders_output = flatten(orders);

    constexpr auto order_inputs_vals = order_inputs(diff_op_inputs, inputs);
    return check_included(flat_orders_output, order_inputs_vals);
}

} // namespace adhoc4

#endif // MONOMIAL_INCLUDED_HPP
