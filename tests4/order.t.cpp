
#include <adhoc.hpp>
#include <constants_type.hpp>
#include <order.hpp>

#include <differential_operator.hpp>

#include <gtest/gtest.h>

#include "type_name.hpp"

namespace adhoc4 {

TEST(SaturationArithmetic, Inverse) {

    constexpr auto max = std::numeric_limits<std::size_t>::max();
    constexpr auto min = std::numeric_limits<std::size_t>::min();

    constexpr auto a = detail::add_sat(static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(4));
    static_assert(a == 7);

    constexpr auto b = detail::add_sat(max, static_cast<std::size_t>(4));
    static_assert(b == max);

    constexpr auto c = detail::add_sat(max, max, max);
    static_assert(c == max);

    constexpr auto d = detail::add_sat(static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(3));
    static_assert(d == 9);
}

TEST(Order, Const) {
    using constants::CD;

    ADHOC(var);
    auto constant = CD<-0.5>();

    // static_assert(order<decltype(constant), decltype(var)>::call() == 0);
    static_assert(order(constant, var) == 0);
}

TEST(Order, VarType) {
    ADHOC(var);
    ADHOC(var2);

    // static_assert(order<decltype(var), decltype(var)>::call() == 1);
    // static_assert(order<decltype(var), decltype(var2)>::call() == 0);
    static_assert(order(var, var) == 1);
    static_assert(order(var, var2) == 0);
}

TEST(Order, UnivariateType) {
    ADHOC(var);

    auto var2 = exp(var);
    // static_assert(order<decltype(var2), decltype(var)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var2, var) == std::numeric_limits<std::size_t>::max());
}

TEST(Order, Add) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 + var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);
}

TEST(Order, Sub) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 - var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);
}

TEST(Order, Mul) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 * var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);

    auto var4 = var1 * var1;
    // static_assert(order<decltype(var4), decltype(var1)>::call() == 2);
    static_assert(order(var4, var1) == 2);
}

TEST(Order, Div) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 / var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var2) == std::numeric_limits<std::size_t>::max());
}

TEST(Order, First) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto var = exp(T);
    // static_assert(order<decltype(var), decltype(T)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var, T) == std::numeric_limits<std::size_t>::max());

    auto sumv = K + v;
    // static_assert(order<decltype(sumv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(sumv), decltype(K)>::call() == 1);
    // static_assert(order<decltype(sumv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(sumv), decltype(T)>::call() == 0);
    static_assert(order(sumv, S) == 0);
    static_assert(order(sumv, K) == 1);
    static_assert(order(sumv, v) == 1);
    static_assert(order(sumv, T) == 0);

    auto mulv = K * v * K;
    // static_assert(order<decltype(mulv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(mulv), decltype(K)>::call() == 2);
    // static_assert(order<decltype(mulv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(mulv), decltype(T)>::call() == 0);
    static_assert(order(mulv, S) == 0);
    static_assert(order(mulv, K) == 2);
    static_assert(order(mulv, v) == 1);
    static_assert(order(mulv, T) == 0);

    auto mulexpv = K * v * exp(K);
    // static_assert(order<decltype(mulexpv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(mulexpv), decltype(K)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    // static_assert(order<decltype(mulexpv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(mulexpv), decltype(T)>::call() == 0);
    static_assert(order(mulexpv, S) == 0);
    static_assert(order(mulexpv, K) == std::numeric_limits<std::size_t>::max());
    static_assert(order(mulexpv, v) == 1);
    static_assert(order(mulexpv, T) == 0);

    auto divv = v / K;
    // static_assert(order<decltype(divv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(divv), decltype(K)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    // static_assert(order<decltype(divv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(divv), decltype(T)>::call() == 0);
    static_assert(order(divv, S) == 0);
    static_assert(order(divv, K) == std::numeric_limits<std::size_t>::max());
    static_assert(order(divv, v) == 1);
    static_assert(order(divv, T) == 0);

    using constants::CD;

    auto poly = K * (K + CD<3.>()) +
                (K + CD<0.5>()) * (K + CD<1.5>()) * (K + CD<2.5>());
    // static_assert(order<decltype(poly), decltype(K)>::call() == 3);
    static_assert(order(poly, K) == 3);
}

TEST(Order, AddSub) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 + var2;
    auto var4 = var3 - var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var4, var1) == 1);
    static_assert(order(var4, var2) == 1);
}

template <detail::StringLiteral literal1, detail::StringLiteral literal2>
constexpr auto operator==(double_t<literal1> /* lhs */,
                          double_t<literal2> /* rhs */) -> bool {
    constexpr auto var1 = std::to_array(literal1.value);
    constexpr auto var2 = std::to_array(literal2.value);

    if constexpr (var1.size() != var2.size()) {
        return false;
    } else {
        return var1 == var2;
    }
}

template <detail::StringLiteral literal1, detail::StringLiteral literal2>
constexpr auto operator!=(double_t<literal1> lhs, double_t<literal2> rhs)
    -> bool {
    return !(lhs == rhs);
}

template <constants::ArgType N, class Denominator>
constexpr auto depends(constants::CD<N> /* num */, Denominator /* den */)
    -> bool {
    // why false? Because a constant dos not depend on any variable
    return false;
};

template <detail::StringLiteral literal1, class Denominator>
constexpr auto depends(double_t<literal1> num, Denominator den) -> bool {
    return num == den;
};

template <template <class> class Univariate, class Input, class Denominator>
constexpr auto depends(Univariate<Input> /* num */, Denominator den) -> bool {
    return depends(Input{}, den);
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Denominator>
constexpr auto depends(Bivariate<Input1, Input2> /* num */, Denominator den)
    -> bool {
    return depends(Input1{}, den) || depends(Input2{}, den);
};

template <constants::ArgType N, class Denominator>
constexpr auto order_var(constants::CD<N> /* num */, Denominator /* den */)
    -> std::size_t {
    // why false? Because a constant dos not depend on any variable
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

template <class Input1, class Input2>
constexpr auto expand(add_t<Input1, Input2> /* num */) {

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

TEST(Order, Vars) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var33);

    static_assert(var1 == var1);
    static_assert(var1 != var2);
    static_assert(var1 != var33);

    static_assert(depends(var1, var1));
    static_assert(!depends(var1, var2));
    static_assert(!depends(var1, var33));

    auto ret1 = exp(var1);
    static_assert(depends(ret1, var1));
    static_assert(!depends(ret1, var2));

    auto ret2 = var1 + var2;
    static_assert(depends(ret2, var1));
    static_assert(depends(ret2, var2));
    static_assert(!depends(ret2, var33));

    auto ret3 = exp(ret2);

    auto result = order_vars(ret3, var1, var2, var33);
    std::cout << type_name2<decltype(result)>() << std::endl;

    auto ret4 = var1 * var2;
    auto result2 = order_vars(var1 * var2, var1, var2, var33);
    // auto result2 = order_vars(ret4, var1);
    std::cout << type_name2<decltype(result2)>() << std::endl;

    auto ret5 = var1 + var2;
    auto result3 = expand(ret5);
    std::cout << type_name2<decltype(result3)>() << std::endl;
}

TEST(Order, Expand) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    ADHOC(var4);

    // {
    //     auto var = var3 * (var1 + var2);
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * var3;
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * (var3 + var4);
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * (var3 + var4) * exp(var1);

    //     auto result2 = order_vars2(var, var3, var2);
    //     std::cout << type_name2<decltype(result2)>() << std::endl;
    // }

    {
        auto var = (var1 + (var3 * var4)) * (var3 + (var3 * var4)) * exp(var1);
        // auto var = var1 + (var1 + var3) * (var3);

        // auto var = var1 + var1 * (var3);

        auto vars = std::tuple(var1, var2, var3, var4);
        auto result2 = convex_hull_orders(var, vars);
        std::cout << type_name2<decltype(result2)>() << std::endl;
    }

    auto v1 = std::index_sequence<18446744073709551615UL, 0, 1, 0>{};
    auto v2 = std::index_sequence<18446744073709551615UL, 0, 1, 1>{};
    auto v3 = std::index_sequence<18446744073709551615UL, 0, 0, 1>{};

    // bool res = v1 <= v2;
    std::cout << (v1 <= v2) << std::endl;
    std::cout << (v2 <= v1) << std::endl;
    std::cout << (v1 <= v1) << std::endl;
    std::cout << (v2 <= v2) << std::endl;

    std::cout << (v1 <= v3) << std::endl;
    std::cout << (v3 <= v1) << std::endl;

    auto vartoreduce = std::make_tuple(v1, v2, v3);
    auto varreduced = filter_smaller(vartoreduce);
    std::cout << type_name2<decltype(vartoreduce)>() << std::endl;
    std::cout << type_name2<decltype(varreduced)>() << std::endl;
}

template <class... Ids, std::size_t... Orders>
constexpr auto
get_vars_single(std::tuple<der::d<Orders, Ids>...> /* diff_ops */) {
    return std::tuple<Ids...>{};
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

template <class... DifferentialOperators>
constexpr auto get_vars(std::tuple<DifferentialOperators...> /* diff_ops */) {
    constexpr auto all_vars =
        std::tuple_cat(get_vars_single(DifferentialOperators{})...);
    return filter_equal(all_vars);
}

TEST(Order, GetVars) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    auto res = d<2>(var1) * d(var2);
    auto res2 = d<2>(var1) * d(var3);

    auto vars = get_vars_single(res);
    auto vars2 = get_vars_single(res2);
    std::cout << type_name2<decltype(vars)>() << std::endl;

    auto res3 = std::tuple(res, res2);
    auto vars3 = get_vars(res3);
    std::cout << type_name2<decltype(vars3)>() << std::endl;
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

template <class Res, class Id, class... Ids, std::size_t Order,
          std::size_t... Orders, class Inputs>
constexpr auto expand_aux(
    Res res,
    std::tuple<der::d<Order, Id>, der::d<Orders, Ids>...> /* diff_op_out */,
    Inputs inputs) {
    constexpr auto new_res =
        std::tuple_cat(res, expand_aux2(convex_hull_orders(Id{}, inputs),
                                        std::make_index_sequence<Order>()));

    return expand_aux(new_res, std::tuple<der::d<Orders, Ids>...>{}, inputs);
}

template <class... Ids, std::size_t... Orders, class Inputs>
constexpr auto expand(std::tuple<der::d<Orders, Ids>...> diff_op_out,
                      Inputs inputs) {
    return expand_aux(std::tuple<>{}, diff_op_out, inputs);
}

template <class Res>
constexpr auto flatten_aux(Res res, std::tuple<> /* idx_seqs */) {
    return res;
}

template <class TupleRes, class... IndexSequence>
constexpr auto flatten_aux3(TupleRes in1,
                            std::tuple<IndexSequence...> /* in2 */) {
    return std::make_tuple(
        std::tuple_cat(in1, std::make_tuple(IndexSequence{}))...);
}

template <class... TupleRes, class... IndexSequence>
constexpr auto flatten_aux2(std::tuple<TupleRes...> /* res */,
                            std::tuple<IndexSequence...> newinputs) {
    return std::tuple_cat(flatten_aux3(TupleRes{}, newinputs)...);
}

template <class Res, class... IndexSequences, class... TupleIndexSequences>
constexpr auto flatten_aux(Res res,
                           std::tuple<std::tuple<IndexSequences...>,
                                      TupleIndexSequences...> /* idx_seqs */) {
    if constexpr (std::tuple_size_v<Res> == 0) {
        constexpr auto new_res =
            std::make_tuple(std::tuple<IndexSequences...>{});
        return flatten_aux(new_res, std::tuple<TupleIndexSequences...>{});
    } else {
        constexpr auto new_res =
            flatten_aux2(res, std::tuple<IndexSequences...>{});
        return flatten_aux(new_res, std::tuple<TupleIndexSequences...>{});
    }
}

template <class... IndexSequences, class... TupleIndexSequences>
constexpr auto
flatten(std::tuple<std::tuple<IndexSequences...>, TupleIndexSequences...>
            idx_seqs) {
    return flatten_aux(std::tuple<>{}, idx_seqs);
}

template <class Input>
constexpr auto order_input_aux2(std::tuple<> /* diff_op */, Input /* input */)
    -> std::size_t {
    return 0;
}

template <class Id, class... Ids, std::size_t Order, std::size_t... Orders,
          class Input>
constexpr auto order_input_aux2(
    std::tuple<der::d<Order, Id>, der::d<Orders, Ids>...> /* diff_op */,
    Input input) -> std::size_t {
    if constexpr (Id{} == input) {
        return Order;
    } else {
        return order_input_aux2(std::tuple<der::d<Orders, Ids>...>{}, input);
    }
}

template <class... Ids, std::size_t... Orders, class... Inputs>
constexpr auto order_input(std::tuple<der::d<Orders, Ids>...> diff_op,
                           std::tuple<Inputs...> /* input */) -> auto {
    return std::index_sequence<order_input_aux2(diff_op, Inputs{})...>{};
}

template <class... DiffOpInputs, class... Inputs>
constexpr auto order_inputs(std::tuple<DiffOpInputs...> /* diff_op_inputs */,
                            std::tuple<Inputs...> input) -> auto {
    return std::make_tuple(order_input(DiffOpInputs{}, input)...);
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
    return element_wise_sum_aux(std::index_sequence<I + IS...>{},
                                std::tuple<ISs...>{});
}

template <class... IS>
constexpr auto element_wise_sum(std::tuple<> /* idx_seqs */) {
    return std::make_index_sequence<0>{};
}

template <class First, class... IS>
constexpr auto element_wise_sum(std::tuple<First, IS...> idx_seqs) {
    return element_wise_sum_aux(AllZeroIS<First{}.size()>(), idx_seqs);
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

template <class Res>
constexpr auto get_minimum_polys_aux(Res res, std::tuple<> /* poly */) {
    return res;
}

template <class Res, class IndexSequence, class... IndexSequences>
constexpr auto get_minimum_polys_aux(
    Res res,
    std::tuple<IndexSequence, IndexSequences...> /* diff_op_output */) {

    constexpr auto separated_IS =
        std::make_tuple(get_minimum_polys_single(IndexSequence{}));

    return res;
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
constexpr auto compare_aginst_input(std::tuple<IndexSequences...> /* poly */,
                                    IndexSequence diff_op_input) {

    return ((IndexSequences{} <= diff_op_input) || ...);
}

TEST(Order, MinimumPolys) {
    std::tuple<std::index_sequence<1, 2, 0>, std::index_sequence<2, 0, 0>,
               std::index_sequence<0, 1, 1>>
        val;
    auto result_exp = get_minimum_polys(val);
    std::cout << type_name2<decltype(result_exp)>() << std::endl;
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
        auto constexpr element_sum = element_wise_sum(poly);
        bool constexpr smallerthan = diff_op_input <= element_sum;
        // if any of the monomial order is larger than all polynomials, it can't
        // be included
        if constexpr (!smallerthan) {
            return false;
        } else {
            // std::cout << type_name2<decltype(poly)>() << std::endl;
            constexpr auto min_poly = get_minimum_polys(poly);
            // std::cout << type_name2<decltype(min_poly)>() << std::endl;
            // std::cout << type_name2<decltype(diff_op_input)>() << std::endl;

            return compare_aginst_input(min_poly, diff_op_input);
            // analyse all polynomials have at least one
            // return true;
        }
    }
}

template <class TupleIndexSequence, class... IndexSequences>
constexpr auto check_included_single_output(TupleIndexSequence /* poly */,
                                            std::tuple<> /* diff_op_input */)
    -> bool {
    return false;
}

template <class TupleIndexSequence, class IndexSequence,
          class... IndexSequences>
constexpr auto check_included_single_output(
    TupleIndexSequence poly,
    std::tuple<IndexSequence, IndexSequences...> /* diff_op_input */) -> bool {
    constexpr bool thisres =
        check_included_final(TupleIndexSequence{}, IndexSequence{});
    if constexpr (thisres) {
        return true;
    } else {
        return check_included_single_output(poly,
                                            std::tuple<IndexSequences...>{});
    }
}

template <class... IndexSequences>
constexpr auto check_included(std::tuple<> /* poly */,
                              std::tuple<IndexSequences...> /* diff_op_input */)
    -> bool {
    return false;
}

template <class TupleIndexSequence, class... TupleIndexSequences,
          class... IndexSequences>
constexpr auto check_included(
    std::tuple<TupleIndexSequence, TupleIndexSequences...> /* poly */,
    std::tuple<IndexSequences...> diff_op_input) -> bool {
    constexpr bool thisres =
        check_included_single_output(TupleIndexSequence{}, diff_op_input);
    if constexpr (thisres) {
        return true;
    } else {
        return check_included(std::tuple<TupleIndexSequences...>{},
                              diff_op_input);
    }
}

template <class... Ids, std::size_t... Orders, class... DiffOpInputs>
constexpr auto monomial_included(std::tuple<der::d<Orders, Ids>...> diff_op_out,
                                 std::tuple<DiffOpInputs...> diff_op_inputs)
    -> bool {
    constexpr auto inputs = get_vars(diff_op_inputs);
    constexpr auto orders = expand(diff_op_out, inputs);
    constexpr auto flat_orders_output = flatten(orders);

    constexpr auto order_inputs_vals = order_inputs(diff_op_inputs, inputs);
    return check_included(flat_orders_output, order_inputs_vals);
}

TEST(Order, MonomialIncluded) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    ADHOC(var4);
    auto res = d<2>(var1) * d(var2);
    auto res2 = d<2>(var1) * d(var3);

    auto res3 = std::tuple(res, res2);

    // auto var = d((var1 + (var3 * var4)) * (var3 + (var3 * var4)) *
    // exp(var1));
    // auto var = d(var1 + (var3 * var4)) * d<2>(var1);
    auto var = d<2>(var1) * d(var1 + (var3 * var4));

    auto resbool = monomial_included(var, res3);
    std::cout << resbool << std::endl;
    std::cout << type_name2<decltype(resbool)>() << std::endl;
}

TEST(Order, MonomialNotIncluded) {
    ADHOC(var1);
    ADHOC(var2);

    auto res = d<2>(var1);
    auto res2 = d<2>(var2);
    auto res_t = std::tuple(res, res2);
    auto out1 = d(exp(var1)) * d(exp(var2));
    constexpr auto resbool = monomial_included(out1, res_t);
    static_assert(!resbool);

    ADHOC(var3);

    auto res3 = d<2>(var3) * d<2>(var1);
    auto res2_t = std::tuple(res2, res3);
    auto out2 = d(exp(var1)) * d(exp(var3));
    constexpr auto resbool2 = monomial_included(out2, res2_t);
    static_assert(resbool2);
}

} // namespace adhoc4
