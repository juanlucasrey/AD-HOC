// #include <dependency.hpp>
#include <adhoc.hpp>
#include <ders.hpp>
#include <init.hpp>
#include <tape.hpp>
#include <tuple_utils.hpp>
// #include <strict_order.hpp>

#include "../public/init.hpp"
#include "../public/tape.hpp"
#include "type_name.hpp"

#include <partition/binomial_coefficient.hpp>
#include <partition/combinations.hpp>
#include <partition/integer_partition.hpp>
#include <partition/multinomial_coefficient.hpp>
#include <partition/partition_function.hpp>

#include <gtest/gtest.h>

#include <numeric>

namespace adhoc3 {

TEST(Raw, First) {
    auto [x, y] = adhoc::Init<2>();
    auto res = (x * y) * (x * cos(y));

    adhoc::Tape t(res);

    t.set(x) = 0.5;
    t.set(y) = 1.5;

    t.evaluate_fwd();
    std::cout << t.val(res) << std::endl;
    std::cout << t.val(x * cos(y)) << std::endl;
    std::cout << t.val(x * y) << std::endl;

    auto [x2, y2] = Init<2>();
    auto res2 = (x2 * y2) * (x2 * cos(y2));

    auto cross1 =
        der_non_null((x2 * y2), var(x)) && der_non_null((x2 * cos(y2)), var(y));
    auto cross2 =
        der_non_null((x2 * y2), var(y)) && der_non_null((x2 * cos(y2)), var(x));
}

// template <class Id, std::size_t Order = 1> struct d {};
// template <class Id, std::size_t Power = 1> struct p {};
// template <class... Ids> struct m {};

// TEST(Raw, dID) {
//     auto [x, y] = Init<2>();
//     m<p<d<decltype(x), 1>, 3>> der1{};
//     m<p<d<decltype(x), 1>, 2>, p<d<decltype(y), 1>, 1>> der2{};
//     m<p<d<decltype(x), 1>, 1>, p<d<decltype(y), 1>, 2>> der3{};
//     m<p<d<decltype(y), 1>, 3>> der4{};
// }

// template <class In, class... TupleInputs>
// inline auto der(std::tuple<TupleInputs...> /* idxs */, In /* in */,
//                 std::array<double, sizeof...(TupleInputs)> vals) -> double &
//                 {
//     static_assert(has_type<In, TupleInputs...>(),
//                   "Only nodes or leafs can have derivatives set or read");
//     // constexpr auto idx = get_idx_first2<In>(std::tuple<TupleInputs...>{});
//     constexpr auto idx = idx_type2<In, TupleInputs...>();
//     return vals[idx];
// }

// template <class... Ders, std::size_t... Orders>
// constexpr auto d_order(m<p<d<Ders, 1>, Orders>...> /* in */) -> std::size_t {
//     return (Orders + ...);
// }

// namespace detail {

// template <typename... OutputsAndDerivativesAlive> struct outputs_t {
//     template <typename... Outputs> constexpr static auto call() noexcept;
// };

// template <> struct outputs_t<> {
//     template <typename... Outputs> constexpr static auto call() noexcept {
//         return std::tuple<Outputs...>{};
//     }
// };

// template <typename... D, typename... OutputsAndDerivativesAlive>
// struct outputs_t<m<D...>, OutputsAndDerivativesAlive...> {
//     template <typename... Outputs> constexpr static auto call() noexcept {
//         // we don't add anything because a derivative is not an output
//         return outputs_t<OutputsAndDerivativesAlive...>::template call<
//             Outputs...>();
//     }
// };

// template <class D, typename... OutputsAndDerivativesAlive>
// struct outputs_t<D, OutputsAndDerivativesAlive...> {
//     template <typename... Outputs> constexpr static auto call() noexcept {
//         return outputs_t<OutputsAndDerivativesAlive...>::template call<
//             Outputs..., D>();
//     }
// };

// template <typename... OutputsAndDerivativesAlive> struct derivatives_t {
//     template <typename... Derivatives> constexpr static auto call() noexcept;
// };

// template <> struct derivatives_t<> {
//     template <typename... Derivatives> constexpr static auto call() noexcept
//     {
//         return std::tuple<Derivatives...>{};
//     }
// };

// template <typename... D, typename... OutputsAndDerivativesAlive>
// struct derivatives_t<m<D...>, OutputsAndDerivativesAlive...> {
//     template <typename... Derivatives> constexpr static auto call() noexcept
//     {
//         return derivatives_t<OutputsAndDerivativesAlive...>::template call<
//             Derivatives..., m<D...>>();
//     }
// };

// template <class D, typename... OutputsAndDerivativesAlive>
// struct derivatives_t<D, OutputsAndDerivativesAlive...> {
//     template <typename... Derivatives> constexpr static auto call() noexcept
//     {
//         // we don't add anything because an output is not a derivative
//         return derivatives_t<OutputsAndDerivativesAlive...>::template call<
//             Derivatives...>();
//     }
// };

// template <typename... NodesAlive> struct inputs_t {
//     template <typename... Inputs> constexpr static auto call() noexcept;
// };

// template <> struct inputs_t<> {
//     template <typename... Inputs> constexpr static auto call() noexcept {
//         return std::tuple<Inputs...>{};
//     }
// };

// template <constants::ArgType D, typename... NodesAlive>
// struct inputs_t<constants::CD<D>, NodesAlive...> {
//     template <typename... Inputs> constexpr static auto call() noexcept {
//         // we don't add anything because a constant is not an input
//         return inputs_t<NodesAlive...>::template call<Inputs...>();
//     }
// };

// template <std::size_t N, typename... NodesAlive>
// struct inputs_t<double_t<N>, NodesAlive...> {
//     template <typename... Inputs> constexpr static auto call() noexcept {
//         using this_type = double_t<N>;

//         constexpr bool already_included = has_type<this_type, Inputs...>();

//         if constexpr (already_included) {
//             // It's already been added, so we don't add the leaf
//             return inputs_t<NodesAlive...>::template call<Inputs...>();
//         } else {
//             // we add the leaf
//             return inputs_t<NodesAlive...>::template call<Inputs...,
//                                                           this_type>();
//         }
//     }
// };

// template <template <class...> class Xvariate, class... Node,
//           typename... NodesAlive>
// struct inputs_t<Xvariate<Node...>, NodesAlive...> {
//     template <typename... Inputs> constexpr static auto call() noexcept {
//         return inputs_t<Node..., NodesAlive...>::template call<Inputs...>();
//     }
// };

// template <typename... Derivatives> struct orders_t {
//     template <std::size_t... Orders> constexpr static auto call() noexcept;
// };

// template <> struct orders_t<> {
//     template <std::size_t... Orders> constexpr static auto call() noexcept {
//         return std::index_sequence<Orders...>{};
//     }
// };

// template <typename D, typename... Derivatives>
// struct orders_t<D, Derivatives...> {
//     template <std::size_t... Orders> constexpr static auto call() noexcept {
//         auto constexpr this_order = d_order(D{});

//         constexpr bool already_included = ((this_order == Orders) || ...);

//         if constexpr (already_included) {
//             // It's already been added, so we don't add the new order
//             return orders_t<Derivatives...>::template call<Orders...>();
//         } else {
//             // we add the order
//             return orders_t<Derivatives...>::template call<Orders...,
//                                                            this_order>();
//         }
//     }
// };

// template <class Output, std::size_t... Orders> constexpr auto expand() {
//     return std::tuple<d<Output, Orders>...>{};
// }

// } // namespace detail

// template <class... OutputsAndDerivatives>
// constexpr auto Outputs(OutputsAndDerivatives... /* o */) {
//     return detail::outputs_t<OutputsAndDerivatives...>::template call<>();
// }

// template <class... OutputsAndDerivatives>
// constexpr auto Derivatives(OutputsAndDerivatives... /* o */) {
//     return detail::derivatives_t<OutputsAndDerivatives...>::template
//     call<>();
// }

// template <class... Outputs>
// constexpr auto Inputs(std::tuple<Outputs...> /* o */) {
//     return detail::inputs_t<Outputs...>::template call<>();
// }

// template <class... Derivatives>
// constexpr auto Orders(std::tuple<Derivatives...> /* o */) {
//     return detail::orders_t<Derivatives...>::template call<>();
// }

// template <class... Outputs, std::size_t... Orders>
// constexpr auto ExpandOutputs(std::tuple<Outputs...> /* o */,
//                              std::index_sequence<Orders...> /* orders */) {
//     return std::tuple_cat((detail::expand<Outputs, Orders...>())...);
// }

// template <class... OutputsAndDerivatives> class Tape2 {
//   private:
//     using outputs =
//         decltype(detail::outputs_t<OutputsAndDerivatives...>::template
//         call());

//     using inputs = decltype(Inputs(outputs{}));

//     using derivatives_inputs =
//         decltype(detail::derivatives_t<
//                  OutputsAndDerivatives...>::template call());

//     using orders = decltype(Orders(derivatives_inputs{}));

//     using derivatives_outputs = decltype(ExpandOutputs(outputs{}, orders{}));

//     using derivatives_outputs_and_inputs =
//         decltype(std::tuple_cat(derivatives_inputs{},
//         derivatives_outputs{}));

//     std::array<double, std::tuple_size_v<derivatives_outputs_and_inputs>>
//         m_derivatives{};

//   public:
//     explicit Tape2(OutputsAndDerivatives... /* out */) {}

//     template <class D> auto inline der(D var) const -> double;
//     template <class D> auto inline der(D var) -> double &;
//     void inline reset_der() {
//         std::fill(std::begin(this->m_derivatives),
//                   std::end(this->m_derivatives), 0.);
//     }
// };

// template <class... OutputsAndDerivatives>
// template <class D>
// auto Tape2<OutputsAndDerivatives...>::der(D /* in */) const -> double {
//     static_assert(has_type<D, OutputsAndDerivatives...>(),
//                   "derivative not on tape");
//     constexpr auto idx =
//         get_idx_first2<D>(std::tuple<OutputsAndDerivatives...>{});
//     return this->m_derivatives[idx];
// }

// template <class... OutputsAndDerivatives>
// template <class D>
// auto Tape2<OutputsAndDerivatives...>::der(D /* in */) -> double & {
//     static_assert(has_type<D, OutputsAndDerivatives...>(),
//                   "derivative not on tape");
//     constexpr auto idx =
//         get_idx_first2<D>(std::tuple<OutputsAndDerivatives...>{});
//     return this->m_derivatives[idx];
// }

TEST(Raw, SeparateInputs) {
    auto [x_t, y_t] = Init<2>();
    auto ly_t = log(y_t);
    auto cx_t = cos(x_t);
    auto lytcx_t = ly_t * cx_t;
    auto slytcx_t = sin(lytcx_t);
    auto clytcx_t = cos(lytcx_t);

    auto d1 = d<3>(x_t);
    auto d2 = d<2>(x_t) * d(y_t);
    auto d3 = d(x_t) * d<2>(y_t);
    auto d4 = d<3>(y_t);

    auto d11 = d(x_t);

    // auto dseed = d<3>(clytcx_t);
    der::d<decltype(clytcx_t), 3> dseed{};

    auto tape = Tape(clytcx_t, slytcx_t, d1, d2, d3, d4, d11);

    std::cout << type_name2<decltype(dseed)>() << std::endl;
    tape.der(dseed) = 1.;
}

TEST(Raw, Der3) {
    auto [x_t, y_t] = Init<2>();
    auto ly_t = log(y_t);
    auto cx_t = cos(x_t);
    auto lytcx_t = ly_t * cx_t;
    auto slytcx_t = sin(lytcx_t);

    // std::tuple<m<p<d<decltype(slytcx_t), 3>, 1>>, m<p<d<decltype(x_t), 1>,
    // 3>>,
    //            m<p<d<decltype(x_t), 1>, 2>, p<d<decltype(y_t), 1>, 1>>,
    //            m<p<d<decltype(x_t), 1>, 1>, p<d<decltype(y_t), 1>, 2>>,
    //            m<p<d<decltype(y_t), 1>, 3>>>
    //     r_t;

    // auto constexpr resx1 = d_order(m<p<d<decltype(x_t), 1>, 3>>{});
    // std::cout << resx1 << std::endl;

    // auto constexpr resx2 =
    //     d_order(m<p<d<decltype(x_t), 1>, 2>, p<d<decltype(y_t), 1>, 1>>{});
    // std::cout << resx2 << std::endl;

    // auto constexpr resx3 =
    //     d_order(m<p<d<decltype(x_t), 1>, 1>, p<d<decltype(y_t), 1>, 2>>{});
    // std::cout << resx3 << std::endl;

    // auto constexpr resx4 = d_order(m<p<d<decltype(y_t), 1>, 3>>{});
    // std::cout << resx4 << std::endl;

    // x3, x2y1, x1y2, y3
    std::array<double, 5> results{};
    results.fill(0);

    results[0] = 1.0;

    // m<p<d<decltype(slytcx_t), 3>, 1>> ind{};
    // std::cout << der(r_t, ind, results) << std::endl;

    // sin( log(y) * cos(x) )
    // y = 0.5 x = 0.3
    double y = 0.5;
    double x = 0.3;
    double ly = std::log(y);
    double cx = std::cos(x);
    double lytcx = ly * cx;
    double slytcx = std::sin(lytcx);

    auto [whatever] = Init<1>();
    auto ders = sin_t<decltype(whatever)>::d2<3>(slytcx, lytcx);

    std::array<std::array<std::size_t, 3>, partition_function(3)> parts{};

    constexpr auto last = FirstPartition<3>();
    constexpr auto first = LastPartition<3>();

    parts[0] = LastPartition<3>();
    parts[1] = PrevPartition(first);
    parts[2] = PrevPartition(parts[1]);

    constexpr auto rnext1 = PrevPartition(first);
    constexpr auto rnext2 = PrevPartition(rnext1);
    static_assert(last == rnext2);

    constexpr std::size_t order1 =
        std::accumulate(first.begin(), first.end(), 0UL) - 1;
    constexpr std::size_t order2 =
        std::accumulate(rnext1.begin(), rnext1.end(), 0UL) - 1;
    constexpr std::size_t order3 =
        std::accumulate(rnext2.begin(), rnext2.end(), 0UL) - 1;

    std::vector<double> buffer(11);

    // std::tuple<m<p<d<decltype(slytcx_t), 3>, 1>>> c1_t;

    // 0: <ln(y)*cos(x), 1, 3>
    // 1: <ln(y)*cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    // 2: <ln(y)*cos(x), 3, 1>
    buffer[0] = BellCoeff(first);
    // buffer[0] = 1.;
    buffer[0] *= ders[order1];
    buffer[1] = BellCoeff(rnext1);
    // buffer[1] = 1.;
    buffer[1] *= ders[order2];
    buffer[2] = BellCoeff(rnext2);
    // buffer[2] = 1.;
    buffer[2] *= ders[order3];

    static_assert(combinations(2, 1) == 2);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 1)>
        multparts1{{{1, 0}, {0, 1}}};

    static_assert(combinations(2, 2) == 3);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
        multparts2{{{2, 0}, {1, 1}, {0, 2}}};

    static_assert(combinations(2, 3) == 4);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 3)>
        multparts3{{{3, 0}, {2, 1}, {1, 2}, {0, 3}}};

    constexpr auto factorials = factorialarr<3>();
    // https://calculus.subwiki.org/wiki/Product_rule_for_higher_derivatives
    // <ln(y) * cos(x), 1, 1> ->
    // cos(x) * <ln(y), 1, 1>
    // ln(y) * < cos(x), 1, 1>
    std::array<double, 2> lncs1{};
    lncs1[0] = MultinomialCoeff<1, 2>(factorials, multparts1[0]);
    lncs1[1] = MultinomialCoeff<1, 2>(factorials, multparts1[1]);
    lncs1.front() *= cx;
    lncs1.back() *= ly;
    // <ln(y) * cos(x), 2, 1> ->
    // cos(x) * <ln(y), 2, 1>,
    // <cos(x), 1, 1> * <ln(y), 1, 1>
    // ln(y) * <cos(x), 2, 1>
    std::array<double, 3> lncs2{};
    lncs2[0] = MultinomialCoeff<2, 2>(factorials, multparts2[0]);
    lncs2[1] = MultinomialCoeff<2, 2>(factorials, multparts2[1]);
    lncs2[2] = MultinomialCoeff<2, 2>(factorials, multparts2[2]);
    lncs2.front() *= cx;
    lncs2.back() *= ly;

    // <ln(y) * cos(x), 3, 1> ->
    // cos(x) * <ln(y), 3, 1>
    // <ln(y), 2, 1> <cos(x), 1, 1>
    // <ln(y), 1, 1> <cos(x), 2, 1>
    // ln(y) * <cos(x), 3, 1>
    std::array<double, 4> lncs3{};
    lncs3[0] = MultinomialCoeff<3, 2>(factorials, multparts3[0]);
    lncs3[1] = MultinomialCoeff<3, 2>(factorials, multparts3[1]);
    lncs3[2] = MultinomialCoeff<3, 2>(factorials, multparts3[2]);
    lncs3[3] = MultinomialCoeff<3, 2>(factorials, multparts3[3]);
    lncs3.front() *= cx;
    lncs3.back() *= ly;

    std::array<double, 4> cxarr{};
    cxarr[0] = 1.;
    cxarr[1] = lncs1[0];
    cxarr[2] = cxarr[1] * lncs1[0];
    cxarr[3] = cxarr[2] * lncs1[0];

    std::array<double, 4> lyarr{};
    lyarr[0] = 1.;
    lyarr[1] = lncs1[1];
    lyarr[2] = lyarr[1] * lncs1[1];
    lyarr[3] = lyarr[2] * lncs1[1];

    // <a * <ln(y), 1, 1> +
    // b * < cos(x), 1, 1>, 3> ->
    std::array<double, combinations(2, 3)> lncs13{};
    lncs13[0] = MultinomialCoeff<3, 2>(factorials, multparts3[0]);
    lncs13[0] *= cxarr[multparts3[0][0]] * lyarr[multparts3[0][1]];
    lncs13[1] = MultinomialCoeff<3, 2>(factorials, multparts3[1]);
    lncs13[1] *= cxarr[multparts3[1][0]] * lyarr[multparts3[1][1]];
    lncs13[2] = MultinomialCoeff<3, 2>(factorials, multparts3[2]);
    lncs13[2] *= cxarr[multparts3[2][0]] * lyarr[multparts3[2][1]];
    lncs13[3] = MultinomialCoeff<3, 2>(factorials, multparts3[3]);
    lncs13[3] = cxarr[multparts3[3][0]] * lyarr[multparts3[3][1]];

    // 0: <cos(x), 1, 3>
    // 1: <ln(y)*cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    // 2: <ln(y)*cos(x), 3, 1>
    // 3: <ln(y), 1, 3>
    // 4: <ln(y), 1, 2> < cos(x), 1, 1>
    // 5: <ln(y), 1, 1> < cos(x), 1, 2>
    buffer[3] = buffer[0] * lncs13[0];
    buffer[4] = buffer[0] * lncs13[1] / 3.; // xxx
    buffer[5] = buffer[0] * lncs13[2] / 3.; // xxx
    buffer[0] *= lncs13[3];

    // 0: <cos(x), 1, 3>
    // 1: <ln(y), 1, 1> <ln(y)*cos(x), 2, 1>
    // 2: <ln(y)*cos(x), 3, 1>
    // 3: <ln(y), 1, 3>
    // 4: <ln(y), 1, 2> < cos(x), 1, 1>
    // 5: <ln(y), 1, 1> < cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    buffer[6] = buffer[1] * lncs1[1];
    buffer[1] *= lncs1[0];

    // <ln(y), 1, 1> <ln(y)*cos(x), 2, 1>
    //->
    // <ln(y), 1, 1> <ln(y), 2, 1>
    // <ln(y), 1, 2> <cos(x), 1, 1>
    // <ln(y), 1, 1> <cos(x), 2, 1>

    // 0: <cos(x), 1, 3>
    // 1: <ln(y), 1, 1> <cos(x), 2, 1>
    // 2: <ln(y)*cos(x), 3, 1>
    // 3: <ln(y), 1, 3>
    // 4: <ln(y), 1, 2> < cos(x), 1, 1>
    // 5: <ln(y), 1, 1> < cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    // 7: <ln(y), 1, 1> <ln(y), 2, 1>
    buffer[7] = buffer[1] * lncs2[0];
    buffer[4] += buffer[1] * lncs2[1] / 3.; // xxx
    buffer[1] *= lncs2[2] / 3.;             // xxx

    // <ln(y) * cos(x), 3, 1> ->
    // cos(x) * <ln(y), 3, 1>
    // <ln(y), 2, 1> <cos(x), 1, 1>
    // <ln(y), 1, 1> <cos(x), 2, 1>
    // ln(y) * <cos(x), 3, 1>

    // 0: <cos(x), 1, 3>
    // 1: <ln(y), 1, 1> <cos(x), 2, 1>
    // 2: <cos(x), 3, 1>
    // 3: <ln(y), 1, 3>
    // 4: <ln(y), 1, 2> <cos(x), 1, 1>
    // 5: <ln(y), 1, 1> <cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    // 7: <ln(y), 1, 1> <ln(y), 2, 1>
    // 8: <ln(y), 3, 1>
    // 9: <ln(y), 2, 1> <cos(x), 1, 1>
    buffer[8] = buffer[2] * lncs3[0];
    buffer[9] = buffer[2] * lncs3[1];
    buffer[1] += buffer[2] * lncs3[2] / 3.; // xxx
    buffer[2] *= lncs3[3];

    // <cos(x), 1, 1> <ln(y)*cos(x), 2, 1>
    //->
    // <ln(y), 2, 1> <cos(x), 1, 1>
    // <ln(y), 1, 1> <cos(x), 1, 2>
    // <cos(x), 1, 1> <cos(x), 2, 1>

    // 0: <cos(x), 1, 3>
    // 1: <ln(y), 1, 1> <cos(x), 2, 1>
    // 2: <cos(x), 3, 1>
    // 3: <ln(y), 1, 3>
    // 4: <ln(y), 1, 2> <cos(x), 1, 1>
    // 5: <ln(y), 1, 1> <cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <cos(x), 2, 1>
    // 7: <ln(y), 1, 1> <ln(y), 2, 1>
    // 8: <ln(y), 3, 1>
    // 9: <ln(y), 2, 1> <cos(x), 1, 1>
    buffer[9] += buffer[6] * lncs2[0];
    buffer[5] += buffer[6] * lncs2[1] / 3; // xxx
    buffer[6] *= lncs2[2];

    auto derslog = log_t<decltype(whatever)>::d2<3>(ly, y);
    // 0: <cos(x), 1, 3>
    // 1: <y, 1, 1> <cos(x), 2, 1>
    // 2: <cos(x), 3, 1>
    // 3: <y, 1, 3> in result
    // 4: <y, 1, 2> <cos(x), 1, 1>
    // 5: <y, 1, 1> <cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <cos(x), 2, 1>
    // 7: <y, 1, 1> <y, 2, 1>
    // 8: <y, 3, 1>
    // 9: <y, 2, 1> <cos(x), 1, 1>
    buffer[1] *= derslog[0];
    results[3] +=
        buffer[3] * derslog[0] * derslog[0] * derslog[0]; // 3 is free now
    buffer[4] *= derslog[0] * derslog[0];
    buffer[5] *= derslog[0];

    results[3] += buffer[7] * derslog[0] * derslog[1]; // 7 is free now
    results[3] += buffer[8] * derslog[2];              // 8 is free now
    buffer[9] *= derslog[1] / 3.;                      // xxx

    auto derscos = cos_t<decltype(whatever)>::d2<3>(cx, x);
    // 0: <cos(x), 1, 3>
    // 1: <y, 1, 1> <cos(x), 2, 1>
    // 2: <cos(x), 3, 1>
    // 4: <y, 1, 2> <cos(x), 1, 1>
    // 5: <y, 1, 1> <cos(x), 1, 2>
    // 6: <cos(x), 1, 1> <cos(x), 2, 1>
    // 9: <y, 2, 1> <cos(x), 1, 1>

    results[0] += buffer[0] * derscos[0] * derscos[0] * derscos[0];
    results[1] += buffer[1] * derscos[1];
    results[0] += buffer[2] * derscos[2];
    results[2] += buffer[4] * derscos[0];
    results[1] += buffer[5] * derscos[0] * derscos[0];
    results[0] += buffer[6] * derscos[0] * derscos[1];
    results[2] += buffer[9] * derscos[0];

    constexpr double third = 1. / 3.;
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << results[0] / 0.081872135144674396 << std::endl;
    std::cout << results[1] / -0.9410337416235649 << std::endl;
    std::cout << results[2] / -1.5274504231928869 << std::endl;
    std::cout << results[3] / -6.9137888891600721 << std::endl;
    // EXPECT_NEAR(results[0], 0.081872135144674396, 1e-14);
    // EXPECT_NEAR(results[1] * third, -0.9410337416235649, 1e-14);
    // EXPECT_NEAR(results[2] * third, -1.5274504231928869, 1e-14);
    // EXPECT_NEAR(results[3], -6.9137888891600721, 1e-14);

    // finite der x3
    double xoriginal = x;
    double yoriginal = y;
    double epsilon = 1e-4;

    x = xoriginal - 2 * epsilon;
    double valuex3m2 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal - epsilon;
    double valuex3m1 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal + epsilon;
    double valuex3p1 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal + 2 * epsilon;
    double valuex3p2 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal;

    double derx3fd =
        (-0.5 * valuex3m2 + valuex3m1 - valuex3p1 + 0.5 * valuex3p2) /
        (epsilon * epsilon * epsilon);
    std::cout << derx3fd << std::endl;

    y = yoriginal - epsilon;
    x = xoriginal - epsilon;
    double valuey1x2v1 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal;
    double valuey1x2v2 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal + epsilon;
    double valuey1x2v3 = std::sin(std::log(y) * std::cos(x));
    double der21 =
        (valuey1x2v1 - 2 * valuey1x2v2 + valuey1x2v3) / (epsilon * epsilon);

    y = yoriginal + epsilon;
    x = xoriginal - epsilon;
    double valuey1x2v4 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal;
    double valuey1x2v5 = std::sin(std::log(y) * std::cos(x));
    x = xoriginal + epsilon;
    double valuey1x2v6 = std::sin(std::log(y) * std::cos(x));
    double der22 =
        (valuey1x2v4 - 2 * valuey1x2v5 + valuey1x2v6) / (epsilon * epsilon);
    double dery1x2 = (-0.5 * der21 + 0.5 * der22) / epsilon;
    y = yoriginal;
    x = xoriginal;
    std::cout << dery1x2 << std::endl;

    x = xoriginal - epsilon;
    y = yoriginal - epsilon;
    double valuey2x1v1 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal;
    double valuey2x1v2 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal + epsilon;
    double valuey2x1v3 = std::sin(std::log(y) * std::cos(x));
    double der11 =
        (valuey2x1v1 - 2 * valuey2x1v2 + valuey2x1v3) / (epsilon * epsilon);

    x = xoriginal + epsilon;
    y = yoriginal - epsilon;
    double valuey2x1v4 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal;
    double valuey2x1v5 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal + epsilon;
    double valuey2x1v6 = std::sin(std::log(y) * std::cos(x));
    double der12 =
        (valuey2x1v4 - 2 * valuey2x1v5 + valuey2x1v6) / (epsilon * epsilon);
    y = yoriginal;
    x = xoriginal;
    double dery2x1 = (-0.5 * der11 + 0.5 * der12) / epsilon;
    std::cout << dery2x1 << std::endl;

    y = yoriginal - 2 * epsilon;
    double valuey3m2 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal - epsilon;
    double valuey3m1 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal + epsilon;
    double valuey3p1 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal + 2 * epsilon;
    double valuey3p2 = std::sin(std::log(y) * std::cos(x));
    y = yoriginal;

    double dery3fd =
        (-0.5 * valuey3m2 + valuey3m1 - valuey3p1 + 0.5 * valuey3p2) /
        (epsilon * epsilon * epsilon);
    std::cout << dery3fd << std::endl;
}

TEST(Raw, Der3Sin) {
    // sin( y * x )
    // y = 0.5 x = 0.3
    double y = 0.5;
    double x = 0.3;
    double yx = y * x;
    double syx = std::sin(yx);

    auto [whatever] = Init<1>();
    auto ders_sin = sin_t<decltype(whatever)>::d2<3>(syx, yx);

    std::array<std::array<std::size_t, 3>, partition_function(3)> parts{};

    constexpr auto last = FirstPartition<3>();
    constexpr auto first = LastPartition<3>();

    parts[0] = LastPartition<3>();
    parts[1] = PrevPartition(first);
    parts[2] = PrevPartition(parts[1]);

    constexpr auto rnext1 = PrevPartition(first);
    constexpr auto rnext2 = PrevPartition(rnext1);
    static_assert(last == rnext2);

    constexpr std::size_t order1 =
        std::accumulate(first.begin(), first.end(), 0UL) - 1;
    constexpr std::size_t order2 =
        std::accumulate(rnext1.begin(), rnext1.end(), 0UL) - 1;
    constexpr std::size_t order3 =
        std::accumulate(rnext2.begin(), rnext2.end(), 0UL) - 1;

    std::vector<double> buffer(2);

    // <y*x, 1, 3>
    // <y*x, 1, 1> <y*x, 2, 1>
    // <y*x, 3, 1> der 0!
    buffer[0] = BellCoeff(first);
    buffer[0] *= ders_sin[order1];
    buffer[1] = BellCoeff(rnext1);
    // buffer[1] = 1.;
    buffer[1] *= ders_sin[order2];
    // buffer[2] = ders_sin[order3] * BellCoeff(rnext2.arr);

    static_assert(combinations(2, 1) == 2);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 1)>
        multparts1{{{1, 0}, {0, 1}}};

    static_assert(combinations(2, 2) == 3);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
        multparts2{{{2, 0}, {1, 1}, {0, 2}}};

    static_assert(combinations(2, 3) == 4);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 3)>
        multparts3{{{3, 0}, {2, 1}, {1, 2}, {0, 3}}};

    constexpr auto factorials = factorialarr<3>();
    // https://calculus.subwiki.org/wiki/Product_rule_for_higher_derivatives
    // <y * x, 1, 1> ->
    // x * <y, 1, 1>
    // y * <x, 1, 1>
    std::array<double, 2> lncs1{};
    lncs1[0] = MultinomialCoeff<1, 2>(factorials, multparts1[0]);
    lncs1[1] = MultinomialCoeff<1, 2>(factorials, multparts1[1]);
    lncs1.front() *= x;
    lncs1.back() *= y;
    // <y * x, 2, 1> ->
    // x * <y, 2, 1>,
    // <x, 1, 1> * <y, 1, 1>
    // y * <x, 2, 1>
    std::array<double, 3> lncs2{};
    lncs2[0] = MultinomialCoeff<2, 2>(factorials, multparts2[0]);
    lncs2[1] = MultinomialCoeff<2, 2>(factorials, multparts2[1]);
    lncs2[2] = MultinomialCoeff<2, 2>(factorials, multparts2[2]);
    lncs2.front() *= x;
    lncs2.back() *= y;

    // <y * x, 3, 1> ->
    // x * <y, 3, 1>
    // <y, 2, 1> <x, 1, 1>
    // <y, 1, 1> <x, 2, 1>
    // y * <x, 3, 1>
    // std::array<double, 4> lncs3{};
    // lncs3[0] = MultinomialCoeff<3, 2>(factorials, multparts3[0]);
    // lncs3[1] = MultinomialCoeff<3, 2>(factorials, multparts3[1]);
    // lncs3[2] = MultinomialCoeff<3, 2>(factorials, multparts3[2]);
    // lncs3[3] = MultinomialCoeff<3, 2>(factorials, multparts3[3]);
    // lncs3.front() *= x;
    // lncs3.back() *= y;

    std::array<double, 4> cxarr{};
    cxarr[0] = 1.;
    cxarr[1] = lncs1[0];
    cxarr[2] = cxarr[1] * lncs1[0];
    cxarr[3] = cxarr[2] * lncs1[0];

    std::array<double, 4> lyarr{};
    lyarr[0] = 1.;
    lyarr[1] = lncs1[1];
    lyarr[2] = lyarr[1] * lncs1[1];
    lyarr[3] = lyarr[2] * lncs1[1];

    // <a * <y, 1, 1> +
    // b * <x, 1, 1>, 3> ->
    std::array<double, combinations(2, 3)> lncs13{};
    lncs13[0] = MultinomialCoeff<3, 2>(factorials, multparts3[0]);
    lncs13[0] *= cxarr[multparts3[0][0]] * lyarr[multparts3[0][1]];

    lncs13[1] = MultinomialCoeff<3, 2>(factorials, multparts3[1]);
    // lncs13[1] = 1;
    lncs13[1] *= cxarr[multparts3[1][0]] * lyarr[multparts3[1][1]];

    lncs13[2] = MultinomialCoeff<3, 2>(factorials, multparts3[2]);
    // lncs13[2] = 1;
    lncs13[2] *= cxarr[multparts3[2][0]] * lyarr[multparts3[2][1]];

    lncs13[3] = MultinomialCoeff<3, 2>(factorials, multparts3[3]);
    lncs13[3] = cxarr[multparts3[3][0]] * lyarr[multparts3[3][1]];

    // <y*x, 1, 3>
    // <y*x, 1, 1> <y*x, 2, 1>
    // <y*x, 3, 1>

    // x3, x2y1, x1y2, y3
    std::array<double, 4> results{};
    results.fill(0);

    // 0: empty
    // 1: <y*x, 1, 1> <y*x, 2, 1>
    // 2: <y*x, 3, 1>
    results[0] += buffer[0] * lncs13[3];
    results[3] += buffer[0] * lncs13[0];
    results[1] += buffer[0] * lncs13[2];
    results[2] += buffer[0] * lncs13[1];
    // buffer[0] is free

    // 0: <y, 1, 1> <y*x, 2, 1>
    // 1: <x, 1, 1> <y*x, 2, 1>
    // 2: <y*x, 3, 1>
    buffer[0] = lncs1[0] * buffer[1];
    buffer[1] *= lncs1[1];

    // 2: <y*x, 3, 1>
    results[2] += buffer[0] * lncs2[1];
    results[1] += buffer[1] * lncs2[1];

    constexpr double third = 1. / 3.;
    // constexpr double third = 1.;
    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << results[0] << std::endl;
    // std::cout << results[1] * third << std::endl;
    // std::cout << results[2] * third << std::endl;
    // std::cout << results[3] << std::endl;
    EXPECT_NEAR(results[0], -0.12359638474200528, 1e-14);
    EXPECT_NEAR(results[1] * third, -0.22359596331880238, 1e-14);
    EXPECT_NEAR(results[2] * third, -0.13415757799128142, 1e-14);
    EXPECT_NEAR(results[3], -0.026696819104273139, 1e-14);
}

TEST(Raw, Der4Sin) {
    // sin( y * x )
    // y = 0.5 x = 0.3
    double y = 0.5;
    double x = 0.3;
    double yx = y * x;
    double syx = std::sin(yx);

    auto [whatever] = Init<1>();
    auto ders_sin = sin_t<decltype(whatever)>::d2<4>(syx, yx);

    std::array<std::array<std::size_t, 4>, partition_function(4)> parts{};

    constexpr auto last = FirstPartition<4>();
    constexpr auto first = LastPartition<4>();

    constexpr std::array<std::array<std::size_t, 4>, 5> parts2{
        {{4, 0, 0, 0}, {2, 1, 0, 0}, {0, 2, 0, 0}, {1, 0, 1, 0}, {0, 0, 0, 1}}};

    parts[0] = LastPartition<4>();
    parts[1] = PrevPartition(first);
    parts[2] = PrevPartition(parts[1]);

    constexpr auto rnext1 = PrevPartition(first);
    constexpr auto rnext2 = PrevPartition(rnext1);
    constexpr auto rnext3 = PrevPartition(rnext2);
    constexpr auto rnext4 = PrevPartition(rnext3);
    static_assert(last == rnext4);

    constexpr std::size_t order1 =
        std::accumulate(first.begin(), first.end(), 0UL) - 1;
    constexpr std::size_t order2 =
        std::accumulate(rnext1.begin(), rnext1.end(), 0UL) - 1;
    constexpr std::size_t order3 =
        std::accumulate(rnext2.begin(), rnext2.end(), 0UL) - 1;
    constexpr std::size_t order4 =
        std::accumulate(rnext3.begin(), rnext3.end(), 0UL) - 1;
    constexpr std::size_t order5 =
        std::accumulate(rnext4.begin(), rnext4.end(), 0UL) - 1;

    std::vector<double> buffer(10);

    // 0: <y*x, 1, 4>
    // 1: <y*x, 1, 2> <y*x, 2, 1>
    // 2: <y*x, 2, 2>
    // <y*x, 1, 1> <y*x, 3, 1> = 0
    // <y*x, 4, 1> = 0

    buffer[0] = BellCoeff(first);
    buffer[0] *= ders_sin[order1];
    buffer[1] = BellCoeff(rnext1);
    buffer[1] *= ders_sin[order2];
    buffer[2] = BellCoeff(rnext2);
    buffer[2] *= ders_sin[order3];

    static_assert(combinations(2, 1) == 2);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 1)>
        multparts1{{{1, 0}, {0, 1}}};

    static_assert(combinations(2, 2) == 3);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
        multparts2{{{2, 0}, {1, 1}, {0, 2}}};
    constexpr auto factorials = factorialarr<4>();
    // https://calculus.subwiki.org/wiki/Product_rule_for_higher_derivatives
    // <y * x, 1, 1> ->
    // x * <y, 1, 1>
    // y * <x, 1, 1>
    std::array<double, 2> lncs1{};
    lncs1[0] = MultinomialCoeff<1, 2>(factorials, multparts1[0]);
    lncs1[1] = MultinomialCoeff<1, 2>(factorials, multparts1[1]);
    lncs1.front() *= x;
    lncs1.back() *= y;
    // <y * x, 2, 1> ->
    // x * <y, 2, 1>,
    // <x, 1, 1> * <y, 1, 1>
    // y * <x, 2, 1>
    std::array<double, 3> lncs2{};
    lncs2[0] = MultinomialCoeff<2, 2>(factorials, multparts2[0]);
    lncs2[1] = MultinomialCoeff<2, 2>(factorials, multparts2[1]);
    lncs2[2] = MultinomialCoeff<2, 2>(factorials, multparts2[2]);
    lncs2.front() *= x;
    lncs2.back() *= y;

    std::array<double, 5> cxarr{};
    cxarr[0] = 1.;
    cxarr[1] = lncs1[0];
    cxarr[2] = cxarr[1] * lncs1[0];
    cxarr[3] = cxarr[2] * lncs1[0];
    cxarr[4] = cxarr[3] * lncs1[0];

    std::array<double, 5> lyarr{};
    lyarr[0] = 1.;
    lyarr[1] = lncs1[1];
    lyarr[2] = lyarr[1] * lncs1[1];
    lyarr[3] = lyarr[2] * lncs1[1];
    lyarr[4] = lyarr[3] * lncs1[1];

    static_assert(combinations(2, 4) == 5);
    constexpr std::array<std::array<std::size_t, 2>, combinations(2, 4)>
        multparts4{{{4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4}}};
    // <a * <y, 1, 1> +
    // b * <x, 1, 1>, 4> ->
    std::array<double, combinations(2, 4)> lncs14{};
    lncs14[0] = MultinomialCoeff<4, 2>(factorials, multparts4[0]);
    lncs14[0] *= cxarr[multparts4[0][0]] * lyarr[multparts4[0][1]];

    lncs14[1] = MultinomialCoeff<4, 2>(factorials, multparts4[1]);
    // lncs13[1] = 1;
    lncs14[1] *= cxarr[multparts4[1][0]] * lyarr[multparts4[1][1]];

    lncs14[2] = MultinomialCoeff<4, 2>(factorials, multparts4[2]);
    // lncs13[2] = 1;
    lncs14[2] *= cxarr[multparts4[2][0]] * lyarr[multparts4[2][1]];

    lncs14[3] = MultinomialCoeff<4, 2>(factorials, multparts4[3]);
    lncs14[3] *= cxarr[multparts4[3][0]] * lyarr[multparts4[3][1]];

    lncs14[4] = MultinomialCoeff<4, 2>(factorials, multparts4[4]);
    lncs14[4] *= cxarr[multparts4[4][0]] * lyarr[multparts4[4][1]];

    // y4, y3x1, y2x2, y1x3, x4
    std::array<double, 5> results{};
    results.fill(0);

    // 0: <y*x, 1, 4>
    // 1: <y*x, 1, 2> <y*x, 2, 1>
    // 2: <y*x, 2, 2>

    results[4] += buffer[0] * lncs14[4];
    results[0] += buffer[0] * lncs14[0];
    results[1] += buffer[0] * lncs14[1];
    results[2] += buffer[0] * lncs14[2];
    results[3] += buffer[0] * lncs14[3];

    // 0: free
    // 1: <y*x, 1, 2> <y*x, 2, 1>
    // 2: <y*x, 2, 2>

    // static_assert(combinations(2, 2) == 3);
    // constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
    //     multparts2{{{2, 0}, {1, 1}, {0, 2}}};

    std::array<double, combinations(2, 2)> lncs12{};
    lncs12[0] = MultinomialCoeff<2, 2>(factorials, multparts2[0]);
    lncs12[0] *= cxarr[multparts2[0][0]] * lyarr[multparts2[0][1]];

    lncs12[1] = MultinomialCoeff<2, 2>(factorials, multparts2[1]);
    lncs12[1] *= cxarr[multparts2[1][0]] * lyarr[multparts2[1][1]];

    lncs12[2] = MultinomialCoeff<2, 2>(factorials, multparts2[2]);
    lncs12[2] *= cxarr[multparts2[2][0]] * lyarr[multparts2[2][1]];

    // 0: <y, 1, 2> <y*x, 2, 1>
    // 1: <x, 1, 2> <y*x, 2, 1>
    // 2: <y*x, 2, 2>
    // 3: <y, 1, 1> <x, 1, 1> <y*x, 2, 1>
    buffer[0] = lncs12[0] * buffer[1];
    buffer[3] = lncs12[1] * buffer[1];
    buffer[1] *= lncs12[2];

    // 0: empty
    // 1: <x, 1, 2> <y*x, 2, 1>
    // 2: <y*x, 2, 2>
    // 3: <y, 1, 1> <x, 1, 1> <y*x, 2, 1>
    // results[0] += buffer[0] * lncs2[0];
    results[1] += buffer[0] * lncs2[1];
    // results[2] += buffer[0] * lncs2[2];

    // 0: empty
    // 1: empty
    // 2: <y*x, 2, 2>
    // 3: <y, 1, 1> <x, 1, 1> <y*x, 2, 1>
    // results[2] += buffer[1] * lncs2[0];
    results[3] += buffer[1] * lncs2[1];
    // results[4] += buffer[1] * lncs2[2];

    // 0: empty
    // 1: empty
    // 2: <y*x, 2, 2>
    // 3: empty
    // results[1] += buffer[3] * lncs2[0];
    results[2] += buffer[3] * lncs2[1];
    // results[3] += buffer[3] * lncs2[2];

    results[2] += buffer[2] * lncs2[1] * lncs2[1];

    EXPECT_NEAR(results[0], 0.00121044887303615, 1e-14);
    EXPECT_NEAR(results[1] / 4, -0.26495077625433783, 1e-14);
    EXPECT_NEAR(results[2] / 6, -0.88877655372816777, 1e-14);
    EXPECT_NEAR(results[3] / 4, -0.73597437848427172, 1e-14);
    EXPECT_NEAR(results[4], 0.009339883279599951, 1e-14);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << results[0] << std::endl;
    // std::cout << results[1] / 4 << std::endl;
    // std::cout << results[2] / 6 << std::endl;
    // std::cout << results[3] / 4 << std::endl;
    // std::cout << results[4] << std::endl;
}

} // namespace adhoc3
