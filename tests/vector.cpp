#include <gtest/gtest.h>

#include "type_name.hpp"
#include <constants_type.hpp>
#include <init.hpp>

#include <algorithm>
#include <cmath>
#include <numeric> // transform_reduce

namespace adhoc {

template <class T> class vector {};
template <class T1, class T2> struct pair {};

template <class T> class star : public Base<star<T>> {};
template <class Loc, class Cont> struct iterator {
    auto operator*() { return star<iterator<Loc, Cont>>{}; }
    auto operator*() const { return star<iterator<Loc, Cont>>{}; }
};

class begin_t {};
class end_t {};

template <class Cont> auto begin(Cont const & /* c */) {
    return iterator<begin_t, Cont>{};
};

template <class Cont> auto end(Cont const & /* c */) {
    return iterator<end_t, Cont>{};
};

template <class I> class prev_t {};
template <class I> class next_t {};
template <class I> auto prev(I const & /* c */) { return prev_t<I>{}; };
template <class I> auto next(I const & /* c */) { return next_t<I>{}; };

template <class Loc, class Cont>
auto prev(iterator<Loc, Cont> const & /* c */) {
    return iterator<prev_t<Loc>, Cont>{};
};

template <class Loc, class Cont>
auto next(iterator<Loc, Cont> const & /* c */) {
    return iterator<next_t<Loc>, Cont>{};
};

template <class I1, class I2> class min_t {};
template <class I1, class I2> class max_t {};
template <class I1, class I2>
auto min(I1 const & /* i1 */, I2 const & /* i2 */) {
    return min_t<I1, I2>{};
};
template <class I1, class I2>
auto max(I1 const & /* i1 */, I2 const & /* i2 */) {
    return max_t<I1, I2>{};
};

template <class Loc1, class Loc2, class Cont>
auto min(iterator<Loc1, Cont> const & /* i1 */,
         iterator<Loc2, Cont> const & /* i2 */) {
    return iterator<min_t<Loc1, Loc2>, Cont>{};
};

template <class Loc1, class Loc2, class Cont>
auto max(iterator<Loc1, Cont> const & /* i1 */,
         iterator<Loc2, Cont> const & /* i2 */) {
    return iterator<max_t<Loc1, Loc2>, Cont>{};
};

template <std::size_t P, class Cont>
class get_t : public Base<get_t<P, Cont>> {};
template <std::size_t P, class Cont> auto get(Cont const & /* c */) {
    return get_t<P, Cont>{};
};

template <class T, class Proj> class upper_bound_t {};

template <class, template <class, class> class>
struct is_instance : public std::false_type {};

template <class T1, class T2, template <class, class> class U>
struct is_instance<U<T1, T2>, U> : public std::true_type {};

// namespace impl {
// // our functor, the niebloid
// struct __upper_bound {
//     template <class I, class T, class Comp = std::pair<double, double>,
//               class Proj = std::pair<double, double>,
//               typename = std::enable_if_t<is_instance<I, vector>::value>>
//     auto operator()(iterator<begin_t, I> /* first */,
//                     iterator<end_t, I> /* last */, const T & /* value */,
//                     Comp /* comp */ = {}, Proj proj = {}) const {
//         return iterator<upper_bound_t<T, decltype(proj(I{}))>, I>{};
//         // swap(a, b);
//     }
// };
// } // namespace impl
// inline constexpr impl::__upper_bound upper_bound{};

namespace ranges {
template <class I, class T, class Comp = std::pair<double, double>,
          class Proj = std::pair<double, double>>
auto upper_bound(iterator<begin_t, I> /* first */,
                 iterator<end_t, I> /* last */, const T & /* value */,
                 Comp /* comp */ = {}, Proj proj = {}) {
    return iterator<upper_bound_t<T, decltype(proj(I{}))>, I>{};
};
} // namespace ranges

template <class Arg1, class Arg2, class Arg3, class Arg4, class Arg5>
auto ranges_upper_bound(Arg1 a1, Arg2 a2, Arg3 a3, Arg4 a4, Arg5 a5) {
    if constexpr (is_instance<Arg1, iterator>::value) {
        return ranges::upper_bound(a1, a2, a3, a4, a5);
    } else {
        return std::ranges::upper_bound(a1, a2, a3, a4, a5);
    }
};

template <class A, class B, class T>
class lerp_t : public Base<lerp_t<A, B, T>> {};

template <class A, class B, class T>
auto lerp(A /* a */, B /* b */, T /* t */) {
    return lerp_t<A, B, T>();
}

template <class Cont, class Init, class BinaryReductionOp,
          class UnaryTransformOp>
class transform_reduce_t
    : public Base<
          transform_reduce_t<Cont, Init, BinaryReductionOp, UnaryTransformOp>> {
};

template <template <class> class Cont, typename ContV, class Init,
          class BinaryReductionOp, class UnaryTransformOp>
auto transform_reduce(iterator<begin_t, Cont<ContV>> /* first */,
                      iterator<end_t, Cont<ContV>> /* last */, Init init,
                      BinaryReductionOp reduce, UnaryTransformOp transform) {
    return transform_reduce_t<Cont<ContV>, Init,
                              decltype(reduce(init, transform(ContV{}))),
                              decltype(transform(ContV{}))>();
}

TEST(VECTOR, Deref) {
    vector<pair<double_t<1>, double_t<2>>> payments;
    vector<pair<double_t<0>, double_t<1>>> values;
    auto LinearInterp = [&values]<typename T>(const T in) {
        auto it_middle = ranges_upper_bound(
            begin(values), end(values), in, std::ranges::less{},
            [](const auto &t) { return get<0>(t); });

        auto it_right =
            max(min(it_middle, prev(end(values))), next(begin(values)));
        auto it_left = prev(it_right);

        auto pos =
            (in - get<0>(*it_left)) / (get<0>(*it_right) - get<0>(*it_left));
        return lerp(get<1>(*it_left), get<1>(*it_right), pos);
    };

    auto reducfunc = []<typename T1, typename T2>(
                         T1 const &lhs, T2 const &rhs) { return lhs + rhs; };

    auto lambdafunc = [&]<typename T>(T const &lhs) {
        return get<1>(lhs) * LinearInterp(get<0>(lhs));
    };

    using constants::CD;
    using constants::encode;
    CD<encode(0)> init;
    double initd = init;
    auto res = transform_reduce(begin(payments), end(payments), init, reducfunc,
                                lambdafunc);
    std::cout << type_name<decltype(res)>() << std::endl;
    std::cout << sizeof(res) << std::endl;
}

template <class PaymentVec, class ValuesVec, class Init>
auto ycfunc(PaymentVec payments, ValuesVec values, Init init) {
    using std::begin;
    using std::end;
    using std::get;
    using std::lerp;
    using std::max;
    using std::min;
    using std::next;
    using std::prev;
    using std::transform_reduce;
    auto LinearInterp = [&values]<typename T>(const T in) {
        auto it_middle = ranges_upper_bound(
            begin(values), end(values), in, std::ranges::less{},
            [](const auto &t) { return get<0>(t); });

        auto it_right =
            max(min(it_middle, prev(end(values))), next(begin(values)));
        auto it_left = prev(it_right);

        auto pos =
            (in - get<0>(*it_left)) / (get<0>(*it_right) - get<0>(*it_left));
        return lerp(get<1>(*it_left), get<1>(*it_right), pos);
    };

    auto reducfunc = []<typename T1, typename T2>(
                         T1 const &lhs, T2 const &rhs) { return lhs + rhs; };

    auto lambdafunc = [&]<typename T>(T const &lhs) {
        return get<1>(lhs) * LinearInterp(get<0>(lhs));
    };

    return transform_reduce(begin(payments), end(payments), init, reducfunc,
                            lambdafunc);
}

TEST(VECTOR, testdouble) {
    std::vector<std::pair<double, double>> payments = {
        {1., 0.02}, {2., 0.02}, {3., 0.01}, {4., 0.01}, {5., 0.03}, {6., 1.03}};

    std::vector<std::pair<double, double>> values;
    values.emplace_back(0, 1.0);
    values.emplace_back(1, 0.98);
    values.emplace_back(5, 0.97);
    values.emplace_back(10, 0.94);

    auto res = ycfunc(payments, values, 0.);
    EXPECT_NEAR(res, 1.080645, 1e-14);
    std::cout << type_name<decltype(res)>() << std::endl;
    std::cout << sizeof(res) << std::endl;

    vector<pair<double_t<1>, double_t<2>>> paymentsT;
    vector<pair<double_t<0>, double_t<1>>> valuesT;
    using constants::CD;
    using constants::encode;
    CD<encode(0)> initT;
    auto resT = ycfunc(paymentsT, valuesT, initT);
    std::cout << type_name<decltype(resT)>() << std::endl;
    std::cout << sizeof(resT) << std::endl;
}

} // namespace adhoc
