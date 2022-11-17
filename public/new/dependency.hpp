#ifndef ADHOC_DEPENDENCY_HPP
#define ADHOC_DEPENDENCY_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"

#include <type_traits>

namespace adhoc2 {

template <class First, class Second> struct equal_or_depends;

#if __cplusplus >= 202002L
template <constants::detail::AsTemplateArg<double> F, class Second>
struct equal_or_depends<constants::CD<F>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::CD<F>, Second> ||
               std::is_same_v<constants::CD<F> const, Second>;
    }
};

#else

template <std::uint64_t N, class Second>
struct equal_or_depends<constants::CD<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::CD<N>, Second>;
    }
};

#endif

template <std::size_t N, class Second>
struct equal_or_depends<double_t<N>, Second> {
    constexpr static auto call() noexcept -> bool;
};

template <std::size_t N, class Second>
constexpr auto equal_or_depends<double_t<N>, Second>::call() noexcept -> bool {
    return std::is_same_v<double_t<N>, Second>;
}

template <template <class> class Univariate, class Input, class Second>
struct equal_or_depends<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> bool;
};

template <template <class> class Univariate, class Input, class Second>
constexpr auto equal_or_depends<Univariate<Input>, Second>::call() noexcept
    -> bool {
    static_assert(
        std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
    return std::is_same_v<Univariate<Input>, Second> ||
           equal_or_depends<Input, Second>::call();
}

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct equal_or_depends<Bivariate<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> bool;
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
constexpr auto
equal_or_depends<Bivariate<Input1, Input2>, Second>::call() noexcept -> bool {
    static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                        Base<Bivariate<Input1, Input2>>>);
    return std::is_same_v<Bivariate<Input1, Input2>, Second> ||
           equal_or_depends<Input1, Second>::call() ||
           equal_or_depends<Input2, Second>::call();
}

template <class First, class Second> struct depends;

template <constants::ArgType N, class Second>
struct depends<constants::CD<N>, Second> {
    constexpr static auto call() noexcept -> bool { return false; }
};

template <std::size_t N, class Second> struct depends<double_t<N>, Second> {
    constexpr static auto call() noexcept -> bool { return false; }
};

template <template <class> class Univariate, class Input, class Second>
struct depends<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(
            std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
        return equal_or_depends<Input, Second>::call();
    }
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct depends<Bivariate<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                            Base<Bivariate<Input1, Input2>>>);
        return equal_or_depends<Input1, Second>::call() ||
               equal_or_depends<Input2, Second>::call();
    }
};

template <class First, class... Others>
inline constexpr auto equal_or_depends_many() noexcept -> bool {
    if constexpr (sizeof...(Others) == 0) {
        return false;
    } else {
        return (equal_or_depends<First, Others>::call() || ...);
    }
}

} // namespace adhoc2

#endif // ADHOC_DEPENDENCY_HPP