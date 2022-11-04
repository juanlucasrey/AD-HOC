#ifndef DEPENDENCY_HPP
#define DEPENDENCY_HPP

#include <adhoc2.hpp>

#include <constants_double.hpp>

// to be removed
#include <constants_type.hpp>

#include <type_traits>

namespace adhoc2 {

template <class First, class Second> struct equal_or_depends;

// to be removed from here
template <int N, class Second>
struct equal_or_depends<constants::Const<N> const, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::Const<N>, Second> ||
               std::is_same_v<constants::Const<N> const, Second>;
    }
};

template <int N, class Second>
struct equal_or_depends<constants::Const<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::Const<N>, Second> ||
               std::is_same_v<constants::Const<N> const, Second>;
    }
};
// to be removed to here

template <int N, class Second>
struct equal_or_depends<constants::Double<N> const, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::Double<N>, Second> ||
               std::is_same_v<constants::Double<N> const, Second>;
    }
};

template <int N, class Second>
struct equal_or_depends<constants::Double<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::Double<N>, Second> ||
               std::is_same_v<constants::Double<N> const, Second>;
    }
};

template <int N, class Second>
struct equal_or_depends<detail::adouble_aux<N> const, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<detail::adouble_aux<N>, Second> ||
               std::is_same_v<detail::adouble_aux<N> const, Second>;
    }
};

template <int N, class Second>
struct equal_or_depends<detail::adouble_aux<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<detail::adouble_aux<N>, Second> ||
               std::is_same_v<detail::adouble_aux<N> const, Second>;
    }
};

template <template <class> class Univariate, class Input, class Second>
struct equal_or_depends<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(
            std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
        return std::is_same_v<Univariate<Input>, Second> ||
               std::is_same_v<Univariate<Input> const, Second> ||
               equal_or_depends<Input, Second>::call();
    }
};

template <template <class> class Univariate, class Input, class Second>
struct equal_or_depends<Univariate<Input> const, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(
            std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
        return std::is_same_v<Univariate<Input>, Second> ||
               std::is_same_v<Univariate<Input> const, Second> ||
               equal_or_depends<Input, Second>::call();
    }
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct equal_or_depends<Bivariate<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                            Base<Bivariate<Input1, Input2>>>);
        return std::is_same_v<Bivariate<Input1, Input2>, Second> ||
               std::is_same_v<Bivariate<Input1, Input2> const, Second> ||
               equal_or_depends<Input1, Second>::call() ||
               equal_or_depends<Input2, Second>::call();
    }
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct equal_or_depends<Bivariate<Input1, Input2> const, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                            Base<Bivariate<Input1, Input2>>>);
        return std::is_same_v<Bivariate<Input1, Input2>, Second> ||
               std::is_same_v<Bivariate<Input1, Input2> const, Second> ||
               equal_or_depends<Input1, Second>::call() ||
               equal_or_depends<Input2, Second>::call();
    }
};

template <class First, class Second> struct depends;

template <int N, class Second>
struct depends<detail::adouble_aux<N> const, Second> {
    constexpr static auto call() noexcept -> bool { return false; }
};

template <int N, class Second> struct depends<detail::adouble_aux<N>, Second> {
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

template <template <class> class Univariate, class Input, class Second>
struct depends<Univariate<Input> const, Second> {
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

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct depends<Bivariate<Input1, Input2> const, Second> {
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

#endif // DEPENDENCY_HPP