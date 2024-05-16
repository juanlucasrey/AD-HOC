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

#ifndef ADHOC3_STRICT_ORDER_HPP
#define ADHOC3_STRICT_ORDER_HPP

#include <string>

// got this from
// https://stackoverflow.com/questions/64794649/is-there-a-way-to-consistently-sort-type-template-parameters

namespace adhoc3::detail {

template <class T> constexpr std::string_view type_name() {
    using namespace std;
#ifdef __clang__
    string_view p = __PRETTY_FUNCTION__;
    return string_view(p.data() + 34, p.size() - 34 - 1);
#elif defined(__GNUC__)
    string_view p = __PRETTY_FUNCTION__;
#if __cplusplus < 201402
    return string_view(p.data() + 36, p.size() - 36 - 1);
#else
    // return string_view(p.data() + 57, p.find(';', 57) - 57);
    return string_view(p.data() + 65, p.find(';', 65) - 65);
#endif
#elif defined(_MSC_VER)
    string_view p = __FUNCSIG__;
    return string_view(p.data() + 84, p.size() - 84 - 7);
#endif
}

template <class... Ts> struct list;

template <template <class...> class Ins, class...> struct instantiate;

template <template <class...> class Ins, class... Ts>
struct instantiate<Ins, list<Ts...>> {
    using type = Ins<Ts...>;
};

template <template <class...> class Ins, class... Ts>
using instantiate_t = typename instantiate<Ins, Ts...>::type;

template <class...> struct concat;

template <class... Ts, class... Us> struct concat<list<Ts...>, list<Us...>> {
    using type = list<Ts..., Us...>;
};

template <class... Ts> using concat_t = typename concat<Ts...>::type;

template <int Count, class... Ts> struct take;

template <int Count, class... Ts>
using take_t = typename take<Count, Ts...>::type;

template <class... Ts> struct take<0, list<Ts...>> {
    using type = list<>;
    using rest = list<Ts...>;
};

template <class A, class... Ts> struct take<1, list<A, Ts...>> {
    using type = list<A>;
    using rest = list<Ts...>;
};

template <int Count, class A, class... Ts> struct take<Count, list<A, Ts...>> {
    using type = concat_t<list<A>, take_t<Count - 1, list<Ts...>>>;
    using rest = typename take<Count - 1, list<Ts...>>::rest;
};

template <class... Types> struct sort_list;

template <class... Ts> using sorted_list_t = typename sort_list<Ts...>::type;

template <class A> struct sort_list<list<A>> {
    using type = list<A>;
};

template <class Left, class Right>
static constexpr bool less_than = type_name<Left>() < type_name<Right>();

template <class A, class B> struct sort_list<list<A, B>> {
    using type = std::conditional_t<less_than<A, B>, list<A, B>, list<B, A>>;
};

template <class...> struct merge;

template <class... Ts> using merge_t = typename merge<Ts...>::type;

template <class... Bs> struct merge<list<>, list<Bs...>> {
    using type = list<Bs...>;
};

template <class... As> struct merge<list<As...>, list<>> {
    using type = list<As...>;
};

template <class AHead, class... As, class BHead, class... Bs>
struct merge<list<AHead, As...>, list<BHead, Bs...>> {
    using type = std::conditional_t<
        less_than<AHead, BHead>,
        concat_t<list<AHead>, merge_t<list<As...>, list<BHead, Bs...>>>,
        concat_t<list<BHead>, merge_t<list<AHead, As...>, list<Bs...>>>>;
};

template <class... Types> struct sort_list<list<Types...>> {
    static constexpr auto first_size = sizeof...(Types) / 2;
    using split = take<first_size, list<Types...>>;
    using type = merge_t<sorted_list_t<typename split::type>,
                         sorted_list_t<typename split::rest>>;
};

} // namespace adhoc3::detail

#endif // ADHOC3_STRICT_ORDER_HPP
