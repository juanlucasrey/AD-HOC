#include <adhoc.hpp>
#include <base.hpp>
#include <name.hpp>

#include <string_view>

int main() {
    {
        ADHOC(variablename);
        constexpr auto res = name(variablename);
        constexpr std::string_view res_check{"variablename"};
        static_assert(res_check == std::string_view(res.data(), res.size()));
    }

    {
        ADHOC(var);
        constexpr auto res = name(tgamma(asin(exp(var))));
        constexpr std::string_view res_check{"tgamma(asin(exp(var)))"};
        static_assert(res_check == std::string_view(res.data(), res.size()));
    }

    {
        ADHOC(var);
        ADHOC(var2);
        constexpr auto res = name(d(var) * d(var2 * var));
        constexpr std::string_view res_check{"d(var)*d(var2*var)"};
        static_assert(res_check == std::string_view(res.data(), res.size()));
    }

    return 0;
}
