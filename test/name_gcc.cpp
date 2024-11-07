#include <adhoc.hpp>
#include <base.hpp>
#include <name.hpp>

#include <string_view>

int main() {

    {
        ADHOC(var);
        ADHOC(var2);
        constexpr auto res_int = d(var) * d(var2 * var);
        constexpr auto res = name(res_int * res_int * d(cos(var)));
        constexpr std::string_view res_check{
            "d^2(var)*d^2(var2*var)*d(cos(var))"};
        static_assert(res_check == std::string_view(res.data(), res.size()));
    }

    return 0;
}
