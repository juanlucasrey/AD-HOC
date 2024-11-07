#include <adhoc.hpp>
#include <dependency.hpp>

int main() {
    // is_input check
    {
        ADHOC(vol);
        ADHOC(strike);
        ADHOC(spot);

        auto res = vol * strike + exp(spot);
        static_assert(is_input(vol));
        static_assert(is_input(strike));
        static_assert(is_input(spot));
        static_assert(!is_input(res));
    }

    // is_constant check
    {
        ADHOC(val0);

        using adhoc4::constants::CD;
        auto constant = CD<0.5>();
        auto res = val0 * constant;

        static_assert(is_constant(constant));
        static_assert(!is_constant(res));
    }

    // oeprator *, + and == check
    {
        ADHOC(var1);
        ADHOC(var2);
        ADHOC(var3);
        ADHOC(var4);

        auto res1 = (var1 + var2) * (var3 * var4);
        auto res2 = (var4 * var3) * (var2 + var1);
        static_assert(res1 == res2);

        auto res3 = (var1 * var2) + (var3 + var4);
        auto res4 = (var4 + var3) + (var2 * var1);
        static_assert(res3 == res4);
    }

    return 0;
}
