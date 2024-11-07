#include <adhoc.hpp>
#include <dependency.hpp>

int main() {
    ADHOC(vol);
    ADHOC(strike);
    ADHOC(spot);

    auto res = vol * strike + exp(spot);
    static_assert(is_input(vol));
    static_assert(is_input(strike));
    static_assert(is_input(spot));
    static_assert(!is_input(res));

    return 0;
}
