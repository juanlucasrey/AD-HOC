#include <adhoc.hpp>
#include <differential_operator.hpp>

int main() {
    ADHOC(var);
    ADHOC(var2);
    auto res = d(var) * d(var2 * var);

    static_assert(std::same_as<decltype(res * res),
                               decltype(d<2>(var) * d<2>(var2 * var))>);

    return 0;
}
