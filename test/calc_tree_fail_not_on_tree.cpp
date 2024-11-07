#include <adhoc.hpp>
#include <calc_tree.hpp>

int main() {
    using namespace adhoc4;
    ADHOC(val0);
    ADHOC(val1);
    auto temp = exp(val0);
    CalcTree t(temp);

    t.set(val0) = 1.0;
    t.set(val1) = 1.0; // fails because it's not on the calculation tree
}
