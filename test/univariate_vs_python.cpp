
#include <adhoc.hpp>
#include <base.hpp>

#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

#include "statistics.hpp"
#include "test_tools.hpp"

int main() {
    TEST_START;

    using namespace adhoc;

    std::string name = "data/lgamma.csv";
    std::vector<double> in;
    std::vector<std::array<double, 7>> out;

    get_values(name, in, out);

    std::function<double(double)> functionstd = [](double d) {
        return std::lgamma(d);
    };

    auto stats = statistics(in, out, functionstd);
    EXPECT_LESS_THAN(std::get<0>(stats), 1e-16);
    EXPECT_LESS_THAN(std::get<1>(stats), 1e-16);
    EXPECT_LESS_THAN(std::get<2>(stats), 1e-15);

    std::vector<std::array<double, 3>> errors{{1e-15, 1e-14, 1e-12},
                                              {1e-16, 1e-15, 1e-15},
                                              {1e-15, 1e-15, 1e-13},
                                              {1e-16, 1e-15, 1e-15}};

    for (std::size_t i = 1; i < 5; i++) {

        std::function<double(double)> function = [i](double d) {
            std::array<double, 10> results1;
            results1.fill(0.);
            lgamma_t<double>::d<6>(d, d, results1);
            return results1[i - 1];
        };

        auto stats = statistics(in, out, function, i);
        EXPECT_LESS_THAN(std::get<0>(stats), errors[i - 1][0]);
        EXPECT_LESS_THAN(std::get<1>(stats), errors[i - 1][1]);
        EXPECT_LESS_THAN(std::get<2>(stats), errors[i - 1][2]);
    }

    TEST_END;
}
