
#include <adhoc.hpp>
#include <base.hpp>

#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

#include "test_tools.hpp"

namespace {

auto get_values(std::string const &name, std::vector<double> &in,
                std::vector<std::array<double, 7>> &out) {
    std::ifstream valuesfile;
    valuesfile.open(name);
    std::string line;
    std::string value;
    while (getline(valuesfile, line)) {
        std::istringstream tokenStream(line);
        std::getline(tokenStream, value, ',');
        in.push_back(std::stod(value));

        std::array<double, 7> results;

        std::getline(tokenStream, value, ',');
        results[0] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[1] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[2] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[3] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[4] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[5] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[6] = std::stod(value);

        out.push_back(results);
    }
    valuesfile.close();
}

auto statistics(std::vector<double> const &in,
                std::vector<std::array<double, 7>> const &out, auto &function,
                std::size_t n = 0) {
    double rel_average = 0.;
    double abs_rel_average = 0.;
    double abs_rel_max = 0.;
    std::size_t count = 0;
    for (std::size_t i = 0; i < in.size(); i++) {
        double const reference_out = out[i][n];
        double const reference_in = in[i];
        double const candidate_out = function(reference_in);
        double const error = candidate_out - reference_out;

        if (reference_out != 0) {
            double const rel_error = error / reference_out;
            double const abs_rel_error = std::abs(rel_error);
            rel_average += rel_error;
            abs_rel_average += abs_rel_error;
            abs_rel_max = std::max(abs_rel_max, abs_rel_error);
            count++;
        }
    }

    rel_average /= static_cast<double>(count);
    rel_average = std::abs(rel_average);
    abs_rel_average /= static_cast<double>(count);

    return std::make_tuple(rel_average, abs_rel_average, abs_rel_max);
};

} // namespace

int main() {
    TEST_START;

    using namespace adhoc4;

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
