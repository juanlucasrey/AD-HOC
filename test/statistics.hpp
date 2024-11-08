#ifndef TEST_STATISTICS_HPP
#define TEST_STATISTICS_HPP

#include <fstream>
#include <sstream>
#include <vector>

template <class D>
auto get_values(std::string const &name, std::vector<D> &in,
                std::vector<std::array<D, 7>> &out) {
    std::ifstream valuesfile;
    valuesfile.open(name);
    std::string line;
    std::string value;
    while (getline(valuesfile, line)) {
        std::istringstream tokenStream(line);
        std::getline(tokenStream, value, ',');
        in.push_back(std::stod(value));

        std::array<D, 7> results;

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

template <class D>
auto statistics(std::vector<D> const &in,
                std::vector<std::array<D, 7>> const &out, auto &function,
                std::size_t n = 0) {
    D rel_average = 0.;
    D abs_rel_average = 0.;
    D abs_rel_max = 0.;
    std::size_t count = 0;
    for (std::size_t i = 0; i < in.size(); i++) {
        D const reference_out = out[i][n];
        D const reference_in = in[i];
        D const candidate_out = function(reference_in);
        D const error = candidate_out - reference_out;

        if (reference_out != 0) {
            D const rel_error = error / reference_out;
            D const abs_rel_error = std::abs(rel_error);
            rel_average += rel_error;
            abs_rel_average += abs_rel_error;
            abs_rel_max = std::max(abs_rel_max, abs_rel_error);
            count++;
        }
    }

    rel_average /= static_cast<D>(count);
    rel_average = std::abs(rel_average);
    abs_rel_average /= static_cast<D>(count);

    return std::make_tuple(rel_average, abs_rel_average, abs_rel_max);
};

#endif // TEST_STATISTICS_HPP
