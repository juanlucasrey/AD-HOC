#pragma once

#include "../credit/mdspan.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

inline auto
loadCSVToVectorDate(std::string const& filename) -> std::vector<std::chrono::year_month_day>
{
    const auto dataPath = std::filesystem::path(__FILE__).parent_path() / filename;
    std::ifstream inputFile(dataPath);
    if (!inputFile) {
        throw std::runtime_error("Unable to open CSV data file: " + dataPath.string());
    }

    std::string contents;
    std::getline(inputFile, contents);

    std::vector<std::chrono::year_month_day> values;
    std::istringstream csvParser(contents);
    std::string token;
    while (std::getline(csvParser, token, ',')) {
        if (token.find_first_not_of(" \t\n\r") == std::string::npos) {
            continue;
        }

        std::istringstream dateParser(token);
        int year = 0;
        unsigned month = 0;
        unsigned day = 0;
        char yMarker = '\0';
        char slash1 = '\0';
        char slash2 = '\0';

        if (!(dateParser >> year >> yMarker >> slash1 >> month >> slash2 >> day) || yMarker != 'y' || slash1 != '/' ||
            slash2 != '/') {
            throw std::runtime_error("Invalid date token in CSV: '" + token + "' in " + dataPath.string());
        }

        values.emplace_back(std::chrono::year(year), std::chrono::month(month), std::chrono::day(day));
    }

    if (values.empty()) {
        throw std::runtime_error("No CSV values found in: " + dataPath.string());
    }

    return values;
}

template<class T = double>
inline auto
loadCSVToVector(std::string const& filename) -> std::vector<T>
{
    static_assert(std::is_arithmetic_v<T>, "loadCSVToVector requires an arithmetic type");

    const auto dataPath = std::filesystem::path(__FILE__).parent_path() / filename;
    std::ifstream inputFile(dataPath);
    if (!inputFile) {
        throw std::runtime_error("Unable to open CSV data file: " + dataPath.string());
    }

    std::string contents;
    std::getline(inputFile, contents);
    for (char& character : contents) {
        if (character == ',') {
            character = ' ';
        }
    }

    std::istringstream parser(contents);
    std::vector<T> values;
    T value = T{};
    while (parser >> value) {
        values.push_back(value);
    }

    if (values.empty()) {
        throw std::runtime_error("No CSV values found in: " + dataPath.string());
    }

    return values;
}

template<class T = double, class... DimTs>
inline auto
loadCSVToTensor(std::string const& filename, DimTs... dims)
  -> std::tuple<std::vector<T>, adhoc::mdspan<const T, sizeof...(DimTs)> >
{
    static_assert(sizeof...(DimTs) > 0, "loadCSVToTensor requires at least one dimension");
    static_assert((std::is_integral_v<DimTs> && ...), "loadCSVToTensor dimensions must be integral");

    auto vals = loadCSVToVector<T>(filename);

    std::size_t expectedSize = (std::size_t{ 1 } * ... * static_cast<std::size_t>(dims));

    if (vals.size() != expectedSize) {
        throw std::runtime_error("CSV data size does not match expected dimensions: " + filename);
    }

    auto mdspanView = adhoc::mdspan<const T, sizeof...(DimTs)>(vals.data(), static_cast<std::size_t>(dims)...);

    return std::make_tuple(std::move(vals), mdspanView);
}
