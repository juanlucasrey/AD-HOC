#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

template<typename T = double>
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
