#ifndef BRNG_READ_CSV
#define BRNG_READ_CSV

#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace adhoc {

template <class UIntType>
auto readCSV(const std::string &filename) -> std::vector<UIntType> {
    std::vector<UIntType> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                UIntType number = std::stoull(value);
                data.push_back(number);
            } catch (const std::exception &e) {
                std::cerr << "Warning: Failed to parse '" << value
                          << "' as uint64_t: " << e.what() << std::endl;
            }
        }
    }

    file.close();
    return data;
}

template <class UIntType = std::uint64_t>
auto split_uint64_to_uint32(const std::vector<UIntType> &input,
                            bool invert = false) -> std::vector<std::uint32_t> {
    std::vector<std::uint32_t> output;
    output.reserve(input.size() * 2); // Reserve space to avoid reallocations

    for (std::uint64_t val : input) {
        std::uint32_t low = static_cast<std::uint32_t>(val & 0xFFFFFFFF);
        std::uint32_t high = static_cast<std::uint32_t>(val >> 32);
        if (invert) {
            output.push_back(high);
            output.push_back(low);
        } else {
            output.push_back(low);
            output.push_back(high);
        }
    }

    return output;
}

} // namespace adhoc

#endif // BRNG_READ_CSV
