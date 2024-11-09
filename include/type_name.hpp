#ifndef TESTS4_TYPE_NAME_HPP
#define TESTS4_TYPE_NAME_HPP

#include <string>

template <class T> constexpr std::string_view type_name() {
    using namespace std;
#ifdef __clang__
    string_view p = __PRETTY_FUNCTION__;
    return string_view(p.data() + 34, p.size() - 34 - 1);
#elif defined(__GNUC__)
    string_view p = __PRETTY_FUNCTION__;
#if __cplusplus < 201402
    return string_view(p.data() + 36, p.size() - 36 - 1);
#else
    return string_view(p.data() + 49, p.find(';', 49) - 49);
#endif
#elif defined(_MSC_VER)
    string_view p = __FUNCSIG__;
    return string_view(p.data() + 84, p.size() - 84 - 7);
#endif
}

template <class T> constexpr std::string type_name2() {

    auto result = type_name<T>();

    std::string_view open{"<"};
    std::string_view close{">"};
    std::string_view comma{", "};
    std::string_view tab{"    "};

    std::size_t current_pos = 0;
    std::size_t next_pos = 0;
    std::size_t tabs = 0;

    std::string final_result;

    std::size_t find_open = result.find(open, current_pos + 1);
    std::size_t find_close = result.find(close, current_pos + 1);
    std::size_t find_comma = result.find(comma, current_pos + 1);

    while (find_open != std::string::npos || find_close != std::string::npos ||
           find_comma != std::string::npos) {
        for (std::size_t i = 0; i < tabs; i++) {
            final_result += tab;
        }

        if (find_close < find_open && find_close < find_comma) {
            next_pos = find_close;

            final_result += result.substr(current_pos, next_pos - current_pos);
            final_result += "\n";
            current_pos = next_pos;

            if (find_close < find_open) {
                if (tabs == 0) {
                    throw;
                }
                tabs--;
            } else {
                tabs++;
            }

            find_open = result.find(open, current_pos + 1);
            find_close = result.find(close, current_pos + 1);
            find_comma = result.find(comma, current_pos + 1);
        } else if (find_open < find_close && find_open < find_comma) {
            next_pos = find_open + open.size();

            final_result += result.substr(current_pos, next_pos - current_pos);
            final_result += "\n";
            current_pos = next_pos;

            if (find_close < find_open) {
                if (tabs == 0) {
                    throw;
                }
                tabs--;
            } else {
                tabs++;
            }

            find_open = result.find(open, current_pos);
            find_close = result.find(close, current_pos);
            find_comma = result.find(comma, current_pos);
        } else {
            next_pos = find_comma + comma.size();

            final_result += result.substr(current_pos, next_pos - current_pos);
            final_result += "\n";
            current_pos = next_pos;

            find_open = result.find(open, current_pos);
            find_close = result.find(close, current_pos);
            find_comma = result.find(comma, current_pos);
        }
    }
    final_result += result.substr(current_pos, result.size() - current_pos);

    return final_result;
}

#endif // TESTS4_TYPE_NAME_HPP
