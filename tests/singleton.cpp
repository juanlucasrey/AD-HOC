#include <type_traits>
#include <vector>

#include <adhoc2.hpp>
#include <tape_input.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

class MySingleton {
  private:
    MySingleton() = default;

  public:
    MySingleton(MySingleton const &) = delete;
    auto operator=(MySingleton const &) -> MySingleton & = delete;

    MySingleton(MySingleton &&) = delete;
    auto operator=(MySingleton &&) -> MySingleton & = delete;

    // std::vector<double> buf;
    double buf;
    static auto getInstance() -> MySingleton & {
        // The only instance
        // Guaranteed to be lazy initialized
        // Guaranteed that it will be destroyed correctly
        static MySingleton instance;
        return instance;
    }

    // ~MySingleton() {
    //     // auto &temp = MySingleton::getInstance();
    //     delete MySingleton::getInstance();
    // };
};

TEST(adhoc2, singleton) {
    auto &temp = MySingleton::getInstance();
    temp.buf = 2;
    // temp.buf.size();
    std::cout << sizeof(temp) << std::endl;
    std::cout << temp.buf << std::endl;
}

} // namespace adhoc2
