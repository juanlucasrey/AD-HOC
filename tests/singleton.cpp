#include <vector>

#include <utils.hpp>

#include <gtest/gtest.h>

class MySingleton {
  private:
    MySingleton() = default;

  public:
    MySingleton(MySingleton const &) = delete;
    auto operator=(MySingleton const &) -> MySingleton & = delete;

    MySingleton(MySingleton &&) = delete;
    auto operator=(MySingleton &&) -> MySingleton & = delete;

    // std::vector<double> buf;
    // double buf;
    static auto getInstance() -> MySingleton & {
        // The only instance
        // Guaranteed to be lazy initialized
        // Guaranteed that it will be destroyed correctly
        static MySingleton instance;
        return instance;
    }
};

TEST(adhoc2, singleton) {
    auto &temp = MySingleton::getInstance();
    // temp.buf.size();
    std::cout << sizeof(temp) << std::endl;
}

template <std::size_t I> class somevar {};

template <std::size_t... Idx> class Tape1 {
    std::array<double, sizeof...(Idx)> buff{};

    template <std::size_t FirstIdx, std::size_t... IdxIn>
    constexpr static auto check_all_different() -> bool {
        if constexpr (sizeof...(IdxIn) == 0) {
            return true;
        } else {
            return ((FirstIdx != IdxIn) && ...) &&
                   check_all_different<IdxIn...>();
        }
    }

  public:
    Tape1() = default;

    explicit Tape1(somevar<Idx>...) {
        static_assert(check_all_different<Idx...>());
    }

    template <std::size_t ThisIdx> void set(somevar<ThisIdx>, double in) {
        static_assert(((ThisIdx == Idx) || ...));
        constexpr auto position_on_buff = adhoc2::Get<ThisIdx, Idx...>::value;
        this->buff[position_on_buff] = in;
    }

    template <std::size_t... IdxPrev, std::size_t... IdxNew>
    friend auto createnew(Tape1<IdxPrev...>, somevar<IdxNew>...)
        -> Tape1<IdxPrev..., IdxNew...>;
};

template <std::size_t... IdxPrev, std::size_t... Idx>
auto createnew(Tape1<IdxPrev...> prevtape, somevar<Idx>...)
    -> Tape1<IdxPrev..., Idx...> {
    auto res = Tape1<IdxPrev..., Idx...>();

    for (std::size_t i = 0; i < prevtape.buff.size(); ++i) {
        res.buff[i] = prevtape.buff[i];
    }
    return res;
}

TEST(adhoc2, tapestatic) {
    somevar<0> v1;
    somevar<1> v2;
    somevar<2> v3;
    auto tape = Tape1(v1, v2);
    // auto tape2 = Tape1(v1, v1); // compilation error
    tape.set(v1, 2.);
    tape.set(v2, 3.);
    // tape.set(v3, 2.); // compilation error

    auto tape2 = createnew(tape, v3);
    tape2.set(v3, 4.);
    std::cout << "done" << std::endl;
}