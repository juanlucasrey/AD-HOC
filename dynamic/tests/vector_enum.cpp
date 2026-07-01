#include <adhoc/vector_enum.hpp>
#include <test_simple_include.hpp>

#include <chrono>
#include <random>

namespace {
void
simple_tests()
{
    // generator of integer values between 0 and 7
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 7);

    std::size_t iters = 10000;

    std::vector<std::uint64_t> v0;
    adhoc::vector_enum<std::uint64_t, 8> v1;

    for (std::size_t i = 0; i < iters; ++i) {
        auto const value = dis(gen);
        v0.push_back(value);
        v1.push_back(value);
    }

    for (std::size_t i = 0; i < iters; ++i) {
        EXPECT_EQUAL(v0[i], v1[i]);
    }
}

void
enum_tests()
{
    enum class OpCode : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
    };

    std::size_t iters = 10000;

    // generator of integer values between 0 and 6
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 6);

    std::vector<OpCode> v0;
    adhoc::vector_enum<OpCode, 7> v1;
    adhoc::vector_enum2<OpCode, 7> v2;

    for (std::size_t i = 0; i < iters; ++i) {
        auto const value = static_cast<OpCode>(dis(gen));
        v0.push_back(value);
        v1.push_back(value);
        v2.push_back(value);
    }

    for (std::size_t i = 0; i < iters; ++i) {
        auto val0 = v0[i];
        auto val1 = v1[i];
        auto val2 = v2[i];

        EXPECT_EQUAL(val0, val1);
        EXPECT_EQUAL(val0, val2);
    }
}

void
resize_tests()
{
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 14);

    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
    };

    std::size_t iters = 67000;
    std::vector<OpCode8> v0;
    adhoc::vector_enum<OpCode8, 15, std::uint64_t> v1;
    adhoc::vector_enum2<OpCode8, 15, std::uint64_t> v2;

    for (std::size_t i = 0; i < iters; ++i) {
        auto const value = static_cast<OpCode8>(dis(gen));
        v0.push_back(value);
        v1.push_back(value);
        v2.push_back(value);
    }

    std::size_t iters_end = 65000;
    for (std::size_t i = iters; i-- > iters_end;) {
        v0.resize(i);
        v1.resize(i);
        v2.resize(i);
        EXPECT_EQUAL(v0.size(), v1.size());
        EXPECT_EQUAL(v0.size(), v2.size());
        for (std::size_t j = 0; j < i; ++j) {
            auto val0 = v0[j];
            auto val1 = v1[j];
            auto val2 = v2[j];
            EXPECT_EQUAL(val0, val1);
        }
    }
}

void
iterator_tests_even()
{
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 14);

    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
    };

    std::size_t iters = 100000;
    std::vector<OpCode8> v0;
    adhoc::vector_enum<OpCode8, 15, std::uint64_t> v1;
    adhoc::vector_enum2<OpCode8, 15, std::uint64_t> v2;

    for (std::size_t i = 0; i < (iters); ++i) {
        auto const value = static_cast<OpCode8>(dis(gen));
        v0.push_back(value);
        v1.push_back(value);
        v2.push_back(value);
    }

    auto it0 = v0.begin();
    auto it1 = v1.begin();
    auto it2 = v2.begin();

    for (std::size_t i = 0; i < (iters); ++i) {
        auto val0 = *it0;
        auto val1 = *it1;
        auto val2 = *it2;
        EXPECT_EQUAL(val0, val1);
        EXPECT_EQUAL(val0, val2);
        ++it0;
        ++it1;
        ++it2;
    }

    auto it0_end = v0.end();
    auto it1_end = v1.end();
    auto it2_end = v2.end();

    EXPECT_EQUAL(it0, it0_end);
    EXPECT_EQUAL(it1, it1_end);
    EXPECT_EQUAL(it2, it2_end);

    for (std::size_t i = 0; i < (iters); ++i) {
        --it0;
        --it1;
        --it2;
        auto val0 = *it0;
        auto val1 = *it1;
        auto val2 = *it2;
        EXPECT_EQUAL(val0, val1);
        EXPECT_EQUAL(val0, val2);
    }
}

void
iterator_tests_odd()
{
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 6);

    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
    };

    std::size_t iters = 100000;
    std::vector<OpCode8> v0;
    adhoc::vector_enum<OpCode8, 7, std::uint64_t> v1;
    adhoc::vector_enum2<OpCode8, 7, std::uint64_t> v2;

    for (std::size_t i = 0; i < (iters); ++i) {
        auto const value = static_cast<OpCode8>(dis(gen));
        v0.push_back(value);
        v1.push_back(value);
        v2.push_back(value);
    }

    auto it0 = v0.begin();
    auto it1 = v1.begin();
    auto it2 = v2.begin();

    for (std::size_t i = 0; i < (iters); ++i) {
        auto val0 = *it0;
        auto val1 = *it1;
        auto val2 = *it2;
        EXPECT_EQUAL(val0, val1);
        EXPECT_EQUAL(val0, val2);
        ++it0;
        ++it1;
        ++it2;
    }

    auto it0_end = v0.end();
    auto it1_end = v1.end();
    auto it2_end = v2.end();

    EXPECT_EQUAL(it0, it0_end);
    EXPECT_EQUAL(it1, it1_end);
    EXPECT_EQUAL(it2, it2_end);

    for (std::size_t i = 0; i < (iters); ++i) {
        --it0;
        --it1;
        --it2;
        auto val0 = *it0;
        auto val1 = *it1;
        auto val2 = *it2;
        EXPECT_EQUAL(val0, val1);
        EXPECT_EQUAL(val0, val2);
    }
}

void
timing_tests()
{
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 15);

    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
        COS,
    };

    std::vector<OpCode8> v0;
    adhoc::vector_enum<OpCode8, 16, std::uint64_t> v1;
    adhoc::vector_enum2<OpCode8, 16, std::uint64_t> v2;

    std::size_t iters = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        iters = std::atoi(env_p);
    }

    if (!iters) {
        return;
    }

    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v0.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector filling:" << std::endl;
        std::cout << time << std::endl;
    }

    gen = std::mt19937(1234);
    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v1.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 filling:" << std::endl;
        std::cout << time << std::endl;
    }

    gen = std::mt19937(1234);
    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v2.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 filling:" << std::endl;
        std::cout << time << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v0[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v1[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v2[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v2.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v2.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }
}

void
timing_tests2()
{
    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 17);

    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
        COS,
        REG_INPUT,
        REG_OUTPUT,
    };

    std::vector<OpCode8> v0;
    adhoc::vector_enum<OpCode8, 18, std::uint64_t> v1;
    adhoc::vector_enum2<OpCode8, 18, std::uint64_t> v2;

    std::size_t iters = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        iters = std::atoi(env_p);
    }

    if (!iters) {
        return;
    }

    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v0.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector filling:" << std::endl;
        std::cout << time << std::endl;
    }

    gen = std::mt19937(1234);
    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v1.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 filling:" << std::endl;
        std::cout << time << std::endl;
    }

    gen = std::mt19937(1234);
    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            auto const value = static_cast<OpCode8>(dis(gen));
            v2.push_back(value);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 filling:" << std::endl;
        std::cout << time << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v0[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v1[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(v2[i]);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v2.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += static_cast<std::uint64_t>(*it);
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::uint64_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v2.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += static_cast<std::uint64_t>(*it);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 2 it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }
}

template<class Container>
auto
size_of_vec(Container const& vec, bool capacity = false) -> std::size_t
{
    std::size_t size = 0;
    size += sizeof(typename Container::value_type) * (capacity ? vec.capacity() : vec.size());
    return size;
}

void
compare_sizes()
{
    enum class OpCode8 : std::uint8_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
        COS,
        REG_INPUT,
        REG_OUTPUT,
    };

    enum class OpCode64 : std::uint64_t {
        ADD,
        SUB,
        MUL,
        ADD_C,
        SUB_C,
        MUL_C,
        NORM, // a.k.a. SQUARE
        ABS,
        INV,
        SQRT,
        POW_C,
        EXP,
        LOG,
        ERF,
        ERFC,
        COS,
        REG_INPUT,
        REG_OUTPUT,
    };

    std::vector<OpCode8> v0;
    std::vector<OpCode64> v1;
    adhoc::vector_enum<OpCode8, 18, std::uint8_t> v2;
    adhoc::vector_enum<OpCode64, 18, std::uint64_t> v3;
    adhoc::vector_enum2<OpCode8, 18, std::uint8_t> v4;
    adhoc::vector_enum2<OpCode64, 18, std::uint64_t> v5;

    std::size_t iters = 10000;

    std::mt19937 gen(1234);
    std::uniform_int_distribution<std::uint64_t> dis(0, 17);

    for (std::size_t i = 0; i < iters; ++i) {
        auto const rnd = dis(gen);
        auto const value8 = static_cast<OpCode8>(rnd);
        auto const value64 = static_cast<OpCode64>(rnd);
        v0.push_back(value8);
        v1.push_back(value64);
        v2.push_back(value8);
        v3.push_back(value64);
        v4.push_back(value8);
        v5.push_back(value64);
    }

    for (std::size_t i = 0; i < iters; ++i) {
        EXPECT_EQUAL(v0[i], v2[i]);
        EXPECT_EQUAL(v1[i], v3[i]);
        EXPECT_EQUAL(v0[i], v4[i]);
        EXPECT_EQUAL(v1[i], v5[i]);
    }

    std::cout << "Size of v0: " << size_of_vec(v0) << " bytes\n";
    std::cout << "Size of v1: " << size_of_vec(v1) << " bytes\n";
    std::cout << "Size of v2: " << v2.size_of() << " bytes\n";
    std::cout << "Size of v3: " << v3.size_of() << " bytes\n";
    std::cout << "Size of v4: " << v4.size_of() << " bytes\n";
    std::cout << "Size of v5: " << v5.size_of() << " bytes\n";
}

} // namespace

auto
main() -> int
{
    simple_tests();
    enum_tests();
    compare_sizes();
    resize_tests();
    iterator_tests_even();
    iterator_tests_odd();
    timing_tests();
    timing_tests2();

    TEST_END;
}
