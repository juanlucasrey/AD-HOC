#include <adhoc/vector_idx.hpp>
#include <test_simple_include.hpp>

#include <chrono>

namespace {
void
simple_tests()
{
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;
    v0.push_back(0);
    v1.push_back(0);

    v0.push_back(1);
    v1.push_back(1);

    v0.push_back(2);
    v1.push_back(2);

    v0.push_back(0);
    v1.push_back(0);

    v0.push_back(1);
    v1.push_back(1);

    for (std::size_t i = 0; i < 251; ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    for (std::size_t i = 253; i < 300; ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    {
        auto it1 = v1.begin();
        auto it2 = v0.begin();

        auto itend = v1.end();

        while (it1 != itend) {
            auto const val1 = *it1;
            auto const val2 = *it2;
            EXPECT_EQUAL(val1, val2);
            ++it1;
            ++it2;
        }
    }

    {
        auto it1 = v1.end();
        auto it2 = v0.end();

        auto itend = v1.begin();
        while (it1 != itend) {
            --it1;
            --it2;
            auto const val1 = *it1;
            auto const val2 = *it2;
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        auto it1 = v0.begin();
        auto it12 = v0.begin();
        auto it2 = v1.begin();
        auto it22 = v1.begin();

        v0.resize(220);
        v1.resize(220);

        auto it13 = v0.begin();
        auto it23 = v1.begin();

        EXPECT_EQUAL(it1, it12);
        EXPECT_EQUAL(it1, it13);

        EXPECT_EQUAL(it2, it22);
        EXPECT_EQUAL(it2, it23);
    }

    for (std::size_t i = 0; i < v1.size(); ++i) {
        EXPECT_EQUAL(v1[i], v0[i]);
    }
}

void
beyond16_random_access_tests()
{
    std::size_t iters = 100000;
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;

    for (std::size_t i = 0; i < (iters); ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    EXPECT_EQUAL(v0.size(), v1.size());

    for (std::size_t i = 0; i < (iters); ++i) {
        auto val0 = v0[i];
        auto val1 = v1[i];
        EXPECT_EQUAL(val0, val1);
    }
}

void
beyond16_resize_tests()
{
    std::size_t iters = 67000;
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;

    for (std::size_t i = 0; i < iters; ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    std::size_t iters_end = 65000;
    for (std::size_t i = iters; i-- > iters_end;) {
        v0.resize(i);
        v1.resize(i);
        EXPECT_EQUAL(v0.size(), v1.size());
        for (std::size_t j = 0; j < i; ++j) {
            auto val0 = v0[j];
            auto val1 = v1[j];
            EXPECT_EQUAL(val0, val1);
        }
    }
}

void
iterator_tests()
{
    std::size_t iters = 100000;
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;

    for (std::size_t i = 0; i < (iters); ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    auto it0 = v0.begin();
    auto it1 = v1.begin();

    for (std::size_t i = 0; i < (iters); ++i) {
        auto val0 = *it0;
        auto val1 = *it1;
        EXPECT_EQUAL(val0, val1);
        ++it0;
        ++it1;
    }

    for (std::size_t i = 0; i < (iters); ++i) {
        --it0;
        --it1;
        auto val0 = *it0;
        auto val1 = *it1;
        EXPECT_EQUAL(val0, val1);
    }
}

void
iterator_bwd_tests()
{
    std::size_t iters = 100000;
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;

    for (std::size_t i = 0; i < (iters); ++i) {
        v0.push_back(i);
        v1.push_back(i);
    }

    auto it0 = v0.end();
    auto it1 = v1.end();

    for (std::size_t i = 0; i < (iters); ++i) {
        --it0;
        --it1;
        auto val0 = *it0;
        auto val1 = *it1;
        EXPECT_EQUAL(val0, val1);
    }
}

void
timing_tests()
{
    std::vector<std::size_t> v0;
    adhoc::vector_idx v1;

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
            v0.push_back(i);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector filling:" << std::endl;
        std::cout << time << std::endl;
    }

    {
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            v1.push_back(i);
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 filling:" << std::endl;
        std::cout << time << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += v0[i];
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += v1[i];
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += *it;
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.begin();
        for (std::size_t i = 0; i < iters; ++i) {
            sum += *it;
            ++it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it fwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v0.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += *it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "normal vector it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }

    {
        std::size_t sum = 0;
        auto time1 = std::chrono::high_resolution_clock::now();
        auto it = v1.end();
        for (std::size_t i = 0; i < iters; ++i) {
            --it;
            sum += *it;
        }
        auto time2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
        std::cout << "index vector 1 it bwd sum:" << std::endl;
        std::cout << time << std::endl;
        std::cout << sum << std::endl;
    }
}

} // namespace

auto
main() -> int
{
    simple_tests();
    beyond16_random_access_tests();
    beyond16_resize_tests();
    iterator_tests();
    iterator_bwd_tests();
    timing_tests();

    TEST_END;
}
