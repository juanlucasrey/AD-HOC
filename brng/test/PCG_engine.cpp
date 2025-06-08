#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/pcg_basic.h"

#include "external/pcg-cpp/pcg_random.hpp"

#include <PCG_engine.hpp>

#include <tools/mask.hpp>

auto main() -> int {

    std::size_t const sims = 10;
    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         355248013,  41705475,   3406281715, 4186697710, 483882979,
    //         2766312848, 1713261421, 154902030,  3085534493, 3877580365};
    //     adhoc::PCG_engine rng(11906719003534950848UL,
    //     15726070495360670683UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng();
    //         auto val2 = pcg32_random();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val2);
    //         }
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     adhoc::PCG_engine rng(11906719003534950848UL,
    //     15726070495360670683UL); check_fwd_and_back(rng, 1000000);
    //     adhoc::PCG_engine rng2(11906719003534950848UL,
    //     15726070495360670683UL); EXPECT_EQUAL(rng, rng2);
    // }

    // {
    //     adhoc::PCG_engine rng(11906719003534950848UL,
    //     15726070495360670683UL); check_back_and_fwd(rng, 1000000);
    //     adhoc::PCG_engine rng2(11906719003534950848UL,
    //     15726070495360670683UL); EXPECT_EQUAL(rng, rng2);
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         676697322, 420258633,  3418632178, 3595600211, 3265791279,
    //         257272927, 3607051826, 1330014364, 1691133457, 2692391003};

    //     pcg32 rng1;
    //     adhoc::PCG_engine rng2(16070087622772795188UL,
    //     1442695040888963407UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         676697322, 420258633,  3418632178, 3595600211, 3265791279,
    //         257272927, 3607051826, 1330014364, 1691133457, 2692391003};

    //     pcg32_oneseq rng1;
    //     adhoc::PCG_engine rng2(16070087622772795188UL,
    //     1442695040888963407UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         2951688802, 530767249, 3930644657, 3055241032, 2771306922,
    //         295818707,  462549900, 3822330474, 1885452019, 1306071709};

    //     pcg32_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs>
    //         rng2(7515307391605212643UL, 0UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         14951315693135216709ULL, 1541401459199960700ULL,
    //         3670514919227316241ULL,  11007308355854268502ULL,
    //         523514384104871782ULL,   18368464326245194742ULL,
    //         17856432868564374551ULL, 3312496325571258406ULL,
    //         2035507739627916964ULL,  10012490566475469317ULL};

    //     pcg64 rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         14951315693135216709ULL, 1541401459199960700ULL,
    //         3670514919227316241ULL,  11007308355854268502ULL,
    //         523514384104871782ULL,   18368464326245194742ULL,
    //         17856432868564374551ULL, 3312496325571258406ULL,
    //         2035507739627916964ULL,  10012490566475469317ULL};

    //     pcg64_oneseq rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         8495917065814552292ULL,  7993034796972095977ULL,
    //         17036305730449054549ULL, 1936042420643248620ULL,
    //         11695933135001989394ULL, 66427974694741102ULL,
    //         4766802956136386423ULL,  16424304722249225872ULL,
    //         1611499284811258992ULL,  14682659323449184743ULL};

    //     pcg64_fast rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
    //         rng2{14627392581883831783ULL, 0, 0, 0};

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    {
        std::vector<std::uint32_t> firstvals = {
            3315726463, 3197513375, 2339211044, 2677264775, 3075764312,
            1537135571, 968113393,  4136885025, 2226883328, 1193408306};

        pcg32_once_insecure rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
                          false>
            rng2{2109010154UL, 2891336453UL};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            if (i < firstvals.size()) {
                EXPECT_EQUAL(firstvals[i], val1);
            }
            auto val2 = rng2();
            if (val1 != val2) {
                std::cout << val1 << ", " << val2 << ", " << i << std::endl;
            }
            EXPECT_EQUAL(val1, val2);
        }
    }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         4073575256423186656ULL,  4976277049109602808ULL,
    //         1131369595235443094ULL,  8483047370033139638ULL,
    //         15929478269862897795ULL, 16801080801056684117ULL,
    //         13846753628636920712ULL, 37468615611691582ULL,
    //         15180297470358637079ULL, 14094073075520863173ULL};

    //     pcg64_once_insecure rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::rxs_m_xs,
    //                       false>
    //         rng2{16070087622772795188ULL, 1442695040888963407ULL};

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         14951315693135216709ULL, 17787361309651119173ULL,
    //         1541401459199960700ULL,  14418444932722016366ULL,
    //         3670514919227316241ULL,  13205094727360841602ULL,
    //         11007308355854268502ULL, 9622203189173594313ULL,
    //         523514384104871782ULL,   4081536189067969393ULL,
    //         18368464326245194742ULL, 13907861900781595380ULL,
    //         17856432868564374551ULL, 3064873753271191496ULL,
    //         3312496325571258406ULL,  4812683234354841000ULL,
    //         2035507739627916964ULL,  13413012600773839155ULL,
    //         10012490566475469317ULL, 13902356169348470650ULL};

    //     pcg128_once_insecure rng1;
    //     adhoc::PCG_engine<adhoc::uint128, 128,
    //     adhoc::tempering_type::xsl_rr_rr,
    //                       false>
    //         rng2(245720598905631564143578724636268694099_ULLL,
    //              117397592171526113268558934119004209487_ULLL);

    //     std::size_t j = 0;
    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         auto val2 = rng2();
    //         auto val1_low = static_cast<std::uint64_t>(val1);
    //         auto val2_low = static_cast<std::uint64_t>(val2);
    //         auto val1_high = static_cast<std::uint64_t>(val1 >> 64U);
    //         auto val2_high = static_cast<std::uint64_t>(val2 >> 64U);

    //         if (i < (firstvals.size() / 2)) {
    //             EXPECT_EQUAL(firstvals[j++], val1_low);
    //             EXPECT_EQUAL(firstvals[j++], val1_high);
    //         }

    //         EXPECT_EQUAL(val1_low, val2_low);
    //         EXPECT_EQUAL(val1_high, val2_high);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         3315726463, 3197513375, 2339211044, 2677264775, 3075764312,
    //         1537135571, 968113393,  4136885025, 2226883328, 1193408306};

    //     pcg32_oneseq_once_insecure rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
    //                       false>
    //         rng2(2109010154UL, 2891336453UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         464549341,  3420025356, 266817216, 2345318095, 3184842467,
    //         1693875720, 2361273422, 242334963, 1737642702, 3187197445};

    //     pcg32_k2 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng(16070087622772795188ULL, 1442695040888963407ULL);
    //     check_fwd_and_back(rng, 1000000);
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);
    //     EXPECT_EQUAL(rng, rng2);
    // }

    // {
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng(16070087622772795188ULL, 1442695040888963407ULL);
    //     check_back_and_fwd(rng, 1000000);
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);
    //     EXPECT_EQUAL(rng, rng2);
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         4092648421, 115971742,  1817535833, 1721513619, 959586403,
    //         60487604,   1684177250, 349045251,  2722158446, 1689624819};

    //     pcg32_k2_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 1, 32, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         3848381134, 2517649454, 720818121,  2595249164, 328196342,
    //         3423925721, 2991614995, 1647535439, 118720549,  3673685198};

    //     pcg32_k64 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 6, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         2074664189, 3522704816, 116223856, 3923406470, 3700224756,
    //         4103933789, 835371017,  296546934, 789310485,  3726250596};

    //     pcg32_k64_oneseq rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 6, 32, true, true, true, true>
    //         rng2(7515307391605212643ULL, 0);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         617584258,  4219051537, 2203214769, 2671283187, 1699695288,
    //         2926587264, 1375943208, 1298550895, 3617730191, 2973591187};

    //     pcg32_k64_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 6, 32, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < 1000000; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         1405310097, 2757682074, 105144385,  3985216533, 1417370028,
    //         3840045086, 3471053572, 1383254411, 2422373624, 1085342612};

    //     pcg32_c64 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 6, 16, false, false, false, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         4007654266, 3086871550, 2592107534, 1691134214, 1455892964,
    //         1658877237, 2433822242, 925810987,  724179278,  2986548162};

    //     pcg32_c64_oneseq rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 6, 32, false, false, false, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         3624681459, 1324265195, 4238366959, 4245115015, 2418776848,
    //         1436629646, 1333019554, 2074155484, 2811724518, 2501806264};

    //     pcg32_c64_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 6, 32, false, false, false, true>
    //         rng2(7515307391605212643ULL, 0);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         15039938985044866497ULL, 5683501356464782399ULL,
    //         144308531516313942ULL,   8099733370339180248ULL,
    //         8830420848678742797ULL,  784593829978671100ULL,
    //         13558795047510121175ULL, 4493558840608307740ULL,
    //         4638863308561733549ULL,  7273529665380397613ULL};

    //     pcg64_k32 rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 16>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         15039938985044866497ULL, 5683501356464782399ULL,
    //         144308531516313942ULL,   8099733370339180248ULL,
    //         8830420848678742797ULL,  784593829978671100ULL,
    //         13558795047510121175ULL, 4493558840608307740ULL,
    //         4638863308561733549ULL,  7273529665380397613ULL};

    //     pcg64_k32_oneseq rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 128>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         16965488128385771295ULL, 4035110450430870263ULL,
    //         11015960111203041517ULL, 2031968004906381767ULL,
    //         17032493867346909189ULL, 16859053264698729785ULL,
    //         3479986990094732936ULL,  14965926616117115728ULL,
    //         11012874564184962801ULL, 4921646156023041314ULL};

    //     pcg64_k32_fast rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 128, false, true>
    //         rng2{14627392581883831783ULL, 0, 0, 0};

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         2029614971409230076ULL,  6514341339642902115ULL,
    //         1112310771763100318ULL,  3214315263717706275ULL,
    //         13244405485125545269ULL, 11110384631416820350ULL,
    //         1574150008570857312ULL,  14567275504595343232ULL,
    //         221245292235614828ULL,   5350028896887104664ULL};

    //     pcg64_c32 rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 16, false, true, false>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         2029614971409230076ULL,  6514341339642902115ULL,
    //         1112310771763100318ULL,  3214315263717706275ULL,
    //         13244405485125545269ULL, 11110384631416820350ULL,
    //         1574150008570857312ULL,  14567275504595343232ULL,
    //         221245292235614828ULL,   5350028896887104664ULL};

    //     pcg64_c32_oneseq rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 128, false, true, false>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         12231051446020525062ULL, 10791781393604019997ULL,
    //         12136335159208195635ULL, 3689426444100948635ULL,
    //         14355892346896396477ULL, 2265510788907471068ULL,
    //         4379860483653266196ULL,  3691947383650620179ULL,
    //         13004466837528046427ULL, 14928106554271064862ULL};

    //     pcg64_c32_fast rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 5, 128, false, true, false>
    //         rng2{14627392581883831783ULL, 0, 0, 0};

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         1370724678, 2801319071, 1099474461, 3076628667, 3991166948,
    //         3088540082, 3255103682, 3587442063, 450804495,  278308096};

    //     pcg32_k1024 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 10, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         1701562813, 1156285936, 3156290616, 1679400757, 2391201473,
    //         997106114,  2720981435, 808610782,  881217557,  1529552752};

    //     pcg32_k1024_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 10, 32, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         2362955264, 879640684,  3759564574, 501056242,  1254933031,
    //         1447276272, 2812717641, 3764630258, 4046875901, 1121121318};

    //     pcg32_c1024 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 10, 16, false, false, false, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         748319511,  4135511899, 1533914471, 545970811,  4258627531,
    //         2151094910, 3834663917, 2337593352, 2389664239, 3459212721};

    //     pcg32_c1024_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 10, 32, false, false, false, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         8269303586148613052ULL,  13880561190893172669ULL,
    //         356522898681002235ULL,   15708926673483292903ULL,
    //         13207894129573472761ULL, 1499432260625917297ULL,
    //         2099356712721218521ULL,  2796489629428879846ULL,
    //         8488052876031346006ULL,  18109046441879503905ULL};

    //     pcg64_k1024 rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 10, 16>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         8269303586148613052ULL,  13880561190893172669ULL,
    //         356522898681002235ULL,   15708926673483292903ULL,
    //         13207894129573472761ULL, 1499432260625917297ULL,
    //         2099356712721218521ULL,  2796489629428879846ULL,
    //         8488052876031346006ULL,  18109046441879503905ULL};

    //     pcg64_k1024_fast rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 10, 128>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         8679299735327275590ULL,  15157656879479246069ULL,
    //         13809342237560597771ULL, 16074321457744229168ULL,
    //         14527115835751351749ULL, 6890878235662325439ULL,
    //         3854616527516944048ULL,  11224601051365220461ULL,
    //         14898189496744227855ULL, 2081511346038753580ULL};

    //     pcg64_c1024 rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 10, 16, false, true, false>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint64_t> firstvals = {
    //         8679299735327275590ULL,  15157656879479246069ULL,
    //         13809342237560597771ULL, 16074321457744229168ULL,
    //         14527115835751351749ULL, 6890878235662325439ULL,
    //         3854616527516944048ULL,  11224601051365220461ULL,
    //         14898189496744227855ULL, 2081511346038753580ULL};

    //     pcg64_c1024_fast rng1;
    //     adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
    //                       true, 10, 128, false, true, false>
    //         rng2;

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         1935467043, 3358424126, 3660229836, 4001473811, 966022768,
    //         1098839924, 2231243075, 1877617528, 262469621,  2936257308};

    //     pcg32_k16384 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 14, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     std::vector<std::uint32_t> firstvals = {
    //         747822173, 1078776174, 2974100056, 923354615, 2529832782,
    //         685402359, 3295016468, 1787336710, 394523013, 3567559088};

    //     pcg32_k16384_fast rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
    //                       true, 14, 32, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         if (i < firstvals.size()) {
    //             EXPECT_EQUAL(firstvals[i], val1);
    //         }

    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }
}
