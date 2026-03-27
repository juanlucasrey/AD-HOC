#ifndef BMC_TOOLS_DATES
#define BMC_TOOLS_DATES

#include <chrono>
#include <flat_set>
#include <vector>

using namespace std::chrono;

namespace adhoc {

template<class T>
inline auto
operator+(year_month_day const& lhs, T d) -> year_month_day
{
    return std::chrono::sys_days(lhs) + days{ d };
}

inline auto
operator-(year_month_day const& lhs, year_month_day const& rhs) -> int
{
    auto const diff = std::chrono::sys_days(lhs) - std::chrono::sys_days(rhs);
    return std::chrono::duration_cast<std::chrono::days>(diff).count();
}

template<class T>
inline auto
operator+=(year_month_day& lhs, T d) -> year_month_day&
{
    lhs = lhs + d;
    return lhs;
}

template<class Float = double>
auto
convert_to_dt(year_month_day const& start_date, year_month_day const& end_date)
{
    constexpr Float DAYS_A_YEAR_ACT36525 = 365.25;
    return static_cast<Float>(end_date - start_date) / DAYS_A_YEAR_ACT36525;
}

template<class Float = double>
auto
convert_to_dt(year_month_day const& asofdate, std::vector<year_month_day> const& dates)
{
    std::vector<Float> result(dates.size());

    for (std::size_t i = 0; i < dates.size(); ++i) {
        result[i] = convert_to_dt(asofdate, dates[i]);
    }

    return result;
}

template<class Float = double>
auto
convert_to_dt_360(year_month_day const& start_date, year_month_day const& end_date)
{
    constexpr Float DAYS_A_YEAR_ACT360 = 360.;
    return static_cast<Float>(end_date - start_date) / DAYS_A_YEAR_ACT360;
}

class Schedule final {
  public:
    Schedule(std::flat_set<year_month_day> const& d)
      : dates(d)
    {
    }

    Schedule(std::vector<year_month_day> const& d)
      : dates(d)
    {
    }

    inline void merge(Schedule const& other) { this->dates.insert(other.dates.cbegin(), other.dates.cend()); }

    inline auto get_dates() const -> std::flat_set<year_month_day> const& { return this->dates; }

    inline auto get_dates_vec() const -> std::vector<year_month_day>
    {
        std::vector<year_month_day> res(this->dates.cbegin(), this->dates.cend());
        return res;
    }

    inline auto operator[](std::size_t i) const -> year_month_day const& { return *std::next(this->dates.cbegin(), i); }

    inline auto size() const -> std::size_t { return this->dates.size(); }

    inline auto front() const -> year_month_day const& { return *this->dates.begin(); }

    inline auto back() const -> year_month_day const& { return *this->dates.rbegin(); }

    template<class Float = double>
    auto convert_to_double() const -> std::flat_set<Float>
    {
        std::flat_set<Float> result;

        constexpr Float DAYS_A_YEAR_ACT36525_FRAC = 1.0 / 365.25;

        auto const& init = *this->dates.cbegin();

        for (auto it = this->dates.cbegin(); it != this->dates.cend(); ++it) {
            result.emplace_hint(result.cend(), static_cast<Float>(*it - init) * DAYS_A_YEAR_ACT36525_FRAC);
        }

        return result;
    }

  private:
    std::flat_set<year_month_day> dates;
};

template<class T>
auto
get(std::flat_set<T> in, std::size_t i) -> T const&
{
    return *std::next(in.begin(), i);
}

enum class HolidayName {
    Saturday,
    Sunday,
    MartinLutherKingBirhtday,
    WashingtonBirthday,
};

template<HolidayName Hol>
class HolidayChecker {
  public:
    static auto is_holiday(std::chrono::year_month_day const& date) -> bool
    {
        if constexpr (Hol == HolidayName::Saturday) {
            return weekday{ date } == std::chrono::Saturday;
        }
        else if constexpr (Hol == HolidayName::Sunday) {
            return weekday{ date } == std::chrono::Sunday;
        }
        else if constexpr (Hol == HolidayName::MartinLutherKingBirhtday) {
            // Martin Luther King Birthday is on the third monday of January
            // so it can't be on or before the 14th or on or after the 22st
            if (date.month() != std::chrono::January) {
                return false;
            }

            if (date.day() <= day{ 14 } || date.day() >= day{ 22 }) {
                return false;
            }

            // find the third monday on January, on that date's year
            // and check if our date is that one
            auto third_monday_january = year_month_day{ Monday[3] / January / date.year() };

            return third_monday_january == date;
        }
        else if constexpr (Hol == HolidayName::WashingtonBirthday) {
            // Washington Birthday is on the third monday of February
            // so it can't be on or before the 14th or on or after the 22st
            if (date.month() != std::chrono::February) {
                return false;
            }

            if (date.day() <= day{ 14 } || date.day() >= day{ 22 }) {
                return false;
            }

            // find the third monday on February, on that date's year
            // and check if our date is that one
            auto third_monday_february = year_month_day{ Monday[3] / February / date.year() };

            return third_monday_february == date;
        }
        else {
            return false;
        }
    }
};

template<HolidayName... Hols>
class Calendar {
  public:
    static auto is_holiday(std::chrono::year_month_day const& date) -> bool
    {
        return (HolidayChecker<Hols>::is_holiday(date) || ...);
    }
};

using NoCalendar = Calendar<>;
using FD = Calendar<HolidayName::Saturday,
                    HolidayName::Sunday,
                    HolidayName::MartinLutherKingBirhtday,
                    HolidayName::WashingtonBirthday>;

} // namespace adhoc

#endif // BMC_TOOLS_DATES
