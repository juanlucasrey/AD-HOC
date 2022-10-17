#ifndef ADHOC_HPP
#define ADHOC_HPP

namespace adhoc {

template <class Input> class mul_scalar;

template <class Input> class add_scalar {
    double m_scalar{0.};
    Input &m_active;

  public:
    add_scalar(double scalar, Input &active);
    auto operator+(double rhs) const -> add_scalar<const add_scalar<Input>>;
    auto operator*(double rhs) const -> mul_scalar<const add_scalar<Input>>;
    operator double() const noexcept;

    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

template <class Input>
add_scalar<Input>::add_scalar(double scalar, Input &active)
    : m_scalar(scalar), m_active(active) {}
template <class Input>
auto add_scalar<Input>::operator+(double rhs) const
    -> add_scalar<const add_scalar<Input>> {
    return {rhs, *this};
}
template <class Input>
auto add_scalar<Input>::operator*(double rhs) const
    -> mul_scalar<const add_scalar<Input>> {
    return {rhs, *this};
}
template <class Input> add_scalar<Input>::operator double() const noexcept {
    return this->m_scalar + this->m_active.operator double();
}

template <class Input>
template <class Denom>
inline auto add_scalar<Input>::d(Denom &in) const noexcept -> double {
    return this->m_active.d(in);
}

template <class Input> class mul_scalar {
    double m_scalar{0.};
    Input &m_active;

  public:
    mul_scalar(double scalar, Input &active);

    auto operator+(double rhs) const -> add_scalar<const mul_scalar<Input>>;
    auto operator*(double rhs) const -> mul_scalar<const mul_scalar<Input>>;
    operator double() const noexcept;

    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

template <class Input>
mul_scalar<Input>::mul_scalar(double scalar, Input &active)
    : m_scalar(scalar), m_active(active) {}
template <class Input>
auto mul_scalar<Input>::operator+(double rhs) const
    -> add_scalar<const mul_scalar<Input>> {
    return {rhs, *this};
}
template <class Input>
auto mul_scalar<Input>::operator*(double rhs) const
    -> mul_scalar<const mul_scalar<Input>> {
    return {rhs, *this};
}
template <class Input> mul_scalar<Input>::operator double() const noexcept {
    return this->m_scalar * this->m_active.operator double();
}

template <class Input>
template <class Denom>
inline auto mul_scalar<Input>::d(Denom &in) const noexcept -> double {
    return this->m_scalar * this->m_active.d(in);
}

class adouble {
    double m_value{0};

  public:
    adouble(double value);
    auto operator+(double rhs) const -> add_scalar<const adouble>;
    auto operator*(double rhs) const -> mul_scalar<const adouble>;

    operator double() const noexcept;

    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

inline adouble::adouble(double value) : m_value(value) {}
inline auto adouble::operator+(double rhs) const -> add_scalar<const adouble> {
    return {rhs, *this};
}
inline auto adouble::operator*(double rhs) const -> mul_scalar<const adouble> {
    return {rhs, *this};
}
inline adouble::operator double() const noexcept { return this->m_value; }

template <class Denom>
inline auto adouble::d(Denom & /* in */) const noexcept -> double {
    return 0.;
}

template <> inline auto adouble::d(adouble &in) const noexcept -> double {
    return static_cast<double>(&in == this);
}

} // namespace adhoc

#endif // ADHOC_HPP