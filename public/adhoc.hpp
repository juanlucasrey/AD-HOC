#ifndef ADHOC_HPP
#define ADHOC_HPP

namespace adhoc {

class Base {
  protected:
    double m_value{0};

  public:
    explicit Base(double value);
    [[nodiscard]] auto v() const noexcept -> double;
};

Base::Base(double value) : m_value(value) {}

inline auto Base::v() const noexcept -> double { return this->m_value; }

template <class Input> class mul_scalar;
template <class Input> class add_scalar;
template <class Input1, class Input2> class mul_active;
template <class Input1, class Input2> class add_active;

template <class Input> class add_scalar : public Base {
    double m_scalar{0.};
    Input &m_active;

  public:
    add_scalar(double scalar, Input &active);
    auto operator+(double rhs) const -> add_scalar<const add_scalar<Input>>;
    auto operator*(double rhs) const -> mul_scalar<const add_scalar<Input>>;
    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

template <class Input>
add_scalar<Input>::add_scalar(double scalar, Input &active)
    : Base(scalar + active.v()), m_scalar(scalar), m_active(active) {}

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

template <class Input>
template <class Denom>
inline auto add_scalar<Input>::d(Denom &in) const noexcept -> double {
    return this->m_active.d(in);
}

template <class Input> class mul_scalar : public Base {
    double m_scalar{0.};
    Input &m_active;

  public:
    mul_scalar(double scalar, Input &active);

    auto operator+(double rhs) const -> add_scalar<const mul_scalar<Input>>;
    auto operator*(double rhs) const -> mul_scalar<const mul_scalar<Input>>;
    template <class Arg>
    auto operator+(Arg &rhs) const -> add_active<const mul_scalar<Input>, Arg>;

    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

template <class Input>
mul_scalar<Input>::mul_scalar(double scalar, Input &active)
    : Base(scalar * active.v()), m_scalar(scalar), m_active(active) {}

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

template <class Input>
template <class Arg>
inline auto mul_scalar<Input>::operator+(Arg &rhs) const
    -> add_active<const mul_scalar<Input>, Arg> {
    return {*this, rhs};
}

template <class Input>
template <class Denom>
inline auto mul_scalar<Input>::d(Denom &in) const noexcept -> double {
    return this->m_scalar * this->m_active.d(in);
}

template <class Input1, class Input2> class add_active : public Base {
    Input1 &m_active1;
    Input2 &m_active2;

  public:
    add_active(Input1 &active1, Input2 &active2);
    template <class Denom>
    [[nodiscard]] auto d(Denom &in) const noexcept -> double;
};

template <class Input1, class Input2>
add_active<Input1, Input2>::add_active(Input1 &active1, Input2 &active2)
    : Base(active1.v() + active2.v()), m_active1(active1), m_active2(active2) {}

template <class Input1, class Input2>
template <class Denom>
inline auto add_active<Input1, Input2>::d(Denom &in) const noexcept -> double {
    return this->m_active1.d(in) + this->m_active2.d(in);
}

template <class Input1, class Input2> class mul_active : public Base {
    Input1 &m_active1;
    Input2 &m_active2;

  public:
    mul_active(Input1 &active1, Input2 &active2);
    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

template <class Input1, class Input2>
mul_active<Input1, Input2>::mul_active(Input1 &active1, Input2 &active2)
    : Base(active1.v() * active2.v()), m_active1(active1), m_active2(active2) {}

template <class Input1, class Input2>
template <class Denom>
inline auto mul_active<Input1, Input2>::d(Denom &in) const noexcept -> double {
    return this->m_active1.d(in) * this->m_active2.v() +
           this->m_active2.d(in) * this->m_active1.v();
}

class adouble {
    double m_value{0};

  public:
    explicit adouble(double value);
    auto operator+(double rhs) const -> add_scalar<const adouble>;
    auto operator*(double rhs) const -> mul_scalar<const adouble>;

    template <class Arg>
    auto operator+(Arg &rhs) const -> add_active<const adouble, Arg>;

    template <class Arg>
    auto operator*(Arg &rhs) const -> mul_active<const adouble, Arg>;

    [[nodiscard]] auto v() const noexcept -> double;

    template <class Denom> auto d(Denom &in) const noexcept -> double;
};

inline adouble::adouble(double value) : m_value(value) {}

inline auto adouble::operator+(double rhs) const -> add_scalar<const adouble> {
    return {rhs, *this};
}

inline auto adouble::operator*(double rhs) const -> mul_scalar<const adouble> {
    return {rhs, *this};
}

template <class Arg>
inline auto adouble::operator+(Arg &rhs) const
    -> add_active<const adouble, Arg> {
    return {*this, rhs};
}

template <class Arg>
inline auto adouble::operator*(Arg &rhs) const
    -> mul_active<const adouble, Arg> {
    return {*this, rhs};
}

inline auto adouble::v() const noexcept -> double { return this->m_value; }

template <class Denom>
inline auto adouble::d(Denom & /* in */) const noexcept -> double {
    return 0.;
}

template <> inline auto adouble::d(adouble &in) const noexcept -> double {
    return static_cast<double>(&in == this);
}

template <class Input>
inline mul_scalar<const Input> operator*(double scalar, Input &in) {
    return in * scalar;
}

template <class Input>
inline add_scalar<Input> operator+(double scalar, Input &in) {
    return in + scalar;
}

} // namespace adhoc

#endif // ADHOC_HPP