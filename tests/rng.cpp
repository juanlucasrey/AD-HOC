// #ifndef RNG_HPP
// #define RNG_HPP

#include <iostream>
using namespace std;

namespace adhoc {

template <class Source> class A {
    Source &m_x;

  public:
    A(Source &x) : m_x(x) { cout << "A's constructor called " << endl; }
};

class B {
    static A<int> a;

  public:
    B() { cout << "B's constructor called " << endl; }
    static A<int> getA() { return a; }
};

int temp = 1;
A<int> B::a(temp); // definition of a

int main() {
    B b1, b2, b3;
    A<int> a = b1.getA();

    return 0;
}

} // namespace adhoc
  // #endif // RNG_HPP