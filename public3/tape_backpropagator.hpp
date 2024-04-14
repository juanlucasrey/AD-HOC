#ifndef ADHOC3_TAPE_BACKPROPAGATOR_HPP
#define ADHOC3_TAPE_BACKPROPAGATOR_HPP

// #include "tape2.hpp"
#include <tape_buffer_size.hpp>

#include <array>

namespace adhoc3 {

template <class NodesValue, class NodesDerivative, class OutputDerivative>
class TapeBackpropagator {
  public:
    explicit TapeBackpropagator(NodesValue /* vn */, NodesDerivative /* dn */,
                                OutputDerivative /* od */) {}

    template <class CalcTreeValue, class InterfaceDerivatives,
              class ArrayInterfaceDerivatives>
    void backpropagate(CalcTreeValue /* ct */, InterfaceDerivatives /* id */,
                       ArrayInterfaceDerivatives /* aid */) {}

  private:
    std::array<double, tape_buffer_size(NodesValue{}, NodesDerivative{},
                                        OutputDerivative{})>
        m_buffer{};
};

} // namespace adhoc3

#endif // ADHOC3_TAPE_BACKPROPAGATOR_HPP
