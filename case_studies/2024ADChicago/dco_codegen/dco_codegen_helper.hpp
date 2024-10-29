// It defines aligned allocation and deallocation routines as well as
// the datatypes used for instantiating the generated code.
// We need header guards here, since we potentially have two files
// with identical content; pragma once doesn't work in that case.
#ifndef DCO_CODEGEN_HELPER_HPP
#define DCO_CODEGEN_HELPER_HPP

namespace dco {
namespace internal {
// Allocate 64 byte aligned memory; round size up to next multiple of alignment.
inline void* codegen_alloc(std::size_t size) {
  std::size_t alignment = 64;
  size = static_cast<std::size_t>(std::ceil(static_cast<double>(size) / static_cast<double>(alignment))) * alignment;
#if defined(_MSC_VER)
  return _aligned_malloc(size, alignment);
#else
  // posix_memalign is supported on older glibc versions.
  void* ret;
  if (posix_memalign(&ret, alignment, size) != 0) ret = nullptr;
  return ret;
#endif
}
// Deallocate aligned memory.
inline void codegen_dealloc(void* buffer) {
  if (buffer == NULL) return;
#if defined(_MSC_VER)
  _aligned_free(buffer);
#else
  free(buffer);
#endif
}
}
// Holds (and reuses) buffer for intermediates/stack used in generated code.
struct codegen_buffer_t {
  void* v  = nullptr;
  void* dv = nullptr;
  void* stack = nullptr;
  void* condition = nullptr;
  std::size_t size_v = 0, size_dv = 0, size_stack = 0, size_condition = 0;

  void _allocate(std::size_t size, void *& p, std::size_t & allocated_size) {
    if (size > allocated_size) {
      if (p != nullptr) {
        internal::codegen_dealloc(p);
        allocated_size = 0;
      }
      p = internal::codegen_alloc(size);
      allocated_size = size;
    }
  }

  void allocate_v(std::size_t size)     { _allocate(size, v, size_v); }
  void allocate_dv(std::size_t size)    { _allocate(size, dv, size_dv); }
  void allocate_stack(std::size_t size) { _allocate(size, stack, size_stack); }
  void allocate_condition(std::size_t size) { _allocate(size, condition, size_condition); }

  void deallocate() {
    if (v) {
      internal::codegen_dealloc(v);
      v = nullptr;
    }
    if (dv) {
      internal::codegen_dealloc(dv);
      dv = nullptr;
    }
    if (stack) {
      internal::codegen_dealloc(stack);
      stack = nullptr;
    }
    if (condition) {
      internal::codegen_dealloc(condition);
      condition = nullptr;
    }
  }
  
  ~codegen_buffer_t() { deallocate(); }
};
}
#endif
