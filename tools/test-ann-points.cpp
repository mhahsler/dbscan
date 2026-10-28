// Exercise partial-construction failures without requiring an R runtime.
// Copyright (c) 2026 Michael Hahsler. GPL version 3 or later.
#include "ann_raii.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <new>
#include <utility>

namespace {
int allocations_until_failure = -1;
int live_arrays = 0;
}

void* operator new[](std::size_t size) {
  if (allocations_until_failure == 0) throw std::bad_alloc();
  if (allocations_until_failure > 0) --allocations_until_failure;
  void* result = std::malloc(size ? size : 1);
  if (!result) throw std::bad_alloc();
  ++live_arrays;
  return result;
}

void operator delete[](void* pointer) noexcept {
  if (pointer) {
    --live_arrays;
    std::free(pointer);
  }
}

void operator delete[](void* pointer, std::size_t) noexcept {
  operator delete[](pointer);
}

int main() {
  // Fail each allocation in turn while creating both data and query storage.
  for (int fail_after = 0; fail_after < 4; ++fail_after) {
    allocations_until_failure = fail_after;
    bool failed = false;
    try {
      ANNPoints data(12, 3);
      ANNPoints query(4, 3);
    } catch (const std::bad_alloc&) {
      failed = true;
    }
    assert(failed);
    assert(live_arrays == 0);
  }
  allocations_until_failure = -1;
  {
    ANNPoints data(12, 3);
    data[2][1] = 42;
    ANNpointArray address = data.get();
    ANNPoints moved(std::move(data));
    assert(!data.get());
    assert(moved.get() == address && moved[2][1] == 42);
    ANNPoints assigned(2, 1);
    assigned = std::move(moved);
    assert(!moved.get());
    assert(assigned.get() == address && assigned[2][1] == 42);
  }
  assert(live_arrays == 0);
  try {
    ANNPoints invalid(0, 3);
    return 1;
  } catch (const std::invalid_argument&) {
    assert(live_arrays == 0);
  }
  std::cout << "ANNPoints allocation-failure and move checks passed.\n";
}
