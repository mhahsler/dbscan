// Copyright (c) 2026 Michael Hahsler.
// Distributed under the GNU General Public License, version 3 or later.

#ifndef DBSCAN_ANN_RAII_H
#define DBSCAN_ANN_RAII_H

#include "ANN/ANN.h"
#include <cstddef>
#include <memory>
#include <stdexcept>

/* Resource Acquisition Is Initialization (RAII) ownership Model to make sure
 * memory is unwound in case of interruptions.
 *
 * ANN uses points and trees.
 *
 * ANNPoints owns two arrays through unique_ptr:
 *  - One contiguous array containing all coordinates.
 *  - An array of row pointers pointing into those coordinates.
 * When the ANNPoints gets deleted, point storage is automatically released.
 *
 * The ANN tree borrows access to point pointers without taking ownership.
 * Its destructor frees its own tree nodes and bookkeeping, while
 *
 * Important: Declaration order ensures safe cleanup
 * ANNPoints dataPts(nrow, ncol);
 * std::unique_ptr<ANNpointSet> kdTree = make_ann_tree(dataPts.get(), ...);
 *
 * C++ destroys local objects in reverse declaration order:
 *  1. kdTree destroys the tree.
 *  2. dataPts destroys the row pointers and coordinates.
 *
 * This also happens when a C++ exception or Rcpp interrupt unwinds the call.
 * If point-storage construction fails halfway through, the already-created
 * unique_ptr members release their allocations automatically.
 *
 * Explicit query points have their own ANNPoints owner.
 * All storage stays alive until the search finishes.

 * annClose() remains separate: it frees ANN’s shared global structure at
 * package unload, after the individual trees are gone.
 */
class ANNPoints {
  std::unique_ptr<ANNcoord[]> coordinates_;
  std::unique_ptr<ANNpoint[]> points_;

public:
  ANNPoints() = default;

  ANNPoints(int rows, int columns) {
    if (rows < 1 || columns < 1)
      throw std::invalid_argument("ANN point storage requires positive dimensions");
    const std::size_t n = static_cast<std::size_t>(rows);
    const std::size_t d = static_cast<std::size_t>(columns);
    if (n > static_cast<std::size_t>(-1) / d / sizeof(ANNcoord))
      throw std::length_error("ANN point storage is too large");
    coordinates_ = std::make_unique<ANNcoord[]>(n * d);
    points_ = std::make_unique<ANNpoint[]>(n);
    for (std::size_t i = 0; i < n; ++i)
      points_[i] = coordinates_.get() + i * d;
  }

  ANNPoints(const ANNPoints&) = delete;
  ANNPoints& operator=(const ANNPoints&) = delete;
  ANNPoints(ANNPoints&&) = default;
  ANNPoints& operator=(ANNPoints&&) = default;

  ANNpointArray get() const { return points_.get(); }
  ANNpoint operator[](std::size_t i) const { return points_[i]; }
};

inline std::unique_ptr<ANNpointSet> make_ann_tree(
    ANNpointArray points, int rows, int columns, int type,
    int bucket_size, int split_rule) {
  if (type == 1)
    return std::make_unique<ANNkd_tree>(points, rows, columns, bucket_size,
      static_cast<ANNsplitRule>(split_rule));
  return std::make_unique<ANNbruteForce>(points, rows, columns);
}

#endif
