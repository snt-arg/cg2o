#ifndef CG2O_EIGEN_LINEAR_SOLVER_CHECKS_H
#define CG2O_EIGEN_LINEAR_SOLVER_CHECKS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>

namespace cg2o::internal {
// All dense backends assemble the same symmetric matrix from upper blocks.
// setZero(rows, cols) reuses the allocation when the dimensions are unchanged.
template <class BlockMatrix, class DenseMatrix>
void fillSymmetricDenseMatrix(const BlockMatrix &blocks, DenseMatrix &matrix) {
  matrix.setZero(blocks.rows(), blocks.cols());
  for (std::size_t col = 0; col < blocks.blockCols().size(); ++col) {
    const int c = blocks.colBaseOfBlock(col);
    for (const auto &[row, block] : blocks.blockCols()[col]) {
      if (row > static_cast<int>(col))
        continue;
      const int r = blocks.rowBaseOfBlock(row);
      matrix.block(r, c, block->rows(), block->cols()) = *block;
      if (row != static_cast<int>(col))
        matrix.block(c, r, block->cols(), block->rows()) = block->transpose();
    }
  }
}

// Dense fallback solvers must not report success for NaN answers or an
// inconsistent least-squares solution of a square Newton/KKT system.
template <class Matrix, class Solution, class RightHandSide>
bool acceptableDenseSolution(const Matrix &matrix, const Solution &x,
                             const RightHandSide &b) {
  if (!matrix.allFinite() || !x.allFinite() || !b.allFinite())
    return false;
  if (matrix.rows() == 0)
    return true;
  const double residual = (matrix * x - b).cwiseAbs().maxCoeff();
  const double scale =
      matrix.cwiseAbs().rowwise().sum().maxCoeff() * x.cwiseAbs().maxCoeff() +
      b.cwiseAbs().maxCoeff();
  return std::isfinite(residual) && std::isfinite(scale) &&
         residual <= 1e-8 * scale;
}

template <class SparseMatrix, class Solution, class RightHandSide>
bool acceptableSymmetricSparseSolution(const SparseMatrix &matrix,
                                       const Solution &x,
                                       const RightHandSide &b) {
  if (!x.allFinite() || !b.allFinite() ||
      !Eigen::Map<const Eigen::VectorXd>(matrix.valuePtr(), matrix.nonZeros())
           .allFinite())
    return false;
  Eigen::VectorXd rowSums = Eigen::VectorXd::Zero(matrix.rows());
  for (int col = 0; col < matrix.outerSize(); ++col)
    for (typename SparseMatrix::InnerIterator it(matrix, col); it; ++it) {
      rowSums[it.row()] += std::abs(it.value());
      if (it.row() != it.col())
        rowSums[it.col()] += std::abs(it.value());
    }
  const double residual =
      (matrix.template selfadjointView<Eigen::Upper>() * x - b)
          .cwiseAbs()
          .maxCoeff();
  const double scale =
      rowSums.maxCoeff() * x.cwiseAbs().maxCoeff() + b.cwiseAbs().maxCoeff();
  return std::isfinite(residual) && std::isfinite(scale) &&
         residual <= 1e-8 * scale;
}
} // namespace cg2o::internal
#endif
