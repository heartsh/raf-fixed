/////////////////////////////////////////////////////////////////
// ScoringScheme.hpp
/////////////////////////////////////////////////////////////////

#ifndef SCORINGSCHEME_HPP
#define SCORINGSCHEME_HPP

#include "IO.hpp"
#include "Options.hpp"
#include "SparseMatrix.hpp"
#include "Utilities.hpp"

class ScoringScheme
{
    const std::vector<SparseMatrix<Real> > &base_pairing_posteriors;
    const std::vector<SparseMatrix<Real> > &alignment_posteriors;
    const std::vector<int> &ids_A;
    const std::vector<int> &ids_B;
    const std::vector<Sequence> &projected_alignment;
    const std::vector<Real> &weights;
    
    const int LA;
    const int LB;   
    
    SparseMatrix<Real> folding_A_mask_sparse;
    SparseMatrix<Real> folding_B_mask_sparse;
    SparseMatrix<Real> alignment_mask_sparse;

    std::vector<Real> folding_A_matrix;
    std::vector<Real> folding_B_matrix;
    std::vector<Real> alignment_matrix;

    Real folding_A_offset;
    Real folding_B_offset;
    Real alignment_offset;

    std::vector<Real> Ramp(const std::vector<Real> &x, const Real lower, const Real upper) const;

    void ComputeFoldingMatrix(SparseMatrix<Real> &folding_mask_sparse,
                              std::vector<Real> &folding_matrix,
                              Real &folding_offset,
                              const std::vector<int> &ids);
    
    void ComputeAlignmentMatrix(SparseMatrix<Real> &alignment_mask_sparse,
                                std::vector<Real> &alignment_matrix,
                                Real &alignment_offset,
                                const std::vector<int> &ids_A,
                                const std::vector<int> &ids_B);

    void ExpandFoldingPosteriors(const int id,
                                 const std::vector<int> &positions,
                                 std::vector<Real> &P_paired,
                                 std::vector<Real> &P_unpaired,
                                 std::vector<Real> &unsparse_mask) const;

    void ExpandAlignmentPosteriors(const int id_A,
                                   const int id_B,
                                   const std::vector<int> &positions_A,
                                   const std::vector<int> &positions_B,
                                   std::vector<Real> &P_aligned,
                                   std::vector<Real> &P_unaligned_A,
                                   std::vector<Real> &P_unaligned_B,
                                   std::vector<Real> &unsparse_mask) const;

    std::pair<std::vector<Real>, Real> ComputeFoldingFeature(const int index,
                                                             const std::vector<Real> &P_paired,
                                                             const std::vector<Real> &P_unpaired) const;

    std::pair<std::vector<Real>, Real> ComputeAlignmentFeature(const int index,
                                                               const int id_A,
                                                               const int id_B,
                                                               const std::vector<Real> &P_aligned,
                                                               const std::vector<Real> &P_unaligned_A,
                                                               const std::vector<Real> &P_unaligned_B) const;
    
    void AugmentWithFoldingLoss(std::vector<Real> &folding_matrix,
                                Real &folding_offset,
                                const std::vector<int> &mapping,
                                const std::vector<int> &ids,
                                const Real folding_FP_margin,
                                const Real folding_FN_margin);
    
    void AugmentWithAlignmentLoss(std::vector<Real> &alignment_matrix,
                                  Real &alignment_offset,
                                  const std::vector<int> &aligned_to_A,
                                  const std::vector<int> &aligned_to_B,
                                  const Real alignment_FP_margin,
                                  const Real alignment_FN_margin);

    std::vector<Real> ComputeParseFoldingFeatures(const std::vector<int> &mapping,
                                                  const std::vector<int> &ids) const;
    
    std::vector<Real> ComputeParseAlignmentFeatures(const std::vector<int> &aligned_to_A,
                                                    const std::vector<int> &aligned_to_B) const;
    
    std::vector<int> GetGroupPositions(const std::vector<Sequence> &alignment, 
                                       const std::vector<int> &ids) const;
    
public:
    
    ScoringScheme(const std::vector<SparseMatrix<Real> > &base_pairing_posteriors,
                  const std::vector<SparseMatrix<Real> > &alignment_posteriors,
                  const std::vector<int> &ids_A,
                  const std::vector<int> &ids_B,
                  const std::vector<Sequence> &projected_alignment,
                  const std::vector<Real> &weights);

    ScoringScheme(const ScoringScheme &rhs);

    void AugmentWithLoss(const std::vector<int> &true_mapping_A,
                         const std::vector<int> &true_mapping_B,
                         const std::vector<int> &true_aligned_to_A,
                         const std::vector<int> &true_aligned_to_B,
                         const Real folding_FP_margin,
                         const Real folding_FN_margin,
                         const Real alignment_FP_margin,
                         const Real alignment_FN_margin);
    
    std::vector<Real> ComputeFeatures(const std::vector<int> &mapping_A,
                                      const std::vector<int> &mapping_B,
                                      const std::vector<int> &aligned_to_A,
                                      const std::vector<int> &aligned_to_B) const;
    
    const SparseMatrix<Real> &GetFoldingAMaskSparse() const { return folding_A_mask_sparse; }
    const SparseMatrix<Real> &GetFoldingBMaskSparse() const { return folding_B_mask_sparse; }
    const SparseMatrix<Real> &GetAlignmentMaskSparse() const { return alignment_mask_sparse; }

    const std::vector<Real> &GetFoldingAMatrix() const { return folding_A_matrix; }
    const std::vector<Real> &GetFoldingBMatrix() const { return folding_B_matrix; }
    const std::vector<Real> &GetAlignmentMatrix() const { return alignment_matrix; }

    const std::vector<Real> &GetWeights() const { return weights; }

    Real GetOffset() const { return folding_A_offset + folding_B_offset + alignment_offset; }
};

#endif
