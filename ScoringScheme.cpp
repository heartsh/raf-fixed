/////////////////////////////////////////////////////////////////
// ScoringScheme.cpp
/////////////////////////////////////////////////////////////////

#include "ScoringScheme.hpp"

/////////////////////////////////////////////////////////////////
// ScoringScheme::Ramp()
//
// Compute ramp basis function.
/////////////////////////////////////////////////////////////////

std::vector<Real> ScoringScheme::Ramp(const std::vector<Real> &x, const Real lower, const Real upper) const
{
    const Real diff = upper - lower;
    Assert(diff >= 0, "Bounds in ramp not in correct order.");
    
    std::vector<Real> ret(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        if (x[i] <= lower)
            ret[i] = Real(0);
        else if (x[i] >= upper)
            ret[i] = Real(1);
        else
            ret[i] = (x[i] - lower) / diff;
    }

    return ret;
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ScoringScheme()
//
// Constructor.
/////////////////////////////////////////////////////////////////

ScoringScheme::ScoringScheme(const std::vector<SparseMatrix<Real> > &base_pairing_posteriors,
                             const std::vector<SparseMatrix<Real> > &alignment_posteriors,
                             const std::vector<int> &ids_A,
                             const std::vector<int> &ids_B,
                             const std::vector<Sequence> &projected_alignment,
                             const std::vector<Real> &weights) :
    base_pairing_posteriors(base_pairing_posteriors),
    alignment_posteriors(alignment_posteriors),
    ids_A(ids_A),
    ids_B(ids_B),
    projected_alignment(projected_alignment),
    weights(weights),

    LA(int(projected_alignment[ids_A[0]].data.length()) - 1),
    LB(int(projected_alignment[ids_B[0]].data.length()) - 1)
{
    
#ifndef NDEBUG

    // check that the alignment given is really projected

    const std::vector<int> positions_A = GetGroupPositions(projected_alignment, ids_A);
    const std::vector<int> positions_B = GetGroupPositions(projected_alignment, ids_B);

    // gaps must be already compressed 
    
    Assert(int(positions_A.size()) == LA+1, "Alignment should already be projected.");
    Assert(int(positions_B.size()) == LB+1, "Alignment should already be projected.");

    for (int i = 1; i <= LA; i++)
        Assert(positions_A[i] == i, "Alignment should already be projected.");
    for (int i = 1; i <= LB; i++)
        Assert(positions_B[i] == i, "Alignment should already be projected.");

    // all sequences in each subset must have the same length

    for (size_t i = 1; i < ids_A.size(); i++)
        Assert(LA == int(projected_alignment[ids_A[i]].data.length()) - 1, "Not all sequences in the group have the same length.");
    for (size_t i = 1; i < ids_B.size(); i++)
        Assert(LB == int(projected_alignment[ids_B[i]].data.length()) - 1, "Not all sequences in the group have the same length.");
#endif

    ComputeFoldingMatrix(folding_A_mask_sparse, folding_A_matrix, folding_A_offset, ids_A);
    ComputeFoldingMatrix(folding_B_mask_sparse, folding_B_matrix, folding_B_offset, ids_B);
    ComputeAlignmentMatrix(alignment_mask_sparse, alignment_matrix, alignment_offset, ids_A, ids_B);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ScoringScheme()
//
// Copy constructor.
/////////////////////////////////////////////////////////////////

ScoringScheme::ScoringScheme(const ScoringScheme &rhs) :
    base_pairing_posteriors(rhs.base_pairing_posteriors),
    alignment_posteriors(rhs.alignment_posteriors),
    
    ids_A(rhs.ids_A),
    ids_B(rhs.ids_B),
    
    projected_alignment(rhs.projected_alignment),
    weights(rhs.weights),
    
    LA(rhs.LA),
    LB(rhs.LB),
    
    folding_A_mask_sparse(rhs.folding_A_mask_sparse),
    folding_B_mask_sparse(rhs.folding_B_mask_sparse),
    alignment_mask_sparse(rhs.alignment_mask_sparse),
    
    folding_A_matrix(rhs.folding_A_matrix),
    folding_B_matrix(rhs.folding_B_matrix),
    alignment_matrix(rhs.alignment_matrix),

    folding_A_offset(rhs.folding_A_offset),
    folding_B_offset(rhs.folding_B_offset),
    alignment_offset(rhs.alignment_offset)
{}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeFoldingMatrix()
//
// Compute matrices associated with structure scoring.
/////////////////////////////////////////////////////////////////

void ScoringScheme::ComputeFoldingMatrix(SparseMatrix<Real> &folding_mask_sparse,
                                         std::vector<Real> &folding_matrix,
                                         Real &folding_offset,
                                         const std::vector<int> &ids)
{
    const int L = int(projected_alignment[ids[0]].data.length()) - 1;
    std::vector<Real> unsparse_mask((L+1)*(L+1));
    
    folding_matrix.clear();
    folding_matrix.resize((L+1)*(L+1));
    folding_offset = Real(0);
    
    // consider each sequence in "ids"
    
    for (size_t k = 0; k < ids.size(); k++)
    {
        const std::vector<int> positions = GetSequencePositions(projected_alignment[ids[k]].data);    
        
        // expand posteriors

        std::vector<Real> P_paired;
        std::vector<Real> P_unpaired;
        std::vector<Real> seq_unsparse_mask;
        
        ExpandFoldingPosteriors(ids[k], positions, P_paired, P_unpaired, seq_unsparse_mask);

        unsparse_mask += seq_unsparse_mask;

        // update scoring tables

        for (int i = 0; i < Weight_NUM_WEIGHTS; i++)
        {
            std::pair<std::vector<Real>, Real> feature = ComputeFoldingFeature(i, P_paired, P_unpaired);
            if (feature.first.size() == 0) continue;

            folding_matrix += weights[i] * ExpandMatrix<Real>(feature.first, L+1, L+1, positions, positions);
            folding_offset += weights[i] * feature.second;
        }

    }

    folding_mask_sparse = SparseMatrix<Real>(&unsparse_mask[0], L+1, L+1, Real(0));
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeAlignmentMatrix()
//
// Compute matrices associated with alignment scoring.
/////////////////////////////////////////////////////////////////

void ScoringScheme::ComputeAlignmentMatrix(SparseMatrix<Real> &alignment_mask_sparse,
                                           std::vector<Real> &alignment_matrix,
                                           Real &alignment_offset,
                                           const std::vector<int> &ids_A,
                                           const std::vector<int> &ids_B)
{
    const int LA = int(projected_alignment[ids_A[0]].data.length()) - 1;
    const int LB = int(projected_alignment[ids_B[0]].data.length()) - 1;
    std::vector<Real> unsparse_mask((LA+1)*(LB+1));
    
    alignment_matrix.clear();
    alignment_matrix.resize((LA+1)*(LB+1));
    alignment_offset = Real(0);

    // consider each pair of sequences in "ids_A" and "ids_B"
    
    for (size_t kA = 0; kA < ids_A.size(); kA++)
    {
        const std::vector<int> positions_A = GetSequencePositions(projected_alignment[ids_A[kA]].data);
        for (size_t kB = 0; kB < ids_B.size(); kB++)
        {
            const std::vector<int> positions_B = GetSequencePositions(projected_alignment[ids_B[kB]].data);

            // expand posteriors

            std::vector<Real> P_aligned;
            std::vector<Real> P_unaligned_A;
            std::vector<Real> P_unaligned_B;
            std::vector<Real> pair_unsparse_mask;
            
            ExpandAlignmentPosteriors(ids_A[kA], ids_B[kB], positions_A, positions_B, P_aligned, P_unaligned_A, P_unaligned_B, pair_unsparse_mask);

            unsparse_mask += pair_unsparse_mask;
            
            // update scoring tables
            
            for (int i = 0; i < Weight_NUM_WEIGHTS; i++)
            {
                std::pair<std::vector<Real>, Real> feature = ComputeAlignmentFeature(i, ids_A[kA], ids_B[kB], P_aligned, P_unaligned_A, P_unaligned_B);
                if (feature.first.size() == 0) continue;
                
                alignment_matrix += weights[i] * ExpandMatrix<Real>(feature.first, LA+1, LB+1, positions_A, positions_B);
                alignment_offset += weights[i] * feature.second;
            }

        }
    }

    alignment_mask_sparse = SparseMatrix<Real>(&unsparse_mask[0], LA+1, LB+1, Real(0));
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ExpandFoldingPosteriors()
//
// Expand folding posterior probability matrix.
/////////////////////////////////////////////////////////////////

void ScoringScheme::ExpandFoldingPosteriors(const int id,
                                            const std::vector<int> &positions,
                                            std::vector<Real> &P_paired,
                                            std::vector<Real> &P_unpaired,
                                            std::vector<Real> &unsparse_mask) const
{
    const int L = int(projected_alignment[id].data.length()) - 1;
    const SparseMatrix<Real> &sparse = base_pairing_posteriors[id];
    const int PL = sparse.GetNumRows()-1;
    Assert(PL+1 == int(positions.size()), "Dimension mismatch.");

    P_paired = sparse.GetUnsparse();
    P_unpaired = std::vector<Real>(PL+1, Real(1));

    for (int i = 1; i <= PL; i++)
    {
        for (int j = i+1; j <= PL; j++)
        {
            P_unpaired[i] -= P_paired[i*(PL+1)+j];
            P_unpaired[j] -= P_paired[i*(PL+1)+j];
        }
    }
    
    P_paired = Clip(P_paired, Real(0), Real(1));
    P_unpaired = Clip(P_unpaired, Real(0), Real(1));

    unsparse_mask = ExpandMatrix<Real>(sparse.GetUnsparseMask(), L+1, L+1, positions, positions);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ExpandAlignmentPosteriors()
//
// Expand alignment posterior probability matrix.
/////////////////////////////////////////////////////////////////

void ScoringScheme::ExpandAlignmentPosteriors(const int id_A,
                                              const int id_B,
                                              const std::vector<int> &positions_A,
                                              const std::vector<int> &positions_B,
                                              std::vector<Real> &P_aligned,
                                              std::vector<Real> &P_unaligned_A,
                                              std::vector<Real> &P_unaligned_B,
                                              std::vector<Real> &unsparse_mask) const
{
    const int LA = int(projected_alignment[id_A].data.length()) - 1;
    const int LB = int(projected_alignment[id_B].data.length()) - 1;
    const int K = int(projected_alignment.size());
    SparseMatrix<Real> sparse = alignment_posteriors[std::min(id_A, id_B) * K + std::max(id_A, id_B)];
    if (id_A > id_B) sparse = SparseMatrix<Real>(sparse, SparseMatrix<Real>::TRANSPOSE);
    const int PLA = sparse.GetNumRows()-1;
    const int PLB = sparse.GetNumCols()-1;
    Assert(PLA+1 == int(positions_A.size()), "Dimension mismatch.");
    Assert(PLB+1 == int(positions_B.size()), "Dimension mismatch.");

    P_aligned = sparse.GetUnsparse();
    P_unaligned_A = std::vector<Real>(PLA+1, Real(1));
    P_unaligned_B = std::vector<Real>(PLB+1, Real(1));

    for (int i = 1; i <= PLA; i++)
    {
        for (int j = 1; j <= PLB; j++)
        {
            P_unaligned_A[i] -= P_aligned[i*(PLB+1)+j];
            P_unaligned_B[j] -= P_aligned[i*(PLB+1)+j];
        }
    }
    
    P_aligned = Clip(P_aligned, Real(0), Real(1));
    P_unaligned_A = Clip(P_unaligned_A, Real(0), Real(1));
    P_unaligned_B = Clip(P_unaligned_B, Real(0), Real(1));

    unsparse_mask = ExpandMatrix<Real>(sparse.GetUnsparseMask(), LA+1, LB+1, positions_A, positions_B);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeFoldingFeature()
//
// Compute a particular folding feature matrix.
/////////////////////////////////////////////////////////////////

std::pair<std::vector<Real>, Real> ScoringScheme::ComputeFoldingFeature(const int index,
                                                                        const std::vector<Real> &P_paired,
                                                                        const std::vector<Real> &P_unpaired) const
{
    const int PL = int(P_unpaired.size())-1;
    std::vector<Real> paired((PL+1)*(PL+1));
    std::vector<Real> unpaired(PL+1);

    
    switch (index)
    {
#if PARAMS_PAIRED_LOCARNA_LOG_ODDS
        case Weight_PAIRED_LOCARNA_LOG_ODDS:            paired = Log(P_paired + Real(1e-10)) / Log(Real(2 * PL)) + Real(1); break;
#endif
#if PARAMS_PAIRED_CONSTANT
        case Weight_PAIRED_CONSTANT:                    paired = std::vector<Real>(P_paired.size(), Real(1)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY
        case Weight_PAIRED_PROBABILITY:                 paired = P_paired; break;
#endif
#if PARAMS_PAIRED_PROBABILITY_SQUARED
        case Weight_PAIRED_PROBABILITY_SQUARED:         paired = P_paired * P_paired; break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG
        case Weight_PAIRED_PROBABILITY_LOG:             paired = Log(Max(Real(1e-10), P_paired)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_000001
        case Weight_PAIRED_PROBABILITY_LOG_0_000001:    paired = Log(Real(1e-6) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_00001
        case Weight_PAIRED_PROBABILITY_LOG_0_00001:     paired = Log(Real(1e-5) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_0001
        case Weight_PAIRED_PROBABILITY_LOG_0_0001:      paired = Log(Real(1e-4) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_001
        case Weight_PAIRED_PROBABILITY_LOG_0_001:       paired = Log(Real(1e-3) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_01
        case Weight_PAIRED_PROBABILITY_LOG_0_01:        paired = Log(Real(1e-2) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_0_1
        case Weight_PAIRED_PROBABILITY_LOG_0_1:         paired = Log(Real(1e-1) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_1
        case Weight_PAIRED_PROBABILITY_LOG_1:           paired = Log(Real(1) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_LOG_10
        case Weight_PAIRED_PROBABILITY_LOG_10:          paired = Log(Real(1e1) + P_paired); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_00_05
        case Weight_PAIRED_PROBABILITY_00_05:           paired = Ramp(P_paired, Real(0.00), Real(0.05)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_05_10
        case Weight_PAIRED_PROBABILITY_05_10:           paired = Ramp(P_paired, Real(0.05), Real(0.10)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_10_15
        case Weight_PAIRED_PROBABILITY_10_15:           paired = Ramp(P_paired, Real(0.10), Real(0.15)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_15_20
        case Weight_PAIRED_PROBABILITY_15_20:           paired = Ramp(P_paired, Real(0.15), Real(0.20)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_20_25
        case Weight_PAIRED_PROBABILITY_20_25:           paired = Ramp(P_paired, Real(0.20), Real(0.25)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_25_30
        case Weight_PAIRED_PROBABILITY_25_30:           paired = Ramp(P_paired, Real(0.25), Real(0.30)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_30_35
        case Weight_PAIRED_PROBABILITY_30_35:           paired = Ramp(P_paired, Real(0.30), Real(0.35)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_35_40
        case Weight_PAIRED_PROBABILITY_35_40:           paired = Ramp(P_paired, Real(0.35), Real(0.40)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_40_45
        case Weight_PAIRED_PROBABILITY_40_45:           paired = Ramp(P_paired, Real(0.40), Real(0.45)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_45_50
        case Weight_PAIRED_PROBABILITY_45_50:           paired = Ramp(P_paired, Real(0.45), Real(0.50)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_50_55
        case Weight_PAIRED_PROBABILITY_50_55:           paired = Ramp(P_paired, Real(0.50), Real(0.55)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_55_60
        case Weight_PAIRED_PROBABILITY_55_60:           paired = Ramp(P_paired, Real(0.55), Real(0.60)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_60_65
        case Weight_PAIRED_PROBABILITY_60_65:           paired = Ramp(P_paired, Real(0.60), Real(0.65)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_65_70
        case Weight_PAIRED_PROBABILITY_65_70:           paired = Ramp(P_paired, Real(0.65), Real(0.70)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_70_75
        case Weight_PAIRED_PROBABILITY_70_75:           paired = Ramp(P_paired, Real(0.70), Real(0.75)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_75_80
        case Weight_PAIRED_PROBABILITY_75_80:           paired = Ramp(P_paired, Real(0.75), Real(0.80)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_80_85
        case Weight_PAIRED_PROBABILITY_80_85:           paired = Ramp(P_paired, Real(0.80), Real(0.85)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_85_90
        case Weight_PAIRED_PROBABILITY_85_90:           paired = Ramp(P_paired, Real(0.85), Real(0.90)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_90_95
        case Weight_PAIRED_PROBABILITY_90_95:           paired = Ramp(P_paired, Real(0.90), Real(0.95)); break;
#endif
#if PARAMS_PAIRED_PROBABILITY_95_00
        case Weight_PAIRED_PROBABILITY_95_00:           paired = Ramp(P_paired, Real(0.95), Real(1.00)); break;
#endif
#if PARAMS_UNPAIRED_CONSTANT
        case Weight_UNPAIRED_CONSTANT:                  unpaired = std::vector<Real>(P_unpaired.size(), Real(1)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY
        case Weight_UNPAIRED_PROBABILITY:               unpaired = P_unpaired; break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_SQUARED
        case Weight_UNPAIRED_PROBABILITY_SQUARED:       unpaired = P_unpaired * P_unpaired; break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG
        case Weight_UNPAIRED_PROBABILITY_LOG:           unpaired = Log(Max(Real(1e-10), P_unpaired)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_000001
        case Weight_UNPAIRED_PROBABILITY_LOG_0_000001:  unpaired = Log(Real(1e-6) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_00001
        case Weight_UNPAIRED_PROBABILITY_LOG_0_00001:   unpaired = Log(Real(1e-5) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_0001
        case Weight_UNPAIRED_PROBABILITY_LOG_0_0001:    unpaired = Log(Real(1e-4) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_001
        case Weight_UNPAIRED_PROBABILITY_LOG_0_001:     unpaired = Log(Real(1e-3) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_01
        case Weight_UNPAIRED_PROBABILITY_LOG_0_01:      unpaired = Log(Real(1e-2) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_0_1
        case Weight_UNPAIRED_PROBABILITY_LOG_0_1:       unpaired = Log(Real(1e-1) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_1
        case Weight_UNPAIRED_PROBABILITY_LOG_1:         unpaired = Log(Real(1) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_LOG_10
        case Weight_UNPAIRED_PROBABILITY_LOG_10:        unpaired = Log(Real(1e1) + P_unpaired); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_00_05
        case Weight_UNPAIRED_PROBABILITY_00_05:         unpaired = Ramp(P_unpaired, Real(0.00), Real(0.05)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_05_10
        case Weight_UNPAIRED_PROBABILITY_05_10:         unpaired = Ramp(P_unpaired, Real(0.05), Real(0.10)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_10_15
        case Weight_UNPAIRED_PROBABILITY_10_15:         unpaired = Ramp(P_unpaired, Real(0.10), Real(0.15)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_15_20
        case Weight_UNPAIRED_PROBABILITY_15_20:         unpaired = Ramp(P_unpaired, Real(0.15), Real(0.20)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_20_25
        case Weight_UNPAIRED_PROBABILITY_20_25:         unpaired = Ramp(P_unpaired, Real(0.20), Real(0.25)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_25_30
        case Weight_UNPAIRED_PROBABILITY_25_30:         unpaired = Ramp(P_unpaired, Real(0.25), Real(0.30)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_30_35
        case Weight_UNPAIRED_PROBABILITY_30_35:         unpaired = Ramp(P_unpaired, Real(0.30), Real(0.35)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_35_40
        case Weight_UNPAIRED_PROBABILITY_35_40:         unpaired = Ramp(P_unpaired, Real(0.35), Real(0.40)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_40_45
        case Weight_UNPAIRED_PROBABILITY_40_45:         unpaired = Ramp(P_unpaired, Real(0.40), Real(0.45)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_45_50
        case Weight_UNPAIRED_PROBABILITY_45_50:         unpaired = Ramp(P_unpaired, Real(0.45), Real(0.50)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_50_55
        case Weight_UNPAIRED_PROBABILITY_50_55:         unpaired = Ramp(P_unpaired, Real(0.50), Real(0.55)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_55_60
        case Weight_UNPAIRED_PROBABILITY_55_60:         unpaired = Ramp(P_unpaired, Real(0.55), Real(0.60)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_60_65
        case Weight_UNPAIRED_PROBABILITY_60_65:         unpaired = Ramp(P_unpaired, Real(0.60), Real(0.65)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_65_70
        case Weight_UNPAIRED_PROBABILITY_65_70:         unpaired = Ramp(P_unpaired, Real(0.65), Real(0.70)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_70_75
        case Weight_UNPAIRED_PROBABILITY_70_75:         unpaired = Ramp(P_unpaired, Real(0.70), Real(0.75)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_75_80
        case Weight_UNPAIRED_PROBABILITY_75_80:         unpaired = Ramp(P_unpaired, Real(0.75), Real(0.80)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_80_85
        case Weight_UNPAIRED_PROBABILITY_80_85:         unpaired = Ramp(P_unpaired, Real(0.80), Real(0.85)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_85_90
        case Weight_UNPAIRED_PROBABILITY_85_90:         unpaired = Ramp(P_unpaired, Real(0.85), Real(0.90)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_90_95
        case Weight_UNPAIRED_PROBABILITY_90_95:         unpaired = Ramp(P_unpaired, Real(0.90), Real(0.95)); break;
#endif
#if PARAMS_UNPAIRED_PROBABILITY_95_00
        case Weight_UNPAIRED_PROBABILITY_95_00:         unpaired = Ramp(P_unpaired, Real(0.95), Real(1.00)); break;
#endif
        default:                                        break;
    }

    // convert to final scoring table

    std::vector<Real> ret((PL+1)*(PL+1));
    Real offset = 0;

    for (int i = 1; i <= PL; i++)
    {
        offset += unpaired[i];
        for (int j = 1; j <= PL; j++)
        {
            ret[i*(PL+1)+j] = paired[i*(PL+1)+j] - unpaired[i] - unpaired[j];
        }
    }

    return std::make_pair(ret, offset);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeAlignmentFeature()
//
// Compute a particular alignment feature matrix.
/////////////////////////////////////////////////////////////////

std::pair<std::vector<Real>, Real> ScoringScheme::ComputeAlignmentFeature(const int index,
                                                                          const int id_A,
                                                                          const int id_B,
                                                                          const std::vector<Real> &P_aligned,
                                                                          const std::vector<Real> &P_unaligned_A,
                                                                          const std::vector<Real> &P_unaligned_B) const
{
    const int PLA = int(P_unaligned_A.size())-1;
    const int PLB = int(P_unaligned_B.size())-1;
    std::vector<Real> aligned((PLA+1)*(PLB+1));
    std::vector<Real> unaligned_A(PLA+1);
    std::vector<Real> unaligned_B(PLB+1);
    
    switch (index)
    {
        
#if PARAMS_ALIGNED_DOT_PLOT
        case Weight_ALIGNED_DOT_PLOT:
        {
            const std::string &x = RemoveGaps(projected_alignment[id_A].data);
            const std::string &y = RemoveGaps(projected_alignment[id_B].data);
            for (size_t i = 1; i < x.length(); i++)
                for (size_t j = 1; j < y.length(); j++)
                    aligned[i*(PLB+1)+j] = Ind(isalpha(x[i]) && x[i] == y[j]);
        }
        break;
#endif
        
#if PARAMS_ALIGNED_CONSTANT
        case Weight_ALIGNED_CONSTANT:                   aligned = std::vector<Real>(P_aligned.size(), Real(1)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY
        case Weight_ALIGNED_PROBABILITY:                aligned = P_aligned; break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_SQUARED
        case Weight_ALIGNED_PROBABILITY_SQUARED:        aligned = P_aligned * P_aligned; break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG
        case Weight_ALIGNED_PROBABILITY_LOG:            aligned = Log(Max(Real(1e-10), P_aligned)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_000001
        case Weight_ALIGNED_PROBABILITY_LOG_0_000001:   aligned = Log(Real(1e-6) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_00001
        case Weight_ALIGNED_PROBABILITY_LOG_0_00001:    aligned = Log(Real(1e-5) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_0001
        case Weight_ALIGNED_PROBABILITY_LOG_0_0001:     aligned = Log(Real(1e-4) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_001
        case Weight_ALIGNED_PROBABILITY_LOG_0_001:      aligned = Log(Real(1e-3) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_01
        case Weight_ALIGNED_PROBABILITY_LOG_0_01:       aligned = Log(Real(1e-2) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_0_1
        case Weight_ALIGNED_PROBABILITY_LOG_0_1:        aligned = Log(Real(1e-1) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_1
        case Weight_ALIGNED_PROBABILITY_LOG_1:          aligned = Log(Real(1) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_LOG_10
        case Weight_ALIGNED_PROBABILITY_LOG_10:         aligned = Log(Real(1e1) + P_aligned); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_00_05
        case Weight_ALIGNED_PROBABILITY_00_05:          aligned = Ramp(P_aligned, Real(0.00), Real(0.05)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_05_10
        case Weight_ALIGNED_PROBABILITY_05_10:          aligned = Ramp(P_aligned, Real(0.05), Real(0.10)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_10_15
        case Weight_ALIGNED_PROBABILITY_10_15:          aligned = Ramp(P_aligned, Real(0.10), Real(0.15)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_15_20
        case Weight_ALIGNED_PROBABILITY_15_20:          aligned = Ramp(P_aligned, Real(0.15), Real(0.20)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_20_25
        case Weight_ALIGNED_PROBABILITY_20_25:          aligned = Ramp(P_aligned, Real(0.20), Real(0.25)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_25_30
        case Weight_ALIGNED_PROBABILITY_25_30:          aligned = Ramp(P_aligned, Real(0.25), Real(0.30)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_30_35
        case Weight_ALIGNED_PROBABILITY_30_35:          aligned = Ramp(P_aligned, Real(0.30), Real(0.35)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_35_40
        case Weight_ALIGNED_PROBABILITY_35_40:          aligned = Ramp(P_aligned, Real(0.35), Real(0.40)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_40_45
        case Weight_ALIGNED_PROBABILITY_40_45:          aligned = Ramp(P_aligned, Real(0.40), Real(0.45)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_45_50
        case Weight_ALIGNED_PROBABILITY_45_50:          aligned = Ramp(P_aligned, Real(0.45), Real(0.50)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_50_55
        case Weight_ALIGNED_PROBABILITY_50_55:          aligned = Ramp(P_aligned, Real(0.50), Real(0.55)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_55_60
        case Weight_ALIGNED_PROBABILITY_55_60:          aligned = Ramp(P_aligned, Real(0.55), Real(0.60)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_60_65
        case Weight_ALIGNED_PROBABILITY_60_65:          aligned = Ramp(P_aligned, Real(0.60), Real(0.65)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_65_70
        case Weight_ALIGNED_PROBABILITY_65_70:          aligned = Ramp(P_aligned, Real(0.65), Real(0.70)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_70_75
        case Weight_ALIGNED_PROBABILITY_70_75:          aligned = Ramp(P_aligned, Real(0.70), Real(0.75)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_75_80
        case Weight_ALIGNED_PROBABILITY_75_80:          aligned = Ramp(P_aligned, Real(0.75), Real(0.80)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_80_85
        case Weight_ALIGNED_PROBABILITY_80_85:          aligned = Ramp(P_aligned, Real(0.80), Real(0.85)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_85_90
        case Weight_ALIGNED_PROBABILITY_85_90:          aligned = Ramp(P_aligned, Real(0.85), Real(0.90)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_90_95
        case Weight_ALIGNED_PROBABILITY_90_95:          aligned = Ramp(P_aligned, Real(0.90), Real(0.95)); break;
#endif
#if PARAMS_ALIGNED_PROBABILITY_95_00
        case Weight_ALIGNED_PROBABILITY_95_00:          aligned = Ramp(P_aligned, Real(0.95), Real(1.00)); break;
#endif

#if PARAMS_UNALIGNED_CONSTANT
        case Weight_UNALIGNED_CONSTANT:
            unaligned_A = std::vector<Real>(P_unaligned_A.size(), Real(1));
            unaligned_B = std::vector<Real>(P_unaligned_B.size(), Real(1));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY
        case Weight_UNALIGNED_PROBABILITY:
            unaligned_A = P_unaligned_A;
            unaligned_B = P_unaligned_B;
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_SQUARED
        case Weight_UNALIGNED_PROBABILITY_SQUARED:
            unaligned_A = P_unaligned_A * P_unaligned_A;
            unaligned_B = P_unaligned_B * P_unaligned_B;
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG
        case Weight_UNALIGNED_PROBABILITY_LOG:
            unaligned_A = Log(Max(Real(1e-10), P_unaligned_A));
            unaligned_B = Log(Max(Real(1e-10), P_unaligned_B));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_000001
        case Weight_UNALIGNED_PROBABILITY_LOG_0_000001:
            unaligned_A = Log(Real(1e-6) + P_unaligned_A);
            unaligned_B = Log(Real(1e-6) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_00001
        case Weight_UNALIGNED_PROBABILITY_LOG_0_00001:
            unaligned_A = Log(Real(1e-5) + P_unaligned_A);
            unaligned_B = Log(Real(1e-5) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_0001
        case Weight_UNALIGNED_PROBABILITY_LOG_0_0001:
            unaligned_A = Log(Real(1e-4) + P_unaligned_A);
            unaligned_B = Log(Real(1e-4) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_001
        case Weight_UNALIGNED_PROBABILITY_LOG_0_001:
            unaligned_A = Log(Real(1e-3) + P_unaligned_A);
            unaligned_B = Log(Real(1e-3) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_01
        case Weight_UNALIGNED_PROBABILITY_LOG_0_01:
            unaligned_A = Log(Real(1e-2) + P_unaligned_A);
            unaligned_B = Log(Real(1e-2) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_0_1
        case Weight_UNALIGNED_PROBABILITY_LOG_0_1:
            unaligned_A = Log(Real(1e-1) + P_unaligned_A);
            unaligned_B = Log(Real(1e-1) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_1
        case Weight_UNALIGNED_PROBABILITY_LOG_1:
            unaligned_A = Log(Real(1) + P_unaligned_A);
            unaligned_B = Log(Real(1) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_LOG_10
        case Weight_UNALIGNED_PROBABILITY_LOG_10:
            unaligned_A = Log(Real(1e1) + P_unaligned_A);
            unaligned_B = Log(Real(1e1) + P_unaligned_B);
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_00_05
        case Weight_UNALIGNED_PROBABILITY_00_05:
            unaligned_A = Ramp(P_unaligned_A, Real(0.00), Real(0.05));
            unaligned_B = Ramp(P_unaligned_B, Real(0.00), Real(0.05));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_05_10
        case Weight_UNALIGNED_PROBABILITY_05_10:
            unaligned_A = Ramp(P_unaligned_A, Real(0.05), Real(0.10));
            unaligned_B = Ramp(P_unaligned_B, Real(0.05), Real(0.10));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_10_15
        case Weight_UNALIGNED_PROBABILITY_10_15:
            unaligned_A = Ramp(P_unaligned_A, Real(0.10), Real(0.15));
            unaligned_B = Ramp(P_unaligned_B, Real(0.10), Real(0.15));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_15_20
        case Weight_UNALIGNED_PROBABILITY_15_20:
            unaligned_A = Ramp(P_unaligned_A, Real(0.15), Real(0.20));
            unaligned_B = Ramp(P_unaligned_B, Real(0.15), Real(0.20));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_20_25
        case Weight_UNALIGNED_PROBABILITY_20_25:
            unaligned_A = Ramp(P_unaligned_A, Real(0.20), Real(0.25));
            unaligned_B = Ramp(P_unaligned_B, Real(0.20), Real(0.25));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_25_30
        case Weight_UNALIGNED_PROBABILITY_25_30:
            unaligned_A = Ramp(P_unaligned_A, Real(0.25), Real(0.30));
            unaligned_B = Ramp(P_unaligned_B, Real(0.25), Real(0.30));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_30_35
        case Weight_UNALIGNED_PROBABILITY_30_35:
            unaligned_A = Ramp(P_unaligned_A, Real(0.30), Real(0.35));
            unaligned_B = Ramp(P_unaligned_B, Real(0.30), Real(0.35));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_35_40
        case Weight_UNALIGNED_PROBABILITY_35_40:
            unaligned_A = Ramp(P_unaligned_A, Real(0.35), Real(0.40));
            unaligned_B = Ramp(P_unaligned_B, Real(0.35), Real(0.40));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_40_45
        case Weight_UNALIGNED_PROBABILITY_40_45:
            unaligned_A = Ramp(P_unaligned_A, Real(0.40), Real(0.45));
            unaligned_B = Ramp(P_unaligned_B, Real(0.40), Real(0.45));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_45_50
        case Weight_UNALIGNED_PROBABILITY_45_50:
            unaligned_A = Ramp(P_unaligned_A, Real(0.45), Real(0.50));
            unaligned_B = Ramp(P_unaligned_B, Real(0.45), Real(0.50));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_50_55
        case Weight_UNALIGNED_PROBABILITY_50_55:
            unaligned_A = Ramp(P_unaligned_A, Real(0.50), Real(0.55));
            unaligned_B = Ramp(P_unaligned_B, Real(0.50), Real(0.55));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_55_60
        case Weight_UNALIGNED_PROBABILITY_55_60:
            unaligned_A = Ramp(P_unaligned_A, Real(0.55), Real(0.60));
            unaligned_B = Ramp(P_unaligned_B, Real(0.55), Real(0.60));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_60_65
        case Weight_UNALIGNED_PROBABILITY_60_65:
            unaligned_A = Ramp(P_unaligned_A, Real(0.60), Real(0.65));
            unaligned_B = Ramp(P_unaligned_B, Real(0.60), Real(0.65));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_65_70
        case Weight_UNALIGNED_PROBABILITY_65_70:
            unaligned_A = Ramp(P_unaligned_A, Real(0.65), Real(0.70));
            unaligned_B = Ramp(P_unaligned_B, Real(0.65), Real(0.70));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_70_75
        case Weight_UNALIGNED_PROBABILITY_70_75:
            unaligned_A = Ramp(P_unaligned_A, Real(0.70), Real(0.75));
            unaligned_B = Ramp(P_unaligned_B, Real(0.70), Real(0.75));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_75_80
        case Weight_UNALIGNED_PROBABILITY_75_80:
            unaligned_A = Ramp(P_unaligned_A, Real(0.75), Real(0.80));
            unaligned_B = Ramp(P_unaligned_B, Real(0.75), Real(0.80));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_80_85
        case Weight_UNALIGNED_PROBABILITY_80_85:
            unaligned_A = Ramp(P_unaligned_A, Real(0.80), Real(0.85));
            unaligned_B = Ramp(P_unaligned_B, Real(0.80), Real(0.85));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_85_90
        case Weight_UNALIGNED_PROBABILITY_85_90:
            unaligned_A = Ramp(P_unaligned_A, Real(0.85), Real(0.90));
            unaligned_B = Ramp(P_unaligned_B, Real(0.85), Real(0.90));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_90_95
        case Weight_UNALIGNED_PROBABILITY_90_95:
            unaligned_A = Ramp(P_unaligned_A, Real(0.90), Real(0.95));
            unaligned_B = Ramp(P_unaligned_B, Real(0.90), Real(0.95));
            break;
#endif
#if PARAMS_UNALIGNED_PROBABILITY_95_00
        case Weight_UNALIGNED_PROBABILITY_95_00:
            unaligned_A = Ramp(P_unaligned_A, Real(0.95), Real(1.00));
            unaligned_B = Ramp(P_unaligned_B, Real(0.95), Real(1.00));
            break;
#endif

        default:
            break;
    }

    // convert to final scoring table

    std::vector<Real> ret((PLA+1)*(PLB+1));
    Real offset = 0;

    for (int i = 1; i <= PLA; i++) offset += unaligned_A[i];
    for (int j = 1; j <= PLB; j++) offset += unaligned_B[j];
    
    for (int i = 1; i <= PLA; i++)
    {
        for (int j = 1; j <= PLB; j++)
        {
            ret[i*(PLB+1)+j] = aligned[i*(PLB+1)+j] - unaligned_A[i] - unaligned_B[j];
        }
    }

    return std::make_pair(ret, offset);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::AugmentWithLoss()
//
// Add Hamming loss for folding and alignment scores to the
// scoring tables.
/////////////////////////////////////////////////////////////////

void ScoringScheme::AugmentWithLoss(const std::vector<int> &true_mapping_A,
                                    const std::vector<int> &true_mapping_B,
                                    const std::vector<int> &true_aligned_to_A,
                                    const std::vector<int> &true_aligned_to_B,
                                    const Real folding_FP_margin,
                                    const Real folding_FN_margin,
                                    const Real alignment_FP_margin,
                                    const Real alignment_FN_margin)
{
    AugmentWithFoldingLoss(folding_A_matrix, folding_A_offset, true_mapping_A, ids_A, folding_FP_margin, folding_FN_margin);
    AugmentWithFoldingLoss(folding_B_matrix, folding_B_offset, true_mapping_B, ids_B, folding_FP_margin, folding_FN_margin);
    AugmentWithAlignmentLoss(alignment_matrix, alignment_offset, true_aligned_to_A, true_aligned_to_B, alignment_FP_margin, alignment_FN_margin);
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::AugmentWithFoldingLoss()
//
// Add Hamming distance loss for folding scores.  The key idea is that
// with loss-augmented inference, instead of searching for the maximum
// scoring parse, we boost each parse's score based on the loss.  We
// can accomplish this implicitly by modifying the scoring matrices.
//
// In the case of folding scores,
//     folding_FP = (i,j) not in TRUE, (i,j) in PREDICTED
//     folding_FN = (i,j) in TRUE, (i,j) not in PREDICTED
//
// The "offset" entry in the scoring table is supposed to
// capture the score that we attain if no base-pairings are
// predicted.  The "scores" matrix in the scoring table gives
// the changes in this score for each base-pairing (i,j) that
// we accept.  Under this, it follows that
// 
// * "offset" should include folding_FN_margin for each (i,j) in
//   TRUE.
//
// * "scores(i,j)" should include
//    * folding_FP_margin for each (i,j) not in TRUE
//    * -(folding_FN_margin) for each (i,j) in TRUE
/////////////////////////////////////////////////////////////////

void ScoringScheme::AugmentWithFoldingLoss(std::vector<Real> &folding_matrix,
                                           Real &folding_offset,
                                           const std::vector<int> &mapping,
                                           const std::vector<int> &ids,
                                           const Real folding_FP_margin,
                                           const Real folding_FN_margin)
{
    const int L = int(projected_alignment[ids[0]].data.length()) - 1;
    
    // update offset
    
    for (int i = 1; i <= L; i++)
        folding_offset += folding_FN_margin * Ind(mapping[i] != 0);

    // update score matrix
    
    for (int i = 1; i <= L; i++)
    {
        for (int j = i+1; j <= L; j++)
        {
            folding_matrix[i*(L+1)+j] += 
                folding_FP_margin * Real(Ind(mapping[i] != j) + Ind(mapping[j] != i))
                - folding_FN_margin * Real(Ind(mapping[i] == j) + Ind(mapping[j] == i));
        }
    }
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::AugmentWithAlignmentLoss()
//
// Add Hamming distance loss for alignment scores.  The key idea is
// that with loss-augmented inference, instead of searching for the
// maximum scoring parse, we boost each parse's score based on the
// loss.  We can accomplish this implicitly by modifying the scoring
// matrices.
//
// In the case of alignment scores,
//     alignment_FP = (i,j) not in TRUE, (i,j) in PREDICTED
//     alignment_FN = (i,j) in TRUE, (i,j) not in PREDICTED
//
// The "offset" entry in the scoring table is supposed to
// capture the score that we attain if no aligned pairs are
// predicted.  The "scores" matrix in the scoring table gives
// the changes in this score for each aligned pair (i,j) that
// we accept.  Under this, it follows that
// 
// * "offset" should include alignment_FN_margin for each (i,j) in
//   TRUE.
//
// * "scores(i,j)" should include
//    * alignment_FP_margin for each (i,j) not in TRUE
//    * -(alignment_FN_margin) for each (i,j) in TRUE
/////////////////////////////////////////////////////////////////

void ScoringScheme::AugmentWithAlignmentLoss(std::vector<Real> &alignment_matrix,
                                             Real &alignment_offset,
                                             const std::vector<int> &aligned_to_A,
                                             const std::vector<int> &aligned_to_B,
                                             const Real alignment_FP_margin,
                                             const Real alignment_FN_margin)
{
    const int LA = int(projected_alignment[ids_A[0]].data.length()) - 1;
    const int LB = int(projected_alignment[ids_B[0]].data.length()) - 1;
    
    // update offset
    
    for (int i = 1; i <= LA; i++)
        alignment_offset += alignment_FN_margin * Ind(aligned_to_A[i] != 0);
    
    for (int i = 1; i <= LB; i++)
        alignment_offset += alignment_FN_margin * Ind(aligned_to_B[i] != 0);
    
    // update score matrix

    for (int i = 1; i <= LA; i++)
    {
        for (int j = 1; j <= LB; j++)
        {
            alignment_matrix[i*(LB+1)+j] += 
                alignment_FP_margin * Real(Ind(aligned_to_A[i] != j) + Ind(aligned_to_B[j] != i))
                - alignment_FN_margin * Real(Ind(aligned_to_A[i] == j) + Ind(aligned_to_B[j] == i));
        }
    }
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeFeatures()
//
// Compute features for the alignment of two groups of sequences.
/////////////////////////////////////////////////////////////////

std::vector<Real> ScoringScheme::ComputeFeatures(const std::vector<int> &mapping_A,
                                                 const std::vector<int> &mapping_B,
                                                 const std::vector<int> &aligned_to_A,
                                                 const std::vector<int> &aligned_to_B) const
{
    std::vector<Real> parse_folding_A_features = ComputeParseFoldingFeatures(mapping_A, ids_A);
    std::vector<Real> parse_folding_B_features = ComputeParseFoldingFeatures(mapping_B, ids_B);
    std::vector<Real> parse_alignment_features = ComputeParseAlignmentFeatures(aligned_to_A, aligned_to_B);
    
    return parse_folding_A_features + parse_folding_B_features + parse_alignment_features;
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeParseFoldingFeatures()
//
// Compute folding features for a single parse.
/////////////////////////////////////////////////////////////////

std::vector<Real> ScoringScheme::ComputeParseFoldingFeatures(const std::vector<int> &mapping,
                                                             const std::vector<int> &ids) const
{
    const int L = int(projected_alignment[ids[0]].data.length()) - 1;
    std::vector<Real> unsparse_mask((L+1)*(L+1));
    std::vector<Real> features(Weight_NUM_WEIGHTS);
    
    // consider each sequence in "ids"
    
    for (size_t k = 0; k < ids.size(); k++)
    {
        const std::vector<int> positions = GetSequencePositions(projected_alignment[ids[k]].data);    
        
        // expand posteriors

        std::vector<Real> P_paired;
        std::vector<Real> P_unpaired;
        std::vector<Real> seq_unsparse_mask;

        ExpandFoldingPosteriors(ids[k], positions, P_paired, P_unpaired, seq_unsparse_mask);

        unsparse_mask += seq_unsparse_mask;

        // update features

        for (int i = 0; i < Weight_NUM_WEIGHTS; i++)
        {
            std::pair<std::vector<Real>, Real> feature = ComputeFoldingFeature(i, P_paired, P_unpaired);
            if (feature.first.size() == 0) continue;

            feature.first = ExpandMatrix<Real>(feature.first, L+1, L+1, positions, positions);
            features[i] = feature.second;
            
            for (int j = 1; j <= L; j++)
            {
                if (mapping[j] > j)
                    features[i] += feature.first[j*(L+1)+mapping[j]];
            }
        }
    }
    
    return features;
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::ComputeParseAlignmentFeatures()
//
// Compute alignment features for a single parse.
/////////////////////////////////////////////////////////////////

std::vector<Real> ScoringScheme::ComputeParseAlignmentFeatures(const std::vector<int> &aligned_to_A,
                                                               const std::vector<int> &aligned_to_B) const
{
    const int LA = int(projected_alignment[ids_A[0]].data.length()) - 1;
    const int LB = int(projected_alignment[ids_B[0]].data.length()) - 1;
    std::vector<Real> unsparse_mask((LA+1)*(LB+1));
    std::vector<Real> features(Weight_NUM_WEIGHTS);
    
    // consider each pair of sequences in "ids_A" and "ids_B"
    
    for (size_t kA = 0; kA < ids_A.size(); kA++)
    {
        const std::vector<int> positions_A = GetSequencePositions(projected_alignment[ids_A[kA]].data);
        for (size_t kB = 0; kB < ids_B.size(); kB++)
        {
            const std::vector<int> positions_B = GetSequencePositions(projected_alignment[ids_B[kB]].data);

            // expand posteriors

            std::vector<Real> P_aligned;
            std::vector<Real> P_unaligned_A;
            std::vector<Real> P_unaligned_B;
            std::vector<Real> pair_unsparse_mask;

            ExpandAlignmentPosteriors(ids_A[kA], ids_B[kB], positions_A, positions_B, P_aligned, P_unaligned_A, P_unaligned_B, pair_unsparse_mask);

            unsparse_mask += pair_unsparse_mask;
            
            // update features
            
            for (int i = 0; i < Weight_NUM_WEIGHTS; i++)
            {
                std::pair<std::vector<Real>, Real> feature = ComputeAlignmentFeature(i, ids_A[kA], ids_B[kB], P_aligned, P_unaligned_A, P_unaligned_B);
                if (feature.first.size() == 0) continue;
                
                feature.first = ExpandMatrix<Real>(feature.first, LA+1, LB+1, positions_A, positions_B);
                features[i] = feature.second;
                
                for (int j = 1; j <= LA; j++)
                {
                    if (aligned_to_A[j] != 0)
                        features[i] += feature.first[j*(LB+1)+aligned_to_A[j]];
                }
            }
        }
    }

    return features;
}

/////////////////////////////////////////////////////////////////
// ScoringScheme::GetGroupPositions()
//
// Compute positions of letters in a group.
/////////////////////////////////////////////////////////////////

std::vector<int> ScoringScheme::GetGroupPositions(const std::vector<Sequence> &alignment, 
                                                  const std::vector<int> &ids) const
{
    std::string s = "@";
    const int L = int(alignment[ids[0]].data.length()) - 1;
    for (int i = 1; i <= L; i++)
    {
        bool letter = false;
        for (size_t j = 0; j < ids.size(); j++)
            if (isalpha(alignment[ids[j]].data[i])) letter = true;
        s.push_back(letter ? 'X' : '-');
    }
    
    return GetSequencePositions(s);
}
