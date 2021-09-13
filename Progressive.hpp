/////////////////////////////////////////////////////////////////
// Progressive.hpp
/////////////////////////////////////////////////////////////////

#ifndef PROGRESSIVE_HPP
#define PROGRESSIVE_HPP

#include "SparseMatrix.hpp"
#include "AlignAndFold.hpp"
#include "ScoringScheme.hpp"
#include "Tree.hpp"
#include "Options.hpp"
#include "Utilities.hpp"

class Progressive 
{
    const std::vector<SparseMatrix<Real> > &base_pairing_posteriors;
    const std::vector<SparseMatrix<Real> > &alignment_posteriors;
    const bool toggle_alignment_shell;
    const bool toggle_verbose;

    void ProjectSequences(const std::vector<int> &ids,
                          std::vector<Sequence> &alignment) const;
    
    std::string AlignGroups(const std::vector<int> &ids_A,
                            const std::vector<int> &ids_B,
                            std::vector<Sequence> &projected_alignment,
                            const std::vector<Real> &weights) const;

    int ComputeMargin(const std::vector<int> &parse_A,
                      const std::vector<int> &parse_B) const;
    
public:

    Progressive(const std::vector<SparseMatrix<Real> > &base_pairing_posteriors,
                const std::vector<SparseMatrix<Real> > &alignment_posteriors,
                const bool toggle_alignment_shell,
                const bool toggle_verbose);
    
    std::string Align(const Tree::TreeNode *root,
                      std::vector<Sequence> &alignment,
                      const std::vector<Real> &weights) const;
    
    std::string IterativeRefinement(std::vector<Sequence> &alignment,
                                    const std::vector<Real> &weights,
                                    std::string consensus) const;
    
    void ComputeLossAndSubgradient(Real &loss,
                                   std::vector<Real> &subgradient,
                                   const std::vector<int> &ids_A,
                                   const std::vector<int> &ids_B,
                                   const std::vector<Sequence> &alignment,
                                   const std::string &consensus,
                                   const std::vector<Real> &weights,
                                   const Real folding_FP_margin,
                                   const Real folding_FN_margin,
                                   const Real alignment_FP_margin,
                                   const Real alignment_FN_margin) const;
    
    static void ComputeTrueParse(std::vector<int> &mapping_A,
                                 std::vector<int> &mapping_B,
                                 std::vector<int> &aligned_to_A,
                                 std::vector<int> &aligned_to_B,
                                 const std::vector<int> &ids_A,
                                 const std::vector<int> &ids_B,
                                 const std::vector<Sequence> &alignment,
                                 const std::string &consensus);
};

#endif
