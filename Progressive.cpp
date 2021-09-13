/////////////////////////////////////////////////////////////////
// Progressive.cpp
/////////////////////////////////////////////////////////////////

#include "Progressive.hpp"

/////////////////////////////////////////////////////////////////
// Progressive::Progressive()
//
// Constructor.
/////////////////////////////////////////////////////////////////

Progressive::Progressive(const std::vector<SparseMatrix<Real> > &base_pairing_posteriors,
                         const std::vector<SparseMatrix<Real> > &alignment_posteriors,
                         const bool toggle_alignment_shell,
                         const bool toggle_verbose) :
    base_pairing_posteriors(base_pairing_posteriors), 
    alignment_posteriors(alignment_posteriors),
    toggle_alignment_shell(toggle_alignment_shell),
    toggle_verbose(toggle_verbose)
{}

/////////////////////////////////////////////////////////////////
// Progressive::Align()
//
// Perform progressive alignment.
/////////////////////////////////////////////////////////////////

std::string Progressive::Align(const Tree::TreeNode *root,
                               std::vector<Sequence> &alignment,
                               const std::vector<Real> &weights) const
{
    Assert(root, "Null pointer reached.");
    Assert(root->ids.size() > 0, "Incorrectly formed TreeNode().");

    if (root->ids.size() == 1) return "";
    
    Align(root->left_child, alignment, weights);
    Align(root->right_child, alignment, weights);
    
    return AlignGroups(root->left_child->ids, root->right_child->ids, alignment, weights);
}

/////////////////////////////////////////////////////////////////
// Progressive::IterativeRefinement()
//
// Perform iterative refinement.
/////////////////////////////////////////////////////////////////

std::string Progressive::IterativeRefinement(std::vector<Sequence> &alignment,
                                             const std::vector<Real> &weights,
                                             std::string consensus) const
{
    std::vector<int> ids_A(1);
    std::vector<int> ids_B;  
    
    for (size_t i = 0; i < alignment.size(); i++)
    {
        ids_A[0] = i;
        ids_B.clear();
        for (size_t j = 0; j < alignment.size(); j++)
            if (j != i) ids_B.push_back (j);
        
        ProjectSequences(ids_A, alignment);
        ProjectSequences(ids_B, alignment);
        
        consensus = AlignGroups(ids_A, ids_B, alignment, weights);
    }
    
    return consensus;
}

/////////////////////////////////////////////////////////////////
// Progressive::ProjectSequences()
//
// Project alignment by removing gaps from a sequence subset.
/////////////////////////////////////////////////////////////////

void Progressive::ProjectSequences(const std::vector<int> &ids,
                                   std::vector<Sequence> &alignment) const
{
    // compute remove mask
    
    std::vector<int> remove(alignment[ids[0]].data.length(), 1);
    for (size_t i = 0; i < ids.size(); i++)
        for (size_t j = 1; j < alignment[ids[i]].data.length(); j++)
            if (isalpha(alignment[ids[i]].data[j])) remove[j] = 0;
    
    // perform projection
    
    for (size_t i = 0; i < ids.size(); i++)
    {
        std::string new_string = "@";
        for (size_t j = 1; j < alignment[ids[i]].data.length(); j++)
            if (!remove[j]) new_string.push_back(alignment[ids[i]].data[j]);
        alignment[ids[i]].data = new_string;
    }
}

/////////////////////////////////////////////////////////////////
// Progressive::AlignGroups()
//
// Align two groups of sequences.
/////////////////////////////////////////////////////////////////

std::string Progressive::AlignGroups(const std::vector<int> &ids_A,
                                     const std::vector<int> &ids_B,
                                     std::vector<Sequence> &projected_alignment,
                                     const std::vector<Real> &weights) const
{
    double starting_time = 0;
    if (toggle_verbose)
    {
        starting_time = GetSystemTime();
        std::ostringstream oss;
        oss << "Aligning groups: " << ids_A << " vs. " << ids_B;
        WriteProgressMessage(oss.str());
    }

    // initialize data structure for scoring scheme and masks

    ScoringScheme scoring(base_pairing_posteriors, alignment_posteriors,
                          ids_A, ids_B, projected_alignment, weights);

    // perform alignment
    
    AlignAndFold aligner(scoring,
                         toggle_alignment_shell,
                         toggle_verbose);
    
    Real score = aligner.DoAlignment();

    // save alignment

    std::string edit_string = aligner.GetEditString();

    for (size_t k = 0; k < ids_A.size(); k++)
        projected_alignment[ids_A[k]].data = aligner.FormatString(projected_alignment[ids_A[k]].data, edit_string, 'A', '-');
    for (size_t k = 0; k < ids_B.size(); k++)
        projected_alignment[ids_B[k]].data = aligner.FormatString(projected_alignment[ids_B[k]].data, edit_string, 'B', '-');

    // print feedback

    if (toggle_verbose)
    {
        WriteProgressMessage("");
        std::cerr << "Aligned groups: " << ids_A << " vs. " << ids_B 
                  << ": " << GetSystemTime() - starting_time << " seconds (score = " << score << ")" << std::endl;
    }
        
    return aligner.GetConsensusStructure();
}

/////////////////////////////////////////////////////////////////
// Progressive::ComputeLossAndSubgradient()
//
// Compute loss and subgradient for the alignment of two groups
// of sequences.
/////////////////////////////////////////////////////////////////

void Progressive::ComputeLossAndSubgradient(Real &loss,
                                            std::vector<Real> &subgradient,
                                            const std::vector<int> &ids_A,
                                            const std::vector<int> &ids_B,
                                            const std::vector<Sequence> &alignment,
                                            const std::string &consensus,
                                            const std::vector<Real> &weights,
                                            const Real folding_FP_margin,
                                            const Real folding_FN_margin,
                                            const Real alignment_FP_margin,
                                            const Real alignment_FN_margin) const
{

#ifndef NDEBUG
    const int L = int(alignment[ids_A[0]].data.length()) - 1;
    Assert(L >= 0, "Invalid sequence length.");
    for (size_t i = 1; i < ids_A.size(); i++)
        Assert(L == int(alignment[ids_A[i]].data.length()) - 1, "Not all sequences have the same length.");
    for (size_t i = 0; i < ids_B.size(); i++)
        Assert(L == int(alignment[ids_B[i]].data.length()) - 1, "Not all sequences have the same length.");
    Assert(L == int(consensus.length()) - 1, "Consensus length does not match alignment.");
#endif

    // make projected copy of alignment in which A and B are not aligned
    
    std::vector<Sequence> projected_alignment(alignment);
    ProjectSequences(ids_A, projected_alignment);
    ProjectSequences(ids_B, projected_alignment);

    // get true parse
    
    std::vector<int> true_mapping_A;
    std::vector<int> true_mapping_B;
    std::vector<int> true_aligned_to_A;
    std::vector<int> true_aligned_to_B;
    
    ComputeTrueParse(true_mapping_A, true_mapping_B,
                     true_aligned_to_A, true_aligned_to_B,
                     ids_A, ids_B, alignment, consensus);

    // initialize loss augmented scoring scheme

    ScoringScheme scoring(base_pairing_posteriors, alignment_posteriors,
                          ids_A, ids_B, projected_alignment, weights);

    ScoringScheme loss_augmented_scoring(scoring);
    
    loss_augmented_scoring.AugmentWithLoss(true_mapping_A, true_mapping_B,
                                           true_aligned_to_A, true_aligned_to_B,
                                           folding_FP_margin, folding_FN_margin,
                                           alignment_FP_margin, alignment_FN_margin);
    
    // perform alignment
    
    AlignAndFold aligner(loss_augmented_scoring,
                         toggle_alignment_shell,
                         toggle_verbose);
    
#ifndef NDEBUG
    Real predicted_score =
#endif
        aligner.DoAlignment();
    
    // get predicted parse
    
    std::vector<int> predicted_mapping_A = aligner.GetMappingA();
    std::vector<int> predicted_mapping_B = aligner.GetMappingB();
    std::vector<int> predicted_aligned_to_A = aligner.GetAlignedToA();
    std::vector<int> predicted_aligned_to_B = aligner.GetAlignedToB();

    // convert true and predicted parses into features
    
    std::vector<Real> true_features =
        scoring.ComputeFeatures(true_mapping_A, true_mapping_B,
                                true_aligned_to_A, true_aligned_to_B);
    
    std::vector<Real> predicted_features =
        scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B,
                                predicted_aligned_to_A, predicted_aligned_to_B);

    // compute margin suffered from differences in predicted and true mapping

    Real margin =
        folding_FP_margin * (ComputeMargin(predicted_mapping_A, true_mapping_A) + ComputeMargin(predicted_mapping_B, true_mapping_B)) +
        folding_FN_margin * (ComputeMargin(true_mapping_A, predicted_mapping_A) + ComputeMargin(true_mapping_B, predicted_mapping_B)) +
        alignment_FP_margin * (ComputeMargin(predicted_aligned_to_A, true_aligned_to_A) + ComputeMargin(predicted_aligned_to_B, true_aligned_to_B)) +
        alignment_FN_margin * (ComputeMargin(true_aligned_to_A, predicted_aligned_to_A) + ComputeMargin(true_aligned_to_B, predicted_aligned_to_B));

    // compute subgradient and objective function

    subgradient = predicted_features - true_features;
    loss = DotProduct(weights, subgradient) + margin;

    // check for correctness

#ifndef NDEBUG
    Real verify_predicted_score = margin + DotProduct(predicted_features, weights);

    std::cerr << "predicted_score = " << predicted_score << std::endl;
    std::cerr << "verify_predicted_score = " << verify_predicted_score << std::endl;
    std::cerr << std::endl;
    std::cerr << "weights = " << weights << std::endl;
    std::cerr << std::endl;
    std::cerr << "Loss-augmented:" << std::endl;
    std::cerr << "   reported score = " << predicted_score << std::endl;
    std::cerr << "   verified score = " << aligner.ComputeScore(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B) << std::endl;
    std::cerr << "   offset = " << loss_augmented_scoring.GetOffset() << std::endl;
    std::cerr << "   features = " << loss_augmented_scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B) << std::endl;
    std::cerr << "   w^T F = " << DotProduct(loss_augmented_scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B), weights) << std::endl;
    std::cerr << "   w^T F + offset = " << DotProduct(loss_augmented_scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B), weights) + loss_augmented_scoring.GetOffset() << std::endl;
    std::cerr << std::endl;

    AlignAndFold aligner2(scoring,
                          toggle_alignment_shell,
                          toggle_verbose);
    
    std::cerr << "Regular:" << std::endl;
    std::cerr << "   verified score = " << aligner2.ComputeScore(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B) << std::endl;
    std::cerr << "   offset = " << scoring.GetOffset() << std::endl;
    std::cerr << "   features = " << scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B) << std::endl;
    std::cerr << "   w^T F = " << DotProduct(scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B), weights) << std::endl;
    std::cerr << "   w^T F + offset = " << DotProduct(scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B), weights) + scoring.GetOffset() << std::endl;
    std::cerr << "   margin = " << margin << std::endl;
    std::cerr << "   w^T F + offset + margin = " << DotProduct(scoring.ComputeFeatures(predicted_mapping_A, predicted_mapping_B, predicted_aligned_to_A, predicted_aligned_to_B), weights) + scoring.GetOffset() + margin << std::endl;
    
    Assert(Abs(predicted_score - verify_predicted_score) / std::max(Real(1), Abs(predicted_score) + Abs(verify_predicted_score)) < Real(1e-5), "Scores do not match.");
#endif

    // special case if true parse does not belong to restricted region
    
    if (loss < 0)
    {
        loss = 0;
        std::fill(subgradient.begin(), subgradient.end(), Real(0));
    }
}

/////////////////////////////////////////////////////////////////
// Progressive::ComputeTrueParse()
//
// Compute true parse for the alignment of two groups of
// sequences.
/////////////////////////////////////////////////////////////////

void Progressive::ComputeTrueParse(std::vector<int> &mapping_A,
                                   std::vector<int> &mapping_B,
                                   std::vector<int> &aligned_to_A,
                                   std::vector<int> &aligned_to_B,
                                   const std::vector<int> &ids_A,
                                   const std::vector<int> &ids_B,
                                   const std::vector<Sequence> &alignment,
                                   const std::string &consensus) 
{
    const int L = int(alignment[ids_A[0]].data.length()) - 1;
    
    // compute occupancy
    
    std::vector<int> occupied_A(L+1);
    std::vector<int> occupied_B(L+1);
    
    for (size_t k = 0; k < ids_A.size(); k++)
        for (int i = 1; i <= L; i++) 
            if (alignment[ids_A[k]].data[i] != '-') occupied_A[i] = 1;
    
    for (size_t k = 0; k < ids_B.size(); k++)
        for (int i = 1; i <= L; i++) 
            if (alignment[ids_B[k]].data[i] != '-') occupied_B[i] = 1;
    
    // analyze alignment
    
    const int LA = Sum(occupied_A);
    const int LB = Sum(occupied_B);

    // clear results
    
    aligned_to_A.clear(); aligned_to_A.resize(LA+1);
    aligned_to_B.clear(); aligned_to_B.resize(LB+1);
    mapping_A.clear(); mapping_A.resize(LA+1);
    mapping_B.clear(); mapping_B.resize(LB+1);

    // compute results
    
    std::vector<std::pair<int,int> > stack;
    
    int i = 0;
    int j = 0;
    for (int k = 1; k <= L; k++)
    {
        i += occupied_A[k];
        j += occupied_B[k];

        // aligned nucleotides
        
        if (occupied_A[k] && occupied_B[k])
        {
            aligned_to_A[i] = j;
            aligned_to_B[j] = i;
        }

        // base-pairings
        
        if (consensus[k] == '(')
        {
            stack.push_back(std::make_pair(occupied_A[k] ? i : 0, occupied_B[k] ? j : 0));
        }
        else if (consensus[k] == ')')
        {
            if (stack.size() == 0) Error("Bad structure in consensus file.");
            if (stack.back().first > 0)
            {
                mapping_A[stack.back().first] = i;
                mapping_A[i] = stack.back().first;
            }
            if (stack.back().second > 0)
            {
                mapping_B[stack.back().second] = j;
                mapping_B[j] = stack.back().second;
            }
            stack.pop_back();
        }
    }
    Assert(i == LA && j == LB, "Lengths should match.");
}

/////////////////////////////////////////////////////////////////
// Progressive::ComputeMargin()
//
// Compute number of pairings in A which are not found in B.
/////////////////////////////////////////////////////////////////

int Progressive::ComputeMargin(const std::vector<int> &parse_A,
                               const std::vector<int> &parse_B) const
{
    Assert(parse_A.size() == parse_B.size(), "Parses not of the same length.");
    
    int loss = 0;
    for (size_t i = 0; i < parse_A.size(); i++)
    {
        if (parse_A[i] != 0)
        {
            loss += Ind(parse_B[i] != parse_A[i]);                
        }
    }
    
    return loss;
}
