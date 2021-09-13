/////////////////////////////////////////////////////////////////
// RAF.cpp
/////////////////////////////////////////////////////////////////

#include "IO.hpp"
#include "Progressive.hpp"
#include "SparseMatrix.hpp"
#include "Tree.hpp"
#include "Options.hpp"
#include "SubgradientDescent.hpp"
//#include "ProbabilisticConsistency.hpp"
//#include "PCT.hpp"

/////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////

// required arguments

bool toggle_train = false;
std::vector<std::string> input_filenames;

// miscellaneous arguments

bool toggle_verbose = false;
std::string tmp_dir_name = "";
bool toggle_alignment_shell = true;
bool toggle_noncomplementary = false;
Real alignment_posterior_cutoff = Real(0.01);
Real base_pairing_posterior_cutoff = Real(0.002);
bool toggle_vienna = false;
int min_hairpin_length = 3;
int num_iterative_refinement_reps = 0;

// environment variables

std::string vienna_rna_dir;
std::string contrafold_dir;
std::string contralign_dir;

// training 

int num_epochs = 5000;
int subset_size = -1;
Real eta = 1;
Real C1 = 0;
Real C2 = 1;
Real folding_FP_margin = 1;
Real folding_FN_margin = 10;
Real alignment_FP_margin = 1;
Real alignment_FN_margin = 1;

// default scoring weights
Real default_weights[Weight_NUM_WEIGHTS] = 
{
    5.69414,   // Weight_PAIRED_PROBABILITY
    4.48409,   // Weight_PAIRED_PROBABILITY_SQUARED
    1.67101,   // Weight_PAIRED_PROBABILITY_LOG
    3.0815,    // Weight_UNPAIRED_PROBABILITY_LOG
    2.44419,   // Weight_ALIGNED_PROBABILITY
    2.13970,   // Weight_ALIGNED_PROBABILITY_SQUARED
    0.0543566, // Weight_ALIGNED_PROBABILITY_LOG
    0.921099   // Weight_UNALIGNED_PROBABILITY_LOG
};

std::vector<Real> weights(default_weights, default_weights + Weight_NUM_WEIGHTS);

/////////////////////////////////////////////////////////////////
// ParseParameters()
//
// Parse command-line arguments.
/////////////////////////////////////////////////////////////////

void ParseParameters (int argc, char **argv)
{
    if (argc == 1)
    {
        std::cerr << std::endl
                  << "Usage: " << argv[0] << " [predict|train] INPUT_MFA(s) [OPTIONS...]" << std::endl
                  << std::endl
                  << "where INPUT_MFA is a collection of input sequences in MFA format" << std::endl
                  << std::endl
                  << "Miscellaneous arguments:" << std::endl
                  << "  --verbose                     verbose output" << std::endl
                  << "  --tmp_dir NAME                use a specific temporary directory" << std::endl
                  << "  --noshell                     do not use alignment shell optimization" << std::endl
                  << "  --noncomplementary            allow non-complementary base-pairs in pairwise alignment steps" << std::endl
                  << "  --align_cutoff VALUE          alignment posterior cutoff" << std::endl
                  << "  --bp_cutoff VALUE             base-pairing posterior cutoff" << std::endl
                  << "  --vienna                      use ViennaRNA instead of CONTRAfold" << std::endl
                  << "  --min_hairpin_len VALUE       minimum hairpin length for CONTRAfold (default: " << min_hairpin_length << ")" << std::endl
                  << std::endl
                  << "Additional arguments for 'predict' mode:" << std::endl
                  << "   --num_ir VALUE                number of iterative refinement reps" << std::endl
                  << std::endl
                  << "Additional arguments for 'train' mode:" << std::endl
                  << "   --num_epochs VALUE            number of epochs for training" << std::endl
                  << "   --subset_size VALUE           size of subset used during training" << std::endl
                  << "   --weights VALUE(s)            set all weights simultaneously (default: " << weights << ")" << std::endl
                  << "   --eta VALUE                   learning rate" << std::endl
                  << "   --margin VALUE(s)             set margin requirements" << std::endl
                  << "                                    (default: [folding_FP   = " << folding_FP_margin << "," << std::endl
                  << "                                               folding_FN   = " << folding_FN_margin << "," << std::endl
                  << "                                               alignment_FP = " << alignment_FP_margin << "," << std::endl
                  << "                                               alignment_FN = " << alignment_FN_margin << "])" << std::endl
                  << "   --regularize C1 C2            set regularization constants" << std::endl
                  << std::endl;
        exit(0);
    }

    std::vector<std::string> required_args;

    // read each argument
    
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (argv[i] == std::string("--verbose"))
            {
                toggle_verbose = true;
            }
            else if (argv[i] == std::string("--tmp_dir"))
            {
                if (++i >= argc) Error("Expected NAME after --tmp_dir.");
                tmp_dir_name = argv[i];
            }
            else if (argv[i] == std::string("--noshell"))
            {
                toggle_alignment_shell = false;
            }
            else if (argv[i] == std::string("--noncomplementary"))
            {
                toggle_noncomplementary = true;
            }
            else if (argv[i] == std::string ("--align_cutoff"))
            {
                if (++i >= argc) Error("Expected VALUE after --align_cutoff.");
                if (!ConvertToNumber(argv[i], alignment_posterior_cutoff)) Error("Could not parse value after --align_cutoff.");
            }
            else if (argv[i] == std::string("--bp_cutoff"))
            {
                if (++i >= argc) Error("Expected VALUE after --bp_cutoff.");
                if (!ConvertToNumber(argv[i], base_pairing_posterior_cutoff)) Error("Could not parse value after --bp_cutoff.");
            }
            else if (argv[i] == std::string("--vienna"))
            {
                toggle_vienna = true;
            }
            else if (argv[i] == std::string("--min_hairpin_len"))
            {
                if (++i >= argc) Error("Expected VALUE after --min_hairpin_len.");
                if (!ConvertToNumber(argv[i], min_hairpin_length)) Error("Could not parse value after --min_hairpin_len.");
                if (min_hairpin_length < 0) Error ("Minimum hairpin length should be nonnegative.");
            }
            else if (argv[i] == std::string("--num_ir"))
            {
                if (++i >= argc) Error("Expected VALUE after --num_ir.");
                if (!ConvertToNumber(argv[i], num_iterative_refinement_reps)) Error("Could not parse value after --num_ir.");
                if (num_iterative_refinement_reps < 0) Error ("Number of iterative refinement repetitions should be nonnegative.");
            }
            else if (argv[i] == std::string("--num_epochs"))
            {
                if (++i >= argc) Error("Expected VALUE after --num_epochs.");
                if (!ConvertToNumber(argv[i], num_epochs)) Error("Could not parse value after --num_epochs.");
                if (num_epochs < 0) Error ("Number of epochs should be nonnegative.");
            }
            else if (argv[i] == std::string("--subset_size"))
            {
                if (++i >= argc) Error("Expected VALUE after --subset_size.");
                if (!ConvertToNumber(argv[i], subset_size)) Error("Could not parse value after --subset_size.");
                if (subset_size <= 0) Error ("Subset size should be positive.");
            }
            else if (argv[i] == std::string("--weights"))
            {
                if (i + int(weights.size()) >= argc) Error("Not enough arguments after --weights.");
                for (size_t j = 0; j < weights.size(); j++)
                    if (!ConvertToNumber(argv[++i], weights[j])) Error("Could not parse all values after --weights.");
            }
            else if (argv[i] == std::string("--eta"))
            {
                if (++i >= argc) Error("Expected VALUE after --eta.");
                if (!ConvertToNumber(argv[i], eta)) Error("Could not parse value after --eta.");
                if (eta <= 0) Error ("Learning rate should be positive.");
            }
            else if (argv[i] == std::string("--margin"))
            {
                if (i + 4 >= argc) Error("Not enough arguments after --margin.");
                if (!ConvertToNumber(argv[++i], folding_FP_margin)) Error("Could not parse all values after --margin.");	
                if (!ConvertToNumber(argv[++i], folding_FN_margin)) Error("Could not parse all values after --margin.");	
                if (!ConvertToNumber(argv[++i], alignment_FP_margin)) Error("Could not parse all values after --margin.");	
                if (!ConvertToNumber(argv[++i], alignment_FN_margin)) Error("Could not parse all values after --margin.");	
            }
            else if (argv[i] == std::string("--regularize"))
            {
                if (i + 2 >= argc) Error("Expected VALUE after --regularize.");
                if (!ConvertToNumber(argv[++i], C1)) Error("Could not parse value after --regularize.");
                if (!ConvertToNumber(argv[++i], C2)) Error("Could not parse value after --regularize.");
            }
            else
            { 
                Error("Unknown argument: %s", argv[i]);
            }
        }
        else
        {
            required_args.push_back(argv[i]);
        }
    }

    // required_args[] should include: required_args[0] == "train" or "predict"
    //                                 required_args[1+] == input filenames
    
    if (required_args.size() < 2) Error("Incorrect number of required arguments.");
    
    if (required_args[0] != "train" && required_args[0] != "predict")
        Error("First required argument must be either \"train\" or \"predict\".");
    
    toggle_train = (required_args[0] == "train");
    
    input_filenames = std::vector<std::string>(required_args.begin() + 1, required_args.end());
    
    if (subset_size == -1 || subset_size > int(input_filenames.size()))
        subset_size = int(input_filenames.size());
}

/////////////////////////////////////////////////////////////////
// GetExternalProgramPaths()
//
// Get paths to external programs.
/////////////////////////////////////////////////////////////////

void GetExternalProgramPaths()
{
    char *ptr;

    if (toggle_vienna)
    {
        ptr = getenv("VIENNA_RNA_DIR");
        if (!ptr) Error("Cannot read environment variable 'VIENNA_RNA_DIR'.");
        vienna_rna_dir = ptr;
    }
    else
    {
        ptr = getenv("CONTRAFOLD_DIR");
        if (!ptr) Error("Cannot read environment variable 'CONTRAFOLD_DIR'.");
        contrafold_dir = ptr;
    }

    ptr = getenv("CONTRALIGN_DIR");
    if (!ptr) Error("Cannot read environment variable 'CONTRALIGN_DIR'.");
    contralign_dir = ptr;
}

/////////////////////////////////////////////////////////////////
// ComputePosteriors()
//
// Compute base-pairing and alignment match posterior matrices.
/////////////////////////////////////////////////////////////////

void ComputePosteriors(const std::vector<Sequence> &sequences,
                       std::vector<SparseMatrix<Real> > &base_pairing_posteriors,
                       std::vector<SparseMatrix<Real> > &alignment_posteriors) 
{
    const int K = int(sequences.size());
    base_pairing_posteriors.resize(K);
    alignment_posteriors.resize(K*K);
    
    // write out sequences
    
    for (int j = 0; j < K; j++)
    {
        IO::WriteToFASTA(SPrintF("%s/seq%d.fasta", tmp_dir_name.c_str(), j),
                         SPrintF("seq%d", j),
                         sequences[j].data);
    }

    // precompute base-pairing probabilities

    for (int j = 0; j < K; j++)
    {
        base_pairing_posteriors[j] =
            (toggle_vienna ? 
             IO::ComputeBasePairingPosteriorsViennaRNA(vienna_rna_dir, tmp_dir_name, j, toggle_verbose,
                                                       base_pairing_posterior_cutoff, min_hairpin_length) :
             IO::ComputeBasePairingPosteriorsCONTRAfold(contrafold_dir, tmp_dir_name, j, toggle_noncomplementary, toggle_verbose,
                                                        base_pairing_posterior_cutoff, min_hairpin_length));
    }
    
    // precompute alignment match probabilities

    for (int j = 0; j < K; j++)
    {
        for (int k = j+1; k < K; k++)
        {
            alignment_posteriors[j*K+k] = 
                IO::ComputeAlignmentPosteriorsCONTRAlign(contralign_dir, tmp_dir_name, j, k, toggle_verbose, alignment_posterior_cutoff);
        }
    }

    /*
    // precompute base-pairing probabilities

    std::vector<SparseMatrix<Real> *> temp_base_pairing_posteriors(K, static_cast<SparseMatrix<Real> *>(NULL));
    
    for (int j = 0; j < K; j++)
    {
        temp_base_pairing_posteriors[j] =
            new SparseMatrix<Real>(toggle_vienna ? 
                                   IO::ComputeBasePairingPosteriorsViennaRNA(vienna_rna_dir, tmp_dir_name, j, toggle_verbose,
                                                                             base_pairing_posterior_cutoff, min_hairpin_length) :
                                   IO::ComputeBasePairingPosteriorsCONTRAfold(contrafold_dir, tmp_dir_name, j, toggle_noncomplementary, toggle_verbose,
                                                                              base_pairing_posterior_cutoff, min_hairpin_length));
    }
    
    // precompute alignment match probabilities

    std::vector<SparseMatrix<Real> *> temp_alignment_posteriors(K*K, static_cast<SparseMatrix<Real> *>(NULL));
    
    for (int j = 0; j < K; j++)
    {
        for (int k = j+1; k < K; k++)
        {
            temp_alignment_posteriors[j*K+k] = 
                new SparseMatrix<Real>(IO::ComputeAlignmentPosteriorsCONTRAlign(contralign_dir, tmp_dir_name, j, k, toggle_verbose, alignment_posterior_cutoff));
        }
    }
    */

    /*
    // run match probabilistic consistency

    ProbabilisticConsistency<Real> consistency(temp_alignment_posteriors, toggle_verbose, 0);
    consistency.Transform();

    // perform PCT

    PCT<Real> pct(temp_base_pairing_posteriors, temp_alignment_posteriors, toggle_verbose, 0);
    pct.Transform();
    */
    // convert posteriors

    /*
    for (int j = 0; j < K; j++)
    {
        base_pairing_posteriors[j] = *temp_base_pairing_posteriors[j];
        delete temp_base_pairing_posteriors[j];
        for (int k = j+1; k < K; k++)
        {
            alignment_posteriors[j*K+k] = *temp_alignment_posteriors[j*K+k];
            delete temp_alignment_posteriors[j*K+k];
        }
    }
    */

    
}

/////////////////////////////////////////////////////////////////
// class OptimizationWrapper
//
// Wrapper routines for optimization.
/////////////////////////////////////////////////////////////////

class OptimizationWrapper : public SubgradientDescent<double>
{
    const std::vector<std::string> &consensus;
    const std::vector<std::vector<Sequence> > &alignments;
    const std::vector<std::vector<Sequence> > &sequences;
    const std::vector<std::vector<SparseMatrix<Real> > > &base_pairing_posteriors;
    const std::vector<std::vector<SparseMatrix<Real> > > &alignment_posteriors;

    std::vector<double> cached_x;
    std::vector<double> cached_g;
    double cached_f;

    // compute function and subgradient

    void ComputeFunctionAndSubgradient(const std::vector<double> &x)
    {
        std::vector<Real> converted_x = ConvertVector<Real,double>(x);

        cached_f = 0;
        cached_x = x;
        cached_g.clear();
        cached_g.resize(x.size(), 0.0);

        // pick random sample
        
        std::vector<int> examples(subset_size, 1);
        examples.resize(alignments.size(), 0);
        std::random_shuffle(examples.begin(), examples.end());

        // for each example
        
        for (size_t i = 0; i < examples.size(); i++)
        {
            if (examples[i] == 0) continue;

            double example_loss = 0.0;
            std::vector<double> example_subgradient(x.size());
            const int K = int(alignments[i].size());
            
            Progressive progressive(base_pairing_posteriors[i], alignment_posteriors[i],
                                    toggle_alignment_shell, toggle_verbose);

            if (subset_size == int(examples.size()))
            {
                // for each pair of sequences
                
                for (int j = 0; j < K; j++)
                {
                    for (int k = j+1; k < K; k++)
                    {
                        Real pair_loss;
                        std::vector<Real> pair_subgradient;
                        
                        progressive.ComputeLossAndSubgradient(pair_loss, pair_subgradient,
                                                              std::vector<int>(1, j), std::vector<int>(1, k), 
                                                              alignments[i], consensus[i], converted_x,
                                                              folding_FP_margin, folding_FN_margin,
                                                              alignment_FP_margin, alignment_FN_margin);
                        example_loss += double(pair_loss);
                        example_subgradient += ConvertVector<double,Real>(pair_subgradient);
                    }
                }
                example_loss /= double(K*(K-1)/2);
                example_subgradient /= double(K*(K-1)/2);
            }
            else
            {
                int j = lrand48() % K;
                int k = lrand48() % (K-1);
                if (k >= j) k++;
                if (j > k) std::swap(j, k);

                Real pair_loss;
                std::vector<Real> pair_subgradient;
                
                progressive.ComputeLossAndSubgradient(pair_loss, pair_subgradient,
                                                      std::vector<int>(1, j), std::vector<int>(1, k), 
                                                      alignments[i], consensus[i], converted_x,
                                                      folding_FP_margin, folding_FN_margin,
                                                      alignment_FP_margin, alignment_FN_margin);
                example_loss += double(pair_loss);
                example_subgradient += ConvertVector<double,Real>(pair_subgradient);
            }
            
            cached_f += example_loss;
            cached_g += example_subgradient;
        }
            
        cached_f /= double(subset_size);
        cached_g /= double(subset_size);

        // add regularization penalty

        cached_f += double(C1) * Sum(Abs(x)) + double(0.5 * C2) * DotProduct(x, x);
        cached_g += double(C1) * Sign(x) + double(C2) * x;
    }

    /////////////////////////////////////////////////////////////////
    // ComputeNormBound()
    //
    // Compute a bound on the maximum possible parameter norm.  To do this,
    // we'll use strong duality.  The optimization problem is
    //
    //    min   (1/2) C ||w||^2 + (1/m) \sum_{i=1}^m \xi_i
    //    s.t.  \xi_i \ge 0
    //          w^T F(x[i], y[i]) \ge w^T F(x[i], y) + MARGIN(y[i],y) - \xi_i
    //
    // whose dual is
    //
    //    max   \sum_{i=1}^m \sum_y \alpha_{i,y} MARGIN(y[i],y) - (1/2) C ||w'||^2
    //    s.t.  \alpha_{i,y} \ge 0
    //          \sum_y \alpha_{i,y} = 1/m
    //
    // where w' = (1/C) \sum_{i=1}^m \sum_y \alpha_{i,y} (F(x[i],y[i]) - F(x[i],y)).
    //
    // Using strong duality, the primal and dual objectives are equal (with w = w') so
    //
    // C ||w||^2 = \sum_{i=1}^m \sum_y \alpha_{i,y} MARGIN(y[i],y) - (1/m) \sum_{i=1}^m \xi_i
    //                          +--------------------------------+                      +----+
    //                              <= (1/m) WORSTMARGIN(y[i])                           >= 0
    //
    // so ||w||^2 <= (1/C) * (1/m) \sum_{i=1}^m WORSTMARGIN(y[i])
    /////////////////////////////////////////////////////////////////
    
    double ComputeNormBound(const std::vector<std::vector<Sequence> > &sequences)
    {
        // compute average of the worst margins for each example
        
        double overall_worst_margin = 0.0;
        for (size_t i = 0; i < sequences.size(); i++)
        {
            // compute average of the worst margins for all pairs within current example
            
            const int K = int(sequences[i].size());
            double example_worst_margin = 0.0;
            
            for (int j = 0; j < K; j++)
            {
                for (int k = j+1; k < K; k++)
                {
                    const int LA = int(sequences[i][j].data.length()) - 1;
                    const int LB = int(sequences[i][k].data.length()) - 1;
                    
                    double pair_worst_margin =
                        (folding_FP_margin + folding_FN_margin + alignment_FP_margin + alignment_FN_margin) * (LA + LB);
                    
                    example_worst_margin += pair_worst_margin;
                }
            }
            example_worst_margin /= double(K*(K-1)/2);
            overall_worst_margin += example_worst_margin;
        }
        overall_worst_margin /= double(sequences.size());
        
        // now finish norm bound computation
        
        return Sqrt(overall_worst_margin / C2);
    }

public:
    
    // constructor
    
    OptimizationWrapper(const std::vector<std::string> &consensus,
                        const std::vector<std::vector<Sequence> > &alignments,
                        const std::vector<std::vector<Sequence> > &sequences,
                        const std::vector<std::vector<SparseMatrix<Real> > > &base_pairing_posteriors,
                        const std::vector<std::vector<SparseMatrix<Real> > > &alignment_posteriors) :
        SubgradientDescent<double>(eta / C2, num_epochs, ComputeNormBound(sequences), 0.0),
        consensus(consensus), alignments(alignments), sequences(sequences),
        base_pairing_posteriors(base_pairing_posteriors), alignment_posteriors(alignment_posteriors) {}

    // report results

    void Report(int iteration, double f, const std::vector<double> &x, const std::vector<double> &g,
                double norm_bound, double step_size)
    {
        std::cerr << "Epoch " << iteration
                  << ", w = " << x
                  << ", ||w|| = " << Norm(x)
                  << ", ||g|| = " << Norm(g)
                  << ", norm bound = " << norm_bound
                  << ", f = " << f
                  << ", alpha = " << step_size << std::endl;
    }

    void Report(const std::string &s) 
    {
        std::cerr << s << std::endl;
    }
   
    // compute function

    double ComputeFunction(const std::vector<double> &x)
    {
        if (x != cached_x) ComputeFunctionAndSubgradient(x);
        return cached_f;
    }

    // compute subgradient
    
    void ComputeSubgradient(std::vector<double> &g, const std::vector<double> &x)
    {
        if (x != cached_x) ComputeFunctionAndSubgradient(x);
        g = cached_g;
    }

};

/////////////////////////////////////////////////////////////////
// TrainAligner()
//
// Perform subgradient descent training.
/////////////////////////////////////////////////////////////////

void TrainAligner()
{
    std::vector<std::string> consensus(input_filenames.size());
    std::vector<std::vector<Sequence> > alignments(input_filenames.size());
    std::vector<std::vector<Sequence> > sequences(input_filenames.size());
    std::vector<std::vector<SparseMatrix<Real> > > base_pairing_posteriors(input_filenames.size());
    std::vector<std::vector<SparseMatrix<Real> > > alignment_posteriors(input_filenames.size());

    // read alignments and compute posteriors

    for (size_t i = 0; i < input_filenames.size(); i++)
    {
        std::cerr << "Processing " << input_filenames[i] << " ..." << std::endl;
        alignments[i] = IO::ReadAlignmentWithConsensus(input_filenames[i], consensus[i]);
        sequences[i] = IO::StripGaps(alignments[i]);
        ComputePosteriors(sequences[i], base_pairing_posteriors[i], alignment_posteriors[i]);
    }

    OptimizationWrapper optimizer(consensus, alignments, sequences, base_pairing_posteriors, alignment_posteriors);
    std::vector<double> w = ConvertVector<double, Real>(weights);
    optimizer.Minimize(w);
    weights = ConvertVector<Real, double>(w);
}

/////////////////////////////////////////////////////////////////
// DoAlignment()
//
// Perform simultaneous multiple alignment and folding.
/////////////////////////////////////////////////////////////////

void DoAlignment(const std::string &input_filename)
{
    std::vector<Sequence> sequences = IO::ReadUnalignedSequences(input_filename);
    std::vector<SparseMatrix<Real> > base_pairing_posteriors;
    std::vector<SparseMatrix<Real> > alignment_posteriors;
    const int K = int(sequences.size());
    
    // compute posterior probabilities
    
    ComputePosteriors(sequences, base_pairing_posteriors, alignment_posteriors); 
   
    // compute distances for each pair of sequences
    
    std::vector<std::vector<Real> > distance(K, std::vector<Real>(K));
    for (int i = 0; i < K; i++)
    {
        for (int j = i+1; j < K; j++)
        {
            distance[i][j] = 
                Real(1.0) - alignment_posteriors[i*K+j].GetSum() / 
                std::min(alignment_posteriors[i*K+j].GetNumRows() - 1,
                         alignment_posteriors[i*K+j].GetNumCols() - 1);
        }
    }

    // construct tree
    
    Tree::TreeNode *root = Tree::UPGMA(distance);
    if (toggle_verbose)
    {
        std::cerr << "Tree: ";
        root->Print(std::cerr);
    }

    // progressive alignment
    
    Progressive progressive(base_pairing_posteriors, alignment_posteriors,
                            toggle_alignment_shell, toggle_verbose);
    
    std::vector<Sequence> alignment(sequences);
    std::string consensus = progressive.Align(root, alignment, weights);
    
    delete root;

    // iterative refinement
    
    for (int i = 0; i < num_iterative_refinement_reps; i++)
        consensus = progressive.IterativeRefinement(alignment, weights, consensus);
    
    // print results
    
    for (size_t i = 0; i < alignment.size(); i++)
    {
        std::cout << ">" << alignment[i].name << std::endl;
        std::cout << alignment[i].data.substr(1) << std::endl;
    }
    std::cout << ">consensus" << std::endl;
    std::cout << consensus.substr(1) << std::endl;
}

/////////////////////////////////////////////////////////////////
// AccumulateAlignmentStatistics()
//
// Compute statistics for a single alignment.
/////////////////////////////////////////////////////////////////

void AccumulateAlignmentStatistics(const std::vector<int> &aligned_to_A,
                                   const std::vector<int> &aligned_to_B,
                                   const SparseMatrix<Real> &posteriors,
                                   const std::vector<Real> &epsilon,
                                   std::vector<Real> &alignment_sparsity_ratio,
                                   std::vector<Real> &alignment_coverage,
                                   std::vector<Real> &alignment_accuracy,
                                   Real weight)
{
    Real num_true_aligned_positions = 0;
    for (size_t i = 1; i < aligned_to_A.size(); i++)
        if (aligned_to_A[i] != 0) num_true_aligned_positions++;
    
    std::vector<Real> num_candidate_aligned_positions(epsilon.size());
    std::vector<Real> num_covered_aligned_positions(epsilon.size());

    for (size_t e = 0; e < epsilon.size(); e++)
    {
        int max_count = 0;
        for (size_t i = 1; i < aligned_to_A.size(); i++)
        {
            int first = -1;
            int last = -1;

            // find valid range
            
            for (const SparseMatrixEntry<Real> *ptr = posteriors.GetRowBegin(i); ptr != posteriors.GetRowEnd(i); ++ptr)
            {
                if (ptr->value > epsilon[e])
                {
                    if (first == -1) first = ptr->column;
                    last = ptr->column;
                }
            }

            // check if non-empty

            if (last != -1)
            {
                max_count = std::max(max_count, last - first + 1);
                if (aligned_to_A[i] > 0 && first <= aligned_to_A[i] && aligned_to_A[i] <= last)
                    num_covered_aligned_positions[e]++;
            }
        }
        num_candidate_aligned_positions[e] = Real(max_count);
    }

    alignment_sparsity_ratio += num_candidate_aligned_positions * weight;
    alignment_coverage += num_covered_aligned_positions / num_true_aligned_positions * weight;
    alignment_accuracy = alignment_coverage;

    /*
    for (size_t i = 1; i < aligned_to_A.size(); i++)
    {
        for (const SparseMatrixEntry<Real> *ptr = posteriors.GetRowBegin(i); ptr != posteriors.GetRowEnd(i); ++ptr)
        {
            for (int e = int(epsilon.size()) - 1; e >= 0; e--)
            {
                if (epsilon[e] > ptr->value) break;
                num_candidate_aligned_positions[e]++;
                if (aligned_to_A[i] == ptr->column)
                    num_covered_aligned_positions[e]++;
            }
        }
    }
    
    alignment_sparsity_ratio += num_candidate_aligned_positions * Real(2) / Real(aligned_to_A.size() + aligned_to_B.size() - 2) * weight;
    alignment_coverage += num_covered_aligned_positions / num_true_aligned_positions * weight;
    alignment_accuracy += (num_covered_aligned_positions + Real(1e-7)) / (num_candidate_aligned_positions + Real(1e-7)) * weight;
    */
}

/////////////////////////////////////////////////////////////////
// AccumulatePairingStatistics()
//
// Compute statistics for a single structure.
/////////////////////////////////////////////////////////////////

void AccumulatePairingStatistics(const std::vector<int> &mapping,
                                 const SparseMatrix<Real> &posteriors,
                                 const std::vector<Real> &epsilon,
                                 std::vector<Real> &pairing_sparsity_ratio,
                                 std::vector<Real> &pairing_coverage,
                                 std::vector<Real> &pairing_accuracy,
                                 Real weight)
{
    Real num_true_pairings = 0;
    for (int i = 1; i < int(mapping.size()); i++)
        if (mapping[i] > i) num_true_pairings++;
    
    std::vector<Real> num_candidate_pairings(epsilon.size());
    std::vector<Real> num_covered_pairings(epsilon.size());

    for (size_t e = 0; e < epsilon.size(); e++)
    {
        int max_count = 0;
        for (size_t i = 1; i < mapping.size(); i++)
        {
            int count = 0;
            for (const SparseMatrixEntry<Real> *ptr = posteriors.GetRowBegin(i); ptr != posteriors.GetRowEnd(i); ++ptr)
            {
                if (ptr->value > epsilon[e])
                {
                    count++;
                    if (mapping[i] == ptr->column)
                        num_covered_pairings[e]++;
                }
            }

            max_count = std::max(max_count, count);
        }
        num_candidate_pairings[e] = Real(max_count);
    }

    pairing_sparsity_ratio += num_candidate_pairings * weight;
    pairing_coverage += (num_covered_pairings + Real(1e-7)) / (num_true_pairings + Real(1e-7)) * weight;
    pairing_accuracy = pairing_coverage;

    /*
    for (int i = 1; i < int(mapping.size()); i++)
    {
        for (const SparseMatrixEntry<Real> *ptr = posteriors.GetRowBegin(i); ptr != posteriors.GetRowEnd(i); ++ptr)
        {
            for (int e = int(epsilon.size()) - 1; e >= 0; e--)
            {
                if (epsilon[e] > ptr->value) break;
                num_candidate_pairings[e]++;
                if (mapping[i] == ptr->column)
                    num_covered_pairings[e]++;
            }
        }
    }
    
    pairing_sparsity_ratio += num_candidate_pairings / Real(mapping.size() - 1) * weight;
    pairing_coverage += num_covered_pairings / num_true_pairings * weight;
    pairing_accuracy += (num_covered_pairings + Real(1e-7)) / (num_candidate_pairings + Real(1e-7)) * weight;
    */
}

/////////////////////////////////////////////////////////////////
// ComputeStatistics()
//
// Compute coverage statistics for the alignment and base-pairing
// probability matrices.
/////////////////////////////////////////////////////////////////

void ComputeStatistics()
{
    alignment_posterior_cutoff = 1e-6;
    base_pairing_posterior_cutoff = 1e-6;
    
    std::vector<std::string> consensus(input_filenames.size());
    std::vector<std::vector<Sequence> > alignments(input_filenames.size());
    std::vector<std::vector<Sequence> > sequences(input_filenames.size());
    std::vector<std::vector<SparseMatrix<Real> > > base_pairing_posteriors(input_filenames.size());
    std::vector<std::vector<SparseMatrix<Real> > > alignment_posteriors(input_filenames.size());

    // read alignments and compute posteriors

    for (size_t i = 0; i < input_filenames.size(); i++)
    {
        std::cerr << "Processing " << input_filenames[i] << " ..." << std::endl;
        alignments[i] = IO::ReadAlignmentWithConsensus(input_filenames[i], consensus[i]);
        sequences[i] = IO::StripGaps(alignments[i]);
        ComputePosteriors(sequences[i], base_pairing_posteriors[i], alignment_posteriors[i]);
    }

    // initialization

    const int NUM_POINTS = 100;
    const Real RATIO = 0.9;
    std::vector<Real> epsilon(NUM_POINTS);
    std::vector<Real> alignment_sparsity_ratio(NUM_POINTS);
    std::vector<Real> base_pairing_sparsity_ratio(NUM_POINTS);
    std::vector<Real> alignment_coverage(NUM_POINTS);
    std::vector<Real> base_pairing_coverage(NUM_POINTS);
    std::vector<Real> alignment_accuracy(NUM_POINTS);
    std::vector<Real> base_pairing_accuracy(NUM_POINTS);
    epsilon[0] = Real(1);
    for (size_t i = 1; i < epsilon.size(); i++)
        epsilon[i] = epsilon[i-1] * RATIO;

    // iterate through sets
    
    for (size_t i = 0; i < input_filenames.size(); i++)
    {
        const int K = int(alignments[i].size());

        // compute alignment statistics
        
        for (int j = 0; j < K; j++)
        {
            for (int k = j+1; k < K; k++)
            {
                std::vector<int> mapping_A;
                std::vector<int> mapping_B;
                std::vector<int> aligned_to_A;
                std::vector<int> aligned_to_B;
                                
                Progressive::ComputeTrueParse(mapping_A, mapping_B, aligned_to_A, aligned_to_B,
                                              std::vector<int>(1,j), std::vector<int>(1,k),
                                              alignments[i], consensus[i]);
                
                AccumulateAlignmentStatistics(aligned_to_A, aligned_to_B, alignment_posteriors[i][j*K+k], epsilon,
                                              alignment_sparsity_ratio, alignment_coverage, alignment_accuracy, Real(2) / Real(K*(K-1)));
            }
        }

        // compute base pairing statistics

        for (int j = 0; j < K; j += 2)
        {
            std::vector<int> mapping_A;
            std::vector<int> mapping_B;
            std::vector<int> aligned_to_A;
            std::vector<int> aligned_to_B;
            
            Progressive::ComputeTrueParse(mapping_A, mapping_B, aligned_to_A, aligned_to_B,
                                          std::vector<int>(1,j), std::vector<int>(1,(j+1)%K),
                                          alignments[i], consensus[i]);
            AccumulatePairingStatistics(mapping_A, base_pairing_posteriors[i][j], epsilon,
                                        base_pairing_sparsity_ratio, base_pairing_coverage, base_pairing_accuracy, Real(1) / Real(K));
            if (j+1 < K)
                AccumulatePairingStatistics(mapping_B, base_pairing_posteriors[i][j+1], epsilon,
                                            base_pairing_sparsity_ratio, base_pairing_coverage, base_pairing_accuracy, Real(1) / Real(K));
        }
    }

    alignment_sparsity_ratio /= Real(input_filenames.size());
    alignment_coverage /= Real(input_filenames.size());
    alignment_accuracy /= Real(input_filenames.size());
    base_pairing_sparsity_ratio /= Real(input_filenames.size());
    base_pairing_coverage /= Real(input_filenames.size());
    base_pairing_accuracy /= Real(input_filenames.size());

    // print results

    printf("%20s%20s%20s%20s%20s%20s%20s\n", "EPSILON", "ALIGNMENT_SPARSITY", "ALIGNMENT_COVERAGE", "ALIGNMENT_ACCURACY", "PAIRING_SPARSITY", "PAIRING_COVERAGE", "PAIRING_ACCURACY");
    for (size_t i = 0; i < epsilon.size(); i++)
    {
        printf("%20lf%20lf%20lf%20lf%20lf%20lf%20lf\n", double(epsilon[i]),
               double(alignment_sparsity_ratio[i]), double(alignment_coverage[i]), double(alignment_accuracy[i]),
               double(base_pairing_sparsity_ratio[i]), double(base_pairing_coverage[i]), double(base_pairing_accuracy[i]));
    }
    
}


/////////////////////////////////////////////////////////////////
// main()
//
// Main program.
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    ParseParameters(argc, argv);
    GetExternalProgramPaths();

    // create temporary directory
    
    if (tmp_dir_name == "")
        tmp_dir_name = MakeTempDirectory();
    else
        MakeDirectory(tmp_dir_name);
    
    //ComputeStatistics();

    if (toggle_train)
        TrainAligner();
    else
    {
        for (size_t i = 0; i < input_filenames.size(); i++)
            DoAlignment(input_filenames[i]);
    }
        
    return 0;
}
