/////////////////////////////////////////////////////////////////
// AlignAndFold.hpp
/////////////////////////////////////////////////////////////////

#ifndef ALIGNANDFOLD_HPP
#define ALIGNANDFOLD_HPP

#include "Options.hpp"
#include "AlignmentShell.hpp"
#include "ScoringScheme.hpp"
#include "SparseMatrix.hpp"
#include "Utilities.hpp"

struct OuterEntry
{
    Real score;
    int traceback;
    
    OuterEntry() : score(0), traceback(0) {}
    OuterEntry(Real score) : score(score), traceback(0) {}
    OuterEntry(Real score, int traceback) : score(score), traceback(traceback) {}
};

struct AlignedPair
{
    union {
        int i;
        int j;
    };

    union {
        int k;
        int l;
    };
    
    AlignedPair() : i(0), k(0) {}
    AlignedPair(int i, int k) : i(i), k(k) {}
    bool operator<(const AlignedPair &rhs) const { return i < rhs.i || (i == rhs.i && k < rhs.k); }
};

struct Subrange
{
    int i;
    int j;
    int k;
    int l;

    Subrange() : i(0), j(0), k(0), l(0) {}
    Subrange(int i, int j, int k, int l) : i(i), j(j), k(k), l(l) {}
};

enum TracebackType
{
    TB_NONE,
    TB_INS_A,
    TB_INS_B,
    TB_MATCH,
    TB_BIFURCATION,
    TB_CONS_BP,
    NUM_TB_TYPES
};

/////////////////////////////////////////////////////////////////
// class AlignAndFold
/////////////////////////////////////////////////////////////////

class AlignAndFold
{
    const ScoringScheme &scoring;
    
    const SparseMatrix<Real> &paired_A_mask_sparse;
    const SparseMatrix<Real> &paired_B_mask_sparse;
    const SparseMatrix<Real> &aligned_mask_sparse;
    
    const std::vector<Real> &S_paired_A;
    const std::vector<Real> &S_paired_B;
    const std::vector<Real> &S_aligned;

    const int LA;
    const int LB;
    
    const std::vector<Extents> shell;
    
    const bool toggle_verbose;

    // for dynamic programming

    std::vector<SparseMatrix<OuterEntry> > Djlik;
    std::vector<SparseMatrix<OuterEntry *> > Dikjl;
    std::vector<AlignedPair> Sik_last;
    std::vector<Real> Sik;
    std::vector<int> SikTB;    
    
    int cached_i;
    int cached_k;
    
    // for traceback
    
    std::vector<int> mapping_A;
    std::vector<int> mapping_B;
    std::vector<int> aligned_to_A;
    std::vector<int> aligned_to_B;
    
    // helpful macros
    
    bool ValidNW(int i, int k) const;
    bool ValidNW(int i, int j, int k, int l) const;
    int EncodeTraceback(int type, int i, int j) const;
    void DecodeTraceback(int traceback, int &type, int &i, int &j) const;
    
    // initialization of outer matrix

    void InitializeOuterMatrix();
    void TouchOuterMatrix(bool toggle_allocate);
    bool TouchTransition(const int i, const int j, const int k, const int l, 
                         const int ip, const int jp, const int kp, const int lp);
    void TouchOuter(std::vector<std::vector<int> > &Dikjl_row_sizes,
                    std::vector<std::vector<int> > &Djlik_row_sizes,
                    int i, int j, int k, int l, bool toggle_allocate);
    
    // perform dynamic programming

    void FillInnerMatrix(int i, int k);
    OuterEntry *GetOuter(int i, int j, int k, int l);
    void FillOuterMatrix(const int i1, const int k1);
    void FillOuterMatrix();
    void PerformTraceback();
    void PerformTraceback(std::queue<Subrange> &traceback_queue,
                          int i, int j, int k, int l);

    // get alignment
    
    std::string ComputeSecondaryStructure(const std::vector<int> &mapping) const;
    std::string ComputeEditString(const std::vector<int> &aligned_to_A,
                                  const std::vector<int> &aligned_to_B) const;
    
public:

    AlignAndFold(const ScoringScheme &scoring,
                 const bool toggle_alignment_shell,
                 const bool toggle_verbose);
    
    Real DoAlignment();
    std::string FormatString(const std::string &s, const std::string &edit_string, char letter, char gap) const;
    std::string GetEditString() const;
    std::string GetConsensusStructure() const;
    
    std::vector<int> GetAlignedToA() const;
    std::vector<int> GetAlignedToB() const;
    std::vector<int> GetMappingA() const;
    std::vector<int> GetMappingB() const;  
    
    Real ComputeScore(const std::vector<int> &mapping_A,
                      const std::vector<int> &mapping_B,
                      const std::vector<int> &aligned_to_A,
                      const std::vector<int> &aligned_to_B) const;
};

#endif
