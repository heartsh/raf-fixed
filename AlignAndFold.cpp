//////////////////////////////////////////////////////////////////
// AlignAndFold.cpp
//
// This class performs sparse dynamic programming in order to
// compute the optimal alignment of two sequences.
//
// I.  Recurrences
//
// For 0 <= i <= j <= LA and 0 <= k <= l <= LB, 
// 
//   S(i,j;k,l)
//
//   = score for aligning A[i+1...j] with B[k+1...l]
//
//         [ 0                                    if i == j && k == l
//         [
//         [ Unpaired nucleotides:
//         [   S(i,j-1;k,l)       + SCORE( -      , A[j] ; -      , -    )
//         [   S(i,j;k,l-1)       + SCORE( -      , -    ; -      , B[l] )
//         [   S(i,j-1;k,l-1)     + SCORE( -      , A[j] ; -      , B[l] )
//         [ 
//         [ Base-pairing indels:
//         [   S(i+1,j-1;k,l)     + SCORE( A[i+1] , A[j] ; -      , -    )
//         [   S(i,j;k+1,l-1)     + SCORE( -      , -    ; B[k+1] , B[l] )
//         [ 
//   = max [ Non-conserved base-pairings:
//         [   S(i+1,j;k+1,l-1)   + SCORE( A[i+1] , -    ; B[k+1] , B[l] )
//         [   S(i,j-1;k+1,l-1)   + SCORE( -      , A[j] ; B[k+1] , B[l] )
//         [   S(i+1,j-1;k+1,l)   + SCORE( A[i+1] , A[j] ; B[k+1] , -    )
//         [   S(i+1,j-1;k,l-1)   + SCORE( A[i+1] , A[j] ; -      , B[l] )
//         [
//         [ Conserved base-pairings:
//         [   S(i+1,j-1;k+1,l-1) + SCORE(A[i+1] , A[j] ; B[k+1] , B[l] )
//         [
//         [ Bifurcation:
//         [   S(i,j';k,l') + S(j',j;l',l)      where i <= j' <= j
//         [                                          k <= l' <= l
//
// These recurrences are very similar to the ones used in the
// FOLDALIGN program.  We note the following:
//
// (1) Indices in the FOLDALIGN recurrences refer to letters in
//     the RNA sequences.  Here, indices refer to positions between
//     letters (ranging from 0 to LA or 0 to LB).
//
// (2) FOLDALIGN uses the following extra recurrences for unpaired nucleotides:
//
//         [ Unpaired nucleotides:
//         [
//         [   S(i+1,j;k,l)   + SCORE( A[i+1] , -    ; -      , -    )
//         [   S(i,j;k+1,l)   + SCORE( -      , -    ; B[k+1] , -    )
//         [   S(i+1,j;k+1,l) + SCORE( A[i+1] , -    ; B[k+1] , -    )
//         [
//         [   S(i+1,j;k,l-1) + SCORE( A[i+1] , -    ; -      , B[l] )
//         [   S(i,j-1;k+1,l) + SCORE( -      , A[j] ; B[k+1] , -    )
//
//     The first three recurrences allow for non base-pairing
//     nucleotides on the left sides of A[i+1...j] and B[k+1...l].
//     However, these cases are already accounted for by the three
//     existing recurrences for unpaired nucleotides in conjunction
//     with the bifurcation rule.  For example, the first case (an
//     insertion of the leftmost character in A[i+1....j]) can be
//     represented as:
//
//                           +---- S(i,j;k,l) ----+
//                           |                    |
//                           V                    V
//                     S(i,i+1;k,k)    +     S(i+1,j;k,l)          (bifurcation rule)
//                           |
//                    +------+----------------+
//                    |                       |
//                    V                       V
//               S(i,i;k,k)  +  SCORE( - , A[i+1] ; - , - )      (right insertion rule)
//                    |                       |
//                    V                       V
//                    0         SCORE( A[i+1] , - ; - , - )
//
//
//     Similarly, the last two recurrences allow for non base-pairing
//     nucleotides on opposite sides of A[i+1...j] and B[k+1...l].
//     Again, these cases should be scored the same as an insertion on
//     the left followed by an insertion on the right.
//
//     Thus, we can ignore both cases altogether.
//
// II. Exploiting base-pairing sparsity
// 
// The LocARNA global alignment algorithm disallows base-pair indels
// as well as non-conserved base-pairings; that is, LocARNA deals with
// a simplified version of the above recurrences:
//
// For 0 <= i <= j <= LA and 0 <= k <= l <= LB,
// 
//   S(i,j;k,l)
//       
//         [ 0                                    if i == j && k == l
//         [
//         [ Unpaired nucleotides:
//         [   S(i,j-1;k,l)       + SCORE( -      , A[j] ; -      , -    )
//         [   S(i,j;k,l-1)       + SCORE( -      , -    ; -      , B[l] )
//         [   S(i,j-1;k,l-1)     + SCORE( -      , A[j] ; -      , B[l] )
//         [ 
//   = max [ Conserved base-pairings:
//         [   S(i+1,j-1;k+1,l-1) + SCORE( A[i+1] , A[j] ; B[k+1] , B[l] )
//         [
//         [ Bifurcation:
//         [   S(i,j';k,l') + S(j',j;l',l)      where i <= j' <= j
//         [                                          k <= l' <= l
//     
// To achieve a speed up, LocARNA rewrites these recurrences by noting
// that a bifurcation is only necessary whenever there is some type of
// new conserved base-pairing being created adjacent to an existing
// structure.  For instance, in the following alignment,
//     
//                i         j'j'j'    j
//                |         | | |     |
//                 A C G G U G U G U C
//                 ( . . . ) . . ( . )
//                 G C G G C   U G U C
//                |          |  |     |
//                k          l' l'    l
//     
// a bifurcation is needed in order to allow the two conserved
// base-pairings.  The location of the bifurcation point (j',l'),
// however, is flexible in the sense that any of the three choices of
// j' shown above and any of the two choices of l' shown above could
// give rise to the same structure.  LocARNA removes this ambiguity by
// forcing the right side of the bifurcation to always be a conserved
// base-pairing:
//
// For 0 <= i <= j <= LA and 0 <= k <= l <= LB,
// 
//   S(i,j;k,l)
//       
//         [ 0                                    if i == j && k == l
//         [
//         [ Unpaired nucleotides:
//         [   S(i,j-1;k,l)       + SCORE( -      , A[j] ; -      , -    )
//         [   S(i,j;k,l-1)       + SCORE( -      , -    ; -      , B[l] )
//   = max [   S(i,j-1;k,l-1)     + SCORE( -      , A[j] ; -      , B[l] )
//         [ 
//         [ Bifurcation with conserved base-pairings:
//         [   S(i,j';k,l') + D(j',j;l',l)        where i <= j' <= j
//         [                                            k <= l' <= l
//
// where
//
//   D(i,j;k,l) = S(i+1,j-1;k+1,l-1) + SCORE( A[i+1] , A[j] ; B[k+1] , B[l] ).
//
// By writing it this way, it is clear that the only (j',l') that must
// be considered in the bifurcation step should be limited to only
// those for which (A[j'+1],A[j]) and (B[l'+1],B[l]) are valid
// base-pairs in A and B, respectively.  Even more importantly, if we
// know the values of all D(*,*;*.*); then, we can compute S(i,*;k,*)
// efficiently for any i and k.  In particular, this means that we can
// implement the recurrences above by storing only the values of the D
// matrix explicitly; given the D matrix, any S values needed can be
// computed efficiently.
//
// III. Exploiting alignment sparsity
//
// The second type of sparsity we will exploit in this algorithm is
// the sparsity that arises from restricted the set of positions
// considered to be alignable.  For two sequences A and B, the
// Needleman-Wunsch (NW) matrix for A and B is an (LA+1) x (LB+1)
// matrix used for performing global sequence alignment; let S denote
// a subset of the NW matrix which we call the "NW shell"; we will
// only evaluate S(i,j;k,l) if (i,k) in S and (j,l) in S.  All other
// S(i,j;k,l) are assumed to be -INF.  By restricting our evaluation
// of the S matrix in this way, we can substantially reduce the amount
// of computation that must be performed.
// 
/////////////////////////////////////////////////////////////////

#include "AlignAndFold.hpp"

const Real INF = Real(1e20);

/////////////////////////////////////////////////////////////////
// UPDATE_MAX_TRACEBACK()
//
// Macro for updating a score/traceback pointer which does not
// evaluate t unless an update is needed.  
/////////////////////////////////////////////////////////////////

#define UPDATE_MAX_TRACEBACK(bs,bt,s,t) \
    do                                  \
    {                                   \
        Real work(s);                   \
        if ((work)>(bs))                \
        {                               \
            (bs)=(work);                \
            (bt)=(t);                   \
        }                               \
    }                                   \
    while (0)

/////////////////////////////////////////////////////////////////
// UPDATE_MAX()
//
// Macro for updating a score.
/////////////////////////////////////////////////////////////////

#define UPDATE_MAX(bs,s)                \
    do                                  \
    {                                   \
        Real work(s);                   \
        if ((work)>(bs))                \
        {                               \
            (bs)=(work);                \
        }                               \
    }                                   \
    while (0)

/////////////////////////////////////////////////////////////////
// ValidNW()
// AlignAndFold::EncodeTraceback()
// AlignAndFold::DecodeTraceback()
//
// Helpful macros.
/////////////////////////////////////////////////////////////////

inline bool AlignAndFold::ValidNW(int i, int k) const
{
    return (shell[i].begin <= k && k <= shell[i].end);
}

inline bool AlignAndFold::ValidNW(int i, int j, int k, int l) const
{
    return i <= j && k <= l && ValidNW(i,k) && ValidNW(j,l);
}

inline int AlignAndFold::EncodeTraceback(int type, int i, int j) const
{
    Assert(TB_NONE <= type && type <= TB_BIFURCATION, "Invalid arguments for traceback.");
    Assert(i >= 0 && j >= 0, "Invalid arguments for traceback.");
    return (i*(LB+1)+j) * NUM_TB_TYPES + type;
}

inline void AlignAndFold::DecodeTraceback(int traceback, int &type, int &i, int &j) const
{
    Assert(traceback >= 0, "Invalid arguments for traceback.");
    
    type = traceback % NUM_TB_TYPES;
    traceback /= NUM_TB_TYPES;
    
    j = traceback % (LB+1);
    i = traceback / (LB+1);
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::AlignAndFold()
//
// Constructor.
/////////////////////////////////////////////////////////////////

AlignAndFold::AlignAndFold(const ScoringScheme &scoring,
                           const bool toggle_alignment_shell,
                           const bool toggle_verbose) :

    scoring(scoring),
    
    paired_A_mask_sparse(scoring.GetFoldingAMaskSparse()),
    paired_B_mask_sparse(scoring.GetFoldingBMaskSparse()),
    aligned_mask_sparse(scoring.GetAlignmentMaskSparse()),

    S_paired_A(scoring.GetFoldingAMatrix()),
    S_paired_B(scoring.GetFoldingBMatrix()),
    S_aligned(scoring.GetAlignmentMatrix()),

    LA(paired_A_mask_sparse.GetNumRows()-1),
    LB(paired_B_mask_sparse.GetNumRows()-1),
    
    shell(toggle_alignment_shell ?
          AlignmentShell::MakeNWShell<Real>(aligned_mask_sparse) :
          AlignmentShell::MakeFullShell<Real>(aligned_mask_sparse)),
    
    toggle_verbose(toggle_verbose),
    
    Sik((LA+1)*(LB+1)), 
    SikTB((LA+1)*(LB+1)),
    
    cached_i(-1), 
    cached_k(-1)
{
    InitializeOuterMatrix();
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::InitializeOuterMatrix()
//
// Initialize the outer dynamic programming matrix, D.  Internally,
// D(i,j;k,l) is represented as a matrix (over indices j and l) of
// sparse matrices (over indices i and k).  This representation
// (called Djlik) allows us to iterate over different i's and k's for
// a fixed choice of j and l (as needed in the dynamic programming
// recursions).
//
// To allow iteration over different j's and l's for a fixed
// choice of i and k, a second matrix of pointers to the first
// matrix is kept (called Dikjl).
/////////////////////////////////////////////////////////////////

void AlignAndFold::InitializeOuterMatrix()
{
    if (toggle_verbose) WriteProgressMessage("Initializing outer DP matrix...");

    // The Sik_last array is used in order to speed up the initial
    // fill of the D matrix.  Although the DP matrix is technically
    // defined for all inner(i,j;k,l) such that (i,k) in S and (j,l)
    // in S, not all of these entries are actually needed in order to
    // compute the D matrix.  The Sik_last[i*(LB+1)+j].{j,l} values
    // keep track of the largest j and l values that are actually
    // needed for filling the D matrix.
    //
    // By default, the last j and l values for i=0, k=0 are set to LA
    // and LB, respectively; this is needed for proper global
    // alignment of the two sequences.
    
    Sik_last.clear();
    Sik_last.resize((LA+1)*(LB+1), AlignedPair(-1,-1));
    Sik_last[0*(LB+1)+0] = AlignedPair(LA,LB);

    // In order to ensure that allocation of the D matrix is done in
    // linear time, we allocate the D matrix in two passes.  In the
    // first pass, we determine the dimensions of each of the sparse
    // matrices in Djlik/Dikjl, and then we must allocate memory for
    // these matrices.  In the second pass, we initialize the entries
    // of the Djlik matrix; we also create pointers from the Dikjl
    // matrix to entries of the Djlik matrix.
    
    TouchOuterMatrix(true);
    TouchOuterMatrix(false);    
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::TouchOuterMatrix()
//
// Touch each of the entries of the dynamic programming matrices.  If
// toggle_allocate is true, the dimensions of the sparse matrices are
// determined and the space for the Djlik and Dijkl matrices is
// allocated.  If toggle_allocate is false, then the Djlik entries are
// initialized and pointers from Dijkl to Djlik entries are set.
/////////////////////////////////////////////////////////////////

void AlignAndFold::TouchOuterMatrix(bool toggle_allocate)
{
    // Transpose the pairing matrices in order to be able to access
    // lists of base-pairs indexed by the rightmost element (rather
    // than leftmost element).  Also, transpose the alignment shell
    // for access by columns rather than by rows.
    
    const SparseMatrix<Real> A_basepairs(paired_A_mask_sparse, SparseMatrix<Real>::TRANSPOSE);
    const SparseMatrix<Real> B_basepairs(paired_B_mask_sparse, SparseMatrix<Real>::TRANSPOSE);
    const std::vector<Extents> shell_T = AlignmentShell::TransposeShell(shell, LA, LB);

    // These variables keep track of the number of elements in each
    // row of the sparse matrices used to represent the D matrix.  In
    // particular, Dikjl[i*(LB+1)+k][j] is the size of the jth row of
    // D(i,*;k,*), and Djlik[j*(LB+1)+l][i] is the size of the ith row
    // of D(*,j;*,l).
    
    std::vector<std::vector<int> > Dikjl_row_sizes((LA+1)*(LB+1));
    std::vector<std::vector<int> > Djlik_row_sizes((LA+1)*(LB+1));

    const SparseMatrixEntry<Real> *pA;
    const SparseMatrixEntry<Real> *pB;    
    
    int ct = 0;
    
    // Iterate through positions (j,l) in S.

    for (int j = 0; j <= LA; j++)
    {
        for (int l = shell[j].begin; l <= shell[j].end; l++)
        {
            // For the current (j,l), build a list of all (i,k)'s such
            // that D(i,j;k,l) must be computed and stored.
            
            std::set<AlignedPair> active_iks;

            // Include (i,k)'s that correspond to potential conserved
            // base-pairings:
            // 
            // Iterate through positions i such that (i+1,j) in BP_A
            //             and positions k such that (k+1,l) in BP_B
            //             and (i,k) in S and (i+1,k) in S.
            
            for (pA = A_basepairs.GetRowBegin(j); pA != A_basepairs.GetRowEnd(j); ++pA)
            {
                const int i = pA->column-1;
                for (pB = B_basepairs.GetRowBegin(l); pB != B_basepairs.GetRowEnd(l); ++pB)
                {
                    const int k = pB->column-1;

                    // Check if this transition is valid.  If so, add
                    // (i,k) to the list of active (i,k)'s.
                    
                    if (TouchTransition(i,j,k,l,i+1,j-1,k+1,l-1))
                        active_iks.insert(AlignedPair(i,k));
                }
            }
            
            // Now, given the list of all active (i,k)'s, we'll
            // "touch" each of the D(i,j,k,l) entries.  During the
            // allocation phase, "touch" means that we will keep track
            // only of the amount of memory needed by these cells;
            // during the initialization phase, "touch" means that we
            // can actually populate the entries with values.
            
            for (std::set<AlignedPair>::const_iterator iter = active_iks.begin(); iter != active_iks.end(); ++iter)
                TouchOuter(Dikjl_row_sizes, Djlik_row_sizes, iter->i, j, iter->k, l, toggle_allocate);
            
            ct += int(active_iks.size());
        }
    }
    
    // Allocate memory for D matrix.
    
    if (toggle_allocate)
    {
        // At this point, the Dikjl row_sizes and Djlik row_sizes matrices should be
        // complete, so we will know exactly how much memory is needed for the D
        // matrix.  We can just allocate it directly.
        
        Dikjl = std::vector<SparseMatrix<OuterEntry *> >((LA+1)*(LB+1));
        Djlik = std::vector<SparseMatrix<OuterEntry> >((LA+1)*(LB+1));
        
        for (int r = 0; r <= LA; r++)
        {
            for (int c = 0; c <= LB; c++)
            {
                const int rc = r*(LB+1)+c;
                if (Dikjl_row_sizes[rc].size() > 0) Dikjl[rc] = SparseMatrix<OuterEntry *>(Dikjl_row_sizes[rc], LA+1, LB+1, NULL);
                if (Djlik_row_sizes[rc].size() > 0) Djlik[rc] = SparseMatrix<OuterEntry>(Djlik_row_sizes[rc], LA+1, LB+1, OuterEntry(0,0));
            }
        }
    }
    else
    {
        if (toggle_verbose)
            std::cerr << ct << " total entries in outer DP matrix initialized." << std::endl;
    }
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::TouchTransition()
// 
// Mark that a particular D(i,j;k,l) cell depends on S(ip,jp;kp,lp) if
// all coordinates are valid.  Return true if valid.
/////////////////////////////////////////////////////////////////

bool AlignAndFold::TouchTransition(const int i, const int j, const int k, const int l, 
                                   const int ip, const int jp, const int kp, const int lp)
{
    // Check that both D(i,j,k,l) and S(ip,jp;kp,lp) refer to valid
    // coordinates.
    
    if (ValidNW(i,j,k,l) && ValidNW(ip,jp,kp,lp))
    {
        const int ipkp = ip*(LB+1)+kp;

        // Update the table of last (j,l) coordinates needed for each
        // (i,k).
        
        Sik_last[ipkp].j = std::max(Sik_last[ipkp].j, jp);
        Sik_last[ipkp].l = std::max(Sik_last[ipkp].l, lp);
        
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::TouchOuter()
//
// Touch a single entry of the dynamic programming matrix.
/////////////////////////////////////////////////////////////////

void AlignAndFold::TouchOuter(std::vector<std::vector<int> > &Dikjl_row_sizes,
                              std::vector<std::vector<int> > &Djlik_row_sizes,
                              int i, int j, int k, int l, bool toggle_allocate)
{
    Assert(i <= j && k <= l, "Invalid indices (i=%d,j=%d,k=%d,l=%d).", i, j, k, l);
    
    const int ik = i*(LB+1)+k;
    const int jl = j*(LB+1)+l;
    
    // allocate memory if Dikjl(i,*;k,*) or Djlik(*,j;*,l) has never been touched before
    
    if (Dikjl_row_sizes[ik].size() == 0) Dikjl_row_sizes[ik].resize(LA+1);
    if (Djlik_row_sizes[jl].size() == 0) Djlik_row_sizes[jl].resize(LA+1);
    
    // fill entries if not in memory allocation mode
    
    if (!toggle_allocate)
    {
        // The next two lines are a little tricky.  In the first line,
        // we access the jth row of D(i,*;k,*), and in the second
        // line, we access the ith row of D(*,j;*,l).  The tricky part
        // is to notice that Dikjl_row_sizes[ik][j] and
        // Djlik_row_size[jl][i] gives the "l" and k" indices of the
        // elements that must be filled.  The fact that
        // Dikjl_row_sizes[ik][j] and Djlik_row_size[jl][i] both get
        // incremented during each call of TouchOuter() ensures uniqueness
        // of every element.
        
        SparseMatrixEntry<OuterEntry *> &Dikjl_entry = Dikjl[ik].GetRowBegin(j)[Dikjl_row_sizes[ik][j]];
        SparseMatrixEntry<OuterEntry> &Djlik_entry = Djlik[jl].GetRowBegin(i)[Djlik_row_sizes[jl][i]];
        
        Dikjl_entry.column = l;
        Dikjl_entry.value = &Djlik_entry.value;
        Djlik_entry.column = k;
        Djlik_entry.value.score = -INF;
    }
    
    // Keep track of current position in each row.
    
    ++(Dikjl_row_sizes[ik][j]);
    ++(Djlik_row_sizes[jl][i]);
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::FillInnerMatrix()
//
// Compute Sik(i,*;k,*) matrix. 
/////////////////////////////////////////////////////////////////

void AlignAndFold::FillInnerMatrix(const int i, const int k)
{
    Assert(ValidNW(i,k), "Invalid aligned nucleotide pair.");
    
    // check cache to see if this is already computed
    
    if (cached_i == i && cached_k == k) return;
    cached_i = i;
    cached_k = k;
    
    // iterate through positions j >= i and l >= k such that (j,l) in S
    
    const int ik = i*(LB+1)+k;
    for (int j = i; j <= Sik_last[ik].j; j++)
    {
        for (int l = std::max(k, shell[j].begin); l <= Sik_last[ik].l; l++)
        {
            const int jl = j*(LB+1)+l;
            const int j1l = (j-1)*(LB+1)+l;
            const int jl1 = j*(LB+1)+l-1;
            const int j1l1 = (j-1)*(LB+1)+l-1;

            Sik[jl] = -INF;

            // initialization

            if (j == i && l == k) 
                UPDATE_MAX_TRACEBACK(Sik[jl], SikTB[jl], 0, EncodeTraceback(TB_NONE,j,l));
            
            // unpaired nucleotides

            if (j > i && ValidNW(j-1,l))
                UPDATE_MAX_TRACEBACK(Sik[jl], SikTB[jl], Sik[j1l], EncodeTraceback(TB_INS_A,j-1,l));

            if (l > k && ValidNW(j,l-1))
                UPDATE_MAX_TRACEBACK(Sik[jl], SikTB[jl], Sik[jl1], EncodeTraceback(TB_INS_B,j,l-1));
            
            if (j > i && l > k && ValidNW(j-1,l-1))
                UPDATE_MAX_TRACEBACK(Sik[jl], SikTB[jl], Sik[j1l1] + S_aligned[jl], EncodeTraceback(TB_MATCH,j-1,l-1));
            
            // bifurcation
            
            const SparseMatrix<OuterEntry> &Djl = Djlik[jl];
            if (Djl.GetNumRows() > 0)
            {
                Real Sikjl = Sik[jl];
                int SikjlTB = SikTB[jl];
                
                // iterate through valid (j',l')
                
                for (int jp = j; jp >= i; jp--)
                {
                    for (const SparseMatrixEntry<OuterEntry> *rp = Djl.GetRowRBegin(jp); rp != Djl.GetRowREnd(jp); --rp)
                    {
                        const int lp = rp->column;
                        if (lp < k) break;
                        
                        Assert(i <= jp && jp <= j && k <= lp && lp <= l, "Invalid indices (i=%d,jp=%d,j=%d,k=%d,lp=%d,l=%d)", i, jp, j, k, lp, l);
                        UPDATE_MAX_TRACEBACK(Sikjl, SikjlTB, Sik[jp*(LB+1)+lp] + rp->value.score, EncodeTraceback(TB_BIFURCATION,jp,lp));
                    }
                }
                
                Sik[jl] = Sikjl;
                SikTB[jl] = SikjlTB;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::GetOuter()
//
// Get pointer to entry of the D matrix.
/////////////////////////////////////////////////////////////////

inline OuterEntry *AlignAndFold::GetOuter(const int i, const int j, const int k, const int l)
{
    if (0 <= i && i <= j && j <= LA && 0 <= k && k <= l && l <= LB)
    {
        const int ik = i*(LB+1)+k;
        if (Dikjl[ik].GetNumRows() <= 0) return NULL;
        return Dikjl[ik](j,l);
    }
    return NULL;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::FillOuterMatrix()
//
// Use Si1k1 matrix to fill D(i,*;k,*).
/////////////////////////////////////////////////////////////////

void AlignAndFold::FillOuterMatrix(const int i1, const int k1)
{
    const int i = i1 - 1;
    const int k = k1 - 1;
    if (i < 0 || k < 0) return;
    const int ik = i*(LB+1)+k;
    SparseMatrix<OuterEntry *> &Dik = Dikjl[ik];

    for (int j = 0; j < Dik.GetNumRows(); j++)
    {
        for (SparseMatrixEntry<OuterEntry *> *p = Dik.GetRowBegin(j); p != Dik.GetRowEnd(j); ++p)
        {
            const int l = p->column;
            
            UPDATE_MAX_TRACEBACK(p->value->score, p->value->traceback,
                                 Sik[(j-1)*(LB+1)+(l-1)]
                                 + S_paired_A[(i+1)*(LA+1)+j]
                                 + S_paired_B[(k+1)*(LB+1)+l]
                                 + S_aligned[(i+1)*(LB+1)+(k+1)]
                                 + S_aligned[j*(LB+1)+l],
                                 TB_CONS_BP);
        }
    }
}

 /*
void AlignAndFold::FillOuterMatrix(const int i, const int k)
{
    const int ik = i*(LB+1)+k;

    for (int j = i; j <= Sik_last[ik].j; j++)
    {
        for (int l = std::max(k, shell[j].begin); l <= Sik_last[ik].l; l++)
        {
            const int jl = j*(LB+1)+l;
            const int ij1 = i*(LA+1)+j+1;
            const int kl1 = k*(LB+1)+l+1;
            const int j1l1 = (j+1)*(LB+1)+l+1;

            OuterEntry *Dentry;

            // conserved base-pairings

            if ((Dentry = GetOuter(i-1,j+1,k-1,l+1)) != NULL)
            {
                UPDATE_MAX_TRACEBACK(Dentry->score, Dentry->traceback, Sik[jl] + S_paired_A[ij1] + S_paired_B[kl1] + S_aligned[ik] + S_aligned[j1l1], TB_CONS_BP);
            }
        }
    }
    }*/


/////////////////////////////////////////////////////////////////
// AlignAndFold::FillOuterMatrix()
//
// Compute dynamic programming matrix.
/////////////////////////////////////////////////////////////////

void AlignAndFold::FillOuterMatrix()
{
    int num_needed = 0;
    int num_done = 0;      
    int percent_complete = -1;
    
    // count number of (i,k) cells needed
    
    if (toggle_verbose)
    {
        for (int i = 0; i <= LA; i++)
            for (int k = shell[i].begin; k <= shell[i].end; k++)
                if (Sik_last[i*(LB+1)+k].j >= 0) num_needed++;
    }
    
    // iterate through candidate aligned positions (i,k)
    
    for (int i = LA; i >= 0; i--)
    {
        for (int k = shell[i].end; k >= shell[i].begin; k--)
        {
            if (Sik_last[i*(LB+1)+k].j == -1) continue;
            
            if (toggle_verbose)
            {
                int new_percent_complete = 100 * num_done / num_needed;
                if (new_percent_complete != percent_complete)
                {
                    percent_complete = new_percent_complete;
                    WriteProgressMessage(SPrintF("Performing pairwise alignment...%d%% complete.", percent_complete));
                }
            }
            
            FillInnerMatrix(i,k);
            FillOuterMatrix(i,k);
            
            num_done++;
        }
    }
    
    if (toggle_verbose)
        WriteProgressMessage("");
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::PerformTraceback()
//
// Perform traceback.
/////////////////////////////////////////////////////////////////

void AlignAndFold::PerformTraceback()
{
    mapping_A = std::vector<int>(LA+1,0);
    mapping_B = std::vector<int>(LB+1,0);
    aligned_to_A = std::vector<int>(LA+1,0);
    aligned_to_B = std::vector<int>(LB+1,0);
    
    if (toggle_verbose)
        WriteProgressMessage("Performing traceback...");
    
    // queue contains the current ((i,j),(k,l)) subrange that needs to be processed
    
    std::queue<Subrange> traceback_queue;
    traceback_queue.push(Subrange(0,LA,0,LB));
    
    while (!traceback_queue.empty())
    {
        Subrange s = traceback_queue.front();
        traceback_queue.pop();
        PerformTraceback(traceback_queue, s.i, s.j, s.k, s.l);
    }
    
    if (toggle_verbose)
        WriteProgressMessage("");
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::PerformTraceback()
//
// Perform traceback for the current (i,k), starting at (j,l).
/////////////////////////////////////////////////////////////////

void AlignAndFold::PerformTraceback(std::queue<Subrange> &traceback_queue,
                                    int i, int j, int k, int l)
{
    Assert(ValidNW(i,j,k,l), "Invalid indices (i=%d,j=%d,k=%d,l=%d).", i, j, k, l);
    FillInnerMatrix(i,k);

    int type = 0;
    int nj = 0;
    int nl = 0;

    while (true)
    {
        // process a traceback pointer
        
        DecodeTraceback(SikTB[j*(LB+1)+l], type, nj, nl);
        
        if (type == TB_NONE)
        {
            break;
        }
        else if (type == TB_MATCH)
        {
            aligned_to_A[j] = l;
            aligned_to_B[l] = j;
        }
        else if (type == TB_BIFURCATION)
        {
            Assert(GetOuter(nj,j,nl,l) != NULL, "Bad traceback.");
            int traceback = GetOuter(nj,j,nl,l)->traceback;
            switch (traceback)
            {
                case TB_CONS_BP: 
                    aligned_to_A[j] = l;
                    aligned_to_B[l] = j;
                    aligned_to_A[nj+1] = nl+1;
                    aligned_to_B[nl+1] = nj+1;
                    mapping_A[nj+1] = j;
                    mapping_A[j] = nj+1;
                    mapping_B[nl+1] = l;
                    mapping_B[l] = nl+1;
                    traceback_queue.push(Subrange(nj+1,j-1,nl+1,l-1));
                    break;
                default:
                    Assert(false, "Should not get here.");
            }
        }
        
        j = nj;
        l = nl;
    }
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::ComputeSecondaryStructure()
//
// Compute secondary structure for a given mapping.
/////////////////////////////////////////////////////////////////

std::string AlignAndFold::ComputeSecondaryStructure(const std::vector<int> &mapping) const
{
    std::string ret = "@";
    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] == 0) ret.push_back('.');
        else if (mapping[i] > i) ret.push_back('(');
        else if (mapping[i] < i) ret.push_back(')');
        else
        {
            Assert(false, "Bad secondary structure mapping.");
        }
    }
    return ret;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::ComputeEditString()
//
// Return a sequence containing {A,B,M}.
/////////////////////////////////////////////////////////////////

std::string AlignAndFold::ComputeEditString(const std::vector<int> &aligned_to_A, const std::vector<int> &aligned_to_B) const {
    size_t i = 1;
    size_t j = 1;
    std::string ret = "@";
    
    while (i < aligned_to_A.size() || j < aligned_to_B.size())
    {
        if (i < aligned_to_A.size() && aligned_to_A[i] == 0)
        {
            i++;
            ret.push_back('A');
        }
        else if (j < aligned_to_B.size() && aligned_to_B[j] == 0)
        {
            j++;
            ret.push_back('B');
        }
        else
        {
            Assert(aligned_to_A[i] == int(j) && aligned_to_B[j] == int(i), "Bad alignment.");
            i++;
            j++;
            ret.push_back('M');
        }
    }
    return ret;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::GetEditString()
//
// Return edit string.
/////////////////////////////////////////////////////////////////

std::string AlignAndFold::GetEditString() const
{
    return ComputeEditString(aligned_to_A, aligned_to_B);
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::FormatString()
//
// Apply edit string to a given string s.
/////////////////////////////////////////////////////////////////

std::string AlignAndFold::FormatString(const std::string &s, const std::string &edit_string, char letter, char gap) const
{
    std::string ret = "@";
    size_t j = 0;
    
    for (size_t i = 1; i < edit_string.size(); i++)
    {
        if (edit_string[i] == letter || edit_string[i] == 'M')
        {
            ret.push_back(s[++j]);
            Assert(j < s.length(), "Invalid character.");
        }
        else
        {
            ret.push_back(gap);
        }
    }
    
    return ret;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::ComputeScore()
//
// Compute score for a fixed structure prediction and
// alignment.
/////////////////////////////////////////////////////////////////

Real AlignAndFold::ComputeScore(const std::vector<int> &mapping_A,
                                const std::vector<int> &mapping_B,
                                const std::vector<int> &aligned_to_A,
                                const std::vector<int> &aligned_to_B) const
{
    
    if (int(mapping_A.size()) != LA+1 ||
        int(mapping_B.size()) != LB+1 ||
        int(aligned_to_A.size()) != LA+1 ||
        int(aligned_to_B.size()) != LB+1)
        Error("Dimension mismatch.");
    
    Real score_pairing_A = 0;
    Real score_pairing_B = 0;
    Real score_alignment = 0;
    
    for (int i = 1; i <= LA; i++)
        score_pairing_A += (mapping_A[i] > i ? S_paired_A[i*(LA+1)+mapping_A[i]] : 0);
    for (int i = 1; i <= LB; i++)
        score_pairing_B += (mapping_B[i] > i ? S_paired_B[i*(LB+1)+mapping_B[i]] : 0);
    for (int i = 1; i <= LA; i++)
        score_alignment += (aligned_to_A[i] != 0 ? S_aligned[i*(LB+1)+aligned_to_A[i]] : 0);

    return
        score_pairing_A + score_pairing_B + score_alignment +
        scoring.GetOffset();
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::DoAlignment()
//
// Perform alignment dynamic programming.
/////////////////////////////////////////////////////////////////

Real AlignAndFold::DoAlignment()
{
    FillOuterMatrix();
    FillInnerMatrix(0,0);
    Real score = Sik[LA*(LB+1)+LB];
    PerformTraceback();
    
    return score + scoring.GetOffset();
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::GetConsensusStructure()
//
// Retrieve predicted consensus structure.
/////////////////////////////////////////////////////////////////

std::string AlignAndFold::GetConsensusStructure() const
{
    const std::string edit_string = ComputeEditString(aligned_to_A, aligned_to_B);
    const std::string struct_A = FormatString(ComputeSecondaryStructure(mapping_A), edit_string, 'A', '-');
    const std::string struct_B = FormatString(ComputeSecondaryStructure(mapping_B), edit_string, 'B', '-');
    
    std::string consensus = "@";
    
    for (size_t i = 1; i < struct_A.length(); i++)
    {
        const char ch =
            (struct_A[i] == struct_B[i]) ? struct_A[i] :
            (struct_A[i] == '-') ? struct_B[i] :
            (struct_B[i] == '-') ? struct_A[i] : '.';
        consensus.push_back(ch);
    }
    
    return consensus;
}

/////////////////////////////////////////////////////////////////
// AlignAndFold::GetAlignedToA()
// AlignAndFold::GetAlignedToB()
// AlignAndFold::GetMappingA()
// AlignAndFold::GetMappingB()
//
// Retrieve predicted features.
/////////////////////////////////////////////////////////////////

std::vector<int> AlignAndFold::GetAlignedToA() const { return aligned_to_A; }
std::vector<int> AlignAndFold::GetAlignedToB() const { return aligned_to_B; }
std::vector<int> AlignAndFold::GetMappingA() const { return mapping_A; }
std::vector<int> AlignAndFold::GetMappingB() const { return mapping_B; }
