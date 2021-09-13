/////////////////////////////////////////////////////////////////
// AlignmentShell.cpp
/////////////////////////////////////////////////////////////////

#include "AlignmentShell.hpp"

/////////////////////////////////////////////////////////////////
// Explicit template instantiations
/////////////////////////////////////////////////////////////////

template std::vector<Extents> AlignmentShell::MakeFullShell<Real>(const SparseMatrix<Real> &mask);
template std::vector<Extents> AlignmentShell::MakeNWShell<Real>(const SparseMatrix<Real> &mask);

/////////////////////////////////////////////////////////////////
// AlignmentShell::MakeFullShell()
//
// Build a NW shell which contains every cell in the Needleman-
// Wunsch matrix.
/////////////////////////////////////////////////////////////////

template<class T>
std::vector<Extents> AlignmentShell::MakeFullShell(const SparseMatrix<T> &mask)
{
    const int LA = mask.GetNumRows()-1;
    const int LB = mask.GetNumCols()-1;
    return std::vector<Extents>(LA+1,Extents(0,LB));
}

/////////////////////////////////////////////////////////////////
// AlignmentShell::MakeNWShell()
//
// Build a list of pairs indicating the first and last cells of
// the Needleman-Wunsch matrix that must be computed in order
// for all entries in the input SparseMatrix to be considered
// in a potential alignment.
/////////////////////////////////////////////////////////////////

template<class T>
std::vector<Extents> AlignmentShell::MakeNWShell(const SparseMatrix<T> &mask)
{
    const int LA = mask.GetNumRows()-1;
    const int LB = mask.GetNumCols()-1;
    
    // initially, all of the left endpoints are on the right edge
    // of the matrix, and all of the right endpoints are on the
    // left edge of the matrix
    
    std::vector<Extents> shell(LA+1,Extents(LB,0));
    
    // for each cell in the sparse matrix
    
    for (int i = 1; i <= LA; i++)
    {
        for (const SparseMatrixEntry<Real> *rp = mask.GetRowBegin(i); rp != mask.GetRowEnd(i); ++rp)
        {
            const int j = rp->column;
            
            Assert(0 < j && j <= LB, "Index out-of-range.");
            
            // if there is an (i,j) aligning pair, then both the
            // (i-1,j-1) and (i,j) cells of the Needleman-Wunsch
            // matrix must be accessible
            
            shell[i-1].begin = std::min(shell[i-1].begin, j-1);
            shell[i-1].end = std::max(shell[i-1].end, j-1);
            shell[i].begin = std::min(shell[i].begin, j);
            shell[i].end = std::max(shell[i].end, j);
        }
    }
    
    // now, perform a backward sweep for the left edge ensuring
    // that the shell left edge is monotonically decreasing
    
    for (int i = LA-1; i >= 0; i--)
        shell[i].begin = std::min(shell[i].begin, shell[i+1].begin);
    
    // similarly, perform a forward sweep for the right edge
    // ensuring that the shell right edge is monotonically
    // increasing
    
    for (int i = 1; i <= LA; i++)
        shell[i].end = std::max(shell[i].end, shell[i-1].end);
    
    // finally, push the left edge to the right edge to ensure 
    // that each row has valid endpoints
    
    for (int i = 0; i <= LA; i++)
        shell[i].begin = std::min(shell[i].begin, shell[i].end);
    
    // the top left and bottom right corners of the matrix must be accessible
    
    shell[0].begin = 0;
    shell[LA].end = LB; 
    
    return shell;  
}

/////////////////////////////////////////////////////////////////
// AlignmentShell::GetShellSize()
//
// Return number of entries in a NW shell.
/////////////////////////////////////////////////////////////////

int AlignmentShell::GetShellSize(const std::vector<Extents> &shell)
{
    int size = 0;
    for (size_t i = 0; i <= shell.size(); i++)
    {
        Assert(shell[i].begin <= shell[i].end, "Improper shell.");
        size += shell[i].end - shell[i].begin + 1;
    }
    return size;
}

/////////////////////////////////////////////////////////////////
// AlignmentShell::TransposeShell()
//
// Transpose shell.
/////////////////////////////////////////////////////////////////

std::vector<Extents> AlignmentShell::TransposeShell(const std::vector<Extents> &shell,
                                                    const int LA, const int LB)
{
    Assert(LA+1 == int(shell.size()), "Dimension mismatch in TransposeShell().");
    
    std::vector<Extents> ret(LB+1, Extents(LA,0));
    for (int r = 0; r < int(shell.size()); r++)
    {
        for (int c = shell[r].begin; c <= shell[r].end; c++)
        {
            ret[c].begin = std::min(ret[c].begin, r);
            ret[c].end = std::max(ret[c].end, r);
        }
    }
    return ret;
}

/////////////////////////////////////////////////////////////////
// AlignmentShell::PrintShells()
//
// Print out NW shell.
/////////////////////////////////////////////////////////////////

void AlignmentShell::PrintShells(const std::vector<Extents> &shell,
                                 const int LA, const int LB)
{
    for (int i = 0; i <= LA; i++)
    {
        for (int j = 0; j <= LB; j++)
        {
            bool inNW = (shell[i].begin <= j && j <= shell[i].end);
            std::cerr << (inNW ? '#' : '.');
        }
        std::cerr << std::endl;
    }
    
    std::cerr << std::endl
              << "#  in NW shell" << std::endl
              << ".  not in NW shell" << std::endl
              << std::endl;
}
