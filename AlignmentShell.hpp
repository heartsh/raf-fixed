/////////////////////////////////////////////////////////////////
// AlignmentShell.hpp
/////////////////////////////////////////////////////////////////

#ifndef ALIGNMENTSHELL_HPP
#define ALIGNMENTSHELL_HPP

#include "Utilities.hpp"
#include "SparseMatrix.hpp"
#include "Options.hpp"

struct Extents
{
    int begin, end;
    
    Extents(int begin, int end) : begin(begin), end(end) {}
};
    

/////////////////////////////////////////////////////////////////
// class AlignmentShell
/////////////////////////////////////////////////////////////////

class AlignmentShell
{
public:
    template<class T> static std::vector<Extents> MakeFullShell(const SparseMatrix<T> &mask);
    template<class T> static std::vector<Extents> MakeNWShell(const SparseMatrix<T> &mask);
    
    static std::vector<Extents> TransposeShell(const std::vector<Extents> &shell,
                                               const int LA, const int LB);
    static int GetShellSize(const std::vector<Extents> &shell);
    static void PrintShells(const std::vector<Extents> &shell, 
                            const int LA, const int LB);
};

#endif
