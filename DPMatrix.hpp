/////////////////////////////////////////////////////////////////
// DPMatrix.hpp
/////////////////////////////////////////////////////////////////

#ifndef DPMATRIX_HPP
#define DPMATRIX_HPP

#include "Utilities.hpp"
#include "SparseMatrix.hpp"
#include "AlignmentShell.hpp"
#include "Options.hpp"

struct DPEntry
{
    Real score;
    int traceback;
    
    DPEntry() : score(0), traceback(0) {}
    DPEntry(Real score) : score(score), traceback(0) {}
    DPEntry(Real score, int traceback) : score(score), traceback(traceback) {}
};

struct LastPair
{
    int j;
    int l;
    
    LastPair(int j, int l) : j(j), l(l) {}
};

struct DPMatrix
{
    std::vector<SparseMatrix<DPEntry> > Djlik;
    std::vector<SparseMatrix<DPEntry *> > Dikjl;
    std::vector<LastPair> Mik_last;
    std::vector<Real> Mik;
    std::vector<int> MikTB;
    
    int cached_i;
    int cached_k;
};
