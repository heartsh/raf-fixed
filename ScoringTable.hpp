/////////////////////////////////////////////////////////////////
// ScoringTable.hpp
/////////////////////////////////////////////////////////////////

#ifndef SCORINGTABLE_HPP
#define SCORINGTABLE_HPP

#include "Options.hpp"
#include "SparseMatrix.hpp"
#include "Utilities.hpp"

class ScoringTable
{
public:
    std::vector<Real> scores;
    Real offset;

    
    ScoringTable();
    ScoringTable(const std::vector<Real> &scores, Real offset);
    ScoringTable(const std::vector<Real> &M,
                 const std::vector<Real> &U,
                 const std::vector<Real> &V,
                 const Real offset);
    ScoringTable(const ScoringTable &rhs);
    ScoringTable &operator=(const ScoringTable &rhs);
};

#endif
