/////////////////////////////////////////////////////////////////
// ScoringTable.cpp
/////////////////////////////////////////////////////////////////

#include "ScoringTable.hpp"

/////////////////////////////////////////////////////////////////
// ScoringTable::ScoringTable()
//
// Default constructor.
/////////////////////////////////////////////////////////////////

ScoringTable::ScoringTable() :
    scores(), offset(0)
{}

/////////////////////////////////////////////////////////////////
// ScoringTable::ScoringTable()
//
// Constructor.
/////////////////////////////////////////////////////////////////

ScoringTable::ScoringTable(const std::vector<Real> &scores, Real offset) :
    scores(scores), offset(offset)
{}

/////////////////////////////////////////////////////////////////
// ScoringTable::ScoringTable()
//
// Constructor.  Given M[i,j], U[i], and V[j], compute
// M[i,j] - U[i] - V[j].
/////////////////////////////////////////////////////////////////

ScoringTable::ScoringTable(const std::vector<Real> &M,
                           const std::vector<Real> &U,
                           const std::vector<Real> &V,
                           const Real offset) :
    scores(M.size()), offset(offset)
{
    const int rows = int(U.size());
    const int cols = int(V.size());
    Assert(rows * cols == int(M.size()), "Dimension mismatch.");
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            scores[i*cols+j] = M[i*cols+j] - U[i] - V[j];
        }
    }
}
    
/////////////////////////////////////////////////////////////////
// ScoringTable::ScoringTable()
//
// Copy constructor.
/////////////////////////////////////////////////////////////////

ScoringTable::ScoringTable(const ScoringTable &rhs) :
    scores(rhs.scores), offset(rhs.offset)
{}

/////////////////////////////////////////////////////////////////
// ScoringTable::operator=()
//
// Assignment operator.
/////////////////////////////////////////////////////////////////

ScoringTable &ScoringTable::operator=(const ScoringTable &rhs)
{
    if (this != &rhs)
    {
        scores = rhs.scores;
        offset = rhs.offset;
    }
    return *this;
}

