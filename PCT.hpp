//////////////////////////////////////////////////////////////////////
// PCT.hpp
//////////////////////////////////////////////////////////////////////

#ifndef PCT_HPP
#define PCT_HPP

#include "SparseMatrix.hpp"
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// class PCT
//////////////////////////////////////////////////////////////////////

template<class RealT>
class PCT
{
    std::vector<SparseMatrix<RealT> *> &base_pairing_posteriors;
    std::vector<SparseMatrix<RealT> *> &alignment_posteriors;
    const int m;
    const bool toggle_verbose;
    const int num_iterations;

    void Accumulate(std::vector<RealT> &res, int x, int y);

public:
    PCT(std::vector<SparseMatrix<RealT> *> &base_pairing_posteriors,
        std::vector<SparseMatrix<RealT> *> &alignment_posteriors,
        const bool toggle_verbose,
        const int num_iterations);
    
    void Transform();
};

#include "PCT.ipp"

#endif
