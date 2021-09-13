/////////////////////////////////////////////////////////////////
// SubgradientDescent.hpp
//
// This file contains an implementation of the subgradient
// descent optimization algorithm.
/////////////////////////////////////////////////////////////////

#ifndef SUBGRADIENTDESCENT_HPP
#define SUBGRADIENTDESCENT_HPP

#include <vector>
#include "Utilities.hpp"

/////////////////////////////////////////////////////////////////
// SubgradientDescent()
//
// Implementation of subgradient optimization routine.
/////////////////////////////////////////////////////////////////

template<class RealT>
class SubgradientDescent
{
    const RealT ETA;
    const int MAX_ITERATIONS;
    const RealT NORM_BOUND;
    const RealT TARGET_LOWER_BOUND;
    const RealT TARGET_TERMINATION_RATIO;
    const RealT SUBGRADIENT_TERMINATION_RATIO;
    const int MAX_BAD_ITERATIONS;
    const RealT IMPROVEMENT_DECREASE_RATIO;

public:
    SubgradientDescent
    (
        const RealT   ETA                            = Real(1e-4),    // base learning rate
        const int     MAX_ITERATIONS                 = 1000,          // maximum number of iterations to run subgradient descent
        const RealT   NORM_BOUND                     = RealT(1e-5),   // maximum parameter vector norm
        const RealT   TARGET_LOWER_BOUND             = RealT(-1e20),  // target function value lower bound
        const RealT   TARGET_TERMINATION_RATIO       = RealT(1e-5),   // required ratio of target suboptimality for termination
        const RealT   SUBGRADIENT_TERMINATION_RATIO  = RealT(1e-5),   // required ratio of subgradient norm to parameter norm for termination
        const int     MAX_BAD_ITERATIONS             = 10,            // maximum number of bad iterations before changing target
        const RealT   IMPROVEMENT_DECREASE_RATIO     = RealT(0.5)     // multiplier for target improvement when unable to improve
    );
    
    virtual ~SubgradientDescent() {}
    
    RealT Minimize(std::vector<RealT> &x0);
    
    virtual RealT ComputeFunction(const std::vector<RealT> &x) = 0;
    virtual void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x) = 0;
    virtual void Report(int iteration, RealT f, const std::vector<RealT> &x, const std::vector<RealT> &g,
                        RealT norm_bound, RealT step_size) = 0;
    virtual void Report(const std::string &s) = 0;
};

#include "SubgradientDescent.ipp"

#endif
