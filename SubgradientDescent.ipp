/////////////////////////////////////////////////////////////////
// SubgradientDescent.ipp
//
// This file contains an implementation of the subgradient
// descent optimization algorithm.
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// SubgradientDescent::SubgradientDescent()
//
// Constructor.
/////////////////////////////////////////////////////////////////

template<class RealT>
SubgradientDescent<RealT>::SubgradientDescent
(
    const RealT   ETA,                                        // base learning rate
    const int     MAX_ITERATIONS,                             // maximum number of iterations to run subgradient descent
    const RealT   NORM_BOUND,                                 // maximum parameter vector norm
    const RealT   TARGET_LOWER_BOUND,                         // target function value lower bound
    const RealT   TARGET_TERMINATION_RATIO,                   // required ratio of target suboptimality for termination
    const RealT   SUBGRADIENT_TERMINATION_RATIO,              // required ratio of subgradient norm to parameter norm for termination
    const int     MAX_BAD_ITERATIONS,                         // maximum number of bad iterations before changing target
    const RealT   IMPROVEMENT_DECREASE_RATIO                  // multiplier for target improvement when unable to improve
) :
    ETA(ETA),
    MAX_ITERATIONS(MAX_ITERATIONS),
    NORM_BOUND(NORM_BOUND),
    TARGET_LOWER_BOUND(TARGET_LOWER_BOUND),
    TARGET_TERMINATION_RATIO(TARGET_TERMINATION_RATIO),
    SUBGRADIENT_TERMINATION_RATIO(SUBGRADIENT_TERMINATION_RATIO),
    MAX_BAD_ITERATIONS(MAX_BAD_ITERATIONS),
    IMPROVEMENT_DECREASE_RATIO(IMPROVEMENT_DECREASE_RATIO)
{}

/////////////////////////////////////////////////////////////////
// SubgradientDescent::Minimize()
//
// Implementation of Pegasos algorithm.
/////////////////////////////////////////////////////////////////

template<class RealT>
RealT SubgradientDescent<RealT>::Minimize(std::vector<RealT> &x)
{
    std::vector<RealT> g;
    RealT f = ComputeFunction(x);
    ComputeSubgradient(g, x);

    // check early termination criteria
    
    if (f >= RealT(1e20))
    {
        Report(SPrintF("Termination before optimization: function value too big (%lf > %lf)", f, 1e20));
        return f;
    }

    for (int epoch = 1; epoch <= MAX_ITERATIONS; epoch++)
    {
        // compute learning rate

        RealT eta = ETA / RealT(epoch);
        
        // take a step

        x -= eta * g;

        // project back to ball

        RealT norm = Norm(x);
        if (norm > NORM_BOUND)
        {
            x *= NORM_BOUND / norm;
        }

        // compute new function

        f = ComputeFunction(x);

        // compute new subgradient
        
        ComputeSubgradient(g, x);
        
        // print updates

        const int update_frequency = std::max(1, MAX_ITERATIONS / 100);
        if (epoch % update_frequency == 0)
        {
            Report(epoch, f, x, g, NORM_BOUND, eta);
        }

        // check convergence criteria

        if (epoch >= MAX_ITERATIONS)
        {
            Report("Termination condition: maximum number of iterations reached");
            break; 
        }
    }
    
    return f;
}

#ifdef COMMENT

/////////////////////////////////////////////////////////////////
// SubgradientDescent::Minimize()
//
// Implementation of subgradient descent optimization routine.
/////////////////////////////////////////////////////////////////

template<class RealT>
RealT SubgradientDescent<RealT>::Minimize(std::vector<RealT> &x0)
{
    RealT f;
    std::vector<RealT> g;
    std::vector<RealT> x = x0;
    
    std::vector<RealT> x_avg(x.size());
    int ct = 0;
    
    // check for termination criteria at beginning
    
    f = ComputeFunction(x);
    if (f >= RealT(1e20))
    {
        Report(SPrintF("Termination before optimization: function value too big (%lf > %lf)", f, 1e20));
        return f;
    }

    ComputeSubgradient(g, x);
    RealT subgradient_ratio = Norm(g) / std::max(RealT(1), Norm(x));
    if (subgradient_ratio < SUBGRADIENT_TERMINATION_RATIO)
    {
        Report(SPrintF("Termination before optimization: subgradient vector small (%lf < %lf)", subgradient_ratio, SUBGRADIENT_TERMINATION_RATIO));
        return f;
    }

    // maintain best so far

    RealT fp = f + 1;
    RealT f0 = f;
    std::vector<RealT> g0 = g;
    
    // compute initial target

    RealT target_improvement = std::min(f - TARGET_LOWER_BOUND, RealT(0.5) * DotProduct(g, g));
    //RealT bad_distance = 0;

    RealT lipschitz_constant = 0;

    // subgradient descent

    for (int epoch = 1; epoch <= MAX_ITERATIONS; epoch++)
    {
        // update target improvement value

        /*
        if (f <= f0 - RealT(0.5) * target_improvement)
        {
            bad_distance = 0;
        }
        else
        {
            if (bad_distance >= 1)
            {
                target_improvement *= IMPROVEMENT_DECREASE_RATIO;
                bad_distance = 0;
                f = f0;
                g = g0;             
                x = x0;
            }
            
            bad_distance += (f - f0 + target_improvement) / Norm(g);
        }
        */
        
        if (f < f0)
        {
            target_improvement *= 2;
        }
        else if (f >= fp)
        {
            target_improvement = std::max(TARGET_TERMINATION_RATIO, target_improvement / 2);
        }
        
        // compute step size
        
        const RealT f_target = std::max(TARGET_LOWER_BOUND, f0 - target_improvement);
        //const RealT alpha = (f - f_target) / DotProduct(g, g);
        //const RealT alpha = ETA / epoch;
        //const RealT alpha = ETA / Norm(g);
        const RealT alpha = ETA;

        lipschitz_constant = std::max(lipschitz_constant, Norm(g));
       
        //const RealT alpha = Sqrt(2 * NORM_BOUND) / lipschitz_constant / Sqrt(Real(epoch));

        // save best so far

        if (f <= f0)
        {
            f0 = f;
            g0 = g;
            x0 = x;
        }

        // take step

        x -= alpha * g;
        
        RealT norm = Norm(x);
        if (norm > NORM_BOUND)
        {
            //Report(SPrintF("Rescaling from ||w|| = %lf to ||w|| = %lf...", norm, NORM_BOUND));
            x *= NORM_BOUND / norm;
        }

        // compute average

        x_avg += x;
        ct++;

        const int update_frequency = std::max(1, MAX_ITERATIONS / 100);
        if (epoch % update_frequency == 0)
        {
            Report(epoch, f0, x_avg / RealT(ct), f, x, g, NORM_BOUND, f_target, alpha);
            //Report(epoch, f0, x0, f, x, g, NORM_BOUND, f_target, alpha);
        }

        // compute new function and subgradient

        fp = f;
        f = ComputeFunction(x);
        ComputeSubgradient(g, x);

        // check convergence criteria
        
        if (epoch >= MAX_ITERATIONS)
        {
            Report("Termination condition: maximum number of iterations reached");
            break; 
        }

        subgradient_ratio = Norm(g) / std::max(RealT(1), Norm(x));
        if (subgradient_ratio < SUBGRADIENT_TERMINATION_RATIO)
        {
            Report(SPrintF("Termination before optimization: subgradient vector small (%lf < %lf)", subgradient_ratio, SUBGRADIENT_TERMINATION_RATIO));
            break;
        }

        /*
        if ((f0 - f_target) / std::max(RealT(1), f_target) < TARGET_TERMINATION_RATIO)
        {
            Report("Termination condition: target value reached");
            break;
            }*/
    }
    
    return f0;
}

#endif
