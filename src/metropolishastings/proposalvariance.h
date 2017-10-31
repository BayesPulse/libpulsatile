#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

#include <Rcpp.h>

// proposalvariance.h

// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions.

class ProposalVariance {

  public:
    // Constructor
    virtual ProposalVariance(double initial_pv,       // proposal variance
                             int adjust_iter, // adjust pv on multiples of adjust_iter
                             int max_iters);  // maximum iteration to adjust pv
    virtual ~ProposalVariance(); // Destructor
    virtual adjustpv();
    void addaccept() { ++accept_ct; ++iter_ct; };            // Add to acceptance count
    void addreject() { ++iter_ct; };                         // Add to iters but not accept count
    double getratio() const { return accept_ct / iter_ct; };
    void resetratio() { accept_ct = 0; iter_ct = 0; };
    double getpv() const { return pv; };

  private:
    int accept_ct; // acceptance count
    int iter_ct;   // iteration count
    double pv;     // proposal variance
    //double ratio;  // acceptance ratio (accept_ct/iter_ct) -- // only calculate when needed (on print or adjust)
    int adjust_iter; // iteration to adjust on
    int max_iter;

};

#endif
