#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>

//
// proposalvariance.h
//   Definitions for the ProposalVariance class
//
// Author: Matt Mulvahill
// Created: 10/13/17
//
// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions or virtual functions.
//
// NOTE: ProposalVariance has a constructor only to document the requirements
// for the constructor.  It can't actually be used since its an abstract class.
//

template <class T>
class ProposalVariance {

  public:
    void addaccept() { ++accept_ct; ++iter_ct; }; // Add to acceptance count
    void addreject() { ++iter_ct; };              // Add to iters but not accept count
    double getratio() const { return accept_ct / iter_ct; };
    void resetratio() { accept_ct = 0; iter_ct = 0; };
    virtual void adjustpv() = 0;
     T getpv() const { return pv; };
     T getpsd() const { return psd; };
    //virtual getpv()  const { return pv; };
    //virtual getpsd() const { return psd; };

  protected:
    ProposalVariance(T initial_pv,    // proposal variance
                     int adjust_iter,      // adjust pv on multiples of adjust_iter
                     int max_iters,        // maximum iteration to adjust pv
                     double target_ratio); // target accept ratio
    ~ProposalVariance();

  private:
    int accept_ct;    // acceptance count
    int iter_ct;      // iteration count
    int adjust_iter;  // iteration to adjust on
    int max_iter;     // iteration to stop adjusting
    T pv;       // proposal variance
    T psd;      // proposal variance
    double target_ratio; // target proposal variance

};


class MarginalProposalVariance : public ProposalVariance
{

  public:
    MarginalProposalVariance(double initial_pv, int adjust_iter, int max_iters, double target_ratio);
    ~ProposalVariance(); // Destructor

    void adjustpv();
    double getpv  const { return pv; };
    double getpsd const { return psd; };

  private:
    void initialize_proposals(double this_pv);
    void set_proposals(double this_pv);
    double pv; // proposal var-covar matrix
    double psd; // proposal var-covar matrix
    double target_ratio = 0.35;

};


class JointProposalVariance : public ProposalVariance
{

  public:
    JointProposalVariance(arma::vec initial_pv, int adjust_iter, int max_iters, double target_ratio);
    ~ProposalVariance(); // Destructor

    void adjustpv(double corr = -0.90);
    arma::mat getpv()  const { return pv; };
    arma::mat getpsd() const { return psd; };

  private:
    void initialize_pv();
    void update_pv();
    double set_covariance(double var00, double var11, double corr);
    double set_covariance(double var00, double var11);
    arma::mat::fixed<2, 2> pv; // proposal var-covar matrix 
    arma::mat::fixed<2, 2> psd;  // proposal correlation matrix
    double target_ratio = 0.25;

};



#endif

