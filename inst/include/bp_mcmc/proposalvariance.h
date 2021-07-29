#ifndef GUARD_bp_mcmc_proposalvariance_h
#define GUARD_bp_mcmc_proposalvariance_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/counter.h>


//
// proposalvariance.h
//   Definitions for the ProposalVariance classes
//
// Author: Matt Mulvahill
// Created: 10/13/17
//
// Description:
//   The ProposalVariance and ProposalVariance2p classes handle the storing and
//   adjustment of the proposal variances (PV) used in the Metropolis Hastings
//   algorithms.  The adjustment of the PV occurs automatically since we are
//   using modified MH algorithms.
//
//   These are vanilla C++ classes.  The 2p version handles 2 dimensional matrix
//   of PVs for the bivariate modifed MH case.
//
//   The classes inherit a Counter object that handles tracking the current
//   iteration and acceptance counts. This object is public, so that the
//   iteration/acceptance counts and ratios can be accessed directly by the
//   MetropolisHastings class, which inherits the ProposalVariance[2p] class.
//
// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions or virtual functions.
//


//
// ProposalVariance Class definition - single parameter version
//
class ProposalVariance {

  public:

    // Constructors
    //ProposalVariance()
    //  : psd(sqrt(5))
    //  , count()
    //  , adjust_iter(500)
    //  , max_iter(25000)
    //  , target_ratio(0.35) { };
    ProposalVariance(double in_pv,
                     int in_adjust_iter,   // adjust pv on multiples of adjust_iter
                     int in_max_iter,     // maximum iteration to adjust pv
                     double in_target_ratio) :
      count(),
      pv(in_pv),
      adjust_iter(in_adjust_iter),
      max_iter(in_max_iter),
      target_ratio(in_target_ratio) { }

    // ProposalVariance functions
    double getpv()  { return pv; }
    double getpsd() { return sqrt(pv);         }

    void check_adjust(int iter) {

      if (iter < (max_iter + 1) && iter % adjust_iter == 0 && iter > 0 &&
          getiter_since_reset() > (adjust_iter - 1)) {
        adjustpv(); 
      }

    }

    void adjustpv() {

      double ratio = getratio();
      double currentpv = getpv();
      double y = 1.0 + 1000.0 * pow(ratio - target_ratio, 3);

      if (y < 0.9)      pv = currentpv * 0.9;
      else if (y > 1.1) pv = currentpv * 1.1;
      resetratio();

    }

    // Counter object implementation
    //   works, but keep an eye out for a better option
    void   addreject(int iter)  { check_adjust(iter); count.addreject(); } ;
    void   addaccept(int iter)  { check_adjust(iter); count.addaccept(); } ;
    double getratio() { return count.getratio();  } ;
    void   resetratio() { count.resetratio();       } ;
    double gettargetratio()   { return target_ratio; } ;
    int    getiter_since_reset()     { return count.getiter_since_reset();   } ;
    int getaccept()   { return count.getaccept(); } ;

  private:

    Counter count;
    double pv;           // proposal variance
    int adjust_iter;     // iteration to adjust on
    int max_iter;        // iteration to stop adjusting
    double target_ratio; // target proposal variance

};




//
// ProposalVariance Class definition - two-parameter version
//
class ProposalVariance2p {

  public:

    //ProposalVariance2p()
    //  : count(), adjust_iter(500), max_iter(25000), target_ratio(0.35) {
    //  arma::vec in_pv(2);
    //  in_pv.fill(0);
    //  initialize_proposals(in_pv);
    //}
    ProposalVariance2p(arma::vec in_pv,
                       int in_adjust_iter,
                       int in_max_iter,
                       double in_target_ratio) :
      count(),
      adjust_iter(in_adjust_iter),
      max_iter(in_max_iter),
      target_ratio(in_target_ratio) {
        initialize_proposals(in_pv);
    }


    // ProposalVariance functions
    arma::mat getpv() const  { return pv; };
    arma::mat getpsd() const { return psd; };

    void check_adjust(int iter) {

      if (iter < (max_iter + 1) && iter % adjust_iter == 0 && iter > 0 &&
          getiter_since_reset() > (adjust_iter - 1)) {
        adjustpv(-0.90);
      }

    }


    void adjustpv(double corr = -0.5) {

      // identity matrix
      arma::mat mydiag(2, 2, arma::fill::eye);

      // y - new diagonal elements of proposal variance-covariance matrix based
      // on inputs
      double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);

      if (y < .90) {

        y = .90;
        pv  = pv % mydiag; // set off-diag to 0 with Schur/Hadamard multiplication
        set_proposals(pv * y, corr);

      } else if (y > 1.1) {

        y = 1.1;
        pv  = pv % mydiag; // set off-diag to 0 with Schur/Hadamard multiplication
        set_proposals(pv * y, corr);

      }

      resetratio();

    }

    // Counter implementation
    //   NOTE: works, but keep an eye out for a better option
    void addreject(int iter)  { check_adjust(iter); count.addreject(); } ;
    void addaccept(int iter)  { check_adjust(iter); count.addaccept(); } ;
    double getratio() { return count.getratio();  } ;
    void resetratio() { count.resetratio();       } ;
    int gettargetratio() { return target_ratio; } ;
    int getiter_since_reset() { return count.getiter_since_reset();   } ;
    int getaccept()   { return count.getaccept(); } ;

  private:

    arma::mat::fixed<2, 2> pv;
    arma::mat::fixed<2, 2> psd; // proposal standard deviation
    Counter count;
    int adjust_iter;     // iteration to adjust on
    int max_iter;        // iteration to stop adjusting
    double target_ratio; // target proposal variance

    void initialize_proposals(arma::vec initial_pv) {
      arma::mat this_pv = arma::diagmat(initial_pv);
      set_proposals(this_pv, -0.9);
    }

    void set_proposals(arma::mat this_pv, double corr) {
      pv = this_pv;
      pv(0, 1) = pv(1, 0) = calc_covariance(pv, corr);
      psd = arma::chol(pv);
    }

    double calc_covariance(arma::mat pv, double corr) {
      return corr * sqrt(pv(0, 0) * pv(1, 1));
    }

};



#endif
