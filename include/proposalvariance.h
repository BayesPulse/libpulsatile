#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "counter.h"


//
// proposalvariance.h
//   Definitions for the ProposalVariance classes
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


//
// ProposalVariance Class definition - single parameter version
//
class ProposalVariance {

  public:

    // Constructors
    ProposalVariance()
      : pv(0)
      , count()
      , adjust_iter(10)
      , max_iter(25000)
      , target_ratio(0.35) { };
    ProposalVariance(double in_pv,
                     int in_adjust_iter,   // adjust pv on multiples of adjust_iter
                     int in_max_iter,     // maximum iteration to adjust pv
                     double in_target_ratio) {
      adjust_iter  = in_adjust_iter;
      max_iter     = in_max_iter;
      target_ratio = in_target_ratio;
      Counter count;
      initialize_proposals(in_pv);
      std::cout << "(in constructor) inpv is = " << in_pv << std::endl;
    }

    // ProposalVariance functions
    double getpsd() {
      std::cout << "(in getpsd) psd is = " << psd << std::endl;
      return sqrt(pv);
    }
    double getpv() {
      //check_adjust();
      return pv; 
    }
    void check_adjust() {
      int iter = count.getiter();
      //std::cout << " iter < max_iter = " << iter << " < " << max_iter << " = " << (iter < max_iter) << std::endl;
      //std::cout << " modulo math gives us: iter % adjust_iter = " << iter << " % " << adjust_iter << " = " << iter % adjust_iter << std::endl;
      if (iter < max_iter && iter % adjust_iter == 0) adjustpv(); 
      //std::cout << " now back to getting pv/psd" << std::endl;
    }
    void adjustpv() {
      double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);
      if (y < 0.9)        {
        set_proposals(pv * 0.9);
        //std::cout << "adjusted! pv = " << pv << std::endl;
      } else if (y > 1.1) {
        set_proposals(pv * 1.1);
        //std::cout << "adjusted!" << pv << std::endl;
      }
    }

    // Counter object implementation
    //   works, but keep an eye out for a better option
    void addreject()  { count.addreject();        } ;
    void addaccept()  { count.addaccept();        } ;
    double getratio() { return count.getratio();  } ;
    void resetratio() { count.resetratio();       } ;
    int getiter()     { return count.getiter();   } ;
    int getaccept()   { return count.getaccept(); } ;

  private:

    double pv;
    double psd;          // proposal standard deviation
    Counter count;
    int adjust_iter;     // iteration to adjust on
    int max_iter;        // iteration to stop adjusting
    double target_ratio; // target proposal variance

    // ProposalVariance internal functions
    void initialize_proposals(double initial_pv) {
      std::cout << "(in initialize_prop) inpv is = " << initial_pv << std::endl;

      set_proposals(initial_pv);
    }
    void set_proposals(double this_pv) {
      std::cout << "(in set_prop) inpv is = " << this_pv << std::endl;
      this->pv = this_pv;
      this->psd = sqrt(pv);
      std::cout << "(in set_prop) pv is = " << pv << std::endl;
      std::cout << "(in set_prop) psd is = " << psd << std::endl;
    }

};




//
// ProposalVariance Class definition - two-parameter version
//
class ProposalVariance2p {

  public:

    ProposalVariance2p()
      : count(), adjust_iter(500), max_iter(25000), target_ratio(0.35) {
      arma::vec in_pv(2);
      in_pv.fill(0);
      initialize_proposals(in_pv);
    }

    ProposalVariance2p(arma::vec in_pv,
                       int in_adjust_iter,
                       int in_max_iter,
                       double in_target_ratio) {
      adjust_iter  = in_adjust_iter;
      max_iter     = in_max_iter;
      target_ratio = in_target_ratio;
      initialize_proposals(in_pv);
    }


    // ProposalVariance functions
    arma::mat getpv() const  { return psd*psd; };
    arma::mat getpsd() const { return psd; };

    void adjustpv(double corr = -0.90) {
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
    }

    // Counter implementation
    //   NOTE: works, but keep an eye out for a better option
    void addreject() { count.addreject(); };
    void addaccept() { count.addaccept(); };
    double getratio() { return count.getratio(); };
    void resetratio() { count.resetratio(); };
    int getiter() { return count.getiter(); };
    int getaccept() { return count.getaccept(); };

  private:

    arma::mat::fixed<2, 2> pv;
    arma::mat::fixed<2, 2> psd; // proposal standard deviation
    Counter count;
    int adjust_iter;     // iteration to adjust on
    int max_iter;        // iteration to stop adjusting
    double target_ratio; // target proposal variance

    void initialize_proposals(arma::vec initial_pv) {
      arma::mat this_pv = arma::diagmat(initial_pv);
      set_proposals(this_pv, -0.90);
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
