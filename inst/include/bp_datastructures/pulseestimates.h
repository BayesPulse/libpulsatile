#ifndef GUARD_bp_datastructures_pulseestimates_h
#define GUARD_bp_datastructures_pulseestimates_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
//#include <bp_mcmc/utils.h>

//
// pulseestimates.h
//   defining the object for holding pulse level estimates
//
// Author: Matt Mulvahill
// Notes:
//

using namespace Rcpp;


//
// PulseEstimates structure
//   aka PulseParms -- pulse chain updated by mcmc algorithm called
//   TODO: Check  on lambda, does fsh need separate definition?
//

class PulseEstimates {

  public:

    // Variables used in all models
    double time;
    double mass;
    double width;
    double tvarscale_mass;   // variance scale for mass t-dist (eta)
    double tvarscale_width;  // variance scale for width t-dist (eta)
    //double lambda; // for fsh pulse only, denomsum - NOT SURE WHAT THIS TERM IS
                   //FOR (from Karen's code)

    // Public user-facing access function for mean_contrib. Only calculates
    // mean_contrib if input values have changed.
    // TODO: change this return type to a pointer
    arma::vec get_mean_contribution(const arma::vec & data_time, double decay_rate)
    {

      if ((prev_mass != mass) | (prev_width != width) | (prev_time != time) |
          (prev_decay_rate != decay_rate)) {
        calc_mean_contribution(data_time, decay_rate);
        prev_time = time;
        prev_mass = mass;
        prev_width = width;
        prev_decay_rate = decay_rate;
      }

      return mean_contribution;

    }

    // For returining to chain
    arma::rowvec get_vector_of_values() 
    {
      arma::rowvec out { time, mass, width, tvarscale_mass, tvarscale_width}; //, lambda };
      return out;
    }

    // Constructor for new PulseEstimates objects
    PulseEstimates(double in_time,
                   double in_mass,
                   double in_width,
                   double in_tvarscale_mass,
                   double in_tvarscale_width,
                   //double fshlambda,
                   double patient_decay,
                   const arma::vec &data_time)
        : time            (in_time)
        , mass            (in_mass)
        , width           (in_width)
        , tvarscale_mass  (in_tvarscale_mass)
        , tvarscale_width (in_tvarscale_width)
        //, lambda          (fshlambda)
        , mean_contribution(data_time.n_elem)
      {
        mean_contribution.fill(0.);
        calc_mean_contribution(data_time, patient_decay);
        prev_time       = time;
        prev_mass       = mass;
        prev_width      = width;
        prev_decay_rate = patient_decay;
      }
    // Constructor for empty pulse object
    PulseEstimates()
      : time(0), mass(0), width(0), tvarscale_mass(0),
        tvarscale_width(0), mean_contribution(1) 
      {
        mean_contribution.fill(0);
      }

  private:

    arma::vec mean_contribution;
    double prev_time, prev_mass, prev_width, prev_decay_rate;

    // mean_contribution() of each pulse to the total mean_concentration
    void calc_mean_contribution(const arma::vec &data_time, double decay_rate)
    {

      double y, z, w;
      arma::vec x(data_time.n_elem);
      x.fill(0.);

      z  = width * decay_rate;
      y  = decay_rate * (0.5 * z  + time);
      z += time;
      w  = sqrt(2. * width);
      x = ((data_time - z) / w) * sqrt(2.);

      double N = data_time.n_elem;
      // NOTE: potentially slow piece
      for (int i = 0; i < N; i++) {
        x(i) = Rf_pnorm5(x(i), 0.0, 1.0, 1, 0);
        //x = arma::normpdf(x); // doesn't give same results
      }

      //Rcpp::Rcout << "New Iteration --------------------------\n"
      //            << "z = " << z << "; y = " <<  y << "; w = " << w
      //            << "\nDecay rate: " << decay_rate
      //            << "\ntime: " << time << "\n"
      //            << "Sum X: " << sum(x) << "\n"
      //            << "mass: " << mass << "\n"
      //            << "width: " << width << "\n";

      // Finish calculating mean_contrib w/ vectorized ops
      //    mass * x is a vector, as is exp(), so use element-wise
      //    multiplication via %
      mean_contribution = (mass * x) % trunc_exp(y - data_time * decay_rate);
      //Rcpp::Rcout << "Sum(mean_contribution) 1 = " << sum(mean_contribution) << "\n";
      // Truncate <0 = 0
      mean_contribution.for_each( [](arma::vec::elem_type& val) { val = std::max(val, 0.); } );
      //Rcpp::Rcout << "Sum(mean_contribution) 2 = " << sum(mean_contribution) << "\n\n";

    }


};



#endif
