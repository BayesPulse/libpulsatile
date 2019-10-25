//#ifndef GUARD_bp_mcmc_birthdeath_h
//#define GUARD_bp_mcmc_birthdeath_h
//
//#include <Rcpp.h>
//#include <RcppArmadillo.h>
//#ifndef NORINSIDE
//#include <RInside.h>
//#endif
//#include <bp_mcmc/proposalvariance.h>
////#include <bp_datastructures/patient.h>
////#include <bp_datastructures/datastructures.h>
//#include <bp_mcmc/utils.h>
//
//// birth_death.h
////   Abstract template class for defining birth death samplers
////
//// Author: Matt Mulvahill
//// Created: 09/04/18
//// Notes:
//
//// T - sampling unit (e.g.1. pulse iterator)  (e.g.2 patient)
//// U - Container of sampling unit  (e.g.1. patient iterator) (e.g.2 population)
//// SAMPLETYPE - type of object being sampled (double/int/arma::vec)
//// PV - proposal variance type corresponding to SAMPLETYPE (double/int/arma::mat)
//template <typename T, typename U, typename SAMPLETYPE, typename PV>
//class BirthDeathProcess
//{
//
//  public:
//    //T * sampling_unit; // either patient or population class
//
//
//class BirthDeathProcess
//{
//
//  public:
//
//    void sample(Patient *patient, bool response_hormone, int iter) {
//      
//      total_birth_rate = 
//
//
//
//    };
//
//  private:
//
//    PulseUtils pu;
//
//    void add_new_pulse(Patient *patient, double position);
//    void remove_pulse(Patient *patient, arma::vec death_rate, int pulse_count);
//    double calculate_total_deathrate(arma::vec death_rate, double pulse_count);
//    double calculate_total_deathrate_original(arma::vec death_rate, double pulse_count);
//    arma::vec calc_death_rate_strauss(Patient *patient,
//                                      arma::vec partial_likelihood,
//                                      int pulse_count,
//                                      bool response_hormone);
//    // arma::vec calc_death_rate_os(Patient *patient,
//    //                              double *partial_likelihood,
//    //                              double total_death_rate,
//    //                              int pulse_count,
//    //                              response_hormone);
//
//
//};
