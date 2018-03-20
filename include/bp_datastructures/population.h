#ifndef GUARD_population_h
#define GUARD_population_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include "bp_datastructures/patient.h"
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"

//
// population.h
//   defining the population class and subclasses
//


// In context of population model, patient priors are population estimates.
//typedef struct PatientPriors PopulationEstimates;


// The user defined values for the priors on the population level parameters
// ** This structure contains all the variables that the user sets when setting the priors **
struct PopulationPriors {

  double baseline_mean;
  double baseline_variance;
  double baseline_sd_max;  // just fyi, only relevant in pop models
  double halflife_mean;
  double halflife_variance;
  double halflife_sd_max;  // just fyi, only relevant in pop models

  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;

  // Objects this class takes over from (single-subject) Patient level classes
  double mass_sd_max;
  double width_sd_max;

  double error_alpha;
  double error_beta;

  double num_pulses;              // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_hardcore_range;  // range of hardcore interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction

};



// Parameters in the cox process for the response hormone intensity
struct AssocEstimates {
  double cluster_size;      // cluster size (rho)
  double log_cluster_size;  // log scale cluster size
  double cluster_width;     // cluster width (nu)
  double log_cluster_width;
};



struct Population {

  // Member objects
  std::vector<Patient> *patients;

  PopulationPriors *priors;
  PopulationEstimates *estimates;
  AssocEstimates *associations;

  // Constructor
  Population(std::vector<Patient> *in_patients,
             PopulationPriors *in_priors,
             PopulationEstimates *in_estimates,
             AssocEstimates *in_associations) {

    patients     = in_patients;
    priors       = in_priors;
    estimates    = in_estimates;
    associations = in_associations;

  }

};



#endif

