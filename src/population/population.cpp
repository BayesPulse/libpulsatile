#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bp_population/bp_population.h>
//#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <iostream>
#ifndef NORINSIDE
#include <RInside.h>
#endif

using namespace Rcpp;

// population.cpp
//   Implementation of population model.
//
// Author: Max McGrath
// Notes:
//

// [[Rcpp::export]]
Rcpp::List population_(Rcpp::NumericMatrix concentrations,
                       Rcpp::NumericVector time,
                       Rcpp::CharacterVector location_prior,
                       Rcpp::List inpriors,
                       Rcpp::List proposalvars,
                       Rcpp::List startingvals,
                       int mcmc_iterations,
                       int thin,
                       int burnin,
                       bool verbose,
                       int pv_adjust_iter,
                       int pv_adjust_max_iter,
                       double bivariate_pv_target_ratio,
                       double univariate_pv_target_ratio,
                       bool test_birthdeath,
                       bool test_sd_masses,
                       bool test_sd_width,
                       bool test_s2s_sd_width,
                       bool test_s2s_sd_mass,
                       bool test_s2s_sd_baseline,
                       bool test_s2s_sd_halflife,
                       bool test_pop_means_width,
                       bool test_pop_means_mass,
                       bool test_pop_means_baseline,
                       bool test_pop_means_halflife,
                       bool test_fixeff_mass,
                       bool test_fixeff_width,
                       bool test_blhl,
                       bool test_error,
                       bool test_locations,
                       bool test_masses,
                       bool test_widths,
                       bool test_tvarscale_mass,
                       bool test_tvarscale_width,
                       Rcpp::NumericVector testMassVec,
                       Rcpp::NumericVector testWidthVec,
                       Rcpp::NumericVector testMKappaVec,
                       Rcpp::NumericVector testWKappaVec,
                       Rcpp::NumericVector testLocVec)
{
  // Set to print every nth iteration for verbose screen output  
  int verbose_iter = 5000;

  // Check for valid input
  if ( !inpriors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !proposalvars.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !startingvals.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  if ( concentrations.nrow() != time.size() ) stop("Time and concentration vectors must be the same size");
  if ( location_prior.size() != 1 ) stop("Location prior vector must be length 1");
  
  // Cleanup variable names 
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;
  std::string loc_prior = Rcpp::as<std::string>(location_prior);

  // Test block for testing trunc_exp
  //arma::vec test_vector(3);
  //test_vector.fill(pow(10, 2));
  //Rcpp::Rcout << "10^2" << trunc_exp(test_vector);
  //Rcpp::Rcout << "10^2" << exp(test_vector);
  //test_vector.fill(pow(10, 4));
  //Rcpp::Rcout << "10^4" << trunc_exp(test_vector);
  //Rcpp::Rcout << "10^4" << exp(test_vector);
  //test_vector.fill(pow(10, 8));
  //Rcpp::Rcout << "10^8" << trunc_exp(test_vector);
  //Rcpp::Rcout << "10^8" << exp(test_vector);
  //Rcpp::Rcout << "Test equation: "
  //            << (test_vector * 0) % exp(test_vector) << "\n"
  //            << (test_vector * 0) % trunc_exp(test_vector) << "\n";
  //Rcpp::Rcout << "Test equation 2: "
  //            << (test_vector * 1.27e-254) % trunc_exp(test_vector) << "\n";
  //
  //arma::vec test_vec2 = (test_vector * 0) % exp(test_vector);
  //arma::vec test_vec3 = (test_vector * 0) % trunc_exp(test_vector);

  //Rcpp::Rcout << "Sums: " << sum(test_vec2) << " " << sum(test_vec3) << "\n";

  //stop("test");


  Rcpp::Rcout << "Location prior is: " << loc_prior << "\n\n";
  
  // Create population priors and estimates objects
  PopulationPriors popPriors(inpriors["mass_mean"],
                             inpriors["mass_var"],                
                             inpriors["mass_p2p_sd_var"],     
                             inpriors["mass_s2s_sd_var"],     
                             inpriors["width_mean"],              
                             inpriors["width_var"],               
                             inpriors["width_p2p_sd_var"],    
                             inpriors["width_s2s_sd_var"],
                             inpriors["baseline_mean"],
                             inpriors["baseline_var"],
                             inpriors["baseline_s2s_sd_var"], 
                             inpriors["halflife_mean"],           
                             inpriors["halflife_var"],            
                             inpriors["halflife_s2s_sd_var"],
                             inpriors["error_alpha"],
                             inpriors["error_beta"]);

  PopulationEstimates popEstimates(startingvals["mass_mean"],
                                   startingvals["mass_s2s_sd_var"],
                                   startingvals["mass_p2p_sd_var"],
                                   startingvals["width_mean"],
                                   startingvals["width_s2s_sd_var"],
                                   startingvals["width_p2p_sd_var"],
                                   startingvals["baseline_mean"],
                                   startingvals["baseline_s2s_sd_var"],
                                   startingvals["halflife_mean"],
                                   startingvals["halflife_s2s_sd_var"]);

  Rcpp::Rcout << "Population priors and estimates created\n";

  PatientPriors patientPriors(startingvals["mass_mean"],
                              startingvals["width_mean"],
                              startingvals["baseline_mean"],
                              startingvals["halflife_mean"],
                              startingvals["mass_s2s_sd_var"],
                              startingvals["width_s2s_sd_var"],
                              startingvals["baseline_s2s_sd_var"],
                              startingvals["halflife_s2s_sd_var"],
                              startingvals["mean_pulse_count"],
                              startingvals["strauss_repulsion"],
                              startingvals["strauss_repulsion_range"]);

  // Create subject level estimates object
  PatientEstimates patEstimates(startingvals["baseline_mean"],
                                startingvals["halflife_mean"],
                                startingvals["error_var"],
                                startingvals["mass_mean"],
                                startingvals["width_mean"],
                                startingvals["mass_p2p_sd_var"],
                                startingvals["width_p2p_sd_var"],
                                false);

  // Get number of patients (should this be a function arg?)
  int numPats = concentrations.ncol();
  
  // Create vector of patients
  std::vector<Patient> pats;

  // Fill vector of patients
  for(int i = 0; i < numPats; i++) {
    PatientData tempdata(time, concentrations.column(i));
    Patient temppatient(tempdata, patientPriors, patEstimates);
    pats.push_back(temppatient);
  }

  // Create population object and pointer to it
  Population pop(pats, popPriors, popEstimates);
  Population * population = &pop;
    
  //----------------------------------------
  // Construct MMH Objects
  //----------------------------------------
  BirthDeathProcess birth_death;
  
  // Population Level   
  Pop_DrawSDRandomEffects draw_sd_masses(proposalvars["ind_mass"], adj_iter,
                                         adj_max, univ_target, false,
                                         verbose, verbose_iter);
  Pop_DrawSDRandomEffects draw_sd_width(proposalvars["ind_width"], adj_iter,
                                         adj_max, univ_target, true,
                                         verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_width(proposalvars["sub_width"], adj_iter, adj_max, univ_target,
                                         true, false, false, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_mass(proposalvars["sub_mass"], adj_iter, adj_max, univ_target,
                                         false, true, false, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_baseline(proposalvars["sub_baseline"], adj_iter, adj_max, univ_target,
                                         false, false, true, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_halflife(proposalvars["sub_halflife"], adj_iter, adj_max, univ_target,
                                         false, false, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_width(proposalvars["pop_width_mean"], adj_iter, adj_max, univ_target,
                                         true, false, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_mass(proposalvars["pop_mass_mean"], adj_iter, adj_max, univ_target,
                                         false, true, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_baseline(proposalvars["pop_baseline_mean"], adj_iter, adj_max, univ_target,
                                         false, false, true, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_halflife(proposalvars["pop_halflife_mean"], adj_iter, adj_max, univ_target,
                                         false, false, false, verbose, verbose_iter);

  // Subject Level
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["sub_mass_mean"], adj_iter, adj_max, univ_target,
                                       false, verbose, verbose_iter);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["sub_width_mean"], adj_iter, adj_max, univ_target,
                                        true, verbose, verbose_iter);

  arma::vec blhl_pv = { proposalvars["sub_baseline"], proposalvars["sub_halflife"] };
  SS_DrawBaselineHalflife draw_blhl(blhl_pv, adj_iter, adj_max, biv_target,
                                    verbose, verbose_iter);
  SS_DrawError draw_error;

  // Pulse Level
  SS_DrawLocations * draw_locations;
  if ( loc_prior == "strauss" ) {
    draw_locations = new SS_DrawLocationsStrauss(proposalvars["pulse_location"], adj_iter, adj_max,
                                                univ_target, verbose, verbose_iter);
  } else {
    Rcpp::stop("Order statistic prior is not yet implemented in the birth death process");
  }

  SS_DrawRandomEffects draw_masses(proposalvars["sub_mass_mean"], adj_iter, adj_max, univ_target,
                                   false, verbose, verbose_iter);
  SS_DrawRandomEffects draw_widths(proposalvars["sub_width_mean"], adj_iter, adj_max, univ_target,
                                   true, verbose, verbose_iter);
  SS_DrawTVarScale draw_tvarscale_mass(proposalvars["sdscale_pulse_mass"],
                                       adj_iter, adj_max, univ_target, false,
                                       verbose, verbose_iter);
  SS_DrawTVarScale draw_tvarscale_width(proposalvars["sdscale_pulse_width"],
                                        adj_iter, adj_max, univ_target, true,
                                        verbose, verbose_iter);

  // Initialize object to store chains 
  PopChains chains(mcmc_iterations, thin, burnin, false, verbose, verbose_iter, numPats);

  //---------------------------------------------------------
  // Set population parameters (for testing)
  //---------------------------------------------------------
  if(!test_pop_means_mass) population->estimates.mass_mean = 4.36547;
  if(!test_pop_means_width) population->estimates.width_mean = 35.23362;
  if(!test_pop_means_halflife) population->estimates.halflife_mean = 43.13908;
  if(!test_pop_means_baseline) population->estimates.baseline_mean = 2.542657;

  if(!test_s2s_sd_mass) population->estimates.mass_s2s_sd = 1.0;
  if(!test_s2s_sd_width) population->estimates.width_s2s_sd = 1.4;
  if(!test_s2s_sd_halflife) population->estimates.halflife_sd = 5;
  if(!test_s2s_sd_baseline) population->estimates.baseline_sd = 0.55;

  if(!test_sd_width) population->estimates.mass_p2p_sd = 1.6;
  if(!test_sd_masses) population->estimates.width_p2p_sd = 5;

  //---------------------------------------------------------
  // Set individual mass means (for testing)
  //---------------------------------------------------------
  
  double j = 0;
  double l = 0;

  if(!test_fixeff_mass) {
    int i = 0;
    j = 0;
    for(auto &pat : population->patients) {

      switch(i) {
        case 0: j = 5.0156; break;
        case 1: j = 4.2086; break;
        case 2: j = 4.1829; break;
        case 3: j = 4.1347; break;
        case 4: j = 6.6892; break;
        case 5: j = 4.1205; break;
        case 6: j = 3.1993; break;
        case 7: j = 3.1900; break;
        case 8: j = 4.6209; break;
        case 9: j = 4.2930; break;
        default: std::cout << "Problem setting mass means\n";
      }

      pat.estimates.mass_mean = j;
      i++;
    }
  
  std::cout << "Patient mass means set\n";
  }

  //---------------------------------------------------------
  // Set individual width means (for testing)
  //---------------------------------------------------------
  if(!test_fixeff_width) {
    int i = 0;
    j = 0;
    for(auto &pat : population->patients) {

      switch(i) {
        case 0: j = 32.8066; break;
        case 1: j = 35.9662; break;
        case 2: j = 35.1099; break;
        case 3: j = 33.0451; break;
        case 4: j = 36.1311; break;
        case 5: j = 34.6148; break;
        case 6: j = 36.0120; break;
        case 7: j = 37.2471; break;
        case 8: j = 35.2640; break;
        case 9: j = 36.0495; break;
        default: std::cout << "Problem setting width means\n";
      }

      pat.estimates.width_mean = j;
      i++;
    }
  
  }
  
  //---------------------------------------------------------
  // Set individual baseline (for testing)
  //---------------------------------------------------------
  if(!test_blhl) {
    int i = 0;
    j = 0;
    for(auto &pat : population->patients) {

      switch(i) {
        case 0: j = 2.8690; break;
        case 1: j = 3.1813; break;
        case 2: j = 2.6171; break;
        case 3: j = 2.4126; break;
        case 4: j = 2.4754; break;
        case 5: j = 2.1407; break;
        case 6: j = 1.4205; break;
        case 7: j = 3.3326; break;
        case 8: j = 2.3291; break;
        case 9: j = 2.6481; break;
        default: std::cout << "Problem setting baselines\n";
      }

      pat.estimates.baseline = j;
      pat.estimates.baseline_halflife(0) = j;
      i++;
    }
  
  }

  //---------------------------------------------------------
  // Set individual halflives (for testing)
  //---------------------------------------------------------
  if(!test_blhl) {
    int i = 0;
    j = 0;
    for(auto &pat : population->patients) {

      switch(i) {
        case 0: j = 41.4493; break;
        case 1: j = 34.6105; break;
        case 2: j = 42.9464; break;
        case 3: j = 50.9084; break;
        case 4: j = 49.6113; break;
        case 5: j = 46.7311; break;
        case 6: j = 44.2629; break;
        case 7: j = 41.9211; break;
        case 8: j = 38.0089; break;
        case 9: j = 40.9405; break;
        default: std::cout << "Problem setting halflives\n";
      }

      pat.estimates.halflife = j;
      pat.estimates.baseline_halflife(1) = j;
      i++;
    }
  
  }
  
  //---------------------------------------------------------
  // Add pulses (for testing)
  //---------------------------------------------------------
  double position, new_mass, new_width, new_tvarscale_mass, new_tvarscale_width,
         new_t_sd_mass, new_t_sd_width;
  l = 0;

  if(!test_birthdeath) {

    for(int i = 0; i < 10; i++) {

      Patient * patient = &population->patients[i];

      switch(i) {
        case 0: j = 14; break;
        case 1: j = 16; break;
        case 2: j = 15; break;
        case 3: j = 16; break;
        case 4: j = 12; break;
        case 5: j = 15; break;
        case 6: j = 11; break;
        case 7: j = 16; break;
        case 8: j = 13; break;
        case 9: j = 14; break;
        default: std::cout << "Problem adding pulses\n";
      }

      Rcpp::RNGScope rng_scope;
      position = (!test_locations) ? testLocVec(l) : Rf_runif(patient->data.fitstart, patient->data.fitend);
      //new_tvarscale_mass = (!test_tvarscale_mass) ? testMKappaVec(l) : Rf_rgamma(2, 0.5);
      new_tvarscale_mass = (!test_tvarscale_mass) ? testMKappaVec(l) : 10;
      //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : Rf_rgamma(2, 0.5);
      new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;
  
      if(!test_masses) {
        new_mass = testMassVec(l);
      } else {
        new_t_sd_mass = patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
        new_mass = -1.0;
        while (new_mass < 0) {
          //new_mass = Rf_rnorm(patient->estimates.mass_mean, new_t_sd_mass);
          new_mass = Rf_rnorm(10, .1);
        }
      }

      if(!test_widths) {
        new_width = testWidthVec(l);
      } else {
        new_t_sd_width = patient->estimates.width_sd / sqrt(new_tvarscale_width);
        new_width = -1.0;
        while (new_width < 0) {
          //new_width = Rf_rnorm(patient->estimates.width_mean, new_t_sd_width);
          new_width = Rf_rnorm(70, 1);
        }
      }
    
      PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                               new_tvarscale_width, patient->estimates.get_decay(),
                               patient->data.time);

      patient->pulses.front() = new_pulse;

      l++;

      for(int k = 0; k < j-1; k++) {

        Rcpp::RNGScope rng_scope;
        position = (!test_locations) ? testLocVec(l) : Rf_runif(patient->data.fitstart, patient->data.fitend);
        //new_tvarscale_mass = (!test_tvarscale_mass) ? testMKappaVec(l) : Rf_rgamma(2, 0.5);
        new_tvarscale_mass = (!test_tvarscale_mass) ? testMKappaVec(l) : 10;
        //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : Rf_rgamma(2, 0.5);
        new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;

        if(!test_masses) {
          new_mass = testMassVec(l);
        } else {
          new_t_sd_mass = patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
          new_mass = -1.0;
          while (new_mass < 0) {
            //new_mass = Rf_rnorm(patient->estimates.mass_mean, new_t_sd_mass);
            new_mass = Rf_rnorm(10, .1);
          }
        }

        if(!test_widths) {
          new_width = testWidthVec(l);
        } else {
          new_t_sd_width = patient->estimates.width_sd / sqrt(new_tvarscale_width);
          new_width = -1.0;
          while (new_width < 0) {
            //new_width = Rf_rnorm(patient->estimates.width_mean, new_t_sd_width);
            new_width = Rf_rnorm(70, 1);
          }
        }
        
        PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                                 new_tvarscale_width, patient->estimates.get_decay(),
                                 patient->data.time);

        patient->pulses.push_back(new_pulse);

        l++;

      }

    }

  Rcpp::Rcout << "Pulses Added Manually\n";

  Rcpp::Rcout << "Details:\n";
  j = l = 1;
  for(auto pat : population->patients) {
    for(auto pulse : pat.pulses) {
      Rcpp::Rcout << std::setprecision(3) << j << " "
                  << l << " "
                  << pulse.mass << " "
                  << pulse.width << " "
                  << pulse.time << " "
                  << pulse.tvarscale_mass << " "
                  << pulse.tvarscale_width << " " << "\n";
      l++;
    }
    l = 1;
    j++;
  }

  Rcpp::Rcout << "\n";

  }




  //----------------------------------------
  // Sample MMH Objects
  //----------------------------------------

  for(int iteration = 0; iteration < mcmc_iterations; iteration++) {

    chains.print_diagnostic_output(population, iteration);

    if (test_sd_masses) draw_sd_masses.sample(population, &population->estimates.mass_p2p_sd, population, iteration);
    if (test_sd_width) draw_sd_width.sample(population, &population->estimates.width_p2p_sd, population, iteration);

    if (test_pop_means_width) draw_pop_means_width.sample(population, &population->estimates.width_mean, iteration);
    if (test_pop_means_mass) draw_pop_means_mass.sample(population, &population->estimates.mass_mean, iteration);
    if (test_pop_means_baseline) draw_pop_means_baseline.sample(population, &population->estimates.baseline_mean, iteration);
    if (test_pop_means_halflife) draw_pop_means_halflife.sample(population, &population->estimates.halflife_mean, iteration);

    if (test_s2s_sd_width) draw_s2s_sd_width.sample(population, &population->estimates.width_s2s_sd, population, iteration);
    if (test_s2s_sd_mass) draw_s2s_sd_mass.sample(population, &population->estimates.mass_s2s_sd, population, iteration);
    if (test_s2s_sd_baseline) draw_s2s_sd_baseline.sample(population, &population->estimates.baseline_sd, population, iteration);
    if (test_s2s_sd_halflife) draw_s2s_sd_halflife.sample(population, &population->estimates.halflife_sd, population, iteration);

    population->matchPatPriorsToPop();

    int j = 1;

    for(auto &pat : population->patients) {
      Patient * patient = &pat;
      if (test_birthdeath) birth_death.sample(patient, false, iteration);
      if (test_fixeff_mass) draw_fixeff_mass.sample(patient, &patient->estimates.mass_mean, iteration);
      if (test_fixeff_width) draw_fixeff_width.sample(patient, &patient->estimates.width_mean, iteration);
      if (test_blhl){ 
        draw_blhl.sample(patient, &patient->estimates.baseline_halflife, iteration);

        //Rcpp::Rcout << "Patient " << j << " Vector: " << patient->estimates.baseline_halflife(0)
        //            << " " << patient->estimates.baseline_halflife(1) << "\n"
        //            << "Before BL: " << patient->estimates.baseline 
        //            << " HL: " << patient->estimates.halflife
        //            << "\n";

        patient->estimates.matchBLHL();

        //Rcpp::Rcout << "After  BL: " << patient->estimates.baseline
        //            << " HL: " << patient->estimates.halflife
        //            << "\n\n";

      }
      if (test_locations) draw_locations->sample_pulses(patient, iteration);
      if (test_masses) draw_masses.sample_pulses(patient, iteration);
      if (test_widths) draw_widths.sample_pulses(patient, iteration);
      if (test_tvarscale_mass) draw_tvarscale_mass.sample_pulses(patient, iteration);
      if (test_tvarscale_width) draw_tvarscale_width.sample_pulses(patient, iteration);

      j++;

    }

    if (test_error) draw_error.sample(population);

    chains.save_sample(population, iteration);

  }


  delete draw_locations;
  
  return(chains.output(population));

};
