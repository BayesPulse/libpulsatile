
//struct PopulationChain { NumericMatrix population; };


  //void save_sample(Population * pop);
    //arma::mat population_chain;



    // Population constructor
    // Chains(int iterations, int thin, int burnin, bool response_hormone, int num_patients) :
    //   Chains(iterations, thin, burnin, response_hormone) {
    //   association_chain(num_outputs, 4);
    //  model_type = "population";
    //   }



  //// Function for adding attributes to population_chain
  //// TODO: Straighten out variance vs SD terms and why is there a variance AND
  //// SD term for mass/width
  //NumericMatrix addattribs_population_chain(NumericMatrix out);




//void save_sample(Population * pop) {

//};

// Function for adding attributes to population_chain
// TODO: Straighten out variance vs SD terms and why is there a variance AND
// SD term for mass/width
//NumericMatrix Chains::addattribs_population_chain(NumericMatrix out) {
//  colnames(out) = CharacterVector::create("iteration",
//                                          "baseline_mean",
//                                          "baseline_variance",
//                                          "halflife_mean",
//                                          "halflife_variance",
//                                          "mass_mean",
//                                          "mass_variance",
//                                          "mass_mean_sd",
//                                          "width_mean",
//                                          "width_variance",
//                                          "width_mean_sd");
//  out.attr("class") = "population_chain";
//
//  return out;
//
//}
