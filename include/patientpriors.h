#ifndef GUARD_patientpriors_h
#define GUARD_patientpriors_h

//#include <RcppArmadillo.h>
#include <RInside.h>
#include "utils.h"

//
// patientpriors.h
//   defining the PatientPriors classes for 
//
// Author: Matt Mulvahill
// Created: 12/30/17
// Notes:
//

using namespace Rcpp;


//
// PopulationToPatientParms structures
//
//   a. PopulationToPatientParms is the base structure containing the components
//      common to all data structures representing the hierarchical level for
//      patient-level priors/population estimates
//
//   1. PopEst is the structure representing this level for the population
//      model, containing two additional parameters. This object holds
//      parameters that are updated via mcmc at the population level and used as
//      'priors' in patient level samplers.
//
//   2. PatPriors is the structure representing this level in the
//      single-subject model.  These are all user-set parameters to the prior
//      distributions.
//

// a. PatientPriors structure
struct PatientTopParms {

  // Used in all models
  double baseline_mean;
  double baseline_variance;
  double halflife_mean;
  double halflife_variance;
  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;

  // Base constructor
  PatientTopParms(double prior_baseline_mean,
                  double prior_baseline_variance,
                  double prior_halflife_mean,
                  double prior_halflife_variance,
                  double prior_mass_mean,
                  double prior_mass_variance,
                  double prior_width_mean,
                  double prior_width_variance) {

    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_variance;
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_variance;
    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_variance;
    width_mean        = prior_width_mean;
    width_variance    = prior_width_variance;

  }

};


// 1. PopulationEstimates struct
struct PopEst : PatientTopParms {
  // Only used in Population models:
  // Haven't quite wrapped my head around these (is it sd or variance), beleive
  // they correspond to mass_sd and width_sd in patient estimates for
  // single-subj model.
  double mass_mean_sd;
  double mass_mean_variance;
  double width_mean_sd;
  double width_mean_variance;

  // Population constructor (Base + 2 parameters)
  PopEst(double prior_baseline_mean,
         double prior_baseline_variance,
         double prior_halflife_mean,
         double prior_halflife_variance,
         double prior_mass_mean,
         double prior_mass_variance,
         double prior_width_mean,
         double prior_width_variance,
         double prior_mass_mean_sd,
         double prior_width_mean_sd
        ) :
    PopulationToPatientParms(prior_baseline_mean, prior_baseline_variance,
                             prior_halflife_mean, prior_halflife_variance,
                             prior_mass_mean, prior_mass_variance,
                             prior_width_mean, prior_width_variance) {

    mass_mean_sd  = prior_mass_mean_sd;
    width_mean_sd = prior_width_mean_sd;

  }

};


// 2. PatientPriors_Single struct
struct PatPriors : PatientTopParms {

  // Member variables not used in Population model:
  double mass_sd_max;
  double width_sd_max;
  double error_alpha;
  double error_beta;
  double num_orderstat;
  int    pulse_count;             // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction
  //double strauss_hardcore_range;  // range of hardcore interaction (only need
  //                                //one for single subject model and doesn't
  //                                //exist in this struct in pop model)

  // TODO: Will need constructor for Xsingle-subject, Xpopulation, and possibly
  // separate ones for trigger and response hormones (so up to 4)

  // Single-subject constructor (Base + 8 parameters)
  PatPriors(double prior_baseline_mean,
            double prior_baseline_variance,
            double prior_halflife_mean,
            double prior_halflife_variance,
            double prior_mass_mean,
            double prior_mass_variance,
            double prior_width_mean,
            double prior_width_variance,
            double prior_mass_sd_max,
            double prior_width_sd_max,
            double prior_error_alpha,
            double prior_error_beta,
            int    prior_pulse_count,
            double prior_strauss_repulsion,
            //double prior_strauss_hardcore_range, // not in single-subj
            double prior_strauss_repulsion_range
           ) :
    PatientTopParms(prior_baseline_mean, prior_baseline_variance,
                    prior_halflife_mean, prior_halflife_variance,
                    prior_mass_mean, prior_mass_variance, prior_width_mean,
                    prior_width_variance) {

      mass_sd_max  = prior_mass_sd_max;
      width_sd_max = prior_width_sd_max;
      error_alpha  = prior_error_alpha;
      error_beta   = prior_error_beta;

      num_orderstat = pulseutils::orderstat_default();

      pulse_count             = prior_pulse_count;
      strauss_repulsion       = prior_strauss_repulsion;
      //strauss_hardcore_range  = prior_strauss_hardcore_range;
      strauss_repulsion_range = prior_strauss_repulsion_range;

    }

};



//// My version of Handler class
//template <typename T> class PatientEstimates
//{
//public:
//    Handle():baseItem(NULL), refCount(new size_t(0)) {}
//    Handle(T& object):baseItem(object.clone()), refCount(new size_t(1)) {}
//
//    Handle(const Handle<T>& other):baseItem(other.baseItem), refCount(other.refCount) { ++*refCount; }
//
//    Handle& operator= (const Handle<T>& other)
//    {
//        ++*other.refCount;
//        dec_count();
//        baseItem = other.baseItem;
//        refCount = other.refCount;
//        return *this;
//    }
//
//    const T* operator->() const {return baseItem;};
//    const T& operator*() const {return *baseItem;};
//
//    T* operator->() {return baseItem;};
//    T& operator*() {return *baseItem;};
//
//    virtual ~Handle(void)
//    {
//        dec_count();
//    }
//private:
//    T *baseItem;
//    std::size_t* refCount;
//    void dec_count()
//    {
//        if (-- *refCount == 0 && baseItem != NULL)
//        {
//            delete baseItem;
//            delete refCount;
//        }
//    }
//};


//// 
//// https://ideone.com/DB7L9p
//// https://stackoverflow.com/questions/18760731/the-handle-class-in-c
////
//#include <iostream>
//#include <vector>
//using std::cout;
//using std::endl;
//using std::vector;
//
//class BaseItem{
//public:
//    virtual BaseItem* clone()
//    {
//        return new BaseItem(*this);
//    }
//    virtual void sayHello()
//    {
//        cout<<"Hello, I am class BaseItem!"<<endl;
//    }
//    virtual ~BaseItem() {}
//};
//
//class ChildItem:public BaseItem{
//public:
//    ChildItem* clone()
//    {
//        return new ChildItem(*this);
//    }
//    void sayHello(){
//        cout<<"Hello, I am class ChildItem!"<<endl;
//    }
//};
//
//
//template <typename T> class PatientEstimates
//{
//public:
//    Handle():baseItem(NULL), refCount(new size_t(0)) {}
//    Handle(T& object):baseItem(object.clone()), refCount(new size_t(1)) {}
//
//    Handle(const Handle<T>& other):baseItem(other.baseItem), refCount(other.refCount) { ++*refCount; }
//
//    Handle& operator= (const Handle<T>& other)
//    {
//        ++*other.refCount;
//        dec_count();
//        baseItem = other.baseItem;
//        refCount = other.refCount;
//        return *this;
//    }
//
//    const T* operator->() const {return baseItem;};
//    const T& operator*() const {return *baseItem;};
//
//    T* operator->() {return baseItem;};
//    T& operator*() {return *baseItem;};
//
//    virtual ~Handle(void)
//    {
//        dec_count();
//    }
//private:
//    T *baseItem;
//    std::size_t* refCount;
//    void dec_count()
//    {
//        if (-- *refCount == 0 && baseItem != NULL)
//        {
//            delete baseItem;
//            delete refCount;
//        }
//    }
//};
//
//int main()
//{
//    BaseItem item1;
//    ChildItem item2;
//
//    vector<Handle<BaseItem> > vec;
//    vec.push_back(Handle<BaseItem>(item1));
//    vec.push_back(Handle<BaseItem>(item2));
//
//    for (vector<Handle<BaseItem> >::iterator iter = vec.begin();
//      iter != vec.end(); iter++)
//    {
//      (*(*iter)).sayHello();
//    }
//    return 0;
//}
