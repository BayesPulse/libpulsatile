#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <cmath>
#include <bp_mcmc/counter.h>
#include <bp_mcmc/proposalvariance.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/pulseestimates.h>
#include <bp_datastructures/utils.h>
#include <bp_mcmc/utils.h>
#include <bpmod_singlesubject/ss_draw_tvarscale.h>
#include <bpmod_singlesubject/ss_draw_randomeffects.h>
#include <bpmod_singlesubject/ss_draw_fixedeffects.h>
#include <bpmod_singlesubject/ss_draw_sdrandomeffects.h>
#include <bpmod_singlesubject/birthdeath.h>
#include <testing/catch.h>


//
// SS_DrawTVarScale::posterior_function
//
// Regression test for the t-variance-scale (kappa) acceptance ratio. The
// scaled-normal representation has theta ~ N(mu, sigma^2 / kappa), so the
// kappa MH log-ratio includes a quadratic data term (theta - mu)^2 * 0.5 * kappa.
// A prior bug multiplied this term by an uninitialized 0, silently dropping it;
// the kappa samplers then ignored how far the random effect sat from its mean.
//
TEST_CASE( "tvarscale posterior_function includes the quadratic data term",
           "[draw_][tvarscale]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  // Known inputs. A small sigma makes the quadratic data term sizeable, so the
  // dropped-term bug is unmistakable.
  const double mu     = 3.5;   // patient mass_mean
  const double sigma  = 1.0;   // patient mass_sd  (sigma^2 = 1)
  const double theta  = 8.0;   // pulse mass (random effect), far from mu
  const double kappa  = 1.0;   // current tvarscale_mass
  const double prop   = 2.0;   // proposed tvarscale

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;

  PulseEstimates pulse;          // empty pulse; posterior_function only reads
  pulse.mass           = theta;  // mass and tvarscale_mass (no mean_contrib)
  pulse.tvarscale_mass = kappa;

  // for_width = false -> sample the mass tvarscale
  SS_DrawTVarScale sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  double got = sampler.posterior_function(&pulse, prop, &pat);

  // Analytic log posterior ratio (same closed form the sampler should compute).
  const double dev2        = (theta - mu) * (theta - mu);
  const double prior_ratio = std::log(Rf_dgamma(prop, 2, 0.5, 0)) -
                             std::log(Rf_dgamma(kappa, 2, 0.5, 0));
  const double stdold      = theta / (sigma / std::sqrt(kappa));
  const double stdnew      = theta / (sigma / std::sqrt(prop));
  const double quad        = (dev2 * 0.5 * kappa - dev2 * 0.5 * prop) / (sigma * sigma);
  const double re_ratio    = quad +
                             Rf_pnorm5(stdold, 0, 1, 1, 1) -
                             Rf_pnorm5(stdnew, 0, 1, 1, 1) -
                             0.5 * std::log(kappa) + 0.5 * std::log(prop);
  const double expected    = prior_ratio + re_ratio;

  SECTION( "matches the analytic scaled-normal log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "regression guard: the quadratic term is present and non-trivial" ) {
    // The quadratic term here is large (~ -10.1). Under the original bug it was
    // identically 0, so `got` would equal `expected - quad`. Assert it does not.
    REQUIRE( quad != Approx(0.0) );
    REQUIRE( got  != Approx(expected - quad) );
  }
}


//
// SS_DrawRandomEffects::posterior_function
//
// Regression test for the pulse random-effect (mass/width) acceptance ratio.
// The random effect is a scale-mixture-of-normals: theta_k ~ N(mu, sigma^2/kappa_k),
// so the prior quadratic term must be weighted by the per-pulse t-scale kappa_k.
// A prior bug omitted kappa entirely, leaving prior variance sigma^2 and making
// the random-effect draw inconsistent with the tvarscale and SD-of-random-effects
// full conditionals (which do use kappa). The inconsistency hit the weakly
// identified pulse *width* hardest, letting its pulse-to-pulse SD space-fill.
//
// The likelihood term does not depend on kappa, so evaluating posterior_function
// at two different kappa values (holding the random effect, mean, sd, and proposal
// fixed) isolates the kappa-weighted prior term: the likelihood contribution
// cancels in the difference. Under the old bug this difference was identically 0.
//
TEST_CASE( "random-effects posterior_function weights the prior by the t-scale kappa",
           "[draw_][randomeffects]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  // The test patient is constructed with a single initial pulse; use it as the
  // random effect being sampled (it must live inside pat.pulses so the internal
  // likelihood() reflects the proposed value).
  REQUIRE( pat.get_pulsecount() >= 1 );
  PulseEstimates & pulse = pat.pulses.front();

  // Known inputs. Small sigma + a random effect far from its mean make the
  // kappa-weighted prior term sizeable, so the dropped-kappa bug is unmistakable.
  const double mu     = 3.5;   // patient mass_mean
  const double sigma  = 1.0;   // patient mass_sd (sigma^2 = 1)
  const double theta  = 8.0;   // pulse mass (random effect), far from mu
  const double prop   = 2.0;   // proposed pulse mass
  const double kappa1 = 1.0;   // first per-pulse t-scale
  const double kappa2 = 2.0;   // second per-pulse t-scale

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;
  pulse.mass              = theta;

  // for_width = false -> sample the pulse mass random effect
  SS_DrawRandomEffects sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  pulse.tvarscale_mass = kappa1;
  double r1 = sampler.posterior_function(&pulse, prop, &pat);

  pulse.tvarscale_mass = kappa2;
  double r2 = sampler.posterior_function(&pulse, prop, &pat);

  // Only the kappa-weighted prior term differs between the two calls:
  //   r(kappa) = kappa * [0.5(theta-mu)^2 - 0.5(prop-mu)^2] / sigma^2 + (loglik term)
  const double prior_unit = (0.5 * (theta - mu) * (theta - mu) -
                             0.5 * (prop  - mu) * (prop  - mu)) / (sigma * sigma);
  const double expected_diff = (kappa2 - kappa1) * prior_unit;

  SECTION( "kappa scales the prior term as predicted by the scaled-normal model" ) {
    REQUIRE( r2 - r1 == Approx(expected_diff) );
  }

  SECTION( "regression guard: the acceptance ratio actually depends on kappa" ) {
    // prior_unit (= 9.0 here) is non-trivial, so a kappa-dependent ratio must
    // differ between kappa1 and kappa2. Under the original bug r2 - r1 == 0.
    REQUIRE( prior_unit != Approx(0.0) );
    REQUIRE( r2 != Approx(r1) );
  }
}


//
// SS_DrawRandomEffects::posterior_function -- log-normal parameterization
//
// Under the papers' log-normal model, log(theta_k) ~ N(mu, sigma^2/kappa_k),
// where mu/sigma are the mean/SD of the LOG random effect. The pulse stores the
// natural-scale value, and the symmetric random-walk proposal is on that natural
// value, so the target of theta carries the log-normal's 1/theta Jacobian; the
// MH log-ratio therefore adds (log(theta) - log(theta')) to the kappa-weighted
// quadratic term. This test checks the full closed form exactly (the likelihood
// delta is computed the same way the sampler does, so what is verified is the
// prior + Jacobian ratio).
//
TEST_CASE( "random-effects posterior_function: log-normal prior + Jacobian",
           "[draw_][randomeffects][lognormal]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses = true;   // papers' parameterization

  REQUIRE( pat.get_pulsecount() >= 1 );
  PulseEstimates & pulse = pat.pulses.front();

  const double mu     = 3.7;   // mean of LOG width (e.g. log(~40))
  const double sigma  = 0.6;   // SD of LOG width
  const double theta  = 55.0;  // current width (natural scale, > 0)
  const double prop   = 30.0;  // proposed width (natural scale, > 0)
  const double kappa  = 1.3;   // per-pulse t-scale

  pat.estimates.width_mean = mu;
  pat.estimates.width_sd   = sigma;
  pulse.width              = theta;
  pulse.tvarscale_width    = kappa;

  // for_width = true -> sample the pulse width random effect
  SS_DrawRandomEffects sampler(1.0, 500, 25000, 0.35, true, false, 5000);
  double got = sampler.posterior_function(&pulse, prop, &pat);

  // Closed-form prior + Jacobian log-ratio.
  const double logold = std::log(theta);
  const double lognew = std::log(prop);
  const double prior_ratio =
      kappa * (0.5*(logold-mu)*(logold-mu) - 0.5*(lognew-mu)*(lognew-mu))
              / (sigma*sigma)
      + (logold - lognew);

  // Likelihood delta, evaluated exactly as the sampler does.
  const double clike = pat.likelihood(false);
  pulse.width = prop;
  const double plike = pat.likelihood(false);
  pulse.width = theta;  // reset

  const double expected = prior_ratio + (plike - clike);

  SECTION( "matches the analytic log-normal log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // The same inputs under the natural-scale prior give a different ratio,
    // confirming the lognormal branch is actually taken.
    pat.lognormal_pulses = false;
    double got_natural = sampler.posterior_function(&pulse, prop, &pat);
    REQUIRE( got != Approx(got_natural) );
  }
}


//
// SS_DrawFixedEffects::posterior_function -- log-normal parameterization
//
// Draws the pulse-population MEAN mu (for mass/width). Under the papers'
// log-normal model, log(theta_k) ~ N(mu, sigma^2/kappa_k), so the per-pulse
// residuals entering mu's full conditional are formed on the LOG scale,
// (log(theta_k) - mu), weighted by the per-pulse t-scale kappa_k. There is no
// truncation at 0, so the truncated-normal normalizing constant (the pnorm
// terms of the natural-scale branch) is dropped. mu is not transformed, so no
// Jacobian appears; its normal hyperprior contributes the usual quadratic ratio.
//
TEST_CASE( "fixedeffects posterior_function: log-normal mean full conditional",
           "[draw_][fixedeffects][lognormal]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses = true;   // papers' parameterization

  REQUIRE( pat.get_pulsecount() == 1 );  // single pulse -> single-term sum
  PulseEstimates & pulse = pat.pulses.front();

  const double mu       = 1.1;   // current mean of LOG mass
  const double prop     = 0.7;   // proposed mean of LOG mass
  const double sigma    = 0.8;   // SD of LOG mass
  const double theta    = 5.0;   // pulse mass (natural scale, > 0)
  const double kappa    = 1.4;   // per-pulse t-scale
  const double pm       = 0.9;   // prior hypermean on mu
  const double pv       = 25.0;  // prior hypervariance on mu

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;
  pat.priors.mass_mean    = pm;
  pat.priors.mass_variance = pv;
  pulse.mass              = theta;
  pulse.tvarscale_mass    = kappa;

  // for_width = false -> sample the mean mass; lognormal = true
  SS_DrawFixedEffects sampler(1.0, 500, 25000, 0.35, false, false, 5000, true);
  bool notused = false;
  double got = sampler.posterior_function(&pat, prop, &notused);

  // Closed form: normal hyperprior ratio on mu + log-scale quadratic ratio,
  // with NO truncation normalizing constant.
  const double logtheta   = std::log(theta);
  const double prior_ratio = (std::pow(mu - pm, 2) - std::pow(prop - pm, 2)) /
                             (2.0 * pv);
  const double psum_old   = std::pow(logtheta - mu,   2) * kappa;
  const double psum_new   = std::pow(logtheta - prop, 2) * kappa;
  const double prop_ratio = 0.5 / (sigma * sigma) * (psum_old - psum_new);
  const double expected   = prior_ratio + prop_ratio;

  SECTION( "matches the analytic log-normal mean log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // A separate natural-scale sampler over the same patient/pulse uses
    // (theta - mu) residuals and keeps the pnorm truncation terms, so the ratio
    // must differ -- confirming the lognormal branch is actually taken.
    SS_DrawFixedEffects sampler_nat(1.0, 500, 25000, 0.35, false, false, 5000, false);
    pat.lognormal_pulses = false;
    double got_natural = sampler_nat.posterior_function(&pat, prop, &notused);
    REQUIRE( got != Approx(got_natural) );
  }
}


//
// SS_DrawSDRandomEffects::posterior_function -- log-normal parameterization
//
// Draws the pulse-to-pulse SD sigma. Under the log-normal model the pulse
// residuals are formed on the LOG scale, (log(theta_k) - mu), weighted by
// kappa_k; there is no truncation at 0, so the truncated-normal normalizing
// constants (old_int/new_int) are dropped. The Gaussian-density n*log(sigma)
// Jacobian, the precision term, and the half-Cauchy prior ratio are unchanged.
//
TEST_CASE( "sdrandomeffects posterior_function: log-normal SD full conditional",
           "[draw_][sdrandomeffects][lognormal]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses = true;   // papers' parameterization

  REQUIRE( pat.get_pulsecount() == 1 );  // single pulse -> single-term sum
  PulseEstimates & pulse = pat.pulses.front();

  const double mu       = 1.1;   // mean of LOG mass
  const double sigma    = 0.8;   // current SD of LOG mass
  const double prop     = 1.3;   // proposed SD of LOG mass
  const double theta    = 5.0;   // pulse mass (natural scale, > 0)
  const double kappa    = 1.4;   // per-pulse t-scale
  const double cauchy   = 5.0;   // half-Cauchy scale parameter

  pat.estimates.mass_mean  = mu;
  pat.estimates.mass_sd    = sigma;
  pat.priors.mass_sd_param = cauchy;
  pulse.mass               = theta;
  pulse.tvarscale_mass     = kappa;

  // for_width = false -> sample the SD of pulse masses
  SS_DrawSDRandomEffects sampler(1.0, 500, 25000, 0.35, false, false, 5000);
  double got = sampler.posterior_function(&pat, prop, &pat);

  // Closed form (log-normal): no truncation constant.
  const double logtheta    = std::log(theta);
  const double third_part  = kappa * (logtheta - mu) * (logtheta - mu);
  const double first_part  = pat.get_pulsecount() * (std::log(sigma) - std::log(prop));
  const double second_part = 0.5 * ((1.0 / (sigma * sigma)) - (1.0 / (prop * prop)));
  const double fourth_part = std::log(cauchy + sigma * sigma) -
                             std::log(cauchy + prop * prop);
  const double expected    = first_part + second_part * third_part + fourth_part;

  SECTION( "matches the analytic log-normal SD log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // Natural-scale uses (theta - mu) residuals and keeps the pnorm truncation
    // terms, so the ratio must differ -- confirming the lognormal branch is taken.
    pat.lognormal_pulses = false;
    double got_natural = sampler.posterior_function(&pat, prop, &pat);
    REQUIRE( got != Approx(got_natural) );
  }
}


//
// SS_DrawSDRandomEffects -- Uniform(0, max) SD prior (papers' default)
//
// Draws the pulse-to-pulse SD sigma. This axis is ORTHOGONAL to the RE scale
// (lognormal_pulses): switching the SD prior from half-Cauchy to Uniform(0, max)
// changes ONLY the prior ratio (the "fourth part") and the support bound; the
// likelihood/precision/Jacobian terms and the log-vs-natural residual geometry
// are untouched.
//
// Uniform(0, max) has constant density 1/max on its support, and the MH base
// rejects out-of-support proposals via parameter_support before ever calling
// posterior_function, so both the current SD and the proposal are guaranteed in
// (0, max). The prior log-ratio is therefore log(1/max) - log(1/max) = 0, i.e.
// the uniform acceptance ratio equals the half-Cauchy ratio with its
// fourth_part = log(c + sigma^2) - log(c + prop^2) term removed.
//
TEST_CASE( "sdrandomeffects Uniform SD prior: prior ratio cancels, support bounded",
           "[draw_][sdrandomeffects][uniform]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  REQUIRE( pat.get_pulsecount() == 1 );  // single pulse -> single-term sum
  PulseEstimates & pulse = pat.pulses.front();

  const double mu       = 3.5;   // patient mass_mean
  const double sigma    = 0.8;   // current SD of pulse masses
  const double prop     = 1.3;   // proposed SD of pulse masses (in (0, max))
  const double theta    = 5.0;   // pulse mass (natural scale, > 0)
  const double kappa    = 1.4;   // per-pulse t-scale
  const double cauchy   = 5.0;   // half-Cauchy scale parameter
  const double sdmax    = 3.0;   // Uniform(0, max) upper bound

  pat.estimates.mass_mean  = mu;
  pat.estimates.mass_sd    = sigma;
  pat.priors.mass_sd_param = cauchy;
  pat.priors.mass_sd_max   = sdmax;
  pulse.mass               = theta;
  pulse.tvarscale_mass     = kappa;

  // for_width = false -> sample the SD of pulse masses
  SS_DrawSDRandomEffects sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  // The half-Cauchy prior ratio (fourth_part) is what the uniform branch drops.
  const double fourth_part = std::log(cauchy + sigma * sigma) -
                             std::log(cauchy + prop * prop);

  // The uniform ratio must equal the half-Cauchy ratio minus fourth_part, in
  // BOTH RE-scale modes (the two axes are independent).
  for (bool lognormal : { false, true }) {
    pat.lognormal_pulses = lognormal;

    pat.uniform_sd_prior = false;
    const double got_cauchy = sampler.posterior_function(&pat, prop, &pat);

    pat.uniform_sd_prior = true;
    const double got_uniform = sampler.posterior_function(&pat, prop, &pat);

    // Uniform prior ratio is 0, so it equals the half-Cauchy ratio less fourth_part.
    REQUIRE( got_uniform == Approx(got_cauchy - fourth_part) );

    // Guard: the two priors give genuinely different ratios (fourth_part != 0).
    REQUIRE( fourth_part != Approx(0.0) );
    REQUIRE( got_uniform != Approx(got_cauchy) );
  }

  SECTION( "parameter_support enforces the open interval (0, max)" ) {
    pat.uniform_sd_prior = true;
    // Inside the interval: accepted.
    REQUIRE( sampler.parameter_support(0.5 * sdmax, &pat) );
    REQUIRE( sampler.parameter_support(0.999 * sdmax, &pat) );
    // At/above the upper bound: rejected.
    REQUIRE_FALSE( sampler.parameter_support(sdmax, &pat) );
    REQUIRE_FALSE( sampler.parameter_support(sdmax + 1.0, &pat) );
    // Non-positive: rejected.
    REQUIRE_FALSE( sampler.parameter_support(0.0, &pat) );
    REQUIRE_FALSE( sampler.parameter_support(-0.1, &pat) );

    // Half-Cauchy support is positive-only, unbounded above (max not consulted).
    pat.uniform_sd_prior = false;
    REQUIRE( sampler.parameter_support(sdmax + 1.0, &pat) );
    REQUIRE_FALSE( sampler.parameter_support(0.0, &pat) );
  }
}


//
// SS_DrawTVarScale::posterior_function -- log-normal parameterization
//
// The per-pulse t-scale (kappa) update under the papers' log-normal model:
// log(theta) ~ N(mu, sigma^2/kappa), with mu/sigma the mean/SD of the LOG random
// effect. The quadratic deviation is formed on the log scale, (log(theta) - mu),
// and there is no truncation at 0, so the two pnorm normalizing constants of the
// natural-scale branch are dropped. The sqrt(kappa) precision normalizing
// constant (-0.5 log kappa + 0.5 log proposal) is retained, and the 1/theta
// Jacobian is constant in kappa (theta fixed) so it cancels. This checks the
// exact closed form and guards that the lognormal branch is actually taken.
//
TEST_CASE( "tvarscale posterior_function: log-normal kappa update",
           "[draw_][tvarscale][lognormal]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses = true;   // papers' parameterization

  const double mu     = 1.2;   // mean of LOG mass
  const double sigma  = 0.5;   // SD of LOG mass
  const double theta  = 5.0;   // pulse mass (natural scale, > 0)
  const double kappa  = 1.3;   // current tvarscale_mass
  const double prop   = 2.0;   // proposed tvarscale

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;

  PulseEstimates pulse;          // empty pulse; posterior_function only reads
  pulse.mass           = theta;  // mass and tvarscale_mass
  pulse.tvarscale_mass = kappa;

  // for_width = false -> sample the mass tvarscale
  SS_DrawTVarScale sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  double got = sampler.posterior_function(&pulse, prop, &pat);

  // Analytic log-normal log-ratio (same closed form the sampler computes).
  const double prior_ratio = std::log(Rf_dgamma(prop, 2, 0.5, 0)) -
                             std::log(Rf_dgamma(kappa, 2, 0.5, 0));
  const double logre       = std::log(theta);
  const double dev2        = (logre - mu) * (logre - mu);
  const double re_ratio    = (dev2 * 0.5 * kappa - dev2 * 0.5 * prop) /
                             (sigma * sigma)
                             - 0.5 * std::log(kappa) + 0.5 * std::log(prop);
  const double expected    = prior_ratio + re_ratio;

  SECTION( "matches the analytic log-normal log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // The same inputs under the natural-scale branch use (theta - mu) residuals
    // and keep the two pnorm truncation terms, so the ratio must differ --
    // confirming the lognormal branch is actually taken.
    pat.lognormal_pulses = false;
    double got_natural = sampler.posterior_function(&pulse, prop, &pat);
    REQUIRE( got != Approx(got_natural) );
  }
}


//
// BirthDeathProcess::add_new_pulse -- log-normal draw path
//
// Under lognormal_pulses the new pulse's mass/width are drawn as exp(rnorm(mu,
// sigma/sqrt(kappa))). With gaussian_random_effects the per-pulse kappa is fixed
// at 1 (no rgamma draw), so add_new_pulse consumes exactly two rnorm draws --
// letting a separately-seeded reference reproduce them exactly. This verifies
// the draw equals exp() of the underlying normal and is strictly positive.
//
TEST_CASE( "birth-death add_new_pulse: log-normal draw is exp(normal) and positive",
           "[birthdeath][lognormal]" ) {

  DataStructuresUtils utils;
  PulseUtils pu;
  BirthDeathProcess bd;

  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses        = true;   // draw on the log scale, exponentiate
  pat.gaussian_random_effects = true;   // kappa fixed at 1 -> exactly two rnorms
  pat.estimates.mass_mean  = 1.2;       // mean of LOG mass
  pat.estimates.mass_sd    = 0.5;
  pat.estimates.width_mean = 3.0;       // mean of LOG width
  pat.estimates.width_sd   = 0.7;

  const int n0 = pat.get_pulsecount();

  // Reference: reproduce the exact two draws add_new_pulse makes under the
  // gaussian-kappa log-normal branch (kappa = 1, so no rgamma is drawn).
  double n_mass = 0.0, n_width = 0.0;
  {
    pu.set_seed(20240706);
    Rcpp::RNGScope scope;
    n_mass  = Rf_rnorm(1.2, 0.5);
    n_width = Rf_rnorm(3.0, 0.7);
  }

  pu.set_seed(20240706);
  bd.add_new_pulse(&pat, 500.0);

  REQUIRE( pat.get_pulsecount() == n0 + 1 );
  PulseEstimates & p = pat.pulses.back();

  SECTION( "mass/width are exp() of the underlying normal draws" ) {
    REQUIRE( p.mass  == Approx(std::exp(n_mass)) );
    REQUIRE( p.width == Approx(std::exp(n_width)) );
  }

  SECTION( "log-normal draws are strictly positive; gaussian kappa is fixed at 1" ) {
    REQUIRE( p.mass  > 0.0 );
    REQUIRE( p.width > 0.0 );
    REQUIRE( p.tvarscale_mass  == 1.0 );
    REQUIRE( p.tvarscale_width == 1.0 );
  }
}
