// #include <RcppArmadillo.h>
// #ifndef NORINSIDE
// #include <RInside.h>
// #endif
// #include <bp_datastructures/chains.h>
// #include <bp_datastructures/utils.h>
// #include <testing/catch.h>
// 
// 
// //
// // chains_tests.cpp
// //     Test chain class & functions
// //
// 
// //
// 
// TEST_CASE( "first chains test -- single subject, single hormone", "[chains]" ) {
// 
//   // Create patient object -- using the default test dataset/specs
//   DataStructuresUtils utils;
//   Patient pat = utils.create_new_test_patient_obj();
//   Patient * patient = &pat;
//   patient = utils.add_default_pulses(patient);
// 
//   // Single subject constructor
//   Chains chains(100, 2, 10, false);
// 
//   SECTION( "Check constructor creates empty chain objects" ) {
//     REQUIRE( chains.patient_chain.n_rows == 45 );
//     REQUIRE( chains.patient_chain.n_cols == 9 );
//     REQUIRE( std::abs(chains.patient_chain.max() - chains.patient_chain.min()) < 0.000001 );
//     REQUIRE( chains.patient_chain(11, 1) == 0.0 );
//     REQUIRE( chains.patient_chain(12, 1) == 0.0 );
//   }
// 
//   chains.save_sample(&pat, 30);
//   chains.save_sample(&pat, 32);
//   chains.save_sample(&pat, 52);
// 
//   SECTION( "Check save_sample() function populates common chain" ) {
//     REQUIRE( std::abs(chains.patient_chain.max() - chains.patient_chain.min()) > 0.000001 );
//     REQUIRE( chains.patient_chain(10, 0) == 30 );
//     REQUIRE( chains.patient_chain(10, 1) == 11 );
//     REQUIRE( chains.patient_chain(10, 2) == 2.6 );
//     REQUIRE( chains.patient_chain(10, 3) == 3.5 );
//     REQUIRE( chains.patient_chain(10, 4) == 30 );
//     REQUIRE( chains.patient_chain(10, 5) == 45 );
//     REQUIRE( chains.patient_chain(10, 6) == 0.05 );
//     REQUIRE( chains.patient_chain(10, 7) == 10 );
//     REQUIRE( chains.patient_chain(10, 8) == 10 );
//     REQUIRE( chains.patient_chain(11, 0) == 32 );
//     REQUIRE( chains.patient_chain(11, 1) == 11 );
//     REQUIRE( chains.patient_chain(21, 0) == 52 );
//     REQUIRE( chains.patient_chain(21, 1) == 11 );
//   }
// 
//   //std::cout << "chain 2 is not empty " << chains.pulse_chains[2] << std::endl;
//   //std::cout << "chain 3 is empty " << (chains.pulse_chains[3] == chains.pulse_chains.end()) << std::endl;
//   SECTION( "Check save_sample() function populates pulse chain vector of matrices" ) {
//     REQUIRE( chains.pulse_chains[0].max() - chains.pulse_chains[0].min() > 1000 );
//     REQUIRE( chains.pulse_chains[1].max() - chains.pulse_chains[1].min() > 1000);
//     REQUIRE( chains.pulse_chains[2].max() - chains.pulse_chains[2].min() > 1000 );
//     REQUIRE( &chains.pulse_chains[2] == &chains.pulse_chains.back() );
//   }
// 
// }
// 
// 
