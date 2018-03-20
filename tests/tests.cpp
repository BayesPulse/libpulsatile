// tests.cpp
//
// file dedicated to providing catch's main()
//

// Let Catch provide main():
//#define CATCH_CONFIG_MAIN
// Let Catch handle setup/cleanup, but use own main()
#define CATCH_CONFIG_RUNNER
#include <RInside.h>
#include <testing/catch.h>

int main( int argc, char* argv[] ) {

  // Create one R instance for all testing
  RInside R;

  // Let Catch do the rest
  int result = Catch::Session().run( argc, argv );

  return result;
}

