#include "Generator.h"



int main(int argc, char** argv) {

  int nevents = 10;
  Generator generator;
  generator.simulate(nevents);
  return 0;
}
