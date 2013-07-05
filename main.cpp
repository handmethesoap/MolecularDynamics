#include "MolDy.h"
#include <iostream>

int main(int argc, char** argv)
{
  if (argc!=3){
    std::cerr<<"There should be two arguments"<<std::endl;
    return -1;
  }
  
  std::string parameterfile = argv[1];
  std::string datafile = argv[2];
  
  MolDy particleSystem;
  
  particleSystem.readParameters(parameterfile);
  particleSystem.testPrintParameters();
  particleSystem.readData(datafile);
  particleSystem.simulate();
  return 0;
}