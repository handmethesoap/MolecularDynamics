#include "MolDy.h"
#include <iostream>
#include "ParameterReader.h"

int main(int argc, char** argv)
{
  if (argc!=3){
    std::cerr<<"There should be two arguments"<<std::endl;
    return -1;
  }
  
  std::string parameterfile = argv[1];
  std::string datafile = argv[2];
  
  ParameterReader parameters;
  parameters.readParameters(parameterfile);
  
  MolDy particleSystem(parameters);
  
  particleSystem.readData(datafile);
  particleSystem.simulate();
  return 0;
}