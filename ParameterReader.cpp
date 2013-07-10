#include "ParameterReader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

bool ParameterReader::readParameters(const std::string& filename)
{
    std::ifstream infile(filename.c_str());
  if (infile.good()){
      while (!infile.eof()){
	  std::string line, param;
	  std::stringstream ss;
	  getline(infile,line);
	  ss<<line;
	  ss>>param;
	  if (param == "name"){
	      ss>>vtkfile;
	      
	      if (ss.fail()){
		  std::cerr<<"Error reading 'vis_space' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      vtkfile_defined = 1;
	  }
	  if (param == "vis_space"){
	      ss>>vis_space;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'vis_space' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      vis_space_defined = 1;
	  }
	  else if (param == "t_start"){
	      ss>>t_start;
	      if (ss.fail()){
		  std::cerr<<"Error reading 't_start' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      t_start_defined = 1;
	  }
	  else if (param == "t_end"){
	      ss>>t_end;
	      if (ss.fail()){
		  std::cerr<<"Error reading 't_end' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      t_end_defined = 1;
	  }
	  else if (param == "delta_t"){
	      ss>>delta_t;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'delta_t' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      delta_t_defined = 1;
	  }
	  else if (param == "x_min"){
	      ss>>x_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'x_min' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      x_min_defined = 1;
	  }
	  else if (param == "y_min"){
	      ss>>y_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'y_min' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      y_min_defined = 1;
	  }
	  else if (param == "z_min"){
	      ss>>z_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'z_min' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      z_min_defined = 1;
	  }
	  else if (param == "x_max"){
	      ss>>x_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'x_max' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      x_max_defined = 1;
	  }
	  else if (param == "y_max"){
	      ss>>y_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'y_max' in '"<<filename<<"'"<<std::endl;
		  return(0);;
	      }
	      y_max_defined = 1;
	  }
	  else if (param == "z_max"){
	      ss>>z_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'z_max' in '"<<filename<<"'"<<std::endl;
		  return(0);;
	      }
	      z_max_defined = 1;
	  }
	  else if (param == "r_cut"){
	      ss>>r_cut;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'r_cut' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      r_cut_defined = 1;
	  }
	  else if (param == "epsilon"){
	      ss>>epsilon;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'epsilon' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      epsilon_defined = 1;
	  }
	  else if (param == "sigma"){
	      ss>>sigma;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'sigma' in '"<<filename<<"'"<<std::endl;
		  return(0);
	      }
	      sigma_defined = 1;
	  }
  // 	  else if (param == "gravity"){
  // 	      ss>>gravity;
  // 	      if (ss.fail()){
  // 		  std::cerr<<"Error reading 'sigma' in '"<<filename<<"'"<<std::endl;
  // 		  return exit(-1);
  // 	      }
  // 	  }
  // 	  else if (param == "reflective_boundaries"){
  // 	      ss>>reflective_boundaries;
  // 	      if (ss.fail()){
  // 		  std::cerr<<"Error reading 'sigma' in '"<<filename<<"'"<<std::endl;
  // 		  return exit(-1);
  // 	      }
  // 	  }
      }
      infile.close(); 
  }
  else{
      std::cerr<<"Error reading '"<<filename<<"'"<<std::endl;
      return 0;
  }
   
  return(1);
}

