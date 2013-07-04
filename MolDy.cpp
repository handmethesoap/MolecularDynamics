#include "MolDy.h"
#include <fstream>
#include <sstream>
#include <iostream>

void MolDy::readParameters(std::string filename)
{
  std::ifstream infile(filename.c_str());
  if (infile.good()){
      while (!infile.eof()){
	  std::string line, param;
	  std::stringstream ss;
	  getline(infile,line);
	  ss<<line;
	  ss>>param;
	  if (param == "vis_space"){
	      ss>>vis_space;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'vis_space' in '"<<filename<<"'"<<std::endl;
		  exit(-1);
	      }
	  }
	  else if (param == "t_start"){
	      ss>>t_start;
	      if (ss.fail()){
		  std::cerr<<"Error reading 't_start' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "t_end"){
	      ss>>t_end;
	      if (ss.fail()){
		  std::cerr<<"Error reading 't_end' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "delta_t"){
	      ss>>delta_t;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'delta_t' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "x_min"){
	      ss>>x_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'x_min' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "y_min"){
	      ss>>y_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'y_min' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "z_min"){
	      ss>>z_min;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'z_min' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "x_max"){
	      ss>>x_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'x_max' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "y_max"){
	      ss>>y_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'y_max' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "z_max"){
	      ss>>z_max;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'z_max' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "r_cut"){
	      ss>>r_cut;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'r_cut' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "epsilon"){
	      ss>>epsilon;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'epsilon' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
	  else if (param == "sigma"){
	      ss>>sigma;
	      if (ss.fail()){
		  std::cerr<<"Error reading 'sigma' in '"<<filename<<"'"<<std::endl;
		  return exit(-1);
	      }
	  }
      }
      infile.close();
  }
  else{
      std::cerr<<"Error reading '"<<filename<<"'"<<std::endl;
      return exit(-1);
  }
}

void MolDy::testPrintParameters(void)
{
  std::cout<< vis_space << std::endl;
  std::cout<< t_start << std::endl;
  std::cout<< t_end << std::endl;
  std::cout<< delta_t << std::endl;
  std::cout<< x_min << std::endl;
  std::cout<< y_min << std::endl;
  std::cout<< z_min << std::endl;
  std::cout<< x_max << std::endl;
  std::cout<< y_max << std::endl;
  std::cout<< z_max << std::endl;
  std::cout<< r_cut << std::endl;
  std::cout<< epsilon << std::endl;
  std::cout<< sigma << std::endl;
  int temp = x_max/r_cut;
  double temp2 = x_max/temp;
  std::cout<< temp << std::endl;
  std::cout<< temp2 << std::endl;
  std::cout<< temp*temp2 << std::endl;
}

void MolDy::readData(std::string filename)
{
  double _mass;
  double _x;
  double _y;
  double _z;
  double _v_x;
  double _v_y;
  double _v_z;
  
  cell_size_x = (x_max - x_min)/(static_cast<int>((x_max - x_min)/r_cut));
  cell_size_y = (y_max - y_min)/(static_cast<int>((y_max - y_min)/r_cut));
  cell_size_z = (y_max - y_min)/(static_cast<int>((z_max - z_min)/r_cut));
  
#ifdef DEBUG
  std::cout << "cell dimensions (x,y,z) = (" << cell_size_x << "," << cell_size_y << "," << cell_size_z << ")" << std::endl;
#endif
  
  std::ifstream infile(filename.c_str());

  if (infile.is_open()){
      while (infile.good()){
	  std::string line, param;
	  std::stringstream ss;
	  getline(infile,line);
	  
	  if(line == "") continue;
	  
	  ss<<line;
	  ss>>_mass;
	  ss>>_x;
	  ss>>_y;
	  ss>>_z;
	  ss>>_v_x;
	  ss>>_v_y;
	  ss>>_v_z;
	  
#ifdef DEBUG
	  std::cout << _mass << " " << _x << " " << _y << " " << _z << " " << _v_x << " " << _v_y << " " << _v_z << std::endl;
#endif
	  
	  cell.push_back(new Particle(_mass, _x, _y, _z, _v_x, _v_y, _v_z));

      }
      infile.close();
  }
  else{
      std::cerr<<"Error reading '"<<filename<<"'"<<std::endl;
      return exit(-1);
  }  
  
}

void MolDy::testPrintData(void)
{
  for(std::list<Particle*>::iterator it = cell.begin(); it != cell.end(); ++it)
  {
    std::cout << "mass = " << (*it)->mass << std::endl;
  }
}