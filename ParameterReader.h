#ifndef PR
#define PR

#include <sstream>
#include <stdlib.h>
#include <iostream>

class ParameterReader
{
private:
  std::string vtkfile;
  std::string vis_space;
  std::string t_start;
  std::string t_end;
  std::string delta_t;
  std::string x_min;
  std::string y_min;
  std::string z_min;
  std::string x_max;
  std::string y_max;
  std::string z_max;
  std::string r_cut;
  std::string epsilon;
  std::string sigma;
  
  bool vtkfile_defined;
  bool vis_space_defined;
  bool t_start_defined;
  bool t_end_defined;
  bool delta_t_defined;
  bool x_min_defined;
  bool y_min_defined;
  bool z_min_defined;
  bool x_max_defined;
  bool y_max_defined;
  bool z_max_defined;
  bool r_cut_defined;
  bool epsilon_defined;
  bool sigma_defined;
  
public:
  ParameterReader(void){
    vtkfile_defined = 0;
    vis_space_defined = 0;
    t_start_defined = 0;
    t_end_defined = 0;
    delta_t_defined = 0;
    x_min_defined = 0;
    y_min_defined = 0;
    z_min_defined = 0;
    x_max_defined = 0;
    y_max_defined = 0;
    z_max_defined = 0;
    r_cut_defined = 0;
    epsilon_defined = 0;
    sigma_defined = 0;
  }
  ~ParameterReader(void){
    
  }
  bool readParameters(const std::string& filename);
  inline bool IsDefined(const std::string& key) const
  {
  if (key == "name"){
	return vtkfile_defined;
    }
    else if (key == "vis_space"){
	return vis_space_defined;
    }
    else if (key == "t_start"){
	return t_start_defined;
    }
    else if (key == "t_end"){
	return t_end_defined;
    }
    else if (key == "delta_t"){
	return delta_t_defined;
    }
    else if (key == "x_min"){
	return x_min_defined;
    }
    else if (key == "y_min"){
	return y_min_defined;
    }
    else if (key == "z_min"){
	return z_min_defined;
    }
    else if (key == "x_max"){
	return x_max_defined;
    }
    else if (key == "y_max"){
	return y_max_defined;
    }
    else if (key == "z_max"){
	return z_max_defined;
    }
    else if (key == "r_cut"){
	return r_cut_defined;
    }
    else if (key == "epsilon"){
	return epsilon_defined;
    }
    else if (key == "sigma"){
	return sigma_defined;
    }
    else
    {
      return 0;
    }
  
  }

  template<typename Type>
  inline void GetParameter(const std::string& key, Type &value) const
  {
    std::stringstream ss;
  //   ss>>key;
    if (key == "name"){
	ss << vtkfile;
	ss >> value;
    }
    else if (key == "vis_space"){
	ss << vis_space;
	ss >> value;
    }
    else if (key == "t_start"){
	ss << t_start;
	ss >> value;
    }
    else if (key == "t_end"){
	ss << t_end;
	ss >> value;
    }
    else if (key == "delta_t"){
	ss << delta_t;
	ss >> value;
    }
    else if (key == "x_min"){
	ss << x_min;
	ss >> value;
    }
    else if (key == "y_min"){
	ss << y_min;
	ss >> value;
    }
    else if (key == "z_min"){
	ss << z_min;
	ss >> value;
    }
    else if (key == "x_max"){
	ss << x_max;
	ss >> value;
    }
    else if (key == "y_max"){
	ss << y_max;
	ss >> value;
    }
    else if (key == "z_max"){
	ss << z_max;
	ss >> value;
    }
    else if (key == "r_cut"){
	ss << r_cut;
	ss >> value;
    }
    else if (key == "epsilon"){
	ss << epsilon;
	ss >> value;
    }
    else if (key == "sigma"){
	ss << sigma;
	ss >> value;
    }
  }
};

#endif
