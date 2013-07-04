#include <stdlib.h>
#include <iostream>
#include <list>

struct Particle{
 double mass;
 double x;
 double y;
 double z;
 double v_x;
 double v_y;
 double v_z;
 
 Particle(double _mass=0.0, double _x=0.0, double _y=0.0, double _z =0.0, double _v_x = 0.0, double _v_y = 0.0, double _v_z = 0.0):
  mass(_mass),
  x(_x),
  y(_y),
  z(_z),
  v_x(_v_x),
  v_y(_v_y),
  v_z(_v_z){}
 
};

class MolDy
{
  private:
  
  int vis_space;
  double t_start;
  double t_end;
  double delta_t;
  double x_min;
  double y_min;
  double z_min;
  double x_max;
  double y_max;
  double z_max;
  double r_cut;
  double epsilon;
  double sigma;
  
  double cell_size_x;
  double cell_size_y;
  double cell_size_z;
  
  std::list<Particle*> cell;
  
  public:
  
  MolDy(void){}
  
  ~MolDy(void){}
  
  void readParameters(std::string filename);
  void testPrintParameters(void);
  void readData(std::string filename);
  void testPrintData(void);
  
};

