#include <stdlib.h>
#include <iostream>
#include <list>
#include <vector>

struct Particle{
 double mass;
 double x;
 double y;
 double z;
 double v_x;
 double v_y;
 double v_z;
 double f_x;
 double f_y;
 double f_z;
 double f_x_prev;
 double f_y_prev;
 double f_z_prev;
 
 Particle(double _mass=0.0, double _x=0.0, double _y=0.0, double _z =0.0, double _v_x = 0.0, double _v_y = 0.0, double _v_z = 0.0, double _f_x = 0.0, double _f_y = 0.0, double _f_z = 0.0, double _f_x_prev = 0.0, double _f_y_prev = 0.0, double _f_z_prev = 0.0):
  mass(_mass),
  x(_x),
  y(_y),
  z(_z),
  v_x(_v_x),
  v_y(_v_y),
  v_z(_v_z),
  f_x(_f_x),
  f_y(_f_y),
  f_z(_f_z),
  f_x_prev(_f_x_prev),
  f_y_prev(_f_y_prev),
  f_z_prev(_f_z_prev){}
 
};

class MolDy
{
  private:
  
  std::string vtkfile;
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
  
  int cells_x;
  int cells_y;
  int cells_z;
  
  int number_particles;
  
  std::list<Particle*> cell;
  std::vector< std::list<Particle*>* > cells;
  
  public:
  
  MolDy(void){}
  
  ~MolDy(void){
    for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
    {
      if(!(*it)->empty())
      {
	for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
	{
	  delete (*it2);
	}
      }
      delete (*it);
    }
  }
  
  void readParameters(std::string filename);
  void testPrintParameters(void);
  void readData(std::string filename);
  void testPrintParticles(void);
  
  void assignParticle(Particle* particle);
  
  void VTKPrint(std::string filename);
  
  void updatePosition(void);
  void calculateForces(void);
  void updateVelocities(void);
  
  void updateCell(void);
  void checkCells(void);
  int adjacentCell(int x_offset, int y_offset, int z_offset, int vector_index);
  void intercellForces(std::list<Particle*>::iterator particle, std::vector< std::list<Particle*>* >::iterator adjacent_cell);
  
  void simulate(void);
  
  
};

