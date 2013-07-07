#include "MolDy.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

void MolDy::readParameters(std::string filename){
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
		  exit(-1);
	      }
	  }
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

void MolDy::testPrintParameters(void){
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

void MolDy::readData(std::string filename){
  double _mass;
  double _x;
  double _y;
  double _z;
  double _v_x;
  double _v_y;
  double _v_z;
  
  cells_x =(static_cast<int>((x_max - x_min)/r_cut));
  cells_y =(static_cast<int>((y_max - y_min)/r_cut));
  cells_z =(static_cast<int>((z_max - z_min)/r_cut));
  cell_size_x = (x_max - x_min)/cells_x;
  cell_size_y = (y_max - y_min)/cells_y;
  cell_size_z = (y_max - y_min)/cells_z;
  
  number_particles = 0;
  
  for(int i = 0; i < cells_x*cells_y*cells_z; ++i)
  {
    cells.push_back(new std::list<Particle*>);
  }
  
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
	  
	  ++number_particles;
	  
	  assignParticle(new Particle(_mass, _x, _y, _z, _v_x, _v_y, _v_z));

      }
      infile.close();
  }
  else{
      std::cerr<<"Error reading '"<<filename<<"'"<<std::endl;
      return exit(-1);
  }  
  
}

void MolDy::assignParticle(Particle* particle){
  
  //deal with periodic boundaries
  if(particle->x >= x_max)
    particle->x = particle->x - (x_max - x_min);
  if(particle->x < x_min)
    particle->x = particle->x + (x_max-x_min);
  if(particle->y >= y_max)
    particle->y = particle->y - (y_max - y_min);
  if(particle->y < y_min)
    particle->y = particle->y + (y_max-y_min);
  if(particle->z >= z_max)
    particle->z = particle->z - (z_max - z_min);
  if(particle->z < z_min)
    particle->z = particle->z + (z_max-z_min);
  
  //find coordinates of cell for particle
  int x = particle->x/cell_size_x;
  int y = particle->y/cell_size_y;
  int z = particle->z/cell_size_z;
  
  //convert coordinates into vector index
  int cellnumber = x + y*cells_x + z*cells_x*cells_y;
  //std::cout << "in " << particle->x << " " << particle->y << " " << particle->z << std::endl;
  //assign particle to cell
  cells[cellnumber]->push_back(particle);
  //std::cout << "out" << std::endl;
}
   
void MolDy::testPrintParticles(void){
  int distance;
  int x;
  int y;
  int z;
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!(*it)->empty())
    {
      distance =  std::distance(cells.begin(), it);
      x = distance%cells_x;
      y =  (distance%(cells_x*cells_y))/cells_x;
      z = distance/(cells_x*cells_y);
      
      std::cout << "CELL: " << distance << ", (x,y,z): (" << x << "," << y << "," << z << "), " << x*cell_size_x << "<x<" << (x+1)*cell_size_x << ", " << y*cell_size_y << "<y<" << (y+1)*cell_size_y << ", " << z*cell_size_z << "<z<" << (z+1)*cell_size_x << std::endl;
      
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	std::cout << "-------------- mass = " << (*it2)->mass << ", (x,y,z) = (" << (*it2)->x << "," << (*it2)->y << "," << (*it2)->z << "), (vx,vy,vz) = (" << (*it2)->v_x << "," << (*it2)->v_y << "," << (*it2)->v_z << ")" << std::endl;
      }
    }
  }
}

void MolDy::VTKPrint(std::string filename){
    std::ofstream outfile(filename.c_str());
    if (outfile.good()){
        outfile<<"# vtk DataFile Version 3.0"<<std::endl;
        outfile<<"SiwiRVisFile"<<std::endl;
        outfile<<"ASCII"<<std::endl;
        outfile<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
        outfile<<"POINTS "<<number_particles<<" DOUBLE" << std::endl;
	for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
	{
	  if(!(*it)->empty())
	  {
	    for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
	    {
	      outfile<< (*it2)->x << " " << (*it2)->y << " " << (*it2)->z << std::endl;
	    }
	  }
	}
	outfile<<"POINT_DATA "<< number_particles << std::endl;
	outfile<< "SCALARS mass double" << std::endl;
	outfile<< "LOOKUP_TABLE default" << std::endl;
	for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
	{
	  if(!(*it)->empty())
	  {
	    for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
	    {
	      outfile<< (*it2)->mass << std::endl;
	    }
	  }
	}
	outfile << "VECTORS force double" << std::endl;
	for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
	{
	  if(!(*it)->empty())
	  {
	    for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
	    {
	      outfile<< (*it2)->f_x << " " << (*it2)->f_y << " " << (*it2)->f_z << std::endl;
	    }
	  }
	}
	outfile << "VECTORS velocity double" << std::endl;
	for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
	{
	  if(!(*it)->empty())
	  {
	    for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
	    {
	      outfile<< (*it2)->v_x << " " << (*it2)->v_y << " " << (*it2)->v_z << std::endl;
	    }
	  }
	}
        outfile.close();
    }
}

void MolDy::updatePosition(void){
  //loop through cells
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!((*it)->empty()))
    { 
      //if cell not empty loop though contents updating positions
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	(*it2)->x = (*it2)->x + delta_t*(*it2)->v_x + (((*it2)->f_x)*delta_t*delta_t)/(2*((*it2)->mass));
	(*it2)->y = (*it2)->y + delta_t*(*it2)->v_y + (((*it2)->f_y)*delta_t*delta_t)/(2*((*it2)->mass)); 
	(*it2)->z = (*it2)->z + delta_t*(*it2)->v_z + (((*it2)->f_z)*delta_t*delta_t)/(2*((*it2)->mass));
      }
    }
  }
  
  updateCell();
}

void MolDy::updateCell(void){
  int distance;
  int cell_x;
  int cell_y;
  int cell_z;
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!(*it)->empty())
    {
      //determine x, y, z coordinate of cell containing the particle
      distance =  std::distance(cells.begin(), it);
      cell_x = distance%cells_x;
      cell_y =  (distance%(cells_x*cells_y))/cells_x;
      cell_z = distance/(cells_x*cells_y);
      
      //loop thtough cell contents
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end();)
      {
	//if the particle has left the cell reassign it and erase it from the current cell
	if((*it2)->x < cell_x*cell_size_x || (*it2)->x >= (cell_x+1)*cell_size_x || (*it2)->y < cell_y*cell_size_y || (*it2)->y >= (cell_y+1)*cell_size_y || (*it2)->z < cell_z*cell_size_z || (*it2)->z >= (cell_z+1)*cell_size_z)
	{
	  assignParticle(*it2);
	  it2 = (*it)->erase(it2);
	}
	else
	{
	  ++it2;
	}
      }
    } 
  }
}

void MolDy::checkCells(void){
  int distance;
  int cell_x;
  int cell_y;
  int cell_z;
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!(*it)->empty())
    {
      distance =  std::distance(cells.begin(), it);
      cell_x = distance%cells_x;
      cell_y =  (distance%(cells_x*cells_y))/cells_x;
      cell_z = distance/(cells_x*cells_y);
     
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	if(((cell_x*cell_size_x) > ((*it2)->x)) || (((*it2)->x) >= ((cell_x +1)*cell_size_x)))
	  std::cout << "particle at x position (" << (*it2)->x << ") is in cell with bounds [" << cell_x*cell_size_x << "," << (cell_x+1)*cell_size_x << "[" << std::endl;
	if(cell_y*cell_size_y > (*it2)->y || (*it2)->y >= (cell_y +1)*cell_size_y)
	  std::cout << "particle at y position (" << (*it2)->y << ") is in cell with bounds [" << cell_y*cell_size_y << "," << (cell_y+1)*cell_size_y << "[" << std::endl;
	if(cell_z*cell_size_z > (*it2)->z || (*it2)->z >= (cell_z +1)*cell_size_z)
	  std::cout << "particle at z position (" << (*it2)->z << ") is in cell with bounds [" << cell_z*cell_size_z << "," << (cell_z+1)*cell_size_z << "[" << std::endl;   
	
      }
    }
  }
}

void MolDy::calculateForces(void){
  
  double xDistance;
  double yDistance;
  double zDistance;
  double xForce;
  double yForce;
  double zForce;
  double distance;
  double forceDistance;
  double sigmaradius;
  int cell_position;
  
   //store current force as previous and reset current force
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!((*it)->empty()))
    { 
      //if cell not empty loop though contents updating positions
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	//store force
	(*it2)->f_x_prev = (*it2)->f_x;
	(*it2)->f_y_prev = (*it2)->f_y; 
	(*it2)->f_z_prev = (*it2)->f_z;
	
	//reset forces to zero
	(*it2)->f_x = 0.0;
	(*it2)->f_y = 0.0;
	(*it2)->f_z = 0.0;
      }
    }
  }
  
  //loop through cells
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    cell_position = std::distance(cells.begin(), it);
    if(!(*it)->empty())
    { 
      
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	std::list<Particle*>::iterator it3 = it2;
	++it3;
	//calculate forces from particles in immediate cell
	for(; it3 !=(*it)->end(); ++it3)
	{
	  
	  
	  xDistance = (*it2)->x - (*it3)->x;
	  yDistance = (*it2)->y - (*it3)->y;
	  zDistance = (*it2)->z - (*it3)->z;
	  //std::cout << (*it2)->x << std::endl;
	  distance = sqrt(xDistance*xDistance + yDistance*yDistance + zDistance*zDistance);
	  std::cout<<distance<<std::endl;
	  
	  if( distance < r_cut)
	  {
	    sigmaradius = (sigma/distance);
	    sigmaradius = sigmaradius*sigmaradius*sigmaradius;
	    sigmaradius = sigmaradius*sigmaradius;
	    
	    forceDistance = 24.0*epsilon*(1.0/(distance*distance))*(sigmaradius)*(1.0-2.0*sigmaradius);
	    
	    xForce = forceDistance*xDistance;
	    yForce = forceDistance*yDistance;
	    zForce = forceDistance*zDistance;
	    
	    (*it2)->f_x += xForce;
	    (*it2)->f_y += yForce;
	    (*it2)->f_z += zForce;
	    
	    (*it3)->f_x -= xForce;
	    (*it3)->f_y -= yForce;
	    (*it3)->f_z -= zForce;
	  }
	}
	
	//std::cout << cells_x*cells_y*cells_z << " "  << adjacentCell(1,0,0,cell_position) <<  " " << cell_position << std::endl;
	//std::cout << cell_position/(cells_x*cells_y) << " " << (cell_position%(cells_x*cells_y)/cells_x) << " " << (cell_position%(cells_x*cells_y))%cells_x <<  std::endl;	

// 	intercellForces(it2, it + adjacentCell(1,0,0,cell_position));	
// 	intercellForces(it2, it + adjacentCell(1,0,-1,cell_position));	
// 	intercellForces(it2, it + adjacentCell(0,0,-1,cell_position));
// 	intercellForces(it2, it + adjacentCell(-1,0,-1,cell_position));	
// 	intercellForces(it2, it + adjacentCell(-1,-1,1,cell_position));
// 	intercellForces(it2, it + adjacentCell(0,-1,1,cell_position));
// 	intercellForces(it2, it + adjacentCell(1,-1,1,cell_position));
// 	intercellForces(it2, it + adjacentCell(-1,-1,0,cell_position));
// 	intercellForces(it2, it + adjacentCell(0,-1,0,cell_position));
// 	intercellForces(it2, it + adjacentCell(1,-1,0,cell_position));
// 	intercellForces(it2, it + adjacentCell(-1,-1,-1,cell_position));
// 	intercellForces(it2, it + adjacentCell(0,-1,-1,cell_position));	
// //std::cout<< "13. " << (cell_position + adjacentCell(1,-1,-1,cell_position))/(cells_x*cells_y) << " " << ((cell_position + adjacentCell(1,-1,-1,cell_position))%(cells_x*cells_y)/cells_x) << " " << ((cell_position + adjacentCell(1,-1,-1,cell_position))%(cells_x*cells_y))%cells_x <<  std::endl;
// 	intercellForces(it2, it + adjacentCell(1,-1,-1,cell_position));
	//std::cout << "helloend" << std::endl;
      }
    }
  }
}

void MolDy::intercellForces(std::list<Particle*>::iterator particle, std::vector< std::list<Particle*>* >::iterator adjacent_cell){
  
  double xDistance;
  double yDistance;
  double zDistance;
  double xForce;
  double yForce;
  double zForce;
  double distance;
  double forceDistance;
  double sigmaradius;
  
  for(std::list<Particle*>::iterator it = (*adjacent_cell)->begin(); it != (*adjacent_cell)->end(); ++it)
  {
    //std::cout << "hello" << std::endl;  
    //std::cout << "hello2 " << (*particle)->x << std::endl; 
    xDistance = (*particle)->x - (*it)->x;
    yDistance = (*particle)->y - (*it)->y;
    zDistance = (*particle)->z - (*it)->z;
    distance = sqrt(xDistance*xDistance + yDistance*yDistance + zDistance*zDistance);
    
    if( distance < r_cut)
    {
      sigmaradius = (sigma/distance);
      sigmaradius = sigmaradius*sigmaradius*sigmaradius;
      sigmaradius = sigmaradius*sigmaradius;
      
      forceDistance = 24*epsilon*1/(distance*distance)*(sigmaradius)*(1-2*sigmaradius);
      
      xForce = forceDistance*xDistance;
      yForce = forceDistance*yDistance;
      zForce = forceDistance*zDistance;
      
      (*particle)->f_x += xForce;
      (*particle)->f_y += yForce;
      (*particle)->f_z += zForce;
      
      (*it)->f_x -= xForce;
      (*it)->f_y -= yForce;
      (*it)->f_z -= zForce;
    }
  }
}

int MolDy::adjacentCell(int x_offset, int y_offset, int z_offset, int vector_index){
  
  int total_offset = 0;
  
  if((vector_index + x_offset)/cells_x  > (vector_index)/cells_x)
  {
    total_offset += (x_offset - cells_x);
  }
  else if((vector_index + x_offset)/cells_x  < (vector_index)/cells_x)
  {
    total_offset += (x_offset + cells_x);
  }
  else
  {
    total_offset += x_offset;
  }
  
  if((vector_index + y_offset*cells_x)/(cells_x*cells_y) > (vector_index)/(cells_x*cells_y))
  {
    total_offset += y_offset*cells_x - cells_x*cells_y;
  }
  else if((vector_index + y_offset*cells_x)/(cells_x*cells_y) < (vector_index)/(cells_x*cells_y))
  {
    total_offset += y_offset*cells_x + cells_x*cells_y;
  }
  else
  {
    total_offset += y_offset*cells_x;
  }
  
  if((vector_index + z_offset*cells_x*cells_y)/(cells_x*cells_y*cells_z) > vector_index/(cells_x*cells_y*cells_z))
  {
    total_offset += z_offset*cells_x*cells_y - cells_x*cells_y*cells_z;
  }
  else if((vector_index + z_offset*cells_x*cells_y)/(cells_x*cells_y*cells_z) < vector_index/(cells_x*cells_y*cells_z))
  {
    total_offset += z_offset*cells_x*cells_y + cells_x*cells_y*cells_z;
  }
  else
  {
    total_offset += z_offset*cells_x*cells_y;
  }
  
  return total_offset;
}

void MolDy::updateVelocities(void){
  //loop through cells
  for(std::vector< std::list<Particle*>* >::iterator it = cells.begin(); it!= cells.end(); ++it)
  {
    if(!((*it)->empty()))
    { 
      //if cell not empty loop though contents updating positions
      for(std::list<Particle*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
      {
	(*it2)->v_x = (*it2)->v_x + (((*it2)->f_x_prev + (*it2)->f_x)*delta_t)/(2*((*it2)->mass));
	(*it2)->v_y = (*it2)->v_y + (((*it2)->f_y_prev + (*it2)->f_y)*delta_t)/(2*((*it2)->mass)); 
	(*it2)->v_z = (*it2)->v_z + (((*it2)->f_z_prev + (*it2)->f_z)*delta_t)/(2*((*it2)->mass));
      }
    }
  }
}

void MolDy::simulate(void){
  for(int t = 0; t < ((t_end-t_start)/delta_t); ++t)
  {
    updatePosition();
    calculateForces();
    updateVelocities();
    if(t%vis_space == 0)
    {
      std::stringstream filename;
      filename << vtkfile << t/vis_space << ".vtk";
      VTKPrint(filename.str());
      std::cout << filename.str() << std::endl;
    }
    //std::cout << time << std::endl;
  }
}
