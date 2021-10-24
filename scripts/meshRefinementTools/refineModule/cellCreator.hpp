
#include "cells.hpp"
#include "scalarFunction.hpp"
#include <vector>

class cellCreator
{


public:

  int maxlevel;

  cellCreator() : maxlevel(3), guaranteedRecursions(3)
  {}

  std::vector<cells> cellData;

  int guaranteedRecursions;

  void setMaxLevels(int levels){maxlevel=levels;}
  void setGuaranteedRecursions(int recurse){guaranteedRecursions=recurse;}

  //
  // 2D versions
  //

  void addCell(int id, double xmid, double ymid, double xlo, double ylo, 
	       double xhi, double yhi, std::vector<double> fValues);

  std::vector<cells> createCells(double xlo, double xhi, double ylo, double yhi,
				 scalarFunction function);

  void recursiveCreateCells(int idcell, int level, double xlo, double xhi, double ylo, 
			    double yhi, scalarFunction function);


  //
  //3D versions
  //

  void addCell(int id, double xmid, double ymid, 
	       double zmid, double xlo, double ylo, double zlo, 
	       double xhi, double yhi, double zhi, std::vector<double> fValues);

  std::vector<cells> createCells(double xlo, double xhi, double ylo, 
				 double yhi, double zlo, double zhi, scalarFunction function);

  void recursiveCreateCells(int idcell, int level, double xlo, double xhi, double ylo, 
			    double yhi, double zlo, double zhi, scalarFunction function);


};

