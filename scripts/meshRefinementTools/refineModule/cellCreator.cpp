
#include "cellCreator.hpp"
#include <vector>
#include "scalarFunction.hpp"
#include <iostream>

std::vector<cells> cellCreator::createCells(double xlo, double xhi, double ylo, 
					    double yhi, double zlo, double zhi, 
					    scalarFunction function)

{

  int idcell=0;
  int level = 0;

  if(guaranteedRecursions > maxlevel)
    guaranteedRecursions = maxlevel;

  std::cout<<" Doing "<<maxlevel<<"  recursions "<<std::endl;

  recursiveCreateCells(idcell, level, xlo, ylo, zlo, 
		       xhi, yhi, zhi, function);


  return cellData;

}



void cellCreator::addCell(int id, double xmid, double ymid, double zmid, 
			  double xlo, double ylo, double zlo, double xhi, 
			  double yhi, double zhi, std::vector<double> fValues)

{


  cells localCell;

  int node_number[8]={0,3,2,1,4,7,6,5};


  localCell.vertices[node_number[0]][0] = xlo;
  localCell.vertices[node_number[0]][1] = ylo;
  localCell.vertices[node_number[0]][2] = zlo;

  localCell.vertices[node_number[1]][0] = xlo;
  localCell.vertices[node_number[1]][1] = yhi;
  localCell.vertices[node_number[1]][2] = zlo;

  localCell.vertices[node_number[2]][0] = xhi;
  localCell.vertices[node_number[2]][1] = yhi;
  localCell.vertices[node_number[2]][2] = zlo;

  localCell.vertices[node_number[3]][0] = xhi;
  localCell.vertices[node_number[3]][1] = ylo;
  localCell.vertices[node_number[3]][2] = zlo;

  localCell.vertices[node_number[4]][0] = xlo;
  localCell.vertices[node_number[4]][1] = ylo;
  localCell.vertices[node_number[4]][2] = zhi;

  localCell.vertices[node_number[5]][0] = xlo;
  localCell.vertices[node_number[5]][1] = yhi;
  localCell.vertices[node_number[5]][2] = zhi;

  localCell.vertices[node_number[6]][0] = xhi;
  localCell.vertices[node_number[6]][1] = yhi;
  localCell.vertices[node_number[6]][2] = zhi;

  localCell.vertices[node_number[7]][0] = xhi;
  localCell.vertices[node_number[7]][1] = ylo;
  localCell.vertices[node_number[7]][2] = zhi;

  for(int i=0 ; i<8 ; ++i)
    {
      localCell.values[node_number[i]] = fValues[i];
    }

  cellData.push_back(localCell);


  return;

}


void cellCreator::recursiveCreateCells(int idcell, int level, double xlo, double ylo, double zlo, 
				       double xhi, double yhi, double zhi, scalarFunction function)
{


  int createflag;
  double xmid,ymid,zmid;
  double mindist,cellsize;
  double size;

  createflag = 0;
  xmid = 0.5 * (xlo+xhi);
  ymid = 0.5 * (ylo+yhi);
  zmid = 0.5 * (zlo+zhi);

  if (level == maxlevel)
    createflag = 1;

  std::vector<double> x;
  x.resize(3);

  std::vector<double> fValues;
  fValues.resize(8);

  x[0] = xlo;
  x[1] = ylo;
  x[2] = zlo;
  fValues[0] = function.evaluateFunction(x);

  x[0] = xlo;
  x[1] = yhi;
  x[2] = zlo;
  fValues[1] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = yhi;
  x[2] = zlo;
  fValues[2] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = ylo;
  x[2] = zlo;
  fValues[3] = function.evaluateFunction(x);

  x[0] = xlo;
  x[1] = ylo;
  x[2] = zhi;
  fValues[4] = function.evaluateFunction(x);

  x[0] = xlo;
  x[1] = yhi;
  x[2] = zhi;
  fValues[5] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = yhi;
  x[2] = zhi;
  fValues[6] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = ylo;
  x[2] = zhi;
  fValues[7] = function.evaluateFunction(x);

  if(level < guaranteedRecursions)
    createflag = 0;

  if (createflag == 1)
    addCell(idcell,xmid,ymid,zmid,xlo,ylo,zlo,xhi,yhi,zhi,fValues);
  else {
    recursiveCreateCells(10*idcell+1,level+1,xlo,ylo,zlo,xmid,ymid,zmid, 
			 function);
    recursiveCreateCells(10*idcell+2,level+1,xlo,ymid,zlo,xmid,yhi,zmid, 
			 function);
    recursiveCreateCells(10*idcell+3,level+1,xmid,ymid,zlo,xhi,yhi,zmid, 
			 function);
    recursiveCreateCells(10*idcell+4,level+1,xmid,ylo,zlo,xhi,ymid,zmid, 
			 function);
    recursiveCreateCells(10*idcell+5,level+1,xlo,ylo,zmid,xmid,ymid,zhi, 
			 function);
    recursiveCreateCells(10*idcell+6,level+1,xlo,ymid,zmid,xmid,yhi,zhi, 
			 function);
    recursiveCreateCells(10*idcell+7,level+1,xmid,ymid,zmid,xhi,yhi,zhi, 
			 function);
    recursiveCreateCells(10*idcell+8,level+1,xmid,ylo,zmid,xhi,ymid,zhi, 
			 function);
    
  }


  return;
}

//-----------------------------------------------------------------------
// create cells 2D
//-----------------------------------------------------------------------

std::vector<cells> cellCreator::createCells(double xlo, double xhi, double ylo, 
					    double yhi, scalarFunction function)

{

  int idcell=0;
  int level = 0;

  if(guaranteedRecursions > maxlevel)
    guaranteedRecursions = maxlevel;


  recursiveCreateCells(idcell, level, xlo, ylo, xhi, 
		       yhi, function);

  std::cout<<" Created "<<cellData.size()<<" 2D cells"<<std::endl;

  return cellData;

}



void cellCreator::addCell(int id, double xmid, double ymid,
			  double xlo, double ylo, double xhi, double yhi,
			  std::vector<double> fValues)

{


  cells localCell;

  int node_number[4]={0,1,2,3};

  localCell.vertices[node_number[0]][0] = xlo;
  localCell.vertices[node_number[0]][1] = ylo;

  localCell.vertices[node_number[1]][0] = xlo;
  localCell.vertices[node_number[1]][1] = yhi;

  localCell.vertices[node_number[2]][0] = xhi;
  localCell.vertices[node_number[2]][1] = yhi;

  localCell.vertices[node_number[3]][0] = xhi;
  localCell.vertices[node_number[3]][1] = ylo;

  for(int i=0 ; i<4 ; ++i)
    localCell.values[node_number[i]] = fValues[i];

  cellData.push_back(localCell);


  //std::cout<<"create vertex "<<xlo<<"   "<<ylo<<std::endl;
  //std::cout<<"create vertex "<<xlo<<"   "<<yhi<<std::endl;
  //std::cout<<"create vertex "<<xhi<<"   "<<yhi<<std::endl;
  //std::cout<<"create vertex "<<xhi<<"   "<<ylo<<std::endl;


  return;

}


void cellCreator::recursiveCreateCells(int idcell, int level, double xlo, double ylo,
				       double xhi, double yhi, scalarFunction function)
{


  int createflag;
  double xmid,ymid;
  double mindist,cellsize;
  double size;

  createflag = 0;
  xmid = 0.5 * (xlo+xhi);
  ymid = 0.5 * (ylo+yhi);

  std::vector<double> x;
  x.resize(2);

  std::vector<double> fValues;
  fValues.resize(4);

  x[0] = xlo;
  x[1] = ylo;
  fValues[0] = function.evaluateFunction(x);

  x[0] = xlo;
  x[1] = yhi;
  fValues[1] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = yhi;
  fValues[2] = function.evaluateFunction(x);

  x[0] = xhi;
  x[1] = ylo;
  fValues[3] = function.evaluateFunction(x);

  //std::cout<<" Function values "<<fValues[0]<<"    "<<fValues[1]<<"    "<<fValues[2]<<"    "<<fValues[3]<<std::endl;

  //If any one of the fValues differs in sign from the others, createflag = 0
  //If all fValues have the same sign, createflag = 1
  bool allPos=true, allNeg=true;
  for(size_t fV=0 ; fV<fValues.size() ; ++fV)
    {
      if(fValues[fV] < 0)
	{
	  allPos = false;
	  break;
	}
    }
  for(size_t fV=0 ; fV<fValues.size() ; ++fV)
    {
      if(fValues[fV] > 0)
	{
	  allNeg = false;
	  break;
	}
    }

  if(allPos || allNeg && level)
    createflag = 1;
  if (level == maxlevel)
    createflag = 1;


  if(level < guaranteedRecursions)
    createflag = 0;


  if (createflag == 1)
    {
      addCell(idcell,xmid,ymid,xlo,ylo,xhi,yhi,fValues);
    }
  else {
    recursiveCreateCells(10*idcell+1,level+1,xlo,ylo,xmid,ymid, 
			 function);
    recursiveCreateCells(10*idcell+2,level+1,xlo,ymid,xmid,yhi, 
			 function);
    recursiveCreateCells(10*idcell+3,level+1,xmid,ymid,xhi,yhi, 
			 function);
    recursiveCreateCells(10*idcell+4,level+1,xmid,ylo,xhi,ymid, 
			 function);
  }


  return;
}
