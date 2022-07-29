#from meshRefine cimport fillCoordinates
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "refine_external.hpp":
     cdef cppclass meshRefine:
          meshRefine()
          void fillCoordinates(vector[double])
          void helloWorld()
          void init(int flag3D)
          void printSurfs()
          void testX(vector[double])
          void printTestX(vector[double])
          void sizeCellNodes(int numCellNodes_)
          void resetNodeCounter()
          void printCoords()
          void freeCoords()
          void setRefinementDistance(double dist_)
          void setRefinementDistanceN(double dist_)
          void setRefinementDistanceP(double dist_)
          void setRefinementDistanceIL(double dist_)
          double getRefinementDistance()
          bool doIRefine()
          bool doIrefineXplane(double x)
          bool doIrefineCentroid()
          bool doIrefineCentroidSided()
          bool doIrefine3D()
          bool doIrefine2D()
          double doIrefineCentroidSidedFinishMetric()
          int getNumSurfs()
          double getMaxSide()
          double getMinSide()
          void setCellMinimum(double min)
          void setCellMinimumN(double min)
          void setCellMinimumP(double min)
          void setCellMinimumIL(double min)
          void setDimension(int dim)
          void setRefinementFactor(double rF)
          int getDimension()
          void readSurfaces(string)
          void setSizeMeasure(string)
          void createSurfaces()
          void journalSurfaces()
          void setDopingFunction(string name, string functionType, string dopingType)
          void setXBounds(string name, string axis, double min, double max)
          void setDopingBounds(string name, double min, double max)
          void setDirection(string name, string direction)
          void listFunctions()
          void printXMLDopingFunctions()

          void setDopingWidth(string name, double xwidth, double ywidth, double zwidth)
          void setDopingWidth(string name, double xwidth, double ywidth)
          void setDopingLocation(string name, double xlocation, double ylocation, double zlocation)
          void setDopingLocation(string name, double xlocation, double ylocation)
          void addRefineToLine(double xmin, double ymin, double xmax, double ymax, double ilThick)
          void addRefineToSurface(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double ilThick)
          void setRefinementBoundingBox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
          void setMaxSurfaceRecursions(int level)
          void setGuaranteedSurfaceRecursions(int level)
          void setWriteJunctions(bool setJunctions)
          void setTetNum(int tn)
          int getTetNum()
          void setAutoRefine()


cdef class PyMeshRefine:
     cdef meshRefine *thisptr     
     def __cinit__(self):
         self.thisptr = new meshRefine()
     def __dealloc__(self):
         del self.thisptr
     def fillCoordinates(self, vector[double] x):
         return self.thisptr.fillCoordinates(x)
     def helloWorld(self):
         return self.thisptr.helloWorld()
     def init(self, int flag3D):
         return self.thisptr.init(flag3D)
     def printSurfs(self):
         return self.thisptr.printSurfs()
     def testX(self, vector[double] x):
         return self.thisptr.testX(x)
     def printTestX(self, vector[double] x):
         return self.thisptr.printTestX(x)
     def sizeCellNodes(self, int numCellNodes_):
         return self.thisptr.sizeCellNodes(numCellNodes_)
     def resetNodeCounter(self):
         return self.thisptr.resetNodeCounter()
     def printCoords(self):
         return self.thisptr.printCoords()
     def freeCoords(self):
         return self.thisptr.freeCoords()
     def setRefinementDistance(self, double dist_):
         return self.thisptr.setRefinementDistance(dist_)
     def setRefinementDistanceN(self, double dist_):
         return self.thisptr.setRefinementDistanceN(dist_)
     def setRefinementDistanceP(self, double dist_):
         return self.thisptr.setRefinementDistanceP(dist_)
     def setRefinementDistanceIL(self, double dist_):
         return self.thisptr.setRefinementDistanceIL(dist_)
     def getRefinementDistance(self):
         return self.thisptr.getRefinementDistance()
     def doIRefine(self):
         return self.thisptr.doIRefine()
     def doIrefineXplane(self, double x):
         return self.thisptr.doIrefineXplane(x)
     def doIrefineCentroid(self):
         return self.thisptr.doIrefineCentroid()
     def doIrefineCentroidSided(self):
         return self.thisptr.doIrefineCentroidSided()
     def doIrefine3D(self):
         return self.thisptr.doIrefine3D()
     def doIrefine2D(self):
         return self.thisptr.doIrefine2D()
     def doIrefineCentroidSidedFinishMetric(self):
         return self.thisptr.doIrefineCentroidSidedFinishMetric()
     def getNumSurfs(self):
         return self.thisptr.getNumSurfs()
     def getMaxSide(self):
         return self.thisptr.getMaxSide()
     def getMinSide(self):
         return self.thisptr.getMinSide()
     def setCellMinimum(self, double min):
         return self.thisptr.setCellMinimum(min)
     def setCellMinimumN(self, double min):
         return self.thisptr.setCellMinimumN(min)
     def setCellMinimumP(self, double min):
         return self.thisptr.setCellMinimumP(min)
     def setCellMinimumIL(self, double min):
         return self.thisptr.setCellMinimumIL(min)
     def setDimension(self, int dim):
         return self.thisptr.setDimension(dim)
     def setRefinementFactor(self, double rF):
         return self.thisptr.setRefinementFactor(rF)
     def getDimension(self):
         return self.thisptr.getDimension()
     def readSurfaces(self, string filename):
         return self.thisptr.readSurfaces(filename)
     def setSizeMeasure(self, string sM):
         return self.thisptr.setSizeMeasure(sM)
     def createSurfaces(self):
         return self.thisptr.createSurfaces()
     def journalSurfaces(self):
         return self.thisptr.journalSurfaces()

     def setDopingFunction(self, string name, string functionType, string dopingType):
         return self.thisptr.setDopingFunction(name, functionType, dopingType)
     def setXBounds(self, string name, string axis, double min, double max):
         return self.thisptr.setXBounds(name, axis, min, max)    
     def setDopingBounds(self, string name, double min, double max):
         return self.thisptr.setDopingBounds(name, min, max)
     def setDirection(self, string name, string direction):
         return self.thisptr.setDirection(name, direction)
     def listFunctions(self):
         return self.thisptr.listFunctions()
     def printXMLDopingFunctions(self):
         return self.thisptr.printXMLDopingFunctions()

     def setDopingWidth(self, string name, double xwidth, double ywidth, double zwidth=0):
         return self.thisptr.setDopingWidth(name, xwidth, ywidth, zwidth)

     def setDopingLocation(self, string name, double xwidth, double ywidth, double zwidth=0):
         return self.thisptr.setDopingLocation(name, xwidth, ywidth, zwidth)

     def addRefineToLine(self, double xmin, double ymin, double xmax, double ymax, double ilThick):
         return self.thisptr.addRefineToLine(xmin, ymin, xmax, ymax, ilThick)

     def addRefineToSurface(self, double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double ilThick):
         return self.thisptr.addRefineToSurface(x0,  y0,  z0,  x1,  y1,  z1,  x2,  y2,  z2,  x3,  y3,  z3,  ilThick)

     def setRefinementBoundingBox(self, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax):
         return self.thisptr.setRefinementBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
     def setMaxSurfaceRecursions(self, int level):
         return self.thisptr.setMaxSurfaceRecursions(level)
     def setGuaranteedSurfaceRecursions(self, int level):
         return self.thisptr.setGuaranteedSurfaceRecursions(level)
     def setWriteJunctions(self, bool setJunctions):
         return self.thisptr.setWriteJunctions(setJunctions)
     def setTetNum(self, int tn):
         return self.thisptr.setTetNum(tn)
     def getTetNum(self):
         return self.thisptr.getTetNum()
     def setAutoRefine(self):
         return self.thisptr.setAutoRefine()
