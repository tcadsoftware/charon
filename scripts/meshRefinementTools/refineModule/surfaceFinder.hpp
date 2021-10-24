
#ifndef SURFACEFINDER_H
#define SURFACEFINDER_H


#ifndef FILE_DEFS
#define FILE_DEFS 
#define FAILURE   0
#define SUCCESS   1
#define UNDEFINED -1
#endif /* !FILE_DEFS */

/* cell types */
#define TETRA                 4
#define HEX                   8

#ifndef Vec
typedef double Vec [3];
#endif /* Vec */


class edgeStruct 
{
public:
  /* two indices of the parent edge */
  int e1;
  int e2;
  
  /* coordinates */
  double x [3];
};


class surfaceFinder
{


  bool colinear_test (Vec p1, Vec p2, Vec p3, double epsilon);
  void select_subcase (int *index, double *values, int inverse_flag, 
		       double threshold);
  void calculate_interior_vertex (int index, Vec *coords, double threshold, 
				  edgeStruct *eptr);
  void calculate_center_vertex (int index, Vec *coords, double *values, 
				double threshold, edgeStruct *eptr);
  char *get_subcase_entry (int index, double *values, int inverse_flag,
			   double threshold, int *triangle_cnt);
  void generate_center_vertex (int index, Vec *coords, double *values, 
			       double threshold, edgeStruct *eptr);
  void tesselate_hex_cell (Vec *coords, double *values,  double threshold, 
			   double epsilon, int *total_vertices, edgeStruct **edge_list, 
			   int *total_triangles, int **t);
  void tesselate_tet_cell (Vec *coords, double *values,  double threshold, 
			   double epsilon, int *total_vertices, edgeStruct **edge_list, 
			   int *total_triangles, int **t);


public:

  void draw_cell (int type, Vec *vertices , double *values, double threshold, 
		  double epsilon, int *total_verts, edgeStruct **edge_list, int *total_tris, 
		  int **triangle_list);
  
};



#endif


