

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "table.h"
#include "tet_table.h"
#include "surfaceFinder.hpp"
#include "vector_macros.h"

/*******************************************************************************/
/* Private variables */

double f_center [6]; /* for ambiguous cases, center vals for each face in cell */
int ambiguous_case; /* remember which ambiguous case we are handling for when
                       we generate extra vertices */ 

/* node renumbering to convert from Exodus node order to table node order */
int renumber_node [8] = {0, 1, 3, 2, 4, 5, 7, 6};  /* circular face ordering */

/*******************************************************************************/
/* Module:
      colinear_test - tests if three points are colinear to determine if
         a zero area triangle has been created
   Inputs:
      p1 - point one as an (x,y,z) vector
      p2 - point two as an (x,y,z) vector
      p3 - point three as an (x,y,z) vector
      epsilon - delta term for determining proximity to zero
   Returns
      true/false - true if they are colinear, false if they are not
*/

bool surfaceFinder::colinear_test (Vec p1, Vec p2, Vec p3, double epsilon)
{
   Vec v1, v2;
   Vec v3;
   double length;

   /* create vectors from the points */
   VecSub (p1, p2, v1);
   VecSub (p3, p2, v2);

   /* get the plane normal by taking the cross product of the vectors */
   VecCross (v1, v2, v3);

   /* get the length of the normal, skip sqrt since only looking for 0
      length */
   length = v3 [0] * v3 [0] + v3 [1] * v3 [1] + v3 [2] * v3 [2];
/*   if (length < epsilon)      * if v1 and v2 are parallel *   */
   if (length < epsilon*epsilon)      /* if close to zero area: RCS mod. Apr06  */
      return (true);
   else return (false);
} /* colinear_test */

/*******************************************************************************/
/* Module:
      select_subcase - uses the center value of each ambiguous face to
         determine which subcase table entry to use
   Inputs:
      index          - the index into the marching cubes table
      values         - array of node values 
      inverse_flag   - set to -1 if the index was >= 128 and was complemented;
                       set to 1 if the index is < 128 and has not been changed
      threshold      - surface value for the isosurface being generated
   Outputs:
      index          - modified to point to the correct subcase in the subcase
                       table
*/

void surfaceFinder::select_subcase (int *index, double *values, int inverse_flag, 
   double threshold)
{
   int num_faces, face;
   int i;
   char *p;
   int offset;
   int v0, v1, v2, v3;

   /* skip over case number and default subtable location */
   p = (char *) &cell_info [*index - 1].triangle [0][2];
   num_faces = (int) (*p++);
   offset = 0;

   for (i = 0; i < num_faces; i++)
   {
      face = (int) (*p++);
      v0 = renumber_node [(int) face_vertices [face][0]];
      v1 = renumber_node [(int) face_vertices [face][1]];
      v2 = renumber_node [(int) face_vertices [face][2]];
      v3 = renumber_node [(int) face_vertices [face][3]];

      /* the center is found by averaging the four nodes of the face */
      f_center [face] = (values [v0] + values [v1] + values [v2] + values [v3]) 
         / 4.0;

      if (inverse_flag == 1)
      {
         if (f_center [face] > threshold)
            offset += 1 << i;
      }
      else
      {
         if (f_center [face] <= threshold)
            offset += 1 << i;
      }
   }
   *index = (int) (*p) + offset;
} /* select_subcase */

/*******************************************************************************/
/* Module:
      calculate_interior_vertex - This is only called for case 10 which has 
         only 2 triangles in its tesselation.
   Inputs:
      index     - the index into the marching cubes table
      coords    - array of node coords 
      threshold - surface value for the isosurface being generated
   Outputs:
      eptr      - pointer in which to return the vertex location vector 
*/

void surfaceFinder::calculate_interior_vertex (int index, Vec *coords, double threshold, 
   edgeStruct *eptr)
{
   int num_faces, face;
   int i;
   int v0, v1, v2, v3;
   double distance;
   double center_val [2];
   Vec center_pt [2];
   char *p;

   /* skip over case number and default subtable location */
   p = (char *) &cell_info [index - 1].triangle [0][2];
   num_faces = (int) (*p++);
   for (i = 0; i < num_faces; i++)
   {
      face = (int) (*p++);
      v0 = renumber_node [(int) face_vertices [face][0]];
      v1 = renumber_node [(int) face_vertices [face][1]];
      v2 = renumber_node [(int) face_vertices [face][2]];
      v3 = renumber_node [(int) face_vertices [face][3]];

      center_pt [i][0] = (coords [v0][0] + coords [v1][0] 
         + coords [v2][0] + coords [v3][0]) / 4.0;
      center_pt [i][1] = (coords [v0][1] + coords [v1][1] 
         + coords [v2][1] + coords [v3][1]) / 4.0;
      center_pt [i][2] = (coords [v0][2] + coords [v1][2] 
         + coords [v2][2] + coords [v3][2]) / 4.0;

      center_val [i] = f_center [face];
   } /* end for i */

   if (center_val [0] != center_val [1])
      distance = (double) ((threshold - center_val [0]) / (center_val [1] - 
         center_val [0]));
   else distance = .5;

   eptr->e1 = -1;
   eptr->e2 = -1;
   eptr->x [0] = center_pt [0][0] + distance * (center_pt [1][0] - 
      center_pt [0][0]);
   eptr->x [1] = center_pt [0][1] + distance * (center_pt [1][1] - 
      center_pt [0][1]);
   eptr->x [2] = center_pt [0][2] + distance * (center_pt [1][2] - 
      center_pt [0][2]);
} /* calculate_interior_vertex */

/*******************************************************************************/
/* Module:
      calculate_center_vertex - finds the coordinates of the center vertex
         for all the cases that are interpolations of two diagonally
         located nodes within the cell
   Inputs:
      coords    - array of node coordinates for the cell
      values    - array of node values for the cell
      index     - the index into the marching cubes table
      threshold - surface value for the isosurface being generated
   Outputs:
      eptr      - pointer in which to return the vertex location vector
*/

void surfaceFinder::calculate_center_vertex (int index, Vec *coords, double *values, 
   double threshold, edgeStruct *eptr)
{
   double distance;
   int v1, v2;

   /* For each case, set v1 and v2 to the node indices for the diagonally
      opposite sides of the cell. */
   switch (index)
   {
      case 22:
      case 104:
         v1 = renumber_node [2];
         v2 = renumber_node [5];
         break;
      case 41:
      case 107:
         v1 = renumber_node [3];
         v2 = renumber_node [4];
         break;
      case 73:
      case 105:
      case 109:
         v1 = renumber_node [0];
         v2 = renumber_node [7];
         break;
      case 97:
      case 121:
         v1 = renumber_node [1];
         v2 = renumber_node [6];
         break;
   }

   distance = (double) ((threshold - values [v1]) / (values [v2] - values [v1]));

   /* Bilinearly interpolate between the values at the nodes on diagonally
      opposite sides of the cell */
   eptr->e1 = -1;
   eptr->e2 = -1;
   eptr->x [0] = coords [v1][0] + distance * (coords [v2][0] - coords [v1][0]);
   eptr->x [1] = coords [v1][1] + distance * (coords [v2][1] - coords [v1][1]);
   eptr->x [2] = coords [v1][2] + distance * (coords [v2][2] - coords [v1][2]);
} /* calculate_center_vertex */

/*******************************************************************************/
/* Module:
      get_subcase_entry - sets the table_ptr to the correct subtable entry
         for ambiguous cases
   Inputs:
      index          - the index into the marching cubes table
      values         - array of node values for the cell
      inverse_flag   - set to -1 if the index was >= 128 and was complemented;
                       set to 1 if the index is < 128 and has not been changed
      threshold      - surface value for the isosurface being generated
   Outputs:
      index          - case table index 
      triangle_cnt   - the number of triangles needed to tesselate a cell with
                       the configuration of nodes >= the surface given by
                       the value of index
   Returns:
      table_ptr      - points to the correct subtable entry for this case
*/

char *surfaceFinder::get_subcase_entry (int index, double *values, int inverse_flag,
   double threshold, int *triangle_cnt)
{
   char *table_ptr;
   int subcase_index; /* modified to pt to correct subcase in subcase table */ 

   ambiguous_case = (int) (cell_info [index - 1].triangle [0][0]);
   subcase_index = index;
   select_subcase (&subcase_index, values, inverse_flag, threshold);
   switch (ambiguous_case)
   {
      case 3:
         table_ptr = &subtable3 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable3 [subcase_index].n_triangles;
         break;
      case 6:
         table_ptr = &subtable6 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable6 [subcase_index].n_triangles;
         break;
      case 7:
         table_ptr = &subtable7 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable7 [subcase_index].n_triangles;
         break;
      case 10:
         table_ptr = &subtable10 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable10 [subcase_index].n_triangles;
         break;
      case 12:
         table_ptr = &subtable12 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable12 [subcase_index].n_triangles;
         break;
      case 13:
         table_ptr = &subtable13 [subcase_index].triangle [0][0];
         *triangle_cnt = (int) subtable13 [subcase_index].n_triangles;
         break;
   }
   return (table_ptr);
} /* get_subcase_entry */

/*******************************************************************************/
/* Module:
      generate_center_vertex - This module creates a vertex that is not along 
         edges zero through eleven.  The vertex is either in the center of the 
         cell or the center of a face, depending on the ambiguous case of the 
         cell.
   Inputs:
      index     - the index into the marching cubes table
      coords    - array of node coordinates for the cell
      values    - array of node values for the cell
      threshold - surface value for the isosurface being generated
   Outputs:
      eptr      - pointer to the edge list element in which to insert the 
                  vertex location being calculated
*/

void surfaceFinder::generate_center_vertex (int index, Vec *coords, double *values, 
   double threshold, edgeStruct *eptr)
{
   if (ambiguous_case == 7)
      calculate_center_vertex (index, coords, values, threshold, eptr);
   else if (ambiguous_case == 10)
      calculate_interior_vertex (index, coords, threshold, eptr);
   else if (ambiguous_case == 13)
      calculate_center_vertex (index, coords, values, threshold, eptr);
} /* generate_center_vertex */

/*******************************************************************************/
/* Module:
      tesselate_hex_cell
   Inputs:
      coords    - array of node coordinates for the cell
      values    - array of node values for the cell
      threshold - surface value for the isosurface being generated
      epsilon   - error difference below which 2 vertex coords are the same
   Outputs:
      total_vertices  - the number of vertices created for the cell
      edge_list       - a struct with the 2 node indices defining the
                        edge used to interpolate the vertex & the vertex
                        coordinates
      total_triangles - the number of triangles created for the cell
      t               - a list of the triangles created
*/
void surfaceFinder::tesselate_hex_cell (Vec *coords, double *values,  double threshold, 
   double epsilon, int *total_vertices, edgeStruct **edge_list, 
   int *total_triangles, int **t)
{
   int e,               /* current edge being interpolated from table lookup */
       v0, v1,          /* endpts of current edge */
       i, j, m, a;      /* loop indices */
   int index;           /* index into table created by categorizing nodes */
   double distance;      /* distance to interpolated pt in pixels */
   int *tptr;           /* running ptr for outputting tris in current cell */
   int triangle_cnt;    /* #tris for configuration defined by index in table */
   char *table_ptr;     /* running ptr to edges in case table */
   double *t1, *t2, *t3;  /* vertex ptrs for checking for 0 area triangles */
   edgeStruct *eptr;    /* running ptr for filling edge_list with new verts */
   edgeStruct *adjptr;  /* ptr for finding duplicate verts in edge_list */
   int *ciptr,          /* triangle list ptr used when filling in gaps */
       inc;             /* increment for tptr (negative if complement) */
   int degenerate;      /* count of zero area triangles found & removed */

   /* look up table for previously created vertices */
   int LUT [13] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

   /* Create Marching Cubes table index by comparing each node's value
      to the surface value given in threshold. */
   index = 0;

   /* Lorensen & Cline test for >= threshold, but this creates triangles of
      zero area for cases 1, 3, 4, 6, 7, 12, and 13, so I switched to testing
      for being above the surface, but not on it. */
   if (values [renumber_node[0]] > threshold)
      index += 1;
   if (values [renumber_node[1]] > threshold)
      index += 2;
   if (values [renumber_node[2]] > threshold)
      index += 4;
   if (values [renumber_node[3]] > threshold)
      index += 8;
   if (values [renumber_node[4]] > threshold)
      index += 16;
   if (values [renumber_node[5]] > threshold)
      index += 32;
   if (values [renumber_node[6]] > threshold)
      index += 64;
   if (values [renumber_node[7]] > threshold)
      index += 128;

   if ((index > 0) && (index < 255)) /* if surface passes through this cell */
   {
      /* The Marching Cubes table only contains half the entries since 
         the tesselation for the complementary index would be the same, only
         the order of enumerating the vertices would be reversed (since
         the right-hand rule would point in the opposite direction).  So
         if the index is for the second half of the table, we complement
         it and then build the triangles in the reverse order.  We set
         the variable 'inc' to +/- 1 and use it to update the triangle_list
         pointer.  It is also a flag for the complementary case when it is
         set to -1. 
      */
      if (index >= 128)
      {
         inc = -1;
         index = (~index) & 0377;
      }
      else inc = 1;

      degenerate = 0; /* initialize counter for # of degenerate tris found */

      /* Get triangle count for this configuration from table */
      triangle_cnt = (int) cell_info [index - 1].n_triangles;

      if (triangle_cnt == 0) /* if = 0, this is an ambiguous case */
      {
         /* Get a pointer to the proper entry in the proper subcase table */
         table_ptr = get_subcase_entry (index, values, inc, threshold,
            &triangle_cnt);
      }
      else /* this is not an ambiguous case, so use normal table entry */ 
         table_ptr = &cell_info [index - 1].triangle [0][0];

      /* Update the triangle count by the triangles in this cell */
      *total_triangles += triangle_cnt;

      /* Initialize running pointers to the edge list and the triangle list. */
      eptr = *edge_list;
      tptr = *t;

      /* if complemented index, set the triangle list pointer to the
         end of the space it will use for this triangle count. 
         The triangle vertices will be added in the reverse order to flip
         the direction of the normal. */
      if (inc == -1)    
         tptr = *t + triangle_cnt * 3 - 1;

      /* For each triangle (the i loop), generate each vertex (the j loop)
         by getting the edge number from the table.  If the vertex has
         already been created, its index within the edge_list will be
         given by the LUT value for that edge.  Values of -1 in the LUT
         represent edges that have not been interpolated yet.
      */
      for (i = 0; i < triangle_cnt; i++)
      {
         for (j = 0; j < 3; j++)
         {
            e = (int) *table_ptr++; /* get next edge from table */
            if (LUT [e] == -1)    /* if this edge's vertex has not been done */
            {
               LUT [e] = *total_vertices; /* set LUT to new vertex index */
               *total_vertices += 1;      /* update vertex count/index */

               /* e = 12 is used as a flag to create a vertex off the edges */
               if (e == 12) /* if vertex is not interpolated on an edge */ 
               {
                  /* create a vertex either at the cell center or in the
                     center of an ambiguous face */
                  generate_center_vertex (index, coords, values, threshold, 
                     eptr);
               }
               else /* the vertex is interpolated along an edge */
               {
                  /* convert from the table numbering for the edge endpoints
                     to exodus numbering for the nodes */
                  v0 = renumber_node [(int) edges [e][0]];
                  v1 = renumber_node [(int) edges [e][1]];

                  /* avoid a divide by zero if the values at the endpoints
                     are the same */
                  if (values [v1] - values [v0] != 0.0) 
                  {
                     /* calculate the interpolative distance for the surface
                        based on the values at the endpoints of the edge */
                     distance = (double) (threshold - values [v0])
                        / (double) (values [v1] - values [v0]);
                  } /* end if (values [v1] - values [v0] != 0.0) */
                  else /* if =, then there is some precision problem */ 
                     distance = 0.0;

                  eptr->e1 = v0; /* return the node index for endpoint 1 */
                  eptr->e2 = v1; /* return the node index for endpoint 2 */

                  /* Calculate the vertex position from the distance and
                     the endpoint node coordinates */
                  eptr->x [0] = coords [v0][0] + distance * 
                     (coords [v1][0] - coords [v0][0]);
                  eptr->x [1] = coords [v0][1] + distance * 
                     (coords [v1][1] - coords [v0][1]);
                  eptr->x [2] = coords [v0][2] + distance * 
                     (coords [v1][2] - coords [v0][2]);

                  /* For vertices generated very near the endpoints, the
                     LUT will not prevent 3 different vertices from being
                     created on each of the edges using that endpoint.
                     To find these duplicate vertices, we do a comparison
                     of vertex coordinates for any vertices already
                     created on edges that share an endpoint with the
                     current edge.  We use an epsilon error term to allow 
                     us to find those vertices that are essentially the
                     same, but are not exactly equal.  The static array
                     adj_edges contains each of the 4 other edges that
                     share endpoints with the edge e.

                            /                       /
                           1                       2
                          /                       /
                         /                       /
                        x-----------e-----------x
                        |                       |
                        |                       |
                        3                       4
                        |                       |
                  */
                  for (m = 0; m < 4; m++)
                  {
                     a = adj_edges [e][m];
                     if (LUT [a] != -1) /* if an earlier vertex was created */
                     {
                        adjptr = *edge_list + LUT [a]; /* pt to it */
                        /* compare its coords with the new vertex just 
                           created, plus or minus epsilon */
                        if ((adjptr->x [0] >= ((eptr->x [0]) - epsilon)) &&
                           (adjptr->x [0] <= ((eptr->x [0]) + epsilon)) &&
                           (adjptr->x [1] >= (eptr->x [1] - epsilon)) &&
                           (adjptr->x [1] <= (eptr->x [1] + epsilon)) &&
                           (adjptr->x [2] >= (eptr->x [2] - epsilon)) &&
                           (adjptr->x [2] <= (eptr->x [2] + epsilon)))
                        {
                           /* if they are essentially equal, remove the
                              current vertex from the output list and put
                              the earlier vertex's index into the LUT */
                           *total_vertices -= 1;
                           LUT [e] = LUT [a];
                           eptr--; /* since this ptr is advanced below, this 
                                      cancels out the update later */
                           break;
                        } /* end if */
                     } /* end if (LUT [a] != -1) */
                  } /* end for m */
               }  /* end else e != 12 */
               eptr++; /* update edge_list pointer */
            } /* end if (LUT [e] == -1) */
            *tptr = LUT [e]; /* set triangle to index new vertex */
            tptr += inc;     /* update triangle list pointer by +/- 1 */
         } /* end j loop - looping once per triangle vertex */

         /* Since the triangle just created could have zero area, check
            for colinear vertices to decide if this triangle could be
            eliminated - no point in taking up space with it otherwise.
            Start by setting up pointers to the vertices. 
         */ 
         t1 = &((*edge_list + *(tptr - inc))->x [0]);
         t2 = &((*edge_list + *(tptr - 2*inc))->x [0]);
         t3 = &((*edge_list + *(tptr - 3*inc))->x [0]);
         if (colinear_test (t1, t2, t3, epsilon))
         {
            /* Since the vertices are colinear, remove the triangle */
            tptr -= 3 * inc;
            *total_triangles -= 1;
            degenerate++; /* needed since could be building list from back */
         }
      } /* end i loop - looping once per triangle */

      /* If we are building the triangle list in reverse order, removing
         triangles (represented by the degenerate flag being > 1) will
         create gaps.  In this case the list must be compressed to fill
         in the gaps.
      */
      if ((degenerate) && (inc == -1))/* if gaps and reverse order */
      {
         tptr -= inc; /* since tptr pts to next location to fill, backup 1 */
         ciptr = *t; /* set to starting triangle list location */

         /* for all the real triangles (cnt - deleted tris), move them */
         for (i = 0; i < (triangle_cnt - degenerate) * 3; i++)
            *ciptr++ = *tptr++;
      }
   } /* end if ((index > 0) && (index < 255)) */
} /* tesselate_hex_cell */

/*******************************************************************************/
/* Module:
      tesselate_tet_cell
   Inputs:
      coords    - array of node coordinates for the cell
      values    - array of node values for the cell
      threshold - surface value for the isosurface being generated
      epsilon   - error difference below which 2 vertex coords are the same
   Outputs:
      total_vertices  - the number of vertices created for the cell
      edge_list       - a struct with the 2 node indices defining the
                        edge used to interpolate the vertex & the vertex
                        coordinates
      total_triangles - the number of triangles created for the cell
      t               - a list of the triangles created
*/
void surfaceFinder::tesselate_tet_cell (Vec *coords, double *values,  double threshold, 
   double epsilon, int *total_vertices, edgeStruct **edge_list, 
   int *total_triangles, int **t)
{
   int e,               /* current edge being interpolated from table lookup */ 
       v0, v1,          /* endpts of current edge */ 
       i, j, m, a;      /* loop indices */
   int index;           /* index into table created by categorizing nodes */
   double distance;      /* distance to interpolated pt in pixels */
   int *tptr;           /* running ptr for outputting tris in current cell */
   int triangle_cnt;    /* #tris for configuration defined by index in table */
   char *table_ptr;     /* running ptr to edges in case table */
   double *t1, *t2, *t3;  /* vertex ptrs for checking for 0 area triangles */
   edgeStruct *eptr;    /* running ptr for filling edge_list with new verts */
   edgeStruct *adjptr;  /* ptr for finding duplicate verts in edge_list */

   /* look up table for previously created vertices */
   int LUT [NUM_TET_EDGES];    

   /* initialize vertex lookup table to undefined for all edges */
   for (i = 0; i < NUM_TET_EDGES; i++)
      LUT [i] = -1;

   /* Create tet table index by comparing each node's value to the surface 
      value given in 'threshold'. */
   index = 0;

   /* I test for being above the surface, but not on it. */
   if (values [renumber_node [0]] > threshold)
      index += 1;
   if (values [renumber_node [1]] > threshold)
      index += 2;
   if (values [renumber_node [2]] > threshold)
      index += 4;
   if (values [renumber_node [3]] > threshold)
      index += 8;

   if ((index > 0) && (index < 15)) /* if surface passes through this cell */
   {
      /* Get triangle count for this configuration from table */
      triangle_cnt = (int) tet_cell_info [index - 1].n_triangles;

      /* There are no ambiguous cases with tets, so use normal table entry */
      table_ptr = &tet_cell_info [index - 1].triangle [0][0];

      /* Update the triangle count by the triangles in this cell */
      *total_triangles += triangle_cnt;

      /* Initialize running pointers to the edge list and the triangle list. */
      eptr = *edge_list;
      tptr = *t;

      /* For each triangle (the i loop), generate each vertex (the j loop)
         by getting the edge number from the table.  If the vertex has
         already been created, its index within the edge_list will be
         given by the LUT value for that edge.  Values of -1 in the LUT
         represent edges that have not been interpolated yet.
      */
      for (i = 0; i < triangle_cnt; i++)
      {
         for (j = 0; j < 3; j++)
         {
            e = (int) *table_ptr++; /* get next edge from table */
            if (LUT [e] == -1)    /* if this edge's vertex has not been done */
            {
               LUT [e] = *total_vertices; /* set LUT to new vertex index */
               *total_vertices += 1;      /* update vertex count/index */

               /* convert from the table numbering for the edge endpoints
                  to exodus numbering for the nodes */
               v0 = renumber_node [(int) tet_edges [e][0]];
               v1 = renumber_node [(int) tet_edges [e][1]];

               /* avoid a divide by zero if the values at the endpoints
                  are the same */
               if ((values [v1] - values [v0]) != 0.0)
               {
                  /* calculate the interpolative distance for the surface
                     based on the values at the endpoints of the edge */
                  distance = (double) (threshold - values [v0])
                     / (double) (values [v1] - values [v0]);
               } /* end if (values [v1] - values [v0] != 0.0) */
               else /* if =, then there is some precision problem */ 
                  distance = 1.0;

               eptr->e1 = v0; /* return the node index for endpoint 1 */
               eptr->e2 = v1; /* return the node index for endpoint 2 */

               /* Calculate the vertex position from the distance and
                  the endpoint node coordinates */
               eptr->x [0] = coords [v0][0] + distance * 
                  (coords [v1][0] - coords [v0][0]);
               eptr->x [1] = coords [v0][1] + distance * 
                  (coords [v1][1] - coords [v0][1]);
               eptr->x [2] = coords [v0][2] + distance * 
                  (coords [v1][2] - coords [v0][2]);

               /* For vertices generated very near the endpoints, the
                  LUT will not prevent 3 different vertices from being
                  created on each of the edges using that endpoint.
                  To find these duplicate vertices, we do a comparison
                  of vertex coordinates for any vertices already
                  created on edges that share an endpoint with the
                  current edge.  We use an epsilon error term to allow
                  us to find those vertices that are essentially the
                  same, but are not exactly equal.  The static array
                  adj_edges contains each of the 4 other edges that
                  share endpoints with the edge e.

                            /               \
                           1                 2
                          /                   \ 
                         /                     \ 
                        x-----------e-----------x
                        \                      /
                         \                    /
                          3                  4
                           \                /
               */
               for (m = 0; m < 4; m++)
               {
                  a = tet_adj_edges [e][m];
                  if (LUT [a] != -1) /* if an earlier vertex was created */
                  {
                     adjptr = *edge_list + LUT [a];
                     /* compare its coords with the new vertex just
                        created, plus or minus epsilon */
                     if ((adjptr->x [0] >= ((eptr->x [0]) - epsilon)) &&
                        (adjptr->x [0] <= ((eptr->x [0]) + epsilon)) &&
                        (adjptr->x [1] >= (eptr->x [1] - epsilon)) &&
                        (adjptr->x [1] <= (eptr->x [1] + epsilon)) &&
                        (adjptr->x [2] >= (eptr->x [2] - epsilon)) &&
                        (adjptr->x [2] <= (eptr->x [2] + epsilon)))
                     {
                        /* if they are essentially equal, remove the
                           current vertex from the output list and put
                           the earlier vertex's index into the LUT */
                        *total_vertices -= 1;
                        LUT [e] = LUT [a];
                        eptr--; /* since this ptr is advanced below, this
                                   cancels out the update later */
                        break;
                     } /* end if */
                  } /* end if (LUT [a] != -1) */
               } /* end for m */
               eptr++; /* update edge_list pointer */
            } /* end if (LUT [e] == -1) */
            *tptr = LUT [e]; /* set triangle to index new vertex */
            tptr += 1;       /* update triangle list pointer */
         } /* end j loop - looping once per triangle vertex */

         /* Since the triangle just created could have zero area, check
            for colinear vertices to decide if this triangle could be
            eliminated - no point in taking up space with it otherwise.
            Start by setting up pointers to the vertices.
         */
         t1 = &((*edge_list + *(tptr - 1))->x [0]);
         t2 = &((*edge_list + *(tptr - 2))->x [0]);
         t3 = &((*edge_list + *(tptr - 3))->x [0]);
         if (colinear_test (t1, t2, t3, epsilon))
         {
            /* Since the vertices are colinear, remove the triangle */
            tptr -= 3;
            *total_triangles -= 1;
         }
      } /* end i loop - looping once per triangle */
   } /* end if ((index > 0) && (index < 15)) */
} /* tesselate_tet_cell */

/*******************************************************************************/
/* Module:
      draw_cell - tesselates hex or tet cells with triangles 
   Inputs:
      type      - cell type (HEX or TETRA)
      vertices  - cell vertices, coordinates consisting of an (x,y,z) triplet
      values    - cell values corresponding to each of the vertices passed
      threshold - surface value
      epsilon   - error difference below which 2 vertex coords are the same
   Outputs:
      total_verts   - the number of vertices generated
      edge_list     - a struct with the 2 node indices defining the
                      edge used to interpolate the vertex & the vertex
                      coordinates
      total_tris    - the number of triangles generated
      triangle_list - a list of the triangles generated
*/
void surfaceFinder::draw_cell (int type, Vec *vertices , double *values, double threshold, 
   double epsilon, int *total_verts, edgeStruct **edge_list, int *total_tris, 
   int **triangle_list)
{
   /* This initialization of the vertex and triangle counts can be
      eliminated if you want to pass in pointers to running counts
      that just get incremented with each cell. 
   */

   *total_verts = 0;
   *total_tris = 0;

   /* These list initializations for the edge list and the triangle
      list can be eliminated if you want to have the vertices and triangles
      added onto global lists.  They are just a ballpark guess as to size.
      I'm assuming that you'd never create more than 10 triangles per
      cell.  This is multiplied by 3 for each of the vertices in a triangle. 
   */

   if (*edge_list == NULL)
      *edge_list = (edgeStruct *) malloc (3 * 10 * sizeof (edgeStruct));
   if (*triangle_list == NULL)
      *triangle_list = (int *) malloc (3 * 10 * sizeof (int));

   if (type == HEX)
   {
      tesselate_hex_cell (vertices, values, threshold, epsilon, total_verts, 
         edge_list,  total_tris, triangle_list);
   }
   else if (type == TETRA)
   {
      tesselate_tet_cell (vertices, values, threshold, epsilon, total_verts, 
         edge_list,  total_tris, triangle_list);
   }
} /* draw_cell */

