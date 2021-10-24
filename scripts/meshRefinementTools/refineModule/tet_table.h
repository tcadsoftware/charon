
#ifndef TET_TABLE_H
#define TET_TABLE_H

/* This table is based on the exodusII vertex numbering scheme.  I take

0: (0,0,0)
1: (1,0,0)
2: (0,1,0)
3: (0,0,1)
as the vertices of a unit tetrahedron (ist binary digit for x, 2nd for y, 3rd
for z)

 How to relate the zero based edge number to the zero based vertex number

    Edge from to  zyx  index

            0  1  001    0
            0  2  010    1
            1  2  011    2
            0  3  100    3
            1  3  101    4
            2  3  110    5

 */
#define NUM_TRIS 				2
#define NUM_SIDES 			3
#define NUM_TET_EDGES 		6

typedef struct {
  char n_triangles;
  char triangle[NUM_TRIS][NUM_SIDES];
} triangulated_tets;

/* for each of the 6 edges, the vertex numbers of the endpts */
char tet_edges[NUM_TET_EDGES][2] = {
  {0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}
};

/* for each of the 4 vertices, the edges that use that vertex as an endpt */
char tet_vertex_edges [4][3] = {
   {0, 1, 3}, {0, 2, 4}, {1, 2, 5}, {3, 4, 5}
};

/* for each of the 6 edges, what other edges share its two endpts */
char tet_adj_edges [NUM_TET_EDGES][4] = {
   {1, 3, 2, 4}, {0, 3, 2, 5}, {0, 4, 1, 5}, {0, 1, 4, 5},
   {0, 2, 3, 5}, {1, 2, 3, 4}
};

char tet_face_vertices [4][3] = { /* vertices in each face */
   {1, 0, 2}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}
};

/* for each edge gives the corresponding face and the corresponding
   vertex on the other side of the face. */
char edge_face_vertex [NUM_TET_EDGES][2][2] = {
   {{0, 2}, {1, 3}},
   {{0, 1}, {2, 3}},
   {{0, 0}, {3, 3}},
   {{1, 1}, {2, 2}},
   {{1, 0}, {3, 2}},
   {{2, 0}, {3, 1}}
};

/*
 * TETRAHEDRON SECTION
 *
 * January 4, 1989
 * This file allows one to do direct bitmaps
 * with the triangles always in right-hand-rule order pointing in the 
 * direction of increasing scalar function, instead of pattern and 
 * permutation tables.
 *
 */

/* The first value is the number of triangles to tesselate the surface 
   within the cube.  The remaining values are the edge numbers of the 
   vertices of each of the triangles.  Edges are numbered starting from zero.   
*/

#define NUM_ENTRIES 	14 
triangulated_tets tet_cell_info [NUM_ENTRIES] = {
 { 1, {{ 0, 3, 1}, { 0, 0, 0}}}, /* 1:   1-*/
 { 1, {{ 0, 2, 4}, { 0, 0, 0}}}, /* 1:   2-*/
 { 2, {{ 1, 2, 3}, { 2, 4, 3}}}, /* 2:   3-*/
 { 1, {{ 1, 5, 2}, { 0, 0, 0}}}, /* 1:   4-*/
 { 2, {{ 0, 3, 2}, { 2, 3, 5}}}, /* 2:   5-*/
 { 2, {{ 0, 5, 4}, { 0, 1, 5}}}, /* 1:   6-*/
 { 1, {{ 3, 5, 4}, { 0, 0, 0}}}, /* 1:   7-*/
 { 1, {{ 4, 5, 3}, { 0, 0, 0}}}, /* 1:   8-*/
 { 2, {{ 0, 4, 1}, { 1, 4, 5}}}, /* 1:   9-*/
 { 2, {{ 0, 2, 5}, { 0, 5, 3}}}, /* 1:  10-*/
 { 1, {{ 2, 5, 1}, { 0, 0, 0}}}, /* 1:  11-*/
 { 2, {{ 1, 4, 2}, { 3, 4, 1}}}, /* 1:  12-*/
 { 1, {{ 4, 2, 0}, { 0, 0, 0}}}, /* 1:  13-*/
 { 1, {{ 1, 3, 0}, { 0, 0, 0}}}, /* 1:  14-*/
};

#endif
