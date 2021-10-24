
#ifndef VECTOR_MACROS_H
#define VECTOR_MACROS_H

#define EPS      .1e-10

#define THREE_D				1

#ifdef THREE_D
#define DIMENSION          3
#else
#define DIMENSION          2
#endif /* THREE_D */

/*
#ifndef Vec
#ifdef DOUBLE
typedef double Vec [DIMENSION];
#else
typedef float Vec [DIMENSION];
#endif 
#endif *//* Vec */

#ifdef THREE_D
/* macro to print out vector A */
#define VecPrint(A) printf ("%g %g %g\n", (A)[0], (A)[1], (A)[2])

/* macro testing if vector A is within EPS of 0.0 */
#define VecNearZero(A,EPS) (((A)[0] < (EPS)) && ((A)[0] > -(EPS)) && \
((A)[1] < (EPS)) && ((A)[1] > -(EPS)) && ((A)[2] < (EPS)) && \
((A)[2] > -(EPS)))

/* macro testing if vector A is inside MN and MX bounds */
#define VecBoundCk(A,MN,MX) (((A)[0] <= (MX)[0]) && \
((A)[0] >= (MN)[0]) && ((A)[1] <= (MX)[1]) && ((A)[1] >= (MN)[1])\
&& ((A)[2] <= (MX)[2]) && ((A)[2] >= (MN)[2]))

/* macro assigning scalar to vector A, A = a */
#define VecAssign(a,A) {(A)[0] = (a); (A)[1] = (a); (A)[2] = (a);}

/* macro computing the length of vector A */
#define VecLen(A) sqrt ((A)[0] * (A)[0] + (A)[1] * (A)[1] + (A)[2] * (A)[2])

/* macro computing the dot product of vectors A and B */
#define VecDot(A,B) ((A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] * (B)[2])

/* macro computing the dot product of vectors A and B scaled by a and b */
#define VecScaleDot(a,A,b,B) (((A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] \
* (B)[2]) / ((a) * (b)))

/* macro copying vector A to vector B */
#define VecCopy(A,B) {(B)[0] = (A)[0]; (B)[1] = (A)[1]; (B)[2] = (A)[2];}

/* macro negating vector A into vector B */
#define VecNegate(A,B) {(B)[0] = -(A)[0]; (B)[1] = -(A)[1]; (B)[2] = -(A)[2];}

/* vector addition macro, C = A + B */
#define VecAdd(A,B,C) {(C)[0] = (A)[0] + (B)[0];\
(C)[1] = (A)[1] + (B)[1]; (C)[2] = (A)[2] + (B)[2];}

/* vector subtraction macro, C = A -  B */
#define VecSub(A,B,C) {(C)[0] = (A)[0] - (B)[0];\
(C)[1] = (A)[1] - (B)[1]; (C)[2] = (A)[2] - (B)[2];}

/* linear combination macro, C = aA + bB */
#define VecComb(a,A,b,B,C) {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
(C)[1] = (a)*(A)[1] + (b)*(B)[1]; (C)[2] = (a)*(A)[2] + (b)*(B)[2];}

/* linear combination summation macro, C += aA + bB */
#define VecCombSum(a,A,b,B,C) {(C)[0] += (a)*(A)[0] + (b)*(B)[0];\
(C)[1] += (a)*(A)[1] + (b)*(B)[1]; (C)[2] += (a)*(A)[2] + (b)*(B)[2];}

/* vector scale macro, B = aA */
#define VecScale(a,A,B) {(B)[0] = (a)*(A)[0];\
(B)[1] = (a)*(A)[1]; (B)[2] = (a)*(A)[2];}

/* vector scale sum macro, B += aA */
#define VecScaleSum(a,A,B) {(B)[0] += (a)*(A)[0];\
(B)[1] += (a)*(A)[1]; (B)[2] += (a)*(A)[2];}

/* add scalar multiple macro, C = aA + B */
#define VecAddS(a,A,B,C) {(C)[0] = (a)*(A)[0] + (B)[0];\
(C)[1] = (a)*(A)[1] + (B)[1]; (C)[2] = (a)*(A)[2] + (B)[2];}

/* sub scalar multiple macro, C = A - bB */
#define VecSubS(A,b,B,C) {(C)[0] = (A)[0] - (b)*(B)[0];\
(C)[1] = (A)[1] - (b)*(B)[1]; (C)[2] = (A)[2] - (b)*(B)[2];}

/* cross product macro, C = A X B */
#define VecCross(A,B,C) {(C)[0] = (A)[1] * (B)[2] - (A)[2] * (B)[1];\
(C)[1] = (A)[2] * (B)[0] - (A)[0] * (B)[2];\
(C)[2] = (A)[0] * (B)[1] - (A)[1] * (B)[0];}

/* add scalar multiple macro, D = (A + B + C) / 3.0 */
#define VecAvg(A,B,C,D) {(D)[0] = ((A)[0] + (B)[0] + (C)[0]) / 3.0;\
(D)[1] = ((A)[1] + (B)[1] + (C)[1]) / 3.0; \
(D)[2] = ((A)[2] + (B)[2] + (C)[2]) / 3.0;}

/* triple scalar product macro, D = Det [A,B,C] (each Vec is row)
   definition taken from "Introduction to Vector Analysis" p45 */
#define VecDet3x3(A,B,C) ((A)[0] * ((B)[1] * (C)[2] - (B)[2] * (C)[1])\
+ (A)[1] * ((B)[2] * (C)[0] - (B)[0] * (C)[2]) + \
(A)[2] * ((B)[0] * (C)[1] - (B)[1] * (C)[0]))

/**********************************************************************************/
#else /* 2D */
/* macro to print out vector A */
#define VecPrint(A) printf ("%f %f\n", (A)[0], (A)[1])

/* macro testing if vector A is near 0.0 */
#define VecNearZero(A,EPS) (((A)[0] < (EPS)) && ((A)[0] > -(EPS)) && \
((A)[1] < (EPS)) && ((A)[1] > -(EPS)))

/* macro testing if vector A is inside MN and MX bounds */
#define VecBoundCk(A,MN,MX) (((A)[0] <= (MX)[0]) && ((A)[0] >= (MN)[0]) && \
((A)[1] <= (MX)[1]) && ((A)[1] >= (MN)[1]))

/* macro assigning scalar to vector A, A = a */
#define VecAssign(a,A) {(A)[0] = (a); (A)[1] = (a);}

/* macro computing the length of vector A */
#define VecLen(A) sqrt ((A)[0] * (A)[0] + (A)[1] * (A)[1])

/* macro computing the dot product of vectors A and B */
#define VecDot(A,B) ((A)[0] * (B)[0] + (A)[1] * (B)[1])

/* macro copying vector A to vector B */
#define VecCopy(A,B) {(B)[0] = (A)[0]; (B)[1] = (A)[1];}

/* macro negating vector A into vector B */
#define VecNegate(A,B) {(B)[0] = -(A)[0]; (B)[1] = -(A)[1];}
/* vector addition macro, C = A + B */
#define VecAdd(A,B,C) {(C)[0] = (A)[0] + (B)[0]; (C)[1] = (A)[1] + (B)[1];}

/* vector subtraction macro, C = A -  B */
#define VecSub(A,B,C) {(C)[0] = (A)[0] - (B)[0]; (C)[1] = (A)[1] - (B)[1];}

/* linear combination macro, C = aA + bB */
#define VecComb(a,A,b,B,C) {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
(C)[1] = (a)*(A)[1] + (b)*(B)[1];}

/* linear combination summation macro, C += aA + bB */
#define VecCombSum(a,A,b,B,C) {(C)[0] += (a)*(A)[0] + (b)*(B)[0];\
(C)[1] += (a)*(A)[1] + (b)*(B)[1];}

/* vector scale macro, B = aA */
#define VecScale(a,A,B) {(B)[0] = (a)*(A)[0]; (B)[1] = (a)*(A)[1];}

/* vector scale sum macro, B += aA */
#define VecScaleSum(a,A,B) {(B)[0] += (a)*(A)[0]; (B)[1] += (a)*(A)[1];}

/* add scalar multiple macro, C = aA + B */
#define VecAddS(a,A,B,C) {(C)[0] = (a)*(A)[0] + (B)[0];\
(C)[1] = (a)*(A)[1] + (B)[1];}

/* sub scalar multiple macro, C = A - bB */
#define VecSubS(A,b,B,C) {(C)[0] = (A)[0] - (b)*(B)[0];\
(C)[1] = (A)[1] - (b)*(B)[1];}

#endif /* 3D */


#endif
