#ifndef LT
#define LT

/* LT_POS to access a lower triangle matrix by C. Buchta
 * modified by M. Hahsler
 * n ... number of rows/columns
 * i,j ... column and row index (starts with 1)
 *
 * LT_POS1 ... 1-based indexing
 * LT_POS0 ... 0-based indexing
 */

/* All arguments are cast to a signed type of pointer size (same as R_xlen_t)
 * so the index arithmetic cannot overflow int for long dist objects
 * (i.e., more than 46,340 objects).
 */
#include <Rcpp.h>
#define LT_XLEN R_xlen_t

/* #define LT_POS1(n, i, j)					\
 * (i)==(j) ? 0 : (i)<(j) ? (n) * ((i) - 1) - (i)*((i)-1)/2 + (j)-(i) -1	\
 *     : (n)*((j)-1) - (j)*((j)-1)/2 + (i)-(j) -1
 */

#define LT_POS1(n, i, j)					\
  (((LT_XLEN)(i)) == ((LT_XLEN)(j)) ? (LT_XLEN) 0 :		\
   ((LT_XLEN)(i)) < ((LT_XLEN)(j)) ?				\
     ((LT_XLEN)(n)) * (((LT_XLEN)(i)) - 1) -			\
       ((LT_XLEN)(i)) * (((LT_XLEN)(i)) - 1) / 2 +		\
       ((LT_XLEN)(j)) - ((LT_XLEN)(i)) - 1 :			\
     ((LT_XLEN)(n)) * (((LT_XLEN)(j)) - 1) -			\
       ((LT_XLEN)(j)) * (((LT_XLEN)(j)) - 1) / 2 +		\
       ((LT_XLEN)(i)) - ((LT_XLEN)(j)) - 1)


/* #define LT_POS0(n, i, j)					\
 *  (i)==(j) ? 0 : (i)<(j) ? (n) * (i) - ((i) + 1)*(i)/2 + (j)-(i) -1	\
 *     : (n)*(j) - ((j) + 1)*(j)/2 + (i)-(j) -1
 */

#define LT_POS0(n, i, j)					\
  (((LT_XLEN)(i)) == ((LT_XLEN)(j)) ? (LT_XLEN) 0 :		\
   ((LT_XLEN)(i)) < ((LT_XLEN)(j)) ?				\
     ((LT_XLEN)(n)) * ((LT_XLEN)(i)) -				\
       (((LT_XLEN)(i)) + 1) * ((LT_XLEN)(i)) / 2 +		\
       ((LT_XLEN)(j)) - ((LT_XLEN)(i)) - 1 :			\
     ((LT_XLEN)(n)) * ((LT_XLEN)(j)) -				\
       (((LT_XLEN)(j)) + 1) * ((LT_XLEN)(j)) / 2 +		\
       ((LT_XLEN)(i)) - ((LT_XLEN)(j)) - 1)

/* M_POS to access matrix column-major order by i and j index (starts with 1)
 * n is the number of rows
 */
#define M_POS(n, i, j) ((i)+(n)*(j))


/*
 * MIN/MAX
 */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


#endif
