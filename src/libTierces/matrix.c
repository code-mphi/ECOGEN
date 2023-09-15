/*
 * HISTORY:
 * This is a matrix C file for the SVD library slightly modified by Kevin Schmidmayer.
 * The original version was aquired from
 * https://www.gb.nrao.edu/~rcreager/Zpectrometer/.
*/

/*
 *  some elementary (2d) matrix operations - also contain
 *  some #define'd sections of code to benchmark
 *  matrix acccess AA[][] vs. A[]
 *
 *  On a P600 the typical cost to operate on a 100x100 matrix is
 *    0.2-0.3ms on elemental functions (allocate, transpose, dup)
 *    1.2ms     on a laplacian
 *    8-12ms    on a matrix multiplication
 *    18ms      on a matrix inversion
 *
 *
 * Some random references:
 *  http://math.nist.gov/javanumerics/jama/         a JAVA matrix package
 */

#include <nemo.h>
#include "matrix.h"

extern double xrandom(double, double);

static char *matrix_version = "$Id: matrix.c,v 1.16 2003/05/09 05:00:31 pteuben Exp $";

Matrix *allocate_matrix(int nr, int nc)
{
  Matrix *m;

  if (nr<0) error("allocate_matrix: illegal size nr=%d",nr);
  if (nc<0) error("allocate_matrix: illegal size nr=%d",nc);
  m = (Matrix *) allocate(sizeof(Matrix));
  m->dat = (double *) NULL;
  m->p = (double **) NULL;
  m->nc = m->nr = 0;
  m->ac = m->ar = (Axis2 *) NULL;
  reallocate_matrix(m,nr,nc);    /* actually allocate the data */

  return m;
}

void reallocate_matrix(Matrix *m, int nr, int nc)
{
  double *d, **p;
  int i;

  if (m->nr == nr && m->nc == nc) return;      /* nothing to do */

  if (m->p)   free(m->p);       /* free up old data */
  if (m->dat) free(m->dat);

  m->nr = nr;
  m->nc = nc;
  if (nr==0 && nc==0) return;             /* return with dummy empty matrix */
  d = m->dat = (double *) allocate(nr*nc*sizeof(double));   /* contigues block of data */
  p = m->p = (double **) allocate(sizeof(double *)*nr);   /* array of pointers */
  nemo_dprintf(1,"reallocate_matrix: @0x%x 0x%x 0x%x\n",m,d,p);
  for (i=0; i<nr; i++) {
    p[i] = d;
    d += nc;
  }
  zero_matrix(m);           /* all data elements are zeroed */
}


void free_matrix(Matrix *m)
{
  if (m == NULL) return;
  nemo_dprintf(1,"free_matrix:     @0x%x 0x%x 0x%x\n",m,m->dat,m->p);
  free(m->p);
  free(m->dat);
  m->nc = -1;
  m->nr = -1;
  free(m->ar);
  free(m->ac);
  free(m);
}

Axis2   *
new_axis(double crpix, double crval, double cdelt, string ctype, string cunit)
{
  Axis2 *a = (Axis2 *) allocate(sizeof(Axis2));
  a->crpix = crpix;
  a->crval = crval;
  a->cdelt = cdelt;
  a->ctype = (ctype ? strdup(ctype) : NULL);
  a->cunit = (cunit ? strdup(cunit) : NULL);
  return a;

}

Axis2   *
new_axis1(Axis2 *o)
{
  Axis2 *a;

  if (o == NULL) return NULL;

  a = (Axis2 *) allocate(sizeof(Axis2));
  a->crpix = o->crpix;
  a->crval = o->crval;
  a->cdelt = o->cdelt;
  a->ctype = (o->ctype ? strdup(o->ctype) : NULL);
  a->cunit = (o->cunit ? strdup(o->cunit) : NULL);
  return a;
}

void free_axis(Axis2 *a)
{
  if (a == NULL) return;
  if (a->ctype) free(a->ctype);
  if (a->cunit) free(a->cunit);
  free(a);
}


void check_matrix(const Matrix *m, int nr, int nc, const char *msg)
{
  if (m->nr != nr || m->nc != nc)
    error("%s: found matrix[%d][%d], expected [%d][%d]",
          msg,m->nr,m->nc,nr,nc);
}


void zero_matrix(Matrix *z)
{
  double *dp;

  int n = z->nr * z->nc;
  dp = z->dat;
  while (n--)
    *dp++ = 0.0;
}


void random_matrix(Matrix *z)
{
  double *dp;

  int n = z->nr * z->nc;
  dp = z->dat;
  while (n--)
    *dp++ = xrandom(0.0,1.0);
}

/*
 * create a unit matrix, or something close to it.
 *      columns > nr or rows > nc are left at 0.0
 */

void unit_matrix(Matrix *u)
{
  int i, n=MIN(u->nr,u->nc);

  double **p = u->p;
  zero_matrix(u);
  for (i=0; i<n; i++)
    p[i][i] = 1.0;
  
}

Matrix *dup_matrix(const Matrix *a)
{
  int i, j;
  Matrix *c;

  c = allocate_matrix(a->nr,a->nc);
  for (i=0; i<a->nr; i++)
    for (j=0; j<a->nc; j++)
      c->p[i][j] = a->p[i][j];

  return c;
}

  

/*
 *    add two matrices
 *           c = a + b
 *    if c==0, create a new one and return it
 */

Matrix *add_matrix(const Matrix *a, const Matrix *b, Matrix *c)
{
  int i, j;

  if (a->nr != b->nr) 
    error("add_matrix: nr not equal: a[%d][%d] b[%d][%d]",a->nr,a->nc,b->nr,b->nc);
  if (a->nc != b->nc) 
    error("add_matrix: nc not equal: a[%d][%d] b[%d][%d]",a->nr,a->nc,b->nr,b->nc);

  if (c == NULL)
    c = allocate_matrix(a->nr,a->nc);
  else
    if (a->nr != c->nr || a->nc != c->nc)
      error("add_matrix: A and C not equal: a[%d][%d] c[%d][%d]",a->nr,a->nc,c->nr,c->nc);
  for (i=0; i<a->nr; i++)
    for (j=0; j<a->nc; j++)
      c->p[i][j] = a->p[i][j] + b->p[i][j]; 

  return c;
}

/* 
 *     a += b (scalar)
 */

void adds_matrix(Matrix *a, const double b)
{
  int n = a->nr * a->nc;

  double *dp = a->dat;
  while (n--)
    *dp++ += b;
}


/*
 *  COST: 100x100 at about 8-12ms/call (depending on access method)
 *
 *   c = a * b
 *
 *   if c==0, allocate new and return it 
 */

#if 1
#define _use_p  1
#endif

Matrix *mul_matrix(const Matrix *a, const Matrix *b, Matrix *c)
{
  int i,j,k;
  double s;
#ifndef _use_p
  double **ap, **bp, **cp;
#else
  double *ad, *bd, *cd;
  double *d1, *d2;
#endif

    
  if (a->nc != b->nr) 
    error("mul_matrix: cannot a[%d][%d] x b[%d][%d]",a->nr,a->nc,b->nr,b->nc);
  if (c == NULL)
    c = allocate_matrix(a->nr,b->nc);

  /* 1000 times 100x100 matrix :  11.9" at -O2 */
#ifndef _use_p
  ap = a->p;
  bp = b->p;
  cp = c->p;
  for (i=0; i<c->nr; i++) {
    for (j=0; j<c->nc; j++) {
      s = 0;
      for (k=0; k<a->nc; k++) {
        s+= ap[i][k] * bp[k][j];
      }
      cp[i][j] = s;
    }
  }
#else
  /*  1000 times 100x100 matrix:  8.0" at -02 */
  ad = a->dat;
  bd = b->dat;
  cd = c->dat;
  for (j=0; j<c->nc; j++) {
    for (i=0; i<c->nr; i++) {
      s = 0;
      d1 = &ad[i*a->nc]; 
      d2 = &bd[j];
      for (k=0; k<a->nc; k++) {
        s+= (*d1) * (*d2);
        d1++;
        d2 += b->nc;
      }
      cd[i*c->nc+j] = s;
    }
  }
#endif
  c->ar = new_axis1(a->ar);
  c->ac = new_axis1(b->ac);
  return c;
}

/*
 * IMUL:  C = A * B, where B is really a diagonal matrix (useful for SVD work)
 *        B[1][n] is treated as a diagonal square matrix B[n][n]
 *
 *        Only useful to speed up SVD type multiplications - an interface hack
 *      
 */

Matrix *imul_matrix(const Matrix *a, const Matrix *b, Matrix *c)
{
  int i,j;
  double **ap, **bp, **cp;
    
  if (b->nr != 1) 
    error("imul_matrix: expect b[1][%d], got b[%d][%d]",b->nr,b->nr,b->nc);
  if (c == NULL)
    c = allocate_matrix(a->nr,b->nc);

  ap = a->p;
  bp = b->p;
  cp = c->p;
  for (i=0; i<c->nr; i++) {
    for (j=0; j<c->nc; j++) {
      cp[i][j] = ap[i][j] * bp[0][j];
    }
  }
  c->ar = new_axis1(a->ar);
  c->ac = new_axis1(b->ac);
  return c;
}

void vmul_matrix(const double *v, const Matrix *a, double *b)
{
  int i,j;
  double s;
  double **ap = a->p;

  for (i=0; i<a->nc; i++) {
    for (j=0, s=0.0; j<a->nr; j++)
      s+= v[j]*ap[j][i];
    b[i] = s;
  }
}

/* 
 *  multiply a matrix with a scalar
 *     a *= b (scalar)
 */

void muls_matrix(Matrix *a, const double b)
{
  int n = a->nr * a->nc;

  double *dp = a->dat;
  while (n--)
    *dp++ *= b;
}

/*
 * print a matrix, with a label
 *     a[0][0]  a[0][1] .....
 *     a[1][0]  ...
 *     ....
 */

void print_matrix(const Matrix *a, const char *msg)
{
  int i,j;
  double **p = a->p;

  printf("============================\n%s[%d][%d]\n",msg,a->nr,a->nc);
  for (j=0; j<a->nr; j++) {
    for (i=0; i<a->nc; i++) {
      printf("%g ",p[j][i]);
    }
    printf("\n");
  }
  printf("============================\n");
}


/*  
 *  COST:  100x100 at 0.26ms/call
 *
 *  100x100  ::  1000 in  0.8"  wo/ free_matrix
 *             100000 in 26.4"  w/ free_matrix   (-O2)  
 *
 */

Matrix *transpose_matrix(const Matrix *a)
{
  Matrix *mp;
  int i,j;
  double **p, **q;

  mp = allocate_matrix(a->nc, a->nr);
  p = a->p;
  q = mp->p;

  for (j=0; j<a->nc; j++)
    for (i=0; i<a->nr; i++)
      q[j][i] = p[i][j];
  mp->ac = a->ar;
  mp->ar = a->ac;
  return mp;
}

/*
 *  an example matrix operation to show when pointer arithmetic can
 *       a) look ugly
 *       b) be slow
 *  compared to a[i][j] notation
 */

Matrix *laplace_matrix(const Matrix *a)
{
  Matrix *mp;
  int i,j,k,nc;
  double **p, **q, d1,d2,d3,d4;
  double *pd, *qd;

  mp = allocate_matrix(a->nr, a->nc);

#if 1
  /* COST: 100x100 at about 1.2ms/call */
  p = a->p;
  q = mp->p;
  for (j=1; j<a->nc-1; j++)
    for (i=1; i<a->nr-1; i++) {
      d1 = p[i][j]-p[i-1][j];
      d2 = p[i][j]-p[i+1][j];
      d3 = p[i][j]-p[i][j-1];
      d4 = p[i][j]-p[i][j+1];
      q[i][j] = sqrt(d1*d1+d2*d2+d3*d3+d4*d4);
    }
#else
  /* COST: 100x100 at about 1.2ms/call */
  pd = a->dat;
  qd = mp->dat;
  nc = a->nc;
  for (j=1; j<a->nc-1; j++)
    for (i=1; i<a->nr-1; i++) {
      k = i*nc + j;
      d1 = pd[k]-pd[k-1];
      d2 = pd[k]-pd[k+1];
      d3 = pd[k]-pd[k+nc];
      d4 = pd[k]-pd[k-nc];
      qd[k] = sqrt(d1*d1+d2*d2+d3*d3+d4*d4);
    }
#endif
  mp->ac = a->ac;
  mp->ar = a->ar;
  return mp;
}

/*
 *  return a matrix that can be used for binning by right-hand 
 *  matrix multiplication, i.e.
 *      M[nx][ny] * B[ny][nbin] = MB[nx][nbin]
 * 
 */

Matrix *bin_matrix(const Matrix *a, const int nbin, const int *bin0, const int *bin1)
{
  Matrix *mp;
  double **p;
  int i, j;

  if (nbin > a->nr) error("Cannot bin this way, too many");
  mp = allocate_matrix(a->nr, nbin);

  zero_matrix(mp);
  p = mp->p;
  for (j=0; j<nbin; j++) {
    for (i=bin0[j]; i<=bin1[j]; i++) {
      p[i][j] = 1.0/(bin1[j]-bin0[j]);
    }
  }
  /* note no axes here yet */
  return mp;
}

/* 
 * Gauss-Jordan matrix inversion (Stoer, I), see NEMO's matinv.c code
 * the words column and row can be interchanged, as this code (and
 * comments) was taken from a C/Fortran-style code.
 *
 * This code only works well for well conditioned matrices
 *
 * COST:  100x100 about 18.2 ms/call
 */
double invert_matrix(Matrix *a)
{
  int    i,j,row,k,evin,size;
  double max,even,mjk,det;
  double **matrix;
  static int maxcol = -1;    /* static members to keep some local arrays */
  static int *per;           /* permutation order array */
  static double *hv;         /* temp column for column permutations */

  if (a == NULL) {           /* option to free static local memory  */
    if (maxcol > 0) {
      free(per);
      free(hv);
      maxcol = -1;
      return -1.0;
    }
    error("invert_matrix: null matrix");
  }

  if (a->nc != a->nr) 
    error("matrix[%d][%d] not square",a->nr,a->nc);

  size = a->nc;
  matrix = a->p;

  if (maxcol < 0) {
    nemo_dprintf(1,"invert_matrix: allocating %d\n",size);
    per = (int *) allocate(sizeof(int)*size);
    hv = (double *) allocate(sizeof(double)*size);
    maxcol = size;
  } else if (size > maxcol) {
    nemo_dprintf(1,"invert_matrix: incrementing allocation to %d\n",size);
    per = (int *) reallocate(per, sizeof(int)*size);
    hv = (double *) reallocate(hv, sizeof(double)*size);
    maxcol = size;
  }

  det=1.0;    /* in case of normal end */
    
  for (i=0; i<size; i++)    /* set permutation array */
    per[i]=i;

  for (j=0; j<size; j++) { /* in j-th column, set row with largest element */
    max=ABS(matrix[j][j]);
    row=j;
    i=j+1;
    while (i < size) {
      if (ABS(matrix[j][i]) > max) {
        max=ABS(matrix[j][i]);
        row=i;
      }
      i++;
    }
    if (matrix[j][row] == 0.0)      /* determinant zero ? */
      return 0.0;                   /* no solution possible */
       
    if (row > j) {  /* if largest element not on diagonal: */
      /* then permute rows */
      for (k=0; k<size; k++) {  /*    permutation loop */
        even=matrix[k][j];
        matrix[k][j]=matrix[k][row];
        matrix[k][row]=even;
      }
      evin=per[j];        /* keep track of permutation */
      per[j]=per[row];
      per[row]=evin;
    }
    even=1.0/matrix[j][j];       /* modify column */
    for(i=0; i<size; i++) 
      matrix[j][i]=even*matrix[j][i];

    matrix[j][j]=even;
    k=0;            /* modify rest of matrix */
    while (k < j) {
      mjk=matrix[k][j];
      i=0;
      while (i < j) {
        matrix[k][i] -= matrix[j][i]*mjk;
        i++;
      }
      i=j+1;
      while (i < size) {
        matrix[k][i] -= matrix[j][i]*mjk;
        i++;
      }
      matrix[k][j] = -even*mjk;
      k++;
    }
    k=j+1;
    while (k < size) {
      mjk=matrix[k][j];
      i=0;
      while (i < j) {
        matrix[k][i] -= matrix[j][i]*mjk;
        i++;
      }
      i=j+1;
      while (i < size) {
        matrix[k][i] -= matrix[j][i]*mjk;
        i++;
      }
      matrix[k][j] = -even*mjk;
      k++;
    }
  }  /*    end of loop through columns */

  for (i=0; i<size; i++) {  /* finally, repermute columns (rows really here) */
    for(k=0; k<size; k++)
      hv[per[k]]= matrix[k][i];
    for(k=0; k<size; k++)
      matrix[k][i] = hv[k];
  }
  if (1) {
    Axis2 *tmp;
    tmp = a->ac;
    a->ac = a->ar;
    a->ar = tmp;
  }

  return det;    /* but note, it's not even computed yet !!!! */
}

/*
 *  refactor matrix A[m][n]  into   U[m][n] * W[n][n] * V^T[n][n]
 *  (W is actually stored as vector W[n] for compactness)
 *  such that
 *            V[n][n] * W[n][n] * U ^T[n][m] * A[m][n] = I[n][n]
 */

void svd_matrix(const Matrix *a, Matrix *u, Matrix *v, Matrix *w)
{
  int nr = a->nr;
  int nc = a->nc;
  int i;
  double **wp;

#if 0
  print_matrix(a,"a-before"); 
#endif
  
  reallocate_matrix(u,nr,nc);
  reallocate_matrix(v,nc,nc);
#if 1
  /* production code ??? */
  reallocate_matrix(w,1,nc);
  /* Kevin's NR-emulator */
  svdcmp(u->dat, w->dat, v->dat, a->dat, a->nr, a->nc);
#else
  warning("svd_matrix: using W as matrix");
  /* only for testing if W is needed as matrix */
  /* various other routines will fail if W is a matrix, and not a vector */
  /* because of the NR interfaces and our inability to do left hand matrix */
  /* multiplications */
  reallocate_matrix(w,nc,nc);  /* !! note !! */

  /* Kevin's NR-emulator */
  svdcmp(u->dat, w->dat, v->dat, a->dat, a->nr, a->nc);
  /* the first row of the W matrix needs to become its diagonal */
  wp = w->p; 
  for (i=1; i<nc; i++) {    /* but skip 1st element !! */
    wp[i][i] = wp[0][i];
    wp[0][i] = 0.0;
  }
#endif
  
#if 0
  print_matrix(u,"U");
  print_matrix(w,"W");
  print_matrix(v,"V");
#endif
  u->ar = new_axis1(a->ar);
  u->ac = new_axis1(a->ac);
  v->ar = new_axis1(a->ac);
  v->ac = new_axis1(a->ac);
  w->ac = new_axis1(a->ac);
}

/*
 *  solve Ax=b, assuming you have factored A into their U,V,W svd matrices
 *              x and b can contain multiple columns (p)
 */
void svd_solve(const Matrix *a, const Matrix *u, const Matrix *v, const Matrix *w,
               Matrix *b, Matrix *x)
{
  int m = a->nr;
  int n = a->nc;
  int p = b->nc;

  check_matrix(a,m,n,"svd_solve: A");
  check_matrix(u,m,n,"svd_solve: U");
  check_matrix(v,n,n,"svd_solve: V");
  check_matrix(w,n,n,"svd_solve: W");
  check_matrix(b,n,p,"svd_solve: B");
  check_matrix(x,n,p,"svd_solve: X");
  

  /* call Kevin's NR routine */
  svbksb(x->dat, u->dat, w->dat, v->dat, b->dat, m, n, p);
}