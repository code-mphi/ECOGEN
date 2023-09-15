/*
 * HISTORY:
 * This is an SVD library slightly modified by Kevin Schmidmayer.
 * The original version was aquired from
 * https://www.gb.nrao.edu/~rcreager/Zpectrometer/.
*/

/*
  $Id: svd.c,v 1.4 2004/03/24 04:58:59 rauch Exp $

  Routines related to singular value decomposition of matrices, based on
  routines in Numerical Recipes which have been modified: to accept
  standard 0-based, row-order C arrays as input; to return error codes if
  memory could not be allocated; to use double precision instead of floats.
  svsolve() and svinverse() functions have also been added.

  Note that all matrix dimensions are specified as the number of rows
  followed by the number of columns.
*/
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"

#define MAXIT 32
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Elements of static temporary storage to allocate. */
#define NSTAT 16

/* Default value for smallest meaningful singular value. */
#define WMIN 1e-13

inline int svdcmp(double u[/*m x n*/], double w[/*n*/], double v[/*n x n*/],
 const double a[/*m x n*/], int m, int n)
{
/*
  Compute the SVD A = U x W x V^T of the m (rows) by n (columns) matrix A;
  on output, u[] holds the m by n matrix U, w[] holds the n diagonal 
  elements of the (diagonal) n by n matrix W, and v[] holds the transpose of
  the n by n matrix V^T; note that u[] can be identical to a[].
  Returns 0 on success, or a system error code otherwise
  (ERANGE if MAXIT iterations are reached without converging).
*/
  static double srv1[NSTAT];
  int flag, i, its, j, jj, k, l=0, nm=0;
  double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

  if (n<=NSTAT) rv1=srv1; else rv1=(double *)malloc(n*sizeof(double));
  if (!rv1) return(errno);

  if (u!=a) memcpy(u, a, m*n*sizeof(double));
  g=scale=anorm=0;

  for (i=0; i<n; i++) {
    l=i+1; rv1[i]=scale*g; g=s=scale=0;
    if (i<m) {
      double *uki=u+i*(n+1);

      for (k=i; k<m; k++, uki+=n) scale+=fabs(*uki);
      if (scale) {
        for (uki=u+i*(n+1), k=i; k<m; k++, uki+=n) {
          *uki/=scale; s+=(*uki)*(*uki);
        }
        f=u[i*n+i]; g=-SIGN(sqrt(s),f); h=f*g-s; u[i*n+i]=f-g;
        for (j=l; j<n; j++) {
          double *ukj;

          for (s=0, uki=u+i*(n+1), ukj=u+i*n+j, k=i; k<m; k++, uki+=n, ukj+=n)
            s+=(*uki)*(*ukj);
          f=s/h;
          for (uki=u+i*(n+1), ukj=u+i*n+j, k=i; k<m; k++, uki+=n, ukj+=n)
            *ukj+=f*(*uki);
        }
        for (uki=u+i*(n+1), k=i; k<m; k++, uki+=n) (*uki)*=scale;
      }
    }
    w[i]=scale*g; g=s=scale=0;
    if (i<=m && i!=n) {
      double *ui=u+i*n;

      for (k=l; k<n; k++) scale+=fabs(ui[k]);
      if (scale) {
        for (ui=u+i*n+l, k=l; k<n; k++, ui++) { (*ui)/=scale; s+=(*ui)*(*ui); }
        f=u[i*n+l]; g=-SIGN(sqrt(s),f); h=f*g-s; u[i*n+l]=f-g;
        for (ui=u+i*n, k=l; k<n; k++) rv1[k]=ui[k]/h;
        for (j=l; j<m; j++) {
          double *uj=u+j*n;

          for (s=0, k=l; k<n; k++) s+=uj[k]*ui[k];
          for (uj=u+j*n, k=l; k<n; k++) uj[k]+=s*rv1[k];
        }
        for (ui=u+i*n, k=l; k<n; k++) ui[k]*=scale;
      }
    }
    { double tmp=fabs(w[i])+fabs(rv1[i]); if (tmp>anorm) anorm=tmp; }
  }

  for (i=n-1; i>=0; i--) {
    if (i<n-1) {
      if (g) {
        double *vji=v+l*n+i, *ui=u+i*n, g1=ui[l]*g;

        if (g1) {
          g1=1/g1; for (j=l; j<n; j++, vji+=n) *vji=g1*ui[j];
        } else {
          g1=1/ui[l]; for (j=l; j<n; j++, vji+=n) *vji=(ui[j]*g1)/g;
        }
        for (j=l; j<n; j++) {
          double *vkj=v+l*n+j, *vki=v+l*n+i;

          for (s=0, ui=u+i*n, k=l; k<n; k++, vkj+=n) s+=ui[k]*(*vkj);
          for (vkj=v+l*n+j, k=l; k<n; k++, vkj+=n, vki+=n) (*vkj)+=s*(*vki);
        }
      }
      {
        double *vi=v+i*n, *vji=v+l*n+i;

        for (j=l; j<n; j++, vji+=n) vi[j]=(*vji)=0;
      }
    }
    v[i*n+i]=1; g=rv1[i]; l=i;
  }

  for (i=(m>n ? n : m)-1; i>=0; i--) {
    l=i+1; g=w[i];
    { double *ui=u+i*n;  for (j=l; j<n; j++) ui[j]=0; }
    if (g) {
      g=1/g;
      for (j=l; j<n; j++) {
        double *uki=u+l*n+i, *ukj=u+l*n+j;

        for (s=0, k=l; k<m; k++, uki+=n, ukj+=n) s+=(*uki)*(*ukj);
        f=(s/u[i*n+i])*g;
        for (uki=u+i*n+i, ukj=u+i*n+j, k=i; k<m; k++, uki+=n, ukj+=n)
          (*ukj)+=f*(*uki);
      }
      { double *uji=u+i*n+i;  for (j=i; j<m; j++, uji+=n) (*uji)*=g; }
    } else {
      double *uji=u+i*n+i; 
      
      for (j=i; j<m; j++, uji+=n) *uji=0;
    }
    ++u[i*n+i];
  }

  for (k=n-1; k>=0; k--) {
    for (its=0; its<=MAXIT; its++) {
      flag=1;
      for (l=k; l>=0; l--) {
        nm=l-1; if (fabs(rv1[l])+anorm==anorm) { flag=0; break; }
        if (fabs(w[nm])+anorm==anorm) break;
      }
      if (flag) {
        c=0; s=1;
        for (i=l; i<=k; i++) {
          f=s*rv1[i]; rv1[i]=c*rv1[i];
          if (fabs(f)+anorm==anorm) break;
          g=w[i]; h=hypot(f, g); w[i]=h; h=1/h; c=g*h; s=-f*h;
          {
            double *ujnm=u+nm, *uji=u+i;

            for (j=0; j<m; j++, ujnm+=n, uji+=n) {
              y=(*ujnm); z=(*uji); (*ujnm)=y*c+z*s; (*uji)=z*c-y*s;
            }
          }
        }
      }
      z=w[k];
      if (l==k) {
        if (z<0) {
          double *vjk=v+k;

          w[k]= -z; for (j=0; j<n; j++, vjk+=n) (*vjk)= -(*vjk);
        }
        break;
      }

      if (its==MAXIT) { if (n>NSTAT) free(rv1);   return(ERANGE); }

      x=w[l]; nm=k-1; y=w[nm]; g=rv1[nm]; h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
      g=hypot(f, 1); f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x; c=s=1;
      for (j=l; j<=nm; j++) {
        i=j+1; g=rv1[i]; y=w[i]; h=s*g; g=c*g; z=hypot(f,h);
        rv1[j]=z; 
        { double z1=1/z;  c=f*z1; s=h*z1; }
        f=x*c+g*s; g=g*c-x*s; h=y*s; y*=c;
        {
          double *vjjj=v+j, *vjji=v+i;

          for (jj=0; jj<n; jj++, vjjj+=n, vjji+=n) {
            x=(*vjjj); z=(*vjji); (*vjjj)=x*c+z*s; (*vjji)=z*c-x*s;
          }
        }
        z=hypot(f,h); w[j]=z;
        if (z) { z=1/z; c=f*z; s=h*z; }
        f=c*g+s*y; x=c*y-s*g;
        {
          double *ujjj=u+j, *ujji=u+i;

          for (jj=0; jj<m; jj++, ujjj+=n, ujji+=n) {
            y=(*ujjj); z=(*ujji); (*ujjj)=y*c+z*s; (*ujji)=z*c-y*s;
          }
        }
      }
      rv1[l]=0; rv1[k]=f; w[k]=x;
    }
  }

  if (n>NSTAT) free(rv1);
  return(0);
}

inline int svbksb(double x[/*n x p*/], const double u[/*m x n*/],
  const double w[/*n*/], const double v[/*n x n*/],
  const double b[/*m x p*/], int m, int n, int p)
{
/*
   Solves A X = B for X, where A = U x W x V^T is the SVD of the m by n
   matrix A, B is an m by p matrix, and X is an n by p matrix.
   Returns 0 on success, or a system error code otherwise.
*/
  int jj, j, i, k;
  static double stmp[NSTAT];
  double s, *tmp=(n<=NSTAT ? stmp : (double *)malloc(n*sizeof(double))), *xjk;

  if (!tmp) return(errno);

  for (k=0; k<p; k++) {
    for (j=0; j<n; j++) {
      s=0;
      if (w[j]) {
        const double *uij=u+j, *bik=b+k;

        for (i=0; i<m; i++, uij+=n, bik+=p) s+=(*uij)*(*bik);
        s/=w[j];
      }
      tmp[j]=s;
    }
    for (xjk=x+k, j=0; j<n; j++, xjk+=p) {
      const double *vj=v+j*n;

      s=0; for (jj=0; jj<n; jj++) s+=vj[jj]*tmp[jj];
      (*xjk)=s;
    }
  }

  if (n>NSTAT) free(tmp);
  return(0);
}

inline int svsolve(double x[/*n x p*/], double wmin, double a[/*m x n*/],
  const double b[/*m x p*/], int m, int n, int p, int preserve)
{
/*
   Solve a (least squares) system of linear equations A X = B by singular value 
   decomposition. Zero all singular values smaller than wmin*max(w[]), using
   wmin=WMIN if wmin<0. If preserve is non-zero, A is copied before
   decomposition; A is m x n, B is m x p, and X will be n x p.
   Returns 0 on success, or a system error code otherwise.
*/
  double *utmp, *w, *v, wmax;
  int i, status;

  w=(double *)malloc(n*sizeof(double));
  v=(double *)malloc(n*n*sizeof(double));
  utmp=(preserve ? (double *)malloc(m*n*sizeof(double)) : a);
  if (!w || !v || !utmp) {
    if (w) free(w);
    if (v) free(v);
    if (utmp && preserve) free(utmp);
    return(errno);
  }

  status=svdcmp(utmp, w, v, a, m, n);
  if (status) {
    free(w); free(v); if (preserve) free(utmp);
    return(status);
  }

  for (wmax=0, i=0; i<n; i++) if (wmax<w[i]) wmax=w[i];
  if (wmin<0) wmax*=WMIN; else wmax*=wmin;
  for (i=0; i<n; i++) if (w[i]<wmax) w[i]=0;

  status=svbksb(x, utmp, w, v, b, m, n, p);
  if (status) {
    free(w); free(v); if (preserve) free(utmp);
    return(status);
  }

  free(w); free(v); if (preserve) free(utmp);
  return(0);
}

inline int svinverse(double x[/*n x m*/], double wmin, const double a[/*m x n*/],
  int m, int n)
{
/* 
   Find the pseudoinverse of the m x n matrix A, zeroing singular values below
   wmin*max(w[]). X should have space for an n x m matrix and can be the
   same as A.  Returns 0 on success, or a system error code otherwise.
*/
  double *u=(double *)malloc(n*m*sizeof(double)),
    *w=(double *)malloc(n*sizeof(double)),
    *v=(double *)malloc(n*n*sizeof(double)), wmax;
  int i, j, k, status;

  if (!u || !w || !v) {
    if (u) free(u);
    if (w) free(w);
    if (v) free(v);
    return(errno);
  }

  status=svdcmp(u, w, v, a, m, n);
  if (status) { free(u); free(w); free(v); return(status); }

  for (wmax=0, i=0; i<n; i++) if (wmax<w[i]) wmax=w[i];
  if (wmin<0) wmax*=WMIN; else wmax*=wmin;
  for (i=0; i<n; i++) if (w[i]<wmax) w[i]=0; else w[i]=1/w[i];
  {
    double *vi=v, *uj;

    for (i=0; i<n; i++, vi+=n)
      for (uj=u, j=0; j<m; j++, uj+=n, x++)
        for (*x=0, k=0; k<n; k++) (*x)+=vi[k]*w[k]*uj[k];
  }

  free(u); free(w); free(v);
  return(0);
}