/*
 * HISTORY:
 * This is a matrix header for the SVD library slightly modified by Kevin Schmidmayer.
 * The original version was aquired from
 * https://www.gb.nrao.edu/~rcreager/Zpectrometer/.
*/

typedef struct _axis {
  double crpix;
  double crval;
  double cdelt;
  char  *ctype;
  char  *cunit;
} Axis2;

typedef struct _vector {
  int n;
  double *dat;
  Axis2 *a;
} Vector;

typedef struct _matrix {
  int nr;                   
  int nc;            
  double *dat;       
  double **p;               
  Axis2 *ar, *ac;
} Matrix;


/* svd.c routines. */
extern int
svdcmp(double u[/*m x n*/], double w[/*n*/], double v[/*n x n*/],
         const double a[/*m x n*/], int m, int n);
extern int
svbksb(double x[/*n x p*/], const double u[/*m x n*/],
         const double w[/*n*/], const double v[/*n x n*/],
         const double b[/*m x p*/], int m, int n, int p);
extern int
svsolve(double x[/*n x p*/], double wmin, double a[/*m x n*/],
          const double b[/*m x p*/], int m, int n, int p, int preserve);
extern int
svinverse(double x[/*n x m*/], double wmin, const double a[/*m x n*/],
            int m, int n);

/* cproto -I../include -I$NEMOLIB -I$NEMOINC matrix.c */

Matrix *allocate_matrix(int nr, int nc);
void    reallocate_matrix(Matrix *m, int nr, int nc);
void    free_matrix(Matrix *m);
Axis2   *new_axis(double crpix, double crval, double cdelt, char *ct, char *cu);
Axis2   *new_axis1(Axis2 *a);
void    check_matrix(const Matrix *m, int nr, int nc, const char *msg);
void    zero_matrix(Matrix *z);
void    random_matrix(Matrix *z);
void    unit_matrix(Matrix *u);
Matrix *dup_matrix(const Matrix *a);
Matrix *add_matrix(const Matrix *a, const Matrix *b, Matrix *c);
void    adds_matrix(Matrix *a, const double b);
Matrix *mul_matrix(const Matrix *a, const Matrix *b, Matrix *c);
Matrix *imul_matrix(const Matrix *a, const Matrix *b, Matrix *c);
void    vmul_matrix(const double *v, const Matrix *a, double *b);
void    muls_matrix(Matrix *a, const double b);
void    print_matrix(const Matrix *a, const char *msg);
Matrix *transpose_matrix(const Matrix *a);
Matrix *laplace_matrix(const Matrix *a);
Matrix *bin_matrix(const Matrix *a, const int nbin, const int *bin0, const int *bin1);
double  invert_matrix(Matrix *a);
void    svd_matrix(const Matrix *a, Matrix *u, Matrix *v, Matrix *w);
void    svd_solve(const Matrix *a, const Matrix *u, const Matrix *v, const Matrix *w, Matrix *b, Matrix *x);

Matrix *simple_image_read(const char *fname);
void    simple_image_write(const char *fname, Matrix *a);