#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double PowerSpec(double k);
double GrowthFactor(double astart, double aend);
double growth(double a);
inline double growth_int(double);
double qromb(double (*func)(double), double a, double b);
double sigma2_int(double k);
inline double TopHatSigma2(double R);
inline double PowerSpec_Efstathiou(double k);
inline double PowerSpec_EH(double k);
inline double PowerSpec_BBKS(double k);
inline double PowerSpec_CMBFAST(double k);

int initialize_powerspectrum(int Spectrum);

double tk_eh(double k);
double transfunc_cmbfast(double k);
inline double Hubble_a(double a);

void readCMB_and_do_spline();

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);
void splie3(double x1a[], double x2a[], double x3a[], double ***ya, int l, int m, int n, double ***y2a);
void splin3(double x1a[], double x2a[], double x3a[], double ***ya, double ***y2a, int l, int m, int n, double x1, double x2, double x3, double *y);

double Omega_de(double a);
double coeff1(double a);
double coeff2(double a);//u"+coeff1(a)u'+coeff2(a)u=0, u=D/a, '=d/dlna
double RungeKutta(double a_in,double a); //calculate linear density growth eq.
double growth_for_any_w(double a);
void growth_spline();

//non-linear matter Pk
double P_nonlinear(double z, double k);
void set_halofit_param(double z, double *param);
double solver(double z);
double get_delta_k(double k);
double sigma2_gauss_int(double lnk, void *params);
double sigma2_gauss(double lnR, void *params);
double dsig_dR_int(double lnk, void *params);
double dsig_dR(double R);
double d2sig_dR2_int(double lnk, void *params);
double d2sig_dR2(double R);
double neff(double R);
double C_halofit(double R);

void stack_data_and_spline_Pk();

double dsig_dR_fast(double R);
double d2sig_dR2_fast(double R);
double sigma2_gauss_fast(double lnR);

double TopHatSigma2_NL(double R);
double sigma2_NL_int(double k, double z);
