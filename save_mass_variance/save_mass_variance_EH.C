#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "proto.h"
#include "allvars.h"
#include "nrutil.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

static char TransTable[2000];
static double xp[1000],yp[1000],yp2[1000];
static double tab_z[NPOINTS],tab_chi[NPOINTS],err_chi[NPOINTS];
static double scale_f[NPOINTS],GF[NPOINTS],GF2[NPOINTS];

static double tab_R[NPOINTS],tab_dsdr[NPOINTS],tab_ds2dr2[NPOINTS],err_dsdr[NPOINTS],err_ds2dr2[NPOINTS];
static double tab_sig2[NPOINTS],err_sig2[NPOINTS];

static int WhichSpectrum, np, WhichWindow, OPT, OPT_fit, WhichGF;

static double bling;

static double r_tophat,Rsmooth,Delta_c,fdelta_c,Const_MF,st_norm;

static double AA,BB,CC;
static double B1,B2,B3,B4,B5;
static double nu, sigma, Omega_z;

static double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
static double G, Hubble;

static double Norm, InitTime, calc_z, a_exp;
static double Dplus; /* growth factor */

static double thre;

static double zhalo;

int main(int argc, char **argv){
  FILE *fin, *fout;
  int i,j;
  
  //input parameters======================================
  WhichSpectrum = 1;  // 1 Eisenstein-Hu, 2 Bond-Efstathiou, 3 CMBFAST table, 4 BBKS
  WhichGF = 2;        // 1 apploximation formula, 2 caluculate linear density growth eq.	
  OPT_fit = 2;        // 1 smith et al 03, 2 takahashi et al 12

  if(argc!=7){
    printf( "usage: %s output Om h0 ns sigma8 z\n", argv[0]);
    return 1;
  }

  char oname[256];
  sprintf(oname, "%s", argv[1]);
			
  //set cosmological model
  w  = -1.;	
  Omega = atof(argv[2]);
  OmegaLambda = 1.-Omega;
  HubbleParam = atof(argv[3]);
  ns = atof(argv[4]);
  Sigma8 = atof(argv[5]);
				
  zhalo = atof(argv[6]);
	
  if(WhichGF==2){growth_spline();}
	
  i=initialize_powerspectrum(WhichSpectrum);	
	
  stack_data_and_spline_Pk();

  fprintf(stdout, "Omega_m = %g, sigma8 (here) = %g\n", Omega*pow(1+zhalo, 3)/(Omega*pow(1+zhalo,3)+OmegaLambda), Sigma8 * GrowthFactor(1.0, 1./(1+zhalo)));

  int Nrbin = 1000;
  double rmin = 0.0;
  double rmax = 100.0;
  double dr = (rmax-rmin)/Nrbin;
      		
  fout = fopen(oname, "w");
  if(fout == NULL){
    printf("you can not make %s\n", oname);
    exit(1);
  }


  for(i=0;i<Nrbin;i++){
    double rbin = rmin+(i+0.5)*dr;
    double sig2 = TopHatSigma2_NL(rbin);
    fprintf(fout, "%e %e\n", rbin, sig2);      
  }
      
  fclose(fout);

  return 0;
}

 //==================================================
int initialize_powerspectrum(int Spectrum)
{
  double res;
  int i;

  for(i=0;i<1000;i++)
    xp[i]=yp[i]=yp2[i]=0;

  if(WhichSpectrum==3){
    //fprintf(stdout,"initialising...\n");
    readCMB_and_do_spline();
  }

  a_exp=1/(1+calc_z);

  AA=6.4/Gamma; 
  BB=3.0/Gamma; 
  CC=1.7/Gamma;  
  nu=1.13;


  B1=2.34/Gamma; 
  B2=3.89/Gamma; 
  B3=16.1/Gamma; 
  B4=5.46/Gamma; 
  B5=6.71/Gamma; 


  Norm = 1.0;
  //Sigma8 = sqrt(TopHatSigma2(8.));

  res = TopHatSigma2(8.);
  Norm=Sigma8*Sigma8/res;
  //fprintf(stdout,"Sigma8 = %g \n",Sigma8);

  Dplus= GrowthFactor(a_exp, 1.0);

  return i;
}

 //==================================================
double PowerSpec(double kmag)
{

  switch(WhichSpectrum)
    {
    case 1:
      return PowerSpec_EH(kmag);
      break;
    case 2:
      return PowerSpec_Efstathiou(kmag);
      break;
    case 3:
      return PowerSpec_CMBFAST(kmag);
      break;
    case 4:
      return PowerSpec_BBKS(kmag);
      break;
    default:
      fprintf(stdout,"Not supported\n");  
    }

}

//==================================================
inline double PowerSpec_Efstathiou(double k)
{
  return Norm*pow(k,ns) / pow(1+pow(AA*k+pow(BB*k,1.5)+CC*CC*k*k,nu),2/nu);
}

//==================================================
inline double PowerSpec_BBKS(double k)
{
  return Norm*pow(k,ns) * pow(log(1.0+B1*k)/(B1*k),2)/ 
    pow(1+ B2*k + B3*B3*k*k + pow(B4*k,3) + pow(B5*k,4),0.5);
}

//==================================================
double   tk_eh(double k)  
{
  double   q,theta,ommh2,a,s,gamma,L0,C0;
  double   tmp;
  double   omegam, ombh2, hubble;

  /* other input parameters */
  hubble= HubbleParam;

  omegam= Omega;
  ombh2=  OmegaBaryon*HubbleParam*HubbleParam;

  //k*= 1000.0;    /* convert to h/Mpc */
  /*k*= HubbleParam;*/  /* convert to 1/Mpc */

  theta = 2.728/2.7;
  ommh2 = omegam*hubble*hubble;
  s     = 44.5*log(9.83/ommh2)/sqrt( 1.+10.*exp(0.75*log(ombh2)) )*hubble;      
  a     = 1.-0.328*log(431.*ommh2)*ombh2/ommh2
    +0.380*log(22.3*ommh2)*(ombh2/ommh2)*(ombh2/ommh2);
  gamma = a+(1.-a)/(1.+exp(4*log(0.43*k*s)) );
  gamma*= omegam*hubble;
  q     = k*theta*theta/gamma;
  L0    = log( 2.*exp(1.)+1.8*q );
  C0    = 14.2 + 731./(1.+62.5*q);
  tmp   = L0/(L0+C0*q*q);

  return tmp;
}

//==================================================
inline double PowerSpec_EH(double k)
{
  return Norm*pow(k,ns)*pow( tk_eh(k), 2);
}

//==================================================
inline double PowerSpec_CMBFAST(double k)
{
  //return Norm*pow(k,ns)*pow(transfunc_cmbfast(k), 2);
  return Norm *transfunc_cmbfast(k);
}


//==================================================
//==================================================
double transfunc_cmbfast(double k)
{
  int i;
  double lk;
  double pow_index;

  //k *= 1000.0; /* convert to h/Mpc */

  lk=log10(k);

  if(lk < xp[0]){
    double dummy = (yp[2]-yp[1])/(xp[2]-xp[1])*(lk-xp[1])+yp[1];
    return pow(10.,dummy);
  } /* usually should not be executed */
  if(lk > xp[np-1]){
    double dummy = (yp[np-2]-yp[np-1])/(xp[np-2]-xp[np-1])*(lk-xp[np-1])+yp[np-1];
    return pow(10.,dummy);
    //return pow(10.,yp[np-1])*pow(k/pow(10.,xp[np-1]),-2.);
    //return pow(10.,yp[np-1]);
  }

  splint(xp-1, yp-1, yp2-1, np, lk, &pow_index);

  return pow(10.0,pow_index);
}

//==================================================
void readCMB_and_do_spline()
{
  int i,iline;
  double yp1,ypn;
  char tfunc_table[100];
  int  errorFlag=0;
  FILE *fd;

  sprintf(tfunc_table,"%s",TransTable);

  fprintf(stdout,"Reading %s .\n",tfunc_table);
  iline=0;
  double dummy1,dummy2,dummy3,dummy4,dummy5;
  fd=fopen(tfunc_table,"r");
  if(fd != NULL){
    while(!feof(fd)){
      fscanf(fd,"%lf %lf\n",&xp[iline],&yp[iline]);
      //fprintf(stdout,"%g %g \n",xp[iline],yp[iline]);
      iline++; 
    }
    fclose(fd);
  }
  else{
    fprintf(stdout,"transfer function file %s not found.\n",tfunc_table);
    exit(1);
  }

  fprintf(stdout,"read in %d data points \n",iline);

  np=iline;

  for(i=0;i<iline;i++){
    xp[i]=log10(xp[i]);
    yp[i]=log10(yp[i]);
  }

  yp1 = 1.e31;
  ypn = 1.e31;

  spline(xp-1, yp-1, iline, yp1, ypn, yp2-1);

  //for(i=0;i<iline;i++)
  //fprintf(stdout,"%g %g \n",xp[i],yp[i]);

}

//==================================================
inline double TopHatSigma2(double R)
{
  r_tophat= R;

  return qromb(sigma2_int, 0, 1000/R);
}

//==================================================
double sigma2_int(double k)
{
  double kr,kr3,kr2,wf,x;

  kr=r_tophat*k;
  kr2=kr*kr;
  kr3=kr2*kr;

  if(kr<1e-8) return 0;

  wf=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
  x=4*PI*k*k*wf*wf*PowerSpec(k)/pow(2*PI,3.);

  return x;
}


//==================================================
double GrowthFactor(double astart, double aend)
{
  return growth(aend)/growth(astart);
}


inline double Hubble_a(double a)
{
  double res;
	 
  res= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));

  return res;
}

double Omega_de(double a){
  double res;
  res = OmegaLambda/(Omega*pow(a,3*w)+OmegaLambda);
  return res;
}	

double coeff1(double a){
  double res;
  res = 0.5*(5-3*w*Omega_de(a));
  return res;
}

double coeff2(double a){
  double res;
  res = 1.5*(1-w)*Omega_de(a);
  return res;
}

double growth(double a)
{
  double hubble_a;

  if(w != -1){
    hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));
  }else{
    hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda);
  }

  switch(WhichGF){
  case 1:
    return hubble_a*qromb(growth_int, 0, a);
    break;
  case 2:
    return growth_for_any_w(a);
    break;
  default:
    fprintf(stdout,"Not supported\n"); 
  }
}

inline double growth_int(double a)
{
  return pow(a / (Omega + (1-Omega-OmegaLambda)*a + OmegaLambda*a*a*a), 1.5);
}

double RungeKutta(double a_in,double a){
  // u=D/a ,initial condition---> du/dlna=0,u=1
  int i,j;
  double h=(log(a)-log(a_in))/10000;
  double x=log(a_in);
  double u=1;
  double dudlna=0;
  double k0[2],k1[2],k2[2],k3[2];

  if(a_in==0){
    printf("you cannot solve calculate linear density growth eq.");
  }if(a == a_in){
    u=1;
  }else{
    for(i=0;i<10000;i++){

      k0[0]=h*dudlna;
      k0[1]=h*(-coeff1(exp(x))*dudlna-coeff2(exp(x))*(u));

      k1[0]=h*(dudlna+k0[1]/2);
      k1[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k0[1]/2)-coeff2(exp(x+h/2))*(u+k0[0]/2));

      k2[0]=h*(dudlna+k1[1]/2);
      k2[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k1[1]/2)-coeff2(exp(x+h/2))*(u+k1[0]/2));

      k3[0]=h*(dudlna+k2[1]);
      k3[1]=h*(-coeff1(exp(x+h))*(dudlna+k2[1])-coeff2(exp(x+h))*(u+k2[0]));

      u = u + (k0[0]+2*k1[0]+2*k2[0]+k3[0])/6;
      dudlna = dudlna + (k0[1]+2*k1[1]+2*k2[1]+k3[1])/6;
      x = x+h;
    }
  }

  return a*u;
}

void growth_spline(){
  int i;
  double yp1,ypn;
  double da = (1.0-(1.0/(1.0+1088.2)))/NPOINTS;

  for(i=0;i<NPOINTS;i++){
    scale_f[i] = (double)(i+1)*da + 1.0/(1.0+1088.2);
    GF[i] = RungeKutta(1.0/(1.0+1088.2),scale_f[i]);
  }

  for(i=0;i<NPOINTS;i++){
    scale_f[i] = log10(scale_f[i]);
    GF[i] = log10(GF[i]);
  }

  yp1 = 1.e31;
  ypn = 1.e31;

  spline(scale_f-1, GF-1, NPOINTS, yp1, ypn, GF2-1);

}

double growth_for_any_w(double a){
  double la;
  double pow_index;

  la=log10(a);

  if(la < scale_f[0])
    return a; /* usually should not be executed */
  if(la > scale_f[NPOINTS-1]){
    double dummy = (GF[NPOINTS-2]-GF[NPOINTS-1])/(GF[NPOINTS-2]-GF[NPOINTS-1])*(la-scale_f[NPOINTS-1])+GF[NPOINTS-1];
    return pow(10.,dummy);
    //return pow(10.,yp[np-1]);
  }

  splint(scale_f-1, GF-1, GF2-1, NPOINTS, la, &pow_index);
  return pow(10.0,pow_index);
}

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a){
  int j;
  for(j=1;j<=m;j++){
    spline(x2a,ya[j],n,1.e31,1.e31,y2a[j]);
  }
}
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y){
  int j;
  double *ytmp,*yytmp;

  ytmp=dvector(1,n);
  yytmp=dvector(1,n);

  for(j=1;j<=m;j++){
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  }
  spline(x1a,yytmp,m,1.e31,1.e31,ytmp);
  splint(x1a,yytmp,ytmp,m,x1,y);
  free_dvector(yytmp,1,n);
  free_dvector(ytmp,1,n);
}

void splie3(double x1a[], double x2a[], double x3a[], double ***ya, int l, int m, int n, double ***y2a){
  int j;
  for(j=1;j<=l;j++){
    splie2(x2a,x3a,ya[j],m,n,y2a[j]);
  }
}
void splin3(double x1a[], double x2a[], double x3a[], double ***ya, double ***y2a, int l, int m, int n, double x1, double x2, double x3, double *y){
  int j;
  double *ytmp,*yytmp;

  ytmp=dvector(1,l);
  yytmp=dvector(1,l);

  for(j=1;j<=l;j++){
    splin2(x2a,x3a,ya[j],y2a[j],m,n,x2,x3,&yytmp[j]);
  }
  spline(x1a,yytmp,l,1.e31,1.e31,ytmp);
  splint(x1a,yytmp,ytmp,l,x1,y);
  free_dvector(yytmp,1,l);
  free_dvector(ytmp,1,l);
}

double P_nonlinear(double z, double k){

  double ksig ,n_eff, C_halo;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  double param[3];
  set_halofit_param(z, param);
  ksig = param[0];
  n_eff = param[1];
  C_halo = param[2];
  double a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n;
  if(OPT_fit==1){
    a_n = pow(10., 1.4861 + 1.8369*n_eff + 1.6762*n_eff*n_eff + 0.7940*n_eff*n_eff*n_eff + 0.1670*n_eff*n_eff*n_eff*n_eff -0.6206*C_halo);
    b_n = pow(10.,0.9463 + 0.9466*n_eff + 0.3084*n_eff*n_eff - 0.9400*C_halo);
    c_n = pow(10.,-0.2807 + 0.6669*n_eff + 0.3214*n_eff*n_eff - 0.0793*C_halo);
    gamma_n = 0.8649 + 0.2989*n_eff + 0.1631*C_halo;
    alpha_n = 1.3884 + 0.3700*n_eff - 0.1452*n_eff*n_eff;
    beta_n = 0.8291 + 0.9854*n_eff + 0.3401*n_eff*n_eff;
    mu_n = pow(10.,-3.5442 + 0.1908*n_eff);
    nu_n = pow(10., 0.9589 + 1.2857*n_eff);
  }
  if(OPT_fit==2){
    double om_w = OmegaLambda*pow(1+z,3*w)/(Omega+OmegaLambda*pow(1+z,3*w)); 
    a_n = pow(10.,1.5222 + 2.8553*n_eff + 2.3706*n_eff*n_eff + 0.9903*n_eff*n_eff*n_eff + 0.2250*n_eff*n_eff*n_eff*n_eff -0.6038*C_halo +0.1749*om_w*(1.+w));
    b_n = pow(10.,-0.5642 + 0.5864*n_eff + 0.5716*n_eff*n_eff - 1.5474*C_halo +0.2279*om_w*(1.+w));
    c_n = pow(10.,0.3698 + 2.0404*n_eff + 0.8161*n_eff*n_eff + 0.5869*C_halo);
    gamma_n = 0.1971 - 0.0843*n_eff + 0.8460*C_halo;
    alpha_n = fabs(6.0835 + 1.3373*n_eff - 0.1959*n_eff*n_eff - 5.5274*C_halo);
    beta_n = 2.0379 - 0.7354*n_eff + 0.3157*n_eff*n_eff + 1.2490*n_eff*n_eff*n_eff + 0.3980*n_eff*n_eff*n_eff*n_eff - 0.1682*C_halo;
    mu_n = 0;
    nu_n = pow(10., 5.2105 + 3.6902*n_eff);
  }
  //printf("%e %e %e %e %e %e %e\n",a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n);


  double delta_L,delta_H,delta_Q;
  //only support flat universe
  double omz;
  omz = Omega/(Omega+OmegaLambda*pow(1+z,3*w));

  double f1 = pow(omz,-0.0307);
  double f2 = pow(omz,-0.0585);
  double f3 = pow(omz,0.0743);

  double fac =4*PI*k*k*k/(2*PI)/(2*PI)/(2*PI);
  delta_L = fac*PowerSpec(k)/inv_D;

  double y = k/ksig;
  delta_Q = delta_L*(pow(1.+delta_L,beta_n)/(1.+alpha_n*delta_L))*exp(-y/4-y*y/8);

  delta_H = a_n*pow(y,3*f1)/(1+b_n*pow(y,f2)+pow(c_n*f3*y,3-gamma_n));
  delta_H = delta_H/(1+mu_n/y+nu_n/y/y);

  return (delta_Q+delta_H)/fac;
}

void set_halofit_param(double z, double *param){
  int i;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  thre = inv_D;
  //double init = pow(10.,-1.065+4.332e-1*(1.+z)-2.516e-2*pow(1.+z,2)+9.069e-4*pow(1.+z,3));
  //double R0=solver(z);
  //printf("z=%e\n",z);
  double R0=solver(z);
  //check the solver
  //double dummy=0;
  //printf("check_solver %e\n",sigma2_gauss(log(R0),&dummy));

  param[0] = 1./R0; //k_sig
  param[1] = neff(R0);
  param[2] = C_halofit(R0);

}

double solver(double z){
  int status; 
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0.;
  double x_lo = -5.0, x_hi = 5.0; 
  gsl_function F; 
  double params = z;
  F.function = &sigma2_gauss; 
  F.params = &params; 
  T = gsl_root_fsolver_brent; 
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, x_lo, x_hi); 

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s); 
    r = gsl_root_fsolver_root(s); 
    x_lo = gsl_root_fsolver_x_lower(s); 
    x_hi = gsl_root_fsolver_x_upper(s); 
    status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7); 
    //if (status == GSL_SUCCESS) printf ("Converged:\n");
    //printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);

  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  double R0 = exp(r);

  //printf("Rsig = %e\n",R0);
  return R0;
}

double get_delta_k(double k){
  double lk;
  double pow_index;
  double fac = k*k*k*4*PI/(2*PI)/(2*PI)/(2*PI);

  return fac*PowerSpec(k);

}

double sigma2_gauss_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);

  res = get_delta_k(k)*exp(-k*k*R*R);
  return res;
}

double sigma2_gauss(double lnR, void *params){
  double result,abserr;
  double params_int=exp(lnR);
  double z = *(double *)params;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&sigma2_gauss_int;
  //F.params=&params_int;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  /*int Nint = 5000;

    double dlnk = 10.0/Nint;
    result = 0;
    for(int i=0;i<Nint;i++){
    result += dlnk*sigma2_gauss_int(i*dlnk-5.0,&params_int);
    }*/
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));

  return sigma2_gauss_fast(params_int)-inv_D;
}

double dsig_dR_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);

  res = get_delta_k(k)*exp(-k*k*R*R)*(-2*k*k*R);
  return res;
}

double dsig_dR(double R){
  double result,abserr;
  double params=R;
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&dsig_dR_int;
  //F.params=&params;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  result = 0;
  for(int i=0;i<Nint;i++){
    result += dlnk*dsig_dR_int(i*dlnk-10.0,&params);
  }
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));

  return result;
}
double d2sig_dR2_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);
  double y = k*R;

  res = get_delta_k(k)*exp(-y*y)*(y*y-y*y*y*y);

  return res;
}
double d2sig_dR2(double R){
  double result,abserr;
  double params=R;
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&d2sig_dR2_int;
  //F.params=&params;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  result = 0;
  for(int i=0;i<Nint;i++){
    result += dlnk*d2sig_dR2_int(i*dlnk-10.0,&params);
  }
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));
  //if(result < 0){printf("R=%e d^2sig/dR^2 = %e\n",R,result);}

  return result;
}

double neff(double R){
  double sig2 = thre;
  double res = -R*dsig_dR_fast(R)/sig2;

  return res -3.0;
}

double C_halofit(double R){
  double sig2 = thre;
  double n_eff = neff(R);

  double res = (3.+n_eff)*(3+n_eff) +4./thre*d2sig_dR2_fast(R); //cf .smith et al

  return res;

}

void stack_data_and_spline_Pk(){
  double logR = (3.-(-3.))/NPOINTS;
  int i;

  for(i=0;i<NPOINTS;i++){
    tab_R[i] = (double)(i)*logR -3.;
    tab_dsdr[i] = dsig_dR(pow(10.,tab_R[i]));
    tab_ds2dr2[i] = d2sig_dR2(pow(10.,tab_R[i]));
    //printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
    tab_dsdr[i] = log10(-tab_dsdr[i]);
    tab_ds2dr2[i] = log10(10+tab_ds2dr2[i]);
    //printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
  }

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  double params_int;
  for(i=0;i<NPOINTS;i++){
    tab_sig2[i]=0.;
    params_int = pow(10.,tab_R[i]);
    for(int j=0;j<Nint;j++){
      tab_sig2[i] += dlnk*sigma2_gauss_int(j*dlnk-10.0,&params_int);
    }
    tab_sig2[i] = log10(tab_sig2[i]);
  }

  double yp1 = 1.e31;
  double ypn = 1.e31;

  spline(tab_R-1, tab_sig2-1, NPOINTS, yp1, ypn, err_sig2-1);
  printf("spline sigma2 ... \n");
  spline(tab_R-1, tab_dsdr-1, NPOINTS, yp1, ypn, err_dsdr-1);
  printf("spline dsig/dr ... \n");
  spline(tab_R-1, tab_ds2dr2-1, NPOINTS, yp1, ypn, err_ds2dr2-1);
  printf("spline d2sig/dr2 ... \n");
}

double dsig_dR_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    return dsig_dR(R);
  }else{
    splint(tab_R-1, tab_dsdr-1, err_dsdr-1, NPOINTS, lR, &pow_index);
    return -pow(10.,pow_index);
  }
}

double d2sig_dR2_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    return d2sig_dR2(R);
  }else{
    splint(tab_R-1, tab_ds2dr2-1, err_ds2dr2-1, NPOINTS, lR, &pow_index);
    return pow(10.,pow_index)-10.;
  }
}

double sigma2_gauss_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    int Nint = 5000;
    double dlnk = 20.0/Nint;
    double result=0.;

    for(int i=0;i<Nint;i++){
      result += dlnk*sigma2_gauss_int(i*dlnk-10.0,&R);
    }
    return result;
  }else{
    splint(tab_R-1, tab_sig2-1, err_sig2-1, NPOINTS, lR, &pow_index);
    return pow(10.,pow_index);
  }
}

//==================================================
double TopHatSigma2_NL(double R)
{
  r_tophat= R;

  double res,abserr;
  double x_lo = 0.0;
  double x_hi = 1000.0/R;

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(sigma2_NL_int(x_hi, zhalo)+sigma2_NL_int(x_lo, zhalo));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*sigma2_NL_int(x_lo+(2*j-1)*h[i-1], zhalo);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
  
}

//==================================================
double sigma2_NL_int(double k, double z)
{
  double kr,kr3,kr2,wf,x;

  kr=r_tophat*k;
  kr2=kr*kr;
  kr3=kr2*kr;

  if(kr<1e-8) return 0;

  wf=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
  x=4*PI*k*k*wf*wf*P_nonlinear(z, k)/pow(2*PI,3.);

  return x;
}
