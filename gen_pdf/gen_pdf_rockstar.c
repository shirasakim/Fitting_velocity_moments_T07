#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h> 

const double PI=3.14159265358979;
const double G=4.299e-9; //km^2 Mpc/M_sun/s^2
const double C=2.99792458e5; //light of speed  km/s
const double TINY=1.0e-48;

static double Lbox; // in Mpc/h

struct halo_data{
  int ID;
  int DescID;
  float Mvir;
  float Vmax;
  float Vrms; 
  float Rvir;
  float Rs;
  int Np;
  //
  float X;
  float Y;
  float Z;
  float VX;
  float VY;
  //
  float VZ;
  float JX;
  float JY;
  float JZ;
  float Spin;
  //
  float rs_klypin;
  float Mvir_all;
  float M200b;
  float M200c;
  float M500c;
  //
  float M2500c;
  float Xoff;
  float Voff;
  float spin_bullock;
  float b_to_a;
  //
  float c_to_a;
  float Ax;
  float Ay; 
  float Az;
  float b_to_a_500c;
  //
  float c_to_a_500c;
  float Ax_500c;
  float Ay_500c;
  float Az_500c;
  float TtoU;
  //
  float M_pe_Behroozi;
  float M_pe_Diemer;
  float Halfmass_Radius;
  int PID;
};

struct myhalo_small
{
  float pos[3], vel[3];
  float mass;
  int central;
};

static int *D_cen1, *D_cen2;

static struct myhalo_small *halos;
static int Nhalo;

void read_halos(char *fname);

int main(int argc, char *argv[])
{
  int i,j;
	
  if(argc!=4){
    printf("usage: %s input outputfile Lbox[Mpc/h]\n", argv[0]);
    return 1;
  }
  
  char prof_name[256];
  sprintf(prof_name, "%s", argv[1]);

  Lbox = atof(argv[3]);

  char oname[256];

  read_halos(prof_name);

  int Nrbin = 200;  
  double rmin = 0.0;
  double rmax = 40.0;
  double dr = (rmax-rmin)/Nrbin;  

  int Ndiv = 50;
  double Ldiv = Lbox/Ndiv;

  FILE *fout;

  int rsearch=(int)(rmax/Ldiv)+1;

  double vmin = -2500.0;
  double vmax = +2500.0;
  int Nvbin = 100;
  double delta = (vmax-vmin)/Nvbin;

  double *tab_vr, *tab_vt, *norm, *tab_vlos;
  
  tab_vr = malloc(sizeof(double) * Nvbin * Nrbin );
  tab_vt = malloc(sizeof(double) * Nvbin * Nrbin );
  tab_vlos = malloc(sizeof(double) * Nvbin * Nrbin );
  norm = malloc(sizeof(double) * Nrbin);

  if(tab_vr == NULL || tab_vt == NULL || norm == NULL || tab_vlos == NULL){
    fprintf(stderr, "cannot malloc tab_pdf\n");
    exit(1);
  }
 
  double M1_min, M1_max;
  double M2_min, M2_max;
  double logMmin = 11.0;
  double logMmax = 15.0;
  int iMbin, jMbin, nMbin = 8;
  double dlogM = (logMmax-logMmin)/nMbin;

  for(iMbin=0;iMbin<nMbin;iMbin++){

    M1_min = pow(10., logMmin + (iMbin+0.0)*dlogM);
    M1_max = pow(10., logMmin + (iMbin+1.0)*dlogM);

    int Ncen1=0;
    for(i=0;i<Nhalo;i++){
      if(halos[i].central == -1 && halos[i].mass >= M1_min && halos[i].mass <= M1_max){
	Ncen1 +=1 ;
      }
    }
    fprintf(stdout, "sample1 : %d halos with their mass range of %5.4e - %5.4e\n", Ncen1, M1_min, M1_max);
  
    D_cen1 = malloc(sizeof(int)*Ncen1);

    if(D_cen1 == NULL){
      fprintf(stderr, "failed to malloc D_cen1\n");
      exit(1);
    }

    int icen1=0;
    for(i=0;i<Nhalo;i++){
      if(halos[i].central == -1 && halos[i].mass >= M1_min && halos[i].mass <= M1_max){
	D_cen1[icen1] = i;
	icen1+=1;
      }
    }

    int *Nsample1_in_subv;
    Nsample1_in_subv = malloc(sizeof(int)*Ndiv*Ndiv*Ndiv);

    if(Nsample1_in_subv == NULL){
      fprintf(stderr, "failed to malloc Nsample1_in_subv\n");
      exit(1);
    }

    for(i=0;i<Ndiv*Ndiv*Ndiv;i++){
      Nsample1_in_subv[i] =0;
    }
    for(i=0;i<Ncen1;i++){
      int ix = (int)(halos[D_cen1[i]].pos[0]/Ldiv);
      int iy = (int)(halos[D_cen1[i]].pos[1]/Ldiv);
      int iz = (int)(halos[D_cen1[i]].pos[2]/Ldiv);

      if(ix >= 0 && ix < Ndiv){
	if(iy >= 0 && iy < Ndiv){
	  if(iz >=0 && iz < Ndiv) Nsample1_in_subv[iz+Ndiv*(iy+Ndiv*ix)] += 1;
	}
      }    

    }

    int NPmax1 = -1;
    for(i=0;i<Ndiv*Ndiv*Ndiv;i++){
      if(Nsample1_in_subv[i] > NPmax1){
	NPmax1 = Nsample1_in_subv[i];
      }
    }

    if(NPmax1 <= 0){
      fprintf(stderr, "found negative value of NPmax1. Something is wrong.\n");
      exit(1);
    }

    int *id_sample1_in_subv;
    id_sample1_in_subv = malloc(sizeof(int)*Ndiv*Ndiv*Ndiv*NPmax1);
    if(id_sample1_in_subv == NULL){
      fprintf(stderr, "faield to malloc id_sample1_in_subv\n");
      exit(1);
    }

    for(i=0;i<Ndiv*Ndiv*Ndiv;i++){
      Nsample1_in_subv[i] =0;
    }

    for(i=0;i<Ncen1;i++){
      int ix = (int)(halos[D_cen1[i]].pos[0]/Ldiv);
      int iy = (int)(halos[D_cen1[i]].pos[1]/Ldiv);    
      int iz = (int)(halos[D_cen1[i]].pos[2]/Ldiv);

      if(ix >= 0 && ix < Ndiv){
	if(iy >= 0 && iy < Ndiv){
	  if(iz >=0 && iz < Ndiv){
	    id_sample1_in_subv[Nsample1_in_subv[iz+Ndiv*(iy+Ndiv*ix)] + NPmax1 * (iz+Ndiv*(iy+Ndiv*ix))] = i;
	    Nsample1_in_subv[iz+Ndiv*(iy+Ndiv*ix)] += 1;
	  }
	}
      }    
    }

    int jMstart = (iMbin <= 2) ? 2 : iMbin;
    for(jMbin=jMstart;jMbin<nMbin;jMbin++){

      M2_min = pow(10., logMmin + (jMbin+0.0)*dlogM);
      M2_max = pow(10., logMmin + (jMbin+1.0)*dlogM);

      int Ncen2=0;
      for(i=0;i<Nhalo;i++){
	if(halos[i].central == -1 && halos[i].mass >= M2_min && halos[i].mass <= M2_max){
	  Ncen2 +=1 ;
	}
      }
      fprintf(stdout, "sample2 : %d halos with their mass range of %5.4e - %5.4e\n", Ncen2, M2_min, M2_max);
  
      D_cen2 = malloc(sizeof(int)*Ncen2);

      if(D_cen2 == NULL){
	fprintf(stderr, "failed to malloc D_cen2\n");
	exit(1);
      }

      int icen2=0;
      for(i=0;i<Nhalo;i++){
	if(halos[i].central == -1 && halos[i].mass >= M2_min && halos[i].mass <= M2_max){
	  D_cen2[icen2] = i;
	  icen2+=1;
	}
      }

      for(i=0;i<Nrbin;i++){
	norm[i] = 0;
      }

      for(i=0;i<Nrbin*Nvbin;i++){
	tab_vr[i] = 0;
	tab_vt[i] = 0;
	tab_vlos[i] =0;
      }

      long long int Npairs = 0;
 
      fprintf(stdout, "start pair counting.\n");
       
      int ic;
      for(ic=0;ic<Ncen2;ic++){
      
	int icx = (int)(halos[D_cen2[ic]].pos[0]/Ldiv);
	int icy = (int)(halos[D_cen2[ic]].pos[1]/Ldiv);
	int icz = (int)(halos[D_cen2[ic]].pos[2]/Ldiv);

	//if(ic > 0 && ic%(1000000) == 0)fprintf(stdout, ".");fflush(stdout);	
	
	int ix, iy, iz;
	for(ix=-rsearch;ix<=rsearch;ix++){
	  for(iy=-rsearch;iy<=rsearch;iy++){
	    for(iz=-rsearch;iz<=rsearch;iz++){
	    
	      int itx = icx + ix;
	      int ity = icy + iy;
	      int itz = icz + iz;

	      if(itx >= 0 && itx < Ndiv){
		if(ity >= 0 && ity < Ndiv){
		  if(itz >= 0 && itz < Ndiv){
		  
		    int Ncen1_here = Nsample1_in_subv[itz+Ndiv*(ity+Ndiv*itx)];
		    int it;
		    for(it=0;it<Ncen1_here;it++){
		      int id = id_sample1_in_subv[it+NPmax1*(itz+Ndiv*(ity+Ndiv*itx))];
		      
		      double dpos[3];
		      for(j=0;j<3;j++){dpos[j] = halos[D_cen1[id]].pos[j] - halos[D_cen2[ic]].pos[j];}
		      
		      double dvel[3];
		      for(j=0;j<3;j++){dvel[j] = halos[D_cen1[id]].vel[j] - halos[D_cen2[ic]].vel[j];}

		      double pos_mid[3];
		      for(j=0;j<3;j++){pos_mid[j] = 0.5 * (halos[D_cen1[id]].pos[j] + halos[D_cen2[ic]].pos[j]); }

		      double radius = sqrt(dpos[0]*dpos[0]+dpos[1]*dpos[1]+dpos[2]*dpos[2]);
		      int irbin = (int)((radius - rmin)/dr);
		      		      
		      if(irbin >= 0 && irbin < Nrbin){
					       		       
			double v_r = (dvel[0] * dpos[0] + dvel[1] * dpos[1] + dvel[2] * dpos[2])/radius;
						
			double cost = dpos[2]/radius;
			double sint = sqrt(dpos[0]*dpos[0]+dpos[1]*dpos[1])/radius;
			double cosphi = dpos[0]/sqrt(dpos[0]*dpos[0]+dpos[1]*dpos[1]);
			double sinphi = dpos[1]/sqrt(dpos[0]*dpos[0]+dpos[1]*dpos[1]);

			if(sqrt(dpos[0]*dpos[0]+dpos[1]*dpos[1]) < 1e-10){
			  cosphi = 1.0;
			  sinphi = 0.0;
			}

			double v_t = dvel[0] * cost * cosphi + dvel[1] * cost * sinphi - dvel[2] * sint;

			double vlos = dvel[2];

			int ivbin;
			ivbin = (int)((v_t - vmin) / (delta));
			if(ivbin >= 0 && ivbin < Nvbin) tab_vt[ivbin + (irbin)*Nvbin] += 1;

			ivbin = (int)((v_r - vmin) / (delta));
			if(ivbin >= 0 && ivbin < Nvbin) tab_vr[ivbin + (irbin)*Nvbin] += 1;

			ivbin = (int)((vlos - vmin) / (delta));
			if(ivbin >= 0 && ivbin < Nvbin) tab_vlos[ivbin + (irbin)*Nvbin] += 1;
			
			norm[irbin] += 1.0;
		      }
		     
		      Npairs +=1; 					
		    }

		  }
		}
	      }
	    }
	  }
	}    
      }    

      fprintf(stdout, "found %lld pairs.\n", Npairs);

      sprintf(oname, "%s_Mbin%d_%d.dat", argv[2], iMbin, jMbin);

      fout = fopen(oname, "w");
      if(fout == NULL){
	fprintf(stderr, "cannot make %s\n", oname);
	exit(1);
      }
      fprintf(fout, "# %d %e %e\n", Ncen1, M1_min, M1_max);
      fprintf(fout, "# %d %e %e\n", Ncen2, M2_min, M2_max);
      
      for(i=0;i<Nrbin;i++){	

	fprintf(fout, "# %e %e %e %e\n", 
		(rmin) + (i+0)*dr, (rmin) + (i+1)*dr, norm[i], delta);

	for(j=0;j<Nvbin;j++){
	  fprintf(fout, "%e %e %e %e\n",  vmin + (j+0.5)*delta, 
		  tab_vt[j+i*Nvbin], tab_vr[j+i*Nvbin], tab_vlos[j+i*Nvbin]);
	}
      }
      fclose(fout);
      
      free(D_cen2);
    }
  
    if(iMbin==nMbin-1)fprintf(stdout, "DONE\n");
    
    free(D_cen1);
    free(Nsample1_in_subv);
    free(id_sample1_in_subv);
  }

  free(norm);
  free(tab_vr);
  free(tab_vt);
  free(tab_vlos);

  free(halos);

  return 0;
}

void read_halos(char *fname){

  int i;
  FILE *fp;

  char readline[2048];
  int row_header = 16;

  fp = fopen(fname, "r");
  if(fp == NULL){

    fprintf(stderr, "can not find %s\n", fname);
    exit(1);

  }

  for(i=0; i<row_header; i++) fgets(readline, 2048, fp);

  Nhalo = 0;
  while(fgets(readline, 2048, fp)) Nhalo += 1;

  fclose(fp);

  if(Nhalo==0){
    fprintf(stderr, "there are no halos in %s.\n", fname);
    exit(1);
  }

  halos = malloc(sizeof(struct myhalo_small)*Nhalo);
  if(halos==NULL){
    fprintf(stderr, "failed to malloc sizeof(struct myhalo) x %d\n", Nhalo);
    exit(1);
  }

  fp = fopen(fname, "r");
  for(i=0; i<row_header; i++) fgets(readline, 2048, fp);

  for(i=0;i<Nhalo;i++){

    if(i > 0 && i % 5000000 == 0) fprintf(stdout, "#"); fflush(stdout);

    struct halo_data h;

    float rvmax;
    fscanf(fp, "%d %d %e %e %e %e %e %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
	   &h.ID, &h.DescID,
	   &h.Mvir, &h.Vmax, &h.Vrms, &h.Rvir, &h.Rs, 
	   &h.Np, 
	   &h.X, &h.Y, &h.Z,
	   &h.VX, &h.VY, &h.VZ, 
	   &h.JX, &h.JY, &h.JZ, 
	   &h.Spin, &h.rs_klypin, &h.Mvir_all, 
	   &h.M200b, &h.M200c, &h.M500c, &h.M2500c, 
	   &h.Xoff, &h.Voff, &h.spin_bullock, &h.b_to_a, &h.c_to_a,
	   &h.Ax, &h.Ay, &h.Az,
	   &h.b_to_a_500c, &h.c_to_a_500c,
	   &h.Ax_500c, &h.Ay_500c, &h.Az_500c,
	   &h.TtoU, &h.M_pe_Behroozi, &h.M_pe_Diemer, &h.Halfmass_Radius, &rvmax, &h.PID);    

    halos[i].pos[0] = h.X; halos[i].pos[1] = h.Y; halos[i].pos[2] = h.Z;
    halos[i].vel[0] = h.VX; halos[i].vel[1] = h.VY; halos[i].vel[2] = h.VZ;
    halos[i].central = h.PID;
    halos[i].mass = h.M200b;

  }
  fclose(fp);

  fprintf(stdout, " Done reading %s.\n", fname);
  
}
