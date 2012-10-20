#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_sf_bessel.h>

# define DIM 3
# define KMAX 20

int *conf[DIM];

/* next holds the nearest neighbours. If DIM=2 then -2,-1,1,2   *
 * denotes the nearest -t,-x,x,t dirn. 0 refers the site itself */
int *next[2*DIM+1];

double *corr,*avgcorr;

double MU, BETA;

int LX,LY,LT,VOL;

int SEED,RNDFLAG;

int START;

int EQUIL,SWEEP,MEAS;

double IB[2*KMAX+1];

double W[2*DIM+1],WT;

double sus, avgsus;
double rho, rhosus,avgrho,avgrhosus;
double wsus,avgwsus;
double spsus,avgspsus;

int main(argc,argv)
     int argc;
     char *argv[];
{
  int d,k1,k2,t;
  int readval;
  char string[100];
  char beta[5];
  char file1[20]="SUS";
  char file2[20]="RHO";
  char file3[20]="RHOSUS";
  char file4[20]="WSUS";
  char file5[20]="SPSUS";
  char file6[20]="CORR";
  FILE *fptr;

  extern void initrandom(int,int);
  extern void initneighbor(void);
  extern void initbessel(void);
  extern void initconf(void);
  extern double update(void);
  extern void checkconf(void);
  extern void printconf(void);
  extern void currsus(void);

  if(argc != 1)
    {
      printf("usage xy (needs a QUEUE file)\n");
      exit(1);
    }

  fptr = fopen("QUEUE","r");
  if(fptr == NULL)
    {
      printf("could not open QUEUE FILE to open\n");
      exit(1);
    }
  readval=fscanf(fptr,"%s %d\n",string,&LX);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&LY);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&LT);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %lf\n",string,&MU);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %lf\n",string,&BETA);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&START);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&EQUIL);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&SWEEP);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&MEAS);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d",string,&SEED);
  if(readval<0) printf("Eror\n");
  readval=fscanf(fptr,"%s %d",string,&RNDFLAG);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s",beta);
  if(readval<0) printf("Error\n");
  fclose(fptr);

  VOL = LX*LY*LT;
  strcat(file1,beta);
  strcat(file2,beta);
  strcat(file3,beta);
  strcat(file4,beta);
  strcat(file5,beta);
  strcat(file6,beta);
  

  for(d=0;d<DIM;d++)
    conf[d] = (int *)malloc(VOL*sizeof(int));

  for(d=0;d<=2*DIM;d++)
    next[d] = (int *)malloc(VOL*sizeof(int));

  corr = (double *)malloc(LT*sizeof(double));
  avgcorr = (double *)malloc(LT*sizeof(double));

  initrandom(SEED,RNDFLAG);
  initneighbor();
  initbessel();
  initconf();

  for(k1=0;k1<EQUIL;k1++)
    {
      sus = update();
      checkconf();
    }

  fptr = fopen(file1,"w");
  fclose(fptr);
  fptr = fopen(file2,"w");
  fclose(fptr);
  fptr = fopen(file3,"w");
  fclose(fptr);
  fptr = fopen(file4,"w");
  fclose(fptr);
  fptr = fopen(file5,"w");
  fclose(fptr);
  fptr = fopen(file6,"w");
  fclose(fptr);


  for(k1=0;k1<MEAS;k1++)
    {
      avgsus = 0;
      avgrho = 0;
      avgrhosus = 0;
      avgwsus = 0;
      avgspsus = 0;
      for(t=0;t<LT;t++) avgcorr[t] = 0.0;
      for(k2=0;k2<SWEEP;k2++)
	{
	  sus = update();
	  avgsus += sus;
          avgspsus += spsus;
	  currsus();
	  avgrho += rho;
	  avgrhosus += rhosus;
	  avgwsus += wsus;
          for(t=0;t<LT;t++) avgcorr[t] += corr[t];
	  checkconf();
	}
      avgsus /= (double)SWEEP;
      avgrho /= (double)SWEEP;
      avgrhosus /= (double)SWEEP;
      avgwsus /= (double)SWEEP;
      avgspsus /= (double)SWEEP;
      for(t=0;t<LT;t++) avgcorr[t] /= (double)SWEEP;
      fptr = fopen(file1,"a");
      fprintf(fptr,"%16.8e\n",avgsus);
      fclose(fptr);
      fptr = fopen(file2,"a");
      fprintf(fptr,"%16.8e\n",avgrho);
      fclose(fptr);
      fptr = fopen(file3,"a");
      fprintf(fptr,"%16.8e\n",avgrhosus);
      fclose(fptr);
      fptr = fopen(file4,"a");
      fprintf(fptr,"%16.8e\n",avgwsus);
      fclose(fptr);
      fptr = fopen(file5,"a");
      fprintf(fptr,"%16.8e\n",avgspsus);
      fclose(fptr);
      fptr = fopen(file6,"a");
      for(t=0;t<LT;t++) fprintf(fptr,"%d %16.8e\n",t,avgcorr[t]);
      fclose(fptr);

    }

  for(d=0;d<DIM;d++) free(conf[d]);
  for(d=0;d<=2*DIM;d++) free(next[d]);
  free(corr); free(avgcorr);

  return 0;
}


void initneighbor(void)
{
  int p,x,y,t;

  for(p=0;p<VOL;p++)
    {
      x = (p%LX);
      y = (p/LX)%LY;
      t = (p/(LX*LY))%LT;

      next[DIM][p] = p;

      next[DIM+1][p] = t*LX*LY + y*LX + ((x+1)%LX);
      next[DIM+2][p] = t*LX*LY + ((y+1)%LY)*LX + x;
      next[DIM+3][p] = ((t+1)%LT)*LX*LY + y*LX + x;

      next[DIM-1][p] = t*LX*LY + y*LX + ((x-1+LX)%LX);
      next[DIM-2][p] = t*LX*LY + ((y-1+LY)%LY)*LX + x;
      next[DIM-3][p] = ((t-1+LT)%LT)*LX*LY + y*LX + x;
    }
}

void initbessel(void)
{
  int k1;

  for(k1=0;k1<=KMAX;k1++) 
    {
      IB[KMAX+k1] = gsl_sf_bessel_In(k1,BETA);
      IB[KMAX-k1] = IB[KMAX+k1];
    }
}

void initconf(void)
{
  int p,d;

  for(p=0;p<VOL;p++)
    {
      for(d=0;d<DIM;d++) conf[d][p] = 0;
    }

}

void checkconf(void)
{
  int p,d,k;

  for(p=0;p<VOL;p++)
    {
      k=0;
      for(d=0;d<DIM;d++)
	{
	  k += conf[d][p];
	  k -= conf[d][next[DIM-(d+1)][p]];
	}

      if(k != 0)
	{
	  printf("checkconf error 1: %6d %6d\n",p,k);
	  for(d=0;d<DIM;d++)
	    {
	      printf(" %4d %4d %4d\n",d+1,
		     conf[d][p],conf[d][next[DIM-(d+1)][p]]);
	    }
	  exit(1);
	}
    }
}

double update(void)
{
  int p,p1,d1;
  int ti,tf,t;


  extern int randint(int);
  extern int nextdir(int);

  p = randint(VOL);
  ti= (p/(LX*LY))%LT;

  for(t=0;t<LT;t++) corr[t]=0.0;

  sus = 0;
  spsus = 0;
  p1 = p;
  while(1)
    {
      tf = (p1/(LX*LY))%LT;
      d1 = nextdir(p1);
      sus += 1.0;
      if(ti==tf) spsus++;

      corr[(tf-ti+LT)%LT] += 0.5;
      corr[(tf-ti+LT)%LT] += 0.5;
      if(d1 > 0) conf[d1-1][p1] += 1;
      else if(d1 < 0) conf[-d1-1][next[DIM+d1][p1]] -= 1;

      p1 = next[DIM+d1][p1];

      if(p1 == p) break;
    }

  return sus;
}

int nextdir(int p)
{
  int d1;
  double prob;

  extern double bndwt(int,int);
  extern int randint(int);
  extern double ran3(void);

  d1 = randint(2*DIM) + 1;
  if(d1 > DIM) d1 = DIM-d1;

  WT = bndwt(p,d1);

  prob = ran3();
  if(prob < WT) return d1;
  else return 0;
}


double bndwt(int p, int d)
{
  int k;
  double w;

  if(d == 0) w = 1;
  else if(d > 0)
    {
      k = conf[d-1][p];
      w = IB[KMAX+k+1]/IB[KMAX+k];
      if(d == DIM) w *= exp(MU);
    }
  else
    {
      k = conf[-d-1][next[DIM+d][p]];
      w = IB[KMAX+k-1]/IB[KMAX+k];
      if(d == -DIM) w *= exp(-MU);
    }

  return w;
}

void printconf(void)
{
  int p,d;

  for(p=0;p<VOL;p++)
    {
      printf("%4d :",p);
      for(d=0;d<DIM;d++)
	{
	  printf("\t (%4d,%4d)",d+1,conf[d][p]);
	}
      printf("\n");
    } 
}

void currsus(void)
{
  int p,d;
  double WIND[DIM];

  for(d=0;d<DIM;d++) WIND[d] = 0;

  for(p=0;p<VOL;p++)
    {
      for(d=0;d<DIM;d++) WIND[d] += conf[d][p];
    }

  rho = WIND[DIM-1]/((double)VOL);

  rhosus = WIND[DIM-1]*WIND[DIM-1]/((double)VOL);

  wsus = 0;
  for(d=0;d<DIM-1;d++)
    wsus += WIND[d]*WIND[d]/((double)(DIM-1)*(double)VOL);

}
