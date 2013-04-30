#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_sf_bessel.h>

# define DIM 2
# define KMAX 20

int *conf[DIM];

/* next holds the nearest neighbours. If DIM=2 then -2,-1,1,2   *
 * denotes the nearest -t,-x,x,t dirn. 0 refers the site itself */
int *next[2*DIM+1];
double *corr,*avgcorr;
double ***pol,***avgpol;
double **rat;
double sus, avgsus;

double BETA,BETAS;
int LX,LT,VOL;
int SEED,RNDFLAG;
int START;
int EQUIL,SWEEP,MEAS;
int NCH,Q1,Q2;
int TMAX,TMIN;
double IB[2*KMAX+1],IBT[2*KMAX+1];
double WT;
int main(argc,argv)
     int argc;
     char *argv[];
{
  int i,j,d,k1,k2,t;
  int delta;
  char string[100];
  char beta[5];
  int readval;
  FILE *fptr;
  char file1[20]="SUS";
  char file2[20]="CORR";
  char file3[20]="POL";

  extern void initrandom(int,int);
  extern void initneighbor(void);
  extern void initbessel(void);
  extern void initconf(void);
  extern double update(void);
  extern void checkconf(void);
  extern void measure(void);
  extern void printconf(void);
  extern void currsus(void);
  extern double ***allocate3d(int, int, int);
  extern void deallocate3d(double ***, int, int, int);
  double **allocatedouble2d(int,int);
  void deallocatedouble2d(double**,int,int);
  if(argc != 1)
    {
      printf("usage xy (needs a QUEUE file)\n");
      exit(1);
    }

  fptr = fopen("QUEUE","r");
  if(fptr == NULL)
    {
      printf("ERROR opening file\n");
      exit(1);
    }
  readval=fscanf(fptr,"%s %d\n",string,&LX);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&LT);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %lf\n",string,&BETA);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %lf\n",string,&BETAS);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&NCH);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&TMAX);
  if(readval<0) printf("Error\n");
  readval=fscanf(fptr,"%s %d\n",string,&TMIN);
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
  fclose(fptr);

  VOL = LX*LT;
  delta=TMAX-TMIN;
  strcat(file1,beta);
  strcat(file2,beta);
  strcat(file3,beta);
  for(d=0;d<DIM;d++)
    conf[d] = (int *)malloc(VOL*sizeof(int));

  for(d=0;d<=2*DIM;d++)
    next[d] = (int *)malloc(VOL*sizeof(int));

  corr = (double *)malloc(LT*sizeof(double));
  avgcorr = (double *)malloc(LT*sizeof(double));
  pol = allocate3d(delta,NCH,NCH);
  avgpol = allocate3d(delta,NCH,NCH);
  rat = allocatedouble2d(NCH,NCH);

  initrandom(SEED,RNDFLAG);
  initneighbor();
  initbessel();
  initconf();

  for(k1=0;k1<EQUIL;k1++)
    {
      sus = update();
      checkconf();
    }
  printf("Equil done\n");
  fptr = fopen(file1,"w");
  fclose(fptr);
  fptr = fopen(file2,"w");
  fclose(fptr);
  fptr = fopen(file3,"w");
  fclose(fptr);
  for(k1=0;k1<MEAS;k1++)
    {
      avgsus = 0;
      for(t=0;t<LT;t++){
         avgcorr[t] = 0.0;
         if(t<delta){ for(i=0;i<NCH;i++) for(j=0;j<NCH;j++) avgpol[t][i][j]=0.0;
     }}
      for(k2=0;k2<SWEEP;k2++)
	{
	  sus = update();
          measure();
	  avgsus += sus;
          for(t=0;t<LT;t++) {
           avgcorr[t] += corr[t];
           if(t<delta){
           for(i=0;i<NCH;i++) for(j=0;j<NCH;j++)
             avgpol[t][i][j] += pol[t][i][j]; 
           }}
	  checkconf();
         }
    avgsus /= (double)SWEEP;
    for(t=0;t<LT;t++) {
      avgcorr[t] /= (double)SWEEP;
      if(t<delta){
       for(i=0;i<NCH;i++) for(j=0;j<NCH;j++)
           avgpol[t][i][j] /= (double)SWEEP;
       }}
    fptr = fopen(file1,"a");
    fprintf(fptr,"%16.8e\n",avgsus);
    fclose(fptr);
    fptr = fopen(file2,"a");
    for(t=0;t<LT;t++) fprintf(fptr,"%d %16.8e\n",t,avgcorr[t]);
    fclose(fptr);
    fptr = fopen(file3,"a");
    for(t=0;t<delta;t++) {
    fprintf(fptr,"%d",t+1);
    for(i=0;i<NCH;i++) for(j=0;j<NCH;j++)
       fprintf(fptr," %2.4le",avgpol[t][i][j]);
    fprintf(fptr,"\n");
    } 
    fclose(fptr);
    }
 //printconf();  
 for(d=0;d<DIM;d++) free(conf[d]);
 for(d=0;d<=2*DIM;d++) free(next[d]);
 free(corr); free(avgcorr);
 deallocate3d(pol,delta,NCH,NCH);
 deallocate3d(avgpol,delta,NCH,NCH);
 return 0;
}

void initneighbor(void)
{
  int p,x,t;

  for(p=0;p<VOL;p++)
    {
      x = (p%LX);
      t = (p/(LX))%LT;
      next[DIM][p] = p;
      next[DIM+1][p] = t*LX + ((x+1)%LX);
      next[DIM-1][p] = t*LX + ((x-1+LX)%LX);
    }
}

void initbessel(void)
{
  int k1;
  for(k1=0;k1<=KMAX;k1++){
      IB[KMAX+k1] = gsl_sf_bessel_In(k1,BETAS);
      IB[KMAX-k1] = IB[KMAX+k1];
      IBT[KMAX+k1]= gsl_sf_bessel_In(k1,BETA);
      IBT[KMAX-k1]= IBT[KMAX+k1];
    }
}

void initconf(void)
{
  int p,d;
  for(p=0;p<VOL;p++){
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
       printf(" %4d %4d %4d\n",d+1,conf[d][p],conf[d][next[DIM-(d+1)][p]]);
	  exit(1);
	}
    }
}

double update(void)
{
  int p,p1,d1;
  int ti,tf,t;

  extern int randint(int);
  extern int nextdir(int,int,int);

  p = randint(VOL);
  ti= (p/LX)%LT;
  for(t=0;t<LT;t++) corr[t]=0.0;

  sus = 0;
  p1 = p;
  while(1)
    {
      tf = (p1/LX)%LT;
      d1 = nextdir(p1,tf,ti);
      sus += 1.0;
    //  if(ti==tf) spsus++;

      corr[(tf-ti+LT)%LT] += 0.5;
      corr[(tf-ti+LT)%LT] += 0.5;
      if(d1 > 0) conf[d1-1][p1] += 1;
      else if(d1 < 0) conf[-d1-1][next[DIM+d1][p1]] -= 1;

      p1 = next[DIM+d1][p1];

      if(p1 == p) break;
    }

  return sus;
}

int nextdir(int p,int tf,int ti)
{
  int d1;
  double prob;

  extern double bndwt(int,int);
  extern int randint(int);
  extern double ran3(void);

  d1 = randint(2*DIM) + 1;
  if(d1 > DIM) d1 = DIM-d1;

  if((tf == ti) && (d1 == DIM)) WT = 0.0;
  else if((tf == ti) && (d1 == -DIM)) WT = 0.0;
  else WT = bndwt(p,d1);
  
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
      if(abs(k) > KMAX) {printf("Out of range\n"); exit(1); }
      if(d == DIM) w = IBT[KMAX+k+1]/IBT[KMAX+k];
      else w = IB[KMAX+k+1]/IB[KMAX+k];
    }
  else
    {
      k = conf[-d-1][next[DIM+d][p]];
      if(abs(k) > KMAX) {printf("Out of range\n"); exit(1); }
      if(d == -DIM) w = IBT[KMAX+k-1]/IBT[KMAX+k];
      else w = IB[KMAX+k-1]/IB[KMAX+k];
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

/*  The distribution of charges in the observable */
/*
q1>q2  -q2 x             x q2    q2>q1 q2 x------>------x -q2
           |             |                |   q2-q1     |       
           |             |                |             |
        q2 ^             v q2          q1 v             ^ q1
           |             |                |             |
           |   q1-q2     |                |             |
        q1 x----->-------x -q1        -q1 x             x q1  
 */
void measure(void)
{
  int q1,q2,p,k,t,tt,site;
  /* Note that the t-link across the boundary has to be neglected */
  for(tt=TMAX;tt>TMIN;tt--){
  for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++)
      pol[tt-TMIN-1][q1-1][q2-1]=0.0;
  for(p=0;p<LX;p++){
    /* calculate the staple shaped observable */    
    site = tt*LX + p;
    //printf("%d %d %d\n",tt,p,site);
    for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++){
      if(q2>q1){
         k = conf[0][site];
         rat[q1-1][q2-1] = IB[KMAX+k+q2-q1]/IB[KMAX+k];
      }
      else rat[q1-1][q2-1]=1.0; 
    }
    t=tt-1; site=next[DIM-2][site];
    //printf("%le\n",ratx[0][0]);
    //printf("%d %d %d\n",t,p,site);
    while(t>=TMIN){
    k = conf[DIM-1][site];
    for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++){
     if(q2>q1){
       rat[q1-1][q2-1] *= IBT[KMAX+k+q1]/IBT[KMAX+k];
     }
     else{
       rat[q1-1][q2-1] *= IBT[KMAX+k+q2]/IBT[KMAX+k];
     }}

    site = next[DIM-2][site];t--; 
    //printf("%le\n",ratx[0][0]);
    //printf("%d %d %d\n",t,p,site);
    }
    site = next[DIM+2][site];t++;
    //printf("%d %d %d\n",t,p,site);

    for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++){
    if(q1>q2){
    /* spatial x-link */
    k = conf[0][site];
    rat[q1-1][q2-1] *= IB[KMAX+k+q1-q2]/IB[KMAX+k];
    } }
    //printf("%le\n",ratx[0][0]);
    //printf("%le %le\n",ratx[0][0],ratx[1][1]);
    site = next[DIM+1][site];
    //printf("%d %d %d %d %d\n",t,p,sitex,sitey,sitez);
    while(t < tt){
    k = conf[DIM-1][site];
    for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++){
     if(q2>q1){
       rat[q1-1][q2-1] *= IBT[KMAX+k-q1]/IBT[KMAX+k];
     }
     else{
       rat[q1-1][q2-1] *= IBT[KMAX+k-q2]/IBT[KMAX+k];
     }}
    //printf("%le %le\n",ratx[0][0],ratx[1][1]);
    //printf("%le\n",ratx[0][0]);
    site = next[DIM+2][site];
    t++; 
    //printf("%d %d %d %d %d\n",t,p,sitex,sitey,sitez);
    } 
    //printf("%le\n",ratx[0][0]);
    for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++){
      pol[tt-TMIN-1][q1-1][q2-1] += rat[q1-1][q2-1];;
    }
   }
   for(q1=1;q1<=NCH;q1++) for(q2=1;q2<=NCH;q2++)
       pol[tt-TMIN-1][q1-1][q2-1] /= (double)LX;
  }
}

double ***allocate3d(int size1, int size2, int size3)
{
   int i,j,k;
   double ***mat;
   mat = (double ***)malloc(size1*sizeof(double**));

   for(i=0;i<size1;i++)
     mat[i]= (double **)malloc(size2*sizeof(double*));

   for(i=0;i<size1;i++)
   for(j=0;j<size2;j++)
     mat[i][j] = (double *)malloc(size3*sizeof(double));

   for(i=0;i<size1;i++)
   for(j=0;j<size2;j++)
   for(k=0;k<size3;k++)
     mat[i][j][k]=0.0;

   return mat;
}

void deallocate3d(double ***mat, int size1, int size2, int size3)
{
  int i,j;
  for(i=0;i<size1;i++)
  for(j=0;j<size2;j++)
    free(mat[i][j]);

  for(i=0;i<size1;i++)
    free(mat[i]);

  free(mat);
}

double **allocatedouble2d(int row, int col)
{
  int i,j;
  double **mat;
  mat = (double **)malloc(row*sizeof(double*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(double *)malloc(col*sizeof(double));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++)
 for(j=0;j<col;j++)
  mat[i][j]=0.0;

 return mat;

}

void deallocatedouble2d(double **mat, int row, int col)
{
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}


