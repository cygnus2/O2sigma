#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

const gsl_rng_type *T;
gsl_rng *r;

void initrandom(int SEED, int rndflag)
{
  if(rndflag == 0) T = gsl_rng_mt19937;
  else T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc(T);

  gsl_rng_set(r,SEED);
}

double ran3(void)
{
  double rnd;

  rnd = gsl_rng_uniform(r);

  return rnd;
}


int randint(int b)
{
  int k;

  k = gsl_rng_uniform_int(r,b);

  return k;
}

void readrand(char *rng)
{
  FILE *fptr;

  fptr = fopen(rng,"r");
  gsl_rng_fread(fptr,r);
  fclose(fptr);
}

void writerand(char *rng)
{
  FILE *fptr;

  fptr = fopen(rng,"w");
  gsl_rng_fwrite(fptr,r);
  fclose(fptr);
}

void randfree(void)
{
  gsl_rng_free(r);
}
