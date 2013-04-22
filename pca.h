#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#define RV_PLOT_CENTRE 1
#define RV_PLOT_DEFAULT 0
#define GNUPLOT_FNAME "dfile.dat"

typedef struct
{
    double* data;
    int len;
} dblarr;

typedef struct
{
    dblarr** data;
    int datalen;
    int len;
} dblarr_2d;

dblarr* generate_rv(int realisations, double mean, double stdev);
dblarr* generate_rv_simcov(double* cov_var, int len, double mean, double stdev, double cov_strength);
double covar(double* data1, double* data2, int len);
double var(double* data, int len);
double expval(double* data, int len);
dblarr* centre(double* data, int len);
void print_arr(double* data, int len);
void cleanup();
void init_rng(int seed);
dblarr* init_dblarr(int len);
dblarr_2d* init_dblarr_2d(int num_vars, int num_arrs);
void free_dblarr(dblarr* arr);
void gnuplot_2d(double* data1, double* data2, int len, int centre);
double arr_max(double* data, int len);
double arr_min(double* data, int len);
double min(double a, double b);
double max(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);
dblarr_2d* generate_rv_set(double* means, double* stdevs, int rv_num, int num_re);
dblarr* random_vector(int len, double multiplier);
void print_dblarr(dblarr_2d* arr);
