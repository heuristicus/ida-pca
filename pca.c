#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#define RV_PLOT_CENTRE 1
#define RV_PLOT_DEFAULT 0

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
dblarr_2d* init_dblarr_2d(int datalen, int len);
void free_dblarr(dblarr* arr);
void gnuplot_2d(double* data1, double* data2, int len, int centre);
double arr_max(double* data, int len);
double arr_min(double* data, int len);
double min(double a, double b);
double max(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);

int rng_init; // Track whether rng initialised
gsl_rng* rng;

int main(int argc, char *argv[])
{
    atexit(cleanup);
    /* double data1[] = {1,2,3,4,5,6,7,8,9,10}; */
    /* double data2[] = {10,8,5,5,6,3,2,5,3,1}; */
    /* int len = sizeof(data1)/sizeof(double); */
    
    /* printf("E[data1]: %lf, Var[data1]: %lf\n", expval(data1, len), var(data1, len)); */
    /* printf("E[data2]: %lf, Var[data2]: %lf\n", expval(data2, len), var(data2, len)); */
    /* printf("covar[data1,data1] %lf\n", covar(data1, data2, len)); */
    /* printf("covar[data1,data2] %lf\n", covar(data1, data2, len)); */

    dblarr* x1 = generate_rv(50, 5, 10);
    dblarr* x2 = generate_rv_simcov(x1->data, x1->len, 10, 5, -0.5);
//    dblarr* x2 = generate_rv(50, 10, 1);
    
    /* dblarr* cx1 = centre(x1->data, x1->len); */
    /* dblarr* cx2 = centre(x2->data, x2->len); */
    
    gnuplot_2d(x1->data, x2->data, x1->len, RV_PLOT_CENTRE);
//    print_arr(x1->data, x1->len);
//    print_arr(x2->data, x2->len);
    
    free_dblarr(x1);
    free_dblarr(x2);
    
    return 0;
}

// Initialises the GSL random number generator with the given seed, or
// the current time if the seed is zero
void init_rng(int seed)
{
    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, seed == 0 ? time(NULL) : seed);
    	
    rng_init = 1;
}

// Frees the random number generator on exit
void cleanup()
{
    gsl_rng_free(rng);
}

// Generate the given number of realisations of a random variable from gaussians with
// the given mean and stdev
dblarr* generate_rv(int realisations, double mean, double stdev)
{
    int i;
    dblarr* randvar = init_dblarr(realisations);
    
    if (!rng_init){
	init_rng(0);
    }

    for (i = 0; i < realisations; ++i) {
	randvar->data[i] = mean + gsl_ran_gaussian(rng, stdev);
    }

    return randvar;
}

// Generates a random variable which covaries with the given variable based on the
// cov_strength parameter. A large negative value will result in a strongly negatively
// correlated value, and a large positive value will result in a strongly positively
// correlated value. Passing a value of zero is equivalent to calling the generate_rv
// function with the given mean and standard deviation.
dblarr* generate_rv_simcov(double* cov_var, int len, double mean, double stdev,
			   double cov_strength)
{
    int i;
    dblarr* randvar = init_dblarr(len);
    
    if (!rng_init){
	init_rng(0);
    }

    // Generate a value for the RV in each dimension
    for (i = 0; i < len; ++i) {
	randvar->data[i] = mean + gsl_ran_gaussian(rng, stdev) + cov_strength * cov_var[i];
    }

    return randvar;    
}

// Calculate an estimate of the expected value for a given set of data
double expval(double* data, int len)
{
    int i;
    double sum = 0;
    
    for (i = 0; i < len; ++i) {
	sum += data[i];
    }
    return sum/len;
}

// Calculate the sample variance for the given set of data
double var(double* data, int len)
{
    int i;
    double sum = 0;
    double exp_val = expval(data, len);
        
    for (i = 0; i < len; ++i) {
	double res = pow(data[i] - exp_val, 2);
	sum += res;
    }

    return sum/(len-1);
}

// Calculate an estimate of the covariance of two sets of samples
double covar(double* data1, double* data2, int len)
{
    // If X1 = X2, then CoVar[X1,X2] = Var[X1]
    if (data1 == data2){
	return var(data1, len);
    }
    
    int i;
    double sum = 0;
    double exp_val_1 = expval(data1, len);
    double exp_val_2 = expval(data2, len);
    
    for (i = 0; i < len; ++i) {
	sum += (data1[i] - exp_val_1) * (data2[i] - exp_val_2);
    }

    return sum/(len-1);
}

// Centres a set of samples by subtracting the mean of the sample from each value.
dblarr* centre(double* data, int len)
{
    dblarr* ret = init_dblarr(len);
    
    int i;
    double mean = expval(data, len);
    
    for (i = 0; i < len; ++i) {
	ret->data[i] = data[i] - mean;
    }
    
    return ret;
}

// Initialises a double array with the given length
dblarr* init_dblarr(int len)
{
    dblarr* ret = malloc(sizeof(dblarr) * len);
    ret->data = malloc(sizeof(double) * len);
    ret->len = len;
    return ret;
}

// Initialises a dblarr_2d struct, which contains an array of pointers
// to dblarr structs. datalen indicates the length of each internal dblarr
// struct. All of them must have the same length. The len parameter indicates
// the number of dblarr structs this dblarr_2d struct will contain.
dblarr_2d* init_dblarr_2d(int datalen, int len)
{
    int i;
    dblarr_2d* ret = malloc(sizeof(dblarr_2d));
    ret->data = malloc(sizeof(dblarr*) * len);
    ret->datalen = datalen;
    ret->len = len;
    
    for (i = 0; i < len; ++i) {
	ret->data[i] = init_dblarr(datalen);
    }
    return ret;
}

// Frees a dblarr struct
void free_dblarr(dblarr* arr)
{
    free(arr->data);
    free(arr);
}

// Plots the values of two random variables as a scatter plot in gnuplot.
// The value of centred defines whether to print centred data or not.
// RV_PLOT_CENTRE will centre values, RV_PLOT_DEFAULT will not. Passing
// 1 or 0 has the same effect.
void gnuplot_2d(double* data1, double* data2, int len, int centred)
{
    FILE *fp = fopen("dfile.dat", "w");

    int i;

    double* d1 = data1;
    double* d2 = data2;
    
    dblarr* cd1 = NULL;
    dblarr* cd2 = NULL;

    if (centred){
	dblarr* cd1 = centre(data1, len);
	dblarr* cd2 = centre(data2, len);
	d1 = cd1->data;
	d2 = cd2->data;
    }
	
    for (i = 0; i < len; ++i) {
	fprintf(fp, "%lf %lf\n", d1[i], d2[i]);
    }

    fclose(fp);

    double xmax = arr_max(d1, len);
    double xmin = arr_min(d1, len);
    double ymax = arr_max(d2, len);
    double ymin = arr_min(d2, len);

    double maxy = absmax(ymax, ymin);
    double maxx = absmax(xmax, xmin);

    double bothmax = max(maxy, maxx);

    double mval = bothmax + bothmax/10;

    printf("xmax %lf xmin %lf ymax %lf ymin %lf\n", xmax, xmin, ymax, ymin);
    printf("maxx %lf, maxy %lf\n", maxx, maxy);
    FILE *pipe = popen("gnuplot -persist","w");
    /* fprintf(pipe, "set xrange [%lf:%lf]\n", xmin + xmin/10, xmax + xmax/10); */
    /* fprintf(pipe, "set yrange [%lf:%lf]\n", ymin + ymin/10, ymax + ymax/10); */
    /* fprintf(pipe, "set xrange [%lf:%lf]\n", -(maxx + maxx/10), maxx + maxx/10); */
    /* fprintf(pipe, "set yrange [%lf:%lf]\n", -(maxy + maxy/10), maxy + maxy/10); */
    fprintf(pipe, "set xrange [%lf:%lf]\n", -mval, mval);
    fprintf(pipe, "set yrange [%lf:%lf]\n", -mval, mval);
    fprintf(pipe, "plot 'dfile.dat' using 1:2 with points\n");
    fclose(pipe); 

    if (centred){
	free_dblarr(cd2);
	free_dblarr(cd1);
    }
}

// Find the largest value in an array
double arr_max(double* data, int len)
{
    int i;
    double max = data[0];
    
    for (i = 0; i < len; ++i) {
	if (data[i] > max)
	    max = data[i];
    }

    return max;
}

// Find the smallest value in an array
double arr_min(double* data, int len)
{
    int i;
    double min = data[0];
    
    for (i = 0; i < len; ++i) {
	if (data[i] < min)
	    min = data[i];
    }

    return min;    
}

double min(double a, double b)
{
    return a < b ? a : b;
}

double max(double a, double b)
{
    return a > b ? a : b;
}

double absmax(double a, double b)
{
    a = fabs(a);
    b = fabs(b);
    
    return max(a, b);
}

double absmin(double a, double b)
{
    a = fabs(a);
    b = fabs(b);
    
    return min(a, b);
}


// Prints an array
void print_arr(double* data, int len)
{
    int i;
    
    for (i = 0; i < len; ++i) {
	printf("%lf\n", data[i]);
    }

    printf("\n");
}
