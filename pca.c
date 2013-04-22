#include "pca.h"

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
    
    dblarr* cx1 = centre(x1->data, x1->len);
    dblarr* cx2 = centre(x2->data, x2->len);

    free_dblarr(cx1);
    free_dblarr(cx2);
    
//    gnuplot_2d(x1->data, x2->data, x1->len, RV_PLOT_CENTRE);
    int rv_num = 2; // Number of random variables
    int rs_num = 10; // Number of realisations for each
    dblarr* means = random_vector(rv_num, 10);
    dblarr* stdevs = random_vector(rv_num, 20);
    
    dblarr_2d* rv_set = generate_rv_set(means->data, stdevs->data, rv_num, rs_num);
    
    print_dblarr(rv_set);

    /* print_arr(means->data, means->len); */
    /* print_arr(stdevs->data, stdevs->len); */
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

// Generates a set of rv_num random variables, each with num_re realisations which
// are defined by the means and stdevs arrays, which should be of length rv_num.
dblarr_2d* generate_rv_set(double* means, double* stdevs, int rv_num, int num_re)
{
    dblarr_2d* rvset = init_dblarr_2d(num_re, rv_num);
    
    int i;
    
    for (i = 0; i < rvset->len; ++i) {
	rvset->data[i] = generate_rv(num_re, means[i], stdevs[i]);
    }
    
    return rvset;
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

// Generates a vector of the requested length containing random values, sampled
// from U[-1,1). Each value is multiplied by the given multiplier value, giving
// an effective distribution of U[-multiplier, +multiplier)
dblarr* random_vector(int len, double multiplier)
{
    dblarr* ret = init_dblarr(len);
    
    int i;
    
    for (i = 0; i < len; ++i) {
	ret->data[i] = (gsl_rng_uniform(rng) * 2 - 1) * multiplier;
    }
    
    return ret;
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
// to dblarr structs. num_vars indicates the length of each internal dblarr
// struct. All of them must have the same length. The num_arrs parameter indicates
// the number of dblarr structs this dblarr_2d struct will contain.
dblarr_2d* init_dblarr_2d(int num_vars, int num_arrs)
{
    int i;
    dblarr_2d* ret = malloc(sizeof(dblarr_2d));
    ret->data = malloc(sizeof(dblarr*) * num_arrs);
    ret->datalen = num_vars;
    ret->len = num_arrs;
    
    for (i = 0; i < num_arrs; ++i) {
	ret->data[i] = init_dblarr(num_vars);
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
    char* fname = GNUPLOT_FNAME;
    FILE *fp = fopen(fname, "w");

    int i;

    double* d1 = data1;
    double* d2 = data2;
    
    dblarr* cd1 = NULL;
    dblarr* cd2 = NULL;

    if (centred){
	cd1 = centre(data1, len);
	cd2 = centre(data2, len);
	d1 = cd1->data;
	d2 = cd2->data;
    }
	
    for (i = 0; i < len; ++i) {
	fprintf(fp, "%lf %lf\n", d1[i], d2[i]);
    }

    fclose(fp);

    printf("RV data output to %s\n", fname);

    double xmax = arr_max(d1, len);
    double xmin = arr_min(d1, len);
    double ymax = arr_max(d2, len);
    double ymin = arr_min(d2, len);

    double maxy = absmax(ymax, ymin);
    double maxx = absmax(xmax, xmin);

    double bothmax = max(maxy, maxx);

    double mval = bothmax + bothmax/10;

    /* printf("xmax %lf xmin %lf ymax %lf ymin %lf\n", xmax, xmin, ymax, ymin); */
    /* printf("maxx %lf, maxy %lf\n", maxx, maxy); */
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

// Prints a dblarr struct
void print_dblarr(dblarr_2d* arr)
{
    int i;
    
    for (i = 0; i < arr->len; ++i) {
	printf("Array %d:\n", i);
	print_arr(arr->data[i]->data, arr->data[i]->len);
    }
}
