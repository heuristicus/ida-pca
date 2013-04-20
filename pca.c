#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <math.h>

double* generate_rv(int dimensionality);
double covar(double* data1, double* data2, int len);
double var(double* data, int len);
double expval(double* data, int len);
double* centre(double* data, int len);
void print_arr(double* data, int len);

int main(int argc, char *argv[])
{
    double data1[] = {1,2,3,4,5,6,7,8,9,10};
    double data2[] = {10,8,5,5,6,3,2,5,3,1};
    int len = sizeof(data1)/sizeof(double);
    
    printf("E[data1]: %lf, Var[data1]: %lf\n", expval(data1, len), var(data1, len));
    printf("E[data2]: %lf, Var[data2]: %lf\n", expval(data2, len), var(data2, len));
    printf("covar[data1,data1] %lf\n", covar(data1, data2, len));
    printf("covar[data1,data2] %lf\n", covar(data1, data2, len));
    double* cd1 = centre(data1, len);
    double* cd2 = centre(data2, len);
    print_arr(cd1, len);
    print_arr(cd2, len);
    
    return 0;
}


// Generate a random variable with the given dimensionality from gaussians with
// randomly initialised mean and standard deviation
double* generate_rv(int dimensionality)
{
    // Store gaussians from which RV values will be generated
    double* randvar = NULL;
    
    int i;
    
    // Generate a value for the RV in each dimension
    for (i = 0; i < dimensionality; ++i) {
	
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
    double sum, exp_val_1, exp_val_2;
    exp_val_1 = expval(data1, len);
    exp_val_2 = expval(data2, len);
    
    for (i = 0; i < len; ++i) {
	sum += (data1[i] - exp_val_1) * (data2[i] - exp_val_2);
    }

    return sum/(len-1);
}

// Centres a set of samples by subtracting the mean of the sample from each value.
double* centre(double* data, int len)
{
    double* ret = malloc(sizeof(double) * len);
    
    int i;
    double mean = expval(data, len);
    
    for (i = 0; i < len; ++i) {
	ret[i] = data[i] - mean;
    }
    
    return ret;
}

// Prints an array
void print_arr(double* data, int len)
{
    int i;
    
    for (i = 0; i < len; ++i) {
	printf("%lf\n", data[i]);
    }
}
