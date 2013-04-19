#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <math.h>

double* generate_rv(int dimensionality);
double covar(double* data1, double* data2, int lenn);
double var(double* data, int len);
double expval(double* data, int len);

int main(int argc, char *argv[])
{
    double data[] = {1,2,3,4,5,6,7,8,9,10};
    printf("ev %lf\n", expval(data, sizeof(data)/sizeof(double)));
    printf("var %lf\n", var(data, sizeof(data)/sizeof(double)));
    printf("expected value: %lf, variance: %lf\n", expval(data, sizeof(data)/sizeof(double)), var(data, sizeof(data)/sizeof(double)));
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
    double sum;
    
    for (i = 0; i < len; ++i) {
	sum += data[i];
	printf("sum is %lf\n", sum);
    }
    printf("sum over len %lf\n", sum/len);
    return sum/len;
}

// Calculate an estimate of the variance for a given set of data
double var(double* data, int len)
{
    int i;
    double sum, exp_val;
    exp_val = expval(data, len);
        
    for (i = 0; i < len; ++i) {
	sum += pow(data[i] - exp_val, 2);
    }
    return sum/len;
}

// Calculate an estimate of the covariance of a 2-dimensional 
double covar(double* data1, double* data2, int len)
{
    // If X1 = X2, then CoVar[X1,X2] = Var[X1]
    if (data1 == data2){
	return var(data1, len);
    }
    
    int i;
    double sum, exp_val_1, exp_val_2;
    exp_val_1 = var(data1, len);
    exp_val_2 = var(data2, len);
    
    for (i = 0; i < len; ++i) {
	sum += (data1[i] - exp_val_1) * (data2[i] - exp_val_2);
    }

    return sum/len;
}
