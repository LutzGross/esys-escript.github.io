#include <math.h>

#include "ext_math.h"

double acosh(double x)
{
	return log(x+sqrt(x*x-1.0));
}

double asinh(double x)
{
	return log(x+sqrt(x*x+1.0));
}

double atanh(double x)
{
	return log((1.0+x)/(1.0-x))/2.0;
}
