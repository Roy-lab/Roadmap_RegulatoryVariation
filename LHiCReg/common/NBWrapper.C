#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
#include "gsl/gsl_randist.h"
#include "gsl/gsl_multimin.h"
#include "NBWrapper.H"
//First need to define a multi-min function which will be computing the likelihood for given set of parameters

double
negBinomLL(const gsl_vector* v, void* params)
{
	double n=gsl_vector_get(v,0);	

}
NBWrapper::NBWrapper()
{
}

NBWrapper::~NBWrapper()
{
}

int 
NBWrapper::estimateParams(double & mean, double& prob, vector<int>& countValues)
{
	gsl_vector* v=gsl_vector_alloc(1);
	double m=0;
	for(int i=0;i<countValues.size();i++)
	{
		m=m+(double)countValues[i];
	}
	m=m/countValues.size();
	//Use MATLAB's  initialization to set s2
	double s2=0;
	for(int i=0;i<countValues.size();i++)
	{
		double diff=m-countValues[i];
		s2=s2+(diff*diff);
	}
	s2=s2/(countValues.size()-1);
	double n_hat=(m*m)/(s2-m);
	//Set initial parameters in params
	gsl_vector_set(v,0,n_hat);
	//gsl_vector_set(v,0,10);
	//gsl_vector_set(v,0,1);
	//gsl_vector_set(v,1,0.5);
	//Set the minimizer function. This is two dimensional because we have two params
	gsl_multimin_function minex_func;
	//minex_func.n=2;
	minex_func.n=1;
	minex_func.f=&negBinomLL;
	minex_func.params=(void*)(&countValues);
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	//Need to set the step size
	//gsl_vector* ss = gsl_vector_alloc (2);
	gsl_vector* ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, 0.01);

	//Set up the minimization object
 	gsl_multimin_fminimizer *s = NULL;
	//s=gsl_multimin_fminimizer_alloc (T, 2);
	s=gsl_multimin_fminimizer_alloc (T, 1);
	gsl_multimin_fminimizer_set (s, &minex_func,v, ss);
	int iter=0;
	int status=0;
	double size=0;
	double pest=0;
	double nest=0;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) 
			break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-5);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
        	}
		//printf ("%d %f %f f() = %f size = %f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval, size);
		nest=gsl_vector_get(s->x,0);
		pest=nest/(nest+m);
		printf ("%d %f  prob=%f f() = %f size = %f\n", iter, nest,pest,  s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);
	gsl_vector_free(v);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	mean=nest;
	prob=pest;
	return status;
}
	
