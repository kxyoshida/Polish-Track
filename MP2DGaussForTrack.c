/* KY modified on 14 Feb 2012 */
/* Adaoted from Craig Markwardt   */
/* Test routines for mpfit library
   $Id: testmpfit.c,v 1.6 2010/11/13 08:18:02 craigm Exp $
*/

#include "Polish_Track.h"

#include "mpfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define PI 3.14159265


/* This is the private data structure which contains the data points
   and their uncertainties */
struct vars_struct {
  double *x;
  double *y;
  double *g;
  double *eg;
};

/* Simple routine to print the fit results */
void printresult(double *x, mp_result *result) 
{
  int i;

  if ((x == 0) || (result == 0)) return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	 result->bestnorm, result->nfunc-result->nfree);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  for (i=0; i<result->npar; i++) {
    printf("  P[%d] = %f +/- %f\n", 
	   i, x[i], result->xerror[i]);
  }
}

/*
def twodgaussian(inpars, circle=False, rotate=True, vheight=True, shape=None):
    """Returns a 2d gaussian function of the form:
        x' = numpy.cos(rota) * x - numpy.sin(rota) * y
        y' = numpy.sin(rota) * x + numpy.cos(rota) * y
        (rota should be in degrees)
        g = b + a * numpy.exp ( - ( ((x-center_x)/width_x)**2 + ((y-center_y)/width_y)**2 ) / 2 )
*/

int gaussfunc(int m, int n, double *p, double *dg, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *g, *eg;
  double xc, yc, sigx2, sigy2, rota;
  double xr, yr, cxr, cyr;

  x = v->x;
  y = v->y;
  g = v->g;
  eg = v->eg;

  sigx2 = p[4] * p[4];
  sigy2 = p[5] * p[5];
  rota = p[6]/180.0 * PI;   // degree translated into radian

  cxr = p[2] * cos(rota) - p[3] * sin(rota);
  cyr = p[2] * sin(rota) + p[3] * cos(rota);

  for (i=0; i<m; i++) {
    xr = x[i] * cos(rota) - y[i] * sin(rota);
    yr = x[i] * sin(rota) + y[i] * cos(rota);
    xc = xr - cxr;
    yc = yr - cyr;
    dg[i] = ( g[i] - (p[0] + p[1] * exp ( - ( xc*xc/sigx2 + yc*yc/sigy2 ) / 2.0)) )/eg[i];
  }
  return 0;
}


/* Test harness routine, which contains test gaussian-peak data 

   Example of fixing two parameter

   Commented example of how to put boundary constraints
*/
int testgaussfix(double* p, int m, double* data)
{
  double perror[7];			   /* Returned parameter errors */
  mp_par pars[7];			   /* Parameter constraints */         
  int i;
  struct vars_struct v;
  int status = 0;
  mp_result result;
  double x[m], y[m], g[m], eg[m];

  memset(&result,0,sizeof(result));      /* Zero results structure */
  result.xerror = perror;

  memset(pars,0,sizeof(pars));        /* Initialize constraint structure */

  for (i=0; i<7; i++) {
    pars[i].fixed = ((int)p[i+7] == 0) ? 0 : 1 ;                   /* Fix parameter 1 */
    pars[i].limits[0] = p[14+i];
    pars[i].limited[0] = ((int)p[21+i] == 0) ? 0 : 1 ;
    pars[i].limits[1] = p[28+i];
    pars[i].limited[1] = ((int)p[35+i] == 0) ? 0 : 1 ;
  }

  for (i = 0; i < m; i++) {
    x[i] = data[i];
    y[i] = data[m+i];
    g[i] = data[2*m+i];
    eg[i] = 1.0;
  }

  v.x = x;
  v.y = y;
  v.g = g;
  v.eg = eg;

  /* Call fitting function for m data points and 7 parameters (no
     parameters fixed) */
  status = mpfit(gaussfunc, m, 7, p, pars, 0, (void *) &v, &result);

  //  printf("*** testgaussfit status = %d\n", status);
  //  printresult(p, &result);

  p[7] = result.xerror[0]; // p0 error
  p[8] = result.xerror[1]; // p1 error
  p[9] = result.xerror[2]; // p2 error
  p[10] = result.xerror[3]; // p3 error
  p[11] = result.xerror[4]; // p4 error
  p[12] = result.xerror[5]; // p5 error
  p[13] = result.xerror[6]; // p6 error
  p[14] = result.bestnorm; // CHI-SQUARE
  p[15] = result.nfunc-result.nfree; // DOF
  p[16] = result.npar; // NPAR
  p[17] = result.nfree; // NFREE
  p[18] = result.npegged; // NPEGGED
  p[19] = result.niter; // NITER
  p[20] = result.nfev; //NFEV

  return status;
}

JNIEXPORT jint JNICALL Java_Polish_1Track_calc (JNIEnv *env, jobject obj, jdoubleArray jparams, jdoubleArray jdata)
{
  int err = 0;
  int i;
  jint nparams = (*env)->GetArrayLength(env,jparams);
  jint ndata = (*env)->GetArrayLength(env,jdata);
  jdouble outparams[nparams];
  if (nparams>0 && ndata>0) {
    jdouble cparams[nparams];
    jdouble cdata[ndata];
    (*env)->GetDoubleArrayRegion(env, jparams, 0, nparams, cparams);
    (*env)->GetDoubleArrayRegion(env, jdata, 0, ndata, cdata);
    err = testgaussfix((double*)cparams, (int)(ndata/3), (double*)cdata);
    printf("MP2DGaussFit: err=%d", err);
    //    for (i = 0; i< nparams; i++)
    //      printf("cparams=%f\n", (double)cparams[i]);
    (*env)->SetDoubleArrayRegion(env, jparams, 0, nparams, cparams);
  }
  return err;
}
