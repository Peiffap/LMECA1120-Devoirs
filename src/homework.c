#include <stdio.h>
#include <math.h>

#ifdef graphic
#include "glfem.h"
#endif



double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
	double xLoc[3] = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
	double yLoc[3] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
	double wLoc[3] = { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };
	int i;
	for (i = 0; i < 3; i++)
	{
		double X = x[0] * (1 - xLoc[i] - yLoc[i]) + x[1] * (xLoc[i]) + x[2] * (yLoc[i]);
		double Y = y[0] * (1 - xLoc[i] - yLoc[i]) + y[1] * (xLoc[i]) + y[2] * (yLoc[i]);
		double W = wLoc[0] * (1 - xLoc[i] - yLoc[i]) + wLoc[1] * (xLoc[i]) + wLoc[2] * (yLoc[i]);
		I = I + W * f(X, Y);
	}

//
// ... A modifier :-)
//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//


#ifdef graphic

    glColor3f (1.0,1.0,1.0); glfemDrawElement(x,y,3);
    glColor3f (1.0,0.0,0.0); glfemDrawNodes(x,y,3);
    glColor3f (0.0,0.0,1.0); glfemDrawNodes(xLoc,yLoc,3);
    
#endif

    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
	double I = integrate(x, y, f);
//
//
//    
     
    return I;
}


