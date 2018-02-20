#include <stdio.h>
#include <math.h>

#ifdef graphic
#include "glfem.h"
#endif

// Code pour la résolution du devoir 1 - Integrate du cours d'introduction aux
// méthodes d'éléments finis donné par le Professeur Vincent Legat à l'EPL,
// faculté de l'UCL.

// @author Gilles Peiffer & Thomas Reniers.
// Février 2018.


double integrate(double x[3], double y[3], double (*f) (double, double))
{
  double I = 0;
	double xLocal[3] = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
	double yLocal[3] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
	double wLocal[3] = { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };
  double xLoc[3];
  double yLoc[3];
	int i;
	for (i = 0; i < 3; i++)
	{
		double X = x[0] * (1 - xLocal[i] - yLocal[i]) + x[1] * (xLocal[i]) + x[2] * (yLocal[i]);
		double Y = y[0] * (1 - xLocal[i] - yLocal[i]) + y[1] * (xLocal[i]) + y[2] * (yLocal[i]);
		double W = wLocal[0] * (1 - xLocal[i] - yLocal[i]) + wLocal[1] * (xLocal[i]) + wLocal[2] * (yLocal[i]);
		I = I + W * f(X, Y);
    xLoc[i] = X;
    yLoc[i] = Y;
	}

//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//


#ifdef graphic

    glColor3f (1.0,1.0,1.0); glfemDrawElement(x,y,3);
    glColor3f (1.0,0.0,0.0); glfemDrawNodes(x,y,3);
    glColor3f (0.0,0.0,1.0); glfemDrawNodes(xLoc,yLoc,3);

#endif

    return I*fabs((x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]));
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if (n == 0)
    {
      return integrate(x,y,f);
    }
    else
    {
      int n1 = n - 1;
      double x0x1 = (x[0] + x[1]) / 2;
      double x1x2 = (x[1] + x[2]) / 2;
      double x2x0 = (x[2] + x[0]) / 2;
      double y0y1 = (y[0] + y[1]) / 2;
      double y1y2 = (y[1] + y[2]) / 2;
      double y2y0 = (y[2] + y[0]) / 2;
      double xa[3] = {x[0],x0x1,x2x0};
      double ya[3] = {y[0],y0y1,y2y0};
      double xb[3] = {x0x1,x[1],x1x2};
      double yb[3] = {y0y1,y[1],y1y2};
      double xc[3] = {x1x2,x2x0,x0x1};
      double yc[3] = {y1y2,y2y0,y0y1};
      double xd[3] = {x2x0,x1x2,x[2]};
      double yd[3] = {y2y0,y1y2,y[2]};
      return integrateRecursive(xa,ya,f,n1) +
      integrateRecursive(xb,yb,f,n1) +
      integrateRecursive(xc,yc,f,n1) +
      integrateRecursive(xd,yd,f,n1);
    }
}
