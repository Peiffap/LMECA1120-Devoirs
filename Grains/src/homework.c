#include"fem.h"


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
    int n = myGrains->n;
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;

    double zeta = 0.0;

//
//  A FAIRE.... :-)    Difficile, difficile :-)
//

	int i, j;
	int k = 0;
	double norm[2];
	double gamma, rCentre, deltax, deltay, deltav, vn, miSum, mjSum, gammaInner, gammaOuter, deltavb, deltavc;

	for (i = 0; i < n; i++)
	{
		// deltavb = 0;

		// Distances
		rCentre = sqrt((x[i] * x[i]) + (y[i] * y[i]));
		gammaOuter = rOut - rCentre - r[i];
		gammaInner = rCentre - rIn - r[i];

		// Normale et vitesse normale
		norm[0] = x[i] / rCentre;
		norm[1] = y[i] / rCentre;
		vn = vx[i] * norm[0] + vy[i] * norm[1];

		// Increments de vitesse
		deltav = fmax(fmax(0, vn + deltavb - (gammaOuter / dt)), -vn - deltavb - (gammaInner / dt)) - deltavb;
		vx[i] -= deltav * norm[0];
		vy[i] -= deltav * norm[1];
		deltavb += deltav;

		// Autres parametres
		zeta = fmax(zeta, abs(deltav));
		dvBoundary[i] += deltav;

		for (j = i + 1; j < n; j++)
		{
			deltavc = 0;

			// Differences de position
			deltax = x[j] - x[i];
			deltay = y[j] - y[i];

			// Distances
			rCentre = sqrt(deltax * deltax + deltay * deltay);
			gamma = rCentre - r[i] - r[j];

			// Normale et vitesse normale
			norm[0] = deltax / rCentre;
			norm[1] = deltay / rCentre;
			vn = (vx[i] - vx[j]) * norm[0] + (vy[i] - vy[j]) * norm[1];

			// Increment de vitesse
			deltav = fmax(0, (vn + deltavc - gamma / dt)) - deltavc;

			// Calcul des masses reduites
			miSum = m[i] / (m[i] + m[j]);
			mjSum = m[j] / (m[i] + m[j]);

			// Calcul des changements de vitesse
			vx[i] -= deltav * norm[0] * miSum;
			vx[j] += deltav * norm[0] * mjSum;
			vy[i] -= deltav * norm[1] * miSum;
			vy[j] += deltav * norm[1] * mjSum;

			// Paramètres supplementaires
			deltavc += deltav;
			zeta = fmax(zeta, fabs(deltav));
			dvContacts[k++] += deltav;
		}
	}

    return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
    int n = myGrains->n;
    int i,iter = 0;
    double zeta;
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];


//
// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
//

//
//  A FAIRE.... :-)    La loi du grand Newton : so easy !
//

	for (i = 0; i < n; i++)
	{
		vx[i] += (gx - (gamma * vx[i] / m[i])) * dt;
		vy[i] += (gy - (gamma * vy[i] / m[i])) * dt;
	}

//
// -2- Correction des vitesses pour tenir compte des contacts
//
    do
	{
        zeta = femGrainsContactIterate(myGrains,dt,iter);
        iter++;
	}
    while ((zeta > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,zeta);

//
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//
    for (i = 0; i < n; ++i)
	{
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
	}
}


#endif
