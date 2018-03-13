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

	if (iter == 0)
	{
		int ii;
		for (ii = 0; ii < n*(n-1)/2.0; ii++)
		{
			dvContacts[ii] = 0.0;
		}
		for (ii = 0; ii < n; ii++)
		{
			dvBoundary[ii] = 0.0;
		}
		return 0.0;
	}

    double zeta = 0.0;

//
//  A FAIRE.... :-)    Difficile, difficile :-)
//

	int i, j;
	int k = 0;
	double *norm = (double *) malloc(2 * sizeof(double));
	double gamma, rCentre, deltax, deltay, deltav, vn, miSum, mjSum, gammaInner, gammaOuter, deltavInner, deltavOuter;

	for (i = 0; i < n; i++)
	{
		// Distances
		rCentre = sqrt(pow(x[i], 2) + pow(y[i], 2));
		gammaOuter = rOut - rCentre - r[i];
		gammaInner = rCentre - rIn - r[i];

		// Normale et vitesse normale
		norm[0] = x[i] / rCentre;
		norm[1] = y[i] / rCentre;
		vn = vx[i] * norm[0] + vy[i] * norm[1];

		// Increments/changements de vitesse
		deltavOuter = fmax(0, vn + dvBoundary[i] - (gammaOuter / dt));
		deltavInner = fmax(0, -vn - dvBoundary[i] - (gammaInner / dt));
		deltav = deltavOuter - deltavInner - dvBoundary[i];
		vx[i] -= deltav * norm[0];
		vy[i] -= deltav * norm[1];

		// Autres parametres
		dvBoundary[i] += deltav;
		zeta = fmax(zeta, abs(deltav));

		for (j = i + 1; j < n; j++)
		{
			// Differences de position
			deltax = x[j] - x[i];
			deltay = y[j] - y[i];

			// Distances
			rCentre = sqrt(pow(deltax,2) + pow(deltay,2));
			gamma = rCentre - r[i] - r[j];

			// Normale et vitesse normale
			norm[0] = deltax / rCentre;
			norm[1] = deltay / rCentre;
			vn = (vx[i] - vx[j]) * norm[0] + (vy[i] - vy[j]) * norm[1];

			// Increment de vitesse
			deltav = fmax(0.0, (vn + dvContacts[k] - gamma / dt)) - dvContacts[k];

			// Calcul des masses reduites
			miSum = m[i] / (m[i] + m[j]);
			mjSum = m[j] / (m[i] + m[j]);

			// Calcul des changements de vitesse
			vx[i] -= deltav * norm[0] * mjSum;
			vx[j] += deltav * norm[0] * miSum;
			vy[i] -= deltav * norm[1] * mjSum;
			vy[j] += deltav * norm[1] * miSum;

			// Paramètres supplementaires
			dvContacts[k++] += deltav;
			zeta = fmax(zeta, fabs(deltav));
		}
	}
	free(norm);
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
