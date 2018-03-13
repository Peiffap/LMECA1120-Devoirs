#include"fem.h"


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
	int n = myGrains->n;
	double *dvBoundary = myGrains->dvBoundary;
	double *dvContacts = myGrains->dvContacts;

	if (iter == 0)
	{
		int ii;
		for (ii = 0; ii < n; ii++)
		{
			dvBoundary[ii] = 0.0;
			dvContacts[ii] = 0.0;
		}
		for (ii = n; ii < n*(n-1)/2.0; ii++)
		{
			dvContacts[ii] = 0.0;
		}
		return 0.0;
	}

	//
	//  A FAIRE.... :-)    Difficile, difficile :-)
	//

	double *x          = myGrains->x;
	double *y          = myGrains->y;
	double *m          = myGrains->m;
	double *r          = myGrains->r;
	double *vy         = myGrains->vy;
	double *vx         = myGrains->vx;
	double rIn         = myGrains->radiusIn;
	double rOut        = myGrains->radiusOut;

	int i, j;
	int k = 0;
	double gamma, rCentre, deltax, deltay, deltav, vn, miSum, mjSum, gammaInner, gammaOuter, deltavInner, deltavOuter, nx, ny, mSum;
	double zeta = 0.0;

	for (i = 0; i < n; i++)
	{
		// Distances
		rCentre = sqrt(pow(x[i], 2) + pow(y[i], 2));
		gammaOuter = rOut - rCentre - r[i];
		gammaInner = rCentre - rIn - r[i];

		// Normale et vitesse normale
		nx = x[i] / rCentre;
		ny = y[i] / rCentre;
		vn = vx[i] * nx + vy[i] * ny;

		// Increments/changements de vitesse
		deltavOuter = fmax(0, vn + dvBoundary[i] - (gammaOuter / dt));
		deltavInner = fmax(0, -vn - dvBoundary[i] - (gammaInner / dt));
		deltav = deltavOuter - deltavInner - dvBoundary[i];
		vx[i] -= deltav * nx;
		vy[i] -= deltav * ny;

		// Autres parametres
		dvBoundary[i] += deltav;
		zeta = fmax(zeta, fabs(deltav));

		for (j = i + 1; j < n; j++)
		{
			// Differences de position
			deltax = x[j] - x[i];
			deltay = y[j] - y[i];

			// Distances
			rCentre = sqrt(pow(deltax,2) + pow(deltay,2));
			gamma = rCentre - r[i] - r[j];

			// Normale et vitesse normale
			nx = deltax / rCentre;
			ny = deltay / rCentre;
			vn = (vx[i] - vx[j]) * nx + (vy[i] - vy[j]) * ny;

			// Increment de vitesse
			deltav = fmax(0.0, (vn + dvContacts[k] - gamma / dt)) - dvContacts[k];

			// Calcul des masses reduites
			mSum = m[i] + m[j];
			miSum = m[i] / mSum;
			mjSum = m[j] / mSum;

			// Calcul des changements de vitesse
			vx[i] -= deltav * nx * mjSum;
			vx[j] += deltav * nx * miSum;
			vy[i] -= deltav * ny * mjSum;
			vy[j] += deltav * ny * miSum;

			// Paramètres supplementaires
			dvContacts[k++] += deltav;
			zeta = fmax(zeta, fabs(deltav));
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
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
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
