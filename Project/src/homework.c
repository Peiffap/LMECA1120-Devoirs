#include"fem.h"
#include"limits.h"
#include <stdlib.h>
#include <stdio.h>


#ifndef NOCONTACTITERATE

/////START OF GRAINS

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
		for (ii = n; ii < n*(n - 1) / 2.0; ii++)
		{
			dvContacts[ii] = 0.0;
		}
		return 0.0;
	}


	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *r = myGrains->r;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double rIn = myGrains->radiusIn;
	double rOut = myGrains->radiusOut;

	int i, j;
	int k = 0;
	double gamma, rCentre, deltax, deltay, deltav, vn, miSum, mjSum, gammaInner, gammaOuter, deltavInner, deltavOuter, nx, ny, mSum;
	double zeta = 0.0;

	for (i = 0; i < n; ++i)
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

		for (j = i + 1; j < n; ++j)
		{
			// Differences de position
			deltax = x[j] - x[i];
			deltay = y[j] - y[i];

			// Distances
			rCentre = sqrt(pow(deltax, 2) + pow(deltay, 2));
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

			// Parametres supplementaires
			dvContacts[k++] += deltav;
			zeta = fmax(zeta, fabs(deltav));
		}
	}
	return zeta;
}

#endif

int ptInTriangle(double x, double y, double x0, double y0, double x1, double y1, double x2, double y2)
{
    double dX = x - x2;
    double dY = y - y2;
    double dX21 = x2 - x1;
    double dY12 = y1 - y2;
    double D = dY12 * (x0 - x2) + dX21 * (y0 - y2);
    double s = dY12 * dX + dX21 * dY;
    double t = (y2 - y0) * dX + (x0 - x2) * dY;
    if (D < 0)
	{
		if (s <= 0 && t <= 0 && s + t >= D)
		{
			return 1;
		}
		return 0;
	}
    else
	{
		if (s >= 0 && t >= 0 && s + t <= D)
		{
			return 1;
		}
		return 0;
	}
}

void vEval(double X[3], double Y[3], double x[2], double ret[2], femDiffusionProblem *theProblemX, femDiffusionProblem *theProblemY, int j)
{
	double phi[3];
	double xi = (x[0] * (Y[2] - Y[0]) - (X[2] - X[0]) * x[1] + (X[2] - X[0]) * Y[0] - X[0] * (Y[2] - Y[0])) / ((Y[2] - Y[0]) * (X[1] - X[0]) - (X[2] - X[0]) * (Y[1] - Y[0]));
	double eta = (x[1] - (Y[1] - Y[0]) * xi - Y[0]) / (Y[2] - Y[0]);

	// printf("xi %f || eta %f \n",xi,eta);

	femDiscretePhi2(theProblemX->space, xi, eta, phi);

	// printf("xi %f eta %f  phi0 %f phi1 %f phi2 %f\n", xi, eta, phi[0], phi[1], phi[2]);

	femMesh *mesh = theProblemX->mesh;
	ret[0] = 0.0;
	ret[1] = 0.0;
	int i;
	for (i = 0; i < 3; ++i)
	{
		ret[0] += theProblemX->soluce[mesh->elem[j + i]] * phi[i];
		ret[1] += theProblemY->soluce[mesh->elem[j + i]] * phi[i];
	}
}

#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, femDiffusionProblem *theProblemX, femDiffusionProblem *theProblemY, double dt, double tol, double iterMax)
{
	int n = myGrains->n;
	int i, iter = 0;
	double zeta;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double gamma = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];

	//
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
	//

	for (i = 0; i < n; i++)
	{
	 	int *elem = theProblemX->mesh->elem;
		femMesh *mesh = theProblemX->mesh;
		double vInc[2];
		int j;
		for (j = 0; j < 3 * mesh->nElem; j += mesh->nLocalNode)
		{
			double X[3];
			double Y[3];
			X[0] = mesh->X[elem[j]];
			Y[0] = mesh->Y[elem[j]];
			X[1] = mesh->X[elem[j + 1]];
			Y[1] = mesh->Y[elem[j + 1]];
			X[2] = mesh->X[elem[j + 2]];
			Y[2] = mesh->Y[elem[j + 2]];
			if (ptInTriangle(x[i], y[i], X[0], Y[0], X[1], Y[1], X[2], Y[2]) == 1)
			{
				double xp[2];
				xp[0] = x[i];
				xp[1] = y[i];
				vEval(X, Y, xp, vInc, theProblemX, theProblemY, j);
				break;
			}
		}
		// printf("vx %f\n", vInc[0]);
		// printf("vy %f\n", vInc[1]);

		// Faire l'isomorphisme (vx - u)
		vx[i] += (gx - (gamma * (vx[i] - vInc[0]) / m[i])) * dt;
		vy[i] += (gy - (gamma * (vy[i] - vInc[1]) / m[i])) * dt;
	}

	//
	// -2- Correction des vitesses pour tenir compte des contacts
	//
	do
	{
		zeta = femGrainsContactIterate(myGrains, dt, iter);
		iter++;
	} while ((zeta > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, zeta);

	//
	// -3- Calcul des nouvelles positions sans penetrations de points entre eux
	//
	for (i = 0; i < n; ++i)
	{
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
	}
}

/////END OF GRAINS

/////START OF ITERATIVESOLVER

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
	int i, j;
	int myRow;
	if (mySolver->iter == 0)
	{
		for (i = 0; i < nLoc; i++) {
			myRow = map[i];
			mySolver->R[myRow] -= Bloc[i];
			for (j = 0; j < nLoc; j++) {
				mySolver->R[myRow] += Aloc[i*nLoc + j] * Uloc[j];
			}
		}
	}
	///*
	for (i = 0; i < nLoc; i++) {
		myRow = map[i];
		for (j = 0; j < nLoc; j++) {
			int myCol = map[j];
			mySolver->S[myRow] += Aloc[i*nLoc + j] * mySolver->D[myCol];
		}
	}
	//*/
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue)
{
	mySolver->R[myNode] = 0.0;//myValue; car pas d'autres possibilites ?
							  //mySolver->D[myNode] = 0.0;
	mySolver->S[myNode] = 0.0;
	//mySolver->X[myNode] = 0.0;
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
	double alpha, beta;
	mySolver->iter++;
	double error = 0.0; int i;
	for (i = 0; i < mySolver->size; i++)
	{
		error += mySolver->R[i] * mySolver->R[i];
	}
	if (mySolver->iter == 1)
	{
		for (i = 0; i < mySolver->size; i++)
		{
			mySolver->X[i] = 0;
			mySolver->D[i] = mySolver->R[i];
		}
	}
	else
	{
		alpha = 0.0;
		beta = 0.0;
		for (i = 0; i < mySolver->size; i++)
		{
			alpha += (mySolver->S[i] * mySolver->R[i]);
		}
		alpha = -error / alpha;
		for (i = 0; i < mySolver->size; i++)
		{
			mySolver->R[i] += alpha * mySolver->S[i];
			beta += mySolver->R[i] * mySolver->R[i];
		}
		beta = beta / error;
		for (i = 0; i < mySolver->size; i++)
		{
			mySolver->X[i] = alpha * mySolver->D[i];
			mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i];
			mySolver->S[i] = 0.0;
		}
	}
	mySolver->error = sqrt(error);
	return(mySolver->X);
}

femMesh *MeshWithGlobalAccess;

int cmpY(const void * x, const void * y)
{
	if (MeshWithGlobalAccess->Y[*(int*)x] < MeshWithGlobalAccess->Y[*(int*)y])
	{
		return 1;
	}
	if (MeshWithGlobalAccess->Y[*(int*)x] > MeshWithGlobalAccess->Y[*(int*)y])
	{
		return -1;
	}
	return 0;
}

int cmpX(const void * x, const void * y)
{
	if (MeshWithGlobalAccess->X[*(int*)x] < MeshWithGlobalAccess->X[*(int*)y])
	{
		return 1;
	}
	if (MeshWithGlobalAccess->X[*(int*)x] > MeshWithGlobalAccess->X[*(int*)y])
	{
		return -1;
	}
	return 0;
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
	int i;
	int *tab;
	switch (renumType) {
	case FEM_NO:
		for (i = 0; i < theProblem->mesh->nNode; i++)
			theProblem->number[i] = i;
		break;
	case FEM_XNUM:
		tab = malloc(sizeof(int)*theProblem->mesh->nNode);
		for (i = 0; i < theProblem->mesh->nNode; i++)
		{
			tab[i] = i;
		}
		MeshWithGlobalAccess = theProblem->mesh;
		qsort(tab, theProblem->mesh->nNode, sizeof(int), cmpX);
		for (i = 0; i < theProblem->mesh->nNode; i++)
		{
			theProblem->number[tab[i]] = i;
		}
		free(tab);
		break;
	case FEM_YNUM:
		tab = malloc(sizeof(int)*theProblem->mesh->nNode);
		for (i = 0; i < theProblem->mesh->nNode; i++)
		{
			tab[i] = i;
		}
		MeshWithGlobalAccess = theProblem->mesh;
		qsort(tab, theProblem->mesh->nNode, sizeof(int), cmpY);
		for (i = 0; i < theProblem->mesh->nNode; i++)
		{
			theProblem->number[tab[i]] = i;
		}
		free(tab);
		break;
	default: Error("Unexpected renumbering option");
	}
}


int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
	int myBand = 0;
	int nodes[4];
	femMesh *theMesh = theProblem->mesh;
	int i, j;
	for (i = 0; i < theMesh->nElem; i++)
	{
		for (j = 0; j < theMesh->nLocalNode; j++)
		{
			nodes[j] = theProblem->number[theMesh->elem[i*theMesh->nLocalNode + j]];
		}
		int max = 0;
		int min = INT_MAX;
		for (j = 0; j < theMesh->nLocalNode; j++)
		{
			max = (int)fmax(max, nodes[j]);
			min = (int)fmin(min, nodes[j]);
		}
		myBand = (int)fmax(myBand, (max - min));
	}
	return(myBand + 1);
}

/*
void femDiffusionUpdate(femDiffusionProblem *theProblem, femGrains *g)
{
	int i,j;
    for (i = 0; i < theProblem->size; i++)
	{
		theProblem->soluce[i] = 0.0;
		theProblem->solver->solver->B[i] = 0.0;
		for (j = 0; j < theProblem->size; j++)
		{
			theProblem->solver->solver->A[i][j] = 0.0;
		}
	}
}
*/

#endif
