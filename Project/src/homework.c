#include"fem.h"
#include"limits.h"


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

			// Param�tres supplementaires
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
		vx[i] += (gx - (gamma * vx[i] / m[i])) * dt;
		vy[i] += (gy - (gamma * vy[i] / m[i])) * dt;
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

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc) //eu besoin de l'aide du drive car no idea of what we were supposed to do
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

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}


# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
	int j;
	for (j = 0; j < theMesh->nLocalNode; j++)
	{
		int a = theMesh->elem[j + i * theMesh->nLocalNode];
		*(map + j) = a;
		*(x + j) = theMesh->X[a];
		*(y + j) = theMesh->Y[a];
	}
}

# endif
# ifndef NOPOISSONSOLVE

// Inspire du code disponible sur le drive
void femPoissonSolve(femPoissonProblem *theProblem)
{
	/*
	PRIMARY CHECK
	femDiscretePrint(theProblem->space);
	printf("edgs\n");
	femEdgesPrint(theProblem->edges);
	printf("sys 1\n");
	femFullSystemPrint(theProblem->system);
	*/
	double **A = theProblem->system->A;
	double *B = theProblem->system->B;
	int i;
	for (i = 0; i < theProblem->mesh->nElem; i++)
	{
		int *map = malloc(theProblem->mesh->nLocalNode * sizeof(int));
		double *x = malloc(theProblem->mesh->nLocalNode * sizeof(double));
		double *y = malloc(theProblem->mesh->nLocalNode * sizeof(double));
		femMeshLocal(theProblem->mesh, i, map, x, y);
		double *dphidxsi = malloc(sizeof(double) * theProblem->mesh->nLocalNode);
		double *dphideta = malloc(sizeof(double) * theProblem->mesh->nLocalNode);
		for (int ii = 0; ii < theProblem->mesh->nLocalNode; ii++)
		{
			double xsi = theProblem->rule->xsi[ii];
			double eta = theProblem->rule->eta[ii];
			double weight = theProblem->rule->weight[ii];
			theProblem->space->dphi2dx(xsi, eta, dphidxsi, dphideta);
			double dxdxsi = 0;
			double dxdeta = 0;
			double dydxsi = 0;
			double dydeta = 0;
			int a;
			for (a = 0; a < theProblem->mesh->nLocalNode; a++)
			{
				dxdxsi += x[a] * dphidxsi[a];
				dxdeta += x[a] * dphideta[a];
				dydxsi += y[a] * dphidxsi[a];
				dydeta += y[a] * dphideta[a];
			}
			double *phi = malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			theProblem->space->phi2(xsi, eta, phi);
			double J = fabs(dxdxsi*dydeta - dydxsi * dxdeta);
			double *dphidx = malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			double *dphidy = malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			for (a = 0; a < theProblem->mesh->nLocalNode; a++)
			{
				dphidx[a] = 1 / J * (dphidxsi[a] * dydeta - dphideta[a] * dydxsi);
				dphidy[a] = 1 / J * (dphideta[a] * dxdxsi - dphidxsi[a] * dxdeta);
			}
			int j, k;
			for (j = 0; j < theProblem->mesh->nLocalNode; j++)
			{
				for (k = 0; k < theProblem->mesh->nLocalNode; k++)
				{
					A[map[j]][map[k]] += weight * J * (dphidx[j] * dphidx[k] + dphidy[j] * dphidy[k]);
				}
				B[map[j]] += phi[j] * J * weight;
			}
			free(phi);
			free(dphidx);
			free(dphidy);
		}
		free(map);
		free(x);
		free(y);
		free(dphidxsi);
		free(dphideta);
	}
	// printf("sys 2\n");
	// femFullSystemPrint(theProblem->system);
	for (i = 0; i < theProblem->edges->nEdge; i++)
	{
		if (theProblem->edges->edges[i].elem[1] == -1)
		{
			femFullSystemConstrain(theProblem->system, theProblem->edges->edges[i].node[0], 0);
			femFullSystemConstrain(theProblem->system, theProblem->edges->edges[i].node[1], 0);
			// printf("n1 : %d | n2 : %d\n", theProblem->edges->edges[i].node[0], theProblem->edges->edges[i].node[1]);
		}
	}
	// printf("sys 3\n");
	// femFullSystemPrint(theProblem->system);
	femFullSystemEliminate(theProblem->system);
	// printf("resolu : \n");
	// femFullSystemPrint(theProblem->system);
}

#endif
