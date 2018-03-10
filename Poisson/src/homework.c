#include"fem.h"

// Thomas

# ifndef NOPOISSONCREATE

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
		int *map = (int *) malloc(theProblem->mesh->nLocalNode * sizeof(int));
		double *x = (double *)malloc(theProblem->mesh->nLocalNode * sizeof(double));
		double *y = (double *)malloc(theProblem->mesh->nLocalNode * sizeof(double));
		femMeshLocal(theProblem->mesh, i, map, x, y);
		double *dphidxsi = (double *) malloc(sizeof(double) * theProblem->mesh->nLocalNode);
		double *dphideta = (double *) malloc(sizeof(double) * theProblem->mesh->nLocalNode);
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
			double *phi = (double *)malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			theProblem->space->phi2(xsi, eta, phi);
			double J = fabs(dxdxsi*dydeta - dydxsi * dxdeta);
			double *dphidx = (double *)malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			double *dphidy = (double *)malloc(sizeof(double) * theProblem->mesh->nLocalNode);
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

# endif
