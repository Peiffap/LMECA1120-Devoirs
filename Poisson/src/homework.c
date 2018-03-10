
#include"fem.h"



# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
	femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
	theProblem->mesh  = femMeshRead(filename);
	theProblem->edges = femEdgesCreate(theProblem->mesh);
	if (theProblem->mesh->nLocalNode == 4) {
		theProblem->space = femDiscreteCreate(4,FEM_QUAD);
		theProblem->rule = femIntegrationCreate(4,FEM_QUAD);
	}
	else if (theProblem->mesh->nLocalNode == 3)
	{
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


	void femPoissonSolve(femPoissonProblem *theProblem)
	{
		//PRIMARY CHECK
		femDiscretePrint(theProblem->space);
		printf("edgs\n");
		femEdgesPrint(theProblem->edges);
		printf("sys\n");
		femFullSystemPrint(theProblem->system);
		//
		double **A = theProblem->system->A;
		double *B = theProblem->system->B;
		// /!\ uniquement pour element parent, à fix
		double *xi = (double *) malloc(theProblem->mesh->nLocalNode * sizeof(double));
		double *eta = (double *) malloc(theProblem->mesh->nLocalNode * sizeof(double));
		theProblem->space->x2(xi, eta);
		int i;
		for (i = 0; i < theProblem->mesh->nElem; i++)
		{
			int *map = (int *) malloc(theProblem->mesh->nLocalNode * sizeof(int));
			double *x = (double *) malloc(theProblem->mesh->nLocalNode * sizeof(double));
			double *y = (double *) malloc(theProblem->mesh->nLocalNode * sizeof(double));
			femMeshLocal(theProblem->mesh, i, map, x, y);
			double *dphidx = (double *) malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			double *dphidy = (double *) malloc(sizeof(double) * theProblem->mesh->nLocalNode);
			theProblem->space->dphi2dx(0, 0, dphidx, dphidy);
			int j,k;
			for (j = 0; j < theProblem->mesh->nLocalNode; j++)
			{
				for (k = 0; k < theProblem->mesh->nLocalNode; k++)
				{
					A[map[j]][map[k]] += dphidx[j] * dphidx[k] + dphidy[j] * dphidy[k];
				}
			}
			free(map);
			free(x);
			free(y);
			free(dphidx);
			free(dphidy);
		}
		/* A changer
		for (i = 0; i < theProblem->edges->nEdge; i++)
		{
				if(theProblem->edges->edges[i].elem[1]==-1)
			{
				femFullSystemConstrain(theProblem->system, i, 1);
			}
		}
		*/
		free(xi);
		free(eta);
		printf("sys\n");
		femFullSystemPrint(theProblem->system);
		femFullSystemEliminate(theProblem->system);
		printf("sysfinal\n");
		femFullSystemPrint(theProblem->system);

	}

	# endif
