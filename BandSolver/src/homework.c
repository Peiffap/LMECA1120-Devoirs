#include"fem.h"
#include <limits.h>

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
        case FEM_NO :
            for (i = 0; i < theProblem->mesh->nNode; i++)
                theProblem->number[i] = i;
            break;
//
// A modifier :-)
// debut
//
        case FEM_XNUM :
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
				//printf("%d:\t%d\n", i, tab[i]);
			}
			free(tab);
			break;
        case FEM_YNUM :
			/* Avons eu besoin de l'aide du drive afin de trouver la maniere generale de faire, nos multiples essais etant infructueux
			for (i = 0; i < theProblem->mesh->nNode; i++)
			{
				ind = i;
				//printf("%d : \tx : %f\t y : %f\n", i, theProblem->mesh->X[i], theProblem->mesh->Y[i]);
				for (j = 0; j < theProblem->mesh->nNode; j++)
				{
					for (k = 0; k < i && isAlreadyIn == 0; k++)
					{
						if (theProblem->number[k] == j)
						{
							isAlreadyIn = 1;
							printf("j = %d\t number[%d] = %d\t bool = %d\n", j, k, theProblem->number[k], isAlreadyIn);
						}
					}
					if (theProblem->mesh->Y[j] < theProblem->mesh->Y[ind] && isAlreadyIn == 0)
					{
						ind = j;
					}
					printf("ind = %d\n", ind);
					isAlreadyIn = 0;
				}
				theProblem->number[i] = ind;
				printf("number[%d] : \t%d\n", i, ind);
			}
			*/
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
				//printf("%d:\t%d\n", i, tab[i]);
			}
			free(tab);
            break;
//
// end
//

        default : Error("Unexpected renumbering option"); }
}


int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
	int myBand = 0;
	int nodes[4];
    femMesh *theMesh = theProblem->mesh;
	int i, j;
	for (i = 0; i < theMesh->nElem; i++)
	{
		//printf("%d : \t%d\t%d\t%d\n", i, theMesh->elem[i%theMesh->nElem], theMesh->elem[i%theMesh->nElem + 1], theMesh->elem[i%theMesh->nElem + 2]);
		for (j = 0; j < theMesh->nLocalNode; j++)
		{
			nodes[j] = theProblem->number[theMesh->elem[i*theMesh->nLocalNode + j]];
			//printf("node : %d\n", nodes[j]);
		}
		int max = 0;
		int min = INT_MAX;
		for (j = 0; j < theMesh->nLocalNode; j++)
		{
			max = (int)fmax(max, nodes[j]);
			min = (int)fmin(min, nodes[j]);
			//printf("max : %d\t min : %d\n", max, min);
		}
		myBand = (int)fmax(myBand, (max - min));
		//printf("\tband : %d\n", myBand);
	}
    return(myBand+1);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc) //eu besoin de l'aide du drive car no idea of what we were supposed to do
{
    int i,j;
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
