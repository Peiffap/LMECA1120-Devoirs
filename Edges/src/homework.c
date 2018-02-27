#include "fem.h"
#include <time.h>


# ifndef NOEXPAND

void edgesExpand(femEdges *theEdges)
{
	femMesh *mesh = theEdges->mesh;
	int i;
	for (i = 0; i < ((theEdges->nEdge) - 2); i+=3)
	{
		theEdges->edges[i].elem[0] = i / 3;
		theEdges->edges[i].elem[1] = -1;
		theEdges->edges[i+1].elem[0] = i / 3;
		theEdges->edges[i+1].elem[1] = -1;
		theEdges->edges[i+2].elem[0] = i / 3;
		theEdges->edges[i+2].elem[1] = -1;
		theEdges->edges[i].node[0] = *(mesh->elem + i);
		theEdges->edges[i].node[1] = *(mesh->elem + i + 1);
		theEdges->edges[i+1].node[0] = *(mesh->elem + i + 1);
		theEdges->edges[i+1].node[1] = *(mesh->elem + i + 2);
		theEdges->edges[i+2].node[0] = *(mesh->elem + i + 2);
		theEdges->edges[i+2].node[1] = *(mesh->elem + i);
	}

	//SECTION TEST
	/*
	for (i = 0; i < theEdges->nEdge; i++)
	{
		printf("%d :	%d	%d :	%d	%d\n", i, theEdges->edges[i].node[0], theEdges->edges[i].node[1], theEdges->edges[i].elem[0], theEdges->edges[i].elem[1]);
	}
	*/
}

# endif
# ifndef NOSORT

void edgesSort(femEdges *theEdges)
{
	qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

int edgesCompare(const void* e0, const void *e1)
{
	if (((femEdge *) e0)->node[0] == ((femEdge *) e1)->node[1] && ((femEdge *) e0)->node[1] == ((femEdge *) e1)->node[0])
	{
		return 0;
	}
	else
	{
		double d = fmin(((femEdge *) e0)->node[0], ((femEdge *) e0)->node[1]);
		double f = fmin(((femEdge *) e1)->node[0], ((femEdge *) e1)->node[1]);
		if (d > f)
		{
			return -1;
		}
		else if(d == f)
		{
			double d = fmax(((femEdge *) e0)->node[0], ((femEdge *) e0)->node[1]);
			double f = fmax(((femEdge *) e1)->node[0], ((femEdge *) e1)->node[1]);
			if (d > f)
			{
				return -1;
			}
			return 1;
		}
		return 1;
	}
}

# endif
# ifndef NOSHRINK

void edgesShrink(femEdges *theEdges)
{
	// Ici commence votre contribution a faire :-)

	// int n = 0;          // Nouveau nombre total de segments
	int nBoundary = 0;  // Nombre de segments frontieres

	/*
	int limit = 0;
	int k;
	for (k = 0; k < (theEdges->nEdge) - 1; k++)
	{
		if (theEdges->edges[k].node[0] == theEdges->edges[k + 1].node[1] && theEdges->edges[k].node[1] == theEdges->edges[k + 1].node[0])
		{
			limit++;
		}
	}
	n = theEdges->nEdge - limit;
	*/

	int i;

	/*
	for (i = 0; i < theEdges->nEdge; i++)//TEST
	{
		printf("%d :	%d	%d :	%d	%d\n", i, theEdges->edges[i].node[0], theEdges->edges[i].node[1], theEdges->edges[i].elem[0], theEdges->edges[i].elem[1]);
	}
	printf("\n");
	clock_t startThomas = clock();
	for (i = 0; i < (theEdges->nEdge) - 1; i++)
	{
		if (theEdges->edges[i].node[0] == theEdges->edges[i + 1].node[1] && theEdges->edges[i].node[1] == theEdges->edges[i + 1].node[0])
		{
			theEdges->edges[i].elem[1] += theEdges->edges[i + 1].elem[0] + 1;
			int j;
			for (j = i + 1; j < (theEdges->nEdge) - 1; j++)
			{
				theEdges->edges[j].elem[0] = theEdges->edges[j + 1].elem[0];
				theEdges->edges[j].elem[1] = theEdges->edges[j + 1].elem[1];
				theEdges->edges[j].node[0] = theEdges->edges[j + 1].node[0];
				theEdges->edges[j].node[1] = theEdges->edges[j + 1].node[1];
			}
		}
	}
	clock_t endThomas = clock();
	float secondsThomas = (float)(endThomas - startThomas) / CLOCKS_PER_SEC;
	printf("Temps : %f \n", secondsThomas);

	for (i = 0; i < theEdges->nEdge; i++)//TEST
	{
		printf("%d :	%d	%d :	%d	%d\n", i, theEdges->edges[i].node[0], theEdges->edges[i].node[1], theEdges->edges[i].elem[0], theEdges->edges[i].elem[1]);
	}
	*/

	// clock_t startGilles = clock();
	int j = 0;
	for (i = 0; i < (theEdges->nEdge) - 1; i++)
	{
		theEdges->edges[j].elem[0] = theEdges->edges[i].elem[0];
		theEdges->edges[j].elem[1] = theEdges->edges[i].elem[1];
		theEdges->edges[j].node[0] = theEdges->edges[i].node[0];
		theEdges->edges[j].node[1] = theEdges->edges[i].node[1];
		if (theEdges->edges[i].node[0] == theEdges->edges[i + 1].node[1] && theEdges->edges[i].node[1] == theEdges->edges[i + 1].node[0])
		{
			theEdges->edges[j].elem[1] += theEdges->edges[i + 1].elem[0] + 1;
			i++;
		}
		j++;
	}
	theEdges->edges[j].elem[0] = theEdges->edges[theEdges->nEdge - 1].elem[0];
	theEdges->edges[j].elem[1] = theEdges->edges[theEdges->nEdge - 1].elem[1];
	theEdges->edges[j].node[0] = theEdges->edges[theEdges->nEdge - 1].node[0];
	theEdges->edges[j].node[1] = theEdges->edges[theEdges->nEdge - 1].node[1];
	/*
	clock_t endGilles = clock();
	float secondsGilles = (float)(endGilles - startGilles) / CLOCKS_PER_SEC;
	printf("Temps %f \n", secondsGilles);
	*/

	int n = j + 1;

	for (i = 0; i < n; i++)
	{
		if (theEdges->edges[i].elem[1] == -1)
		{
			nBoundary++;
		}
	}

	// Ici, finit votre contribution

	// Reallocation du tableau des edges

	theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
	theEdges->nEdge = n;
	theEdges->nBoundary = nBoundary;

	//JUST SOME TEST TO BE SURE
	/*
	for (i = 0; i < theEdges->nEdge; i++)
	{
		printf("%d :	%d	%d :	%d	%d\n", i, theEdges->edges[i].node[0], theEdges->edges[i].node[1], theEdges->edges[i].elem[0], theEdges->edges[i].elem[1]);
	}
	*/

}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{
	double L = 0;
	int i;
	for (i = 0; i < theEdges->nEdge; i++)
	{
		if (theEdges->edges[i].elem[1] == -1)
		{
			int n0 = theEdges->edges[i].node[0];
			int n1 = theEdges->edges[i].node[1];
			double X0 = theEdges->mesh->X[n0];
			double Y0 = theEdges->mesh->Y[n0];
			double X1 = theEdges->mesh->X[n1];
			double Y1 = theEdges->mesh->Y[n1];
			L += sqrt(pow(X0 - X1, 2) + pow(Y0 - Y1, 2));
		}
	}
	return L;
}

# endif
