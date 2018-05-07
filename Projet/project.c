#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define FALSE 0 
#define TRUE  1

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    double radiusIn;
    double radiusOut;
    double gravity[2];
    double gamma;
    double *x;
    double *y;
    double *vx;
    double *vy;
    double *r;
    double *m;
    double *dvBoundary;
    double *dvContacts;
} femGrains;

typedef struct 
{
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femIterativeSolver *solver;
    int size;
    int *number;
    double *soluce;
} femDiffusionProblem;


