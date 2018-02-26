/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}


femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }
    
    char str[256]; if (fgets(str, sizeof(str), file) == NULL) Error("Corrupted mesh file !");

    if (!strncmp(str,"Number of triangles",19))  { 
        ErrorScan(sscanf(str,"Number of triangles %d \n", &theMesh->nElem));
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            ErrorScan(fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2])); }}
  
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->elem);
    free(theMesh);
}

void femMeshWrite(const femMesh *theMesh, const char *filename)
{
    int i,*elem;
    
    FILE* file = fopen(filename,"w");
    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]); }
    
    if (theMesh->nLocalNode == 3) {
        fprintf(file, "Number of triangles %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fprintf(file,"%6d : %6d %6d %6d \n", i,elem[0],elem[1],elem[2]);   }}
    
    fclose(file);
}
                                                                                                          
femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int n    = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    return theEdges;
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;    
    for (i = 0; i < theEdges->nEdge; ++i) {
        printf("%6d : %4d %4d : %4d %4d \n",i,
               theEdges->edges[i].node[0],theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0],theEdges->edges[i].elem[1]); }
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
