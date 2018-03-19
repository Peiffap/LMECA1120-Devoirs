/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2017 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
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

void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}


femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
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
    else if (!strncmp(str,"Number of quads",15))  { 
        printf("%s \n",str);
        ErrorScan(sscanf(str,"Number of quads %d \n", &theMesh->nElem));  
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            ErrorScan(fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3])); }}
  
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
    
    if (theMesh->nLocalNode == 4) {
        fprintf(file, "Number of quads %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[0],elem[1],elem[2],elem[3]);   }}
    else if (theMesh->nLocalNode == 3) {
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
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;          
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }
      
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
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

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
                        return  0;
}

femSolver *femSolverFullCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *)femFullSystemCreate(size);
    return(mySolver);
}

femSolver *femSolverBandCreate(int size, int band)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femSolver *femSolverIterativeCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *)femIterativeSolverCreate(size);
    return(mySolver);
}

void femSolverFree(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    free(mySolver);
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}


double femSolverGet(femSolver *mySolver,int i,int j)
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrint((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

  
void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver,Aloc,Bloc,Uloc,map,nLoc); break;
        default : Error("Unexpected solver type"); }
}
    
void femSolverConstrain(femSolver *mySolver, int myNode, double myValue)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_ITER : femIterativeSolverConstrain((femIterativeSolver *)mySolver->solver,myNode,myValue); break;
        default : Error("Unexpected solver type"); }
}
  
double *femSolverEliminate(femSolver *mySolver)
{
    double *soluce;
    switch (mySolver->type) {
        case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(soluce);
}



int femSolverConverged(femSolver *mySolver)
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}


femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    theSystem->A = malloc(sizeof(double*) * size); 
    theSystem->B = malloc(sizeof(double) * size * (size+1));
    theSystem->A[0] = theSystem->B + size;  
    theSystem->size = size;
    int i;
    for (i=1 ; i < size ; i++) 
        theSystem->A[i] = theSystem->A[i-1] + size;
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}


void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

void femFullSystemPrintInfos(femFullSystem *mySystem)
{
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));     
}

void  femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    mySystem->B[map[i]] += Bloc[i]; }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]); 
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}


 
void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size, band, ifirst, iend;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    ifirst = fmax(0,myNode - band + 1);
    iend   = myNode;
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }

    ifirst = myNode+1;
    iend = fmin(myNode + band,size);
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[myNode][i];
        A[myNode][i] = 0; }
        
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}
 
double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}



femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*4);      
    mySolver->D = mySolver->R + size;       
    mySolver->S = mySolver->R + size*2;       
    mySolver->X = mySolver->R + size*3;       
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++) 
        mySolver->R[i] = 0;        
}
 
void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

femDiffusionProblem *femDiffusionCreate(const char *filename, femSolverType solverType, femRenumType renumType)
{
    int i,band;
    femDiffusionProblem *theProblem = malloc(sizeof(femDiffusionProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->size = theProblem->mesh->nNode;
    theProblem->number  = malloc(sizeof(int)*theProblem->size);
    femDiffusionRenumber(theProblem,renumType);
         
    switch (solverType) {
        case FEM_FULL : 
                theProblem->solver = femSolverFullCreate(theProblem->size); break;
        case FEM_BAND : 
                band = femDiffusionComputeBand(theProblem);
                theProblem->solver = femSolverBandCreate(theProblem->size,band); break;
        case FEM_ITER : 
               theProblem->solver = femSolverIterativeCreate(theProblem->size); break;
        default : Error("Unexpected solver option"); }
        
    theProblem->soluce = malloc(sizeof(double)*theProblem->size);
    for (i = 0; i < theProblem->size; i++)      
        theProblem->soluce[i] = 0;
 
    return theProblem;
}

void femDiffusionFree(femDiffusionProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    femSolverFree(theProblem->solver);
    free(theProblem->number);
    free(theProblem->soluce);
    free(theProblem);
}
    

void femDiffusionMeshLocal(const femDiffusionProblem *theProblem, const int iElem, int *map, double *x, double *y, double *u)
{
    femMesh *theMesh = theProblem->mesh;
    int j,nLocal = theMesh->nLocalNode;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]]; 
        u[j]   = theProblem->soluce[map[j]];
        map[j] = theProblem->number[map[j]];
        }   
}

void femDiffusionCompute(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femSolver *theSolver = theProblem->solver;
    femEdges *theEdges = theProblem->edges;
    int *number = theProblem->number;
       
    if (theSpace->n > 4) Error("Unexpected discrete space size !"); 
    
    double Xloc[4],Yloc[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4],Aloc[16],Bloc[4],Uloc[4];
    int iEdge,iElem,iInteg,i,j,map[4];
   
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(theProblem,iElem,map,Xloc,Yloc,Uloc);  
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {    
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    Aloc[i*(theSpace->n)+j] += (dphidx[i] * dphidx[j] 
                                            + dphidy[i] * dphidy[j]) * jac * weight; }}                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                Bloc[i] += phi[i] * jac *weight; }}
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); } 

    for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) {      
        if (theEdges->edges[iEdge].elem[1] < 0) {       
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[0]],0.0); 
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[1]],0.0); }}
  
    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->mesh->nNode; i++)
        theProblem->soluce[i] += soluce[number[i]];
  
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
