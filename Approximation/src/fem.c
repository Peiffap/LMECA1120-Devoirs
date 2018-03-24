/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussQuad4Xsi[4]     = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Eta[4]     = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Weight[4]  = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000 };
static const double _gaussTri3Xsi[3]      = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3]      = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3]   = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };
static const double _gaussTri12Xsi[12]    = { 0.249286745170910,0.249286745170910,0.501426509658179, 
										                          0.063089014491502,0.063089014491502,0.873821971016996, 
										                          0.310352451033785,0.636502499121399,0.053145049844816, 
										                          0.310352451033785,0.636502499121399,0.053145049844816 };
static const double _gaussTri12Eta[12]    = { 0.249286745170910,0.501426509658179,0.249286745170910,
										                          0.063089014491502,0.873821971016996,0.063089014491502,
										                          0.636502499121399,0.053145049844816,0.310352451033785,
										                          0.053145049844816,0.310352451033785,0.636502499121399 };
static const double _gaussTri12Weight[12] = { 0.058393137863189,0.058393137863189,0.058393137863189,
										                          0.025422453185104,0.025422453185104,0.025422453185104,
										                          0.041425537809187,0.041425537809187,0.041425537809187,
										                          0.041425537809187,0.041425537809187,0.041425537809187 };

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
    else if (type == FEM_TRIANGLE && n == 12) {
        theRule->n      = 12;
        theRule->xsi    = _gaussTri12Xsi;
        theRule->eta    = _gaussTri12Eta;
        theRule->weight = _gaussTri12Weight; }
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


void _p3c0_x(double *xsi, double *eta) 
{
    xsi[0] = 0.0;     eta[0] = 0.0;
    xsi[1] = 1.0;     eta[1] = 0.0;
    xsi[2] = 0.0;     eta[2] = 1.0;
    xsi[3] = 1.0/3.0; eta[3] = 0.0;
    xsi[4] = 2.0/3.0; eta[4] = 0.0;
    xsi[5] = 2.0/3.0; eta[5] = 1.0/3.0;
    xsi[6] = 1.0/3.0; eta[6] = 2.0/3.0;
    xsi[7] = 0.0;     eta[7] = 2.0/3.0;
    xsi[8] = 0.0;     eta[8] = 1.0/3.0;
    xsi[9] = 1.0/3.0; eta[9] = 1.0/3.0;
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
    else if (type == FEM_TRIANGLE && n == 10) {
        theSpace->n       = 10;
        theSpace->x2      = _p3c0_x;
        theSpace->phi2    = femApproxPhi;
        theSpace->dphi2dx = femApproxDphi; }
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
    double xsi[10], eta[10], phi[10], dphidxsi[10], dphideta[10];
    
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

void femMeshClean(femMesh *theMesh)
{
    int i,j,*elem;
    
     
    int *check = malloc(sizeof(int)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i)
      	check[i] = 0;
    int *map = malloc(sizeof(int)*theMesh->nNode);
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        for (j = 0; j < theMesh->nLocalNode; ++j) {
            check[elem[j]] = 1; }}
    int iGlo = 0;
    for (i = 0; i < theMesh->nNode; ++i)  {
      	if (check[i] != 0) {
            map[i] = iGlo;
            theMesh->X[iGlo] = theMesh->X[i];
            theMesh->Y[iGlo] = theMesh->Y[i];
            iGlo++; }}
    theMesh->nNode = iGlo;
    
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        for (j = 0; j < theMesh->nLocalNode; ++j) {
        	elem[j] = map[elem[j]]; }}
            
    free(check);
    free(map);
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

double femApproxStommel(double x, double y)
{
    //
    // Solution analytique de Stommel dans un carre [0,1]x[0,1]
    // Modelisation de l'elevation de l'ocean Atlantique dans un carre adimensionnel
    // Ce modele que l'on attribue generalement au grand oceanographe Henry M.
    // Stommel (1920-1992), est considere comme le premier modele qualitativement correct du Gulf Stream
    //
    
    const double tau0 = 0.1;
    const double L = 1e6;
    const double gamm = 1e-6;
    const double rho = 1000;
    const double delta = 1;
    const double g = 9.81;
    const double h = 1000;
    const double f0 = 1e-4;
    const double beta = 2e-11;
    
    double Y = y - 0.5;
    double epsilon = gamm / (L * beta);
    double Z1 = (-1 + sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double Z2 = (-1 - sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double D = ((exp(Z2) - 1) * Z1 + (1 - exp(Z1)) * Z2) / (exp(Z1) - exp(Z2));
    double f1 = M_PI / D * (1 + ((exp(Z2) - 1) * exp(x * Z1) + (1 - exp(Z1)) * exp(x * Z2)) / (exp(Z1) - exp(Z2)));
    double f2 = 1 / D* (((exp(Z2) - 1) * Z1 * exp(x * Z1) + (1 - exp(Z1)) * Z2 * exp(x * Z2)) / (exp(Z1) - exp(Z2)));
    double eta = D * tau0 * f0 * L / (M_PI * gamm * rho * delta * g * h) *
    (  -gamm / (f0 * delta * M_PI) * f2 * sin(M_PI * Y) 
     + 1 / M_PI * f1 * (cos(M_PI * Y) * (1 + beta * Y) 
                        - beta / M_PI * sin(M_PI * Y)));
    return eta;
    
}


femApproxProblem *femApproxCreate(const char *filename)
{
    femApproxProblem *theProblem = malloc(sizeof(femApproxProblem));
    theProblem->mesh  = femMeshRead(filename);   
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    
    // Pour une approximation P1-CO, remplacer les 3 lignes 
    // qui suivent par la version en commentaire :-)
    // 
    //  theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
    //  theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE);   
    //  int n = theProblem->mesh->nNode ;
    //
    
    theProblem->space = femDiscreteCreate(10,FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(12,FEM_TRIANGLE);   
    int n = 2*(theProblem->edges->nEdge) + theProblem->mesh->nNode + theProblem->mesh->nElem; 
    theProblem->system = femFullSystemCreate(n);
    theProblem->rule = femIntegrationCreate(12,FEM_TRIANGLE); 
    theProblem->minValue = 0.0;
    theProblem->maxValue = 0.0;
    return theProblem;
}


void femApproxFree(femApproxProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}
  
  
double femApproxEval(const femApproxProblem *theProblem, int iElem, double xsi, double eta)
{
    double phi[20];
    int map[20];
    
    femDiscretePhi2(theProblem->space,xsi,eta,phi);
    femApproxLocal(theProblem,iElem,map);     
    double *U = theProblem->system->B;
    
    double u = 0.0;    
    for (int i = 0; i < theProblem->space->n; i++)    
        u += U[map[i]]*phi[i];        
    return u;   
}  


femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;  
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    mySystem->A = malloc(sizeof(double*) * size); 
    mySystem->B = elem;
    mySystem->A[0] = elem + size;  
    mySystem->size = size;
    for (i=1 ; i < size ; i++) 
        mySystem->A[i] = mySystem->A[i-1] + size;
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

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
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

void  femFullSystemConstrain(femFullSystem *mySystem, 
                             int myNode, double myValue) 
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

