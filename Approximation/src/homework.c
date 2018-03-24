
#include "fem.h"


  

# ifndef NOAPPROXPHI

void femApproxPhi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) *  9.0/2.0; 
    phi[1] = 0.0;
    phi[2] = 0.0;  
    phi[3] = 0.0;
    phi[4] = 0.0;
    phi[5] = 0.0;
    phi[6] = 0.0;  
    phi[7] = 0.0;
    phi[8] = 0.0;
    phi[9] = 0.0;   
    
    // La premiere fonction vous est donnee gracieusement....
    // Il vous reste les 9 autres a obtenir !
    // A modifier :-)

}

# endif
# ifndef NOAPPROXDPHI

void femApproxDphi(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    int i;
    for (i=0;i < 10; i++) {
      dphidxsi[i] = 0;
      dphideta[i] = 0;  }
    
     
    // A faire a titre de bonus, car on n'a pas vraiment besoin de ces derivees...
    // Il n'a qu'un unique point a gagner avec ceci
}

# endif
# ifndef NOAPPROXLOCAL

void femApproxLocal(const femApproxProblem *theProblem, const int iElem, int *map)
{
    femMesh *theMesh = theProblem->mesh;
    
    
    // Pour une approximation P1-CO, remplacer ce qui suit par :
    //
    //  for (int j=0; j < 3; j++) {  
    //      map[j] = theMesh->elem[iElem*3 + j]; }
    // 
    if (iElem > 1)  Error("It only works with meshes one.txt and two.txt ");
    int mapElem[2][10] = {{0,1,2,3,4,5,6,7,8,9},{1,10,2,11,12,13,14,6,5,15}};
    for (int j=0; j < 10; j++) { 
        map[j] = mapElem[iElem][j]; }  
        
    // A generaliser pour un maillage quelconque avec une interpolation P3-C0
    // Avec la version actuelle, le code ne fonctionne qu'avec les maillages one.txt et two.text
    //
  
       
}

# endif
# ifndef NOAPPROXSOLVE

void femApproxSolve(femApproxProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    double x[3],y[3],phi[10];
    int iElem,iInteg,iEdge,i,j,map[10];

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femApproxLocal(theProblem,iElem,map); 
        for (j=0; j < 3; ++j) {
            int *elem = &theMesh->elem[iElem*3];
            x[j]   = theMesh->X[elem[j]];
            y[j]   = theMesh->Y[elem[j]]; } 
          
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
   
     // A completer et a modifier ce qui suit 
     // pour obtenir l'approximation de la fonction de Stommel !
     // Ici on construit une matrice diagonale et 
     // et un menbre de droite nulle : la solution sera donc nulle :-)
   
            for (i = 0; i < theSpace->n; i++) { 
                theSystem->A[map[i]][map[i]] = 1.0; }                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                theSystem->B[map[i]] = 0.0; }}} 
                
    femFullSystemEliminate(theSystem);
}


# endif

