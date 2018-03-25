
#include "fem.h"




# ifndef NOAPPROXPHI

void femApproxPhi(double xsi, double eta, double *phi)
{
	double third = 1.0 / 3.0;
	double tthird = 2.0 * third;

    phi[0] = (1.0 - xsi - eta) * (third - xsi - eta) * (tthird - xsi - eta) *  4.5;
    phi[1] = xsi * (xsi - third) * (xsi - tthird) * 4.5;
    phi[2] = eta * (eta - third) * (eta - tthird) * 4.5;
    phi[3] = (1.0 - xsi - eta) * (tthird - xsi - eta) * xsi * 13.5;
    phi[4] = (1.0 - xsi - eta) * (xsi - third) * xsi * 13.5;
    phi[5] = (xsi - third) * xsi * eta * 13.5;
    phi[6] = (eta - third) * xsi * eta * 13.5;
    phi[7] = (1.0 - xsi - eta) * (eta - third) * eta * 13.5;
    phi[8] = (1.0 - xsi - eta) * (tthird - xsi - eta) * eta * 13.5;
    phi[9] = (1.0 - xsi - eta) * xsi * eta * 27.0;

    // La premiere fonction vous est donnee gracieusement....
    // Il vous reste les 9 autres a obtenir !
    // A modifier :-)

}

# endif
# ifndef NOAPPROXDPHI

void femApproxDphi(double xsi, double eta, double *dphidxsi, double *dphideta)
{
	double eta2 = pow(eta,2);
	double xsi2 = pow(xsi,2);

    dphidxsi[0] = 0.5 * (-27.0 * eta2 + eta * (36.0 - 54.0 * xsi) - 27.0 * xsi2 + 36.0 * xsi - 11.0);
	dphideta[0] = 0.5 * (-27.0 * eta2 + eta * (36.0 - 54.0 * xsi) - 27.0 * xsi2 + 36.0 * xsi - 11.0);
	dphidxsi[1] = 13.5 * xsi2 - 9.0 * xsi + 1.0;
	dphideta[1] = 0.0;
	dphidxsi[2] = 0.0;
	dphideta[2] = 13.5 * eta2 - 9.0 * eta + 1.0;
	dphidxsi[3] = 4.5 * (3.0 * eta2 + eta * (12.0 * xsi - 5.0) + 9.0 * xsi2 - 10.0 * xsi + 2.0);
	dphideta[3] = 4.5 * xsi * (6.0 * eta + 6.0 * xsi - 5.0);
	dphidxsi[4] = -4.5 * (eta * (6.0 * xsi - 1.0) + 9.0 * xsi2 - 8.0 * xsi + 1.0);
	dphideta[4] = -4.5 * xsi * (3.0 * xsi - 1.0);
	dphidxsi[5] = 4.5 * eta * (6.0 * xsi - 1.0);
	dphideta[5] = 13.5 * (xsi - 1.0 / 3.0) * xsi;
	dphidxsi[6] = 13.5 * (eta - 1.0 / 3.0) * eta;
	dphideta[6] = 4.5 * (6.0 * eta - 1.0) * xsi;
	dphidxsi[7] = -4.5 * eta * (3.0 * eta - 1.0);
	dphideta[7] = -4.5 * (9.0 * eta2 + eta * (6.0 * xsi - 8.0) - xsi + 1.0);
	dphidxsi[8] = 4.5 * eta * (6.0 * eta + 6.0 * xsi - 5.0);
	dphideta[8] = 4.5 * (9.0 * eta2 + 2.0 * eta * (6.0 * xsi - 5.0) + 3.0 * xsi2 - 5.0 * xsi + 2.0);
	dphidxsi[9] = -27.0 * eta * (eta + 2.0 * xsi - 1.0);
	dphideta[9] = -27.0 * xsi * (2.0 * eta + xsi - 1.0);

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
    // Avec la version actuelle, le code ne fonctionne qu'avec les maillages one.txt et two.txt
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
     // et un membre de droite nuls : la solution sera donc nulle :-)

            for (i = 0; i < theSpace->n; i++) {
                theSystem->A[map[i]][map[i]] = 1.0; }
            for (i = 0; i < theSpace->n; i++) {
                theSystem->B[map[i]] = femApproxStommel(theRule->xsi[iInteg], theRule->eta[iInteg]); }}}

    femFullSystemEliminate(theSystem);
}


# endif
