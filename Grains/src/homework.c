#include"fem.h"


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
    int n = myGrains->n; 
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;
    
    double zeta = 0.0;
 
//
//  A FAIRE.... :-)    Difficile, difficile :-)
//

    return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
    int n = myGrains->n;
    int i,iter = 0;
    double zeta;
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

 
// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravité et de la trainee
//

//
//  A FAIRE.... :-)    La loi du grand Newton : so easy !
//

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       
    do {
        zeta = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((zeta > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,zeta);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//
    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}


#endif