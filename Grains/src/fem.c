
/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGrains *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut)
{
    int i,nContact = n*(n-1)/2;
    
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = n;
    theGrains->radiusIn = radiusIn;
    theGrains->radiusOut = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -9.81;
    theGrains->gamma = 0.5;
    
       
    theGrains->x  = malloc(n*sizeof(double));
    theGrains->y  = malloc(n*sizeof(double));
    theGrains->vx = malloc(n*sizeof(double));
    theGrains->vy = malloc(n*sizeof(double));
    theGrains->r  = malloc(n*sizeof(double));
    theGrains->m  = malloc(n*sizeof(double));       
    theGrains->dvBoundary = malloc(n * sizeof(double));
    theGrains->dvContacts = malloc(nContact * sizeof(double));
   
    for(i = 0; i < n; i++) {
        theGrains->r[i] = r;
        theGrains->m[i] = m;
        theGrains->x[i] = (i%5) * r * 2.5 - 5 * r + 1e-8; 
        theGrains->y[i] = (i/5) * r * 2.5 + 2 * r + radiusIn;
        theGrains->vx[i] = 0.0;
        theGrains->vy[i] = 0.0; 
        theGrains->dvBoundary[i] = 0.0; }
 
    for(i = 0; i < nContact; i++)  
        theGrains->dvContacts[i] = 0.0;

  
    return theGrains;
}

femGrains *femGrainsCreateTiny(double radiusIn, double radiusOut)
{
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = 2;
    theGrains->radiusIn   = radiusIn;
    theGrains->radiusOut  = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -10;
    theGrains->gamma = 0.0;
           
    theGrains->x  = malloc(2 * sizeof(double));
    theGrains->y  = malloc(2 * sizeof(double));
    theGrains->vx = malloc(2 * sizeof(double));
    theGrains->vy = malloc(2 * sizeof(double));
    theGrains->r  = malloc(2 * sizeof(double));
    theGrains->m  = malloc(2 * sizeof(double));       
    theGrains->dvBoundary = malloc(2 * sizeof(double));
    theGrains->dvContacts = malloc(sizeof(double));
   
    theGrains->r[0] = 0.1;
    theGrains->r[1] = 0.1;  
    theGrains->m[0] = 1.0;
    theGrains->m[1] = 1.0; 
    theGrains->x[0] = 0.0; 
    theGrains->x[1] = 0.0;
    theGrains->y[0] = -radiusOut + 0.3; 
    theGrains->y[1] = -radiusOut + 0.1; 
    theGrains->vx[0] = 0.0; 
    theGrains->vx[1] = 0.0;            
    theGrains->vy[0] = 0.0; 
    theGrains->vy[1] = 0.0; 
    theGrains->dvBoundary[0] = 0.0;
    theGrains->dvBoundary[1] = 0.0;

    theGrains->dvContacts[0] = 0.0;
     
    return theGrains;
}

void femGrainsFree(femGrains *theGrains)
{
    free(theGrains->x);
    free(theGrains->y);
    free(theGrains->vx);
    free(theGrains->vy);
    free(theGrains->r);
    free(theGrains->m);
    free(theGrains->dvBoundary);
    free(theGrains->dvContacts);
    free(theGrains);
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
