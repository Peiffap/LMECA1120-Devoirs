
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1


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


femGrains  *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut);
femGrains  *femGrainsCreateTiny(double radiusIn, double radiusOut);
void        femGrainsFree(femGrains *myGrains);
void        femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax);
double      femGrainsContactIterate(femGrains *myGrains, double dt, int iter); 

double      femMin(double *x, int n);
double      femMax(double *x, int n);
void        femError(char *text, int line, char *file);
void        femErrorScan(int test, int line, char *file);
void        femWarning(char *text, int line, char *file);


#endif
