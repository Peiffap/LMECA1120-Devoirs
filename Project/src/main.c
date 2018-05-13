
/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *  Project for 17-18
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */


#include "glfem.h"

int main(void)
{


    int    n = 15;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.5;
    double radiusOut = 2.0;
    double dt      = 1e-1;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);

	femPoissonProblem* theProblem = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");


    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);

    femPoissonSolve(theProblem);

    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);

  //  A decommenter pour obtenir l'exemple de la seance d'exercice :-)
  //  femGrains* theGrains = femGrainsCreateTiny(radiusIn,radiusOut);;

    GLFWwindow* window = glfemInit("MECA1120 : Project 17-18");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(radiusOut,w,h);
        char theMessageTime[256];
        sprintf(theMessageTime,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessageTime);
		glfemPlotField(theProblem->mesh,theProblem->system->B);
		for (i=0 ;i < theGrains->n; i++) {
            glColor3f(1,1,1);
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);
  //
  // A decommenter pour pouvoir progresser pas par pas
  //          printf("press CR to compute the next time step >>");
  //          char c= getchar();
  //
            femGrainsUpdate(theGrains,dt,tol,iterMax);
            t += dt; }

        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS)
        	    theRunningMode = 1;
          if (glfwGetKey(window,'S') == GLFW_PRESS)
        	    theRunningMode = 0; }

    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	        (!glfwWindowShouldClose(window)));


    glfwTerminate();
    femGrainsFree(theGrains);
	femPoissonFree(theProblem);
    exit(EXIT_SUCCESS);
}
