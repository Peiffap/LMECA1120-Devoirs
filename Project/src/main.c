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
#include <time.h>
#include <math.h>

int main(void)
{


    int    n = 42;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;
    double dt      = 5e-2;
    double tEnd    = 20.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;

    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);

	// femPoissonProblem* theProblem = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");

	femDiffusionProblem* theProblemX = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);
	femDiffusionProblem* theProblemY = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);
	femDiffusionProblem* theProblem = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);

	femGrainsUpdate(theGrains, theProblemX, theProblemY, dt, tol, iterMax);

    // femPoissonSolve(theProblem);
	femDiffusionCompute(theProblem, theGrains, 0);
	femDiffusionCompute(theProblemX, theGrains, 0);
	femDiffusionCompute(theProblemY, theGrains, 1);
	femDiffusionCombine(theProblem, theProblemX, theProblemY);
	femSolverPrintInfos(theProblem->solver);

    // printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->size));
    // fflush(stdout);

  //  A decommenter pour obtenir l'exemple de la seance d'exercice :-)
  //  femGrains* theGrains = femGrainsCreateTiny(radiusIn,radiusOut);;

    GLFWwindow* window = glfemInit("MECA1120 : Project 17-18");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.05;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(radiusOut,w,h);
        char theMessageTime[256];
        sprintf(theMessageTime,"Time = %g sec",t);
		/*
		clock_t tic = clock();
		int testConvergence;
		do {
			femDiffusionCompute(theProblem);
			femSolverPrintInfos(theProblem->solver);
			testConvergence = femSolverConverged(theProblem->solver);
			printf("Thomas. \n");
		}
		while (testConvergence == 0);
		if (testConvergence == -1)  printf("    Iterative solver stopped after a maximum number of iterations\n");
		printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
		*/
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessageTime);
		glfemPlotField(theProblem->mesh, theProblem->soluce);
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
            femGrainsUpdate(theGrains, theProblemX, theProblemY, dt, tol, iterMax);
            t += dt; }

			theProblemX = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);
			theProblemY = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);
			theProblem = femDiffusionCreate("../data/meca1120-projet-meshMedium.txt", FEM_BAND, FEM_XNUM, theGrains);

			femDiffusionCompute(theProblem, theGrains, 0);
			femDiffusionCompute(theProblemX, theGrains, 0);
			femDiffusionCompute(theProblemY, theGrains, 1);
			femDiffusionCombine(theProblem, theProblemX, theProblemY);

        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS)
        	    theRunningMode = 1;
          if (glfwGetKey(window,'S') == GLFW_PRESS)
        	    theRunningMode = 0; }



    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	        (!glfwWindowShouldClose(window)));


    glfwTerminate();
	femDiffusionFree(theProblem);
	femDiffusionFree(theProblemX);
	femDiffusionFree(theProblemY);
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
}
