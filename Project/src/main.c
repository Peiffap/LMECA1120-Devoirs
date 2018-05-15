
/* C EST NOTRE FONCTIOOOOOOOOOOOOOOOOOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *  Homework 4 for 17-18 : Discrete Grains
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
 
#include "glfem.h"
#include <time.h>


int main(void)
{  
    
    
    int    n = 15;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt      = 1e-1;
    double tEnd    = 20.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;
    //femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);
	femGrains* theGrains = femGrainsCreateTiny(radiusIn, radiusOut);
   
  //  A decommenter pour obtenir l'exemple de la seance d'exercice :-)
	//femGrains* theGrains = femGrainsCreateTiny(radiusIn, radiusOut);

	femSolverType solverType = FEM_BAND;
	femRenumType  renumType = FEM_XNUM;
	char meshFileName[] = "..\\data\\tiny.txt";

	femDiffusionProblem* theProblemX = NULL;
	femDiffusionProblem* theProblemY = NULL;
	femDiffusionProblem* theProblem = NULL;

	femGrainsUpdate(theGrains, dt, tol, iterMax, theProblemX, theProblemY, t);

	theProblemX = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
	theProblemY = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
	theProblem = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
	clock_t tic = clock();
	femDiffusionCompute(theProblem, theGrains);
	femDiffusionCompute(theProblemX, theGrains);
	femDiffusionComputeY(theProblemY, theGrains);
	femDiffusionCombo(theProblem, theProblemX, theProblemY);
	t += dt;
	
	femSolverPrintInfos(theProblem->solver);
	printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 / CLOCKS_PER_SEC);
	printf("    Maximum value : %.4f\n", femMax(theProblem->soluce, theProblem->size));
	fflush(stdout);


	int option = 1;
	femSolverType newSolverType = solverType;
	femRenumType  newRenumType = renumType;


    GLFWwindow* window = glfemInit("MECA1120 : Projet");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1;
    float theVelocityFactor = 0.25;
     

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

		glfwGetFramebufferSize(window, &w, &h);
		glfemReshapeWindowsMesh(theProblem->mesh, w, h);
		glfemPlotField(theProblem->mesh, theProblem->soluce);
		glColor3f(1.0, 0.0, 0.0); glfemPlotBndIn(theProblem->edges);
		glColor3f(0.0, 0.0, 1.0); glfemPlotBndOut(theProblem->edges);
        for (i=0 ;i < theGrains->n; i++) {
            glColor3f(0.184313,0.5,0.309804); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); }       
        //glColor3f(1,0,1); glfemDrawCircle(0,0,radiusOut);
        //glColor3f(1,0,1); glfemDrawCircle(0,0,radiusIn); 
        char theMessage[256];
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);    
        glfwSwapBuffers(window);
        glfwPollEvents();
 
        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);  
  //
  // A decommenter pour pouvoir progresser pas par pas
          //  printf("press CR to compute the next time step >>");
           // char c= getchar();
  //
            femGrainsUpdate(theGrains,dt,tol,iterMax, theProblemX, theProblemY, t);

			theProblemX = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
			theProblemY = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
			theProblem = femDiffusionCreate(meshFileName, solverType, renumType, theGrains);
			femDiffusionCompute(theProblem, theGrains);
			femDiffusionCompute(theProblemX, theGrains);
			femDiffusionComputeY(theProblemY, theGrains);
			femDiffusionCombo(theProblem, theProblemX, theProblemY);
            t += dt; 
		}
         
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
    exit(EXIT_SUCCESS);
}



