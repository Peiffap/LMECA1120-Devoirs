/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{

    femPoissonProblem* theProblem = femPoissonCreate("../data/tiny.txt");

    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt")
    // par :
    // ("..\\data\\triangles_166.txt")
    //
    // Sorry for the inconvenience :-)
    // On r�fl�chit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !


    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);

    femPoissonSolve(theProblem);

    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);

    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));

    GLFWwindow* window = glfemInit("MECA1120 : homework 3 ");
    glfwMakeContextCurrent(window);
    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotField(theProblem->mesh,theProblem->system->B);
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );

    // Check if the ESC key was pressed or the window was closed

    glfwTerminate();
    femPoissonFree(theProblem);
    exit(EXIT_SUCCESS);



}
