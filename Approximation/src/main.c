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
 
    femApproxProblem* theProblem = femApproxCreate("../data/two.txt");
  //  femMeshWrite(theProblem->mesh, "../data/triangles_22_bis.txt");
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
     
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);
    femApproxSolve(theProblem);   
 
    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);
    
    char theMessage[256];
    theProblem->maxValue = femMax(theProblem->system->B,theProblem->system->size);
    sprintf(theMessage, "Max : %.4f", theProblem->maxValue );
    
    GLFWwindow* window = glfemInit("MECA1120 : homework 6 ");
    glfwMakeContextCurrent(window);
    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotFieldRecursiveTriangle(theProblem,5);
        glColor3f(1.0,1.0,0.0); glfemDrawIntegrationNodes(theProblem);
        glColor3f(1.0,0.0,0.0); glfemDrawValueNodes(theProblem);
        glColor3f(0.0,0.0,0.0); glfemDrawNodes(theProblem->mesh->X,theProblem->mesh->Y,theProblem->mesh->nNode);
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femApproxFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

