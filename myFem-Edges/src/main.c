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
    femMesh  *theMesh  = femMeshRead("../data/conge.txt");    
    femEdges *theEdges = femEdgesCreate(theMesh);    
    
    edgesExpand(theEdges);               //   femEdgesPrint(theEdges);
    edgesSort(theEdges);                 //   femEdgesPrint(theEdges);
    edgesShrink(theEdges);               //   femEdgesPrint(theEdges);
    printf("Boundary edges  : %i \n", theEdges->nBoundary);
    printf("Boundary length : %14.7e \n", edgesBoundaryLength(theEdges));

    char theMessage[256];
    sprintf(theMessage, "Boundary edges : %i", theEdges->nBoundary);
      
//
//  On superpose le maillage (en bleu), 
//  tous les segments frontieres (en noir),
//  et la frontiere (en rouge)
//
//  Au depart de votre travail, vous devriez obtenir un maillage bleu....
//  et a la fin de l'exercice un maillage noir avec bord rouge :-)
//

    GLFWwindow* window = glfemInit("MECA1120 : homework 2 ");
    glfwMakeContextCurrent(window);
    do
    {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theMesh,w,h);
        glColor3f(0.4,0.4,1.0); glfemPlotMesh(theMesh);
        glColor3f(0.0,0.0,0.0); glfemPlotEdges(theEdges);  
        glColor3f(1.0,0.0,0.0); glfemPlotBnd(theEdges);          
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);              
        glfwSwapBuffers(window);
        glfwPollEvents();
    } 
   while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           (!glfwWindowShouldClose(window)) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femEdgesFree(theEdges);
    femMeshFree(theMesh);
    exit(EXIT_SUCCESS);
}



