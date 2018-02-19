#include <stdio.h>
#include <math.h>


double integrate(double x[3], double y[3], double(*f)(double,double));
double integrateRecursive(double x[3], double y[3], double(*f)(double,double), int n);

double fun(double x, double y)  { return cos(x) + y * y; }
double stupid(double x, double y)  { return 1.0; }


// -2- Version non-graphique du code
//     Mettre le flag graphic dans l'instruction de compilation



#ifndef graphic

int main()
{
    double x[3] = {-1, 1, 2};
    double y[3] = {0, 5, 1};
    int n;
    
    printf("Surface integration    : %14.7e \n", integrate(x,y,stupid));
    printf("More funny integration : %14.7e \n", integrate(x,y,fun));
    for (n=0;  n <= 4; n++) 
        printf("Recursive integration (n = %2d) : %14.7e \n", n, integrateRecursive(x,y,fun,n));      
    return 0;
}

#endif

// -2- Version graphique du code
//     Mettre le flag graphic dans l'instruction de compilation

#ifdef graphic


#include "glfem.h"

int main()
{
    double x[3] = {0, 1, 0};
    double y[3] = {0, 0, 1};
    char theMessage[256];

  //  double x[3] = {-1, 1, 2};
  //  double y[3] = {0, 5, 1};
    
    GLFWwindow* window = glfemInit("MECA1120 : homework 1 ");
    glfwMakeContextCurrent(window);
    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(x,y,3,w,h);

       
        
        double I = integrateRecursive(x,y,fun,2);
        sprintf(theMessage,"Integral = %14.7e",I);
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);    

          
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );

    glfwTerminate();    
    exit(EXIT_SUCCESS);
    
    return 0;
}

#endif