/*
 *  glfem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2015 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1.2)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fem.h"



#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>



void        glfemDrawElement(double *x, double *y, int n);
void        glfemDrawNodes(double *x, double *y, int n);

void        glfemMessage(char *message);

void        glfemMakeRasterFont(void);
void        glfemDrawMessage(int h, int v, char *message);
void        glfemSetRasterSize(int width, int height);

void        glfemReshapeWindows(double *x, double *y, int n, int w, int h);
void        glfemDraw(void (*glFunction)(int width, int height));
GLFWwindow* glfemInit(char *windowName);

#endif


