#include"fem.h"


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
	int n = myGrains->n;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *r = myGrains->r;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double *dvBoundary = myGrains->dvBoundary;
	double *dvContacts = myGrains->dvContacts;
	double rIn = myGrains->radiusIn;
	double rOut = myGrains->radiusOut;
	double tol = 1e-6;

	double zeta = 0.0;

	int i;
	int j;
	int coll = 0;
	double distG = 0.0;
	double gamma = 0.0;
	double vn = 0.0;
	double nx = 0.0;
	double ny = 0.0;
	double rx = 0.0;
	double ry = 0.0;
	double dv = 0.0;
	double dvx = 0.0;
	double dvy = 0.0;
	double delta_dvContacts = 0.0;

	if (iter == 0) {
		for (i = 0; i < n*(n - 1) / 2; i++) {
			dvContacts[i] = 0;
		}
		for (i = 0; i < n; i++) {
			dvBoundary[i] = 0;
		}

		return zeta;
	}

	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			rx = x[j] - x[i];
			ry = y[j] - y[i];
			distG = sqrt(pow(rx, 2) + pow(ry, 2));
			if (distG < (rOut-rIn) /4) {
				gamma = distG - r[i] - r[j];
				nx = rx / distG;
				ny = ry / distG;
				vn = (vx[i] * nx + vy[i] * ny) - (vx[j] * nx + vy[j] * ny);

				if (vn + dvContacts[coll] - (gamma / dt) >= 0) {
					dv = fmax(0.0, vn + dvContacts[coll] - (gamma / dt)) - dvContacts[coll];

					dvx = -dv * nx*m[j] / (m[j] + m[i]);
					dvy = -dv * ny*m[j] / (m[j] + m[i]);
					vx[i] += dvx;
					vy[i] += dvy;
					dvx = dv * nx*m[i] / (m[j] + m[i]);
					dvy = dv * ny*m[i] / (m[j] + m[i]);
					vx[j] += dvx;
					vy[j] += dvy;

					delta_dvContacts = dv;
					dvContacts[coll] += delta_dvContacts;

					zeta = fmax(zeta, fabs(dv));
				}
			}
			

			coll++;
		}
	}

	double distC = 0.0;
	double gammaIn = 0.0;
	double gammaOut = 0.0;
	double delta_dvBoundary = 0.0;
	double dvo = 0.0;
	double dvi = 0.0;


	for (i = 0; i < n; i++) {
		distC = sqrt(pow(x[i], 2) + pow(y[i], 2));
		nx = x[i] / distC;
		ny = y[i] / distC;
		gammaOut = rOut - distC - r[i];
		gammaIn = distC - rIn - r[i];
		vn = vx[i] * nx + vy[i] * ny;

		dvo = fmax(0, vn + dvBoundary[i] - (gammaOut / dt));
		dvi = fmax(0, -vn - dvBoundary[i] - (gammaIn / dt));
		dv = dvo - dvi - dvBoundary[i];

		dvx = -dv * nx;
		vx[i] += dvx;
		dvy = -dv * ny;
		vy[i] += dvy;

		delta_dvBoundary = dv;
		dvBoundary[i] += delta_dvBoundary;

		zeta = fmax(zeta, fabs(dv));
	}

	return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femDiffusionProblem *theProblemX, femDiffusionProblem *theProblemY, double t)
//void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femDiffusionProblem *theProblemX, double t)
{
	int n = myGrains->n;
	int i, j, k, iter = 0;
	double zeta;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double gamma = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];
	double vxfluide = 1.0;
	double vyfluide = 1.0;
	double X1, X2, X3, Y1, Y2, Y3, Xb, Yb, J1, J2, J3;
	double phi[4];
	int currElem = 0;


	// 
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
	//
	double dvx = 0;
	double dvy = 0;
	//printf("gy = %f", gy);
	for (i = 0; i < n; i++) {
		
		vxfluide = 0;
		vyfluide = 0;
		
		if (t != 0) {
			femMesh *theMesh = theProblemX->mesh;
			for (j = 0; j < theMesh->nElem; j++) {

				//printf("caribou");
				X1 = theMesh->X[theMesh->elem[j * 3 + 0]];
				X2 = theMesh->X[theMesh->elem[j * 3 + 1]];
				X3 = theMesh->X[theMesh->elem[j * 3 + 2]];
				Xb = myGrains->x[i];
				Y1 = theMesh->Y[theMesh->elem[j * 3 + 0]];
				Y2 = theMesh->Y[theMesh->elem[j * 3 + 1]];
				Y3 = theMesh->Y[theMesh->elem[j * 3 + 2]];
				Yb = myGrains->y[i];
				J1 = (X2 - X1)*(Yb - Y1) - (Xb - X1)*(Y2 - Y1);
				J2 = (X2 - Xb)*(Y3 - Yb) - (X3 - Xb)*(Y2 - Yb);
				J3 = (Xb - X1)*(Y3 - Y1) - (X3 - X1)*(Yb - Y1);
				if (J1 > 0 && J2 > 0 && J3 > 0) {
					currElem = j;
					j = theMesh->nElem;
				}
				//printf("carotte");
			}


			//printf("gy = %f",gy);

			double x = Xb;
			double y = Yb;

			double xsi = (x*(Y3 - Y1) - (X3 - X1)*y + (X3 - X1)*Y1 - X1 * (Y3 - Y1)) / ((Y3 - Y1)*(X2 - X1) - (X3 - X1)*(Y2 - Y1));
			double eta = (y - (Y2 - Y1)*xsi - Y1) / (Y3 - Y1);
			femDiscretePhi2(theProblemX->space, xsi, eta, phi);


			for (k = 0; k < 3; k++) {
				vxfluide += theProblemX->soluce[theMesh->elem[currElem * 3 + k]] * phi[k];
				vyfluide += theProblemY->soluce[theMesh->elem[currElem * 3 + k]] * phi[k];
			}
		}
		
		
		//vxfluide = vyfluide = 0;
		
		dvx = gx * dt - (gamma * dt / m[i])* (vx[i] - vxfluide);
		vx[i] += dvx;
		dvy = gy * dt - (gamma * dt / m[i])* (vy[i] - vyfluide);
		vy[i] += dvy;
		//printf("vy = %f\n", vy[0]);

		//printf("dindon");
	}
	//printf("goéland");
	//
	// -2- Correction des vitesses pour tenir compte des contacts        
	//       
	
	do {
		zeta = femGrainsContactIterate(myGrains, dt, iter); //zeta = delta_ v
		iter++;
	} while ((zeta > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, zeta);
	
	//printf("vy = %f\n", vy[0]);

	//  
	// -3- Calcul des nouvelles positions sans penetrations de points entre eux
	//
	for (i = 0; i < n; ++i) {
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
		//printf("vy = %f\n", vy[i]);
		//printf("y = %f",y[i]);
	}
}


#endif
