/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2021 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License fo`r more details.                        *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef DRAG_FORCE_MODEL_H
#define DRAG_FORCE_MODEL_H

// Drag force models

void DragForceDarcy(double rho, double nu, double phi, double dp, Vector3d vf, Vector3d vs, Vector3d& drag, Matrix3d& dfdu)
{
	double n = 1.-phi;								// porosity
	Vector3d vsf = vf-vs;							// Relative velocity

	// vsf *= n;

	double fe = 1.75/sqrt(150.*n*n*n);
	double k = n*n*n*dp*dp/(150.*phi*phi);

	double a = -n*nu/k;
	double b = -n*fe/sqrt(k);
	double vn = vsf.norm();

	drag = rho*(a + b*vsf.norm())*vsf;

	if (vn>0.)
	{
		dfdu(0,0) = a+b*vn+b*vsf(0)*vsf(0)/vn;
		dfdu(1,1) = a+b*vn+b*vsf(1)*vsf(1)/vn;
		dfdu(2,2) = a+b*vn+b*vsf(2)*vsf(2)/vn;

		dfdu(0,1) = b*vsf(0)*vsf(1)/vn;
		dfdu(0,2) = b*vsf(0)*vsf(2)/vn;
		dfdu(1,0) = b*vsf(1)*vsf(0)/vn;
		dfdu(1,2) = b*vsf(1)*vsf(2)/vn;
		dfdu(2,0) = b*vsf(2)*vsf(0)/vn;
		dfdu(2,1) = b*vsf(2)*vsf(1)/vn;
	}


	// cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	// cout << "rho: " << rho << endl;
	// cout << "nu: " << nu << endl;
	// cout << "phi: " << phi << endl;
	// cout << "dp: " << dp << endl;
	// cout << "vf: " << vf.transpose() << endl;
	// cout << "vs: " << vs.transpose() << endl;

	// cout << "in fe: " << fe << endl;
	// cout << "in k: " << k << endl;
	// cout << "in a: " << a << endl;
	// cout << "in b: " << b << endl;
	// cout << "rho*(a + b*vsf.norm()): " << rho*(a + b*vsf.norm()) << endl;
	// cout << "in drag: " << drag.transpose() << endl;
	// // abort();
	// if (phi==0.51)	abort();
	// cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

	// return drag;
}

// Drag Force of Intermediate Reynolds Number Flow Past Mono- and Bidisperse Arrays of Spheres
void DragForceHoef(double rho, double nu, double phi, double dp, Vector3d vf, Vector3d vs, Vector3d& drag, Matrix3d& dfdu)
{
	double n = 1.-phi;								// porosity
	Vector3d vsf = vf-vs;							// Relative velocity

	double n2 = n*n;
	double re = n*dp*vsf.norm()/nu;					// Reynolds number

	double fHat0 = 10.*phi/n2 + n2*(1.+1.5*sqrt(phi));
	double fHat = fHat0;

	double aa = -18.*phi*n*nu*rho/(dp*dp);
	double bb = 0.01720833/n2;
	double pp = 1./n+3.*n*phi+8.4*pow(re, -0.343);
	double qq = 1.+pow(10.,3.*phi)*pow(re,-0.5*(1.+4.*phi));
	double dpp = -2.8812*pow(re, -1.343);
	double dqq = -pow(10.,3.*phi)*0.5*(1.+4.*phi)*pow(re,-0.5*(3.+4.*phi));

	if (re>1.e-12)	fHat += bb*re*pp/qq;
	drag = aa*fHat*vsf;	// volume averaged drag force (N/m^3)

	double dfHdre = bb*pp/qq + bb*re*(dpp*qq - pp*dqq)/(qq*qq);
	
	double vn = vsf.norm();

	double dreC = n*dp/nu/vn;

	if (vn>0.)
	{
		double coef = aa*dfHdre*dreC;
		// coef = 0.;
		dfdu(0,0) = aa*fHat + coef*vsf(0)*vsf(0);
		dfdu(1,1) = aa*fHat + coef*vsf(1)*vsf(1);
		dfdu(2,2) = aa*fHat + coef*vsf(2)*vsf(2);

		dfdu(0,1) = coef*vsf(0)*vsf(1);
		dfdu(0,2) = coef*vsf(0)*vsf(2);
		dfdu(1,0) = coef*vsf(1)*vsf(0);
		dfdu(1,2) = coef*vsf(1)*vsf(2);
		dfdu(2,0) = coef*vsf(2)*vsf(0);
		dfdu(2,1) = coef*vsf(2)*vsf(1);
	}

	// cout << "re: " << re << endl;
	// cout << "drag: " << drag.transpose() << endl;
	// abort();
}

// vfd: darcy velocity of fluid, vs: solid velocity, acc: acceleration, ve: the velocity from sum fi*ei, vc: corrected velocity, fe: external force
void VelDragForceHoef(double rho, double nu, double phi, double dp, Vector3d vfd, Vector3d vs, Vector3d acc, Vector3d ve, Vector3d& vc, Vector3d& fe)
{
	double n = 1.-phi;								// porosity
	Vector3d vf = vfd/n;
	Vector3d vsf = 1.0e-5*(vf-vs);				    // Relative velocity

	double n2 = n*n;
	double re = n*dp*vsf.norm()/nu;					// Reynolds number

	double fHat0 = 10.*phi/n2 + n2*(1.+1.5*sqrt(phi));
	double fHat = fHat0;

	double aa = -18.*phi*n*nu*rho/(dp*dp);
	double bb = 0.01720833/n2;
	double pp = 1./n+3.*n*phi+8.4*pow(re, -0.343);
	double qq = 1.+pow(10.,3.*phi)*pow(re,-0.5*(1.+4.*phi));
	double dpp = -2.8812*pow(re, -1.343);
	double dqq = -pow(10.,3.*phi)*0.5*(1.+4.*phi)*pow(re,-0.5*(3.+4.*phi));

	if (re>1.e-12)	fHat += bb*re*pp/qq;
	Vector3d drag = aa*fHat*vsf;	// volume averaged drag force (N/m^3)

	double dfHdre = bb*pp/qq + bb*re*(dpp*qq - pp*dqq)/(qq*qq);
	
	double vn = vsf.norm();

	double dreC = n*dp/nu/vn;

	Matrix3d dfdu = Matrix3d::Zero();

	if (vn>0.)
	{
		double coef = aa*dfHdre*dreC;
		// coef = 0.;
		dfdu(0,0) = aa*fHat + coef*vsf(0)*vsf(0);
		dfdu(1,1) = aa*fHat + coef*vsf(1)*vsf(1);
		dfdu(2,2) = aa*fHat + coef*vsf(2)*vsf(2);

		dfdu(0,1) = coef*vsf(0)*vsf(1);
		dfdu(0,2) = coef*vsf(0)*vsf(2);
		dfdu(1,0) = coef*vsf(1)*vsf(0);
		dfdu(1,2) = coef*vsf(1)*vsf(2);
		dfdu(2,0) = coef*vsf(2)*vsf(0);
		dfdu(2,1) = coef*vsf(2)*vsf(1);
	}

	Vector3d v = ve + 0.5*n*acc;
	Vector3d rh = rho*(v-vfd) + 0.5*drag;

	Vector3d dudt = (n*rho*Matrix3d::Identity() - 0.5*dfdu).inverse()*rh;

	vc = (vf + dudt)*n;
	fe = 2.*rho*(vc-v);
}

// Vector3d DragForceTenneti(double rho, double nu, double phi, double dp, Vector3d vf, Vector3d vs)
// {
// 	double n = 1.-phi;								// porosity
// 	Vector3d vsf = vs-vf;							// Relative velocity	
// 	double re = n*dp*vsf.norm()/nu;					// Reynolds number

// 	double cd0 = 1.+0.15*pow(re, 0.687);
// 	double a = 5.81*phi/pow(n, 3) + 0.48*pow(phi, 1./3.)/pow(n, 4);
// 	double b = phi*re*(0.95 + 0.61*pow(phi, 3)/(n*n));

// 	double cd = n*(cd0/pow(n,3) + a + b);

// 	Vector3d drag = 3.*M_PI*dp*n*cd*vsf;
// 	return drag;
// }

#endif