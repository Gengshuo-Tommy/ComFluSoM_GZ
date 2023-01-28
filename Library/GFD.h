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

#include "../HEADER.h"
#include <WEIGHT.h>

// ======================================================================================
// only for 2d now
VectorXd GFD (vector<Vector3d> l, vector<Vector3d> lu, Vector3d x, Vector3d u, size_t wtype, double h, bool order2)
{
	MatrixXd A(5,5);
	// A.setZero();

	VectorXd b0(5);
	VectorXd b1(5);

	size_t n = l.size();

	for (size_t i=n; ++i)
	{
		double h = l[i](0)-x(0);
		double k = l[i](1)-x(1);
		double hh = h*h;
		double hk = h*k;
		double kk = k*k;
		double hhh = hh*h;
		double hhk = hh*k;
		double hkk = hk*k;
		double kkk = kk*k;
		double hhhh= hh*hh;
		double hhhk= hh*hk;
		double hhkk= hh*kk;
		double hkkk= hk*kk;
		double kkkk= kk*kk;

		double s = sqrt(hh+kk)/h;
		double w = WeightQ(s);
		double w2 = w*w;

		A(0,0) += hh*w2;
		A(0,1) += hk*w2;
		A(0,2) += hhh*w2;
		A(0,3) += hkk*w2;
		A(0,4) += hhk*w2;

		A(1,1) += kk*w2;
		A(1,3) += kkk*w2;

		A(2,2) += hhhh*w2;
		A(2,3) += hhkk*w2;
		A(2,4) += hhhk*w2;

		A(3,3) += kkkk*w2;
		A(3,4) += hkkk*w2;

		double b0i = lu(i)(0)-u(0);
		b0(0) += b0i*h*w2;
		b0(1) += b0i*k*w2;
		b0(2) += b0i*hh*w2;
		b0(3) += b0i*kk*w2;
		b0(4) += b0i*hk*w2;

		double b1i = lu(i)(1)-u(1);
		b1(0) += b0i*h*w2;
		b1(1) += b0i*k*w2;
		b1(2) += b0i*hh*w2;
		b1(3) += b0i*kk*w2;
		b1(4) += b0i*hk*w2;
	}

	A(0,2) *= 0.5;
	A(1,3) *= 0.5;
	A(1,2) = 0.5*A(0,4);
	A(1,4) = A(0,3);
	A(0,3) *= 0.5;
	A(2,2) *= 0.25;
	A(4,4) = A(2,3);
	A(2,3) *= 0.25;
	A(2,4) *= 0.5;
	A(3,3) *= 0.25;
	A(3,4) *= 0.5;

	b0(2) *= 0.5;
	b0(3) *= 0.5;

	b1(2) *= 0.5;
	b1(3) *= 0.5;

	for (size_t m=0; m<5; ++m)
	for (size_t n=0; n<m; ++n)
	{
		A(m,n) = A(n,m);
	}
}