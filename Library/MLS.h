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
#include <../Interpolations/WEIGHT.h>
// ======================================================================
// // polynomials
template <class T>
VectorXd PolyVec (T& x, T& xc, size_t d, size_t np, bool order2, double h)
{
	VectorXd p(np);
	T xr = (x-xc)/h;
	// T xr = x;

	p(0) = 1.;
	size_t count = 0;
	for (size_t i=0; i<d; ++i)
	{
		count++;
		p(count) = xr(i);
	}
	if (order2)
	{
		for (size_t i=0; i<d; ++i)
		for (size_t j=i; j<d; ++j)
		{
			count++;
			p(count) = xr(i)*xr(j);
		}
	}
	return p;
}

// // ======================================================================================
// // moving least square
template <class T>
VectorXd MLS (vector<T> l, T x, size_t d, size_t wtype, double h, double mu, bool order2)
{
	size_t np = 1+d;
	size_t nm = 0;
	if (order2)
	{
		switch(d)
		{
			case 1: 	nm = 1;
			break;
			case 2:		nm = 3;
			break;
			case 3:		nm = 6;
			break;
		}
		np += nm;
	}
	size_t n = l.size();		// number of data points 
	MatrixXd P(n,np);
	MatrixXd W(n,n);
	W.setZero();

	for (size_t i=0; i<n; ++i)
	{
		VectorXd pi = PolyVec(l[i],x,d,np,order2,h);
		P.row(i) = pi.transpose();
		double s = (l[i]-x).norm()/h;
		switch(wtype)
		{
			case 1:	W(i,i) = WeightQ(s);
			break;
			case 2: W(i,i) = WeightC(s);
			break;			
			case 3: W(i,i) = 1.;
			break;
		}
	}

	MatrixXd PW = P.transpose()*W;
	MatrixXd M = PW*P;	
	MatrixXd Mi = M.inverse();
	if (Mi.hasNaN())
	{
		for (size_t i=0; i<nm; ++i)
		{
			size_t j = np-i-1;
			M(j,j) += mu;
		}
		Mi = M.inverse();	
	}
	// VectorXd px = PolyVec(x,x,d,np,order2,h);
	// VectorXd phi = px.transpose()*Mi*PW;
	VectorXd phi = Mi.row(0)*PW;
	cout << "phi size: " << phi.size() << endl;
	return phi;
}

// // moving least square
template <class T>
MatrixXd LS (vector<T> l, size_t d, double mu, bool order2)
{
	size_t np = 1+d;
	size_t nm = 0;
	if (order2)
	{
		switch(d)
		{
			case 1: 	nm = 1;
			break;
			case 2:		nm = 3;
			break;
			case 3:		nm = 6;
			break;
		}
		np += nm;
	}
	size_t n = l.size();		// number of data points 
	MatrixXd P(n,np);
	MatrixXd W(n,n);
	W.setZero();

	T x;
	x.setZero();

	for (size_t i=0; i<n; ++i)
	{
		VectorXd pi = PolyVec(l[i],x,d,np,order2,1);
		P.row(i) = pi.transpose();
	}

	MatrixXd PW = P.transpose();
	MatrixXd M = PW*P;	
	MatrixXd Mi = M.inverse();
	if (Mi.hasNaN())
	{
		for (size_t i=0; i<nm; ++i)
		{
			size_t j = np-i-1;
			M(j,j) += mu;
		}
		Mi = M.inverse();	
	}
	// VectorXd px = PolyVec(x,x,d,np,order2,h);
	// VectorXd phi = px.transpose()*Mi*PW;
	MatrixXd A = Mi*PW;
	return A;
}