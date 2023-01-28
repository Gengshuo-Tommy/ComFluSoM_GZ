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

#ifndef WEIGHT_H
#define WEIGHT_H

// Weight functions
// quadratic spline weight function
// double WeightQ(double s)
// {
// 	double w = 0.;
// 	if (s<0.5)		w = 0.75-s*s;
// 	else if (s<1.5)	w = 0.5*pow(1.5-s,2);
// 	return w;
// }
double WeightQ(double s)
{
	double w = 0.;
	if (s<1.)
	{
		double ss = s*s;
		w = 1.-6.*ss + 8.*ss*s - 3.*ss*ss;
	}
	return w;
}
// // cubic spline weight function
double WeightC(double s)
{
	double w = 0.;
	if (s<0.5)		w = 0.66666666666+4.*s*s*(s-1.);
	else if (s<1.)	w = 1.33333333333+4.*s*(-1.+s-s*s/3.);
	return w;
}

double Delta3P(double s)
{
	double w = 0.;
	if (s<0.5)		w = 0.33333333333*(1.+sqrt(-3.*s*s+1.));
	else if (s<1.5)	w = 0.16666666666*(5.-3.*s-sqrt(-3.*(1.-s)*(1.-s)+1.));
	return w;
}

double Delta4P(double s)
{
	double w = 0.;
	if (s<1.)		w = 0.125*(3.-2.*s+sqrt(1.+4.*s-4.*s*s));
	else if (s<2.)	w = 0.125*(5.-2.*s-sqrt(-7.+12.*s-4.*s*s));
	return w;
}

#endif