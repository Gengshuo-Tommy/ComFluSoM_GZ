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

#ifndef EFFECTIVE_VISCOSITY_MODEL_H
#define EFFECTIVE_VISCOSITY_MODEL_H

// A coupled lattice Boltzmann method and discrete element method for discrete particle simulations of particulate flows
// miu0: viscosity of fluid, phi: solid fraction, phiC: critical solid fraction, iMiu: intrinsic viscosity, 
double CalEffectiveViscosity(double phiC, double iMiu, double phi)
{
    double a = 1.+0.5*iMiu*phi*phiC/(1.- phi);
    double miu = a*a;
    return miu;
}

#endif