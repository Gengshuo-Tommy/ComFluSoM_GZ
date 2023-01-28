/***********************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
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
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MPM_NODE_H
#define MPM_NODE_H

class MPM_NODE
{
public:
	MPM_NODE();
	MPM_NODE(size_t level, const Vector3d& x);
	void Reset();
	void ResetwithDEM(size_t nproc);
	void CombineDPs();
	void NonSlippingBC();
	void SlippingBC(Vector3d& norm);
	void FrictionBC(double dt, Vector3d& n);

	bool 							Actived;					// Flag of actived node
	bool  							IsBoundary;					// Flag of boundary ndoe

	size_t 							Type;						// 0 for nodes, 1 for gauss quadrature points
	size_t 							ID; 				    	// Index of gird in the list 
	size_t 							Level; 				    	// Level of the node, 0 for basic node, 1 for level 1 refined node

	double							M;							// Mass
	double                          Mp;                         // Porous mass
	double  						Mbar;						// Smoothed mass for surface tension
	double							Vol;						// Volume
	double                          Poro;                       // Porosity
	double 							Mu;							// Friction coefficient

	Vector3d						X;							// Position
	Vector3d						V;							// Velocity
	Vector3d						Mv;							// Momentum
	Vector3d						F;							// Total force
	Vector3d						Fi;							// Internal force
	Vector3d						Fe;							// External force
	Vector3d  						SurfNorm;					// Normal direction for surface tension
	Vector3d  						GradM;						// Mbar gradient for surface tension

	Matrix3d						Stress;						// Stress

	Matrix3d						VGrad;						// Velocity gradient

	vector<size_t>					GPs;						// Gauss quadrature points
	vector<double> 					Ng;							// Shape function
	vector<Vector3d>				GNg;						// Gradient of shape function
	vector<size_t>					MPs;						// Material points
	vector<size_t>					DPs;						// DEM particles
	vector<vector<size_t>>			DPs_proc;					// DEM particles for openMP to avoid race condition
	vector<size_t> 					BCTypes;					// Boundary condition type
	vector<Vector3d> 				Norms;						// Normal direction for boundaries

	// unordered_map<size_t, bool> 	PMap;						// Map for all particles which have influence on this node
};

inline MPM_NODE::MPM_NODE()
{
	Actived = false;
	IsBoundary = false;
	Type 	= 0;
	ID 		= 0;
	Level 	= 0;
	M 		= 0.;
	Mp      = 0.;
	Poro    = 0.;
	Mbar 	= 0.;
	Vol 	= 1.;
	X 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	SurfNorm = Vector3d::Zero();
	GradM 	= Vector3d::Zero();
	Stress.setZero();
	VGrad.setZero();
	GPs.resize(0);
	MPs.resize(0);
	DPs.resize(0);
	DPs_proc.resize(0);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline MPM_NODE::MPM_NODE(size_t level, const Vector3d& x)
{
	Actived = false;
	IsBoundary = false;
	Type 	= 0;
	ID 		= 0;
	Level 	= level;
	M 		= 0.;
	Mp      = 0.;
	Poro    = 0.;
	Mbar 	= 0.;
	Vol 	= 1.;
	X 		= x;
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	Stress.setZero();
	GPs.resize(0);
	MPs.resize(0);
	DPs.resize(0);
	DPs_proc.resize(0);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline void MPM_NODE::Reset()
{
	Actived = false;
	IsBoundary = false;
	M = 0.;
	Mbar = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
	Fi.setZero();
	Fe.setZero();
	SurfNorm.setZero();
	GradM.setZero();
	Stress.setZero();
	VGrad.setZero();
	// MPs.resize(0);
	MPs.clear();
	DPs.clear();
	// DPs_proc.resize(nproc);
	for (size_t i=0; i<DPs_proc.size(); ++i)
	{
		DPs_proc[i].clear();
	}
}

inline void MPM_NODE::ResetwithDEM(size_t nproc)
{
	Actived = false;
	IsBoundary = false;
	M = 0.;
	Mbar = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
	Fi.setZero();
	Fe.setZero();
	SurfNorm.setZero();
	GradM.setZero();
	Stress.setZero();
	MPs.clear();
	DPs.clear();
	// DPs_proc.resize(nproc);
	for (size_t i=0; i<DPs_proc.size(); ++i)
	{
		DPs_proc[i].clear();
	}
}

inline void MPM_NODE::CombineDPs()
{
	// no need to sort the list since DEM is parallized base on particles
	for (size_t i=0; i<DPs_proc.size(); ++i)
	{
		DPs.insert( DPs.end(), DPs_proc[i].begin(), DPs_proc[i].end() );
	}
	// sort( DPs.begin(), DPs.end() );
	// DPs.erase(unique(DPs.begin(), DPs.end()), DPs.end());
}

inline void MPM_NODE::NonSlippingBC()
{
	F.setZero();
	Mv.setZero();
}

inline void MPM_NODE::SlippingBC(Vector3d& n)
{	
	// make sure it's compressing the bc then set norm force and momentum to zero
	if (Mv.dot(n)>0.)
	{
		F  = F-F.dot(n)*n;
		Mv = Mv-Mv.dot(n)*n;
	}
}

inline void MPM_NODE::FrictionBC(double dt, Vector3d& n)
{
	double mvn = Mv.dot(n);							
	if (mvn>0.)										// 176
	{
		double fn = -mvn/dt;						// 182
		Vector3d mvtV = Mv - mvn*n;
		Vector3d t = mvtV.normalized();				// Tangential vector
		double mvt = mvtV.dot(t);
		double fs = -mvt/dt;						// 184
		double ff = Mu*abs(fn);
		double ft = -min(abs(ff), abs(fs));			// 185

		Vector3d mvb = Mv - F*dt;
		Vector3d mva = (mvt+ft*dt)*t;
		F = (mva - mvb)/dt;
		Mv = mva;
	}
}

#endif