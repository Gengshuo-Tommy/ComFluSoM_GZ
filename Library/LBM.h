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

#ifndef LBM_H
#define LBM_H

#include <HEADER.h>
#include <WEIGHT.h>
#include <IBM.h>
#include <DRAG_FORCE_MODEL.h>
#include <EFFECTIVE_VISCOSITY_MODEL.h>

// discrete model
enum DnQm
{
	D2Q5 ,
	D2Q9 ,
	D3Q7 ,
	D3Q15,
    D3Q19,
    D3Q27
};

enum CollisionModel
{
	SRT,
	MRT,
    CM
};

struct VelBC
{
	Vector3i Pos;
	Vector3i Nei;
	Vector3d Vel;
};

struct NoGradBC
{
	Vector3i Pos;
	Vector3i Nei;
};

class LBM
{
public:
	LBM();
	~LBM();
	LBM(DnQm dnqm, CollisionModel cmodel, size_t nx, size_t ny, size_t nz, double nu);
	void LESModel(VectorXd fneq, double omega, double& omegaT);
	void Init(double rho0, Vector3d initV);
	void InitSolidFraction();
	void InitPressureGradient();
	void InitCorrectMomentum();
	void Stream();
	void StreamGLBM();
	void SBounceBack(int i, int j, int k);
	void CalRhoVLocal(size_t i, size_t j, size_t k);
	void CalRhoVLocalNoForce(size_t i, size_t j, size_t k);
	void CalRhoVLocalNoForce(VectorXd f, double& rho, Vector3d& v);
	void CalRhoV();
	void CalRhoVNoForce();
	void SolveOneStep();
	void SolveOneStepPorous();
	void WriteFileH5(string fname, int n);
	Vector3d InterpolateVLiniear(const Vector3d& x);
	Vector3d InterpolateV3PD(const Vector3d& x);
	void SetEffectiveViscosityModel(double phiC, double iMiu);
	void UpdateEffectiveViscosity();
	/*===================================Functions for SRT=====================================================*/
	void  (LBM::*CalFeq)(VectorXd& feq, double rho, Vector3d v);							// Function pointer to calculate equilibrium distribution 
	double CalFeqCi(size_t q, double rho, Vector3d v);
	void CalFeqC(VectorXd& feq, double rho, Vector3d v);									// Function to calculate equilibrium distribution for compressible fluid
	void CalFeqPorous(VectorXd& feq, double rho, double eps, Vector3d v);
	void CollideSRTLocal(size_t i, size_t j, size_t k);
	void CollideSRTPSMLocal(size_t i, size_t j, size_t k, double g, Vector3d vw, Vector3d& fh);
	void CollideSRTPorousLocal(size_t i, size_t j, size_t k);
	void BodyForceSRTLocal(size_t i, size_t j, size_t k, Vector3d v, double omega, Vector3d force);
	void CollideSRT();
	void CollideSRTPorous();
	void UpdatePorousForce();
	void UpdatePorousForceHoef();
	void CalPGradLocal(size_t i, size_t j, size_t k);
	void CalPGrad();
	/*===================================Functions for MRT=====================================================*/
	void  (LBM::*CalMeq)(VectorXd& meq, double rho, Vector3d v);							// Function pointer to calculate equilibrium distribution 
	void CalMeqD2Q5(VectorXd& meq, double rho, Vector3d v);
	void CalMeqD2Q9(VectorXd& meq, double rho, Vector3d v);									// Function to calculate equilibrium distribution for compressible fluid
	void CalMeqD3Q7(VectorXd& meq, double rho, Vector3d v);
	void CalMeqD3Q15(VectorXd& meq, double rho, Vector3d v);
	void CalMeqD3Q19(VectorXd& meq, double rho, Vector3d v);
	void CollideMRTLocal(size_t i, size_t j, size_t k);
	void CollideMRTPorousLocal(size_t i, size_t j, size_t k);
	void BodyForceMRTLocal(size_t i, size_t j, size_t k, Vector3d v, MatrixXd mf, Vector3d force);
	void CollideMRT();
	void CollideMRTPorous();
	/*===================================Functions for Boundaries=====================================================*/
	void IBBYu(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh/*, bool refill*/);
	void HBB(size_t i, size_t j, size_t k, size_t q, Vector3d vw, Vector3d& fh/*, bool refill*/);
	void MBB(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh/*, bool refill*/);
	void VIBB(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh);
	void VIBB2(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh);
	void ApplyIBBFix();
	void SetFixedWall(string str);
	void ApplyWall();
	void ApplyVelBC();
	void ApplyNoGradBC();

	int 							Nproc;													// Number of processors which used

	double*** 						Rho;													// Fluid density
	double*** 						SolidFraction;											// Solid volume fraction
    double*** 						G;		             									// Gray coefficient
	double*** 						Fl;														// Liquid volume fraction
	double*** 						Source;													// Source term for CDE
	double***						Omega;													// Reversed Tua for variable viscosities
	int***		 					Flag0;													// Flag of lattice type
	int***		 					Flag;													// Flag of lattice type
	Vector3d***						V;														// Fluid velocity
	Vector3d***						ExForce;												// External force
	Vector3d***						PGrad;													// Pressure gradient
	Vector3d***						CMV;													// Correction of momentum

	Vector3d						Acc;													// Global Acceleration
	CollisionModel 					Cmodel;
    size_t 							Nx;														// Domain size
    size_t 							Ny;
    size_t 							Nz;
	int 							DomSize[3]; 

	bool  							UseLES;
	bool  							UseSolidFraction;
	bool  							ResetSolidFraction;
	bool  							UsePhaseChange;
	bool							Periodic[3];

	double  						Dx;														// Only used to convert to physical unit
	double  						Dt;

    double 							Nu0;													// Vsicosity
    double 							Tau0;													// Relaxation time                                                     
    double 							Omega0;													// Reciprocal of relaxation time
    double 							Rho0;
    double 							Sc2;													// Square of Smagorinsky constant
    // paramters for effective viscosity
    double  						PhiC;													// Critical solid volume fraction
    double  						IMiu;													// Intrinsic viscosity

    VectorXd						S;														// Relaxation matrix
    MatrixXd						M;														// Transform matrix
    MatrixXd						Mi;														// Inverse of transform matrix
    MatrixXd						Ms;														// Inverse of M multiply by S, for speed up.
    MatrixXd						Mf;														// Inverse of M multiply by (I-0.5*S) for force term, for speed up

    vector<size_t>					ViscoRT;												// Index of viscosity related component in S

    double							Dp;														// Particle diameter of porous media
    double							Rhop;													// Particle density of porous media

    size_t 							D;														// Dimension
    size_t 							Q;														// Number of discrete velocity vectors
    vector<Vector3d> 				E;														// Discrete velocities
    vector<Vector3i> 				Ne;														// Relative location of neighbor cells
    vector <double> 				W;														// Weights
    vector <int> 					Op;														// Opposite directions

    VectorXd*** 					F;														// Distribution function
    VectorXd*** 					Ft;														// Distribution function

    vector<tuple<size_t, size_t, size_t, size_t, double>>	Libb;							// List for ibb (lattice ID, q, delta)
    vector<Vector3i>				Lwall;													// List of wall nodes
    vector<VelBC>					Lvel;													// List of velocity BC nodes
    vector<NoGradBC>				Lnog;													// List of zero gradient BC nodes
};

inline LBM::LBM(DnQm dnqm, CollisionModel cmodel, size_t nx, size_t ny, size_t nz, double nu)
{
	Cmodel = cmodel;
	Nx	= nx;
	Ny	= ny;
	Nz	= nz;
	DomSize[0] = Nx;
	DomSize[1] = Ny;
	DomSize[2] = Nz;

	Periodic[0] = true;
	Periodic[1] = true;
	Periodic[2] = true;

	Nu0  = nu;
    Tau0 = 3.*Nu0 + 0.5;
    Omega0 = 1./Tau0;
    Sc2 = 0.;
    PhiC = 0.64;
    IMiu = 2.5;
    cout << "Relaxation time = " << Tau0 << endl;
    Nproc = 1;
	Acc = Vector3d::Zero();
	
	UseLES = false;
	UseSolidFraction = false;
	ResetSolidFraction = false;
	UsePhaseChange = false;
	ViscoRT.resize(0);
	Libb.resize(0);

	CalFeq = &LBM::CalFeqC;
	if (dnqm==D2Q5)
	{
		E   = {     { 0, 0, 0}, { 1, 0, 0}, { -1, 0, 0}, 
                	{0, 1, 0}, { 0,-1, 0}};
		W   = {     1./3. , 1./6. , 1./6. ,
                	1./6. , 1./6.};
		Op  = {     0, 2, 1, 
                	4, 3};
		D   = 2;
		Q   = 5;

		if (Cmodel==SRT)
		{
			cout << "Please use MRT for CDE." << endl;
			abort();
		}

	    if (Cmodel==MRT)
	    {
	    	S.resize(Q);
		    M.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0) <<  1,  1,  1,  1,  1;
		    M.row(1) <<  0,  1, -1,  0,  0;
		    M.row(2) <<  0,  0,  0,  1, -1;
		    M.row(3) <<  4, -1, -1, -1, -1;
		    M.row(4) <<  0,  1,  1, -1, -1;
		    // S values from "Lattice Boltzmann models for the convection-diffusion equation: D2Q5 vs D2Q9"
		    S(1) = S(2) = Omega0;
		    S(0) = S(3) = S(4) = 1.;

		    ViscoRT.push_back(1);
		    ViscoRT.push_back(2);

		    Mi = M.inverse();
		    Ms = Mi*S.asDiagonal();
		    VectorXd one(Q);
		    for (size_t q=0; q<Q; ++q) 	one(q) = 1.;
		    Mf = Mi*(one-0.5*S).asDiagonal()*M;
		    CalMeq = &LBM::CalMeqD2Q5;
		}
	}

	else if (dnqm==D2Q9)
	{
		E   = {     { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0}, 
                	{-1, 0, 0}, { 0,-1, 0}, { 1, 1, 0}, 
                	{-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
		W   = {     4./9. , 1./9. , 1./9. ,
                	1./9. , 1./9. , 1./36., 
                	1./36., 1./36., 1./36. };
		Op  = {     0, 3, 4, 
                	1, 2, 7, 
                	8, 5, 6 };
		Ne  = {     { 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
		D   = 2;
		Q   = 9;

	    if (Cmodel==MRT)
	    {
		    M.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0) <<  1,  1,  1,  1,  1,  1,  1,  1,  1;
		    M.row(1) << -4, -1, -1, -1, -1,  2,  2,  2,  2;
		    M.row(2) <<  4, -2, -2, -2, -2,  1,  1,  1,  1;
		    M.row(3) <<  0,  1,  0, -1,  0,  1, -1, -1,  1;
		    M.row(4) <<  0, -2,  0,  2,  0,  1, -1, -1,  1;
		    M.row(5) <<  0,  0,  1,  0, -1,  1,  1, -1, -1;
		    M.row(6) <<  0,  0, -2,  0,  2,  1,  1, -1, -1;
		    M.row(7) <<  0,  1, -1,  1, -1,  0,  0,  0,  0;
		    M.row(8) <<  0,  0,  0,  0,  0,  1, -1,  1, -1;

		    S.resize(Q);
		    S(0) = 1.;
		    S(1) = 1.4;
		    S(2) = 1.4;
		    S(3) = 1.;
		    S(4) = S(6) = 1.2;
		    S(5) = 1.;
		    S(7) = S(8) = Omega0;

		    ViscoRT.push_back(7);
		    ViscoRT.push_back(8);

		    Mi = M.inverse();
		    Ms = Mi*S.asDiagonal();
		    VectorXd one(Q);
		    for (size_t q=0; q<Q; ++q) 	one(q) = 1.;
		    Mf = Mi*(one-0.5*S).asDiagonal()*M;
		    CalMeq = &LBM::CalMeqD2Q9;
	    }
	}

	else if (dnqm==D3Q7)
	{
		E   = {     { 0, 0, 0}, 
					{ 1, 0, 0}, {-1, 0, 0}, 
                	{ 0, 1, 0}, { 0,-1, 0},
                	{ 0, 0, 1}, { 0, 0,-1}};
		W   = {     1./4. , 1./8. , 1./8. ,
                	1./8. , 1./8., 1./8. , 1./8.};
		Op  = {     0, 2, 1, 
                	4, 3, 6, 5};
		D   = 3;
		Q   = 7;

		if (Cmodel==SRT)
		{
			cout << "Please use MRT for CDE." << endl;
			abort();
		}
		// Modify alpha2 for D3Q7 "Multiple-relaxation-time lattice Boltzmann model for the convectionand anisotropic diffusion equation"
		// CalFeq = &LBM::CalGeq;
		Tau0 = 4.*Nu0+0.5;
		Omega0 = 1./Tau0;
		cout << "Tau0: " << Tau0 << endl;

	    if (Cmodel==MRT)
	    {
	    	S.resize(Q);
		    M.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0) <<  1,  1,  1,  1,  1,  1,  1;
		    M.row(1) <<  0,  1, -1,  0,  0,  0,  0;
		    M.row(2) <<  0,  0,  0,  1, -1,  0,  0;
		    M.row(3) <<  0,  0,  0,  0,  0,  1, -1;
		    M.row(4) <<  6, -1, -1, -1, -1, -1, -1;
		    M.row(5) <<  0,  2,  2, -1, -1, -1, -1;
		    M.row(6) <<  0,  0,  0,  1,  1, -1, -1;
		    // S values from "Multiple-relaxation-time lattice Boltzmann model for the convectionand anisotropic diffusion equation"
		    S(1) = S(2) = S(3)= Omega0;
		    S(0) = S(4) = S(5) = S(6) = 1.;

		    ViscoRT.push_back(1);
		    ViscoRT.push_back(2);
		    ViscoRT.push_back(3);

		    Mi = M.inverse();
		    Ms = Mi*S.asDiagonal();
		    VectorXd one(Q);
		    for (size_t q=0; q<Q; ++q) 	one(q) = 1.;
		    Mf = Mi*(one-0.5*S).asDiagonal()*M;
		    CalMeq = &LBM::CalMeqD3Q7;
		}		
	}

	else if (dnqm==D3Q15)
	{
		E   = {	   { 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
			       { 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
			       {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, {-1, 1, 1}, { 1,-1,-1} };
		W   = {	   2./9. , 1./9. , 1./9. , 1./9. , 1./9. ,
			       1./9. , 1./9. , 1./72., 1./72., 1./72.,
			       1./72., 1./72., 1./72., 1./72., 1./72. };
		Op  = {	   0,  2,  1,  4,  3, 
                   6,  5,  8,  7, 10, 
                   9, 12, 11, 14, 13 };
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
		D   = 3;
		Q   = 15;

	    if (Cmodel==MRT)
	    {
		    M.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0 ) <<  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(1 ) << -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(2 ) << 16.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(3 ) <<  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0;
		    M.row(4 ) <<  0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0;
		    M.row(5 ) <<  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0;
		    M.row(6 ) <<  0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0;
		    M.row(7 ) <<  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0;
		    M.row(8 ) <<  0.0,  0.0,  0.0,  0.0,  0.0, -4.0,  4.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0;
		    M.row(9 ) <<  0.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0;
		    M.row(10) <<  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0;
		    M.row(11) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0;
		    M.row(12) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0;
		    M.row(13) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0;
		    M.row(14) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0;

		    S.resize(Q);
		    S(0) = 0.;
		    S(1) = 1.6;
		    S(2) = 1.2;
		    S(3) = S(5) = S(7) = 0.;
		    S(4) = S(6) = S(8) = 1.6;
		    S(9) = S(10) = S(11) = S(12) = S(13) = Omega0;
		    S(14) = 1.2;

		    ViscoRT.push_back(9);
		    ViscoRT.push_back(10);
		    ViscoRT.push_back(11);
		    ViscoRT.push_back(12);
		    ViscoRT.push_back(13);

		    Mi = M.inverse();
		    Ms = Mi*S.asDiagonal();
		    VectorXd one(Q);
		    for (size_t q=0; q<Q; ++q) 	one(q) = 1.;
		    Mf = Mi*(one-0.5*S).asDiagonal()*M;

		    CalMeq = &LBM::CalMeqD3Q15;
	    }
	}

    else if (dnqm==D3Q19)
    {
        E   = {    { 0, 0, 0}, 
                   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
                   { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},
                   { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1} };
        W   = {    1./3. ,
                   1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 
                   1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 
                   1./36., 1./36., 1./36., 1./36., 1./36., 1./36. };
        Op  = {    0 , 
                   2 , 1 , 4 , 3 , 6 , 5 , 
                   8 , 7 , 10, 9 , 12, 11, 
                   14, 13, 16, 15, 18, 17 };
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
        D   = 3;
        Q   = 19;
        if (Cmodel==MRT)
        {
 		    M.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0 ) <<  1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.;
		    M.row(1 ) << -30., -11., -11., -11., -11., -11., -11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8.;
		    M.row(2 ) <<  12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.;
		    M.row(3 ) << 0., 1., -1., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0.;
		    M.row(4 ) << 0., -4., 4., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0.;
		    M.row(5 ) << 0., 0., 0., 1., -1., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1.;
		    M.row(6 ) << 0., 0., 0., -4., 4., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1.;
		    M.row(7 ) << 0., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1.;
		    M.row(8 ) << 0., 0., 0., 0., 0.,-4.,  4., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1.;
		    M.row(9 ) << 0., 2., 2., -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2.;
		    M.row(10) << 0., -4., -4., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2.;
		    M.row(11) << 0., 0., 0., 1., 1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0.;
		    M.row(12) << 0., 0., 0.,-2.,-2.,  2.,  2., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0.;
		    M.row(13) << 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., 0., 0., 0., 0.;
		    M.row(14) << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1.;
		    M.row(15) << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0.;
		    M.row(16) << 0., 0., 0., 0., 0., 0., 0., 1., -1., 1., -1., -1., 1., -1., 1., 0., 0., 0., 0.;
		    M.row(17) << 0., 0., 0., 0., 0., 0., 0., -1., -1., 1., 1., 0., 0., 0., 0., 1., -1., 1., -1.;
		    M.row(18) << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., -1., -1., 1., 1.;

		    S.resize(Q);
		    S(0) = 0.;
		    S(1) = 1.19;
		    S(2) = 1.4;
		    S(3) = S(5) = S(7) = 0.;
		    S(4) = S(6) = S(8) = 1.2;
		    S(9) = S(11) = Omega0;
		    S(10) = S(12) = 1.4;
		    S(13) = S(14) = S(15) = Omega0;	// related to bulk viscosity
		    S(16) = S(17) = S(18) = 1.98;

		    ViscoRT.push_back(9);
		    ViscoRT.push_back(11);
		    ViscoRT.push_back(13);
		    ViscoRT.push_back(14);
		    ViscoRT.push_back(15);

		    Mi = M.inverse();
		    Ms = Mi*S.asDiagonal();
		    VectorXd one(Q);
		    for (size_t q=0; q<Q; ++q) 	one(q) = 1.;
		    Mf = Mi*(one-0.5*S).asDiagonal()*M;

		    CalMeq = &LBM::CalMeqD3Q19;       	
        }
    }
    else if (dnqm==D3Q27)
    {
        E   = {    { 0, 0, 0},

                   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},

                   { 1, 1, 0}, {-1, 1, 0}, { 1,-1, 0}, {-1,-1, 0}, 
                   { 1, 0, 1}, {-1, 0, 1}, { 1, 0,-1}, {-1, 0,-1},
                   { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1},

                   { 1, 1, 1}, {-1, 1, 1}, { 1,-1, 1}, {-1,-1, 1},
                   { 1, 1,-1}, {-1, 1,-1}, { 1,-1,-1}, {-1,-1,-1} };
        W   = {    8./27., 
        		   2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 
        		   1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 
        		   1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.};
        Op  = {    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19};
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
        D   = 3;
        Q   = 27;
    }
    else
    {
        cout << "Undefined DnQm model: " << dnqm << endl;
        abort();
    }
}

inline double LBM::CalFeqCi(size_t q, double rho, Vector3d v)
{
    double vv  = v.dot(v);
    double ev  = E[q].dot(v);
    double fi = W[q]*rho*(1. + 3.*ev + 4.5*ev*ev - 1.5*vv);
    return fi;
}

inline void LBM::CalFeqC(VectorXd& feq, double rho, Vector3d v)
{
    double vv  = v.dot(v);
    for (size_t q=0; q<Q; ++q)
    {
    	double ev  = E[q].dot(v);
    	feq(q) = W[q]*rho*(1. + 3.*ev + 4.5*ev*ev - 1.5*vv);
	}
}

inline void LBM::CalFeqPorous(VectorXd& feq, double rho, double eps, Vector3d v)
{
	double vv  = v.dot(v);
	double epsi = 1./eps;
	for (size_t q=0; q<Q; ++q)
	{
		double ev  = E[q].dot(v);
		feq(q) = W[q]*rho*(1. + 3.*ev + (4.5*ev*ev - 1.5*vv)*epsi);
	}
}

inline void LBM::CalMeqD2Q5(VectorXd& meq, double rho, Vector3d v)
{
	meq(0) = rho;
	meq(1) = rho*v(0);
	meq(2) = rho*v(1);
	meq(3) = 2./3.*rho;
	meq(4) = 0.;
}

inline void LBM::CalMeqD2Q9(VectorXd& meq, double rho, Vector3d v)
{
	meq(0) = rho;
	meq(1) = rho*(-2 + 3*(v(0)*v(0) + v(1)*v(1)));
	meq(2) = -rho - meq(1);
	meq(3) = rho*v(0);
	meq(4) = -meq(3);
	meq(5) = rho*v(1);
	meq(6) = -meq(5);
	meq(7) = rho*(v(0)*v(0) - v(1)*v(1));
	meq(8) = rho*v(0)*v(1);
}

// "Lattice Boltzmann models for the convection-diffusion equation: D2Q5 vs D2Q9"
inline void LBM::CalMeqD3Q7(VectorXd& meq, double rho, Vector3d v)
{
	meq(0) = rho;
	meq(1) = rho*v(0);
	meq(2) = rho*v(1);
	meq(3) = rho*v(2);
	meq(4) = 3./4.*rho;
	meq(5) = meq(6) = 0.;
}

inline void LBM::CalMeqD3Q15(VectorXd& meq, double rho, Vector3d v)
{
	double v2 = v.squaredNorm();
	meq(0) = rho;
	meq(1) = rho*(v2-1.);
	meq(2) = rho*(1.-5*v2);
	meq(3) = rho*v(0);
	meq(4) = -7.*rho*v(0)/3.;
	meq(5) = rho*v(1);
	meq(6) = -7.*rho*v(1)/3.;
	meq(7) = rho*v(2);
	meq(8) = -7.*rho*v(2)/3.;
	meq(9) = rho*(3.*v(0)*v(0)-v2);
	meq(10)= rho*(v(1)*v(1)-v(2)*v(2));
	meq(11)= rho*v(0)*v(1);
	meq(12)= rho*v(1)*v(2);
	meq(13)= rho*v(0)*v(2);
	meq(14)= 0.;
}

inline void LBM::CalMeqD3Q19(VectorXd& meq, double rho, Vector3d v)
{
	double v2 = v.squaredNorm();
	meq(0) = rho;
	meq(1) = rho*(-11.+19.*v2);
	meq(2) = rho*(3.+-5.5*v2);
	meq(3) = rho*v(0);
	meq(4) = meq(3)*(-2./3.);
	meq(5) = rho*v(1);
	meq(6) = meq(5)*(-2./3.);
	meq(7) = rho*v(2);
	meq(8) = meq(7)*(-2./3.);
	meq(9) = rho*(3.*v(0)*v(0)-v2);
	meq(10)= -0.5*meq(9);
	meq(11)= rho*(v(1)*v(1)-v(2)*v(2));
	meq(12)= -0.5*meq(11);
	meq(13)= rho*v(0)*v(1);
	meq(14)= rho*v(1)*v(2);
	meq(15)= rho*v(0)*v(2);
	meq(16)=meq(17)=meq(18)=0.;
}

inline void LBM::Init(double rho0, Vector3d initV)
{
    cout << "================ Start init. ================" << endl;
    Rho0 	= rho0;
	Rho		= new double**  [Nx+1];
	Omega	= new double**  [Nx+1];
	G		= new double**	[Nx+1];
	Flag0 	= new int** 	[Nx+1];
	Flag 	= new int** 	[Nx+1];
	V		= new Vector3d**[Nx+1];
	ExForce	= new Vector3d**[Nx+1];
	F		= new VectorXd**[Nx+1];
	Ft		= new VectorXd**[Nx+1];

	for (size_t i=0; i<=Nx; ++i)
	{
		Rho    [i]	= new double*  [Ny+1];
		Omega  [i]	= new double*  [Ny+1];
		G	   [i]  = new double*  [Ny+1];
		Flag0  [i]  = new int* 	   [Ny+1];
		Flag   [i]  = new int* 	   [Ny+1];
		V      [i]	= new Vector3d*[Ny+1];
		ExForce[i]	= new Vector3d*[Ny+1];
		F      [i]	= new VectorXd*[Ny+1];
		Ft     [i]	= new VectorXd*[Ny+1];

		for (size_t j=0; j<=Ny; j++)
		{
			Rho    [i][j]	= new double  [Nz+1];
			Omega  [i][j]	= new double  [Nz+1];
			G      [i][j]	= new double  [Nz+1];
			Flag0  [i][j]	= new int 	  [Nz+1];
			Flag   [i][j]	= new int 	  [Nz+1];
			V      [i][j]	= new Vector3d[Nz+1];
			ExForce[i][j]	= new Vector3d[Nz+1];
			F      [i][j]	= new VectorXd[Nz+1];
			Ft     [i][j]	= new VectorXd[Nz+1];

			for (size_t k=0; k<=Nz; k++)
			{
				Rho    [i][j][k]	= Rho0;
				Omega  [i][j][k]	= Omega0;
				G      [i][j][k]	= 0;
				Flag0  [i][j][k]    = -1;
				Flag   [i][j][k]    = -1;
				V      [i][j][k]	= initV;
				ExForce[i][j][k]	= Vector3d::Zero();

				VectorXd feq(Q);
				if (Cmodel==SRT)
				{
					(this->*CalFeq)(feq, rho0, initV);
				}
				else if (Cmodel==MRT)
				{
					VectorXd meq(Q);
					(this->*CalMeq)(meq, rho0, initV);
					feq = Mi*meq;
				}

				F [i][j][k]	= feq;
				Ft[i][j][k]	= feq;
			}
		}
	}
	cout << "================ Finish init. ================" << endl;
}

inline void LBM::InitSolidFraction()
{
    cout << "================ Start init G. ================" << endl;
	SolidFraction		= new double**  [Nx+1];

	for (size_t i=0; i<=Nx; ++i)
	{
		SolidFraction[i]	= new double*  [Ny+1];

		for (size_t j=0; j<=Ny; j++)
		{
			SolidFraction[i][j]	= new double  [Nz+1];

			for (size_t k=0; k<=Nz; k++)
			{
				SolidFraction[i][j][k]	= 0.;
			}
		}
	}
	UseSolidFraction = true;
	cout << "================ Finish init G. ================" << endl;
}

inline void LBM::InitPressureGradient()
{
    cout << "================ Start init Pressure Gradient. ================" << endl;
	PGrad		= new Vector3d**  [Nx+1];

	for (size_t i=0; i<=Nx; ++i)
	{
		PGrad[i]	= new Vector3d*  [Ny+1];

		for (size_t j=0; j<=Ny; j++)
		{
			PGrad[i][j]	= new Vector3d  [Nz+1];

			for (size_t k=0; k<=Nz; k++)
			{
				PGrad[i][j][k] = Vector3d::Zero();
			}
		}
	}
	cout << "================ Finish init Pressure Gradient. ================" << endl;
}

inline void LBM::InitCorrectMomentum()
{
    cout << "================ Start init Pressure Gradient. ================" << endl;
	CMV		= new Vector3d**  [Nx+1];

	for (size_t i=0; i<=Nx; ++i)
	{
		CMV[i]	= new Vector3d*  [Ny+1];

		for (size_t j=0; j<=Ny; j++)
		{
			CMV[i][j]	= new Vector3d  [Nz+1];

			for (size_t k=0; k<=Nz; k++)
			{
				CMV[i][j][k] = Vector3d::Zero();
			}
		}
	}
	cout << "================ Finish init Pressure Gradient. ================" << endl;
}

void LBM::SetEffectiveViscosityModel(double phiC, double iMiu)
{
	PhiC = phiC;
	IMiu = iMiu;
}

inline void LBM::LESModel(VectorXd fneq, double omega, double& omegaT)
{
	double qHat2 = 0.;
	for (size_t m=0; m<D; ++m)
	for (size_t n=0; n<D; ++n)
	{
		double qij = 0.;
		for (size_t q=0; q<Q; ++q)	qij += E[q](m)*E[q](n)*fneq(q);
		qHat2 += qij*qij;
	}
	double qHat = sqrt(2.*qHat2);
	double tau = 1./omega;
	double tauEddy = 0.5*(sqrt(tau*tau+18.*Sc2*qHat)-tau);
	omegaT = 1./(tau+tauEddy);
}

inline void LBM::CollideSRTLocal(size_t i, size_t j, size_t k)
{
	double rho = Rho[i][j][k];
	Vector3d v = V[i][j][k];
	VectorXd feq(Q), fneq(Q);
	(this->*CalFeq)(feq, rho, v);
	fneq = F[i][j][k] - feq;
	double omega = Omega[i][j][k];
	double omegaT = omega;
	if (UseLES)	LESModel(fneq, omega, omegaT);
	Ft[i][j][k] = F[i][j][k] - omegaT*fneq;
	// apply external force
	Vector3d fb = rho*Acc + ExForce[i][j][k];
	BodyForceSRTLocal(i, j, k, v, omegaT, fb);
	ExForce[i][j][k].setZero();
}

inline void LBM::CollideSRTPorousLocal(size_t i, size_t j, size_t k)
{
	double rho = Rho[i][j][k];
	Vector3d v = V[i][j][k];
	VectorXd feq(Q), fneq(Q);
	double poro = 1.-SolidFraction[i][j][k];
	CalFeqPorous(feq, rho, poro, v);
	fneq = F[i][j][k] - feq;
	double omega = Omega[i][j][k];
	Ft[i][j][k] = F[i][j][k] - omega*fneq;
	// apply external force
	Vector3d fb = ExForce[i][j][k] + poro*rho*Acc;
	Vector3d vf = v/poro;
	BodyForceSRTLocal(i, j, k, vf, omega, fb);
	ExForce[i][j][k].setZero();
}

inline void LBM::CollideSRTPSMLocal(size_t i, size_t j, size_t k, double g, Vector3d vw, Vector3d& fh)
{
	double rho = Rho[i][j][k];
	Vector3d v = V[i][j][k];
	VectorXd feq(Q), fneq(Q), feqv(Q), collideS(Q);
	(this->*CalFeq)(feq, rho, v);
	(this->*CalFeq)(feqv, rho, vw);
	fneq = F[i][j][k] - feq;
	double omega = Omega[i][j][k];
	double omegaT = omega;
	if (UseLES)	LESModel(fneq, omega, omegaT);
	double bn = g*(1./omegaT-0.5)/((1.-g) + (1./omegaT-0.5));

	for (size_t q=0; q<Q; ++q)
	{
		size_t op = Op[q];
		collideS(q) = F[i][j][k](op) - feqv(op) - F[i][j][k](q) + feqv(q);
	}

	Ft[i][j][k] = F[i][j][k] - (1.-bn)*omegaT*fneq + bn*collideS;
	// apply external force
	Vector3d fb = (1.-bn)*(rho*Acc + ExForce[i][j][k]);
	BodyForceSRTLocal(i, j, k, v, omegaT, fb);
	ExForce[i][j][k].setZero();
}

// Lattice Boltzmann model for incompressible flows through porous media
// Generalized lattice Boltzmann model for flow through tight porous media with Klinkenberg’s effect
inline void LBM::BodyForceSRTLocal(size_t i, size_t j, size_t k, Vector3d v, double omega, Vector3d force)
{
	double a = (1.-0.5*omega);
	double vf = 3.*v.dot(force);
	for (size_t q=0; q<Q; ++q)
	{
		double ef = E[q].dot(force);
		double ev = E[q].dot(v);
		#pragma omp atomic
		Ft[i][j][k](q) += W[q]*a*(3.*ef*(1.+3.*ev)-vf);

		if (Ft[i][j][k](q)<0.)
		{
			cout << "negative after force term!" << endl;
			cout << "force: " << force.transpose() << endl;
			cout << "i: " << i << " j: " << j << " k: " << k << endl;
			abort();
		}
	}
}

inline void LBM::CalPGradLocal(size_t i, size_t j, size_t k)
{
	Vector3d pg = Vector3d::Zero();
	for (size_t q=0; q<Q; ++q)
	{
		int in = (i+ (int) E[q][0]+Nx+1)%(Nx+1);
		int jn = (j+ (int) E[q][1]+Ny+1)%(Ny+1);
		int kn = (k+ (int) E[q][2]+Nz+1)%(Nz+1);

		// pg += W[q]*E[q]*(Rho[in][jn][kn]/*-Rho0*/)/*/(1.-SolidFraction[in][jn][kn])*/;
		pg += W[q]*E[q]*(Rho[in][jn][kn]- Rho0)/(1.-SolidFraction[in][jn][kn]);
	}
	PGrad[i][j][k] = pg;
}

inline void LBM::CollideMRTLocal(size_t i, size_t j, size_t k)
{
	double rho = Rho[i][j][k];
	Vector3d v = V[i][j][k];
	VectorXd meq(Q), m(Q), mneq(Q);
	(this->*CalMeq)(meq, rho, v);
	m = M*F[i][j][k];
	mneq = m-meq;
	Vector3d fb = rho*Acc + ExForce[i][j][k];
	// For LES subgrid model
	if (UseLES)
	{
		VectorXd fneq(Q);
		fneq = Mi*mneq;
		double omega = Omega[i][j][k];
		double omegaT = omega;
		LESModel(fneq, omega, omegaT);

		VectorXd s = S;
		for (size_t n=0; n<ViscoRT.size(); ++n)		s(ViscoRT[n]) = omegaT;
		Ft[i][j][k] = F[i][j][k] - Mi*s.asDiagonal()*mneq;
		MatrixXd mf(Q,Q);
		mf.setIdentity();
		mf -= 0.5*Mi*s.asDiagonal()*M;
		BodyForceMRTLocal(i, j, k, v, mf, fb);
	}
	else
	{
		Ft[i][j][k] = F[i][j][k] - Ms*mneq;
		BodyForceMRTLocal(i, j, k, v, Mf, fb);
	}
	ExForce[i][j][k].setZero();
}

inline void LBM::CollideMRTPorousLocal(size_t i, size_t j, size_t k)
{
	double rho = Rho[i][j][k];
	Vector3d v = V[i][j][k];
	VectorXd feq(Q), meq(Q), m(Q), mneq(Q);
	double poro = 1.-SolidFraction[i][j][k];
	CalFeqPorous(feq, rho, poro, v);
	m = M*F[i][j][k];
	meq = M*feq;
	mneq = m-meq;

	double omegaT = Omega[i][j][k];
	VectorXd s = S;
	for (size_t n=0; n<ViscoRT.size(); ++n)		s(ViscoRT[n]) = omegaT;
	Ft[i][j][k] = F[i][j][k] - Mi*s.asDiagonal()*mneq;

	Vector3d fb = ExForce[i][j][k] + poro*rho*Acc;
	MatrixXd mf(Q,Q);
	mf.setIdentity();
	mf -= 0.5*Mi*s.asDiagonal()*M;
	Vector3d vf = v/poro;
	BodyForceMRTLocal(i, j, k, vf, mf, fb);

	ExForce[i][j][k].setZero();
}

inline void LBM::BodyForceMRTLocal(size_t i, size_t j, size_t k, Vector3d v, MatrixXd mf, Vector3d force)
{
	VectorXd fb(Q), fbm(Q);
	double vf = 3.*v.dot(force);
	for (size_t q=0; q<Q; ++q)
	{
		double ef = E[q].dot(force);
		double ev = E[q].dot(v);
		fb(q) = W[q]*(3.*ef*(1.+3.*ev)-vf);
	}
	fbm = mf*fb;
	for (size_t q=0; q<Q; ++q)
	{
		#pragma omp atomic
		Ft[i][j][k](q) += fbm(q);
	}
}

inline void LBM::CalRhoVLocal(size_t i, size_t j, size_t k)
{
    double rho = 0.;
    Vector3d v;
    v.setZero();
    for (size_t q = 0; q < Q; ++q)
    {
    	rho	+= F[i][j][k](q);
    	v 	+= F[i][j][k](q)*E[q];
    }
    Rho[i][j][k] = rho;
    V[i][j][k] = (v+0.5*ExForce[i][j][k])/rho + 0.5*Acc;
}

inline void LBM::CalRhoVLocalNoForce(size_t i, size_t j, size_t k)
{
    double rho = 0.;
    Vector3d v;
    v.setZero();
    for (size_t q = 0; q < Q; ++q)
    {
    	rho	+= F[i][j][k](q);
    	v 	+= F[i][j][k](q)*E[q];
    }
    Rho[i][j][k] = rho;
    V[i][j][k] = v/rho;
}

inline void LBM::CalRhoVLocalNoForce(VectorXd f, double& rho, Vector3d& v)
{
	rho = 0.;
	v.setZero();
    for (size_t q = 0; q < Q; ++q)
    {
    	rho	+= f(q);
    	v 	+= f(q)*E[q];
    }
    v /= rho;
}

inline void LBM::CalRhoV()
{
	// double err = 0.;
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CalRhoVLocal(i,j,k);
		Flag[i][j][k] = Flag0[i][j][k];
		Flag0[i][j][k] = -1;
    }
}

inline void LBM::CalRhoVNoForce()
{
	// double err = 0.;
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CalRhoVLocalNoForce(i,j,k);
		Flag[i][j][k] = Flag0[i][j][k];
		Flag0[i][j][k] = -1;
    }
}

inline void LBM::CollideSRT()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CollideSRTLocal(i, j, k);
    }
}

inline void LBM::UpdateEffectiveViscosity()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
    	double phi = SolidFraction[i][j][k];
		double visRatio = CalEffectiveViscosity(PhiC, IMiu, phi);
		double tau = 3.*Nu0*visRatio+0.5;
		Omega[i][j][k] = 1./tau;
    }	
}

inline void LBM::CalPGrad()
{
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
    	CalPGradLocal(i, j, k);
    }
}

inline void LBM::UpdatePorousForce()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
    	CalRhoVLocalNoForce(i,j,k);

    	double rho = Rho[i][j][k];
    	double eps = 1.-SolidFraction[i][j][k];
		if (eps==1.)
		{
			// ExForce[i][j][k] = rho*Acc;
			V[i][j][k] += 0.5*Acc;
		}
		else
		{
			Vector3d v0 = V[i][j][k];
			Vector3d v = v0 + 0.5*eps*Acc;							// Eq.17
			double fe = 1.75/(sqrt(150.*eps)*eps);					// Eq. 3	geometric function
			double perm = Dp*Dp*eps/(150.*pow(1./eps-1.,2));		// Eq. 4	permeability
			double c0 = 0.5*(1.+0.5*eps*Nu0/perm);					// Eq.18
			double c1 = 0.5*eps*fe/sqrt(perm);						// Eq.18
			Vector3d vc = v/(c0+sqrt(c0*c0 + c1*v.norm()));			// Eq.16
			V[i][j][k] = vc;
			ExForce[i][j][k] = 2.*rho*(vc-v0) - rho*eps*Acc;		// Eq.15
		}
		// if (ResetSolidFraction)	SolidFraction[i][j][k] = 0.;
    }
}

inline void LBM::UpdatePorousForceHoef()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
    	Vector3d vt = V[i][j][k];

    	CalRhoVLocalNoForce(i,j,k);

    	double rho = Rho[i][j][k];
    	double phi = SolidFraction[i][j][k];
    	double eps = 1.-phi;
		if (eps==1.)
		{
			V[i][j][k] += 0.5*Acc;
		}
		else
		{
			double nu = (1./Omega[i][j][k]-0.5)/3.;
			Vector3d vs (0., 0.,0.);
			Vector3d ve = V[i][j][k];
			Vector3d vc, fe;
			VelDragForceHoef(rho, nu, phi, Dp, vt, vs, Acc, ve, vc, fe);
			V[i][j][k] = vc;
			ExForce[i][j][k] = fe;		// Eq.15
		}
		// if (ResetSolidFraction)	SolidFraction[i][j][k] = 0.;
    }
}

inline void LBM::CollideSRTPorous()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CollideSRTPorousLocal(i, j, k);
    }
}

inline void LBM::CollideMRT()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CollideMRTLocal(i, j, k);
    }
}

inline void LBM::CollideMRTPorous()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
    {
		CollideMRTPorousLocal(i, j, k);
    }
}

inline void LBM::Stream()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
	{
		for (size_t q=0; q< Q;  ++q)
		{
			int ip = (i- (int) E[q][0]+Nx+1)%(Nx+1);
			int jp = (j- (int) E[q][1]+Ny+1)%(Ny+1);
			int kp = (k- (int) E[q][2]+Nz+1)%(Nz+1);

			F[i][j][k][q] = Ft[ip][jp][kp][q];
		}
		if (ResetSolidFraction)	SolidFraction[i][j][k] = 0.;
	}
}

inline void LBM::StreamGLBM()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i = 0; i <= Nx; ++i)
    for (size_t j = 0; j <= Ny; ++j)
    for (size_t k = 0; k <= Nz; ++k)
	{
		double ns =  G[i][j][k];
		for (size_t q=0; q< Q;  ++q)
		{
			int ip = (i- (int) E[q][0]+Nx+1)%(Nx+1);
			int jp = (j- (int) E[q][1]+Ny+1)%(Ny+1);
			int kp = (k- (int) E[q][2]+Nz+1)%(Nz+1);

			F[i][j][k](q) = (1.-ns)*Ft[ip][jp][kp][q] + ns*Ft[ip][jp][kp](Op[q]);
		}
		if (ResetSolidFraction)	SolidFraction[i][j][k] = 0.;
	}
}

inline void LBM::SBounceBack(int i, int j, int k)
{
	VectorXd ft(Q);
	ft = F[i][j][k];
	for (size_t q=0; q<Q; ++q)
	{
		Ft[i][j][k](q) 	= ft(Op[q]);
	}
}

inline void LBM::IBBYu(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh)
{
	size_t op = Op[q];
	double fw = delta*Ft[i][j][k](q) + (1-delta)*F[i][j][k](q) + 6.*W[q]*Rho[i][j][k]*E[op].dot(vw);
	F[i][j][k](op) = fw/(1+delta) + delta*Ft[i][j][k](op)/(1+delta);

	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	fh = Ft[i][j][k](q)*(E[q]-vw) + F[i][j][k](op)*(E[q]+vw);
}

inline void LBM::HBB(size_t i, size_t j, size_t k, size_t q, Vector3d vw, Vector3d& fh)
{
	size_t op = Op[q];
	F[i][j][k](op) = Ft[i][j][k](q) + 6.*W[q]*Rho0*E[op].dot(vw);
	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	fh = Ft[i][j][k](q)*(E[q]-vw) + F[i][j][k](op)*(E[q]+vw);
}

inline void LBM::VIBB(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh)
{
	size_t op = Op[q];
	size_t in = (i + (int) E[op](0)+Nx+1)%(Nx+1);
	size_t jn = (j + (int) E[op](1)+Ny+1)%(Ny+1);
	size_t kn = (k + (int) E[op](2)+Nz+1)%(Nz+1);

	double rhof = Rho[i][j][k];
	double rhoff = Rho[in][jn][kn];
	Vector3d vf = V[i][j][k];
	Vector3d vff = V[in][jn][kn];

	double rhod = 2.*delta*rhof + (1.-2.*delta)*rhoff;

	if (delta>0.5) rhod = rhof;

	double fneqf = Ft[i][j][k](q)-CalFeqCi(q, rhof, vf);
	double fneqff = Ft[in][jn][kn](q)-CalFeqCi(q, rhoff, vff);
	double fneqd = 2.*delta*fneqf*/*rhod/rhof*/ + (1.-2.*delta)*(fneqff)/**rhod/rhoff*/;

	Vector3d vd0;
	if (delta<0.5)	vd0 = 2.*delta*vf + (1.-2.*delta)*vff;
	else  			vd0 = ((1.-delta)*vf + (2.*delta-1.)*vw)/delta;
	Vector3d vd1 = ((1.-delta)*vff + 2.*delta*vw)/(1.+delta);
	Vector3d vd = (vd0 + 2.*vd1)/3.;

	// Vector3d vd = vd0;

	double feqd = CalFeqCi(q, rhod, vd);

	F[i][j][k](op) = feqd +fneqd + 6.*W[q]*Rho[i][j][k]*E[op].dot(vw);

	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	// fh = F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);
	// fh = (Ft[i][j][k](q) + F[i][j][k](op))*E[q];
	// if (refill)	fh = (F[i][j][k](q) + F[i][j][k](op))*E[q];
	fh = Ft[i][j][k](q)*(E[q]-vw) + F[i][j][k](op)*(E[q]+vw);
	// if (refill)	fh = F[i][j][k](q)*(E[q]-vw) + F[i][j][k](op)*(E[q]+vw);
	// if (refill)	fh.setZero();
}

inline void LBM::VIBB2(size_t i, size_t j, size_t k, size_t q, double delta, Vector3d vw, Vector3d& fh)
{
	size_t op = Op[q];
	size_t in = (i + (int) E[op](0)+Nx+1)%(Nx+1);
	size_t jn = (j + (int) E[op](1)+Ny+1)%(Ny+1);
	size_t kn = (k + (int) E[op](2)+Nz+1)%(Nz+1);

	double rhof = Rho[i][j][k];
	double rhoff = Rho[in][jn][kn];
	Vector3d vf = V[i][j][k];
	Vector3d vff = V[in][jn][kn];

	double rhod = 2.*delta*rhof + (1.-2.*delta)*rhoff;

	double fneqf = Ft[i][j][k](q)-CalFeqCi(q, rhof, vf);
	double fneqff = Ft[in][jn][kn](q)-CalFeqCi(q, rhoff, vff);
	double fneqd = 2.*delta*fneqf*rhod/rhof + (1.-2.*delta)*(fneqff)*rhod/rhoff;

	Vector3d vd0;
	if (delta<0.5)	vd0 = 2.*delta*vf + (1.-2.*delta)*vff;
	else  			vd0 = ((1.-delta)*vf + (2.*delta-1.)*vw)/delta;
	Vector3d vd1 = ((1.-delta)*vff + 2.*delta*vw)/(1.+delta);
	Vector3d vd = (vd0 + 2.*vd1)/3.;

	double feqd = CalFeqCi(q, rhod, vd);

	F[i][j][k](op) = feqd +fneqd + 6.*W[q]*Rho[i][j][k]*E[op].dot(vw);

	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	fh = Ft[i][j][k](q)*(E[q]-vw) + F[i][j][k](op)*(E[q]+vw);
}

inline void LBM::SetFixedWall(string str)
{
	if (str=="Z_bot")
	{
		size_t k = 0;
		for (size_t i=0; i<=Nx; ++i)
		for (size_t j=0; j<=Ny; ++j)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](2)<0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}
	else if (str=="Z_top")
	{
		size_t k = Nz;
		for (size_t i=0; i<=Nx; ++i)
		for (size_t j=0; j<=Ny; ++j)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](2)>0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}

	else if (str=="Y_bot")
	{
		size_t j = 0;
		for (size_t i=0; i<=Nx; ++i)
		for (size_t k=0; k<=Nz; ++k)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](1)<0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}
	else if (str=="Y_top")
	{
		size_t j = Ny;
		for (size_t i=0; i<=Nx; ++i)
		for (size_t k=0; k<=Nz; ++k)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](1)>0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}

	else if (str=="X_bot")
	{
		size_t i = 0;
		for (size_t j=0; j<=Ny; ++j)
		for (size_t k=0; k<=Nz; ++k)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](0)<0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}
	else if (str=="X_top")
	{
		size_t i = Nx;
		for (size_t j=0; j<=Ny; ++j)
		for (size_t k=0; k<=Nz; ++k)
		{
			for (size_t q=1; q<Q; ++q)
			{
				if (E[q](0)>0.)	Libb.push_back({i,j,k,q, 0.5});
			}
		}
	}
}

inline void LBM::ApplyWall()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<Lwall.size(); ++c)
	{
		SBounceBack(Lwall[c](0),Lwall[c](1),Lwall[c](2));
	}
}

inline void LBM::ApplyIBBFix()
{
	Vector3d vw (0., 0., 0.);
	Vector3d fh (0., 0., 0.);	
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t m=0; m<Libb.size(); ++m)
	{
		Vector3d fhi;
		size_t i = get<0>(Libb[m]);
		size_t j = get<1>(Libb[m]);
		size_t k = get<2>(Libb[m]);
		size_t q = get<3>(Libb[m]);
		double delta = get<4>(Libb[m]);
		IBBYu(i, j, k, q, delta, vw, fhi);
		// HBB(i, j, k, q, vw, fhi);
		fh += fhi;
	}
}

inline void LBM::ApplyVelBC()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<Lvel.size(); ++c)
	{
		int i = Lvel[c].Pos(0);
		int j = Lvel[c].Pos(1);
		int k = Lvel[c].Pos(2);

		int in = Lvel[c].Nei(0);
		int jn = Lvel[c].Nei(1);
		int kn = Lvel[c].Nei(2);
		double rhon = Rho[in][jn][kn];
		Vector3d vn = V[in][jn][kn];

		VectorXd feq(Q), feqn(Q);
		(this->*CalFeq)(feq, rhon, Lvel[c].Vel);
		(this->*CalFeq)(feqn, rhon, vn);

		Ft[i][j][k] = feq + (F[in][jn][kn]-feqn);
		CalRhoVLocal(i,j,k);
	}
}

inline void LBM::ApplyNoGradBC()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<Lnog.size(); ++c)
	{
		int i = Lnog[c].Pos(0);
		int j = Lnog[c].Pos(1);
		int k = Lnog[c].Pos(2);

		int in = Lnog[c].Nei(0);
		int jn = Lnog[c].Nei(1);
		int kn = Lnog[c].Nei(2);

		F[i][j][k] = F[in][jn][kn];
	}	
}

// Linear interpolate velocity at any poit x
inline Vector3d LBM::InterpolateVLiniear(const Vector3d& x)
{
    Vector3d vx (0., 0., 0.);
    // Linear interpolation
    // Find node
    Vector3i min0, max0;
    min0(0) = floor(x(0));
    max0(0) = min0(0)+1;
    min0(1) = floor(x(1));
    max0(1) = min0(1)+1;
    min0(2) = floor(x(2));
    max0(2) = min0(2)+1;

    vector <Vector3d> ver;
    ver.resize(8);

    ver[0] << min0(0), min0(1), min0(2);
    ver[1] << max0(0), min0(1), min0(2);
    ver[2] << max0(0), max0(1), min0(2);
    ver[3] << min0(0), max0(1), min0(2);

    ver[4] << max0(0), min0(1), max0(2);
    ver[5] << min0(0), min0(1), max0(2);
    ver[6] << min0(0), max0(1), max0(2);
    ver[7] << max0(0), max0(1), max0(2);

    for (size_t l=0; l<pow(2,D); ++l)
    {
        int i = (int) ver[l](0);
        int j = (int) ver[l](1);
        int k = (int) ver[l](2);
        i = (i+Nx+1)%(Nx+1);
        j = (j+Ny+1)%(Ny+1);
        k = (k+Nz+1)%(Nz+1);
        Vector3d s = x-ver[7-l];
        double vol = abs(s(0)*s(1)*s(2));
        vx += vol*V[i][j][k];
    }
    return vx;
}

inline Vector3d LBM::InterpolateV3PD(const Vector3d& x)
{
    Vector3d vx (0., 0., 0.);
    // Find node
    Vector3i min0 (0,0,0);
    Vector3i max0 (0,0,0);

    for(size_t d=0; d<D; ++d)
    {
    	min0(d) = floor(x(d)-1.55);
    	max0(d) = floor(x(d)+1.55);
    }

    for (int k=min0(2); k<=max0(2); ++k)
    for (int j=min0(1); j<=max0(1); ++j)
    for (int i=min0(0); i<=max0(0); ++i)
    {
    	Vector3d xn (i,j,k);
    	double s = (x-xn).norm();
    	vx += Delta3P(s)*V[i][j][k];
    }
    return vx;
}

// inline void LBM::DistributeForce3PD(const Vector3d& x, const Vector3d& force)
// {

// }

inline void LBM::SolveOneStep()
{
	CalRhoV();
	if (Cmodel==SRT)		CollideSRT();
	else if (Cmodel==MRT)	CollideMRT();
	ApplyWall();
	ApplyVelBC();
	Stream();
	ApplyIBBFix();
	// ApplyVelBC();
	ApplyNoGradBC();
}

inline void LBM::SolveOneStepPorous()
{
	UpdatePorousForce();
	// UpdatePorousForceHoef();
	// CalRhoV();
	// ApplyVelBC();
	// ApplyNoGradBC();
	// CollideSRTPorous();
	if (Cmodel==SRT)		CollideSRTPorous();
	else if (Cmodel==MRT)	CollideMRTPorous();
	ApplyWall();
	Stream();
	ApplyIBBFix();
	ApplyVelBC();
	ApplyNoGradBC();
}

inline void LBM::WriteFileH5(string fname, int n)
{
	stringstream	out;					//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = fname+"_LBM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.

	hsize_t nx = Nx+1;
	hsize_t ny = Ny+1;
	hsize_t nz = Nz+1;

	int numLat = nx*ny*nz;

	hsize_t	dims_scalar[3] = {nx, ny, nz};			//create data space.
	hsize_t	dims_vector[4] = {nx, ny, nz, 3};		//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* rho_h5 	= new double[numLat];
	double* g_h5 	= new double[numLat];
	double* vel_h5 	= new double[3*numLat];

	string rho_name = "Density";
	string g_name = "Gamma";

	int len = 0;
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	if (UseSolidFraction)
	{
		g_name = "SolidFraction";
		for (size_t k=0; k<nz; k++)
		for (size_t j=0; j<ny; j++)
		for (size_t i=0; i<nx; i++)
		{
			rho_h5[len] 	= Rho[i][j][k];
			g_h5[len] 		= SolidFraction[i][j][k];
			double phi = 1.-g_h5[len];
			vel_h5[3*len  ] = V[i][j][k](0)/phi;
			vel_h5[3*len+1] = V[i][j][k](1)/phi;
			vel_h5[3*len+2] = V[i][j][k](2)/phi;
	        len++;
		}
	}
	else if (UsePhaseChange)
	{
		rho_name = "Enthalpy";
		g_name = "LiquidFraction";
		for (size_t k=0; k<nz; k++)
		for (size_t j=0; j<ny; j++)
		for (size_t i=0; i<nx; i++)
		{
			rho_h5[len] 	= Rho[i][j][k];
			g_h5[len] 		= Fl[i][j][k];
			vel_h5[3*len  ] = V[i][j][k](0);
			vel_h5[3*len+1] = V[i][j][k](1);
			vel_h5[3*len+2] = V[i][j][k](2);
	        len++;
		}
	}
	else
	{
		for (size_t k=0; k<nz; k++)
		for (size_t j=0; j<ny; j++)
		for (size_t i=0; i<nx; i++)
		{
			rho_h5[len] 	= Rho[i][j][k];
			g_h5[len] 		= Flag0[i][j][k];
			vel_h5[3*len  ] = V[i][j][k](0);
			vel_h5[3*len+1] = V[i][j][k](1);
			vel_h5[3*len+2] = V[i][j][k](2);
	        len++;
		}
	}

	DataSet	*dataset_rho = new DataSet(file.createDataSet(rho_name, PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_g	 = new DataSet(file.createDataSet(g_name, PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete dataset_rho;
	delete dataset_g;
	delete dataset_vel;

	delete rho_h5;
	delete g_h5;
	delete vel_h5;

	file.close();

	string file_name_xmf = fname+"_LBM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"LBM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << nz << " " << ny << " " << nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\""<< rho_name << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/" << rho_name << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\""<< g_name << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/" << g_name << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}

#endif
