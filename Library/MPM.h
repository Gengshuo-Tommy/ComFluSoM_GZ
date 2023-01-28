/************************************************************************
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

#ifndef MPM_H
#define MPM_H

#include <HEADER.h>
#include <SHAPE.h>
#include <MPM_PARTICLE.h>
#include <MPM_NODE.h>
#include <MPM_CELL.h>
// #include "../Mesh/SURFACE_MESH.h"
// #include "../Mesh/READOBJ.h"

class MPM
{
public:
	MPM(): gen(std::random_device()()) {};
	MPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx);
	void Init(bool useFbar);
	void UpdateLn(MPM_PARTICLE* p0);
	void CalNGN(MPM_PARTICLE* p0);
	void CalNGN_MLS(MPM_PARTICLE* p0);
	void UpdateLAn();
	void ParticleToNode();
	void ParticleToNodePoro();
	void ParticleToCell();
	void ParticleToNodeMLS();
	void CalFOnNode(bool firstStep);
	void UpdateParticle(MPM_PARTICLE* p0);
	void UpdateStress(MPM_PARTICLE* p0);
	void NodeToParticle();
	void NodeToParticleDoubleMapping();
	void CalVGradLocal(int p);
	void CalPSizeCP(MPM_PARTICLE* p0);
	void CalPSizeR(MPM_PARTICLE* p0);
	void CalVOnNode();
	void CalSurfaceTension();
	void CalVGradOnNode();
	void CalVOnNodeDoubleMapping();
	void SetNonSlippingBC(size_t n);
	void SetNonSlippingBC(size_t i, size_t j, size_t k);
	void SetSlippingBC(size_t n, Vector3d& norm);
	void SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm);
	void SetFrictionBC(size_t n, double mu, Vector3d& norm);
	void SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm);
	void SolveMUSL(int tt, int ts);
	void SolveUSF(int tt, int ts);
	void SolveUSA(int tt, int ts);
	void AddNode(size_t level, Vector3d& x);
	void AddParticle(int tag, Vector3d& x, double m);
	void DeleteParticles();
	void AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m);
	void WriteFileH5(int n);
	void FindIndex(size_t n, size_t& i, size_t& j, size_t& k);
	double GetUniformD1();
	void LoadMPMFromH5(string fname, double ratio);

	double 		(*N)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d 	(*GN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	void 			(*NGN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);

	vector <size_t>					LAn;														// List of actived nodes
	vector <size_t>					LAc;														// List of actived nodes
	vector <MPM_PARTICLE*>			Lp;														// List of all MPM particles
	vector <MPM_PARTICLE*>			Lbp;														// List of boundary MPM particles
	vector <MPM_NODE*>				Ln;														// List of all MPM nodes
	vector <MPM_CELL*>				Lcell;													// List of all MPM cells

	double  								CoefSurf;												// Coefficient of surface tension
	double  								RhoH;														// Density of heavy phase
	double  								RhoL;														// Density of light phase

	vector <Vector3d>					Nei;														// Relative position of neighbor nodes
	vector <Vector3d>					NeiG;														// Shape function gradient of neighbor nodes

	bool									Periodic[3];
	bool 									MLSv;
	bool									ActiveShift;
	bool  								UseFbar;

   size_t 								Nx;														// Domain size
   size_t 								Ny;
   size_t 								Nz;
   size_t 								Ncz;
   size_t 								Ncy;
   size_t 								Nnode;													// Total number of nodes

   size_t 								Nproc;
   size_t 								D;															// Dimension	
   size_t 								Ntype;													// Type of shape function 0 for Linear, 1 for Quadratic and 2 for Cubic 3 for GIMP

   double 								Nrange;													// Influence range of shape function
   double 								Dt;														// Time step
   double 								Dc;														// Damping coefficient
   double 								Cs;														// Speed of sound
   Vector3d								Dx;														// Space step
   Vector3d 							Shift;													// Shift for Galilean invariant
   mt19937 					   		gen;														// Random number generator
};

MPM::MPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx)
{
	Nproc	= 1;
	Ntype 	= ntype;
	Nx 		= nx;
	Ny 		= ny;
	Nz 		= nz;
	Dx 		= dx;
	D 			= 3;
	Dt 		= 1.;
	Dc 		= 0.;
	Cs 		= 0.;
	Periodic[0] = false;
	Periodic[1] = false;
	Periodic[2] = false;

	MLSv 			= false;
	ActiveShift	= false;
	UseFbar		= false;

	Ncz = (Nx+1)*(Ny+1);
	Ncy = (Nx+1);

	Nnode = (Nx+1)*(Ny+1)*(Nz+1);

	Shift.setZero();

	if (Nz==0)
	{
		D = 2;
		if (Ny==0)	D = 1;
	}
	// Linear
	if 		(Ntype == 0)
	{
		if (D==1) 			NGN =& LS1D;
		else if (D==2)		NGN =& LS2D;
		else if (D==3)		NGN =& LS3D;
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.;
		cout << "Using Linear shape function." << endl;
	}
	else if (Ntype == 2)
	{
		
	}
	// GIMP
	else if (Ntype == 3)
	{
		if (D==1)			NGN =& GIMP1D;
		else if (D==2) 	NGN =& GIMP2D;
		else if (D==3)		NGN =& GIMP3D;
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.;
		cout << "Using GIMP shape function." << endl;
	}
	else
	{
		cout << "Undefined shape function type. Retry 0 for Linear, 1 for Quadratic and 2 for Cubic." << endl;
		abort();
	}

	if (D==3)
	{
		Nei  = {   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
	}
	else if (D==2)
	{
		Nei  = {	{ 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
	}

	// for surface tension
	Vector3d x0 (0.,0.,0.);
	for (size_t i=0; i<Nei.size(); ++i)
	{
		Vector3d xn = Nei[i];
		double n = 0.;
		Vector3d gn;
		Vector3d psize (0.5, 0.5, 0.5);
		if (D==2)	psize(2) = 0.; 
		NGN(x0, xn, Dx, psize, n, gn);
		NeiG.push_back(gn);
	}
	CoefSurf = 0.;
	RhoH = 0.;
	RhoL = 0.;
}

// Uniform distribution from -0.5 to 0.5
inline double MPM::GetUniformD1()
{
	uniform_real_distribution<double>	dis(-0.5,0.5);
	return dis(gen);
}

void MPM::Init(bool useFbar)
{
	cout << "================ Start init.  ================" << endl;
	Lp.resize(0);
	Ln.resize(0);
	Lcell.resize(0);
	// Add basic nodes
	for (size_t n=0; n<(Nx+1)*(Ny+1)*(Nz+1); ++n)
	{
    	size_t i, j, k;
    	FindIndex(n, i, j, k);
		Vector3d x (i, j, k);
		Ln.push_back(new MPM_NODE(0,x));
    	Ln[Ln.size()-1]->ID = Ln.size()-1;
	}

	UseFbar = useFbar;
	if (UseFbar)
	{
		for (size_t n=0; n<(Nx+1)*(Ny+1)*(Nz+1); ++n)
		{
	    	size_t i, j, k;
	    	FindIndex(n, i, j, k);
			Vector3d x (i+0.5, j+0.5, k+0.5);
			Lcell.push_back(new MPM_CELL(0,x));
	    	Lcell[Lcell.size()-1]->ID = Lcell.size()-1;
		}
	}
	cout << "=============== Finish init.  ================" << endl;
}
// Find index for grid
inline void MPM::FindIndex(size_t n, size_t& i, size_t& j, size_t& k)
{
	k = n/Ncz;
	j = (n%Ncz)/Ncy;
	i = (n%Ncz)%Ncy;
}

void MPM::CalNGN(MPM_PARTICLE* p0)
{
	// Reset shape function (N) and gradient of shape function (GN)
	p0->Lni.resize(0);
	p0->LnN.resize(0);
	p0->LnGN.resize(0);
	// Find min position of nodes which is infuenced by this particle
	Vector3i minx 	= Vector3i::Zero();
	Vector3i maxx 	= Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		maxx(d) = (int) trunc(p0->X(d) + p0->PSize(d) + 1.);
		minx(d) = (int) ceil(p0->X(d) - p0->PSize(d) - 1.);
	}
	// Find nodes within the influence range
	for (int i=minx(0); i<=maxx(0); ++i)
	for (int j=minx(1); j<=maxx(1); ++j)
	for (int k=minx(2); k<=maxx(2); ++k)
	{
		double n;
		Vector3d gn;
		Vector3d xn (i,j,k);
		NGN(p0->X, xn, Dx, p0->PSize, n, gn);
		// Find id of current node
		size_t ii = (i+Nx+1)%(Nx+1);
		size_t jj = (j+Ny+1)%(Ny+1);
		size_t kk = (k+Nz+1)%(Nz+1);
		size_t id = ii+jj*Ncy+kk*Ncz;
		if (n>0.)
		{
			p0->Lni.push_back(id);
			p0->LnN.push_back(n);
			p0->LnGN.push_back(gn);
		}
	}
}

void MPM::ParticleToNode()
{
	// reset mass internal force velocity for nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Ln[id]->Reset();
    }
    LAn.resize(0);

    // Update shape function and grad
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	MPM_PARTICLE* p0 = Lp[p];
    	// Apply linear spring to boudnary particles and pull them back to intial position
    	if (p0->FixBySpring)		p0->Fc = p0->Kn*(p0->X0-p0->X)/**sqrt((Lp[p]->X0-Lp[p]->X).norm())*/;
    	
    	p0->X += Shift;
    	CalNGN(p0);
    	Matrix3d vsp = -p0->Vol*p0->Stress;
		Vector3d fex = p0->M*p0->B + p0->Fh + p0->Fc;

		for (size_t l=0; l<p0->Lni.size(); ++l)
		{
			// Grid id
			size_t id = p0->Lni[l];
			// weight
			double 		n 		= p0->LnN[l];
			Vector3d 	gn 	= p0->LnGN[l];
			Vector3d 	df 	= n*fex + vsp*gn;
			// weigthed mass contribution
			double nm = n*p0->M;

			Vector3d ve = p0->V;
			// Vector3d ve = p0->V-Lp[p]->L*(Lp[p]->X-Ln[id]->X);
			#pragma omp atomic
			Ln[id]->M += nm;
			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->Mv(d) += nm*ve(d);
				#pragma omp atomic
				Ln[id]->F(d) += df(d);
/*				// smooth stress on node
				for (size_t c=0; c<D; ++c)
				{
					#pragma omp atomic
					Ln[id]->Stress(d,c) += nm*p0->Stress(d,c);
				}*/
			}
		}
    }
    UpdateLAn();
}

void MPM::ParticleToNodePoro()
{
	// reset mass internal force velocity for nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Ln[id]->Reset();
   }
   LAn.resize(0);

   // Update shape function and grad
   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t p=0; p<Lp.size(); ++p)
   {
   	CalNGN(Lp[p]);
    	double   poro = Lp[p]->Poro;
    	Matrix3d vsp  = -Lp[p]->Vol*Lp[p]->Stress;
		Vector3d fex  = (1-poro)*Lp[p]->M*Lp[p]->B + Lp[p]->Fh;
		double   pex  = (1-poro)*Lp[p]->Vol*Lp[p]->P;

		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 		n 	= Lp[p]->LnN[l];
			Vector3d 	gn = Lp[p]->LnGN[l];
			Vector3d 	df = n*fex + pex*gn + vsp*gn;
			// weigthed mass contribution
			double nm  = n*Lp[p]->M*(1-Lp[p]->Poro);
			double npm = n*Lp[p]->M;
			if (nm<0.)
			{
				cout << "mass problem" << endl;
				cout << "nm= " << nm << endl;
				cout << "p= " << p << endl;
				cout << "id= " << id << endl;
				cout << "Lp[p]->X= " << Lp[p]->X.transpose() << endl;
				cout << "Ln[id]->X= " << Ln[id]->X.transpose() << endl;
				abort();
			}
			Vector3d ve = Lp[p]->V-Lp[p]->L*(Lp[p]->X-Ln[id]->X);
			#pragma omp atomic
			Ln[id]->M  += nm;
			#pragma omp atomic
			Ln[id]->Mp += npm;
			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->Mv(d) += nm*ve(d);
				#pragma omp atomic
				Ln[id]->F(d) += df(d);
				// #pragma omp atomic
				// Ln[id]->Mv(d) += df(d)*Dt;
				// smooth stress on node
				for (size_t c=0; c<D; ++c)
				{
					#pragma omp atomic
					Ln[id]->Stress(d,c) += nm*Lp[p]->Stress(d,c);
				}
			}
		}
	}

   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t p=0; p<Lp.size(); ++p)
   {
   	for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
   	{
   		size_t id = Lp[p]->Lni[l];
   		Ln[id]->Poro = 1 - Ln[id]->M/Ln[id]->Mp;
   	}
   }

   UpdateLAn();
}


// need to improve efficiency
void MPM::ParticleToCell()
{
   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t n=0; n<Lcell.size(); ++n)
   {
   	Lcell[n]->Reset();
   }

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		MPM_PARTICLE* p0 = Lp[p];
		int i = ceil(p0->X(0))-1;
		int j = ceil(p0->X(1))-1;
		int k = ceil(p0->X(2))-1;
		// Find id of current cell
		size_t ii = (i+Nx+1)%(Nx+1);
		size_t jj = (j+Ny+1)%(Ny+1);
		size_t kk = (k+Nz+1)%(Nz+1);
		size_t cid = ii+jj*Ncy+kk*Ncz;
		p0->CID = cid;

		Vector3d xcp = Lcell[cid]->X-p0->X;
		if (abs(xcp(0))>0.5 || abs(xcp(1))>0.5)
		{
			cout << "Lcell[cid]->X: " << Lcell[cid]->X.transpose() << endl;
			cout << "p0->X: " << p0->X.transpose() << endl;
			abort();
		}

		#pragma omp atomic
		Lcell[cid]->Vol0 += p0->Vol0;
		#pragma omp atomic
		Lcell[cid]->Vol += p0->Vol0*p0->Td.determinant();

		// #pragma omp atomic
		// Lcell[cid]->Np++;

		// for (size_t d=0; d<D; ++d)
		// for (size_t c=0; c<D; ++c)
		// {
		// 	#pragma omp atomic
		// 	Lcell[cid]->Stress(d,c) += p0->Stress(d,c);
		// }		
	}

   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t n=0; n<Lcell.size(); ++n)
   {
   	Lcell[n]->J = Lcell[n]->Vol/Lcell[n]->Vol0;
   	// cout << "j: " << Lcell[n]->J << endl;
   	// Lcell[n]->Stress /= (double) Lcell[n]->Np;
   }
}

void MPM::CalNGN_MLS(MPM_PARTICLE* p0)
{
	// Reset shape function (N) and gradient of shape function (GN)
	p0->Lni.resize(0);
	p0->Lgi.resize(0);
	p0->LnN.resize(0);
	p0->LnGN.resize(0);
	// Find min position of nodes which is infuenced by this particle
	// for nodes
	Vector3i minn 	= Vector3i::Zero();
	Vector3i maxn 	= Vector3i::Zero();
	// for gauss points
	Vector3i ming 	= Vector3i::Zero();
	Vector3i maxg 	= Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		// minn(d) = ceil(p0->X(d) -1.);
		// maxn(d) = trunc(p0->X(d) +1.);
		minn(d) = ceil(p0->X(d) -1.5);
		maxn(d) = trunc(p0->X(d) +1.5);
		ming(d) = ceil(p0->X(d)-2.);
		maxg(d) = trunc(p0->X(d)+1.);
	}

	// Find nodes within the influence range
	for (int i=minn(0); i<=maxn(0); ++i)
	for (int j=minn(1); j<=maxn(1); ++j)
	for (int k=minn(2); k<=maxn(2); ++k)
	{
		double n;
		Vector3d gn;
		// Find id of current node
		int id = i+j*Ncy+k*Ncz;
		// if (id>Ln.size())
		// {
		// 	cout << "id too large" << endl;
		// 	cout << p0->X.transpose() << endl;
		// 	cout << i << " " << j << " " << k << endl;
		// 	cout << id << endl;
		// 	cout << p0->ID << endl;
		// 	cout << Lp.size() << endl;
		// }
		NGN(p0->X, Ln[id]->X, Dx, p0->PSize, n, gn);
		if (n>0.)
		{
			p0->Lni.push_back(id);
			p0->LnN.push_back(n);
			p0->LnGN.push_back(gn);
			Ln[id]->Actived = true;
		}
		// store particle ID for MLS
		#pragma omp critical
		{
			Ln[id]->MPs.push_back(p0->ID);
		}
	}
	// Find nodes within the influence range
	// for (int i=ming(0); i<=maxg(0); ++i)
	// for (int j=ming(1); j<=maxg(1); ++j)
	// for (int k=ming(2); k<=maxg(2); ++k)
	// {
	// 	// Find id of current node
	// 	int id = i+j*Ncy+k*Ncz;
	// 	if (id>Ln.size())
	// 	{
	// 		cout << "id too large GP" << endl;
	// 		cout << p0->X.transpose() << endl;
	// 		cout << i << " " << j << " " << k << endl;
	// 		cout << id << endl;
	// 		cout << p0->ID << endl;
	// 	}
	// 	// store particle ID for MLS
	// 	#pragma omp critical
	// 	{
	// 		Lg[id]->MPs.push_back(p0->ID);
	// 		// if (Lg[id]->MPs.size()>50)
	// 		// {
	// 		// 	cout << "Lg[id]->MPs too large" << endl;
	// 		// 	cout << Lg[id]->MPs.size() << endl;
	// 		// 	abort();
	// 		// }
	// 	}
	// 	p0->Lgi.push_back(id);
	// }
}

void MPM::ParticleToNodeMLS()
{
	// reset mass internal force velocity for nodes
	// #pragma omp parallel for schedule(static) num_threads(1)
	// for (size_t n=0; n<LAn.size(); ++n)
	// {
	// 	size_t id = LAn[n];
	// 	Ln[id]->Reset();
 //    }
    for (size_t n=0; n<Ln.size(); ++n)
    {
    	Ln[n]->Reset();
    }
    LAn.resize(0);

    // Update shape function and grad
    // #pragma omp parallel for schedule(static) num_threads(1)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	Lp[p]->X += Shift;
    	// UpdateLn(Lp[p]);
    	CalNGN_MLS(Lp[p]);

    	Matrix3d vsp = -Lp[p]->Vol*Lp[p]->Stress;
		Vector3d fex = Lp[p]->M*Lp[p]->B + Lp[p]->Fh + Lp[p]->Fc;

		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 		n 	= Lp[p]->LnN[l];
			Vector3d 	gn 	= Lp[p]->LnGN[l];
			Vector3d 	df 	= n*fex + vsp*gn;
			// weigthed mass contribution
			double nm = n*Lp[p]->M;

			if (nm<0.)
			{
				cout << "mass problem" << endl;
				cout << "nm= " << nm << endl;
				cout << "p= " << p << endl;
				cout << "id= " << id << endl;
				cout << "Lp[p]->X= " << Lp[p]->X.transpose() << endl;
				cout << "Ln[id]->X= " << Ln[id]->X.transpose() << endl;
				abort();
			}
			#pragma omp atomic
			Ln[id]->M += nm;
			#pragma omp atomic
			Ln[id]->F(0) += df(0);
			#pragma omp atomic
			Ln[id]->F(1) += df(1);
			#pragma omp atomic
			Ln[id]->F(2) += df(2);
		}
    }
    UpdateLAn();
    // cout << "start mls" << endl;
	// #pragma omp parallel for schedule(static) num_threads(1)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		size_t nmp = Ln[id]->MPs.size();
		bool wrong = false;
		if (nmp>800)
		{
			cout << "nmp: " << nmp << endl;
			wrong = true;
		}
		// position of MPs
		vector<Vector3d> Xp;
		// velocity of MPs
		MatrixXd Vp(nmp,3);
		for (size_t i=0; i<nmp; ++i)
		{
			// ID of MP
			size_t p = Ln[id]->MPs[i];
			Xp.push_back(Lp[p]->X);
			Vp.row(i) = Lp[p]->V.transpose();
			if (wrong)
			{
				cout << p << endl;
			}
		}
		// shape functions
		// VectorXd phi = MLS(Xp, Ln[id]->X, 1);
		VectorXd phi = MLS(Xp, Ln[id]->X, 0);
		// // node velocity
		Ln[id]->V = phi.transpose()*Vp;
		Ln[id]->Mv = Ln[id]->M*Ln[id]->V+Ln[id]->F;
	}
	// cout << "end mls" << endl;
}

void MPM::UpdateLAn()
{
	vector<vector <size_t>> lan;
	lan.resize(Nproc);
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		auto id = omp_get_thread_num();
		lan[id].insert( lan[id].end(), Lp[p]->Lni.begin(), Lp[p]->Lni.end() );
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Nproc; ++n)
	{
		sort( lan[n].begin(), lan[n].end() );
		lan[n].erase(unique(lan[n].begin(), lan[n].end()), lan[n].end());
	}
	for (size_t n=0; n<Nproc; ++n)
	{
		LAn.insert( LAn.end(), lan[n].begin(), lan[n].end() );
	}
	sort( LAn.begin(), LAn.end() );
	LAn.erase(unique(LAn.begin(), LAn.end()), LAn.end());
}

void MPM::CalSurfaceTension()
{

	double mav = 0.5*(RhoH+RhoL);

	vector<size_t> lb(0);	// bounary nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id0 = LAn[n];
		Vector3d x0 = Ln[id0]->X;				// current node position
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;

			if ((Ln[id]->M-mav)*(Ln[id0]->M-mav) < 0.)	
			{
				Ln[id0]->IsBoundary = true;
				#pragma omp critical
				{
					lb.push_back(id0);
				}
				break;
			}
		}
	}

	// lb = LAn;

	// calculate mbar for all node, should improve the efficiency latter
	vector<double> lmbar0(Ln.size());
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Ln.size(); ++n)
	{
		Vector3d x0 = Ln[n]->X;
		double mbar = 0.;
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;
			mbar += Ln[id]->M;
		}
		lmbar0[n] = mbar / (double) Nei.size();
	}
	// second iteration for smooth
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Ln.size(); ++n)
	{
		Vector3d x0 = Ln[n]->X;
		double mbar = 0.;
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;
			mbar += lmbar0[id];
		}
		Ln[n]->Mbar = mbar / (double) Nei.size();
	}

	// find the surrounding nodes that need to calcualte their mbar gradient
	vector<size_t> ls(0);
	for (size_t i=0; i<lb.size(); ++i)
	{
		Vector3d x0 = Ln[lb[i]]->X;
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;
			ls.push_back(id);
		}
	}
   sort( ls.begin(), ls.end() );
   ls.erase( unique( ls.begin(), ls.end() ), ls.end() );

   // ls = LAn;

   // calculate gradient of mbar
   for (size_t i=0; i<ls.size(); ++i)
   {
   	Vector3d x0 = Ln[ls[i]]->X;
   	Vector3d gm;
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;
			gm -= NeiG[m]*Ln[id]->Mbar;
		}
		Ln[ls[i]]->GradM = gm;
		gm.normalize();
		Ln[ls[i]]->SurfNorm = gm;
   }

   // calculate 
   for (size_t i=0; i<lb.size(); ++i)
   {
   	size_t id0 = lb[i];
		Vector3d x0 = Ln[id0]->X;
		double kc = 0.;
		for (size_t m=0; m<Nei.size(); ++m)
		{
			Vector3d xn = x0 + Nei[m];
			// Find id of current node
			size_t ii = ((int) xn(0)+Nx+1)%(Nx+1);
			size_t jj = ((int) xn(1)+Ny+1)%(Ny+1);
			size_t kk = ((int) xn(2)+Nz+1)%(Nz+1);
			size_t id = ii+jj*Ncy+kk*Ncz;
			kc += NeiG[m].dot(Ln[id]->SurfNorm);
		}   	
		Vector3d fsur = CoefSurf*kc*Ln[id0]->GradM/(RhoL*RhoL-RhoH*RhoH)/0.5*Ln[id0]->M;
		Ln[id0]->F += fsur;
		// // cout << "fsur: " << fsur.transpose() << endl;
		// cout << "GradM: " << Ln[id0]->GradM.transpose() << endl;
		// cout << "kc: " << kc << endl;
		// cout << "Ln[id0]->SurfNorm: " << Ln[id0]->SurfNorm.transpose() << endl;
		// // cout << "CoefSurf: " << CoefSurf << endl;
		// cout << "Ln[id0]->M: " << Ln[id0]->M << endl;
		// // cout << "RhoH*RhoH-RhoL*RhoL: " << RhoH*RhoH-RhoL*RhoL << endl;
		// cout << "x: " << Ln[id0]->X.transpose() << endl;
		// abort();

		// if (x0(0)==31. && x0(1)==31.)
		// {
		// 	cout << "x: " << Ln[id0]->X.transpose() << endl;
		// 	cout << "fsur: " << fsur.transpose() << endl;
		// 	cout << "!!!!!!!!!!!!!!" << endl;
		// 	abort();
		// }

		// if (x0(0)>20. && x0(1)>20.)
		// {
		// 	cout << "x: " << Ln[id0]->X.transpose() << endl;
		// 	cout << "fsur: " << fsur.transpose() << endl;
		// 	cout << "Ln[id0]->M: " << Ln[id0]->M << endl;
		// 	// abort();
		// }
   }
   // abort();
}

void MPM::CalVOnNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Ln[id]->Stress /= Ln[id]->M;
		if (Ln[id]->M==0.)	Ln[id]->V.setZero();
		else
		{
			Ln[id]->F  *= (1.-Dc*Sign(Ln[id]->Mv.dot(Ln[id]->F)));
			Ln[id]->Mv += Ln[id]->F*Dt;

			if (Ln[id]->BCTypes.size()>0)
			{
				for (size_t i=0; i<Ln[id]->BCTypes.size(); ++i)
				{
					if (Ln[id]->BCTypes[i]==1)				Ln[id]->NonSlippingBC();
					else if (Ln[id]->BCTypes[i]==2)		Ln[id]->SlippingBC(Ln[id]->Norms[i]);
					else if (Ln[id]->BCTypes[i]==3)		Ln[id]->FrictionBC(Dt, Ln[id]->Norms[i]);
				}
			}
			// Ln[id]->Mv += Ln[id]->F*Dt;
			Ln[id]->V = Ln[id]->Mv/Ln[id]->M;
			Ln[id]->F /= Ln[id]->M;
		}
	}
}

void MPM::CalVGradOnNode()
{
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Vector3i vin (1, Ncy, Ncz);
		for (size_t d=0; d<D; ++d)
		{
			size_t idL = id-vin(d);
			size_t idR = id+vin(d);
			Ln[id]->VGrad.row(d) = 0.5*(Ln[idR]->V - Ln[idL]->V);
			cout << "Ln[idR]->V: " << Ln[idR]->V.transpose() << endl;
			cout << "Ln[idL]->V: " << Ln[idL]->V.transpose() << endl;
		}
			cout << "++++++++++++" << endl;
			cout << Ln[id]->VGrad << endl;
			abort();
	}
}

void MPM::SetNonSlippingBC(size_t n)
{
	Ln[n]->BCTypes.push_back(1);
}
void MPM::SetNonSlippingBC(size_t i, size_t j, size_t k)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(1);
}

void MPM::SetSlippingBC(size_t n, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(2);
	Ln[n]->Norms.push_back(norm);
}
void MPM::SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	// cout << "n= " << n << endl;
	// cout << "Nnode= " << Nnode << endl;
	Ln[n]->BCTypes.push_back(2);
	// cout << "push_back 1 " << endl;
	Ln[n]->Norms.push_back(norm);
	// cout << "push_back 2 " << endl;
}

void MPM::SetFrictionBC(size_t n, double mu, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;

}
void MPM::SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}

void MPM::UpdateParticle(MPM_PARTICLE* p0)
{
	p0->StressSmooth.setZero();
	// p0->L = Lfd;
	// Vector3d xp = p0->X;
	p0->X -= Shift;

	if (!p0->FixV)
	{
		// Vector3d xb = p0->X;
		for (size_t l=0; l<p0->Lni.size(); ++l)
		{
			size_t 	id = p0->Lni[l];
			double 	n  = p0->LnN[l];
			Vector3d vb = p0->V;
			// Update velocity of this particle
			p0->V += n*Ln[id]->F*Dt;				// force is divided by nodal mass already

			#pragma omp critical
			{
				if (p0->V.norm()>0.5)
				{
					cout << "vb: " << vb.transpose() << endl;
					cout << "p0->V: " << p0->V.transpose() << endl;
					cout << "p0->X: " << p0->X.transpose() << endl;
					for (size_t l=0; l<p0->Lni.size(); ++l)
					{
						size_t 	id = p0->Lni[l];
						// double 	n  = p0->LnN[l];
						Vector3d an = Ln[id]->F/Ln[id]->M;
						cout << "an: " << an.transpose() << endl;
						cout << "Ln[id]->F: " << Ln[id]->F.transpose() << endl;
						cout << "Ln[id]->Mv: " << Ln[id]->Mv.transpose() << endl;
						cout << "Ln[id]->M: " << Ln[id]->M << endl;			
					}
					abort();
				}
			}

			// p0->X += n*Ln[id]->Mv/Ln[id]->M*Dt;
			p0->X += n*Ln[id]->V*Dt;
			p0->StressSmooth += n*Ln[id]->Stress;
		}
	}
	else
	{
		p0->V = p0->Vf;
		p0->X += p0->V*Dt;
	}

	if (p0->X(0)<0.) 				p0->X(0) += Nx+1.;
	else if (p0->X(0)>Nx+1.) 	p0->X(0) -= Nx+1.;

	if (p0->X(1)<0.) 				p0->X(1) += Ny+1.;
	else if (p0->X(1)>Ny+1.) 	p0->X(1) -= Ny+1.;

	p0->L = Matrix3d::Zero();

	for (size_t l=0; l<p0->Lni.size(); ++l)
	{
		size_t id = p0->Lni[l];
		Vector3d gn = p0->LnGN[l];
		// Calculate velocity gradient tensor
		p0->L += Ln[id]->V*gn.transpose();
		// p0->L += gn*Ln[id]->V.transpose();
	}
	// F-bar method
	if (UseFbar)
	{
		double bar = pow(Lcell[p0->CID]->J/p0->Td.determinant(), 1./(double) D);
		p0->Td *= bar;
		// p0->Stress = Lcell[p0->CID]->Stress;
	}
	p0->Td = (Matrix3d::Identity() + p0->L*Dt)*p0->Td;
}

void MPM::UpdateStress(MPM_PARTICLE* p0)
{
	// p0->Td = (Matrix3d::Identity() + p0->L*Dt)*p0->Td;


	// Update particle length
	if (Ntype==3)
	{
		if (p0->Type==0)	CalPSizeR(p0);
	}
	// else 				CalPSizeCP(p);
	// else if (Lp[p]->Type==2 || Lp[p]->Type==3)	CalPSizeCP(p);
	// else if (Lp[p]->Type==2 || Lp[p]->Type==3)	CalPSizeR(p);
	// Update volume of particles
	p0->Vol = p0->Td.determinant()*p0->Vol0;
	// update porosity of particles
	p0->Poro = 1 - (1 - p0->PoroInit)/p0->Td.determinant();
	// Update strain
	Matrix3d de = 0.5*Dt*(p0->L + p0->L.transpose());
	Matrix3d dw = 0.5*Dt*((p0->L - p0->L.transpose()));
	// Update stress
	p0->Stress += dw*p0->Stress - p0->Stress*dw.transpose();
	// p0->Stress += dw*p0->Stress - p0->Stress*dw;
	// p0->Vol 	*= 1.+de.trace();
	if (p0->Type==0)			p0->Elastic(de);
	else if (p0->Type==1)
	{
		p0->EOSMonaghan(Cs);
		p0->Newtonian(de);
	}
	else if (p0->Type==2)	p0->MohrCoulomb(de);
	else if (p0->Type==3)	p0->DruckerPrager(de);
	else if (p0->Type==5)
	{
		// Lp[p]->EOSMonaghan(Cs);
		p0->EOSMorris(Cs);
		p0->Granular(de,dw);
	}
	// Reset hydro force and contact force
	p0->Fh.setZero();
	p0->Fc.setZero();
	p0->P = 0;
}

// void MPM::NodeToParticle()
// {
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		// cout << "p" << p << endl;
// 		// Reset position increasement of this particle
// 		Lp[p]->DeltaX 	= Vector3d::Zero();
// 		Lp[p]->StressSmooth.setZero();
//  		if (!Lp[p]->FixV)
// 		{
// 			for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
// 			{
// 				size_t 	id = Lp[p]->Lni[l];
// 				double 	n  = Lp[p]->LnN[l];
// 				// Update velocity of this particle
// 				Vector3d an = Ln[id]->F/Ln[id]->M;				
// 				Lp[p]->V += n*an*Dt;
// 				Lp[p]->X += n*Ln[id]->V*Dt;
// 				Lp[p]->StressSmooth += n*Ln[id]->Stress;
// 				// Lp[p]->X += n*Ln[id]->Mv/Ln[id]->M*Dt;
// 			}
// 		}
// 		else
// 		{
// 			Lp[p]->V = Lp[p]->Vf;
// 			Lp[p]->X += Lp[p]->V*Dt;
// 		}
// 		// CalVOnNodeDoubleMapping();
// 		// Velocity gradient tensor
// 		CalVGradLocal(p);
// 		Lp[p]->Td = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->Td;

// 		// Update particle length
// 		if (Lp[p]->Type==0)	CalPSizeR(p);
// 		// else 				CalPSizeCP(p);
// 		// else if (Lp[p]->Type==2 || Lp[p]->Type==3)	CalPSizeCP(p);
// 		// else if (Lp[p]->Type==2 || Lp[p]->Type==3)	CalPSizeR(p);
// 		// Update volume of particles
// 		Lp[p]->Vol 	= Lp[p]->Td.determinant()*Lp[p]->Vol0;
// 		// Update strain
// 		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
// 		// Update stress
// 		// Matrix3d w = 0.5*Dt*((Lp[p]->L - Lp[p]->L.transpose()));
// 		// Lp[p]->Stress += w*Lp[p]->Stress+Lp[p]->Stress*w.transpose();
// 		if (Lp[p]->Type==0)			Lp[p]->Elastic(de);
// 		else if (Lp[p]->Type==1)
// 		{
// 			Lp[p]->EOSMonaghan(Cs);
// 			Lp[p]->Newtonian(de);
// 		}
// 		else if (Lp[p]->Type==2)	Lp[p]->MohrCoulomb(de);
// 		else if (Lp[p]->Type==3)	Lp[p]->DruckerPrager(de);
// 		else if (Lp[p]->Type==5)
// 		{
// 			// Lp[p]->EOSMonaghan(Cs);
// 			Lp[p]->EOSMorris(Cs);
// 			Lp[p]->Granular(de);
// 		}
// 		// Reset hydro force and contact force
// 		Lp[p]->Fh.setZero();
// 		Lp[p]->Fc.setZero();
// 	}
// }

void MPM::NodeToParticle()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// cout << "p: " << p << endl;
		UpdateParticle(Lp[p]);
		// cout << "UpdateParticle" << endl;
		UpdateStress(Lp[p]);
		// cout << "UpdateStress" << endl;

	}
}

void MPM::NodeToParticleDoubleMapping()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateParticle(Lp[p]);
	}
	CalVOnNodeDoubleMapping();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateStress(Lp[p]);
	}		
}

void MPM::CalVOnNodeDoubleMapping()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<LAn.size(); ++c)
	{
		Ln[c]->V = Vector3d::Zero();
	}
	// Double map velocity from particles to nodes to aviod small mass problem
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 	n = Lp[p]->LnN[l];
			// weigthed mass contribution
			double nm = n*Lp[p]->M;
			// #pragma omp atomic
			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->V(d) += nm*Lp[p]->V(d)/Ln[id]->M;
			}
		}
    }
	// #pragma omp parallel for schedule(static) num_threads(1)
	// for (size_t c=0; c<LAn.size(); ++c)
	// {
	// 	Ln[c]->V = Ln[c]->Mv/Ln[c]->M;
	// }
}

// void MPM::CalVGradLocal(int p)
// {
// 	Lp[p]->L = Matrix3d::Zero();
// 	for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
// 	{
// 		size_t	 	id = Lp[p]->Lni[l];
// 		Vector3d 	gn 	= Lp[p]->LnGN[l];
// 		// Calculate velocity gradient tensor
// 		Lp[p]->L += Ln[id]->V*gn.transpose();
// 	}
// }

void MPM::CalVGradLocal(int p)
{
	Lp[p]->L = Matrix3d::Zero();

	size_t xmin = (int) Lp[p]->X(0);
	size_t ymin = (int) Lp[p]->X(1);
	size_t xmax = xmin+1;
	size_t ymax = ymin+1;

	int id0 = xmin+ymin*Ncy+0*Ncz;
	int id1 = xmax+ymin*Ncy+0*Ncz;
	int id2 = xmax+ymax*Ncy+0*Ncz;
	int id3 = xmin+ymax*Ncy+0*Ncz;

	double qx = Lp[p]->X(0)-xmin;
	double qy = Lp[p]->X(1)-ymin;

	Vector3d v0 = (1.-qx)*Ln[id0]->V + qx*Ln[id1]->V;
	Vector3d v1 = (1.-qy)*Ln[id1]->V + qy*Ln[id2]->V;
	Vector3d v2 = (1.-qx)*Ln[id3]->V + qx*Ln[id2]->V;
	Vector3d v3 = (1.-qy)*Ln[id0]->V + qy*Ln[id3]->V;

	Matrix3d Lfd = Matrix3d::Zero();

	Lfd(0,0) = 0.25*(v1(0)-v3(0));
	Lfd(1,0) = 0.25*(v1(1)-v3(1));

	Lfd(0,1) = 0.25*(v2(0)-v0(0));
	Lfd(1,1) = 0.25*(v2(1)-v0(1));


	for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
	{
	size_t id = Lp[p]->Lni[l];
	Vector3d gn = Lp[p]->LnGN[l];
	// Calculate velocity gradient tensor
	Lp[p]->L += Ln[id]->V*gn.transpose();
	}

	Lp[p]->L *= 0.75;
	Lp[p]->L += 0.25*Lfd;
}

// void MPM::CalVGradLocal(int p)
// {
// 	Lp[p]->L = Matrix3d::Zero();
// 	for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
// 	{
// 		size_t	 	id = Lp[p]->Lni[l];
// 		double 	n 	= Lp[p]->LnN[l];
// 		// Calculate velocity gradient tensor
// 		Lp[p]->L += n*Ln[id]->VGrad;
// 	}
// }

void MPM::CalPSizeCP(MPM_PARTICLE* p0)
{
	p0->PSize(0) =  p0->PSize0(0)*p0->Td(0,0);
	p0->PSize(1) =  p0->PSize0(1)*p0->Td(1,1);
	p0->PSize(2) =  p0->PSize0(2)*p0->Td(2,2);
}

// Based on "iGIMP: An implicit generalised interpolation material point method for large deformations"
void MPM::CalPSizeR(MPM_PARTICLE* p0)
{
	p0->PSize(0) =  p0->PSize0(0)*sqrt(p0->Td(0,0)*p0->Td(0,0) + p0->Td(1,0)*p0->Td(1,0) + p0->Td(2,0)*p0->Td(2,0));
	p0->PSize(1) =  p0->PSize0(1)*sqrt(p0->Td(0,1)*p0->Td(0,1) + p0->Td(1,1)*p0->Td(1,1) + p0->Td(2,1)*p0->Td(2,1));
	p0->PSize(2) =  p0->PSize0(2)*sqrt(p0->Td(0,2)*p0->Td(0,2) + p0->Td(1,2)*p0->Td(1,2) + p0->Td(2,2)*p0->Td(2,2));
}

// void MPM::CalStressOnParticleElastic()
// {
// 	// Update stresses on particles
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		// Velocity gradient tensor
// 		CalVGradLocal(p);
// 		// Update deformation tensor
// 		Lp[p]->Td = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->Td;
// 		// Update particle length
// 		CalPSizeR(p);
// 		// CalPSizeCP(p);
// 		// Update volume of particles
// 		Lp[p]->Vol 	= Lp[p]->Td.determinant()*Lp[p]->Vol0;
// 		// Update strain
// 		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
// 		// Update stress
// 		Matrix3d w = 0.5*Dt*((Lp[p]->L - Lp[p]->L.transpose()));
// 		Lp[p]->Stress += w*Lp[p]->Stress-Lp[p]->Stress*w.transpose();
// 		Lp[p]->Elastic(de);
// 	}
// }

// void MPM::CalStressOnParticleMohrCoulomb()
// {
// 	// Update stresses on particles
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		// Velocity gradient tensor
// 		CalVGradLocal(p);
// 		// Update deformation tensor
// 		Lp[p]->Td = (Matrix3d::Identity() + Lp[p]->L)*Lp[p]->Td;
// 		// Update particle length
// 		// CalPSizeR(p);
// 		CalPSizeCP(p);
// 		// Update volume of particles
// 		Lp[p]->Vol 	= Lp[p]->Td.determinant()*Lp[p]->Vol0;
// 		// Update strain
// 		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
// 		// Update stress
// 		Lp[p]->MohrCoulomb(de);
// 	}
// }

// void MPM::CalStressOnParticleNewtonian()
// {
// 	// cout << "start CalStressOnParticleNewtonian " << endl;
// 	// Update stresses on particles
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		// Velocity gradient tensor
// 		CalVGradLocal(p);
// 		// Update deformation tensor
// 		Lp[p]->Td = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->Td;
// 		// Update particle length
// 		// CalPSizeCP(p);
// 		// Update volume of particles
// 		Lp[p]->Vol 	= Lp[p]->Td.determinant()*Lp[p]->Vol0;
// 		// Update strain
// 		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
// 		// Update EOS
// 		// cout << "start EOSMorris " << endl;
// 		// Lp[p]->EOSMorris(Cs);
// 		Lp[p]->EOSMonaghan(Cs);
// 		// Update stress
// 		// cout << "start Newtonian " << endl;
// 		Lp[p]->Newtonian(de);
// 	// cout << "finish CalStressOnParticleNewtonian " << endl;

// 	}
// }

void MPM::SolveMUSL(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (ActiveShift)
		{
			for (size_t d=0; d<D; ++d)	Shift(d) = GetUniformD1();		
		}
		// cout << "Shift: " << Shift.transpose() << endl;
		bool show = false;
		if (t%100==0)	show = true;
		if (show) 	cout << "Time Step = " << t << endl;
		if (t%ts == 0)
		// if (t> 170600)
		// if (t> 40700)
		// if (t> 27000)
		{
			WriteFileH5(t);
		}
		auto t_start = std::chrono::system_clock::now();
		if (MLSv)	ParticleToNodeMLS();
		else 		ParticleToNode();
		auto t_end = std::chrono::system_clock::now();
		if (show)	cout << "ParticleToNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (UseFbar)	ParticleToCell();
		t_start = std::chrono::system_clock::now();
		CalVOnNode();
		// CalVGradOnNode();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalVOnNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		// if (t==2)
		// {
		// 	cout << Lp[4]->P << endl;
		// }
		// cout << "t= " << t << endl;
		// cout << Lp[4]->P << endl;
		NodeToParticle();
		// #pragma omp parallel for schedule(static) num_threads(Nproc)
  //   	for (size_t p=0; p<Lp.size(); ++p)
  //   	{
  //   		if (Lp[p]->V.norm()>0.5)
  //   		{
  //   			cout << "Lp[p]->V.norm()>0.5" << endl;
  //   			cout << "v:" << Lp[p]->V.transpose() << endl;
		// 		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		// 		{
		// 			size_t 	id = Lp[p]->Lni[l];
		// 			double 	n  = Lp[p]->LnN[l];
		// 			Vector3d an = Ln[id]->F/Ln[id]->M;
		// 			cout << "an: " << an.transpose() << endl;
		// 			cout << "Ln[id]->F: " << Ln[id]->F.transpose() << endl;
		// 			cout << "Ln[id]->Mv: " << Ln[id]->Mv.transpose() << endl;
		// 			cout << "Ln[id]->M: " << Ln[id]->M << endl;			
		// 		}
  //   			abort();
  //   		}
  //   	}
		// NodeToParticleDoubleMapping();

/*		size_t p=42396;

			cout << "t: " << t << endl;
			cout << "p: " << p << endl;
			cout << "vol0: " << Lp[p]->Vol0 << endl;
			cout << "vol: " << Lp[p]->Vol << endl;

			cout << "Lp[p]->V: " << Lp[p]->V.transpose() << endl;
			cout << "Lp[p]->X: " << Lp[p]->X.transpose() << endl;
			cout << "Lp[p]->PSize: " << Lp[p]->PSize.transpose() << endl;
			cout << "Lp[p+2]->PSize: " << Lp[p+2]->PSize.transpose() << endl;
			cout << "Td: " << endl;
			cout << Lp[p]->Td << endl;

			cout << "Td p+2: " << endl;
			cout << Lp[p+2]->Td << endl;

			// for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
			// {
			// 	size_t 	id = Lp[p]->Lni[l];
			// 	double 	n  = Lp[p]->LnN[l];
			// 	Vector3d an = Ln[id]->F/Ln[id]->M;
			// 	cout << "an: " << an.transpose() << endl;
			// 	cout << "Ln[id]->F: " << Ln[id]->F.transpose() << endl;
			// 	cout << "Ln[id]->Mv: " << Ln[id]->Mv.transpose() << endl;
			// 	cout << "Ln[id]->M: " << Ln[id]->M << endl;			
			// }
			cout << "==================" << endl;

		if (Lp[p]->V.norm()>0.015)	abort();*/
		// cout << "t= " << t << endl;
		// cout << Lp[4]->P << endl;

		t_end = std::chrono::system_clock::now();
		if (show)	cout << "NodeToParticle= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		// if 		(CMType==0) 	CalStressOnParticleElastic();
		// else if (CMType==1) 	CalStressOnParticleMohrCoulomb();
		// else if (CMType==2) 	CalStressOnParticleNewtonian();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalStressOnParticleMohrCoulomb= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (show) 	cout << "===========================" << endl;
	}
}

void MPM::AddNode(size_t level, Vector3d& x)
{
    Ln.push_back(new MPM_NODE(level,x));
    Ln[Ln.size()-1]->ID = Ln.size()-1;
}

void MPM::AddParticle(int tag, Vector3d& x, double m)
{
    Lp.push_back(new MPM_PARTICLE(tag,x,m));
    Lp[Lp.size()-1]->ID = Lp.size()-1;
}

void MPM::DeleteParticles()
{
	vector <MPM_PARTICLE*>	Lpt;
	Lpt.resize(0);

	for (size_t p=0; p<Lp.size(); ++p)
	{
		if (!Lp[p]->Removed)	Lpt.push_back(Lp[p]);
	}
	Lp = Lpt;

	for (size_t p=0; p<Lp.size(); ++p)
	{
		Lp[p]->ID = p;
	}
}

void MPM::AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m)
{
	Vector3i maxx = Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		maxx(d) = (int) (l(d)/ratio)-1;
	}

	for (int k=0; k<=maxx(2); ++k)
    for (int j=0; j<=maxx(1); ++j)
    for (int i=0; i<=maxx(0); ++i)
    {
    	Vector3d x = Vector3d::Zero();
    					x(0) = ratio*(i + 0.5)+x0(0);
    	if (D>1)		x(1) = ratio*(j + 0.5)+x0(1);
    	if (D>2)		x(2) = ratio*(k + 0.5)+x0(2);

    	AddParticle(tag, x, m);
    }

    for (size_t p=0; p<Lp.size(); ++p)
    {
    	Lp[p]->Vol0 = 1.;
    	for (size_t d=0; d<D; ++d)
    	{
    		Lp[p]->Vol0 *= ratio;
    		if (Ntype==3)	Lp[p]->PSize0(d) = 0.5*ratio;
    		else 			Lp[p]->PSize0(d) = 0.;
    		Lp[p]->PSize(d) = Lp[p]->PSize0(d);
    	}
    	Lp[p]->Vol 	= Lp[p]->Vol0;
    	// Assume a radius of MPM particle for DEMPM 
    	if (D==3)			Lp[p]->R = pow(3.*Lp[p]->Vol0/(4.*M_PI), 1./3.);
    	else if (D==2) 	Lp[p]->R = sqrt(Lp[p]->Vol0/M_PI);
    	Lp[p]->R = 0.5*ratio;
    }
}

inline void MPM::LoadMPMFromH5(string fname, double ratio)
{
	cout << "========= Start loading MPM particles from " << fname << "==============" << endl;
	H5std_string FILE_NAME( fname );
	H5std_string DATASET_NAME_POS( "Position" );
	H5File file_pos( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_pos = file_pos.openDataSet( DATASET_NAME_POS );
	DataSpace dataspace_pos = dataset_pos.getSpace();
    hsize_t dims_pos[2];
    dataspace_pos.getSimpleExtentDims( dims_pos, NULL);
    hsize_t dimsm_pos = dims_pos[0];
    cout <<"Position" << endl;

	H5std_string DATASET_NAME_VEL( "Velocity" );
	H5File file_vel( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_vel = file_pos.openDataSet( DATASET_NAME_VEL );
	DataSpace dataspace_vel = dataset_vel.getSpace();
    hsize_t dims_vel[2];
    dataspace_vel.getSimpleExtentDims( dims_vel, NULL);
    hsize_t dimsm_vel = dims_vel[0];
    cout <<"Velocity" << endl;

	H5std_string DATASET_NAME_PSIZE( "Psize" );
	H5File file_psize( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_psize = file_pos.openDataSet( DATASET_NAME_PSIZE );
	DataSpace dataspace_psize = dataset_psize.getSpace();
    hsize_t dims_psize[2];
    dataspace_psize.getSimpleExtentDims( dims_psize, NULL);
    hsize_t dimsm_psize = dims_psize[0];
    cout <<"Psize" << endl;

	H5std_string DATASET_NAME_S( "Stress" );
	H5File file_s( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_s = file_pos.openDataSet( DATASET_NAME_S );
	DataSpace dataspace_s = dataset_s.getSpace();
    hsize_t dims_s[2];
    dataspace_s.getSimpleExtentDims( dims_s, NULL);
    hsize_t dimsm_s = dims_s[0];
    cout <<"Stress" << endl;

	H5std_string DATASET_NAME_TD( "Derformation_Tensor" );
	H5File file_td( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_td = file_pos.openDataSet( DATASET_NAME_TD );
	DataSpace dataspace_td = dataset_td.getSpace();
    hsize_t dims_td[2];
    dataspace_td.getSimpleExtentDims( dims_td, NULL);
    hsize_t dimsm_td = dims_td[0];
    cout <<"Derformation_Tensor" << endl;

	H5std_string DATASET_NAME_M( "Mass" );
	H5File file_m( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_m = file_m.openDataSet( DATASET_NAME_M );
	DataSpace dataspace_m = dataset_m.getSpace();
    hsize_t dims_m[2];
    dataspace_m.getSimpleExtentDims( dims_m, NULL);
    hsize_t dimsm_m = dims_m[0];
    cout <<"Mass" << endl;

	H5std_string DATASET_NAME_VOL( "Volume" );
	H5File file_vol( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_vol = file_vol.openDataSet( DATASET_NAME_VOL );
	DataSpace dataspace_vol = dataset_vol.getSpace();
    hsize_t dims_vol[2];
    dataspace_vol.getSimpleExtentDims( dims_vol, NULL);
    hsize_t dimsm_vol = dims_vol[0];
    cout <<"Volume" << endl;

	H5std_string DATASET_NAME_TAG( "Tag" );
	H5File file_tag( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_tag = file_tag.openDataSet( DATASET_NAME_TAG );
	DataSpace dataspace_tag = dataset_tag.getSpace();
    hsize_t dims_tag[2];
    dataspace_tag.getSimpleExtentDims( dims_tag, NULL);
    hsize_t dimsm_tag = dims_tag[0];
    cout <<"Tag" << endl;

    double* data_pos = new double[dimsm_pos];
    dataset_pos.read( data_pos, PredType::NATIVE_DOUBLE, dataspace_pos, dataspace_pos );
    cout <<"data_pos" << endl;

    double* data_vel = new double[dimsm_vel];
    dataset_vel.read( data_vel, PredType::NATIVE_DOUBLE, dataspace_vel, dataspace_vel );
    cout <<"data_vel" << endl;

    double* data_s = new double[dimsm_s];
    dataset_s.read( data_s, PredType::NATIVE_DOUBLE, dataspace_s, dataspace_s );
    cout <<"data_s" << endl;

    double* data_td = new double[dimsm_td];
    dataset_td.read( data_td, PredType::NATIVE_DOUBLE, dataspace_td, dataspace_td );
    cout <<"data_td" << endl;

    double* data_m = new double[dimsm_m];
    dataset_m.read( data_m, PredType::NATIVE_DOUBLE, dataspace_m, dataspace_m );
    cout <<"data_m" << endl;

    double* data_tag = new double[dimsm_tag];
    dataset_tag.read( data_tag, PredType::NATIVE_DOUBLE, dataspace_tag, dataspace_tag );
    cout <<"data_tag" << endl;

    double* data_vol = new double[dimsm_vol];
    dataset_vol.read( data_vol, PredType::NATIVE_DOUBLE, dataspace_vol, dataspace_vol );
    cout <<"data_vol" << endl;

    double* data_psize = new double[dimsm_psize];
    dataset_psize.read( data_psize, PredType::NATIVE_DOUBLE, dataspace_psize, dataspace_psize );
    cout <<"data_psize" << endl;

    int np = dimsm_pos/3;
    for (int i=0; i<np; ++i)
    {
    	Vector3d pos (data_pos[3*i], data_pos[3*i+1], data_pos[3*i+2]);
    	Vector3d vel (data_vel[3*i], data_vel[3*i+1], data_vel[3*i+2]);
    	Vector3d psize (data_psize[3*i], data_psize[3*i+1], data_psize[3*i+2]);
    	Matrix3d stress;
		stress(0,0) = data_s[6*i];
		stress(0,1) = stress(1,0) = data_s[6*i+1];
		stress(0,2) = stress(2,0) = data_s[6*i+2];
		stress(1,1) = data_s[6*i+3];
		stress(2,1) = stress(1,2) = data_s[6*i+4];
		stress(2,2) = data_s[6*i+5];
    	Matrix3d td;
		td(0,0) = data_td[9*i];
		td(0,1) = data_td[9*i+1];
		td(0,2) = data_td[9*i+2];
		td(1,0) = data_td[9*i+3];
		td(1,1) = data_td[9*i+4];
		td(1,2) = data_td[9*i+5];
		td(2,0) = data_td[9*i+6];
		td(2,1) = data_td[9*i+7];
		td(2,2) = data_td[9*i+8];

    	double m = data_m[i];
    	double vol = data_vol[i];
    	int tag = (int) data_tag[i];
    	AddParticle(tag, pos, m);
    	Lp[Lp.size()-1]->V = vel;
    	Lp[Lp.size()-1]->PSize = psize;
    	Lp[Lp.size()-1]->Vol = vol;
    	Lp[Lp.size()-1]->Stress = stress;
    	Lp[Lp.size()-1]->Td = td;
    	Lp[Lp.size()-1]->Vol0 = 1.;
    	for (size_t d=0; d<D; ++d)
    	{
    		Lp[Lp.size()-1]->Vol0 *= ratio;
    		if (Ntype==3)	Lp[Lp.size()-1]->PSize0(d) = 0.5*ratio;
    	}
    }

    delete data_pos;
    delete data_vel;
    delete data_s;
    delete data_td;
    delete data_m;
    delete data_tag;
    delete data_vol;
    delete data_psize;

    cout << "========= Loaded "<< Lp.size()<< " MPM particles from " << fname << "==============" << endl;
}

inline void MPM::WriteFileH5(int n)
{
	stringstream	out;							         //convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "MPM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);     //create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};			   //create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.
	hsize_t	dims_tensor[1] = {6*Lp.size()};			//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);
	int rank_tensor = sizeof(dims_tensor) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);
	DataSpace	*space_tensor = new DataSpace(rank_tensor, dims_tensor);

	double* rho_h5 	= new double[  Lp.size()];
	double* vol_h5 	= new double[  Lp.size()];
	double* poro_h5 	= new double[  Lp.size()];
	double* gammap_h5 = new double[  Lp.size()];
	double* pos_h5 	= new double[3*Lp.size()];
	double* vel_h5 	= new double[3*Lp.size()];
	double* s_h5 	   = new double[6*Lp.size()];

	for (size_t i=0; i<Lp.size(); ++i)
	{
		rho_h5[  i  ]  = Lp[i]->M;
		vol_h5[  i  ] 	= Lp[i]->Vol;
		poro_h5[ i  ]  = Lp[i]->Poro;
		//gammap_h5[i ]  = Lp[i]->gammap_dot_tr;
		gammap_h5[i ]  = Lp[i]->gammap_dot;
 		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);

		s_h5  [6*i  ] 	= Lp[i]->Stress(0,0);
		s_h5  [6*i+1] 	= Lp[i]->Stress(0,1);
		s_h5  [6*i+2] 	= Lp[i]->Stress(0,2);
		s_h5  [6*i+3] 	= Lp[i]->Stress(1,1);
		s_h5  [6*i+4] 	= Lp[i]->Stress(1,2);
		s_h5  [6*i+5] 	= Lp[i]->Stress(2,2);
	}

	DataSet	*dataset_rho   = new DataSet(file.createDataSet("Rhos", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_vol	= new DataSet(file.createDataSet("Volume", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet  *dataset_poro  = new DataSet(file.createDataSet("Porosity", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet  *dataset_gammap= new DataSet(file.createDataSet("Gammap_dot", PredType::NATIVE_DOUBLE, *space_scalar));
   DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
   DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
   DataSet	*dataset_s   	= new DataSet(file.createDataSet("Stress", PredType::NATIVE_DOUBLE, *space_tensor));

	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_vol->write(vol_h5, PredType::NATIVE_DOUBLE);
	dataset_poro->write(poro_h5, PredType::NATIVE_DOUBLE);
	dataset_gammap->write(gammap_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
	dataset_s->write(s_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;
	delete dataset_rho;
	delete dataset_vol;
	delete dataset_poro;
	delete dataset_gammap;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_s;

	delete rho_h5;
	delete vol_h5;
	delete poro_h5;
	delete gammap_h5;
	delete pos_h5;
	delete vel_h5;
	delete s_h5;

	file.close();

	string file_name_xmf = "MPM_"+out.str()+".xmf";

	std::ofstream oss;
   oss.open(file_name_xmf);
   oss << "<?xml version=\"1.0\" ?>\n";
   oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
   oss << "<Xdmf Version=\"2.0\">\n";
   oss << " <Domain>\n";
   oss << "   <Grid Name=\"MPM\" GridType=\"Uniform\">\n";
   oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Lp.size() << "\"/>\n";
   oss << "     <Geometry GeometryType=\"XYZ\">\n";
   oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Lp.size() << " 3\" >\n";
   oss << "        " << file_name_h5 <<":/Position \n";
   oss << "       </DataItem>\n";
   oss << "     </Geometry>\n";
   oss << "     <Attribute Name=\"Rhos\" AttributeType=\"Scalar\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Rhos \n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "     <Attribute Name=\"Volume\" AttributeType=\"Scalar\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Volume \n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "     <Attribute Name=\"Porosity\" AttributeType=\"Scalar\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Porosity \n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "     <Attribute Name=\"Gammap_dot\" AttributeType=\"Scalar\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Gammap_dot \n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Velocity\n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
   oss << "       <DataItem Dimensions=\"" << Lp.size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
   oss << "        " << file_name_h5 <<":/Stress\n";
   oss << "       </DataItem>\n";
   oss << "     </Attribute>\n";
   oss << "   </Grid>\n";
   oss << " </Domain>\n";
   oss << "</Xdmf>\n";
   oss.close();
}

#endif
