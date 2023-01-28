#ifndef MPLBM_H
#define MPLBM_H

#include <HEADER.h>
#include <MPM.h>
#include <LBM.h>
// #include <DRAG_FORCE_MODEL.h>

class MPLBM
{
public:
	MPLBM();
	~MPLBM();
	MPLBM(size_t nx, size_t ny, size_t nz);
	void InitMPM(size_t ntype, Vector3d dx, bool useFbar);
	void InitLBM(DnQm dnqm, CollisionModel cmodel, double nu, double rho0, Vector3d v0);
	void WaterSoilCoupling();
	void WaterSoilCouplingOneWay();
	void SolveOneStep();

	MPM*							DomMPM;													// Domain of MPM
	LBM*							DomLBM;													// Domain of LBM
    size_t 							Nx;														// Domain size
    size_t 							Ny;
    size_t 							Nz;	

	size_t 							Nproc;													// Number of processors which used
	int 							D;

    double 							Rhop;													// Density of solid particle
};

inline MPLBM::MPLBM(size_t nx, size_t ny, size_t nz)
{
	Nx = nx;
	Ny = ny;
	Nz = nz;	
}

inline void MPLBM::InitMPM(size_t ntype, Vector3d dx, bool useFbar)
{
	DomMPM = new MPM(ntype, Nx, Ny, Nz, dx);
	DomMPM->Init(useFbar);
	D = DomMPM->D;
}

inline void MPLBM::InitLBM(DnQm dnqm, CollisionModel cmodel, double nu, double rho0, Vector3d v0)
{
	DomLBM = new LBM(dnqm, cmodel, Nx, Ny, Nz, nu);
	DomLBM->Init(rho0, v0);
	DomLBM->InitSolidFraction();
	DomLBM->ResetSolidFraction = true;
	DomLBM->InitPressureGradient();
	DomLBM->InitCorrectMomentum();
}

inline void MPLBM::WaterSoilCoupling()
{
	for (size_t p=0; p<DomMPM->Lp.size(); ++p)
	{
		MPM_PARTICLE* p0 = DomMPM->Lp[p];

		Vector3d xs = DomMPM->Lp[p]->X;								// Position of solid
		double volp = p0->M/DomLBM->Rhop;							// solid particle volume
		double phi = volp/p0->Vol;									// solid volume fraction
		// double poro = 1.-phi;								    // porosity

		size_t kernel = 0;											// kernel type
		vector<StructMeshIndexWeight> lw;							// list of weight and index

		CalWeightStructMesh(DomLBM->D, DomLBM->DomSize, DomLBM->Periodic, kernel, xs, lw);

		Vector3d vf (0.,0.,0.);										// darcy velocity
		InterpolateStructMesh(lw, DomLBM->V, vf);

		double omegaf = 0.;
		InterpolateStructMesh(lw, DomLBM->Omega, omegaf);
		double nu = (1./omegaf-0.5)/3.;

		VectorXd ff (DomLBM->Q);
		ff.setZero();
		InterpolateStructMesh(lw, DomLBM->F, ff);

		double rhof; 
		Vector3d ve;
		DomLBM->CalRhoVLocalNoForce(ff, rhof, ve);

		Vector3d vs = p0->V;										// Solid velocity

		if (phi<0. || phi>1.)
		{
			cout << "wrong phi" << endl;
			cout << phi << endl;
			abort();
		}

		Vector3d vc (0., 0., 0.);
		Vector3d fe (0., 0., 0.);
		VelDragForceHoef(rhof, nu, phi, DomLBM->Dp, vf, vs, DomLBM->Acc, ve, vc, fe);

		fe *= p0->Vol;

		Vector3d pgrad (0., 0., 0.);
		InterpolateStructMesh(lw, DomLBM->PGrad, pgrad);

		// p0->V = pgrad;

		// p0->Fh = -p0->Vol*(fe - phi*pgrad);
		// p0->Fh = -fe;

		// pgrad << 0., -9.8e-8, 0.;

		p0->Fh = -fe-p0->Vol*phi*pgrad;

		// cout << "pgrad: " << pgrad.transpose() << endl;

		// // abort();

		// if (fe.norm()>1.e12)
		// {
		// 	cout << "rhof: " << rhof << endl;
		// 	cout << "nu: " << nu << endl;
		// 	cout << "phi: " << phi << endl;
		// 	cout << "dp: " << DomLBM->Dp << endl;
		// 	cout << "vf: " << vf.transpose() << endl;
		// 	cout << "vs: " << vs.transpose() << endl;
		// 	cout << "acc: " << DomLBM->Acc.transpose() << endl;
		// 	cout << "ve: " << ve.transpose() << endl;
		// 	cout << "vc: " << vc.transpose() << endl;
		// 	cout << "fe: " << fe.transpose() << endl;
		// 	abort();
		// }

		DistributeStructMesh(lw, volp, DomLBM->SolidFraction);
		DistributeStructMesh(lw, fe, DomLBM->ExForce);

		Vector3d cmv = rhof*p0->Vol*(vc-ve);
		DistributeStructMesh(lw, cmv, DomLBM->CMV);
	}

    #pragma omp parallel for schedule(static) num_threads(DomLBM->Nproc)
	for (size_t i = 0; i <= DomLBM->Nx; ++i)
    for (size_t j = 0; j <= DomLBM->Ny; ++j)
    for (size_t k = 0; k <= DomLBM->Nz; ++k)
    {
    	DomLBM->CalRhoVLocalNoForce(i,j,k);

    	double phi = DomLBM->SolidFraction[i][j][k];
		if (phi==0.)
		{
			DomLBM->V[i][j][k] += 0.5*DomLBM->Acc;
		}
		else
		{
			DomLBM->V[i][j][k] += DomLBM->CMV[i][j][k]/DomLBM->Rho[i][j][k];
		}
		DomLBM->CMV[i][j][k].setZero();
    }
}

inline void MPLBM::WaterSoilCouplingOneWay()
{
	for (size_t p=0; p<DomMPM->Lp.size(); ++p)
	{
		MPM_PARTICLE* p0 = DomMPM->Lp[p];

		Vector3d xs = DomMPM->Lp[p]->X;								// Position of solid
		double poro = p0->Poro;						 			    // solid porosity
		double phi  = 1. - poro;								    // solid volume fraction

		size_t kernel = 0;											// kernel type
		vector<StructMeshIndexWeight> lw;							// list of weight and index

		CalWeightStructMesh(DomLBM->D, DomLBM->DomSize, DomLBM->Periodic, kernel, xs, lw);

		// Vector3d vf (0.,0.,0.);										// darcy velocity
		// InterpolateStructMesh(lw, DomLBM->V, vf);

		// double omegaf = 0.;
		// InterpolateStructMesh(lw, DomLBM->Omega, omegaf);
		// double nu = (1./omegaf-0.5)/3.;

		VectorXd ff (DomLBM->Q);
		ff.setZero();
		InterpolateStructMesh(lw, DomLBM->F, ff);

		double rhof; 
		Vector3d ve;
		DomLBM->CalRhoVLocalNoForce(ff, rhof, ve);

		Vector3d vs = p0->V;										// Solid velocity

		if (phi<0. || phi>1.)
		{
			// cout << "wrong phi" << endl;
			// cout << phi << endl;
			phi = 1.00;
			// abort();
		}

		// Vector3d vc (0., 0., 0.);
		Vector3d fe (0., 0., 0.);
		// VelDragForceHoef(rhof, nu, phi, DomLBM->Dp, vf, vs, DomLBM->Acc, ve, vc, fe);
		fe  = 0.00001*(ve-vs);

		fe *= p0->Vol;

		Vector3d pgrad (0., 0., 0.);
		InterpolateStructMesh(lw, DomLBM->PGrad, pgrad);

		p0->Fh = -fe-p0->Vol*phi*pgrad;

	}
}

inline void MPLBM::SolveOneStep()
{
	WaterSoilCoupling();
	DomLBM->UpdateEffectiveViscosity();
	// DomLBM->CalRhoV();
	DomLBM->CalPGrad();
	// DomLBM->CollideSRTPorous();
	DomLBM->CollideMRTPorous();
	DomLBM->Stream();
	DomLBM->ApplyIBBFix();
	DomLBM->ApplyVelBC();
	DomLBM->ApplyNoGradBC();
	// DomLBM->UpdatePorousForce();
	// DomLBM->UpdatePorousForceHoef();

	DomMPM->ParticleToNode();
	if (DomMPM->UseFbar)	DomMPM->ParticleToCell();
	DomMPM->CalVOnNode();
	DomMPM->NodeToParticle();
}

#endif