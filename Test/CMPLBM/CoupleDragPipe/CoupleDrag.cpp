#include <MPLBM.h>

int main(int argc, char const *argv[])
{
	int nx = 200;
	int ny = 10;
	int nz = 0;

	// Solid phase parameters
	double rhosPhysical 	= 2500;				// Physical density, unit [kg/m^3]
	double niuPhysical      = 1.0e-3;           // Niu physica; [m^2/s]
	double YoungPhysical 	= 7.5e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio
	double SolidRatio 		= 0.58;				// solid volume fraction
	double DpPhysical       = 1.0e-3;           // Particle diamter, unit [m]
	// Space time and mass step
	double dx = 1.0e-4;							// unit [m]
	double dt = 1.0e-7;							// unit [s]
	double dm = 1;
	// How many particles in a cell
	double Ratio = 1./2.;
	// Vector3d G 		= (1.-rhofPhysical/rhosPhysical)*GPhysical*dt*dt/dx;
	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,2)/dm;
	double niu      = niuPhysical*dt/dx/dx;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Poisson 	= PoissonPhysical;
	double Rhos 	= rhosPhysical/dm*dx*dx*dx;
	double Rhop  	= Rhos/SolidRatio;
	double Dp       = DpPhysical/dx;

	Vector3d v0 (0., 0., 0.);			   // initial velocity	
	Vector3d G  (1.0e-5/2., 0., 0.);        // external force
	Vector3d gridSize (1,1,1);
	bool useFbar = false;

	MPLBM* a  = new MPLBM(nx, ny, nz);
	a->InitMPM(3, gridSize, useFbar);
	a->InitLBM(D2Q9, MRT, niu, 1.0, v0);

	a->DomLBM->Dp = Dp;
	a->DomLBM->Rhop = Rhop;

	a->DomLBM->Acc  = G;

	a->DomLBM->Nproc = 10;

	// a->DomLBM->SetFixedWall("Y_bot");
	// a->DomLBM->SetFixedWall("Y_top");

	a->DomLBM->Periodic[1] =  true;

	Vector3d x0 (50,  0, 0);
	// demention of the box
	Vector3d l0 (100, 11, 0);
	// Generate a box of particles
	a->DomMPM->AddBoxParticles(-1, x0, l0, Ratio, Mp);

	for (size_t p=0; p<a->DomMPM->Lp.size(); ++p)
	{
		a->DomMPM->Lp[p]->SetElastic(Young, Poisson);
		a->DomMPM->Lp[p]->B  = Vector3d::Zero();
	}

	cout << "Particle Diamter = " << Dp << endl;
	for (size_t t=0; t<=1.0e5; ++t)
	{
		a->WaterSoilCoupling();
		if (t%1000==0)
		{
			a->DomLBM->WriteFileH5("flowinporous",t);
			a->DomMPM->WriteFileH5(t);
			cout << "SolidFraction :"<<a->DomLBM->SolidFraction[100][5][0] << endl;
			cout << "Velocity :"<<a->DomLBM->V[100][5][0].transpose()<< endl;
			cout << "time: " << t << endl;
		}
		a->DomLBM->CalRhoV();
		a->DomLBM->CollideSRTPorous();	
		a->DomLBM->Stream();
	}

	return 0;
}