// 2D simulation of a elastic beam under gravity

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of one grid
	Vector3d gridSize (1,1,1);
	// Domain size
	int nx = 100;
	int ny = 100;
	int nz = 0;
	// Create MPM domain
	MPM* a = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// switch MLS for velocity interpolation
	a->MLSv = false;
	// Initialization
	a->Init();
	// Physcial parameters of particles
	double rhosPhysical 	= 2039.435;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double YoungPhysical 	= 7.5e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio
	// Space time and mass step
	double dx = 0.5;							// unit [m]
	double dt = 1.0e-4;							// unit [s]
	double dm = 1.0e-1;							// unit [kg]
	// How many particles in a cell
	double Ratio = 1./4.;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,2)/dm;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Poisson 	= PoissonPhysical;
	// Start point of the box for generating particles
	Vector3d x0 (20, 80, 0);
	// demention of the box
	Vector3d l0 (20, 4,  0);
	a->Nproc = 1;
	a->Dc = 0.;
	// Generate a box of particles
	a->AddBoxParticles(-1, x0, l0, Ratio, Mp);
	// Define gravity
	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		// a->Lp[p]->SetElastic(Young, Poisson);
		a->Lp[p]->SetGranular(Young, Poisson);
		a->Lp[p]->B  = G;
		a->Lp[p]->Stress(0,0) = -0.000001;
		a->Lp[p]->Stress(0,1) = 0.000001;
		a->Lp[p]->Stress(1,0) = -0.000001;
		a->Lp[p]->Stress(1,1) = 0.000001;
	}
	// Define boundary
	for (int i=19; i<=20; ++i)
	for (int j=0;  j<ny;  ++j)
	{
		a->SetNonSlippingBC(i,j,0);
	}

	// Solve
	for (int step=0; step <= 50000; ++step)
	{
		if (step % 100 == 0) // output the file per 100 times
		{
			cout<<"step == "<<step<<endl;
			cout<<"=============="<<endl;
			a->WriteFileH5(step);
		}

		a->ParticleToNode();

		a->CalVOnNode();

		a->NodeToParticle();
	}

	return 0;
}