// 2D simulation of an elastic beam under gravity

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
	MPM* x = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// switch MLS for velocity interpolation
	x->MLSv = false;
	// Initialization
	x->Init();
	// Physcial parameters of particles
	double rhosPhysical 	= 2600;			    // Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double YoungPhysical 	= 5.0e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio
	// Space time and mass step
	double dx = 1.0e-1;							// unit [m]
	double dt = 1.0e-4;							// unit [s]
	// How many particles in a cell
	double Ratio = 1./4.;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Rhos     = rhosPhysical*pow(dx,3);
	double Young 	= YoungPhysical*dx*dt*dt;
	double Poisson 	= PoissonPhysical;
	// Start point of the box for generating particles
	Vector3d x0 (20, 80, 0);
	// demention of the box
	Vector3d l0 (20, 4,  0);
	x->Nproc = 1;
	x->Dc = 0.;
	// Generate a box of particles
	x->AddBoxParticles(-1, x0, l0, Ratio, Rhos);
	// Define gravity
	for (size_t p=0; p<x->Lp.size(); ++p)
	{
		x->Lp[p]->SetElastic(Young, Poisson);
		// a->Lp[p]->SetGranular(Young, Poisson, 10, 10);
		x->Lp[p]->B  = G;
	}
	// Define boundary
	for (int i=19; i<=20; ++i)
	for (int j=0;  j<ny;  ++j)
	{
		x->SetNonSlippingBC(i,j,0);
	}

	// Solve
	for (int step=0; step <= 50000; ++step)
	{
		if (step % 100 == 0) // output the file per 100 times
		{
			cout<<"step == "<<step<<endl;
			cout<<"=============="<<endl;
			x->WriteFileH5(step);
		}

		// a->ParticleToNodePoro();
		x->ParticleToNode();

		x->CalVOnNode();

		x->NodeToParticle();
	}

	return 0;
}