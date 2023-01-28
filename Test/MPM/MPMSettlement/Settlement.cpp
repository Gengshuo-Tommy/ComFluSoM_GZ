// 2D simulation of settlement under gravity

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
	bool useFbar = false;
	x->Init(useFbar);
	// Physcial parameters of particles
	double RhosPhysical 	= 2650;			    // Physical density, unit [kg/m^3]
	// double RhosDiaPhysical  = 225e-6;           // Diameter [m]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double YoungPhysical 	= 9.8e4;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio
	// double Solidfraction    = 0.585;           

	// Space time and mass step
	double dx = 1.0e-3;							// unit [m]
	double dt = 1.0e-5;					        // unit [s]
	// How many particles in a cell
	double Ratio = 1.0/4.0;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Rhos     = RhosPhysical*pow(dx,3);
	// double Rhos_Dia = RhosDiaPhysical/dx;
	double Young 	= YoungPhysical*dx*dt*dt;
	double Poisson 	= PoissonPhysical;
	// double Kcoe     = Poisson/(1 - Poisson);
	// double Poro     = 1 - Solidfraction;
	// Start point of the box for generating particles
	Vector3d x0 (20, 20, 0);
	// demention of the box
	Vector3d l0 (50, 50,  0);
	x->Nproc = 10;
	x->Dc = 0.2;
	// Generate a box of particles
	x->AddBoxParticles(-1, x0, l0, Ratio, Rhos);
	// Define gravity
	for (size_t p=0; p<x->Lp.size(); ++p)
	{
		x->Lp[p]->SetElastic(Young, Poisson);
		// x->Lp[p]->SetGranular(Young, Poisson, Rhos, Rhos_Dia, Poro);
		x->Lp[p]->B  = G;
	}

	// Define boundary
	for (int i=0; i<=nx; ++i)
	{
		x->SetNonSlippingBC(i,20,0);
	}

	for (int j = 0; j <= ny; ++j)
	{
		Vector3d norm1 (-1., 0. , 0.);
		// x->DomMPM->SetSlippingBC(1,  j,  0, norm1);
		x->SetSlippingBC(20, j, 0, norm1);
	}
	for (int j = 0; j <= ny; ++j)
	{
		Vector3d norm2 (1., 0. , 0.);
		// x->DomMPM->SetSlippingBC(1,  j,  0, norm2);
		x->SetSlippingBC(70,  j,  0, norm2);
	}
	cout << x->Lp.size() << endl;
	// Solve
	for (int step=0; step <= 1.0e6; ++step)
	{
		if (step % 1000 == 0) // output the file per 100 times
		{
			cout<<"step == "<<step<<endl;
			cout<<"=============="<<endl;
			x->WriteFileH5(step);
		}

		// x->ParticleToNodePoro();

		x->ParticleToNode();

		x->CalVOnNode();

		x->NodeToParticle();

	}

	return 0;
}