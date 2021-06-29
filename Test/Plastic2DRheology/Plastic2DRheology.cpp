// 2D simulation of sand

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of one grid
	Vector3d gridSize (1,1,1);
	// Domain size
	int nx = 150;
	int ny = 100;
	int nz = 0;
	// Create MPM domain
	MPM* x = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// switch MLS for velocity interpolation
	x->MLSv = false;
	// Initialization
	x->Init();
	// Physcial parameters of particles
	double RhosPhysical 	= 2500;			    // Physical density, unit [kg/m^3]
	double RhosDiaPhysical  = 225e-6;           // Diameter [m]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double YoungPhysical 	= 7.5e4;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio

	// Space time and mass step
	double dx = 1.0e-2;							// unit [m]
	double dt = 1.0e-5;					        // unit [s]
	// How many particles in a cell
	double Ratio = 1.0/4.0;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Mp 		= RhosPhysical*pow(dx,3)*pow(Ratio,2);
	double Rhos     = RhosPhysical*pow(dx,3);
	double Rhos_Dia = RhosDiaPhysical/dx;
	double Young 	= YoungPhysical*dx*dt*dt;
	double Poisson 	= PoissonPhysical;
	double Kcoe     = Poisson/(1 - Poisson);
	// Start point of the box for generating particles
	Vector3d x0 (1,   1,  0);
	// demention of the box
	Vector3d l0 (40, 80,  0);
	x->Nproc = 10;
	x->Dc = 0.05;
	// Generate a box of particles
	x->AddBoxParticles(-1, x0, l0, Ratio, Mp);
	// Define gravity and parameters
	for (size_t p=0; p<x->Lp.size(); ++p)
	{
		x->Lp[p]->SetGranular(Young, Poisson, Rhos, Rhos_Dia);
		x->Lp[p]->B  = G;
	}
	// // set initial stress
	for (size_t p=0; p<x->Lp.size(); ++p)
	{
		x->Lp[p]->Stress(1,1) = Rhos*G(1)*(41 - x->Lp[p]->X(1));
		x->Lp[p]->Stress(0,0) = Kcoe*x->Lp[p]->Stress(1,1);
		cout << x->Lp[p]->Stress << endl;
	}

	x->SetNonSlippingBC(0,0,0);
	x->SetNonSlippingBC(1,1,0);
    x->SetNonSlippingBC(0,1,0);
    x->SetNonSlippingBC(1,0,0);

	// Define boundary
	for (int i = 0;  i <=  1;  ++i)
	for (int j = 1;  j <= ny;  ++j)
	{
		x->SetNonSlippingBC(i,j,0);
	}

	for (int i = 1;  i <= nx;  ++i)
	for (int j = 0;  j <= 1;   ++j)
	{
		// a->SetNonSlippingBC(i,j,0);
		Vector3d norm (0, -1, 0);
		x->SetSlippingBC(i, j, 0, norm);
	}
	// Solve
	for (int step = 0; step <= 50000; ++step)
	{
		if (step % 1000 == 0) // output the file per 100 times
		{
			cout<<"step == "<<step<<endl;
			cout<<"=============="<<endl;
			x->WriteFileH5(step);
		}

		x->ParticleToNodePoro();

		x->CalVOnNode();

		x->NodeToParticle();
	}

	return 0;
}