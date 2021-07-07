// 2D Simulation Collapse of Rectangular Glass Beads

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of one grid
	Vector3d gridSize (1,1,1);
	// Domain size
	int nx = 200;
	int ny = 100;
	int nz = 0;
	// Glass Beads Parameters
	double rhosPhysical 	= 2500;			    // Soil density, unit [kg/m^3]
	Vector3d GPhysical (0, -9.8, 0);		    // Body force, unit [m/s^2]
	double YoungPhysical 	= 1.5e4;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Possion ratio
	double CPhysical		= 0.0;		        // Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double PhiPhysical		= 25./180.*M_PI;   // Angle of internal friction
	double PsiPhysical		= 0./180.*M_PI;		// Angle of dilatation
	// double Mu               = 0.25;
	// Space time and mass step
	double dx  = 1.0e-2;						// unit [m]
	double dt  = 1.0e-5;					    // unit [s]
	// How many particles in a cell
	double Ratio = 1./2.;
	// Dimentionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Rhos 	= rhosPhysical*dx*dx*dx;
	double Young 	= YoungPhysical*dx*dt*dt;
	double Poisson 	= PoissonPhysical;
	double C 		= CPhysical*dx*dt*dt;
	double Phi 		= PhiPhysical;      
	double Psi 		= PsiPhysical;
	double Kcoe     = Poisson/(1.-Poisson);

	// Start point of the box for generating particles
	Vector3d x0 (1,   1,  0);
	// demention of the box
	Vector3d l0 (39, 74, 0);

	// Create MPM domain
	MPM* x = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// switch MLS for velocity interpolation
	x->MLSv = false;
	// Initialization
	x->Init();
	x->AddBoxParticles(-1, x0, l0, Ratio, Rhos);

	x->Nproc = 10;
	x->Dc = 0.0;

	// Define Mohr Coulomb parameters and gravity
	for (size_t p=0; p<x->Lp.size(); ++p)
	{
		x->Lp[p]->SetDruckerPrager(2, Young, Poisson, Phi, Psi, C);
		x->Lp[p]->B  = G;
		// Init stress
		x->Lp[p]->Stress(1,1) = Rhos*G(1)*(75 - x->Lp[p]->X(1));
		x->Lp[p]->Stress(0,0) = Kcoe*x->Lp[p]->Stress(1,1);
		cout << x->Lp[p]->Stress << endl;
	}

	// Define boundary
	x->SetNonSlippingBC(0,0,0);
	x->SetNonSlippingBC(1,1,0);
    x->SetNonSlippingBC(0,1,0);
    x->SetNonSlippingBC(1,0,0);

	for (int i = 1;  i<=nx;  ++i)
	for (int j = 0;  j<=1; ++j)
	{
		Vector3d norm (0, -1, 0);
		x->SetSlippingBC(i, j, 0, norm);
	}

	for (int i =  0;  i<=1;   ++i)
	for (int j =  1;  j<=ny;  ++j)
	{
		// x->SetNonSlippingBC(i, j, 0);
		Vector3d norm (-1., 0, 0);
		x->SetSlippingBC(i, j, 0, norm);
	}

	for (int step = 0; step<= 200000; step++)
	{
		if (step % 1000 == 0 )
		{
			cout<<"Time Step == "<<step<<endl;
			x->WriteFileH5(step);
			cout<<x->Lp[5]->Stress<<endl;
			cout<<x->Lp[10]->Stress<<endl;
			cout<<x->Lp[160]->Stress<<endl;
			cout<<x->Lp[500]->Poro<<endl;
		}

		x->ParticleToNode();

		x->CalVOnNode();

		x->NodeToParticle();
	}

	return 0;
}