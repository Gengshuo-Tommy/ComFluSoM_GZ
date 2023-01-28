#include <LBM.h>

int main(int argc, char const *argv[])
{
	size_t nx = 50;
	size_t ny = 50;
	size_t nz = 0;

	double rhofPhysical 	= 1000.;						// Physical density, unit [kg/m^3]
	double miuPhysical 		= 1.0e-3;						    // Physical viscosity, unit [N.s/m^2]
	double nuPhysical 		= miuPhysical/rhofPhysical;		// Physical viscosity, unit [m^2/s]
	Vector3d GPhysical (0., -9.8, 0.);						// Body force

	// Space time and mass step
	double dx = 0.01;							// unit [m]
	double dt = 1.0e-5;							// unit [s]
	double dm = 1.0e-3;							// unit [kg]

	double Rhof 	= rhofPhysical/dm*dx*dx*dx;
	double Nu 		= nuPhysical*dt/(dx*dx);
	// double Nu 		= 0.1;
	Vector3d G 		= GPhysical*dt*dt/dx;

	double K  = 1.0e-4*ny*ny;
	double Poro  = 0.5;
	double dp = sqrt(K*180.*pow(1.-Poro,2)/pow(Poro,3));

	Vector3d v0 (0., 0., 0.);			// initial velocity	

	cout << "nu: " << Nu << endl;
	cout << "Rhof: " << Rhof << endl;
	cout << "G: " << G.transpose() << endl;

	// abort();
/*===============================================================================*/
	LBM* lbm = new LBM(D2Q9, MRT, nx, ny, nz, Nu);
	// LBM* lbm = new LBM(D3Q15, MRT, nx, ny, nz, nu);
	lbm->Init(Rhof, v0);
	lbm->InitSolidFraction();
	lbm->InitPressureGradient();
	lbm->Dp =  dp;
	lbm->Acc = G;

	lbm->Nproc = 20;

	lbm->SetFixedWall("Y_bot");
	lbm->SetFixedWall("Y_top");

	for (int i=0;  i<=50; ++i)
	for (int j=0;  j<=25; ++j)
	{
		lbm->SolidFraction[i][j][0] = 0.9;
	}

	size_t tt = 10e4;
	size_t ts = 1000;
	for (size_t t=0; t<tt; ++t)
	{
		if (t%ts==0)
		{
			lbm->WriteFileH5("pressure",t);
		}
		auto t_start = std::chrono::system_clock::now();
		lbm->SolveOneStepPorous();
		// lbm->SolveOneStep();
		auto t_end = std::chrono::system_clock::now();
		double time0 = std::chrono::duration<double, std::milli>(t_end-t_start).count();

		if (t%ts==0)	cout << "time: " << time0 << endl;
	}

	lbm->CalPGrad();

	// lbm->CalPGradLocal(25, 25, 0);
	for (int j = 0; j <= 50; j++)
	{
		cout << "pressure gradient: " << lbm->PGrad[25][j][0].transpose() << "  "<< j << endl;
	}

	// lbm->WriteFileH5("pressure",0);
	// lbm->WriteFileH5("pressure",1);

	return 0;
}