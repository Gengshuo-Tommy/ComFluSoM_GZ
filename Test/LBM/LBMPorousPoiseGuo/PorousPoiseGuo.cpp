#include <LBM.h>

int main(int argc, char const *argv[])
{
	int L = 50;
	int nx = L;
	int ny = L;
	int nz = 0;

	double nu = 0.01;					// viscosity
	double rho0 = 1.;					// density
	Vector3d v0 (0., 0., 0.);			// initial velocity	

/*===============================================================================*/

	LBM* lbm = new LBM(D2Q9, MRT, nx, ny, nz, nu);
	lbm->Init(rho0, v0);
	lbm->InitSolidFraction();
	lbm->Nproc = 1;						// number of proccessors
	
	double Re = 0.01;
	double Da = 1e-4;
	double Poro = 0.1;
	
	double K = Da*L*L;
	double r = sqrt(Poro/K);
	double U0 = Re*nu/L;
	double G = U0*nu/K/(1.-1./cosh(0.5*L*r));
	double dp = sqrt(K*180.*pow(1.-Poro,2)/pow(Poro,3));

	cout << "k: " << K << endl;
	cout << "U0: " << U0 << endl;
	cout << "G" << G << endl;
	cout << "r: " << r << endl;

	lbm->Acc << 1.0e-6, 0., 0.;			// body force
	lbm->Dp = dp;

	lbm->SetFixedWall("Y_bot");
	lbm->SetFixedWall("Y_top");

	for (int i=0;  i<=nx; ++i)
	for (int j=0;  j<=ny; ++j)
	{
		lbm->SolidFraction[i][j][0] = 1.-Poro;
	}

	size_t tt = 1e4;
	size_t ts = 100;
	for (size_t t=0; t<tt; ++t)
	{
		if (t%ts==0)
		{
			lbm->WriteFileH5("PorousPoiseGuo",t);
			cout << lbm->V[25][25][0].transpose() << endl;
		}
		auto t_start = std::chrono::system_clock::now();
		lbm->SolveOneStepPorous();
		auto t_end = std::chrono::system_clock::now();
		double time0 = std::chrono::duration<double, std::milli>(t_end-t_start).count();

		if (t%ts==0)	cout << "time: " << time0 << endl;
	}

	ostringstream info;
	info << "re_" << Re << "_Da_" << Da << ".res";
	ofstream log(info.str().c_str(), ios_base::out | ios_base::app);
	log << "\"Y\"     \"V\"     \"Vt\"\n";
	size_t i = (size_t) 0.5*L;
	for (int j=0; j<=ny; ++j)
	{
		double y = (double) j;
		double vt = G*K/nu*(1.-cosh(r*(y-0.5*L))/cosh(r*0.5*L));
		log << setprecision(9) << fixed << y/L << "     " << lbm->V[i][j][0].norm()/U0 << "     " << vt/U0 << "\n";
	}
	return 0;
}