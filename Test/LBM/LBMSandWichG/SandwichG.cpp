#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 61;
	int ny = 61;
	int nz = 0;

	double nu = 1.2e-2;
	double k  = 0.25;
	double Ns1 = 0.3*0.5*nu/(k + nu);
	double Ns2 = 0.1*0.5*nu/(k + nu);

	Vector3d G (1.0e-6, 0., 0.);
/*===============================================================================*/

	LBM* lbm = new LBM(D2Q9, MRT, nx ,ny ,nz , nu);

	Vector3d v0 (0.,0.,0.);
	lbm->Init(1., v0);
	lbm->Acc = G;

	lbm->Nproc = 1;

	lbm->SetFixedWall("Y_bot");
	lbm->SetFixedWall("Y_top");

	for (int i=0; i<=nx; ++i)
	for (int j=1; j<=20; ++j)
	{
		lbm->G[i][j][0] = Ns1; 
	}

	for (int i=0;  i<=nx; ++i)
	for (int j=21; j<=40; ++j)
	{
		lbm->G[i][j][0] = Ns2; 
	}

	for (int i=0;  i<=nx; ++i)
	for (int j=41; j<=60; ++j)
	{
		lbm->G[i][j][0] = Ns1; 
	}

	for (int t=0; t<=500000; ++t)
	{
		if (t%2000==0)
		{
			cout << "Time= " << t << endl;
			lbm->WriteFileH5("SandwichG",t);
			cout << lbm->V[25][61][0].transpose() << endl;
			cout << lbm->V[25][40][0].transpose() << endl;
			cout << lbm->V[25][25][0].transpose() << endl;
			cout << lbm->V[25][10][0].transpose() << endl;
			cout << lbm->V[25][0][0].transpose() << endl;
		}
		lbm->CollideSRT();
		lbm->StreamGLBM();
		lbm->ApplyIBBFix();
		lbm->CalRhoV();	
	}
	return 0;
}