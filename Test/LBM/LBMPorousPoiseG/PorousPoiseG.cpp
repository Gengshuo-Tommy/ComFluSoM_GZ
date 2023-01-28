#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 50;
	int ny = 50;
	int nz = 0;

	double Rhof 	= 1.;
	double Nu 		= 0.1;
	Vector3d G (1.0e-6, 0., 0.);

	Vector3d v0 (0., 0., 0.);			// initial velocity	

	// v0 -= 0.5*G;

	// cout << "nu: " << Nu << endl;
	// cout << "Rhof: " << Rhof << endl;
	// cout << "G: " << G.transpose() << endl;

	// abort();
/*===============================================================================*/
	LBM* lbm = new LBM(D2Q9, SRT, nx, ny, nz, Nu);
	
	lbm->Acc = G;
	lbm->Init(Rhof, v0);

	lbm->Nproc = 1;

	lbm->SetFixedWall("Y_bot");
	lbm->SetFixedWall("Y_top");

	for (int i=0; i<=50; ++i)
	for (int j=0; j<=50; ++j)
	{
		lbm->G[i][j][0] = 0.1; 
	}

	size_t tt = 1e4;
	size_t ts = 100;
	for (size_t t=0; t<tt; ++t)
	{
		if (t%ts==0)
		{
			cout << "time = " << t << endl;
			lbm->WriteFileH5("PorousPoiseGray",t);
			cout << lbm->V[25][50][0].transpose() << endl;
			cout << lbm->V[25][40][0].transpose() << endl;
			cout << lbm->V[25][25][0].transpose() << endl;
			cout << lbm->V[25][10][0].transpose() << endl;
			cout << lbm->V[25][0][0].transpose() << endl;
		}
		lbm->CalRhoV();
		lbm->CollideSRT();
	    lbm->StreamGLBM();
	    // lbm->Stream();
	    lbm->ApplyIBBFix();
	}

	return 0;
}