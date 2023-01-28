#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 50;
	int ny = 50;
	int nz = 0;

	double nu = 1.2e-2;

	Vector3d Veltop (1.0e-2, 0., 0.);
/*===============================================================================*/

	LBM* lbm = new LBM(D2Q9, SRT, nx, ny, nz, nu);

	Vector3d v0 (0.,0.,0.);
	lbm->Init(1., v0);

	for (int i=0; i<=nx; ++i)
	{
		VelBC nodeVel;
		nodeVel.Pos << i,50,0;
		nodeVel.Nei << i,49,0;
		nodeVel.Vel << Veltop;
		lbm->Lvel.push_back(nodeVel);
	}

	for (int ii=0; ii<=50; ++ii)
	for (int jj=0; jj<=50; ++jj)
	{
		lbm->G[ii][jj][0] = 0.9; 
	}

	lbm->SetFixedWall("Y_bot");

	for (int t=0; t<=200000; ++t)
	{
		if (t%100==0)
		{
			cout << "Time= " << t << endl;
			lbm->WriteFileH5("PorousCouetteG",t);
			cout << lbm->V[25][50][0].transpose() << endl;
			cout << lbm->V[25][40][0].transpose() << endl;
			cout << lbm->V[25][25][0].transpose() << endl;
			cout << lbm->V[25][10][0].transpose() << endl;
			cout << lbm->V[25][0][0].transpose() << endl;
		}
		lbm->CalRhoV();
		lbm->CollideSRT();
		lbm->ApplyVelBC();
		lbm->StreamGLBM();
		lbm->ApplyIBBFix();
	}
	return 0;
}