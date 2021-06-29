// 2D simulation of Darcy flow

#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 50;
	int ny = 120;
	int nz = 0;

	double nu = 0.1;

	Vector3d ac (0., -1.0e-6, 0.);

/*===============================================================================*/

	LBM* lbm = new LBM(D2Q9, SRT, false, nx ,ny ,nz , nu);

	Vector3d v0 (0.,0.,0.);
	lbm->Init(1., v0);

	lbm->SetA(ac);

	for (int i = 0; i <= nx;  ++i)
	for (int j = 0; j <= 39;  ++j)
	{
		lbm->Ns[i][j][0] = 0.5/1.5;
	}

	for (int i = 0;  i <= nx;  ++i)
	for (int j = 40; j <= 80;  ++j)
	{
		lbm->Ns[i][j][0] = 0.1/1.1;
	}

	for (int i = 0;  i <= nx;  ++i)
	for (int j = 81; j <= 120;  ++j)
	{
		lbm->Ns[i][j][0] = 0.5/1.5;
	}

	// for (int i = 0; i <= nx; ++i)
	// {
	// 	Vector3i x_top (i,ny,0);
	// 	Vector3i x_bot (i,0 ,0);
	// 	lbm->Lwall.push_back(x_top);
	// 	lbm->Lwall.push_back(x_bot);
	// }

	for (int j = 0; j <= ny; ++j)
	{
		Vector3i x_left  (0,  ny, 0);
		Vector3i x_right (50, ny ,0);
		lbm->Lwall.push_back(x_left);
		lbm->Lwall.push_back(x_right);
	}

	for (int t=0; t<10000; ++t)
	{
		if (t%100==0)
		{
			cout << "Time= " << t << endl;
			lbm->WriteFileH5(t,1);
		}
		lbm->CollideSRT();
		lbm->ApplyWall();
		lbm->StreamGLBM();
		lbm->CalRhoV();	
	}

	return 0;
}