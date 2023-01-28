#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 10;
	int ny = 100;  
	int nz = 0;

	double niu = 0.5;   // unit [m^2/s]

/*===============================================================================*/

	LBM* x = new LBM(D2Q9, SRT, nx, ny, nz, niu);

	Vector3d v0 (0.,0.,0.);
	x->Init(1.001, v0);

	//Set Solid Fraction
	for (int i=0; i<=nx; ++i)
	for (int j=0; j<=ny; ++j)
	{
		x->G[i][j][0] = 0.49; 
	}

	// Set Non Gradient Boundary Condition
	for (int i=0; i<=nx; ++i)
	{
		NoGradBC nodeNog;
		nodeNog.Pos << i,0,0;
		nodeNog.Nei << i,1,0;
		x->Lnog.push_back(nodeNog);
	}
	
	for (int t=0; t<=10000; ++t)
	{
		for (int i=0;i<=nx;++i)
		{
			x->CalFeqC(x->F[i][100][0], 1., x->V[i][100][0]);
			VectorXd feqn(x->Q), fneq(x->Q);
			x->CalFeqC(feqn, x->Rho[i][99][0], x->V[i][99][0]);
			fneq = x->F[i][99][0]-feqn;
			x->F[i][100][0] += fneq;
		}
		cout << "1" << endl;
		x->CalRhoV();
		if (t%10==0)
		{
			cout << "Time= " << t << endl;
			x->WriteFileH5("ConsolidationWaterPhase",t);
			cout << "top dentsiy "<< x->Rho[5][100][0] << endl;
		}
		cout << "2" << endl;
		x->CollideSRT();
		cout << "3" << endl;
		x->StreamGLBM();
		cout << "4" << endl;
		x->ApplyNoGradBC();
	}
	return 0;
}