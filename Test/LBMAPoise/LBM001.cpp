#include "LBM.h"

int main(int argc, char* argv[])
{
	bool   incompress               = false;
	double nuphysical               = 0.1;                      // (m^2)/s
	double rhophysical              = 1000;                     // kg/m^3
	Vector3d forcephysical (0.1, 0, 0);                         // m/(s^2)
	Vector3d initvphysical (0,0,0);                             // m/s

	int nx = 100;
	int ny = 50;
	int nz = 0;

	// Program works  
	LBM* x =  new LBM (D2Q9, SRT, incompress, nx, ny, nz, nuphysical);
	x->Init(rhophysical,initvphysical);

	for (int i = 0; i<101; ++i)
	for (int j = 1; j<50; ++j)
	{
		x->ABody[i][j][0]=forcephysical;
	}

	for (int i=0; i<101; ++i)
	{
		Vector3i h(i,0,0);
		Vector3i n(i,50,0);
		x->Lwall.push_back(h);
		x->Lwall.push_back(n);
	}

	for (int step = 1; step <=1000; ++step)
	{	
		// cout << "1" << endl;
		x->CalRhoV();
		// cout << "2" << endl;
		x->CollideSRT();
        // cout << "3" << endl;
		x->Stream();
		// cout << "4" << endl;
		x->ApplyWall();
		// cout << "5" << endl;
		cout<<"Time Step "<<step<<endl;
		cout<<x->V[50][0][0].transpose()<<endl;
		cout<<x->V[50][10][0].transpose()<<endl;
		cout<<x->V[50][25][0].transpose()<<endl;
		cout<<x->V[50][40][0].transpose()<<endl;
		cout<<x->V[50][50][0].transpose()<<endl;
		x->WriteFileH5(step,1);
	}

	return 0;
}