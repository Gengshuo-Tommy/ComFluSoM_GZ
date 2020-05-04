#include "LBM.h"

int main(int argc, char* argv[])
{
	bool   incompress               = false;
	double nuphysical               = 0.5;                      // (m^2)/s
	double rhophysical              = 1000;                     // kg/m^3
	Vector3d forcephysical (0.1, 0, 0);                         // m/(s^2)
	Vector3d initvphysical (0,0,0);                             // m/s
 
	int nx = 100;
	int ny = 50;
	int nz = 0;

	// Program works  
	LBM* x =  new LBM (D2Q9, SRT, incompress, nx, ny, nz, nuphysical);
	x->Init(rhophysical,initvphysical);

	// Set Driven Force
	for (int i = 0; i<=100; ++i)
	for (int j = 1; j<=49; ++j)
	{
		x->ABody[i][j][0]=forcephysical;
	}

	//Set Solid Fraction
	for (int i = 0; i<=100; ++i)
	for (int j = 1; j<=49; ++j)
	{
		x->Ns[i][j][0]= 0.5/1.5;
	}

	// Set Boundary Condition
	for (int i=0; i<=100; ++i)
	{
		Vector3i h(i,0,0);
		Vector3i n(i,50,0);
		x->Lwall.push_back(h);
		x->Lwall.push_back(n);
	}

	for (int step = 1; step <=1000; ++step)
	{			

		x->CalRhoVGray();                            // calculate velocity
		
		x->CollideSRTZhu();

		x->StreamGray();

		x->ApplyWall();

		cout<<"Time Step "<<step<<endl;
		cout<<x->V[50][15][0].transpose()<<endl;
		cout<<x->V[50][25][0].transpose()<<endl;
		cout<<x->V[50][50][0].transpose()<<endl;
		x->WriteFileH5(step,1);
	}

	return 0;
}