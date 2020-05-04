#include "LBM.h"

int main(int argc, char* argv[])
{
	bool   incompress               = false;
	double nuphysical               = 0.5;                      // (m^2)/s
	double rhophysical              = 1000;                     // kg/m^3
	Vector3d forcephysical (0, 0.5, 0);                         // m/(s^2)
	Vector3d initvphysical (0,0,0);                             // m/s
 
	int nx = 51;
	int ny = 101;
	int nz = 0;

	// Program works  
	LBM* x =  new LBM (D2Q9, SRT, incompress, nx, ny, nz, nuphysical);
	x->Init(rhophysical,initvphysical);

	// Set Driven Force
	for (int i = 1; i<=50; ++i)
	for (int j = 0; j<=101; ++j)
	{
		x->ABody[i][j][0]=forcephysical;
	}

	//Set Solid Fraction
	for (int i = 1; i<=50; ++i)
	for (int j = 0; j<=51; ++j)
	{
		x->Ns[i][j][0]= 0.8/1.8;
	}

	for (int i = 1; i<=50; ++i)
	for (int j = 52; j<=101; ++j)
	{
		x->Ns[i][j][0]= 0.2/1.2;
	}

	// Set Boundary Condition
	for (int j=0; j<=101; ++j)
	{
		Vector3i h(0,j,0);
		Vector3i n(51,j,0);
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
		cout<<x->V[25][25][0].transpose()<<endl;
		cout<<x->V[25][50][0].transpose()<<endl;
		cout<<x->V[25][75][0].transpose()<<endl;
		x->WriteFileH5(step,1);
	}

	return 0;
}