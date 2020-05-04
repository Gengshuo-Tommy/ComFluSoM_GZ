#include "LBM.h"

int main(int argc, char* argv[])
{
	bool   incompress               = false;
	double nuphysical               = 0.5;                      // (m^2)/s
	double rhophysical              = 1000;                     // kg/m^3
	Vector3d forcephysical (0.1, 0, 0);                         // m/(s^2)
	Vector3d initvphysical (0,0,0);                             // m/s
 
	int nx = 20;
	int ny = 103;
	int nz = 0;

	// Program works  
	LBM* x =  new LBM (D2Q9, SRT, incompress, nx, ny, nz, nuphysical);
	x->Init(rhophysical,initvphysical);

	// Set Driven Force
	for (int i = 0; i<=20; ++i)
	for (int j = 1; j<=102; ++j)
	{
		x->ABody[i][j][0]=forcephysical;
	}

	//Set Solid Fraction
	for (int i = 0; i<=20; ++i)
	for (int j = 1; j<=34; ++j)
	{
		x->Ns[i][j][0]= 0.1;
	}

	for (int i = 0; i<=20; ++i)
	for (int j = 35; j<=68; ++j)
	{
		x->Ns[i][j][0]= 0;
	}

	for (int i = 0; i<=20; ++i)
	for (int j = 69; j<=102; ++j)
	{
		x->Ns[i][j][0]= 0.1;
	}

	// Set Boundary Condition
	for (int i=0; i<=20; ++i)
	{
		Vector3i h(i,0,0);
		Vector3i n(i,103,0);
		x->Lwall.push_back(h);
		x->Lwall.push_back(n);
	}

	for (int step = 1; step <=75; ++step)
	{			

		x->CalRhoVGray();                            // calculate velocity
		
		x->CollideSRTWalsh();

		x->StreamWalsh();

		x->ApplyWall();

		cout<<"Time Step "<<step<<endl;
		cout<<x->V[10][0][0].transpose()<<endl;
		cout<<x->V[10][20][0].transpose()<<endl;
		cout<<x->V[10][51][0].transpose()<<endl;
		cout<<x->V[10][82][0].transpose()<<endl;
		cout<<x->V[10][103][0].transpose()<<endl;
		x->WriteFileH5(step,1);
	}

	return 0;
}