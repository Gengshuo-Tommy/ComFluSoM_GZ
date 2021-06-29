/***********************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

class MPM_PARTICLE
{
public:
	MPM_PARTICLE();
	MPM_PARTICLE(int tag, const Vector3d& x, double rho_input);
	void Elastic(Matrix3d& de);
	void SetElastic(double young, double poisson);
	void SetGranular(double young, double poisson, double rho_input, double d_input);
	void SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c);
	void SetTensionCutoff(double pmax);
	void CalcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out);
	void Granular(Matrix3d& de,  Matrix3d& dw);
	void DruckerPrager(Matrix3d& de);

    int 						Type;                       // Type of particle, 0 for elastic 1 for fluid 2 for soil.
	int 						ID; 				    	// Index of particle in the list 
	int 						Tag;				    	// Tag of particle

	double 						Rhos;				        // Particle density
	double 						Vol;						// Volume
	double 						Vol0;						// Init Volume
	double 						R;							// Radius for DEMPM
	double 						Arc;						// Arc length for boundary nodes
// Parameters for Drucker Prager constitutive laws ============================================
	double 						Mu;							// Shear modulus (Lame's second parameter) or viscosity for fluid
	double						La;							// Lame's first parameter
	double 						K;							// bulk modulus
	double 						H;							// Liniear hardening modulus
	double 						Young;						// Young's modus
	double						Poisson;					// Possion ratio
	double 						C;							// Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double 						Phi;						// Angle of internal friction
	double 						Psi;						// Angle of dilatation
	double 						A_dp;						// Drucker–Prager parameters
	double 						B_dp;						// Drucker–Prager parameters
	double 						Ad_dp;						// Drucker–Prager parameters
	double 						Pmax;						// Drucker–Prager parameters for tension cutoff (max pressure)
	bool 						TensionCut;					// Drucker–Prager parameters for tension cutoff

// Parameters for Inertial Rheology model ======================================================
	double 						grains_d;				    // Density of granular particles (not macro density of granular materials)
	double 						grains_rho;					// Diameter of granular particles
	double                      phi;                        // Solid packing fraction
	double                      phi_m;                      // Maximum packing fraction
	double                      phi_eq;                     // Equilibrium packing fraction
	double 						eta;                        // Water viscosity                   
	double 						I_tr;                       // Inertial number
	double 						I_v_tr;                     // Viscous inertial number
	double                      I_m_tr;                     // Mixed inertial number
	double                      I_0;                        
	double                      mu_1;                       
	double                      mu_2;                       
	double                      mu;                     
	double                      a;
	double                      K_3;      
	double                      gammap_dot_tr;              // Equivalent plastic shear strain      
// ======================================================================================
	double                      Poro;                       // Solid phase porosity
	double                      PoroInit;                   // Initial solid phase porosity
	double                      P;                          // Granular Pressure

	Vector3d					PSize0;						// Vector of half length of particle domain at init
	Vector3d					PSize;						// Vector of half length of particle domain

	Vector3d 					X;				            // Position
	Vector3d 					X0;				            // Init position
	Vector3d 					DeltaX;				        // Increasement of position
	Vector3d					V;							// Velocity
	Vector3d					Vf;							// Fixed velocity
	Vector3d					B;							// Body force acc
	Vector3d					Fh0;						// Hydro force
	Vector3d					Fh;							// Hydro force
	Vector3d					Fc;							// Contact force
	Vector3d					Nor;						// Normal direction (only non-zero for boundary particles)

	Matrix3d					Strain;						// Strain
	Matrix3d					StrainP;					// Plastic strain tensor
	Matrix3d					Stress;						// Stress
	Matrix3d					StressSmooth;				// Smoothed Stress for visualization
	Matrix3d					L;							// Velocity gradient tensor
	Matrix3d					Td;							// Derformation gradient tensor
	Matrix3d					Dp;							// Elastic tensor in principal stress space
	Matrix3d					Dpi;						// Inverse of Dp

	bool						FixV;						// Whether the velocity is fixed
	bool						Removed;					// whether this particle is removed

	vector<int>					Lnei;						// List of neighor nodes indexs, used to calculate arc lengh for FSI problems
	vector<size_t>				Lni;						// List of node indexs
	vector<size_t>				Lgi;						// List of gauss point indexs
	vector<double>				LnN;						// List of shape functions
	vector<Vector3d>			LnGN;						// List of gradient of shape functions
};

inline MPM_PARTICLE::MPM_PARTICLE()
{
    Type	= -1;
	ID		= 0;
	Tag		= 0;
	Rhos    = 0.;
	Vol     = 0.;
	X 		= Vector3d::Zero();
	X0 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fh0 	= Vector3d::Zero();
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

inline MPM_PARTICLE::MPM_PARTICLE(int tag, const Vector3d& x, double rho_input)
{
    Type	= -1;
	ID		= 0;
	Tag		= tag;
	Rhos    = rho_input;
	X 		= x;
	X0 		= x;
	Poro    = 0.0;
	PoroInit= 0.5;
	Vol     = 0.0;
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

inline void MPM_PARTICLE::SetElastic(double young, double poisson)
{
	Type 	= 0;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi = Dp.inverse();
}

inline void MPM_PARTICLE::SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c)
{
	Type 	= 3;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 	= Dp.inverse();
	Phi 	= phi;
	Psi 	= psi;
	C 		= c;
	TensionCut = false;
	// inner
	if (dptype==0)
	{
		double bot = sqrt(3)*(3.+sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.+sin(Psi)));
	}
	// outer
	else if (dptype==1)
	{
		double bot = sqrt(3)*(3.-sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.-sin(Psi)));
	}
	// plane strain
	else if (dptype==2)
	{
		double bot = sqrt(9.+12*tan(Phi)*tan(Phi));
		A_dp = 3.*tan(Phi)/bot;
		B_dp = 3./bot;
		Ad_dp = 3.*tan(Psi)/sqrt(9.+12*tan(Psi)*tan(Psi));
	}
}
// only works with Drucker Prager (12.3.2019)
inline void MPM_PARTICLE::SetTensionCutoff(double pmax)
{
	TensionCut 	= true;
	Pmax 		= pmax;
}

// Elastic model
inline void MPM_PARTICLE::Elastic(Matrix3d& de)
{
	Stress += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();
}

inline void MPM_PARTICLE::DruckerPrager(Matrix3d& de)
{
	// Apply elastic model first
	Elastic(de);
	double p = Stress.trace()/3.;
	Matrix3d ss = Stress - p*Matrix3d::Identity();
	double j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
	double f = j2sqr+A_dp*p-B_dp*C;
	if (f>0.)
	{
		// return to cone
		if (A_dp*(p-j2sqr/Mu*K*Ad_dp)-B_dp*C<0.)
		{
			double dl = (j2sqr+A_dp*p-B_dp*C)/(Mu+A_dp*K*Ad_dp);
			Stress -= dl*(Mu/j2sqr*ss+K*Ad_dp*Matrix3d::Identity());
		}
		// return to apex
		else
		{
			Stress = B_dp*C/A_dp*Matrix3d::Identity();
		}
		p = Stress.trace()/3.;
		ss = Stress - p*Matrix3d::Identity();
		j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
		double fb = f;
		f = j2sqr+A_dp*p-B_dp*C;
		if (f>1.0e-8)
		{
			cout << "f before: " << fb << endl;
			cout << "f after: " << f << endl;
			abort();
		}
	}
}

// Set granular parameters
inline void MPM_PARTICLE::SetGranular(double young, double poisson, double rho_input, double d_input)
{
	Type 	= 5;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;

	phi     = 0.6;
	phi_m   = 0.584;
	phi_eq  = 0.0;
	eta     = 0.0;
	I_tr    = 0.0;
	I_v_tr  = 0.0;
	I_m_tr  = 0.0;
	I_0     = 0.3085;
	mu_1    = 0.35;
	mu_2    = 1.387;
	mu      = 0.0;
	a       = 1.23;
	K_3     = 4.715;
	gammap_dot_tr = 0.0;

	grains_d   = d_input;
	grains_rho = rho_input;

	return;
}

// Calculate beta
inline void MPM_PARTICLE::CalcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out)
{
	I_out  = gdp * grains_d * sqrt(grains_rho/p);
    Iv_out = gdp * eta_in / p;
    Im_out = sqrt(I_out * I_out + 2 * Iv_out);

    if (Im_out == 0){
        mu_out = mu_1;
    } else if (p <= 0){
        mu_out = mu_2;
    } else {
        mu_out = mu_1 + (mu_2 - mu_1)/(1 + I_0/Im_out) + 2.5 * phi_in * Iv_out/(a*Im_out);
    }

    if ( p <= 0) {
        phi_eq = 0;
    } else {
        phi_eq = phi_m / (1 + a * Im_out);
    }

    beta_out = K_3 * (phi_in - phi_eq);

    return;
}


// Granular material plastic deformation
inline void MPM_PARTICLE::Granular(Matrix3d& de, Matrix3d& dw)
{
	Matrix3d T, D, W, T_tr;
	int k, j;
	double G, K;
	double trD, tau_bar_k, tau_bar_tr, tau_bar, p_tr, p_k, p, beta;
	double rp_max, rp_min, rt_max, rt_min, p_max, p_min, tau_max, tau_min, r_p, r_p_1, b_p, r_rt, b_rt;
	double p_zero_strength_f1, phi_eq_final, ABS_TOL, REL_TOL;
	bool is_solved, tau_too_large;

	ABS_TOL = 1.0e-30;
	REL_TOL = 1.0e-30;

	G   = Mu;
	K   = La;

	T   = Stress;

	D   = de;
	W   = dw;

	trD = D.trace();

	// calculate trial stress
	T_tr = T + 2 * Mu * D + La * trD * Matrix3d::Identity();
	tau_bar_tr = (T_tr - (T_tr.trace() / 3.0) * Matrix3d::Identity()).norm() / sqrt(2.0);
	p_tr = -T_tr.trace() / 3.0;

	// state determined
	is_solved     = false;
	tau_too_large = false;

	// initial value input
	tau_bar   = tau_bar_tr;
	p         = p_tr;

	// check elasticity
	if (!is_solved)
	{
		beta = K_3 *(phi - phi_m);

		if ( (tau_bar_tr <= ((mu_1 + beta)*p_tr)) && (p_tr >= 0) && (phi >= phi_m) )
		{
			gammap_dot_tr = 0;
			p = p_tr;
			tau_bar = tau_bar_tr; 
			is_solved = true;  // Check elasticity, f1,f2 and f3 should all be less than 0
		}
	}

	// check f1&&f2 condition
	if (!is_solved)
	{
		beta  = K_3 * (phi - 0.0);
		if ((p_tr + (K*tau_bar_tr/G)*beta) <= 0)
		{
			gammap_dot_tr = tau_bar_tr/ G;
			p = 0;
			tau_bar = 0;
			is_solved = true; // Under f2 solution, f1 $ f2 are both equal to zero, f3 <= 0.
		}
	}

	// Calculate maximaum P for function f1
	if (!is_solved){

	    // set iteration time and initial p_k value
		k   = 0;
		p_k = 0.0;
		// set bisection tolerance
		beta  = K_3 * (phi - phi_m);
		r_p_1 = max(abs(p_tr),abs(p_tr + (K*tau_bar_tr/G)*beta));
		b_p   = r_p_1; 

		// calculate bisection value at zero strength limit
		beta    = K_3 * (phi - 0.0); // beta max
		p_max   = p_tr + (K*tau_bar_tr/G)*beta; 
		p_min   = 0.0;

		// calculate gammap_dot_tr at tau equals to zero
		gammap_dot_tr = tau_bar_tr / G;

		// calculate state for p max
		CalcState(gammap_dot_tr, p_max, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

		// calculate residual for p_max
		rp_max = p_max - p_tr - K*beta*gammap_dot_tr; // must higher than 0

		// calculate state for p min
		CalcState(gammap_dot_tr, p_min, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

		// calculate residual for p_max
		rp_min = p_min - p_tr - K*beta*gammap_dot_tr; // must lower than 0 ..

		// cout << "1  " << rp_max << endl;
		// cout << "2  " << rp_min << endl;

		if (rp_min * rp_max > 0){
			cout << "ERROR: f1 weak residuals in binary search have same sign! That's bad!" << endl;
		}

		// bisection method
		while ( abs(r_p_1) > ABS_TOL && abs(r_p_1) / abs(b_p) > REL_TOL ){
			k += 1;
			if (k >50)
			{
				break;
			}

			//binary search
            p_k = 0.5*(p_max + p_min);

            //calculate state for step
            CalcState(gammap_dot_tr, p_k, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

            //calculate residual
            r_p_1 = p_k - p_tr - K*beta*gammap_dot_tr;

            //check residual sign
            if (r_p_1 * rp_min > 0){
                //r_p replaces r_min
                rp_min = r_p_1;
                p_min = p_k;
            } else {
                //r_p replaces r_max
                rp_max = r_p_1;
                p_max = p_k;
            }
		}

		p_zero_strength_f1 = p_k; // Maximum p satisfied
		// cout << "p_zero_strength_f1 ==" << p_zero_strength_f1 <<endl;
	}

	// check f1 condition
	if (!is_solved){
		// initial residual
		tau_bar_k = 0.0;
		r_rt   = tau_bar_tr;
		b_rt   = r_rt;

		// use bisection for function f1
		tau_max = tau_bar_tr;
		tau_min = 0.0;

		// initilize bisection
		gammap_dot_tr = 0.0;
		CalcState(gammap_dot_tr, p_tr,eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
		rt_max = tau_max - (mu + beta)*p_tr;

		gammap_dot_tr = tau_bar_tr/G;
		CalcState(gammap_dot_tr, p_zero_strength_f1, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
		rt_min = tau_min - (mu + beta)*p_zero_strength_f1;

		if (rt_max*rt_min>0)
		{
			cout << "ERROR: f1 residuals in binary search have same sign! That's bad!" << endl;
			abort();
		}

		tau_too_large = false;
		k = 0;

        while (abs(r_rt) > abs(b_rt)*REL_TOL && abs(r_rt) > ABS_TOL){
        	k += 1;
        	if (k > 50){
        		break;
        	}

            beta = phi - phi_m; //high pressure limit
            r_p  = max(abs(p_tr), abs( p_tr + (K * tau_bar_tr / G)*beta)); //reference
            b_p  = r_p;

            //set tau_bar_k
            tau_bar_k = 0.5 * (tau_max + tau_min);

            //calculate gammap_dot for tau going to zero
            gammap_dot_tr = (tau_bar_tr - tau_bar_k) / G;

            //calculate state for p_max
            CalcState(gammap_dot_tr, p_max, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

            //calculate residual for p_max
            rp_max = p_max - p_tr - K * beta * gammap_dot_tr;

            //calculate state for p_min
            CalcState(gammap_dot_tr, p_min, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

            //calculate residual for r_min
            rp_min = p_min - p_tr - K * beta * gammap_dot_tr;

            tau_too_large = false;
            if (rp_min * rp_max > 0) 
            {
                tau_max = tau_bar_k;
                tau_too_large = true;
                r_p = 0;            
            }

            // bisection method
            j = 0;
            while (std::abs(r_p) > ABS_TOL and std::abs(r_p) / std::abs(b_p) > REL_TOL) {
                j += 1;
                if (j > 50){
                	break;
                }

                //binary search
                p_k = 0.5 * (p_max + p_min);

                //calculate state for step
                CalcState(gammap_dot_tr, p_k, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                //calculate residual
                r_p = p_k - p_tr - K * beta * gammap_dot_tr;

                //check residual sign
                if (r_p * rp_min > 0) {
                    //r_p replaces r_min
                    rp_min = r_p;
                    p_min = p_k;
                } else{
                    //r_p replaces r_max
                    rp_max = r_p;
                    p_max = p_k;
                }
            }

            if (!tau_too_large) {
                //calculate equiv shear rate
                gammap_dot_tr = (tau_bar_tr - tau_bar_k) / G;

                //calculate state
                CalcState(gammap_dot_tr, p_k, eta, phi, I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                //calculate residual
                r_rt = tau_bar_k - (mu + beta) * p_k;

                //check sign of residual
                if (r_rt * rt_min > 0) {
                    rt_min = r_rt;
                    tau_min = tau_bar_k;
                } else {
                    rt_max = r_rt;
                    tau_max = tau_bar_k;
                }
            }
        }

        gammap_dot_tr = (tau_bar_tr - tau_bar_k) / G;

        phi_eq_final = phi_m / (1 + a*I_m_tr);

        if ((phi >= phi_m) || ( phi >= phi_eq_final) ) {
        	p = p_k;
            tau_bar = tau_bar_k;
            is_solved = true;
        }

	}

	// update stress
	if ( p >= 0 && tau_bar >0 )
	{
		T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace() / 3.0 * Matrix3d::Identity());
		T = T - p * Matrix3d::Identity();
	}
	else {
		T.setZero();
	}

	Stress = T;

	return;
}