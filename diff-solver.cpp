#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "phase.h"

using namespace std;

Real alpha = 0.75;

Real K;
Real D1 = 0.01;
Real D2 = 0.01;
Real w1 = 1;
Real w2 = 5;
Real dt = 0.001;
Real r1,r2, r1_new, r2_new;
Real psi1,psi2, psi1_new, psi2_new;

ofstream output ;

void Transform()
{
	int n1 = floor(psi1 / (2*M_PI));
	int n2 = floor(psi2 / (2*M_PI));
	psi1 -= n1*2*M_PI;
	psi2 -= n2*2*M_PI;
}

Real Find_Integrand(const Real& w, const Real& fraction, const Real& r, const Real& psi, const Real& psi_other)
{
	return (w + 0.5*(1-alpha)*K*(r*r*r + 1/r)*r2*sin(psi2 - psi1));
}

void One_Step(Real dt)
{
	Real psi1_integrand = (w1 + 0.5*alpha*K*(r1*r1*r1 + 1/r1)*r2*sin(psi2 - psi1));
	Real psi2_integrand = (w2 + 0.5*(1-alpha)*K*(r2*r2*r2 + 1/r2)*r1*sin(psi1 - psi2));
	Real c = cos(psi2-psi1);
	Real r1_integrand = 0.5*K*(1-r1*r1*r1*r1)*((1-alpha)*r1 + alpha*r2*c) - D1*r1;
	Real r2_integrand = 0.5*K*(1-r2*r2*r2*r2)*(alpha*r2 + (1-alpha)*r1*c) - D2*r2;

	psi1_new = psi1 + psi1_integrand*dt;
	psi2_new = psi2 + psi2_integrand*dt;
	r1_new = r1 + r1_integrand*dt;
	r2_new = r2 + r2_integrand*dt;

	Real psi1_new_integrand = (w1 + 0.5*alpha*K*(r1_new*r1_new*r1_new + 1/r1_new)*r2_new*sin(psi2_new - psi1_new));
	Real psi2_new_integrand = (w2 + 0.5*(1-alpha)*K*(r2_new*r2_new*r2_new + 1/r2_new)*r1*sin(psi1_new - psi2_new));
	Real c_new = cos(psi2-psi1);
	Real r1_new_integrand = 0.5*K*(1-r1_new*r1_new*r1_new*r1_new)*((1-alpha)*r1_new + alpha*r2_new*c_new) - D1*r1_new;
	Real r2_new_integrand = 0.5*K*(1-r2_new*r2_new*r2_new*r2_new)*(alpha*r2_new + (1-alpha)*r1_new*c_new) - D2*r2_new;

	psi1 += 0.5*(psi1_integrand + psi1_new_integrand)*dt;
	psi2 += 0.5*(psi2_integrand + psi2_new_integrand)*dt;
	r1 += 0.5*(r1_integrand + r1_new_integrand)*dt;
	r2 += 0.5*(r2_integrand + r2_new_integrand)*dt;

	Transform();
}

void Eq(const Real& duration, const Real& dt, const Real& interval, bool time_info_flag = false)
{
	int savestep, steps;
	steps = (int) round(duration / dt);
	if (interval == 0)
		savestep = 200000;
	else
		savestep = (int) round(interval / dt);
	for (int i = 0; i < steps; i++)
	{
		if (i % savestep == 0)
			output << i*dt << "\t" << r1 << "\t" << r2 << "\t" << psi1 << "\t" << psi2 << endl;
		if (i % 100000 == 0 && time_info_flag)
			cout << "\r" << "Step is at: " << i << flush;
		One_Step(dt);
	}
	if (time_info_flag)
		cout << "\tDone" << endl;
}

void Ave(const  int& steps, const Real dt, Real& mr1, Real& mr2)
{
	mr1 = mr2 = 0;
	int counter = 0;
	for (int i = 0; i < steps; i++)
	{
		One_Step(dt);
		if (i % 10 == 0)
		{
			mr1 += r1;
			mr2 += r2;
			counter++;
		}
	}
	mr1 /= counter;
	mr2 /= counter;
}

void Set_Output()
{
	output.close();
	stringstream address;
	address.str("");
	address << "alpha_" << alpha << "_w1_" << std::setprecision(10) << w1 << "_w2_" << w2 << "_D1_" << D1 << "_D2_" << D2 << "_K_" << K << ".dat";
	output.open(address.str().c_str());
	cout << address.str() << endl;
}

void Single_Run(int argc, char* argv[])
{

//	psi1 = 2*M_PI*rand() / RAND_MAX;
//	psi2 = 2*M_PI*rand() / RAND_MAX;
	alpha  = atof(argv[1]);
	w1= atof(argv[2]);
	w2 = atof(argv[3]);
	D1 = atof(argv[4]);
	D2 = atof(argv[5]);
	K = atof(argv[6]);
	Real sim_t = atof(argv[7]);



	r1 = r2 = 0.9;

	Eq(sim_t, 0.001, 0.001,true);
}

void Omega_Chnage(int argc, char* argv[])
{
	alpha  = atof(argv[1]);
	Real omega = atof(argv[2]);
	Real omega_end = atof(argv[3]);
	Real dw = atof(argv[4]);
	Real D = atof(argv[5]);

	Real eq_t = atof(argv[6]);
	Real sim_t = atof(argv[7]);

	K = 1;
	w1 = omega;
	w2 = 0;
	D1 = D2 = D;

	r1 = r2 = 0.9;

	Real dt = 0.001;
	Eq(eq_t, dt, 0.1);

	if (omega_end > omega)
		while (omega <= omega_end)
		{
			w1 = omega;
			w2 = 0;

			Set_Output();
			Eq(sim_t, dt, 0.1);
			omega += dw;
		}
	else
		while (omega >= omega_end)
		{
			w1 = omega;
			w2 = 0;

			Set_Output();
			Eq(sim_t, dt, 0.1);
			omega -= dw;
		}
}

Real tmax = 100000;
Real r2c = 0.7;

Real Go_Till_Spike(const Real& dt)
{
	Real t = 0;
	while (r2 > r2c && t < (tmax+1))
	{
		One_Step(dt);
		t += dt;
	}
	while (r2 < r2c)
	{
		One_Step(dt);
		t += dt;
	}
	return t;
}

Real Spike_Interval(const Real& omega, const Real& D, const Real& eq_t)
{
	w1 = omega;
	D1 = D2 = D;
	r1 = r2 = 0.99;

	Real dt = 0.01;
	Eq(eq_t, dt, 1000);

	Real t = Go_Till_Spike(dt);
	if (t > tmax)
		t = Go_Till_Spike(dt);
	if (t < tmax)
		t = Go_Till_Spike(dt);
	return (t);
}

Real Spike_Interval(int argc, char* argv[])
{
	alpha  = atof(argv[1]);
	Real D = atof(argv[2]);
	Real eq_t = atof(argv[3]);
	Real omega = atof(argv[4]);
	Real factor = atof(argv[5]);
	tmax = atof(argv[6]);
	Real omega_c, domega;

	K = 1;
	w2 = 0;

	Real t = tmax+1;

	while (t > tmax)
	{
		cout << "#" << setprecision(20) << omega << endl;
		omega_c = omega;
		omega *= factor;
		t = Spike_Interval(omega, D, eq_t);
	}
	domega = omega - omega_c;

	cout << setprecision(10) << domega << "\t" << t << endl;

	factor = 1.5;
	while (t < tmax && t > 10)
	{
		domega *= factor;
		omega = omega_c + domega;
		t = Spike_Interval(omega, D, eq_t);
		if (t < tmax)
			cout << setprecision(20) << omega - omega_c << "\t" << t << endl;
	}
}

int main(int argc, char* argv[])
{
	srand(time(NULL));

//	Single_Run(argc, argv);
	Omega_Chnage(argc,argv);
//	Real t = Spike_Interval(argc, argv);

//	Trace_k(100000,100000, dt);
}
