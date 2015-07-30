#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <sstream>
#include <fstream>
#include "cluster.h"

struct Simulator{
	Cluster c1;
	Cluster c2;
	Cluster k1_c1, k1_c2;

	Real t;

	stringstream info;
	ofstream output;

	Simulator(int N1, int N2, Real omega1, Real omega2, Real D1, Real D2, Real K);
	~Simulator();
	
	void One_Step(const Real& dt);
	void One_Step_OMP(const Real& dt);
	void Run(Real duration, const Real& dt, const Real& save_interval);

	void Write();
};

Simulator::Simulator(int N1, int N2, Real omega1, Real omega2, Real D1, Real D2, Real K): c1(N1,omega1,D1), c2(N2,omega2,D2), k1_c1(N1,omega1,D1), k1_c2(N2,omega2,D2), t(0)
{
	Oscillator::K = K;
	c1.alpha = N1;
	c1.alpha /= N1 + N2;
	c2.alpha = 1 - c1.alpha;
	k1_c1.alpha = c1.alpha;
	k1_c2.alpha = c2.alpha;

	info.str("");
	info << "N_" << N1+N2 << "_alpha_" << c2.alpha << "_w1_" << omega1 << "_w2_" << omega2 << "_D1_" << D1 << "_D2_" << D2 << "_K_" << K << ".dat";

	output.open(info.str().c_str());
	cout << info.str() << endl;
	c1.Find_Order_Parameter();
	c2.Find_Order_Parameter();
}

Simulator::~Simulator()
{
	output.close();
}



void Simulator::One_Step(const Real& dt)
{
// 	c1.Find_Order_Parameter();
// 	c2.Find_Order_Parameter();
// 
// 	c1.All_Interactions(c2);
// 	c2.All_Interactions(c1);
// 	k1_c1 = c1;
// 	k1_c2 = c2;
// 
// 	k1_c1.Noiseless_Euler_Move(dt);
// 	k1_c2.Noiseless_Euler_Move(dt);
// 	k1_c1.Reset();
// 	k1_c2.Reset();
// 
// 	k1_c1.Find_Order_Parameter();
// 	k1_c2.Find_Order_Parameter();
// 	k1_c1.All_Interactions(k1_c2);
// 	k1_c2.All_Interactions(k1_c1);
// 
// 	c1 += k1_c1;
// 	c2 += k1_c2;
// 
// 	c1.Noiseless_Euler_Move(0.5*dt);
// 	c2.Noiseless_Euler_Move(0.5*dt);
// 	c1.Reset();
// 	c2.Reset();
// 
// 	c1.Add_Noise(dt);
// 	c2.Add_Noise(dt);

	c1.Find_Order_Parameter();
	c2.Find_Order_Parameter();

	c1.All_Interactions(c2);
	k1_c1 = c1;
	k1_c1.Noiseless_Euler_Move(dt);
	k1_c1.Reset();
	k1_c1.Find_Order_Parameter();

	c2.All_Interactions(c1);
	k1_c2 = c2;
	k1_c2.Noiseless_Euler_Move(dt);
	k1_c2.Reset();
	k1_c2.Find_Order_Parameter();

	k1_c1.All_Interactions(k1_c2);
	c1 += k1_c1;
	c1.Noiseless_Euler_Move(0.5*dt);
	c1.Reset();
	c1.Add_Noise(dt);

	k1_c2.All_Interactions(k1_c1);
	c2 += k1_c2;
	c2.Noiseless_Euler_Move(0.5*dt);
	c2.Reset();
	c2.Add_Noise(dt);

	t += dt;
}

void Simulator::One_Step_OMP(const Real& dt)
{
	Real C = 0;
	Real S = 0;
	#pragma omp parallel for default(shared) reduction(+:C,S)
	for (int i = 0; i < c1.N; i++)
	{
		c1.os[i].Interact(c1.alpha*c1.r,c1.psi);
		c1.os[i].Interact(c2.alpha*c2.r,c2.psi);
		c1.os[i].Add_Drift();
		k1_c1.os[i] = c1.os[i];
		k1_c1.os[i].Noiseless_Euler_Move(dt);
		k1_c1.os[i].Reset();
		C += cos(k1_c1.os[i].phi);
		S += sin(k1_c1.os[i].phi);
	}
	C /= c1.N;
	S /= c1.N;
	k1_c1.r = sqrt(C*C + S*S);
	k1_c1.psi = atan2(S,C);

	C = 0;
	S = 0;
	#pragma omp parallel for default(shared) reduction(+:C,S)
	for (int i = 0; i < c2.N; i++)
	{
		c2.os[i].Interact(c2.alpha*c2.r,c2.psi);
		c2.os[i].Interact(c1.alpha*c1.r,c1.psi);
		c2.os[i].Add_Drift();
		k1_c2.os[i] = c2.os[i];
		k1_c2.os[i].Noiseless_Euler_Move(dt);
		k1_c2.os[i].Reset();
		C += cos(k1_c2.os[i].phi);
		S += sin(k1_c2.os[i].phi);
	}
	C /= c2.N;
	S /= c2.N;
	k1_c2.r = sqrt(C*C + S*S);
	k1_c2.psi = atan2(S,C);

	C = 0;
	S = 0;
	#pragma omp parallel for default(shared) reduction(+:C,S)
	for (int i = 0; i < c1.N; i++)
	{
		k1_c1.os[i].Interact(c1.alpha*k1_c1.r,k1_c1.psi);
		k1_c1.os[i].Interact(c2.alpha*k1_c2.r,k1_c2.psi);
		k1_c1.os[i].Add_Drift();
		c1.os[i].int_sum += k1_c1.os[i].int_sum;
		c1.os[i].Noiseless_Euler_Move(0.5*dt);
		c1.os[i].Reset();
		c1.os[i].Add_Noise(dt);
		C += cos(c1.os[i].phi);
		S += sin(c1.os[i].phi);
	}
	C /= c1.N;
	S /= c1.N;
	c1.r = sqrt(C*C + S*S);
	c1.psi = atan2(S,C);

	C = 0;
	S = 0;
	#pragma omp parallel for default(shared) reduction(+:C,S)
	for (int i = 0; i < c2.N; i++)
	{
		k1_c2.os[i].Interact(c2.alpha*k1_c2.r,k1_c2.psi);
		k1_c2.os[i].Interact(c1.alpha*k1_c1.r,k1_c1.psi);
		k1_c2.os[i].Add_Drift();
		c2.os[i].int_sum += k1_c2.os[i].int_sum;
		c2.os[i].Noiseless_Euler_Move(0.5*dt);
		c2.os[i].Reset();
		c2.os[i].Add_Noise(dt);
		C += cos(c2.os[i].phi);
		S += sin(c2.os[i].phi);
	}
	C /= c2.N;
	S /= c2.N;
	c2.r = sqrt(C*C + S*S);
	c2.psi = atan2(S,C);

	t += dt;
}

void Simulator::Run(Real duration, const Real& dt, const Real& save_interval = 0)
{
	long int steps = (int) round(duration / dt);
	int save_steps = (int) round(save_interval / dt);
	if (save_interval == 0)
		save_steps = 0;

	c1.Find_Order_Parameter();
	c2.Find_Order_Parameter();

	for (int i = 0; i < steps; i++)
	{
		if ((i % save_steps == 0) && (save_steps != 0))
		{
			if (i % (1000*save_steps) == 0)
				cout << "time is: " << t << endl;
			Write();
		}
		One_Step_OMP(dt);
//		One_Step(dt);
		}
}

void Simulator::Write()
{
	c1.Find_Order_Parameter();
	c2.Find_Order_Parameter();
	output << t << "\t" << c1.r << "\t" << c2.r << "\t" << c1.psi << "\t" << c2.psi << endl;
}

#endif
