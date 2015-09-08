#ifndef _CLUSTER_
#define _CLUSTER_

#include "oscillator.h"

struct Cluster{
	int N;
	Oscillator* os;
	Real r; //order parameter
	Real psi;
	Real alpha;

	Cluster(int input_N, Real omega, Real D);
	Cluster(int input_N);
	Cluster(const Cluster& c);
	~Cluster();

	void Set_Omega(const Real &);
	void Set_D(const Real &);

	void Reset();
	void Self_Interact();
	void Interact(Cluster& c);
	void All_Interactions(Cluster& c);

	void Add_Noise(const Real& dt);
	void Noiseless_Euler_Move(const Real&);

	void Find_Order_Parameter();

	void Copy_Interaction(Cluster& c);

	void Write_Interactions(ostream& out);

	const Cluster& operator= (const Cluster& c);
	const Cluster& operator+= (const Cluster& c);
};

Cluster::Cluster(int input_N, Real omega = 0, Real D = 0)
{
	N = input_N;
	os = new Oscillator[N];
	Set_Omega(omega);
	Set_D(D);
}

Cluster::Cluster(const Cluster& c)
{
	N = c.N;
	os = new Oscillator[N];
	for (int i = 0; i < N; i++)
		os[i] = c.os[i];
}

Cluster::~Cluster()
{
	if (N != 0)
		delete [] os;
}

void Cluster::Set_Omega(const Real& omega)
{
	for (int i = 0; i < N; i++)
		os[i].omega = omega;
}

void Cluster::Set_D(const Real& D)
{
	for (int i = 0; i < N; i++)
		os[i].D = D;
}

void Cluster::Reset()
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
		os[i].Reset();
}

void Cluster::Self_Interact()
{
	Find_Order_Parameter();
	for (int i = 0; i < N; i++)
	{
		os[i].Interact(alpha*r,psi);
		os[i].Add_Drift();
	}
}

void Cluster::Interact(Cluster& c)
{
	for (int i = 0; i < N; i++)
		os[i].Interact(c.alpha*c.r,c.psi);
}

void Cluster::All_Interactions(Cluster& c)
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
	{
		os[i].Interact(alpha*r,psi);
		os[i].Interact(c.alpha*c.r,c.psi);
		os[i].Add_Drift();
	}
}


void Cluster::Add_Noise(const Real& dt)
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
		os[i].Add_Noise(dt);
}

void Cluster::Noiseless_Euler_Move(const Real& dt)
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
		os[i].Noiseless_Euler_Move(dt);
}

void Cluster::Find_Order_Parameter()
{
	Real C = 0;
	Real S = 0;
	#pragma omp parallel for default(shared) reduction(+:C,S)
	for (int i = 0; i < N; i++)
	{
		C += cos(os[i].phi);
		S += sin(os[i].phi);
	}
	C /= N;
	S /= N;
	r = sqrt(C*C + S*S);
	psi = atan2(S,C);
}

void Cluster::Write_Interactions(ostream& out)
{
	for (int i = 0; i < N; i++)
		out << os[i].int_sum / 2 << endl;
}

const Cluster& Cluster::operator= (const Cluster& c)
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
		os[i] = c.os[i];
	return (*this);
}

const Cluster& Cluster::operator+= (const Cluster& c)
{
	#pragma omp parallel for default(shared)
	for (int i = 0; i < N; i++)
		os[i].int_sum += c.os[i].int_sum;
	return (*this);
}

#endif
