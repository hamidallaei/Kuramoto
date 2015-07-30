#ifndef _OSCILLATOR_
#define _OSCILLATOR_

#include "phase.h"

struct Oscillator{
	Phase phi;
	Real omega;
	Real D;
	Real int_sum; //interaction summation

	static Real K;

	Oscillator(const Real& input_omega, const Real&, const Real&);
	Oscillator(const Real& input_omega, const Real& input_D);
	Oscillator();
	Oscillator(const Oscillator&);

	void Reset();
	void Interact(const Real& factor, const Phase& psi);
	void Add_Drift();
	void Add_Noise(const Real& dt);
	void Noiseless_Euler_Move(const Real&);

	Oscillator& operator= (const Oscillator& c);
};

Oscillator::Oscillator(const Real& input_omega, const Real& input_D, const Real& input_phi) : phi(input_phi), D(input_D), omega(input_omega)
{
	Reset();
}

Oscillator::Oscillator(const Real& input_omega, const Real& input_D) : omega(input_omega), D(input_D)
{
	Reset();
}

Oscillator::Oscillator()
{
	omega = 0;
	D = 0;
	Reset();
}

Oscillator::Oscillator(const Oscillator& os)
{
	phi = os.phi;
	omega = os.omega;
	D = os.D;
	int_sum = os.int_sum;
}

void Oscillator::Reset()
{
	int_sum = 0;
}

void Oscillator::Interact(const Real& factor, const Phase& psi)
{
	int_sum += factor*K*sin(psi - phi);
}

void Oscillator::Add_Drift()
{
	int_sum += omega;
}

void Oscillator::Add_Noise(const Real& dt)
{
	phi.Add_Noise(sqrt(2*D*dt));
}

void Oscillator::Noiseless_Euler_Move(const Real& dt)
{
	phi += int_sum*dt;
}

Oscillator& Oscillator::operator= (const Oscillator& c)
{
	phi = c.phi;
	int_sum = int_sum;
	return (*this);
}

Real Oscillator::K = 1;

#endif