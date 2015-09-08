#include "simulator.h"
#include <cstdlib>

void Single_Run(int argc, char* argv[])
{
	int N1 = atoi(argv[1]);
	int N2  = atoi(argv[2]);
	Real omega1 = atof(argv[3]);
	Real omega2 = atof(argv[4]);
	Real D1 = atof(argv[5]);
	Real D2 = atof(argv[6]);
	Real K = atof(argv[7]);
	Real sim_t = atof(argv[8]);
	Phase::Init_Rand(time(NULL));
	Simulator sim(N1,N2,omega1,omega2,D1,D2,K);
	sim.Run(sim_t, 0.01, 0.05); // simulation time, time steps, saving interval (zero means no save)
}

void Omega_Chnage(int argc, char* argv[])
{
	int N1 = atoi(argv[1]);
	int N2  = atoi(argv[2]);
	Real omega = atof(argv[3]);
	Real omega_end = atof(argv[4]);
	Real dw = atof(argv[5]);
	Real D = atof(argv[6]);
	Real eq_t = atof(argv[7]);
	Real sim_t = atof(argv[8]);
	Phase::Init_Rand(time(NULL));
	Simulator sim(N1,N2,omega,0,D,D,1);
	sim.Run(eq_t, 0.05, 0.1);
	if (omega_end > omega)
		while (omega <= omega_end)
		{
			sim.Set_Omega(omega,0);
			sim.Run(sim_t, 0.05, 0.1);
			omega += dw;
		}
	else
		while (omega >= omega_end)
		{
			sim.Set_Omega(omega,0);
			sim.Run(sim_t, 0.05, 0.1);
			omega -= dw;
		}
		
}

int main(int argc, char* argv[])
{
//	Omega_Chnage(argc, argv);
	Single_Run(argc, argv);
}
