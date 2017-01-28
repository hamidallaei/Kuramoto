#include "simulator.h"
#include <cstdlib>

void Single_Run(int argc, char* argv[])
{
	int N1 = atoi(argv[1]);
	int N2  = atoi(argv[2]);
	Real D = atof(argv[3]);
	Real alpha = atof(argv[4]);
	Real sim_t = atof(argv[5]);
	Phase::Init_Rand(time(NULL));
	Simulator sim(1.0/32,N1,N2,0,0,D,1);
	sim.Set_kappa(0.5);
	sim.Set_alpha(alpha);
	sim.Set_Output();

	sim.Run(sim_t, 2*Oscillator::dt); // simulation time, time steps, saving interval (zero means no save)
}

void Change_Noise(int argc, char* argv[])
{
	int N1 = atoi(argv[1]);
	int N2  = atoi(argv[2]);
	Real alpha = atof(argv[3]);
	Real sim_t = atof(argv[4]);
	Phase::Init_Rand(time(NULL));
	Simulator sim(1.0/32, N1, N2, 0, 0, 0, 1);
	sim.Set_kappa(0.5);
	sim.Set_alpha(alpha);

	Real D_max = 0.5*(1 + sim.c1.kappa)*cos(alpha);
	Real D = 0;
	for (D = 0; D < 1; D += 0.01)
	{
		sim.Set_D(D);
		sim.Set_Output();
		sim.Run(sim_t, 2*Oscillator::dt); // simulation time, time steps, saving interval (zero means no save)
	}
}

int main(int argc, char* argv[])
{
//	Omega_Chnage(argc, argv);
//	Single_Run(argc, argv);
	Change_Noise(argc, argv);
}
