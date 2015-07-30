#include "simulator.h"
#include <cstdlib>

int main(int argc, char* argv[])
{
	int N1 = atoi(argv[1]);
	int N2  = atoi(argv[2]);
	Real omega1 = atof(argv[3]);
	Real omega2 = atof(argv[4]);
	Real D1 = atof(argv[5]);
	Real D2 = atof(argv[6]);
	Real K = atof(argv[7]);
	Real sim_t = atof(argv[8]);
	Phase::Init_Rand(124);
	Simulator sim(N1,N2,omega1,omega2,D1,D2,K);
	sim.Run(sim_t, 0.01, 0.1);
}
