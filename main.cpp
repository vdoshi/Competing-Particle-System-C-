#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/log/trivial.hpp>
#include <vector>
#include <string>
#include "Network.hpp"
#include "Particle.hpp"
#include "CombinatorialKernel.hpp"
#include <cmath>

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, VertexInfo> unDirAdjList;

int main() {
	for (double t = 0; t <= 0; ++t)
	{
		double Tau = t * 20;
		for (double k = 1; k <= 9; ++k)
		{
			double kappa = k;// std::pow(10, k);
			//double kappa = 0;
			std::vector<CombParticle<unDirAdjList>> particles;
			for (int i = 0; i < 200; ++i) //adjust number of particles here
			{

				if (i % 2 == 0) {
					CombParticle<unDirAdjList> temp("dodgerblue", kappa / 100);
					particles.push_back(temp);
				}

				else
				{
					CombParticle<unDirAdjList> temp("orangered", kappa / 100);
					particles.push_back(temp);
				}

			}
			//std::cout << "particles initialized" << std::endl;
			std::vector<Particle<unDirAdjList>*> particle_ptr;
			for (int i = 0; i < 200; ++i) //adjust number of particles here
			{
				particle_ptr.push_back(&particles[i]);
			}

			//std::ofstream RR_out;
			//RR_out.open("kappa" + std::to_string(kappa) + "_RR.txt");

			for (int i = 1; i <= 10; ++i)
			{
				Network<unDirAdjList> test("facebook_combined.txt", particle_ptr);
				test.WalkAndLog_Hist(i, 250, kappa, Tau);
			}

		}
		//RR_out.close();

	}
	system("pause");
	return 0;
}