// Abstract base class
// Sets the interface for derived particle classes differentiated by kernel
// Also reffered to as walker sometimes
// NOTE: The KERNEL is a charachteristic of the derived class. Kappa is set by user.

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/log/trivial.hpp>
#include <time.h>
#include <fstream>
#include <unordered_map>
#include <boost/random/discrete_distribution.hpp>
#include <boost/graph/graph_traits.hpp>

typedef boost::random::mt19937 RNG;

template <class Graph>
class Particle
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter;
	typedef typename boost::graph_traits<Graph>::adjacency_iterator adj_iter;

private:
	std::string m_type;
	double m_KernelRate;	//rate at which this particle leaves the current vertex, uninfluenced by other particles
	double m_RedistRate;	//Redistribution or death rate of the particle. 
	double Kappa_M;			//this is kappa/M, the rate at which each particle 'kills' another. 
							//Please remember to divide intended kappa by M, the number of particles of the type in main().

public:
	//=============Constructors and destructor here===================
	Particle()                              //default constructor
	{
		m_type = "TypeA";
		m_KernelRate = 1;
		m_RedistRate = 0;
		Kappa_M = 0;
	}
	Particle(std::string Type, double kappa_M) :m_type(Type), Kappa_M(kappa_M) //constructor with type and kappa
	{
		m_KernelRate = 1;
		m_RedistRate = 0;
	}
	virtual ~Particle() {};                    //destructor
	//================================================================

	//==================Access for private variables==================
	void SetType(const std::string type);
	std::string GetType();
	void SetKernelRate(double rate);
	double GetKernelRate();
	double GetRedistRate();
	void RedistRateIncrease();	//Increase the redist rate by kappa_M. Use this when updating for an incoming particle of another type
	void RedistRateDecrease();	//Decrease the redist rate by kappa_M. Use this when updating for an outgoing particle of another type
	void NewRedistRate(int Competition);	//Sets the new redist rate as kappa_M*Competition. Use this after every movement of the particle.


	//================================================================

	//Purely virtual function (main algorithm)
	virtual void WalkByKernel(Vertex &CurrentVertex, const Graph &G, RNG &gen) = 0;	//perform the random walk according to kernel
	virtual void UpdateKernelRate(const Vertex &CurrentVertex, const Graph &G) = 0;	//update the redist rate according to current position

};


//==========Source Code===============================================
template<class Graph>
void Particle<Graph>::SetType(const std::string type) { m_type = type; }

template<class Graph>
std::string Particle<Graph>::GetType() { return m_type; }

template<class Graph>
void Particle<Graph>::SetKernelRate(double rate) { m_KernelRate = rate; }

template<class Graph>
double Particle<Graph>::GetKernelRate() { return m_KernelRate; }

template<class Graph>
double Particle<Graph>::GetRedistRate() { return m_RedistRate; }

template<class Graph>
void Particle<Graph>::RedistRateIncrease() { m_RedistRate += Kappa_M; }

template<class Graph>
void Particle<Graph>::RedistRateDecrease() { m_RedistRate -= Kappa_M; }

template<class Graph>
void Particle<Graph>::NewRedistRate(int Competition) { m_RedistRate = Kappa_M * Competition; }

#endif