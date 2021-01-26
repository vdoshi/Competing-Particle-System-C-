// Derived class
// Makes the particle walk according to kenrel L=D-A

#ifndef COMBINATORIALKERNEL_H
#define COMBINATORIALKERNEL_H

#include "Particle.hpp"

template <class Graph>
class CombParticle : public Particle<Graph>
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter;

public:
	CombParticle():Particle<Graph>() {};
	~CombParticle() {};
	CombParticle(std::string Type, double kappa) :Particle<Graph>(Type, kappa) {};

	void WalkByKernel(Vertex &CurrentVertex, const Graph& G, RNG &gen);
	void UpdateKernelRate(const Vertex &CurrentVertex, const Graph &G);

};

template <class Graph>
void CombParticle<Graph>::WalkByKernel(Vertex &CurrentVertex, const Graph& G, RNG &gen)
{
	std::vector<Vertex> adj_vertices;
	adj_iter vi, vi_end;
	for (boost::tie(vi, vi_end) = boost::adjacent_vertices(CurrentVertex, G); vi != vi_end; ++vi)
	{
		adj_vertices.push_back(*vi);
	}
	typename std::vector<Vertex>::iterator it = adj_vertices.begin();
	if (adj_vertices.size() >1)
	{
		boost::random::uniform_int_distribution<> moves(0, adj_vertices.size() - 1);
		std::advance(it, moves(gen));
	}
	CurrentVertex = *it;
}

template <class Graph>
void CombParticle<Graph>::UpdateKernelRate(const Vertex &CurrentVertex, const Graph &G)
{
	Particle<Graph>::SetKernelRate(boost::degree(CurrentVertex, G));
}

#endif