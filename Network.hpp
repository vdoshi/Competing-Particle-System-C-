// Use this to read the graph and prepare the simulations etup for the graph.
// Only use the adjacency list data structure for now.
// Works with the continuous time implementation of the markov chain

#ifndef NETWORK_H
#define NETWORK_H

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <stdio.h>
#include <cmath>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/expressions.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>

#include "Particle.hpp"



//====Convert ~U[0,1] to ~exp(\lambda)
double GenExpTime(double rate, double u)
{
	return -(1 / rate)*log(1 - u);
};

//====Logfile settings=========================================//
void init_trace_logging(std::string &filename)
{
	boost::log::add_file_log(filename, boost::log::keywords::auto_flush = true);
	boost::log::core::get()->set_filter
	(
		boost::log::trivial::severity == boost::log::trivial::trace
	);
};

void stop_trace_logging() {
	boost::log::core::get()->remove_all_sinks();
}

//====Property bundle for node information====================//
struct VertexInfo
{
	unsigned int ID;
	std::string Type;
	int Density;
};

typedef boost::random::mt19937 RNG;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, VertexInfo> unDirAdjList;

//===="Label writer" functor for writing graphviz output======//
//====Edit this to edit the graph output style/info===========//
template <class IDMAP, class TYPEMAP, class DENSITYMAP>
class vertex_writer
{
private:
	IDMAP id_m;
	TYPEMAP type_m;
	DENSITYMAP density_m;

public:
	vertex_writer(IDMAP id, TYPEMAP type, DENSITYMAP density) :id_m(id), type_m(type), density_m(density) {}  //constructor

	template <class Vertex>
	void operator()(std::ostream &out, const Vertex &v) const //main operator overloading
	{//print the node information in dot language
		out << "[label=\"" << id_m[v] << "\", color=\"" << type_m[v] << "\", shape=circle, style=filled,  width=\"" << (double(density_m[v] + 1)) / 100 << "\"]";
	}
};

//Initialize the functor for "global" use
template <class IDMAP, class TYPEMAP, class DENSITYMAP>
inline vertex_writer<IDMAP, TYPEMAP, DENSITYMAP>
make_vertex_writer(IDMAP id, TYPEMAP type, DENSITYMAP density) {
	return vertex_writer<IDMAP, TYPEMAP, DENSITYMAP>(id, type, density);
}

//=========================================================================================================================================================//
//==========The Main Network class=========================================================================================================================//
template <class Graph = unDirAdjList>
class Network
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

	typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter;
	typedef typename boost::graph_traits<Graph>::adjacency_iterator adj_iter;

	typedef typename std::map<double, std::vector<Vertex>> MAP_timeHistory;

private:
	Graph m_Graph;						//This is the adjacency list (or any other Graph data structure)
	std::vector<std::string> m_types;	//Vector of all different particle types
	double m_NetKernelRate;				//Rate at which the next "walk event" occurs
	double m_NetRedistRate;				//Rate at which the next "redistribution event" occurs
	std::vector<double> vecKernelRates;	//vector of all kernel rates, ordered according to ordering of particles
	std::vector<double> vecRedistRates;	//vector of all redist rates, ordered according to ordering of particles
	std::vector<Particle<Graph>*> m_particles;	//vector of all particles

	std::unordered_map<std::string, std::vector<Particle<Graph>*> > type2particles;	//Contains arrays of particles seggregated by their type
	std::unordered_map<int, Vertex> name2vertex;					//Map to retrieve vertex descriptor from its ID
	std::unordered_map<Particle<Graph>*, Vertex> particle2vertex;   //Map to retrieve location of a particle

	std::unordered_map<std::string, MAP_timeHistory> type2history;	//Contains all the history. Updated using the TraceLog_Hist method
	std::unordered_map<std::string, std::vector<double>> type2HistSum;	//contains history weighted by time and summed up


	void Redist(Particle<Graph>* Ptr, RNG &gen);		//Redistributes one of the particles by density distribution. Updates the m_RedistRate for all particles involved. Updates m_KernelRate for redistributed particle.
	void Redist_Hist(Particle<Graph>* Ptr, RNG &gen, double time);	//Same as Redist but redistribution vector is the historical distribution at some past time
	void WalkAccToKernel(Particle<Graph>* Ptr, RNG &gen);	//Calls the walk by kernel procedure for the particle. Updates the m_RedistRate for all the particles involved. Updates m_KernelRate for redistributed particle.
	std::vector<Particle<Graph>*> FindLocalCompetition(Particle<Graph>* Ptr);	//Returns the vector of particles of different type located at the same vertex
	void IncreaseRedistRates(std::vector<Particle<Graph>*> &particles);	//given a vector of particles, increase the redistribution rate for all of them
	void DecreaseRedistRates(std::vector<Particle<Graph>*> &particles);	//given a vector of particles, decrease the redistribution rate for all of them

	void TraceLog(const double time, const double time_previous, const int out, const int event_count);	//Logs the trace of the algorithm
	void TraceLog_Hist(const double time, const double time_previous, const int out, const int event_count, const double TAU);	//Works same as tracelog but also updates the History variable
	void UpdateNetRates();				//Updates the rate at which events occur. Should be called immediately after every step

public:
	//Not yet default constructible.
	//Construct only using filename for reading graph and a vector of particle pointers
	Network(std::string filename, std::vector<Particle<Graph>*> particles);
	~Network() {};   //destructor

	void WalkAndLog(const int run_no, const double T, const double kappa);	//The MAIN function
	void WalkAndLog_Hist(const int run_no, const double T, const double kappa, const double TAU);	//Same as WalkAndLog but uses Historical average for redistribution
};



//=========Source Code=========================================//

//Constructor
template <class Graph>
Network<Graph>::Network(std::string filename, std::vector<Particle<Graph>*> particles)
{
	typedef unsigned int uint;
	BOOST_LOG_TRIVIAL(info) << "Creating graph";
	std::ifstream input(filename);
	std::string line;
	int node_count = 0;
	int edge_count = 0;
	//-----------This while loop will create the graph----------
	while (getline(input, line))
	{

		uint source_name = 0;
		uint target_name = 0;
		if (line.size() && line[0] == '#')
			continue;
		else
		{

			sscanf(line.c_str(), "%d\t%d", &source_name, &target_name);	//read source and target node ids from the text file

			//make sure a vertex for the source has been added to the graph
			Vertex s;
			if (name2vertex.find(source_name) == name2vertex.end())
			{
				s = add_vertex(m_Graph);
				node_count = node_count + 1;
				//BOOST_LOG_TRIVIAL(info) << "Node count: " << node_count;
				name2vertex.insert(std::make_pair(source_name, s));
				m_Graph[s].ID = source_name;
				m_Graph[s].Type = "gainsboro";	//this basically means empty
				m_Graph[s].Density = 0;
			}
			else
				s = name2vertex[source_name];

			//make sure a vertex for the target node has been added to the graph
			Vertex t;
			if (name2vertex.find(target_name) == name2vertex.end())
			{
				t = add_vertex(m_Graph);
				node_count = node_count + 1;
				//BOOST_LOG_TRIVIAL(info) << "Node count: " << node_count;
				name2vertex.insert(std::make_pair(target_name, t));
				m_Graph[t].ID = target_name;
				m_Graph[t].Type = "gainsboro";
				m_Graph[t].Density = 0;
			}
			else
				t = name2vertex[target_name];

			//add the edge to the graph
			boost::add_edge(s, t, m_Graph);
			++edge_count;
			//BOOST_LOG_TRIVIAL(info) << "Edge count: "<<edge_count;
		}
	}
	//--------------------------------------------------------------
	BOOST_LOG_TRIVIAL(info) << "Graph created. Counting particle types and assigning initial positions";
	//counting particle types and assigning initial positions at random
	//Note: particles of different types should not be assigned the same initial position
	//Assuming that the redist rate is always zero intiaiaklly for all particles

	//--------------------------------------------------------------
	RNG gen(std::time(0));	//random number generator
	v_iter vi_start, vi_end, vi;
	boost::tie(vi_start, vi_end) = boost::vertices(m_Graph);

	//iterate through all particles
	for (typename std::vector<Particle<Graph>*>::iterator p_it = particles.begin(); p_it != particles.end(); ++p_it)
	{
		//counting down and recording the different particle types
		std::string particle_type = (*p_it)->GetType();
		if (find(m_types.begin(), m_types.end(), particle_type) == m_types.end())
		{
			m_types.push_back(particle_type);
			std::vector<Particle<Graph>*> temp;
			type2particles.insert(std::make_pair(particle_type, temp));
		}
		type2particles[particle_type].push_back(*p_it);

		//assigning the particle to a vertex
		int numberOfVertices = boost::num_vertices(m_Graph);

		boost::random::uniform_int_distribution<> moves(0, numberOfVertices - 1);
		vi = vi_start;

		std::advance(vi, moves(gen));

		//making sure it not assigned to a vertex containing particles of a different type

		while ((m_Graph[*vi].Type != particle_type) && (m_Graph[*vi].Type != "gainsboro"))
		{
			vi = vi_start;
			std::advance(vi, moves(gen));
		}

		//at the exit of the while loop, vi will be pointing to a feasible vertex
		//add this feasible vertex to the map recording the location of the particle
		//also update the density and type in the vertex property bundle
		//update the Kernel rate of the particle under contention
		particle2vertex.insert(std::make_pair(*p_it, *vi));
		m_Graph[*vi].Type = particle_type;
		m_Graph[*vi].Density = m_Graph[*vi].Density + 1;
		(*p_it)->UpdateKernelRate(*vi, m_Graph);
		//BOOST_LOG_TRIVIAL(info) << "Initial kernel rate of particle = " << (*p_it)->GetKernelRate();
		m_particles.push_back(*p_it);	//add the particle to list of all particles
	}
	//---------------------------------------------------------------
	UpdateNetRates();
	BOOST_LOG_TRIVIAL(info) << "Net initial kernel rate = " << m_NetKernelRate;
	BOOST_LOG_TRIVIAL(info) << "Done. Finished all setups";
}


//Walk and log
template <class Graph>
void Network<Graph>::WalkAndLog(const int run_no, const double T, const double kappa)
{
	double time_count = 0.0;
	//Initialize zero hist sum
	int n_vertices = name2vertex.size();

	for (std::vector<std::string>::iterator type_it = m_types.begin(); type_it != m_types.end(); ++type_it)
	{
		for (int id = 0; id < n_vertices; ++id)
		{
			type2HistSum[*type_it].push_back(0.0);
		}
	}
	//output the initial graph
	std::cout << "Run: " << run_no << std::endl;
	double time_previous = 0;
	double time = 0; //will use this as stopping criteria
					 /*
					 std::ofstream dotfile("InitialGraphRun" + std::to_string(run_no) + ".dot");
					 boost::write_graphviz(dotfile, m_Graph,
					 make_vertex_writer(boost::get(&VertexInfo::ID, m_Graph),
					 boost::get(&VertexInfo::Type, m_Graph),
					 boost::get(&VertexInfo::Density, m_Graph)));
					 */

					 //initialize trace logging
	std::string filename = "DistfacebookParticles" + std::to_string(m_particles.size()) + "_Kappa" + std::to_string(kappa) + "_traceRun" + std::to_string(run_no) + ".log";
	init_trace_logging(filename);
	TraceLog(time, time_previous, 1, 0);

	RNG gen(std::time(0));
	boost::random::uniform_01<> unif;

	//----------------MAIN LOOP----------------
	int event_count = 0;
	while (time < T)
	{
		//first, generate time at which next event occurs
		//time += GenExpTime(unif(gen));
		//std::cout << "finding time for next event...\n";
		boost::random::exponential_distribution<> EXP(m_NetKernelRate + m_NetRedistRate);
		time += EXP(gen);
		//std::cout << "next event @ time " <<time <<" \n";
		//next, generate the type of event
		double p = unif(gen);
		if (p < m_NetRedistRate / (m_NetKernelRate + m_NetRedistRate)) //i.e. if event is a redistribution. Now generate whom to redistribute
		{
			//std::cout << "Redistributing..." << std::endl;
			boost::random::discrete_distribution<> dist(vecRedistRates.begin(), vecRedistRates.end());
			int particle = dist(gen);
			Redist(m_particles[particle], gen);	//redistribute the respective particle. This also updates the rates
		}
		else	//i.e. event is a walk. Now generate whom to walk
		{
			//std::cout << "Walking..." << std::endl;
			boost::random::discrete_distribution<> dist(vecKernelRates.begin(), vecKernelRates.end());
			int particle = dist(gen);
			WalkAccToKernel(m_particles[particle], gen);	//walk the corressponding particle acc to algoroithm. This also updates the kernel and redist rates of particles involved]
		}
		int out;
		if (time > time_count)
		{
			out = 1;
			time_count += 0.1;
		}
		else
			out = 0;

		event_count++;
		TraceLog(time, time_previous, out, event_count);	//log the current particle positions
		time_previous = time;
		//std::cout << "Run: " << run_no << ", Event count = " << event_count << ", RR = "<< m_NetRedistRate << std::endl;
	}
	std::cout << "Run: " << run_no << ", Event count = " << event_count << std::endl;
	stop_trace_logging();
}


//Walk and log with history
template <class Graph>
void Network<Graph>::WalkAndLog_Hist(const int run_no, const double T, const double kappa, const double TAU)
{
	double time_count = 0.0;
	//Initialize zero hist sum
	int n_vertices = name2vertex.size();

	for (std::vector<std::string>::iterator type_it = m_types.begin(); type_it != m_types.end(); ++type_it)
	{
		for (int id = 0; id < n_vertices; ++id)
		{
			type2HistSum[*type_it].push_back(0.0);
		}
	}

	//in this function, will have to keep a record of all the history over all times.
	//************FINISH THIS PART*************
	//output the initial graph
	std::cout << "Run: " << run_no << std::endl;
	double time_previous = 0;
	double time = 0; //will use this as stopping criteria
					 /*
					 std::ofstream dotfile("InitialGraphRun" + std::to_string(run_no) + ".dot");
					 boost::write_graphviz(dotfile, m_Graph,
					 make_vertex_writer(boost::get(&VertexInfo::ID, m_Graph),
					 boost::get(&VertexInfo::Type, m_Graph),
					 boost::get(&VertexInfo::Density, m_Graph)));
					 */

					 //initialize trace logging
	std::string filename = "Hist_Distfacebook_Particles" + std::to_string(m_particles.size()) + "_Kappa" + std::to_string(kappa) + "_Tau" + std::to_string(TAU) + "_traceRun" + std::to_string(run_no) + ".log";
	init_trace_logging(filename);
	TraceLog_Hist(time, time_previous, 1, 0, TAU);

	RNG gen(std::time(0));
	boost::random::uniform_01<> unif;

	//----------------MAIN LOOP----------------
	int event_count = 0;
	while (time < T)
	{

		//first, generate time at which next event occurs
		//time += GenExpTime(unif(gen));
		//std::cout << "finding time for next event...\n";
		boost::random::exponential_distribution<> EXP(m_NetKernelRate + m_NetRedistRate);
		time += EXP(gen);
		//std::cout << "next event @ time " <<time <<" \n";
		//next, generate the type of event
		double p = unif(gen);
		if (p < m_NetRedistRate / (m_NetKernelRate + m_NetRedistRate)) //i.e. if event is a redistribution. Now generate whom to redistribute
		{
			double Redist_time;
			//HERE: Need to first identify the time in history we are redistributing to.
			if (TAU <= time)
				Redist_time = time - unif(gen)*TAU; //redist uniformly over the last TAU time units. CAN BE MODIFIED
			else
				Redist_time = unif(gen)*time;


			//std::cout << "Redistributing..." << std::endl;
			boost::random::discrete_distribution<> dist(vecRedistRates.begin(), vecRedistRates.end());
			int particle = dist(gen);
			Redist_Hist(m_particles[particle], gen, Redist_time);	//redistribute the respective particle. This also updates the rates
		}
		else	//i.e. event is a walk. Now generate whom to walk
		{
			//std::cout << "Walking..." << std::endl;
			boost::random::discrete_distribution<> dist(vecKernelRates.begin(), vecKernelRates.end());
			int particle = dist(gen);
			WalkAccToKernel(m_particles[particle], gen);	//walk the corressponding particle acc to algoroithm. This also updates the kernel and redist rates of particles involved]
		}
		int out;
		if (time > time_count)
		{
			out = 0;
			//--temp code--
			if (time_count == 250)
				out = 1;
			else
				out = 0;
			//--temp code--
			time_count += 1.0;
		}
		else
			out = 0;

		event_count++;
		TraceLog_Hist(time, time_previous, out, event_count, TAU);	//log the current particle positions
		time_previous = time;
	}
	std::cout << "Run: " << run_no << ", Event count = " << event_count << "; KR = " << m_NetKernelRate << ", RR = " << m_NetRedistRate << std::endl;
	stop_trace_logging();
	type2history.erase(type2history.begin(), type2history.end());
}


//Redist
template <class Graph>
void Network<Graph>::Redist(Particle<Graph>* Ptr, RNG &gen)
{

	//Here, it is important for every particle involved to go through update redist rate procedure.
	std::string particle_type = Ptr->GetType();     //particle type
	unsigned int current_node = m_Graph[particle2vertex[Ptr]].ID;   //current position
	std::vector<Particle<Graph>*>::iterator it_s = type2particles[particle_type].begin();
	boost::random::uniform_int_distribution<> DensityDist(0, type2particles[particle_type].size() - 1);
	std::vector<Particle<Graph>*>::iterator it;

	int flag = 1;   //will be 0 when a feasible vertex is found

	while (flag == 1)   //while a feasible vertex is not found, redistribute
	{

		it = it_s;
		std::advance(it, DensityDist(gen));
		if (m_Graph[particle2vertex[*it]].ID != current_node)
			flag = 0;
	}

	DecreaseRedistRates(FindLocalCompetition(Ptr));	//adjust the redist rate of the old competition
	particle2vertex[Ptr] = particle2vertex[*it];	//move particle to feasible vertex
	std::vector<Particle<Graph>*> NewComp = FindLocalCompetition(Ptr);	//find new competition
	IncreaseRedistRates(NewComp);	//adjust the redist rate of the new competition
	Ptr->UpdateKernelRate(particle2vertex[Ptr], m_Graph);	//update the kernel rate
	Ptr->NewRedistRate(NewComp.size());
	//update the rates
	UpdateNetRates();
}


//Redist according to history dist
template <class Graph>
void Network<Graph>::Redist_Hist(Particle<Graph>* Ptr, RNG &gen, double time)
{
	//Here, it is important for every particle involved to go through update redist rate procedure.
	//Also, need to do a (backward) search to find the correct historical distribution
	std::string particle_type = Ptr->GetType();     //particle type
	unsigned int current_node = m_Graph[particle2vertex[Ptr]].ID;   //current position
	MAP_timeHistory::reverse_iterator t_it;
	double inter_time = time;
	for (t_it = type2history[particle_type].rbegin(); t_it != type2history[particle_type].rend(); ++t_it)
	{
		if (time >= t_it->first)
		{	//found the time of the transition just before the time we want to go to
			inter_time = t_it->first;
			break;
		}
	}
	std::vector<Vertex>::iterator it_s = t_it->second.begin();
	boost::random::uniform_int_distribution<> DensityDist_Hist(0, t_it->second.size() - 1);
	std::vector<Vertex>::iterator it;

	int flag = 1;   //will be 0 when a feasible vertex is found

	while (flag == 1)   //while a feasible vertex is not found, redistribute
	{

		it = it_s;
		std::advance(it, DensityDist_Hist(gen));
		if (*it != current_node)
			flag = 0;
	}
	//std::cout << "(Redist_time, inter_time): (" << time << ", " << inter_time << ")\n";
	DecreaseRedistRates(FindLocalCompetition(Ptr));	//adjust the redist rate of the old competition
	particle2vertex[Ptr] = *it;	//move particle to feasible vertex
	std::vector<Particle<Graph>*> NewComp = FindLocalCompetition(Ptr);	//find new competition
	IncreaseRedistRates(NewComp);	//adjust the redist rate of the new competition
	Ptr->UpdateKernelRate(particle2vertex[Ptr], m_Graph);	//update the kernel rate
	Ptr->NewRedistRate(NewComp.size());
	//update the rates
	UpdateNetRates();
}

//Walk according to Kernel
template <class Graph>
void Network<Graph>::WalkAccToKernel(Particle<Graph>* Ptr, RNG &gen)
{

	DecreaseRedistRates(FindLocalCompetition(Ptr));	//adjust the redist rate of the old competition
	Ptr->WalkByKernel(particle2vertex[Ptr], m_Graph, gen);	//walk according to kernel
	std::vector<Particle<Graph>*> NewComp = FindLocalCompetition(Ptr);	//find new competition


	Ptr->UpdateKernelRate(particle2vertex[Ptr], m_Graph);	//update the kernel rate
	int rr = Ptr->GetRedistRate();
	Ptr->NewRedistRate(NewComp.size());
	IncreaseRedistRates(NewComp);	//adjust the redist rate of the new competition
	//update the rates
	UpdateNetRates();
}


//Increase Redistribution Rates
template <class Graph>
void Network<Graph>::IncreaseRedistRates(std::vector<Particle<Graph>*> &particles)
{
	for (typename std::vector<Particle<Graph>*>::iterator it = particles.begin(); it != particles.end(); ++it)
		(*it)->RedistRateIncrease();

}


//Decrease Redistribution Rates
template <class Graph>
void Network<Graph>::DecreaseRedistRates(std::vector<Particle<Graph>*> &particles)
{
	for (typename std::vector<Particle<Graph>*>::iterator it = particles.begin(); it != particles.end(); ++it)
		(*it)->RedistRateDecrease();
}


//Find local competition
template <class Graph>
std::vector<Particle<Graph>*> Network<Graph>::FindLocalCompetition(Particle<Graph>* Ptr)
{
	std::vector<Particle<Graph>*> competition;
	unsigned int CurrentID = m_Graph[particle2vertex[Ptr]].ID;
	std::string particle_type = Ptr->GetType();

	for (std::vector<Particle<Graph>*>::iterator p_it = m_particles.begin(); p_it != m_particles.end(); ++p_it)
	{
		if ((particle_type != (*p_it)->GetType()) && (CurrentID == m_Graph[particle2vertex[*p_it]].ID))
		{
			competition.push_back(*p_it);
		}
	}
	return competition;
}


//Update Net Rate
template <class Graph>
void Network<Graph>::UpdateNetRates()
{
	m_NetKernelRate = 0;
	m_NetRedistRate = 0;
	vecKernelRates.clear();
	vecRedistRates.clear();
	for (typename std::vector<Particle<Graph>*>::iterator p_it = m_particles.begin(); p_it != m_particles.end(); ++p_it)
	{
		double KR = (*p_it)->GetKernelRate();
		double RR = (*p_it)->GetRedistRate();

		m_NetKernelRate += KR;
		m_NetRedistRate += RR;
		vecKernelRates.push_back(KR);
		vecRedistRates.push_back(RR);
	}

}


//Trace log
template <class Graph>
void Network<Graph>::TraceLog(const double time, const double time_previous, const int out, const int event_count)
{
	std::string TraceString = "# Time = " + std::to_string(time) + "\n";
	for (std::vector<std::string>::iterator type_it = m_types.begin(); type_it != m_types.end(); ++type_it)
	{
		std::vector<Vertex> temp;
		TraceString += *type_it + ": [";
		for (typename std::vector<Particle<Graph>*>::iterator p_it = type2particles[*type_it].begin(); p_it != type2particles[*type_it].end(); ++p_it)
		{

			type2HistSum[*type_it][m_Graph[particle2vertex[*p_it]].ID] += 1 * (time - time_previous);
			/*//Comment out if want to log exact trace and not the final history dist

			if (p_it == type2particles[*type_it].end() - 1)
				TraceString += std::to_string(m_Graph[particle2vertex[*p_it]].ID);
			else
				TraceString += std::to_string(m_Graph[particle2vertex[*p_it]].ID) + ", ";

			*/
		}

		if (out != 0)
		{
			for (std::vector<double>::iterator h_it = type2HistSum[*type_it].begin(); h_it != type2HistSum[*type_it].end(); ++h_it)
			{
				if (h_it == type2HistSum[*type_it].end() - 1)
					TraceString += std::to_string(*h_it);
				else
					TraceString += std::to_string(*h_it) + ", ";
			}
		}
		TraceString += "]\n";
		TraceString += "Event Count = " + std::to_string(event_count) + "\n";
	}
	if (out != 0)
		BOOST_LOG_TRIVIAL(trace) << TraceString;
}

//Trace log with history updated
template <class Graph>
void Network<Graph>::TraceLog_Hist(const double time, const double time_previous, const int out, const int event_count, const double TAU)
{
	std::string TraceString = "# Time = " + std::to_string(time) + "\n";
	for (std::vector<std::string>::iterator type_it = m_types.begin(); type_it != m_types.end(); ++type_it)
	{
		std::vector<Vertex> temp;

		TraceString += *type_it + ": [";
		for (typename std::vector<Particle<Graph>*>::iterator p_it = type2particles[*type_it].begin(); p_it != type2particles[*type_it].end(); ++p_it)
		{
			temp.push_back(particle2vertex[*p_it]);
			type2HistSum[*type_it][m_Graph[particle2vertex[*p_it]].ID] += 1 * (time - time_previous);
			/*//Comment out if want to log exact trace and not the final history dist

			if (p_it == type2particles[*type_it].end() - 1)
				TraceString += std::to_string(m_Graph[particle2vertex[*p_it]].ID);
			else
				TraceString += std::to_string(m_Graph[particle2vertex[*p_it]].ID) + ", ";

			*/
		}
		type2history[*type_it].insert(std::make_pair(time, temp));	//record history and map it to time

		if (out != 0)
		{
			for (std::vector<double>::iterator h_it = type2HistSum[*type_it].begin(); h_it != type2HistSum[*type_it].end(); ++h_it)
			{
				if (h_it == type2HistSum[*type_it].end() - 1)
					TraceString += std::to_string(*h_it);
				else
					TraceString += std::to_string(*h_it) + ", ";
			}
		}
		TraceString += "]\n";
		TraceString += "Event Count = " + std::to_string(event_count) + "\n";

		if (time > TAU)
		{	//erasing history tail as time moves forward
			MAP_timeHistory::iterator time_it = type2history[*type_it].begin();
			while (time_it->first < time - TAU)
				++time_it;

			type2history[*type_it].erase(type2history[*type_it].begin(), --time_it);
		}

	}
	if (out != 0)
		BOOST_LOG_TRIVIAL(trace) << TraceString;
}

#endif