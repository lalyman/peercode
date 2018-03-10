/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

//damping constant
static double c;


/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
	double og_length;
	double K;
	EdgeData() : og_length(1.0), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
	
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
	if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0))
    n.position() += n.value().vel * dt;
  }
	
  //Constraints
	constraints(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
	if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0))
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
//struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
 /* template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
	if (n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0)) {
		return Point(0 ,0 ,0);
	} else {
		Point f(0,0,0);
		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
		f += -(*it).value().K*(n.position() - (*it).node2().position())/(*it).length()*((*it).length() - (*it).value().og_length);
		}
		f += Point(0,0,-n.value().mass*grav);
    		return f;
	}
  }
};*/
struct GravityForce {
template <typename NODE>
  Point operator()(NODE n, double t) {
	return Point(0,0,-n.value().mass*grav);
  }
};

struct MassSpringForce {
template <typename NODE>
  	Point operator()(NODE n, double t) {
	Point f(0,0,0);
	for (auto it = n.edge_begin(); it != n.edge_end(); ++it) 
		f += -(*it).value().K*(n.position() - (*it).node2().position())/(*it).length()*((*it).length() - (*it).value().og_length);
    	return f;
	}
};

struct DampingForce {
template <typename NODE>
  	Point operator()(NODE n, double t) {
    	return -c*n.value().vel;
	}
};

struct NullForce {
template <typename NODE>
  	Point operator()(NODE n, double t) {
    	return Point(0,0,0);
	}
};

template <class T1, class T2, class T3 = NullForce>
struct CombinedForce {
	CombinedForce(T1 f1, T2 f2, T3 f3 = NullForce()){
	F1 = f1;
	F2 = f2;
	F3 = f3;
	}

	template <typename NODE>
  	Point operator()(NODE n, double t) {
    	return F1(n,t) + F2(n,t) + F3(n,t);
	}
	
private:
T1 F1;
T2 F2;
T3 F3;
};

template <class T1, class T2, class T3>
CombinedForce<T1,T2,T3> make_combined_force(T1 f1, T2 f2, T3 f3) {
	return CombinedForce<T1,T2,T3>(f1, f2, f3);
	}

template <class T1, class T2>
CombinedForce<T1,T2> make_combined_force(T1 f1, T2 f2) {
	return CombinedForce<T1, T2>(f1, f2);
	}


struct PlaneConstraint {
template <typename Graph>
  	void operator()(Graph &graph, double t) {
    		for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
			if ((*it).position().z < -0.75) {
				(*it).position().z = -0.75;
				(*it).value().vel.z = 0;	
			}
	}
};

struct SphereConstraint {
template <typename Graph>
  	void operator()(Graph &graph, double t) {
    		for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
			if (norm((*it).position() - c) < r) {
				Point R = ((*it).position()-c)/norm((*it).position() - c);
				(*it).position() = c + R*r;
				(*it).value().vel = (*it).value().vel - dot((*it).value().vel, R)*R;
			}
		}
Point c = Point(0.5,0.5,-0.5);
double r = 0.15;
};

struct SphereConstraintDelete {
template <typename Graph>
  	void operator()(Graph &graph, double t) {
    		for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
			if (norm((*it).position() - c) < r) {
				graph.remove_node(it);
			}
		}
Point c = Point(0.5,0.5,-0.5);
double r = 0.15;
};

template <class T1, class T2 >
struct CombinedConstraint {
	CombinedConstraint(T1 f1, T2 f2){
	F1 = f1;
	F2 = f2;
	}

	template <typename Graph>
  	void operator()(Graph& graph, double t) {
    	F1(graph,t);
	F2(graph, t);
	}	
private:
T1 F1;
T2 F2;
};




template <class T1, class T2>
CombinedConstraint<T1,T2> make_combined_constraint(T1 f1, T2 f2) {
	return CombinedConstraint<T1, T2>(f1, f2);
	}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
	c = 1.0/graph.num_nodes();
	for (auto it = graph.node_begin(); it != graph.node_end(); ++it) 
		(*it).value().mass = 1.0 / graph.num_nodes();
	for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) 
		(*it).value().og_length = (*it).length();
	

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()), make_combined_constraint(PlaneConstraint(), SphereConstraintDelete()));
//make_combined_constraint(PlaneConstraint(), SphereConstraint())

        // Update viewer with nodes' new positions
       // Clear the viewer 's nodes and edges
	viewer.clear();
	node_map.clear();
	// Update viewer with nodes ' new positions and new edges
	viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
	viewer.add_edges (graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
