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

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double len;
  double K_ = 100.0;
  EdgeData() {}
  EdgeData(double k, double l) {
    K_ = k;
    len = l;
  }
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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
	
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt

		// if the position of the node is not equal to the two points, update the position of the node
		if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
			n.position() += n.value().vel * dt;
		}
    
		
  }
	
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
		
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
		if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
			n.value().vel += force(n, t) * (dt / n.value().mass);
		}
    
  }
	constraint(g, t);
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point springForce = Point(0, 0, 0);
   
		auto curr = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
  		springForce += -(*it).value().K_ * (curr - (*it).node2().position()) * (norm(curr - (*it).node2().position()) - (*it).value().len) / norm(curr - (*it).node2().position());
		}

    Point gravityForce = n.value().mass * Point(0, 0, -grav);
    Point totalForce = gravityForce + springForce;
		(void) t;
    return totalForce;
  }
};

/** Gravity Force object */
struct GravityForce {
	template <typename NODE>
	Point operator()(NODE n, double t) {
		Point gravityForce = n.value().mass * Point(0, 0, -grav);
		(void) t;
		return gravityForce;
	}
};

/** MassSpring Force object */
struct MassSpringForce {
	template <typename NODE>
	Point operator()(NODE n, double t) {
		Point springForce = Point(0, 0, 0); 
		auto curr = n.position();
 		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
  		springForce += -(*it).value().K_ * (curr - (*it).node2().position()) * (norm(curr - (*it).node2().position()) - (*it).value().len) / norm(curr - (*it).node2().position());
		}
		(void) t;
  	return springForce;
	}
};

/** Damping Force object */
struct DampingForce {
	DampingForce() {
	} 
	DampingForce(double c) {
		c_ = c;	
	}
	
	template <typename NODE>
	Point operator()(NODE n, double t) {
		(void) t;	
		return -c_ * n.value().vel * 1.0;		
	}
	double c_ = 0.0;
};

/** Combined Force function object to combine two kinds of forces */
template <typename F1, typename F2>
struct CombinedForce {
	CombinedForce(F1 f1, F2 f2):f1_(f1), f2_(f2){} 
	template <typename NODE>
	Point operator()(NODE n, double t) {
		(void) t;		
		return f1_(n, t) + f2_(n, t);
	}
	F1 f1_;
	F2 f2_;
};

/** Function that returns the combined of two forces */
template<typename F1, typename F2> 
CombinedForce<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2 = F2()) {
	return CombinedForce<F1, F2>(f1, f2);
}

/** Function that returns the combined of three forces */
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 f1 = F1(), F2 f2 = F2(), F3 f3 = F3()) {
	return CombinedForce<CombinedForce<F1, F2>, F3>(make_combined_force(f1, f2), f3);
}

/* Constant Constraint function object */
struct ConstantConstraint {
	ConstantConstraint() {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
			if ((*itr).position() != Point(0, 0, 0) && (*itr).position() != Point(1, 0, 0)) {
				continue;
			}
			(*itr).value().vel = Point(0, 0, 0);
		}
		(void) t;
	}

};

/* Plane Constraint function object */
struct PlaneConstraint {
	PlaneConstraint(double z): z_(z) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
			if (dot((*itr).position(), Point(0, 0, 1)) >= -0.75) {
				continue;
			}
			(*itr).value().vel.elem[2] = 0;
			(*itr).position().elem[2] = z_;
		}
		(void)t;
	}
	double z_;
};

/* Sphere Constraint function object */
struct SphereConstraint {
	SphereConstraint(Point c, double r):center_(c), radius_(r) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
			if (norm((*itr).position() - center_) >= radius_) {
				continue;
			}
			Point R = ((*itr).position() - center_) / norm((*itr).position() - center_);
			(*itr).value().vel -= (dot((*itr).value().vel, R)) * R;
			(*itr).position() = R * radius_ + center_;
		}
		(void)t;
	}
	
	Point center_;
	double radius_;	
};

/* Another Sphere Constraint function object to remove nodes */
struct SphereConstraint2 {
	SphereConstraint2(Point c, double r):center_(c), radius_(r) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		for (auto itr = g.node_begin(); itr != g.node_end();) {
			if (norm((*itr).position() - center_) >= radius_) {
				++itr;
			} else {
				g.remove_node(*itr);
			}
		}
		(void)t;
	}
	
	Point center_;
	double radius_;	
};



/** Function object to combine two constraints */
template<typename C1, typename C2>
struct CombinedConstraint {
	CombinedConstraint(C1 c1, C2 c2):c1_(c1), c2_(c2) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		c1_(g, t);
		c2_(g, t);
		(void) t;
	}
	C1 c1_;
	C2 c2_;
};

/** Function that returns the combined constraint of two constraints */
template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
	return CombinedConstraint<C1, C2>(c1, c2);
}

/** Function that returns the combined constraint of three constraints */
template<typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
	return CombinedConstraint<CombinedConstraint<C1, C2>, C3>(make_combined_constraint(c1, c2), c3);
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
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  int numOfNodes = graph.num_nodes();
  for (auto nodeItr = graph.node_begin(); nodeItr != graph.node_end(); ++nodeItr) {
    (*nodeItr).value().vel = Point(0, 0, 0);
    (*nodeItr).value().mass = (1.0 / numOfNodes);
  }
  for (auto edgeItr = graph.edge_begin(); edgeItr != graph.edge_end(); ++edgeItr) {
    (*edgeItr).value().K_ = 100.0;
    (*edgeItr).dir().value().K_ = 100.0;
    
    (*edgeItr).value().len = (*edgeItr).length();
		(*edgeItr).dir().value().len = (*edgeItr).length();
		
  }

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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
				
        // symp_euler_step(graph, t, dt, Problem1Force());
				
				// customized force
				auto myForce = make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0 / graph.num_nodes()));
				// customized constraint
			 //	auto myConstraint = make_combined_constraint(ConstantConstraint(), PlaneConstraint(-0.75), SphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
				auto myConstraint2 = make_combined_constraint(ConstantConstraint(), PlaneConstraint(-0.75), SphereConstraint2(Point(0.5, 0.5, -0.5), 0.15));
				
				symp_euler_step(graph, t, dt, myForce, myConstraint2);
				
				viewer.clear();
				node_map.clear();
        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
				viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
				
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
