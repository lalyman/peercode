/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdarg>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to stor
e with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */

struct EdgeData {
	double Kij;
	double Lij;

	EdgeData() :Kij(100), Lij(0.1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned int;


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

struct DefaultConstraint {
	template<typename GRAPH>
	void operator() (GRAPH&, double) {
	}
};

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

	
	/*Problem 3: Constraints handling*/
	/*
	if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
		continue;
	*/

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

	/*Problem 3: Constraints handling*/
	/*
	if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
		continue;
	*/

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
	return symp_euler_step(g, t, dt,force,DefaultConstraint());
}



/** Force function object for HW2 #1. */
//Notes, here if you call Problem1Force, the graph generated is not for Problem 1 because I built my Problem 2 based on this function
//So the results is actually for problem 2 and this Problem1Force does not take any argument input
//L is set to the edge distance between two nodes. K is always kept to 100.
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings

	 if (n.position()==Point(0,0,0)||n.position()==Point(1,0,0))
		return Point(0,0,0);

	 /*GravityForce*/
	 Point force=Point(0, 0, -grav);

	 /*SpringForce*/
	 Point p1 = n.position();
	 for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
		 auto curredge = (*it);
		 Point p2 = curredge.node2().position();
		 double d = norm(p1 - p2);
		 force = force - ((curredge.value().Kij) / d)*(p1 - p2)*(d - curredge.value().Lij);
	 }
	std::cout<<"F= "<<force<<std::endl;
	return force;
  }

};

/** Problem 3 **/

class GravityForce {
public:
	template <typename NODE>
	Point operator()(NODE , double) {
		Point force = Point(0, 0, -grav);
		return force;
	}
};

class MassSpringForce {
public:
	template <typename NODE>
	Point operator()(NODE n, double) {
		Point force = Point(0, 0, 0);
		Point p1 = n.position();
		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
			auto curredge = (*it);
			Point p2 = curredge.node2().position();
			double d = norm(p1 - p2);
			force = force - ((curredge.value().Kij) / d)*(p1 - p2)*(d - curredge.value().Lij);
		}
		return force;
	}
};

class DampingForce {
public:
	DampingForce(GraphType G){
		c = (double) 1 / (G.num_nodes());
	}

	template <typename NODE>
	Point operator()(NODE n, double) {
		return (-c)*n.value().vel;//proportional to the velocity of this node
	}
	
private:
	double c;
};

/** A class that combines all types of forces which 
* has overloading () operator takes arguments on NODE n and double t (time)
*
* Given two forces struct [@a f1, @a f2]
* this class combines the forces generated by these two structs together, which has overloading operator ()
* takes arguments on Node n and double t as well
*/

template <typename T1, typename T2>
class combined_force {
public:
	combined_force(T1 f1, T2 f2) : f1_(f1), f2_(f2) {}

	template <typename NODE>
	Point operator()(NODE n, double t) {
		Point force1 = f1_(n,t);
		Point force2 = f2_(n,t);

		return force1 + force2;
	}

private:
	T1 f1_;
	T2 f2_;
};

/** Helper function for constructing combined_forces. This deduces the type of
* the forces so the user doesn't have to write it.
*
* Usage:
* // Construct a force which has overloading operator () takes arguments of Node n and double t (time)
* auto force= make_combined_force(f1(),f2()); 
* where f1() and f2() are the constructors for these force structs
*/

template <typename T1, typename T2>
combined_force<T1, T2> make_combined_force(T1 f1, T2 f2) {
	return combined_force<T1, T2>(f1,f2);
}

/** Helper function for constructing combined_forces.
* Similar to the functions above, but takes three force structs as the inputs
* Usage:
* // Construct a force which has overloading operator () takes arguments of Node n and double t (time)
* auto force= make_combined_force(f1(),f2(),f3());
* where f1(), f2() and f3() are the constructors for these force structs
*/
template <typename T1, typename T2, typename T3>
combined_force<combined_force<T1,T2>, T3>make_combined_force(T1 f1, T2 f2, T3 f3) {
	auto a= make_combined_force(f1, f2);
	return combined_force<combined_force<T1,T2>, T3>(a, f3);
}

/*Constant Constraint apply only to specific nodes and keep these nodes fixed
* It takes the argument of the index of these specific nodes and keep track of them
*/
struct ConstantConstraint {
	ConstantConstraint(size_type i, size_type j):i_(i),j_(j) {}

	template <typename GRAPH>
	void operator()(GRAPH& g, double) {
		g.node(i_).position() = Point(0, 0, 0);
		g.node(i_).value().vel = Point(0, 0, 0);
		g.node(j_).position() = Point(1, 0, 0);
		g.node(j_).value().vel = Point(0, 0, 0);
	}

private:
	size_type i_;//index for Point(0,0,0)
	size_type j_;//index for Point(1,0,0)
};

/*Constraint applies to all nodes.If a node is not on one side of the plane, 
* we set its position to the nearest point on the plane
* and set the velocity normal to the direction of the plane to be 0
*/

struct PlaneConstraint {
	template <typename GRAPH>
	void operator()(GRAPH& g, double) {
		for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			auto n = *it;
			if (n.position().z < -0.75) {
				n.position().z = -0.75;
				n.value().vel.z = 0;
			}
		}
	}
};

/*Constraint applies to all nodes.If a node is not outside the sphere
* we set its position to the nearest point on the sphere
* and set the velocity normal to the surface of the sphere to be 0
*/
struct SphereConstraint {
	template <typename GRAPH>
	void operator()(GRAPH& g, double) {
		Point c = Point(0.5, 0.5, -0.5);
		for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			auto n = *it;
			Point R = n.position() - c;
			if (norm(R) < 0.15) {
				Point unitR = (1.0 / norm(R))*R;
				n.position() = c + 0.15*unitR;
				n.value().vel -= dot(n.value().vel, unitR)*unitR;
			}
			
		}
	}
};

/*Constraint applies to all nodes.If a node is not outside the sphere
* we remove this node and all its incident edges
*/
struct P5Constraint {
	template <typename GRAPH>
	void operator()(GRAPH& g, double) {
		Point c = Point(0.5, 0.5, -0.5);
		for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			auto n = *it;
			Point R = n.position() - c;
			if (norm(R) < 0.15) {
				g.remove_node(n);
			}

		}
	}
};

/** A class that combines all types of constraints which
* has overloading () operator takes arguments on GRAPH g and double t (time)
*
* Given two forces struct [@a c1, @a c2]
* this class combines the constraintsgenerated by these two structs together, which has overloading operator ()
* takes arguments on GRAPH g and double t as well
*/
template <typename T1, typename T2>
class combined_constraint {
public:
	combined_constraint(T1 c1, T2 c2) : c1_(c1), c2_(c2) {}

	template <typename GRAPH>
	void operator()(GRAPH& g, double t) {
		c1_(g,t);
		c2_(g,t);

	}

private:
	T1 c1_;
	T2 c2_;
};


/** Helper function for constructing combined_constraint. This deduces the type of
* the constraints so the user doesn't have to write it.

*
* Usage:
* // Construct a constraint which has overloading operator () takes arguments of GRAPH g and double t (time)
* auto force= make_combined_constraint(c1(),c2());
* where c1() and c2() are the constructors for these force structs
*/
template <typename T1, typename T2>
combined_constraint<T1, T2> make_combined_constraint(T1 c1, T2 c2) {
	return combined_constraint<T1, T2>(c1, c2);
}

/** Similar as above, but with three input constraints**/
template <typename T1, typename T2, typename T3>
combined_constraint<combined_constraint<T1, T2>, T3>make_combined_constraint(T1 c1, T2 c2, T3 c3) {
	auto a = make_combined_constraint(c1, c2);
	return combined_constraint<combined_constraint<T1, T2>, T3>(a, c3);
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
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary. And find the node index which position is (0,0,0) or (1,0,0)
  size_type id1;
  size_type id2;

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  auto n = *it;
	  n.value().vel = Point(0,0,0);
	  n.value().mass = (double)1 / graph.num_nodes();
	  if (n.position() == Point(0, 0, 0))
		  id1 = n.index();
	  if (n.position() == Point(1, 0, 0))
		  id2 = n.index();
  }

  //P3: Initialize L for each edge
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
	  auto edg = *it;
	  edg.value().Lij = edg.length();
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
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

	  
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
		//symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce(graph)));
		  //symp_euler_step(graph,t,dt,Problem1Force(L));//This function not valid any more
         //symp_euler_step(graph, t, dt, Problem1Force());
		 /* symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(graph)),
			  make_combined_constraint(ConstantConstraint(id1, id2), PlaneConstraint(),SphereConstraint()));
			  */
		 /*symp_euler_step(graph, t, dt, 
			              make_combined_force(GravityForce(), MassSpringForce(), DampingForce(graph)),
						  P5Constraint());
						  */
       symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(graph)),
	                    make_combined_constraint(ConstantConstraint(id1,id2),P5Constraint()));
		 //Problem 5
		 viewer.clear();
		 node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
		//Problem 5
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
