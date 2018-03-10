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

using size_type = unsigned;

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double K;       //< spring constant 
  double L;       //< rest-length
  EdgeData() : K(0), L(0) {}
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
  // std::cout << g.node(0).position() <<std::endl;
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    //if (n.position()!=Point(0 ,0 ,0) && n.position()!=Point(1 ,0 ,0)) {
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  	//} // end if
  }

  // Find violated constraints after the position update;
  // reset positions and velocities before the forces are calculated.
  // combined_constraints constraint3 = make_combined_constraints(constant_const(i1,i2),
  // 																	sphere_const(),plane_const());
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    //if (n.position()!=Point(0 ,0 ,0) && n.position()!=Point(1 ,0 ,0)) {
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    
    n.value().vel += force(n, t) * (dt / n.value().mass);
    //} //end if
  }

  return t + dt;
}

//----------------------------------------------------

struct GravityForce {
	template <typename NODE>
	Point operator()(NODE n, double ) {
		return n.value().mass*Point (0,0,-grav);
	}
};

struct MassSpringForce {
	template <typename NODE>
	Point operator()(NODE n, double ) {
		Point force;

		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
  		Edge e = *it;

  		Node u;
  		if (e.node1() == n) u = e.node2();
      else u = e.node1();

      double e_len = norm_2(n.position()-u.position());
  		force += -e.value().K*(n.position()-u.position())/e_len*(e_len-e.value().L);
  	}
	  return force;
	}
};

struct DampingForce {
	DampingForce(double dconst) 
		: dconst_(dconst){
	}

	template <typename NODE>
	Point operator()(NODE n, double ) {
		return -dconst_*n.value().vel;
	}

private: 
	double dconst_;
};

//----------------------------------------------------

template <typename F1, typename F2>
class combined_force {	
public: 
	combined_force(F1 f1, F2 f2)
		: f1_(f1), f2_(f2){
	}

	Point operator()(Node n, double t) {
		return f1_(n,t)+f2_(n,t);
	}
private:
	F1 f1_;
	F2 f2_;
};

template <typename F1, typename F2>
combined_force<F1,F2> make_combined_force(F1 f1, F2 f2){
	return combined_force<F1,F2>(f1,f2);
}

template <typename F1, typename F2, typename F3>
combined_force<combined_force<F1,F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3){
	return make_combined_force(make_combined_force(f1, f2), f3);
}

//---------PROBLEM 1 CODE-----------

/** Force function object for HW2 #1. */
//struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  /*template <typename NODE>
  Point operator()(NODE n, double t) {
  	// HW2 #1: YOUR CODE HERE
  	if (n.position()==Point(0 ,0 ,0) || n.position()==Point(1 ,0 ,0))
    	return Point (0 ,0 ,0);
    else {
	  	Point force = n.value().mass*Point (0,0,-grav);
	    
	  	for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
	  		Edge e = *it;
	  		Node u;
	  		if (e.node1() == n) u = e.node2();
	      else u = e.node1();

	      double e_len = norm_2(n.position()-u.position());
	  		
	  		force += -e.value().K*(n.position()-u.position())/e_len*(e_len-e.value().L);
	  	}
	  	return force;
  	}
  }
};*/

//----------------------------------------------------

struct constant_const {
	// void operator()(GraphType& n, double t) { 
	// 	for (auto it = n.node_begin(); it!=n.node_end(); ++it){
	// 		auto n = *it;			
	// 		if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
	// 			n.value().vel = Point(0,0,0);
	// 		}
	// 	}
	// }
	constant_const(size_type i1, size_type i2)
		:i1_(i1), i2_(i2) { 
	}

	void operator()(GraphType& g, double ) {
		g.node(i1_).position() = Point(0,0,0);
		g.node(i1_).value().vel = Point(0,0,0);
		g.node(i2_).position() = Point(1,0,0);
		g.node(i2_).value().vel = Point(0,0,0);
		
	}
	size_type i1_;
	size_type i2_;
};

struct sphere_const {
	void operator()(GraphType& n, double ) {
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
		for (auto it = n.node_begin(); it!=n.node_end(); ++it){
			auto n = *it;
			
			Point dist = n.position()-c;
			double dist_n = norm_2(dist);
			if (dist_n < r) {
				n.position() = c+ r*(dist/dist_n);
				double temp = inner_prod(n.value().vel,dist/dist_n);
				n.value().vel = n.value().vel-temp*(dist/dist_n);
			}
		}
	}
};

struct plane_const {
	void operator()(GraphType& n, double ) {
		for (auto it = n.node_begin(); it!=n.node_end(); ++it){
			auto n = *it;
			if (inner_prod(n.position(),Point(0,0,1))<-0.75) {
				n.position().elem[2] = -0.75;
				n.value().vel.elem[2] = 0;
			}
		}
	}
};

struct removal_const {
  void operator()(GraphType& n, double ) {
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    auto it = n.node_begin();
    while (it!=n.node_end()){
      auto node = *it;
      if(norm(node.position()-c) < r){
        //apply the version where we take in an iterator
        it = n.remove_node(it); 
      }
      else ++it;
    }
  }
};

//----------------------------------------------------


template <typename C1, typename C2>
class combined_constraints {	
public: 
	combined_constraints(C1 c1, C2 c2)
		: c1_(c1), c2_(c2){
	}

	void operator()(GraphType& g, double t) {
		c1_(g,t); c2_(g,t);
	}
private:
	C1 c1_;
	C2 c2_;
};

template <typename C1, typename C2>
combined_constraints<C1,C2> make_combined_constraints(C1 c1, C2 c2){
	return combined_constraints<C1,C2>(c1,c2);
}

template <typename C1, typename C2, typename C3>
combined_constraints<combined_constraints<C1,C2>,C3> make_combined_constraints(C1 c1, C2 c2, C3 c3){
	return make_combined_constraints(make_combined_constraints(c1, c2), c3);
}

//----------------------------------------------------


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
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
  	auto n = *it;
  	n.value().vel = Point(0,0,0);
  	n.value().mass = 1/double(graph.num_nodes());
  }

  for (auto it = graph.edge_begin();it != graph.edge_end(); ++it) {
  	auto e = *it;
  	e.value().K = 100;
  	e.value().L = e.length();
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

      size_type i1, i2;
	  	for (auto it = graph.node_begin(); it!=graph.node_end(); ++it){
	  		auto n = *it;
	  		if (n.position()==Point(0,0,0)) i1 = n.index();
	  		if (n.position()==Point(1,0,0)) i2 = n.index();
	  	}

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, 
        	make_combined_force(GravityForce(), 
        									MassSpringForce(),DampingForce(1/graph.num_nodes())),
        	make_combined_constraints(constant_const(i1,i2), removal_const(),
     															plane_const()));

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
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