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

struct EdgeData {
  double spring_const;
  double rest_len;
  EdgeData() : spring_const(100.0), rest_len(1.0) {}
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

/** Force function object for HW2 #1. */

struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double ) {
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
    	return Point(0,0,0);
    }

    //const int K = 100;
    Point springF(0,0,0);
    Point gravF = n.value().mass * Point(0,0,-grav);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){
    	Edge current = *i;
    	Point diff = current.node1().position() - current.node2().position();
    	double dist = current.length();
        //springF += -K * diff / dist * (dist - 0.25);
        springF += (-current.value().spring_const * diff 
                   / dist * (dist - current.value().rest_len));
    }
    return springF + gravF;
  }
}; 

/** Function object for gravity force */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double ) {
  	return n.value().mass * Point(0,0,-grav);
  }
};

/** Function object for mass spring force */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double ) {
  	Point springF(0,0,0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){
    	Edge current = *i;
    	Point diff = current.node1().position() - current.node2().position();
    	double dist = current.length();
        springF += (-current.value().spring_const * diff 
                   / dist * (dist - current.value().rest_len));
    }
    return springF;
  }
};

/** Function object for damping force */
struct DampingForce {
  double degree;
  DampingForce(double d): degree(d) {};
  template <typename NODE>
  Point operator()(NODE n, double ) {
  	return -(1.0/degree)*n.value().vel; 
  }
};

template <typename F1, typename F2>
struct combineF{
  F1 force1;
  F2 force2;
  combineF(F1 f1_ = F1(), F2 f2_ = F2()): force1(f1_), force2(f2_){}

  template <typename NODE>	
  Point operator()(NODE n, double t){
    return force1(n,t)+force2(n,t);
  }
};

template <typename F1, typename F2>
combineF<F1,F2> make_combined_force(F1 f1_ = F1(), F2 f2_ = F2()){
  return combineF<F1,F2>(f1_, f2_);
}

template <typename F1, typename F2, typename F3>
combineF<combineF<F1,F2>, F3> make_combined_force(F1 f1_, F2 f2_, F3 f3_){
  return make_combined_force(make_combined_force(f1_,f2_),f3_);
}

/** Function object for const constraint */
struct ConstConstraint {
  void operator()(GraphType& g, double ){ 	
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
    	if ((*i).position()==Point(0,0,0) || (*i).position()==Point(1,0,0)){
    	  (*i).value().vel = Point(0,0,0);
        }
    }
  } 
};

/** Function object for plane constraint */
struct PlaneConstraint {
	void operator()(GraphType& g, double ){
      for (auto i = g.node_begin(); i != g.node_end(); ++i){
        Node current = *i;
        if (inner_prod(current.position(), Point(0,0,1)) < -0.75){
          current.position().z = -0.75;
          current.value().vel.z = 0.0;
        }
      }
	}
};

/** Function object for sphere constraint */
struct SphereConstraint {
  void operator()(GraphType& g, double ){
  	double r = 0.15;
    Point center(0.5,0.5,-0.5);
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
    	Node current = *i;
    	if (norm(current.position()-center) < r){
    	  /** HW2 #4
           current.position() = center + ((current.position() - center) / 
    	                         norm(current.position() - center) * r);
    	   Point R = ((current.position()-center) /
    	             norm(current.position() - center));
           current.value().vel = (current.value().vel - 
                                 inner_prod(current.value().vel, R) * R);
          */
    	  g.remove_node(current);
    	}
    }
  }
};

template <typename C1, typename C2>
struct combineC{
  C1 constraint1;
  C2 constraint2;
  combineC(C1 c1_ = C1(), C2 c2_ = C2()): constraint1(c1_), constraint2(c2_){}

  void operator()(GraphType& g, double t){
    constraint1(g,t);
    constraint2(g,t);
  }
};

template <typename C1, typename C2>
combineC<C1,C2> make_combinedC(C1 c1_ = C1(), C2 c2_ = C2()){
  return combineC<C1,C2>(c1_, c2_);
}

template <typename C1, typename C2, typename C3>
combineC<combineC<C1,C2>, C3> make_combinedC(C1 c1_, C2 c2_, C3 c3_){
  return make_combinedC(make_combinedC(c1_,c2_),c3_);
}

template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    
    /** HW2: #3
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
    	n.value().vel = Point(0,0,0);
    }
    else{
    	n.position() += n.value().vel * dt;
    } 
    */
  }
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  
  
  auto constraint_one = make_combinedC(ConstConstraint(), PlaneConstraint());
  auto constraint_two = make_combinedC(constraint_one, SphereConstraint());
  constraint_two(g,t); 
  
  return t + dt;
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
  // Set initial conditions for your nodes, if necessary.
  for (auto i = graph.node_begin(); i != graph.node_end(); ++i){
  	Node current = *i;
  	//current.position() = Point(0,0,0);
  	current.value().vel = Point(0,0,0);
  	current.value().mass = 1.0 / graph.num_nodes();
  }

  // HW2 #2: Set initial conditions for edges
  for (auto i = graph.edge_begin(); i != graph.edge_end(); ++i){
  	Edge edge_ = *i;
  	edge_.value().spring_const = 100.0;
    edge_.value().rest_len = edge_.length();
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
        //symp_euler_step(graph, t, dt, Problem1Force());
        
        symp_euler_step(graph, t, dt, 
        	make_combined_force(GravityForce(), 
        	                    DampingForce(graph.num_nodes()), 
        	                    MassSpringForce())); 

        // Clear the viewer’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes’ new positions and new edges
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
