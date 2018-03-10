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

/** HW2 #2
 * Custom structure of data to store with Nodes */
struct EdgeData {
  double K;        //Spring constant, K_ij
  double L;        // Rest-length, L_ij
  EdgeData() : K(0), L(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


template < typename C1, typename C2 >
struct Constraints {
  C1 c1_;
  C2 c2_;
  Constraints(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}

  template < typename GRAPH >
  void operator()(GRAPH& g, double t) {
    c1_(g, t);
    c2_(g, t);
  }
};

template < typename C1, typename C2 >
Constraints < C1, C2 > make_combined_constraint ( C1 c1, C2 c2 ) {
 return Constraints < C1, C2 > (c1, c2);
}

template < typename C1, typename C2, typename C3 >
Constraints< C1, Constraints< C2, C3> > make_combined_constraint ( C1 c1, C2 c2, C3 c3 ) {
  return Constraints< C1, Constraints< C2, C3 > > (c1, c2, c3);
}

struct Constraint1 {
  template < typename GRAPH >
  void operator()( GRAPH& g, double t ) {
    (void) t;
    for( auto it = g.node_begin(); it != g.node_end(); ++it ) {
      auto n = *it;
      if( n.position()[2] < -0.75 ) {
        n.position() = n.position()*Point(1,1,0) + Point(0,0,-0.75);
        n.value().vel = n.value().vel*Point(1,1,0);
      }
    }
  }
};

struct Constraint2 {
  template < typename GRAPH >
  void operator()( GRAPH& g, double t ) {
    (void) t;
    for( auto it = g.node_begin(); it != g.node_end(); ++it ) {
      auto n = *it;
      Point r_dist = n.position() - Point(0.5, 0.5, -0.5);
      if( norm( r_dist ) < 0.15 ) {
        n.position() = Point(0.5, 0.5, -0.5) + 0.15 * (r_dist/norm(r_dist));
        n.value().vel = n.value().vel - dot(n.value().vel, r_dist/norm(r_dist))*(r_dist/norm(r_dist));
      }
    }
  }
};

struct Constraint3 {
  template < typename GRAPH >
  void operator()( GRAPH& g, double t ) {
    (void) t;
    for( auto it = g.node_begin(); it != g.node_end(); ++it ) {
      auto n = *it;
      Point r_dist = n.position() - Point(0.5, 0.5, -0.5);
      if( norm( r_dist ) < 0.15 ) {
        g.remove_node(n);
      }
    }
  }
};

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

// Combine constraints and define user constraints.
auto user_constraints = make_combined_constraint(Constraint1(), Constraint3());

//Constraint3 user_constraints;
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if( n.position() == Point(0) or n.position() == Point(1, 0, 0) )
      continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Impose the constraints before compute the force.
  user_constraints(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if( n.position() == Point(0) or n.position() == Point(1, 0, 0) )
      continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

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
    (void) t;
    // HW2 #1: YOUR CODE HERE
    /** Constrain the two corners of the cloth */
    if ( n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    
    /** Turn on for HW2 #1 */
    #if 0
    /** For the initial condition, when @a t = 0, 
     * Initial position:  Is set from the initial coordinates of the nodes file.
     * Therefore, no action required for initial position
     * Intial position of n = n.position()
     * K = 100.
     * L: initial length of the edges.
     * */
    double K = 100;
    Point f_spring = Point(0);
    //grid0
    //double L = 0.25;
    //grid1
    double L = 0.04166666667; 
    //grid2 
    //double L = 0.02040816327;
    //grid3
    //double L = 0.0101010101;
    #endif
    Point f_spring = Point(0);
    for ( auto ii = n.edge_begin(); ii != n.edge_end(); ++ii ) {
      auto adj_n = (*ii).node1();
      if ( adj_n == n )
        adj_n = (*ii).node2();
      f_spring += (-1 * (*ii).value().K ) * (( n.position() - adj_n.position() ) 
                  / norm( n.position() - adj_n.position() ))
                  * ( norm( n.position() - adj_n.position() ) - (*ii).value().L );
    }
    /** Add gravitaion force */
    return f_spring + ( Point(0, 0, -grav) * n.value().mass );
  }
};

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0, 0, -grav) * n.value().mass;
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point f_spring = Point(0);
    for ( auto ii = n.edge_begin(); ii != n.edge_end(); ++ii ) {
      auto adj_n = (*ii).node1();
      if ( adj_n == n )
        adj_n = (*ii).node2();
      f_spring += (-1 * (*ii).value().K ) * (( n.position() - adj_n.position() ) 
                  / norm( n.position() - adj_n.position() ))
                  * ( norm( n.position() - adj_n.position() ) - (*ii).value().L );
    }
    return f_spring;
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    /** for HW2 #3, c = 1/N, which is same as mass 
     * As such, avoid including extra value, c_ij, in NodeInfo */
    double c = n.value().mass; 
    return (-1)*(c * n.value().vel); 
  }
};
/** Return the sum of the two forces 
 * Inspired by std::pair, using template to inherit the type of the inputs.
 * @param[in] @a f1, @a f2 forces acting on the node 
 * @return @a f1 + @a f2, the summation of the two forces
 *
 * @tparam Force1, Force2 are structure that returns the force acting on certain
 *                        node( NODE n ) at time @a t. Take @a n and the @a t for the
 *                        input.
 * */
template < typename Force1, typename Force2 >
struct ExternalForces {
  Force1 f1_;
  Force2 f2_;
  ExternalForces(Force1 f1, Force2 f2) : f1_(f1), f2_(f2) {}

  template < typename NODE >
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t);
  }
};
/** Return the Point type of the two forces */
template < typename Force1, typename Force2 >
ExternalForces<Force1, Force2> make_combined_force(Force1 f1, Force2 f2) {
  return ExternalForces<Force1, Force2> (f1, f2);
}
/** Return the Point type of the three forces */
template < typename Force1, typename Force2, typename Force3 >
ExternalForces<Force1, ExternalForces<Force2, Force3> > 
              make_combined_force(Force1 f1, Force2 f2, Force3 f3) {
  return ExternalForces<Force1, ExternalForces<Force2, Force3> > (f1, make_combined_force(f2, f3));
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
  /**  
   * Zero initial velocity: set all nodes' velocity to be zero.
   * Mass: m_i = 1/N, where N = num_nodes
   * */
  for ( auto ni = graph.node_begin(); ni != graph.node_end(); ++ni ) {
    (*ni).value().vel = Point(0);
    (*ni).value().mass = float(1)/float(graph.num_nodes());
  }
  for ( auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei ) {
    (*ei).value().K = 100;
    (*ei).value().L = (*ei).length();
  }
    // Print out the stats
  std::cout << graph.size() << " " << graph.num_edges() << std::endl; //total number of the nodes and edges.
  //std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

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
        symp_euler_step(graph, t, dt, 
                        make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));
        // Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions and edges
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
