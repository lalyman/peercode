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
#include <cmath>
#include <thread>
#include <iostream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,double>;
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  //Constrain the points if needed
  ConstrainPlane(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
struct MassSpringForce {
  // Return the force applying to @a n at time @a t.
  
  template <typename NODE>
  Point operator()(NODE n) {
    const double k = 100;
    auto pos = n.position();

    // Spring force start
    Point f = Point(0,0,0);

    for (auto a = n.edge_begin(); a != n.edge_end(); ++a){
      
      //Get the Edge and find the position
      Edge e = *a;
      double L = e.value();
      Point pos2;
      if (e.node2() == n){
        pos2 = e.node1().position();
      }
      else{
        pos2 = e.node2().position();
      }

      //Calculate the Spring force 
      auto sub = pos-pos2;
      auto dist = std::sqrt(norm_1(sub * sub));  
      f += -k * sub / dist * (dist-L);      
    }; 
    return f;
  }
};

struct GravityForce {
  // Return the force applying to @a n at time @a t.
  
  template <typename NODE>
  Point operator()(NODE n) {

    // set original force
    Point f = Point(0,0,0);

    for (auto a = n.edge_begin(); a != n.edge_end(); ++a){
      
      //Calculate the force from gravity
      f +=  n.value().mass * Point(0,0,-grav);      
    }; 
    return f;
  }
};


struct DampingForce {
  // Return the force applying to @a n at time @a t.
  
  template <typename NODE>
  Point operator()(NODE n) {

    // set original force
    Point f = Point(0,0,0);

    for (auto a = n.edge_begin(); a != n.edge_end(); ++a){
      
      //Calculate the force from gravity
      f +=  n.value().vel * c;      
    }; 
    return f;
  }
};


struct ConstrainPlane{

  template <typename G>
  void operator()(G& graph, double t){
     for (auto a = graph.node_begin(); a != graph.node_end(); a++){
      if ((*a).position().z < -0.75){
        (*a).position().z = -0.75;
        (*a).value().vel.z = 0;
      };
     };
  }
};


struct ConstrainSphere{

  template <typename G>
  void operator()(G& graph, double t){
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
     for (auto a = graph.node_begin(); a != graph.node_end(); a++){
      if (((*a).position()-c) < r){
        auto R = ((*a).position()-c)/std::abs((*a).position()-c);
        (*a).position() = c + r * R;
        (*a).value().vel -= dot((*a).value().vel,R)*R;
      };
     };
  }
};


struct emptyFunctor {

  template <typename NODE>
  Point operator()(NODE n) {
    return Point(0,0,0);
  }
};

// Return the combined 
struct make_combined_force {

  template <typename F>
  make_combined_force(F a){
    static const auto first = a;
    static const auto second = emptyFunctor();
    static const auto third = emptyFunctor();
  }
  
  template <typename F>
  make_combined_force(F a, F b){
    static const auto first = a;
    static const auto second = b;
    static const auto third = emptyFunctor();
  } 
  
  template <typename F>
  make_combined_force(F a, F b, F c){
    static const auto first = a;
    static const auto second = b;
    static const auto third = c;
  }

  template <typename NODE>
  Point operator()(NODE n){
    return first(n) + second(n) + third(n);
  }

};





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

  //return 0;

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  // Set initial mass and edge values
  c = 1/double(graph.num_nodes());
  for (auto x = graph.node_begin(); x != graph.node_end(); ++x){
    (*x).value().mass = c;
  }
  for (auto x = graph.edge_begin(); x != graph.edge_end(); ++x){
    (*x).value() = (*x).length();
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
      double dt = 0.0001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce()));

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
