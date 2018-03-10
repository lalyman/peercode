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

// Spring constant
static constexpr double K = 100;

// Original length of edges

//static double L = 0;

// Damping constant

static double c_damping = 0;


/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double rest_len;
  EdgeData() : rest_len(0) {}
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
    n.position() += n.value().vel * dt;
  }
 
  constraint(g, t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) n.value().vel = Point(0, 0, 0);
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
    // HW2 #1: YOUR CODE HERE
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) return Point(0, 0, 0);
    else{
      Point result(0, 0, -1*grav*n.value().mass);
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto nit = *it;
        result += (-1*K*((nit.length() - nit.value().rest_len)/nit.length())*(n.position() - nit.node2().position()));
      }
      return result;
    }
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    //return Point(0);
  }
};

/** Gravity force */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) return Point(0, 0, 0);
    //else{
      return Point(0, 0, -1*grav*n.value().mass);
    //}
  }
};

struct MassSpringForce {
  template<typename NODE> 
  Point operator()(NODE n, double t) { 
    //if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) return Point(0, 0, 0);
    //else{
      Point result(0, 0, 0);  
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto nit = *it;
        result += (-1*K*((nit.length() - nit.value().rest_len)/nit.length())*(n.position() - nit.node2().position()));
      }
      return result;
    //}
  }
};

struct DampingForce {
  template<typename NODE> 
  Point operator()(NODE n, double t) { 
    //if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) return Point(0, 0, 0);
    //else{
      return -1*c_damping*n.value().vel;
    //}
  }
};

/** helper structs for make_combined_force() */
template<typename FORCE1, typename FORCE2>
struct CombinedForce2{
  CombinedForce2(FORCE1 f10, FORCE2 f20): f1(f10), f2(f20){ } 
  
  template<typename NODE>
  Point operator()(NODE n, double t = 0){
    return f1(n, t) + f2(n, t); 
  }

  private:
  FORCE1 f1;
  FORCE2 f2;
};

template<typename FORCE1, typename FORCE2, typename FORCE3>
struct CombinedForce3{
  CombinedForce3(FORCE1 f10, FORCE2 f20, FORCE3 f30): f1(f10), f2(f20), f3(f30){ } 
  
  template<typename NODE>
  Point operator()(NODE n, double t = 0){
    return f1(n, t) + f2(n, t) + f3(n, t); 
  }

  private:
  FORCE1 f1;
  FORCE2 f2;
  FORCE3 f3;
};

/** combine 2 force */
template<typename FORCE1, typename FORCE2>
CombinedForce2<FORCE1, FORCE2> make_combined_force(FORCE1 f1, FORCE2 f2){
  return CombinedForce2<FORCE1, FORCE2>(f1, f2);
}

/** combine 3 force */
template<typename FORCE1, typename FORCE2, typename FORCE3>
CombinedForce3<FORCE1, FORCE2, FORCE3> make_combined_force(FORCE1 f1, FORCE2 f2, FORCE3 f3){
  return CombinedForce3<FORCE1, FORCE2, FORCE3>(f1, f2, f3);
}


/** constraints functor */
template<typename C1, typename C2>
struct ConstraintCombine{
  ConstraintCombine(C1 c10, C2 c20): c1(c10), c2(c20) { }

  template<typename GRAPH>
  void operator()(GRAPH& g, double t = 0){
    c1(g, t);
    c2(g, t);
  }

  private:
  C1 c1;
  C2 c2;
};

/** Plane constraint */
struct PlaneConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      auto n0 = *ni; 
      if(n0.position().z < -0.75){
        n0.position().z = -0.75;
        n0.value().vel.z = 0; 
      }
    }
  } 
};

/** Circle constraint */
struct CircleConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    Point center(0.5, 0.5, -0.5);
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      auto n0 = *ni;
      if(norm(n0.position() - center) < 0.15){
        Point r_norm = (n0.position() - center) / norm(n0.position() - center);
        n0.position() = center + 0.15 * r_norm;
        n0.value().vel = n0.value().vel - (n0.value().vel * r_norm)*r_norm;
      } 
    }
  }
};

/** Sphere constraint(remove node) */
struct BallConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    Point center(0.5, 0.5, -0.5);
    std::vector<typename GRAPH::Node> temp;
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      //auto n0 = *ni;
      //if(norm(n0.position() - center) < 0.15) temp.push_back(n0);
      //if(norm(n0.position() - center) < 0.15) g.remove_node(n0);
      while(ni != g.node_end() && norm((*ni).position() - center) < 0.15) g.remove_node(ni); 
    }
    //for(auto r_n : temp) g.remove_node(r_n);
  }
};

/** Combine constraint function*/
template<typename C1, typename C2>
ConstraintCombine<C1, C2> make_combined_constraint(C1 c10, C2 c20){
  return ConstraintCombine<C1, C2>(c10, c20);
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
#if 1
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for(auto itr_n = graph.node_begin(); itr_n != graph.node_end(); ++itr_n){
    (*itr_n).value().vel = Point(0, 0, 0);
    (*itr_n).value().mass = 1 / double(graph.num_nodes());
  }

  // Set initial condition for edges
  for(auto itr_e = graph.edge_begin(); itr_e != graph.edge_end(); ++itr_e){
    (*itr_e).value().rest_len = (*itr_e).length();
  }
  //auto it = graph.edge_begin();
  c_damping = (*graph.node_begin()).value().mass;

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();
#if 1
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
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()), make_combined_constraint(CircleConstraint(), PlaneConstraint()));
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()), BallConstraint());
  
        // Clear viewer
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        //if (graph.size() < 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        
        
      }

    });  // simulation thread
#endif
  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
