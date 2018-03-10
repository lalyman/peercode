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
#include <iostream>
#include <stdlib.h>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

using size_type = std::size_t;

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double K = 100.0; // spring constant
double c; // damping constant
Point Zero_Vector = Point(0,0,0);

void set_damping_const(size_type i){
   c = 1.0/((double) i);
}


/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double L;
  double K;
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/* Remove sphere constraint */
struct RemoveSphere {
  Point c = Point(0.5,0.5,-0.5); double r = 0.15;
  template<typename G>
  void operator()(G& g){

    for(auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it; auto p = n.position();
      if(norm(p - c) < r){
        it = g.remove_node(it);
      }
    }
  }
};

/* plain constraint functor */
struct zBarrier {
   template <typename G>
   void operator()(G& g){
     for (auto it = g.node_begin(); it != g.node_end(); ++it) {
       auto n = *it;
       if (n.position().z < -0.75){
          n.position().z = -0.74; // nearest point
          n.value().vel.z = 0.0; // reset velocity in z direction
       }
     }
     return;
   }
};

/** Sphere block functor */
struct SphereBlock {
  Point c = Point(0.5,0.5,-0.5); double r = 0.15;
  template<typename G>
  void operator()(G& g){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point x = n.position();
      auto R_i = (x - c)/norm(x - c);
      if(norm(x - c) < r){
        n.value().vel -= (n.value().vel * R_i) * R_i;
        n.position() = c + r * R_i;
      }
    }
    return;
  }
};

/* fix node constraint */
struct FixNode {
   template<typename N>
   void operator()(N& n){
     n.value().vel = Zero_Vector;
     return;
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

template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  auto i = 0;
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() = n.position() + n.value().vel * dt;
    i++;
  }

  /* fixing node constraint */
  FixNode fix_node;

  /* uncomment to test ball barrier */
  //SphereBlock zz;
  //zz(g);

  /* uncomment to test floor effect */
  //zBarrier zz;
  //zz(g);

  RemoveSphere yy;
  yy(g);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
       fix_node(n); // apply constraint fix node
    } else {
       // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
       n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }
  return t + dt;
}


template <typename G>
void initialize (G& g){
  for(auto it=g.node_begin(); it!=g.node_end();++it){
    // initialize mass
    (*it).value().mass = (1.0/((double)g.num_nodes()));
    // initialize velocity vector
    (*it).value().vel = Point(0,0,0);
  }
  for(auto it=g.edge_begin(); it!=g.edge_end();++it){
    // initialize length
    (*it).value().L = (*it).length();
    // initialize spring constant
    (*it).value().K = K;
  }
}

/** Gravity Force function object for HW2 #3. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings
    return  n.value().mass*Point(0,0,-1*grav);
  }
};

/** Mass Spring Force function object for HW2 #3. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings

    Point force_s, x_i;
    force_s = Point(0,0,0);
    x_i = n.position();

    // spring force contribution
    for (auto it=n.edge_begin(); it!=n.edge_end();++it){
       Point x_j;
       if (n.index() == (*it).node1().index()){
          x_j = (*it).node2().position();
       } else {
          x_j = (*it).node1().position();
       }
       Point temp1 = x_i - x_j;

       double temp2 = norm(x_i - x_j);
       double L = (*it).value().L;
       double K = (*it).value().K;

       force_s += (-K*(temp2 - L)  / temp2) *(temp1) ;
    }
    return force_s;  
  }
};

/** Damping Force function object for HW2 #3. */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings
    return -1 * c * n.value().vel;
  }
};

template <typename force1, typename force2, typename force3 >
struct make_c3{
   force1 f1; force2 f2; force3 f3;
   make_c3 ( force1 ff1, force2 ff2, force3 ff3 );
   template <typename NODE>
   Point operator()(NODE n, double t) {
     return f1(n, t) + f2(n, t) + f3(n,t);
   }
};
template <typename force1, typename force2, typename force3>
make_c3<force1,force2,force3>::make_c3(force1 ff1,force2 ff2,force3 ff3){
   f1 = ff1; f2 = ff2; f3 = ff3;
}

template <typename force1, typename force2>
struct make_c2{
   force1 f1; force2 f2;
   make_c2 ( force1 ff1, force2 ff2 );
   template <typename NODE>
   Point operator()(NODE n, double t) {
     return f1(n, t) + f2(n, t);
   } 
};
template <typename force1, typename force2>
make_c2<force1,force2>::make_c2(force1 ff1,force2 ff2){
   f1 = ff1; f2 = ff2;
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
    (void) t;    // silence compiler warnings
    Point force_s, force_g, x_i;

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
       return Point(0,0,0);

    force_s = Point(0,0,0);
    force_g = Point(0,0,0);
    x_i = n.position();

    // gravity force contribution 
    force_g = n.value().mass*Point(0,0,-1*grav);
    // spring force contribution
    for (auto it=n.edge_begin(); it!=n.edge_end();++it){
       Point x_j;
       if (n.index() == (*it).node1().index())
          x_j = (*it).node2().position();
       else 
          x_j = (*it).node1().position();
       Point temp1 = x_i - x_j;
       double temp2 = norm(x_i - x_j);
       double L = (*it).value().L;
       double K = (*it).value().K;
       force_s += (-K*(temp2 - L)  / temp2) *(temp1) ;
    }
    return force_s+force_g;
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
  initialize(graph);
  set_damping_const(graph.num_nodes());

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  //Problem1Force.set_num_nodes(graph.num_nodes());

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

        //symp_euler_step(graph, t, dt, Problem1Force());
                  
        // combine 2 forces case
/*           make_c2<GravityForce, MassSpringForce>
                     make_combined_force( (GravityForce())
                                        , (MassSpringForce()));
*/
        // combine 3 forces case
        make_c3<GravityForce, MassSpringForce, DampingForce> 
                     make_combined_force( (GravityForce())
                                        , (MassSpringForce())
                                        , (DampingForce())); 

        symp_euler_step(graph, t, dt, make_combined_force);

        // clearing viewer nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions and new edges
        viewer.add_nodes(graph.node_begin(),graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(),graph.edge_end(), node_map);

        // Update viewer with nodes' new positions
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
