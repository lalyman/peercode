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
#include <cmath>
#include <iostream>

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

struct  EdgeData {
  double K;
  double L;
  EdgeData(){}
  EdgeData(double K_,double L_): K(K_), L(L_){
  }
};


// Define the Graph type
// using GraphType = Graph<NodeData>;
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

    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      continue;
    }

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      continue;
    }
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  constraint(g, t);

  return t + dt;
}




template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      continue;
    }

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      continue;
    }
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
    // HW2 #1: YOUR CODE HERE
    // (void) n; (void) t; (void) grav;    // silence compiler warnings
    if(n.position() == Point(0.0, 0.0, 0.0) || n.position() == Point(1.0, 0.0, 0.0)){
      return Point(0.0, 0.0, 0.0);
    } else {
      Point spring_force = Point(0.0,0.0,0.0);
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto e = (*it);
        auto n2 = e.node2();
        Point diff = n.position() - n2.position();
        
        spring_force += -e.value().K * (diff) / norm(diff) * (norm(diff) - e.value().L);
      }
      return (spring_force + n.value().mass * Point(0, 0, -grav));
    }
  }
};



/** Function object for gravity force

*/
struct  GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return n.value().mass * Point(0, 0, -grav);
}
};


/** Function object for mass spring force.
*/
struct  MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
      Point spring_force = Point(0.0,0.0,0.0);
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto e = (*it);
        auto n2 = e.node2();
        Point diff = n.position() - n2.position();
        spring_force += -e.value().K * (diff) / norm(diff) * (norm(diff) - e.value().L);
      }
      return spring_force;
    }
};


/** Function object of damping force
*/
struct  DampingForce{
  double c_;
  DampingForce(double c):c_(c){
  };
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return -c_ * n.value().vel;
  }
};

/** Zero force object
*/
struct ZeroForce{
  template<typename NODE>
  Point operator()(NODE n, double t){
    return Point(0, 0, 0);
  }
};

/** Function object of combined force
*/
template<typename F1, typename F2, typename F3>
struct make_combined_force
{
  F1 f1_;
  F2 f2_;
  F3 f3_;
  // combine three forces constructor
  make_combined_force(F1 f1, F2 f2, F3 f3):f1_(f1),\
  f2_(f2), f3_(f3){
  }
  // combine two forces constructor
  make_combined_force(F1 f1, F2 f2): f1_(f1), f2_(f2),\
  f3_(F3()){
  }
  // return combined force
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1_(n, t) + f2_(n, t) + f3_(n, t);
  }
  
};

/** plane Constraint function object.
*/
struct planeConstraint{
  double z_;
  planeConstraint(double z):z_(z){
  }
  template<typename GRAPH>
  void operator()(GRAPH &g, double t){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      // std::cout << "here" << std::endl;
      auto n = *it;
      if(dot(n.position(), Point(0, 0, 1)) < z_){
        n.position().z = z_;
        n.value().vel.z = 0.0;
      }
    }
  }
};

/** sphere constraint function object.
*/
struct sphereConstraint
{
  Point c_;
  double r_;
  sphereConstraint(Point c, double r):c_(c), r_(r){
  }
  template<typename GRAPH>
  void operator()(GRAPH &g, double t){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(norm(n.position() - c_) < r_){
        Point R = (n.position() - c_) / norm(n.position() - c_);
        n.position() = c_ + R * r_;
        n.value().vel = n.value().vel - dot(n.value().vel, R) * R;
      }
    }
  }  
};

/** sphere constraint for question 4.
*/
struct sphereHoleConstraint{
  Point c_;
  double r_;
  sphereHoleConstraint(Point c, double r):c_(c), r_(r){
  }
  template<typename GRAPH>
  void operator()(GRAPH &g, double t){
    auto it = g.node_begin();
    while(it != g.node_end()){
      auto n = *it;
      if(norm(n.position() - c_) < r_){
        g.remove_node(n);
      } else{
        ++it;
      }
    }
  } 
};

/** Zero constrain function object
*/
struct ZeroConstraint
{
  template<typename GRAPH>
  void operator()(GRAPH &g, double t){
  }
};


/** combined constraint function object
* can combine two or three constraints
*/
template<typename C1, typename C2, typename C3>
struct make_combined_constraint
{
  C1 c1_;
  C2 c2_;
  C3 c3_;
  make_combined_constraint(C1 c1, C2 c2, C3 c3):c1_(c1),\
  c2_(c2), c3_(c3){
  }
  make_combined_constraint(C1 c1, C2 c2):c1_(c1),\
  c2_(c2), c3_(C3()){
  }  
  template<typename GRAPH>
  void operator()(GRAPH &g, double t){
    c1_(g,t);
    c2_(g,t);
    c3_(g,t);
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
  double spring_constant = 100.0;
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
    auto node = *it;
    // node.value().vel = Point(0, 0, 0);
    node.value().mass = 1.0 / graph.num_nodes();
  }
  for(auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
    auto e = *it;
    e.value().K = spring_constant;
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
      double t_end = 5;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        //*************** p1, p2, UNCOMMENT this part to see p1, p2 results
        // symp_euler_step(graph, t, dt, Problem1Force());
        // // Update viewer with nodes' new positions
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);
        //*************** p1, p2


        //*********** p3, generalize forces, UNCOMMENT this part to see p3 result
        //// combine two forces
        // auto  myforce = make_combined_force<GravityForce, MassSpringForce, ZeroForce>\
        // (GravityForce(), MassSpringForce());
        // //// combine three forces
        // // auto  myforce = make_combined_force<GravityForce, MassSpringForce, DampingForce>\
        // // (GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));
        // // Constraint is set to be zero.
        // symp_euler_step(graph, t, dt, myforce, ZeroConstraint());
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);
        //************** p3
    


        // p4, generalizing constraints, UNCOMMENT to see p4 result
        //  also did combine constraints. Used combined force.

        // ****** comment to see plane constraint
        // auto myforce = make_combined_force<GravityForce, MassSpringForce, DampingForce>\
        // (GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));

        // double z = -0.75;
        // symp_euler_step(graph, t, dt, myforce, planeConstraint(z));
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);
        // ************** 

        //// sphere constraint
        // // *********** uncomment to see sphere constraint.
        // auto myforce = make_combined_force<GravityForce, MassSpringForce, DampingForce>\
        // (GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));       
        // Point c = Point(0.5, 0.5, -0.5);
        // double r = 0.15;
        // symp_euler_step(graph, t, dt, myforce, sphereConstraint(c, r));
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);
        // // ******

        // combine constraint
        // ******* UNCOMMENT to see combined constraint.
        double z = -0.75;
        // Point c = Point(0.5, 0.5, -0.5);
        // double r = 0.15;

        // auto myforce = make_combined_force<GravityForce, MassSpringForce, DampingForce>\
        // (GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));
        // auto my_constraint = make_combined_constraint<planeConstraint,\
        //  sphereConstraint, ZeroConstraint>\
        // (planeConstraint(z), sphereConstraint(c,r));
        // symp_euler_step(graph, t, dt, myforce, my_constraint);

        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);
        // *******************************

        // ************** p4



        //********* test problem 5, remove nodes, UNCOMMENT to see p5 result

        auto myforce = make_combined_force<GravityForce, MassSpringForce, DampingForce>\
        (GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));

        // auto myforce = make_combined_force<GravityForce, MassSpringForce, ZeroForce>\
        // (GravityForce(), MassSpringForce());       


        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;
        auto my_constraint = make_combined_constraint<planeConstraint,\
         sphereHoleConstraint, ZeroConstraint>(planeConstraint(z), sphereHoleConstraint(c, r));
        symp_euler_step(graph, t, dt, myforce, my_constraint);

                // symp_euler_step(graph, t, dt, myforce, sphereHoleConstraint(c, r));

        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);
        //*******************

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
