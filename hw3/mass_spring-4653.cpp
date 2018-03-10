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
  EdgeData(): spring_const(100.0), rest_len(1.0){}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;





/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double) {
    // HW2 #1: YOUR CODE HERE

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }

    Point springF = Point(0,0,0);
    Point gravity = n.value().mass * Point(0,0,-grav);

    // use incident iterator to go through adjacent vertices
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      auto current = *iter;
      Point vec_diff = current.node1().position() - current.node2().position();
      double current_len = current.length();
      double K  = current.value().spring_const; // spring constant
      double rest_length = current.value().rest_len;
      springF += -K * vec_diff / current_len * (current_len - rest_length);
    }
    return springF + gravity;
  }
};

/** function object for gravity **/
struct Gravity {
  template <typename NODE>
  Point operator()(NODE n, double) {
    return n.value().mass * Point(0,0,-grav);
  }
};

/** fuction object for mass spring force **/
struct SpringForce{
  template <typename NODE>
  Point operator()(NODE n, double){
    Point springF(0,0,0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      auto current = *iter;
      Point vec_diff = current.node1().position() - current.node2().position();
      double current_len = current.length();
      double K  = current.value().spring_const; // spring constant
      double rest_length = current.value().rest_len;
      springF += -K * vec_diff / current_len * (current_len - rest_length);
    }
    return springF;
  }
};

/** function object for damping **/
struct DampingForce {
  double deg;
  DampingForce(double d): deg(d){};
  template <typename NODE>
  Point operator()(NODE n, double){
    return -1.0/deg * n.value().vel;
  }
};

/** function object for combining forces **/
template <typename F1, typename F2>
struct combineForces{
  F1 f1;
  F2 f2;
  combineForces(F1 f1_=F1(), F2 f2_=F2()): f1(f1_), f2(f2_){}
  template <typename NODE>
  Point operator() (NODE n, double t){
    return f1(n,t) + f2(n,t);
  }
};

/** function for combining two forces **/
template <typename F1, typename F2>
combineForces<F1,F2> make_combined_force(F1 f1_=F1(), F2 f2_=F2()){
  return combineForces<F1,F2>(f1_,f2_);
}

/** function for combining three forces **/
template <typename F1, typename F2, typename F3>
combineForces<combineForces<F1,F2>,F3> make_combined_force(F1 f1_=F1(), F2 f2_=F2(), F3 f3_=F3()){
  return make_combined_force(make_combined_force(f1_,f2_),f3_);
}


/** function object for the constant constraint **/
struct constConstraint {
  void operator()(GraphType& g, double){
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter){
      auto current = *iter;
      if (current.position() == Point(0,0,0) || current.position() == Point(1,0,0)){
        current.value().vel = Point(0,0,0);
      }
    }
  }
};

/** function object for the plane constraint **/
struct planeConstraint {
  void operator()(GraphType& g, double){
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter){
      auto current = *iter;
      auto angle = inner_prod(current.position(), Point(0,0,1));
      if (angle < -0.75){
        current.position().z = -0.75;
        current.value().vel.z=0.0;
      }
    }
  }
};

/** function object for the sphere constraint **/
struct sphereConstraint {
  void operator()(GraphType& g, double){
    double radius = 0.15;
    Point center(0.5,0.5,-0.5);
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter){
      Node current = *iter;
      auto diff_vec = current.position()-center;
      if (norm(diff_vec)<radius){
        Point R_pt = diff_vec / norm(diff_vec);
        current.position() = center+ R_pt * radius;
        current.value().vel=(current.value().vel-inner_prod(current.value().vel,R_pt)*R_pt);
        g.remove_node(current);
      }
    }
  }
};

/** function object for combining constrants **/
template <typename C1, typename C2>
struct combineConstraints{
  C1 c1;
  C2 c2;
  combineConstraints(C1 c1_=C1(), C2 c2_=C2()): c1(c1_), c2(c2_){}
  void operator() (GraphType& g, double t){
    c1(g,t);
    c2(g,t);
  }
};

/** function for combining two constraints **/
template <typename C1, typename C2>
combineConstraints<C1,C2> make_combined_constr(C1 c1_=C1(), C2 c2_=C2()){
  return combineConstraints<C1,C2>(c1_,c2_);
}

/** function for combining three forces **/
template <typename C1, typename C2, typename C3>
combineConstraints<combineConstraints<C1,C2>,C3> make_combined_constr(C1 c1_=C1(), C2 c2_=C2(), C3 c3_=C3()){
  return make_combined_constr(make_combined_constr(c1_,c2_),c3_);
}



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

    // code for problem 3
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
    //   n.value().vel = Point(0,0,0);
    // }
    // else {
    //   n.position() += n.value().vel * dt;
    // }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  auto constraint_one = make_combined_constr(constConstraint(),planeConstraint());
  auto total_constraint = make_combined_constr(constraint_one, sphereConstraint());

  total_constraint(g,t);
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
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter){
    Node current = *iter;
    current.value().vel = Point(0,0,0);
    current.value().mass  = 1.0 / graph.num_nodes();
  }

  for (auto iter = graph.edge_begin(); iter != graph.edge_end(); ++iter){
    Edge current = *iter;
    current.value().spring_const =100.0;
    current.value().rest_len  = current.length();
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
        symp_euler_step(graph, t, dt, make_combined_force(
          Gravity(), SpringForce(), DampingForce(graph.num_nodes()))
        );


        // Clear the viewer's nodes and EDGES
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
