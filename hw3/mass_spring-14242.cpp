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

// Define the Graph type
using GraphType = Graph<NodeData, double>;
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
 * @tparam G::node_value_type contains elements @a mass and @a vel
 *           which support +, -, *, / with doubles.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where @a t is the current time.
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

  // Reset nodes that violate constraints
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update velocity of the node according to the force
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
  double k; // Spring constant
  double l; // Spring rest-length

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    //return Point(0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    else{
      Point force = n.value().mass * Point(0,0,-grav);
      for(auto iit = n.edge_begin(); iit != n.edge_end(); ++iit){
        Point diff = n.position()-(*iit).node2().position();
        force += -k * diff * (norm(diff)-l) / norm(diff);
      }
      return force;
    }
  }
};

/** Force function object that represents gravitational force */
struct GravityForce{
  /**
   * Return the gravitational force applying to @a n at time @a t.
   * @tparam NODE Node object that supports value() to obtain internal
   *    values
   * @param n Node on which force acts
   * @param t Time index
   * @return Point representing amount of gravitational force
   *    acting on node @a n
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Force function object that represents spring force */
struct MassSpringForce{
  double k; // Spring constant

  /**
   * Return the spring force applying to @a n at time @a t.
   * @tparam NODE Node object that supports iterators edge_begin()
   *    and edge_end() as well as element position() of type Point
   * @param n Node on which force acts
   * @param t Time index
   * @return Point representing amount of spring force
   *    acting on node @a n
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force = Point(0);
    for(auto iit = n.edge_begin(); iit != n.edge_end(); ++iit){
      Point diff = n.position()-(*iit).node2().position();
      force += -k * diff * (norm(diff)-(*iit).value()) / norm(diff);
    }
    return force;
  }
};

/** Force function object that represents damping force */
struct DampingForce{
  double c; // Damping constant

  /**
   * Return the damping force applying to @a n at time @a t.
   * @tparam NODE Node object that supports value() to obtain internal
   *    values
   * @param n Node on which force acts
   * @param t Time index
   * @return Point representing amount of damping force
   *    acting on node @a n
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return -c*n.value().vel;
  }
};

/** Generic Forces object that combines to Force function objects */
template<class F1, class F2>
struct Forces{
  F1 force1;
  F2 force2;

  /**
   * Return the sum of force1 and force2 acting on node @a n
   *    at time @a t
   * @tparam NODE Node type
   * @param n Node on which forces act
   * @param t Time index
   * @return Point representing amount of combined force acting
   *    on node @a n at time @a t
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return force1(n,t) + force2(n,t);
  }
};

/**
 * Two argument version to make Forces object
 * @tparam F1 Function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam F2 Function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @param f1 First force to combine in Forces object
 * @param f2 Second force to combine in Forces object
 * @return
 */
template<class F1, class F2>
Forces<F1, F2> make_combined_force(F1 f1, F2 f2){
  return Forces<F1, F2>{f1, f2};
}

/**
 * Three argument version to make Forces object
 * @tparam F1 Function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam F2 Function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam F3 Function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @param f1 First force to combine in Forces object
 * @param f2 Second force to combine in Forces object
 * @param f3 Third force to combine in Forces object
 * @return
 */
template<class F1, class F2, class F3>
Forces<F3, Forces<F1, F2>> make_combined_force(F1 f1, F2 f2, F3 f3){
  return Forces<F3, Forces<F1, F2>>{f3, Forces<F1, F2>{f1, f2}};
}

/** Constraint function object that fixes a node's position */
struct ConstConstraint{
  Node n;     // Node to fix
  Point pos;  // Position of node

  /**
   * Reset position of Node n to pos
   * @tparam G Graph object
   * @param g Graph on which to reset node
   * @param t Time index
   */
  template <typename G>
  void operator()(G& g, double t) {
    n.position() = pos;
    n.value().vel = Point(0,0,0);
  }
};

/** Constraint function object that prevents points from falling below
 * the plane z = -0.75*/
struct PlaneConstraint{
  /**
   * Reset position of all nodes in @a g that fall below plane z = -0.75
   * to be on plane z = -0.75.
   * @tparam G Graph object that supports iterators Node iterators
   * node_begin() and node_end()
   * @param g Graph on which to apply constraint
   * @param t Time index
   */
  template <typename G>
  void operator()(G& g, double t) {
    for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
      double z = dot((*nit).position(), Point(0,0,1));
      if(z < -0.75){
        (*nit).position().z = -0.75;
        (*nit).value().vel.z = 0.0;
      }
    }
  }
};

/** Constraint function object that prevents points from falling through
 * the sphere with radius 0.15 and center (0.5, 0.5, -0.5) */
struct SphereConstraint{
  /**
   * Reset position of all nodes in @a g that fall into the sphere
   * with radius 0.15 and center (0.5, 0.5, -0.5)
   * @tparam G Graph object that supports iterators Node iterators
   * node_begin() and node_end()
   * @param g Graph on which to apply constraint
   * @param t Time index
   */
  template <typename G>
  void operator()(G& g, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;

    for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
      Point diff = (*nit).position() - center;
      if(norm(diff) < r){
        Point R = diff / norm(diff);
        // Set position to point nearest on surface of sphere
        (*nit).position() = R*r + center;

        // Set velocity component normal to sphere surface to zero
        (*nit).value().vel = (*nit).value().vel - R*dot((*nit).value().vel,R);
      }
    }
  }
};

/** Constraint function object that removes points that falling through
 * the sphere with radius 0.15 and center (0.5, 0.5, -0.5) */
struct TearConstraint{
  /**
   * Remove all nodes in @a g that fall into the sphere
   * with radius 0.15 and center (0.5, 0.5, -0.5)
   * @tparam G Graph object that supports iterators Node iterators
   * node_begin() and node_end()
   * @param g Graph on which to apply constraint
   * @param t Time index
   */
  template <typename G>
  void operator()(G& g, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;

    for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
      Point diff = (*nit).position() - center;
      if(norm(diff) < r){
        g.remove_node(nit);
      }
    }
  }
};

/** Constraint function object that combines two other constraint function
 * objects
 * @tparam C1 Constraint function object called as @a constraint(@a g, @a t),
 *           where g is the graph and @a t is the current time.
 * @tparam C2 Constraint function object called as @a constraint(@a g, @a t),
 *           where g is the graph and @a t is the current time.
 */
template<class C1, class C2>
struct Constraints{
  C1 constraint1;
  C2 constraint2;

  template <typename G>
  void operator()(G& g, double t) {
    constraint1(g,t);
    constraint2(g,t);
  }
};

/**
 * Combines two constraint function objects to return another constraint object
 * @tparam C1 Constraint function object called as @a constraint(@a g, @a t),
 *           where g is the graph and @a t is the current time.
 * @tparam C2 Constraint function object called as @a constraint(@a g, @a t),
 *           where g is the graph and @a t is the current time.
 * @param c1 First constraint to combine
 * @param c2 Second constraint to combine
 * @return Constraint function object that combines @a c1 and @a c2
 */
template<class C1, class C2>
Constraints<C1, C2> combine_constraints(C1 c1, C2 c2){
  return Constraints<C1, C2>{c1, c2};
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

  // Set node initial velocity to 0 and mass to 1/num_nodes()
  for(auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
    (*nit).value().vel = Point(0);
    (*nit).value().mass = 1 / (double) graph.num_nodes();
  }

  // Set edge initial rest lengths to initial lengths
  for(auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit){
    (*eit).value() = norm((*eit).node1().position()-(*eit).node2().position());
  }

  // Spring constant is 100
  double k = 100;

  // Damping constant is 1/num_nodes()
  double c = 1 / (double) graph.num_nodes();

  // Get nodes at position (0,0,0) and position (1,0,0)
  Node node0, node1;
  for(auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
    if((*nit).position() == Point(0,0,0)) {
      node0 = *nit;
    }else if((*nit).position() == Point(1,0,0)){
      node1 = *nit;
    }
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
        symp_euler_step(graph, t, dt,
                        make_combined_force(GravityForce{},
                                            MassSpringForce{k},
                                            DampingForce{c}),
                        combine_constraints(ConstConstraint{node0, Point(0,0,0)},
                                            combine_constraints(
                                                ConstConstraint{node1, Point(1,0,0)},
                                                TearConstraint{})));

        // Clear the viewer's nodes and edges
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
