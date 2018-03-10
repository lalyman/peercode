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
  double K;       //< Node velocity
  double L;     //< Node mass
  EdgeData() : K(1), L(1) {}
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
       n.value().vel = Point(0,0,0);
    }
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    else
        n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C cons) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;

  }

  cons(g);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
       n.value().vel = Point(0,0,0);
    }
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    else
        n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//double K;
//double L;
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
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        return Point(0,0,0);
    }
    Point force = Point(0,0,0);

    for( auto it = n.edge_begin(); it!=n.edge_end(); ++it){
      auto currEdge = *it;
      auto n2 = currEdge.node2();
      auto pos1 = n.position();
      auto pos2 = n2.position();
      double K = currEdge.value().K;
      double L = currEdge.value().L;
      force += -K*(pos1-pos2)/norm(pos1-pos2)*(norm(pos1-pos2)-L);


    }
    force += n.value().mass*Point(0,0,-grav);
    return force;
  }
};

struct MassSpringForce {
  /** Return the spring force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    Point force = Point(0,0,0);

    for( auto it = n.edge_begin(); it!=n.edge_end(); ++it){
      auto currEdge = *it;
      auto n2 = currEdge.node2();
      auto pos1 = n.position();
      auto pos2 = n2.position();
      double K = currEdge.value().K;
      double L = currEdge.value().L;
      force += -K*(pos1-pos2)/norm(pos1-pos2)*(norm(pos1-pos2)-L);
    }
    //force += n.value().mass*Point(0,0,-grav);
    return force;
  }
};

struct GravityForce {
  /** Return the gravity force acting on node @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    Point force = Point(0,0,-grav);
    return n.value().mass*force;
  }
};

struct DampingForce {
  /** Return the force damping dorce acting on node @a n at time @a t. */

  double c;

  DampingForce(const double coeff) : c(coeff) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;

    Point force = n.value().vel;
    return -c * force;
  }
};

/**
 * @brief The NoForce struct returns a zero force
 */
struct NoForce {
  /** Return zero force */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    (void) n;
    return Point(0,0,0);
  }
};

/**
 * @brief make_combined_force struct combines the two or three genralized forces passed as arguments
 */
template <typename force1, typename force2, typename force3>
struct make_combined_force {
    force1 f1;
    force2 f2;
    force3 f3;

    make_combined_force(force1 for1, force2 for2, force3 for3) : f1(for1), f2(for2), f3(for3) {}

    make_combined_force(force1 for1, force2 for2) : f1(for1), f2(for2), f3(NoForce()) {}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      // HW2 #1: YOUR CODE HERE
      (void) t;

      return f1(n,t) + f2(n,t) + f3(n,t);
    }
};

/**
 * @brief The FlatConstraint struct changes z velocity and position of nodes that violate constraint
 */

struct FlatConstraint {

    double limit = -0.75;
    template<typename G>
    void operator()(G& g){
        for(auto it = g.node_begin(); it!= g.node_end(); ++it){
            auto currNode = *it;
            auto pos = currNode.position();
            if (pos.z < limit){
                currNode.position().z = limit;
                currNode.value().vel.z = 0;
            }
        }
    }
};

/**
 * @brief The SphericalConstraint struct changes velocity and position of nodes that violate constraint
 */
struct SphericalConstraint {

    double r = 0.15;
    Point c = Point(0.5, 0.5, -0.5);
    template<typename G>
    void operator()(G& g){
        for(auto it = g.node_begin(); it!= g.node_end(); ++it){
            auto currNode = *it;
            auto pos = currNode.position();
            auto unitVec = (pos-c)/norm(pos-c);
            if (norm(pos - c) < r){
                auto unitVec = (pos-c)/norm(pos-c);
                currNode.position() = unitVec * r + c;
                currNode.value().vel -= inner_prod(currNode.value().vel, unitVec) * unitVec;
            }
        }
    }
};

/**
 * @brief The SphericalDestroyer struct removes nodes when they violate this constraint
 */
struct SphericalDestroyer {

    double r = 0.15;
    Point c = Point(0.5, 0.5, -0.5);
    template<typename G>
    void operator()(G& g){
        for(auto it = g.node_begin(); it!= g.node_end(); ++it){
            auto pos = (*it).position();
            if (norm(pos - c) < r){
                it = g.remove_node(it);
            }
        }
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
#if 1
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Setting initial conditions
  for (auto it = graph.node_begin(); it!=graph.node_end(); ++it){
      (*it).value().mass = 1/((double)graph.num_nodes());
      (*it).value().vel = Point(0,0,0);
  }

  for (auto it = graph.node_begin(); it!=graph.node_end(); ++it){
      for (auto in = (*it).edge_begin(); in!=(*it).edge_end(); ++in){
          (*in).value().K = 100;
          (*in).value().L = (*in).length();
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
        //std::cout << "t = " << t << std::endl;
        auto forces = make_combined_force<MassSpringForce, GravityForce, NoForce> (MassSpringForce(), GravityForce());
        //auto withDamping = make_combined_force<MassSpringForce, GravityForce, DampingForce> (MassSpringForce(), GravityForce(), DampingForce(1/((double)graph.num_nodes())));
        //symp_euler_step(graph, t, dt, Problem1Force());

        //symp_euler_step(graph, t, dt, forces);
        //symp_euler_step(graph, t, dt, withDamping);
        //symp_euler_step(graph, t, dt, withDamping, FlatConstraint());
        symp_euler_step(graph, t, dt, forces, SphericalDestroyer());
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
