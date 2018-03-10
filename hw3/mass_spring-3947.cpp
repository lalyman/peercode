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

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K;       //< Spring constant
  double L;       //< Rest length
  EdgeData() : K(0), L(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

namespace mass_spring {
  GraphType::size_type const_node1_index {0};
  GraphType::size_type const_node2_index {0};
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
    //if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0))
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F>
void apply_constraints(G& g, double t, F constraints) {
  constraints(g, t);
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
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    else {
      double K {100.0};
      //double L {0.25};  // For grid0
      //double L {0.0416667}; // For grid1
      double L {0.0204082}; // For grid2
      //double L {0.010101};  // For grid3
      Point force {Point(0,0,0)};
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        force += -K*((*it).node1().position() - (*it).node2().position())*(1.0 - L/(*it).length());
      }
      return force + Point(0,0,-n.value().mass*grav);
    }
  }
};

struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    else {
      Point force {Point(0,0,0)};
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        double K {(*it).value().K};
        double L {(*it).value().L};
        force += -K*((*it).node1().position() - (*it).node2().position())*(1.0 - L/(*it).length());
      }
      return force + Point(0,0,-n.value().mass*grav);
    }
  }
};


struct GravityForce {
  /** Return the force of gravity applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    return Point(0,0,-n.value().mass*grav);
  }
};

struct MassSpringForce {
  /** Return the force of mass springs applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    Point force {Point(0,0,0)};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double K {(*it).value().K};
      double L {(*it).value().L};
      force += -K*((*it).node1().position() - (*it).node2().position())*(1.0 - L/(*it).length());
    }
    return force;
  }
};

struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    double c {1.0/ ((double) n.graph()->size())};   // Damping constant
    return -c*n.value().vel;
  }
};

/** A combination of two functors */
template <typename F1, typename F2>
struct function_pair {
  F1 f1;
  F2 f2;

  function_pair(const F1& _f1, const F2& _f2)
  : f1(_f1), f2(_f2) {}

  //template <typename TYPE>
  Point operator()(Node& n, double t) {
    return f1(n,t) + f2(n,t);
  }

  void operator()(GraphType& g, double t) {
    f1(g,t);
    f2(g,t);
  }
};

/** A combination of three functors */
template <typename F1, typename F2, typename F3>
struct function_trio {
  F1 f1;
  F2 f2;
  F3 f3;

  function_trio(const F1& _f1, const F2& _f2, const F3& _f3)
  : f1(_f1), f2(_f2), f3(_f3) {}

  //template <typename TYPE>
  Point operator()(Node& n, double t) {
    return f1(n,t) + f2(n,t) + f3(n,t);
  }

  void operator()(GraphType& g, double t) {
    f1(g,t);
    f2(g,t);
    f3(g,t);
  }
};

/** Combine two functors */
template <typename F1, typename F2>
function_pair <F1, F2> make_combined_function(F1 f1, F2 f2) {
  return function_pair <F1, F2> (f1, f2);
}

/** Combine three functors */
template <typename F1, typename F2, typename F3>
function_trio <F1, F2, F3> make_combined_function(F1 f1, F2 f2, F3 f3) {
  return function_trio <F1, F2, F3> (f1, f2, f3);
}

struct ConstantConstraint {
  /** Enforce constant constraints */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    g.node(mass_spring::const_node1_index).position() = Point(0,0,0);
    g.node(mass_spring::const_node2_index).position() = Point(1,0,0);
    g.node(mass_spring::const_node1_index).value().vel = Point(0,0,0);
    g.node(mass_spring::const_node2_index).value().vel = Point(0,0,0);
  }
};

struct PlaneConstraint {
  /** Find and fix nodes that violate plane constraint */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position().z < -0.75) {
        (*it).position().z = -0.75;
        (*it).value().vel.z = 0.0;
      }
    }
  }
};

struct SphereConstraint {
  /** Find and fix nodes that violate sphere constraint */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    Point c {Point(0.5,0.5,-0.5)};
    double r {0.15};
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      double L {norm((*it).position() - c)};
      if (L < r) {
        Point R {((*it).position() - c)/L};
        (*it).position() = c + R*r;
        (*it).value().vel -= dot((*it).value().vel,R)*R;
      }
    }
  }
};

struct SphereRemoveConstraint {
  /** Find and remove nodes that violate sphere constraint */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    Point c {Point(0.5,0.5,-0.5)};
    double r {0.15};
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      double L {norm((*it).position() - c)};
      if (L < r) {
        g.remove_node(*it);
        --it;
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
    (*it).value().mass = 1.0/((double) graph.size());
    if ((*it).position() == Point(0,0,0))
      mass_spring::const_node1_index = (*it).index();
    else if ((*it).position() == Point(1,0,0))
      mass_spring::const_node2_index = (*it).index();
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    (*it).value().K = 100.0;
    (*it).value().L = (*it).length();
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
        //symp_euler_step(graph, t, dt, Problem2Force());
        //symp_euler_step(graph, t, dt, make_combined_function(GravityForce(), MassSpringForce()));
        symp_euler_step(graph, t, dt, make_combined_function(GravityForce(), MassSpringForce(), DampingForce()));
        //apply_constraints(graph, t, make_combined_function(ConstantConstraint(), PlaneConstraint(), SphereConstraint()));
        apply_constraints(graph, t, make_combined_function(ConstantConstraint(), PlaneConstraint(), SphereRemoveConstraint()));

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
