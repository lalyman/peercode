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
#include <tuple>

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
  double K;
  double L;
  EdgeData() : K(1), L(0) {}
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
double symp_euler_step(G& g, double t, double dt, F force, C c) {
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

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  c(g, t);

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
    // return Point(0);
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    else {
      Point f(0,0,-grav);
      f *= n.value().mass;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        auto n_adj = (*it).node2();
        double dist = (*it).length();
        f -= (n.position()-n_adj.position()) / dist * (dist - (*it).value().L) * (*it).value().K;
      }
      return f;
    }
  }
};

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0,0,-grav) * n.value().mass;
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point f(0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto n_adj = (*it).node2();
      double dist = (*it).length();
      f -= (n.position()-n_adj.position()) / dist * (dist - (*it).value().L) * (*it).value().K;
    }
    return f;
  }
};

struct DampingForce {
  double c_;
  DampingForce(double c) : c_(c) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -n.value().vel * c_;
  }
};

template <typename F1, typename F2>
struct TwoForces {
  F1 force1;
  F2 force2;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return force1(n,t)+force2(n,t);
  }
};

template <typename F1, typename F2>
TwoForces<F1, F2> make_combined_force(F1 force1, F2 force2) {
  return TwoForces<F1, F2> {force1, force2};
};

template <typename F1, typename F2, typename F3>
TwoForces<TwoForces<F1,F2>,F3> make_combined_force(F1 force1, F2 force2, F3 force3) {
  return TwoForces<TwoForces<F1,F2>,F3> {make_combined_force(force1, force2), force3};
}

struct ConstantConstraint {
  Point p1_;
  Point p2_;
  ConstantConstraint(Point p1, Point p2) : p1_(p1), p2_(p2) {}
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position() == p1_ || (*it).position() == p2_)
        (*it).value().vel = Point(0);
    }
  }
};

struct PlaneConstraint {
  double z_;
  Point e_;
  PlaneConstraint(double z, Point e) : z_(z), e_(e) {}
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (dot((*it).position(), e_) < z_) {
        (*it).position().z = z_;
        (*it).value().vel.z = 0.0;
      }
    }
  }
};

struct SphereConstraint {
  Point c_;
  double r_;
  SphereConstraint(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      double dist = norm_2((*it).position()-c_);
      if (dist < r_) {
        Point R = ((*it).position()-c_) / dist;
        (*it).position() = c_ + R * r_;
        (*it).value().vel -= R * dot((*it).value().vel, R);
      }
    }
  }
};

struct SphereConstraintRemoveNode {
  Point c_;
  double r_;
  SphereConstraintRemoveNode(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ) {
      double dist = norm_2((*it).position()-c_);
      if (dist < r_)
        it = g.remove_node(it);
      else
        ++it;
    }
  }
};

template <typename C1, typename C2>
struct TwoConstraint{
  C1 c1;
  C2 c2;
  void operator()(GraphType& g, double t) {
    c1(g,t);
    c2(g,t);
  }
};

template <typename C1, typename C2>
TwoConstraint<C1,C2> make_combined_constraint(C1 c1, C2 c2) {
  return TwoConstraint<C1,C2> {c1, c2};
}

template <typename C1, typename C2, typename C3>
TwoConstraint<TwoConstraint<C1,C2>,C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return TwoConstraint<TwoConstraint<C1,C2>,C3> {make_combined_constraint(c1, c2),c3};
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
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().mass = 1.0 / graph.num_nodes();
    (*it).value().vel = Point(0,0,0);
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
      // double dt = 0.0005; // data/grid3
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());

        // symp_euler_step(graph, t, dt, \
          make_combined_force(GravityForce(), MassSpringForce()), \
          ConstantConstraint(Point(0,0,0), Point(1,0,0)));

        // symp_euler_step(graph, t, dt, \
          make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())), \
          make_combined_constraint(ConstantConstraint(Point(0,0,0), Point(1,0,0)), PlaneConstraint(-0.75, Point(0,0,1)), SphereConstraint(Point(0.5,0.5,-0.5), 0.15)));

        symp_euler_step(graph, t, dt, \
          make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())), \
          make_combined_constraint(ConstantConstraint(Point(0,0,0), Point(1,0,0)), PlaneConstraint(-0.75, Point(0,0,1)), SphereConstraintRemoveNode(Point(0.5,0.5,-0.5), 0.15)));

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
