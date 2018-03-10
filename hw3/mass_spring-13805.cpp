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

/** Custom structure of data to store with Edges */
struct EdgeData {
  double spring_const;
  double rest_length;
  EdgeData() : spring_const(100), rest_length(1) {}
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
    // Enforce constraint that two corners are pinned down so blanket doesn't fall
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      continue;
    } else {
      n.position() += n.value().vel * dt;
    }
  }

  // Add additional constraints (PlaneConstraint, SphereConstraint, etc)
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/******************************************************************************
 ******************************    FORCES     *********************************
 *****************************************************************************/

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
    Point force_total = Point(0, 0, 0);
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return force_total; // these points don't move
    }
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      auto edge = (*ii);
      NODE neighbor = edge.node2();
      if (neighbor.index() == n.index()) {
        neighbor = edge.node1();
      }
      Point pos_diff = n.position() - neighbor.position();
      force_total += -edge.value().spring_const * (pos_diff) / norm(pos_diff) 
                    * (norm(pos_diff) - edge.value().rest_length);
    }
    force_total += Point(0, 0, -grav * n.value().mass);
    (void) t;
    return force_total;
  }
};


/** Gravity force for HW #2, Problem 3*/
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point gravity_force = Point(0, 0, -grav * n.value().mass);
    return gravity_force;
  }
};


/** Mass-spring force for HW #2, Problem 3*/
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force = Point(0, 0, 0);
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      auto edge = (*ii);
      NODE neighbor = edge.node2();
      if (neighbor.index() == n.index()) {
        neighbor = edge.node1();
      }
      Point pos_diff = n.position() - neighbor.position();
      spring_force += -edge.value().spring_const * (pos_diff) / norm(pos_diff) 
                    * (norm(pos_diff) - edge.value().rest_length);
    }
    return spring_force;
  }
};


/** Damping force for HW #2, Problem 3*/
struct DampingForce {
  double c_;
  DampingForce(double c) : c_(c) {};
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point damping_force = -c_ * n.value().vel;
    (void) t;
    return damping_force;
  }
};


/** Combine the different force functors */
template <typename F1, typename F2>
struct CombinedForce {
  F1 f1_;
  F2 f2_;
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t);
  }
};

/** Combine two force functors */
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1, F2>(f1, f2);
} 

/** Combine three force functors */
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return make_combined_force(make_combined_force(f1, f2), f3);
} 


/******************************************************************************
 ***************************     CONSTRAINTS     ******************************
 *****************************************************************************/

/** Constrains nodes from moving past a plane */
struct PlaneConstraint {
  double z_;
  PlaneConstraint(double z) : z_(z) {}
  template <typename G>
  void operator()(G& graph, double t) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      Node curr = (*ni);
      if (dot(curr.position(), Point(0, 0, 1)) < z_) {
        curr.position().z = z_;
        curr.value().vel.z = 0.0;
      }
    }
    (void) t;
  }
};

/** Constrains nodes from entering a sphere */
struct SphereConstraint {
  Point c_;
  double r_;
  SphereConstraint(Point c, double r) : c_(c), r_(r) {}
  template <typename G>
  void operator()(G& graph, double t) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      Node curr = (*ni);
      Point diff_vec = curr.position() - c_;
      double dist = norm(diff_vec);
      if (dist < r_) {
        curr.position() = c_ + diff_vec / dist * r_;
        curr.value().vel = curr.value().vel - 
               dot(curr.value().vel, diff_vec / dist) * diff_vec / dist;
      }
    }
    (void) t;
  }
};


/** Remove nodes and associated edges from entering a sphere */
struct SphereConstraint2 {
  Point c_;
  double r_;
  SphereConstraint2(Point c, double r) : c_(c), r_(r) {}
  template <typename G>
  void operator()(G& graph, double t) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      Node curr = (*ni);
      double dist = norm(curr.position() - c_);
      if (dist < r_) {
        graph.remove_node(curr);
      }
    }
    (void) t;
  }
};


/** Functor to combine two templated constraints */
template <typename C1, typename C2>
struct CombinedConstraint {
  C1 c1_;
  C2 c2_;
  CombinedConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  template <typename G>
  void operator()(G& g, double t) {
    c1_(g, t);
    c2_(g, t);
  }
};

/** Combines two constraint functors */
template <typename C1, typename C2>
CombinedConstraint<C1, C2> combine_constraints(C1 c1, C2 c2) {
  return CombinedConstraint<C1, C2>(c1, c2);
}

/** Combines three constraint functors */
template <typename C1, typename C2, typename C3>
CombinedConstraint<C1, CombinedConstraint<C2, C3>> combine_constraints(C1 c1, C2 c2, C3 c3) {
  return combine_constraints(c1, combine_constraints(c2, c3));
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
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    (*ni).value().mass = 1.0 / graph.num_nodes();
    (*ni).value().vel  = Point(0, 0, 0);
  }
  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    (*ei).value().spring_const = 100;
    (*ei).value().rest_length  = (*ei).length();
  }
  double c = 1.0 / graph.num_nodes();
  auto combined_force = make_combined_force(GravityForce(), MassSpringForce(), 
                                            DampingForce(c));

  // Add constraints
  double z {-0.75};
  Point center = Point(0.5, 0.5, -0.5);
  double r {0.15};
  auto plane_constraint = PlaneConstraint(z);
  auto sphere_constraint = SphereConstraint(center, r);
  auto sphere_constraint2 = SphereConstraint2(center, r);
  auto combined_constraint = combine_constraints(plane_constraint, sphere_constraint2);
  (void) sphere_constraint;

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
        symp_euler_step(graph, t, dt, combined_force, combined_constraint);
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
