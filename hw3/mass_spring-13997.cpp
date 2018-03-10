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
  double K, L;
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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    //n.position() = n.position() + n.value().vel * dt;
  }

  // Apply constraints
  constraints(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    //if (n.position() != Point(0) && n.position() != Point(1,0,0)) {
      n.value().vel += force(n, t) * (dt / n.value().mass);
    //}

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
    (void) t;

    // Create force to be returned.
    Point force_tot = Point(0);

    // Check if fixed point.
    Point pos = n.position();
    if (pos == Point(0) || pos == Point(1,0,0)) {
      return force_tot;
    }

    // Calculate mass spring force force_ms.
    Point force_ms = Point(0);
    double dist_norm(0.0);
    for (auto iedge_=n.edge_begin(); iedge_!=n.edge_end(); ++iedge_) {
      Point pos_j;
      if ((*iedge_).node1() == n) {
        pos_j     = (*iedge_).node2().position();
      } else {
        pos_j     = (*iedge_).node1().position();
      }
      dist_norm = norm(pos-pos_j);
      double K((*iedge_).value().K), L((*iedge_).value().L);
      force_ms += -K * (pos-pos_j) / dist_norm * (dist_norm - L);
    }

    // Calculate gravity force force_g.
    Point force_g = Point(0, 0, -grav*n.value().mass);

    // Add the forces and return the total force.
    force_tot += force_ms + force_g;
    return force_tot;
  }
};


/** Force function object for the gravity. */
struct GravityForce {
  /** Return the force applying to @a n at time @a t. */
  GravityForce() {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    // Point pos = n.position();

    /*// Check if fixed point. If so, return zero force.
    if (pos == Point(0) || pos == Point(1,0,0)) {
      return Point(0);
    }*/

    // Else return gravity force.
    return Point(0, 0, -grav*n.value().mass);
  }
};


/** Force function object for the mass spring. */
struct MassSpringForce {
  /** Return the force applying to @a n at time @a t. */
  MassSpringForce() {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    Point pos = n.position();

    /*// Check if fixed point. If so, return zero force.
    if (pos == Point(0) || pos == Point(1,0,0)) {
      return Point(0);
    }*/

    // Else calculate and return mass spring force.
    Point force_ms = Point(0);
    double dist_norm(0.0);
    for (auto iedge_=n.edge_begin(); iedge_!=n.edge_end(); ++iedge_) {
      Point pos_j;
      if ((*iedge_).node1() == n) {
        pos_j     = (*iedge_).node2().position();
      } else {
        pos_j     = (*iedge_).node1().position();
      }
      dist_norm = norm(pos-pos_j);
      double K((*iedge_).value().K), L((*iedge_).value().L);
      force_ms += -K * (pos-pos_j) / dist_norm * (dist_norm - L);
    }
    return force_ms;
  }
};


/** Force function object for the damping. */
struct DampingForce {
  /** Return the force applying to @a n at time @a t. */
  DampingForce() {}
  DampingForce(unsigned N) : c(1.0/N) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    // Point pos = n.position();

    /*// Check if fixed point. If so, return zero force.
    if (pos == Point(0) || pos == Point(1,0,0)) {
      return Point(0);
    }*/

    // Else return damping force.
    // double c(1.0/25);
    return -c*n.value().vel;
  }

  // Member variable.
  double c;
};


struct make_combined_force {

  // Contructors.
  make_combined_force()
    : damp_(false),         grav_(false),         ms_(false) {}
  make_combined_force(DampingForce f_damp)
    : f_damp_(f_damp),
      damp_(true),          grav_(false),         ms_(false) {}
  make_combined_force(GravityForce f_grav)
    :                       f_grav_(f_grav),
      damp_(false),         grav_(true),          ms_(false) {}
  make_combined_force(MassSpringForce f_ms)
    :                                             f_ms_(f_ms),
      damp_(false),         grav_(false),         ms_(true) {}
  make_combined_force(GravityForce f_grav, MassSpringForce f_ms)
    :                       f_grav_(f_grav),      f_ms_(f_ms),
      damp_(false),         grav_(true),          ms_(true) {}
  make_combined_force(DampingForce f_damp, MassSpringForce f_ms)
    : f_damp_(f_damp),                            f_ms_(f_ms),
      damp_(true),          grav_(false),         ms_(true) {}
  make_combined_force(DampingForce f_damp, GravityForce f_grav)
    : f_damp_(f_damp),      f_grav_(f_grav),
      damp_(true),          grav_(true),          ms_(false) {}
  make_combined_force(DampingForce f_damp, GravityForce f_grav, MassSpringForce f_ms)
    : f_damp_(f_damp),      f_grav_(f_grav),      f_ms_(f_ms),
      damp_(true),          grav_(true),          ms_(true) {}

  make_combined_force(MassSpringForce f_ms, GravityForce f_grav)
    : make_combined_force(f_grav, f_ms) {}
  make_combined_force(MassSpringForce f_ms, DampingForce f_damp)
    : make_combined_force(f_damp, f_ms) {}
  make_combined_force(GravityForce f_grav, DampingForce f_damp)
    : make_combined_force(f_damp, f_grav) {}
  make_combined_force(DampingForce f_damp, MassSpringForce f_ms, GravityForce f_grav)
    : make_combined_force(f_damp, f_grav, f_ms) {}
  make_combined_force(GravityForce f_grav, DampingForce f_damp, MassSpringForce f_ms)
    : make_combined_force(f_damp, f_grav, f_ms) {}
  make_combined_force(GravityForce f_grav, MassSpringForce f_ms, DampingForce f_damp)
    : make_combined_force(f_damp, f_grav, f_ms) {}
  make_combined_force(MassSpringForce f_ms, DampingForce f_damp, GravityForce f_grav)
    : make_combined_force(f_damp, f_grav, f_ms) {}
  make_combined_force(MassSpringForce f_ms, GravityForce f_grav, DampingForce f_damp)
    : make_combined_force(f_damp, f_grav, f_ms) {}


  // Returning the resulting force.
  Point operator()(Node n, double t) {
    Point f_return = Point(0);
    if (damp_) f_return += f_damp_(n, t);
    if (grav_) f_return += f_grav_(n, t);
    if (ms_)   f_return += f_ms_(n, t);
    return f_return;
  }

  // Member variables.
  DampingForce f_damp_;
  GravityForce f_grav_;
  MassSpringForce f_ms_;
  bool damp_, grav_, ms_;
};


struct constraint_plane {
  constraint_plane() : zmin(-0.75) {}

  // Apply constraint.
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;

    for (auto in_=g.node_begin(); in_!=g.node_end(); ++in_) {
      auto n_ = *in_;
      if (n_.position().z < zmin) {
        n_.position().z = zmin;
        n_.value().vel.z = 0.0;
      }
    }
  }

  // Member variable.
  double zmin;
};


struct constraint_sphere {
  constraint_sphere()
    : radius(0.15), center(Point(0.5, 0.5, -0.5)) {}

  // Apply constraint.
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;

    for (auto in_=g.node_begin(); in_!=g.node_end(); ++in_) {
      auto n_ = *in_;
      Point delta = n_.position() - center;
      if (norm(delta) < radius) {
        if (norm(delta) == 0) {
          n_.position().z += radius;
          n_.value().vel.z = 0.0;
        } else {
          n_.position() = center + delta / norm(delta) * radius;
          n_.value().vel -= delta / norm(delta) * inner_prod(n_.value().vel, delta/norm(delta));
        }
      }
    }
  }

  // Member variables.
  double radius;
  Point center;
};


struct constraint_sphRmv {
  constraint_sphRmv()
    : radius(0.15), center(Point(0.5, 0.5, -0.5)) {}


  // Apply constraint.
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;

    for (auto in_=g.node_begin(); in_!=g.node_end(); ++in_) {
      auto n_ = *in_;
      Point delta = n_.position() - center;
      if (norm(delta) < radius) {
        g.remove_node(n_);
      }
    }
  }

  // Member variables.
  double radius;
  Point center;
};


struct make_combined_constraint {

  // Contructors.
  make_combined_constraint()
    : srm_(false),    sph_(false),   pln_(false)    {}
  make_combined_constraint(constraint_sphRmv c_srm)
    : c_srm_(c_srm),
      srm_(true),    sph_(false),   pln_(false)    {}
  make_combined_constraint(constraint_sphere c_sph)
    :                c_sph_(c_sph),
      srm_(false),   sph_(true),    pln_(false)    {}
  make_combined_constraint(constraint_plane c_pln)
    :                               c_pln_(c_pln),
      srm_(false),   sph_(false),   pln_(true)     {}
  make_combined_constraint(constraint_sphere c_sph, constraint_plane c_pln)
    :                c_sph_(c_sph), c_pln_(c_pln),
      srm_(false),   sph_(true),    pln_(true)     {}
  make_combined_constraint(constraint_sphRmv c_srm, constraint_plane c_pln)
    : c_srm_(c_srm),                c_pln_(c_pln),
      srm_(true),    sph_(false),   pln_(true)     {}
  make_combined_constraint(constraint_sphRmv c_srm, constraint_sphere c_sph)
    : c_srm_(c_srm), c_sph_(c_sph),
      srm_(true),    sph_(true),    pln_(false)    {}
  make_combined_constraint(constraint_sphRmv c_srm, constraint_sphere c_sph, constraint_plane c_pln)
    : c_srm_(c_srm), c_sph_(c_sph), c_pln_(c_pln),
      srm_(true),    sph_(true),    pln_(true)     {}

  make_combined_constraint(constraint_plane c_pln, constraint_sphere c_sph)
    : make_combined_constraint(c_sph, c_pln) {}
  make_combined_constraint(constraint_plane c_pln, constraint_sphRmv c_srm)
    : make_combined_constraint(c_srm, c_pln) {}
  make_combined_constraint(constraint_sphere c_sph, constraint_sphRmv c_srm)
    : make_combined_constraint(c_srm, c_sph) {}
  make_combined_constraint(constraint_sphRmv c_srm, constraint_plane c_pln, constraint_sphere c_sph)
    : make_combined_constraint(c_srm, c_sph, c_pln) {}
  make_combined_constraint(constraint_plane c_pln, constraint_sphRmv c_srm, constraint_sphere c_sph)
    : make_combined_constraint(c_srm, c_sph, c_pln) {}
  make_combined_constraint(constraint_plane c_pln, constraint_sphere c_sph, constraint_sphRmv c_srm)
    : make_combined_constraint(c_srm, c_sph, c_pln) {}
  make_combined_constraint(constraint_sphere c_sph, constraint_sphRmv c_srm, constraint_plane c_pln)
    : make_combined_constraint(c_srm, c_sph, c_pln) {}
  make_combined_constraint(constraint_sphere c_sph, constraint_plane c_pln, constraint_sphRmv c_srm)
    : make_combined_constraint(c_srm, c_sph, c_pln) {}

  // Apply constraints.
  void operator()(GraphType& g, double t) {
    if (srm_) c_srm_(g, t);
    if (sph_) c_sph_(g, t);
    if (pln_) c_pln_(g, t);
    return;
  }

  // Member variables.
  constraint_sphRmv c_srm_;
  constraint_sphere c_sph_;
  constraint_plane  c_pln_;
  bool srm_, sph_, pln_;
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  double m = 1.0/nodes.size();
  for (auto inode_=graph.node_begin(); inode_!=graph.node_end(); ++inode_) {
    auto node_ = *inode_;
    node_.value().vel = Point(0);
    node_.value().mass = m;
    //node_.position()
  }

  for (auto iedge_=graph.edge_begin(); iedge_!=graph.edge_end(); ++iedge_) {
    auto edge_ = *iedge_;
    edge_.value().K = 100.0;
    edge_.value().L = edge_.length();
    // edge_.value().L = norm(edge_.node1().position()-edge_.node1().position());
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

        symp_euler_step(graph, t, dt,
          make_combined_force(DampingForce(graph.num_nodes()), GravityForce(), MassSpringForce()),
          make_combined_constraint(constraint_sphRmv(), constraint_plane()));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()));
        //symp_euler_step(graph, t, dt, make_combined_force(DampingForce(), GravityForce(), MassSpringForce()));
        //symp_euler_step(graph, t, dt, Problem1Force());

        // Update viewer with nodes' new positions
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
