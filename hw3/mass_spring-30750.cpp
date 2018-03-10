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
  double rest_len;  //< Rest length
  double spr_const; //< Spring constant
  EdgeData() : rest_len(0), spr_const(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


//
// SYMPLECTIC EULER METHOD
//

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

  // Apply constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


//
// GENERALIZED FORCES
//

/** Force function object. */
struct Problem1Force {
  /** Return the force applied to @a n at time @a t.
   *
   * This is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    // If node is at (0,0,0) or (1,0,0), apply zero force
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    // Compute spring force acting on the node
    Point f_cumm = Point(0,0,-n.value().mass*grav); // Start with force due to gravity
    // Add spring force from each adjacent edge
    for (auto iiter = n.edge_begin(); iiter != n.edge_end(); ++iiter) {
      auto e = *iiter;
      f_cumm += -e.value().spr_const*(e.length() - e.value().rest_len)*
                (e.node1().position() - e.node2().position())/e.length();
    }
    return f_cumm;
  }
};

/** No force function object */
struct ZeroForce {
  /** Return no force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; (void) n;
    return Point(0,0,0);
  }
};

/** Gravity force function object */
struct GravityForce {
  /** Return the gravity force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0,0,-n.value().mass*grav);
  }
};

/** Spring force function object */
struct MassSpringForce {
  /** Return the spring force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; 
    Point f_cumm = Point(0,0,0); // Start with zero force
    // Add spring force from each adjacent edge
    for (auto iiter = n.edge_begin(); iiter != n.edge_end(); ++iiter) {
      auto e = *iiter;
      f_cumm += -e.value().spr_const*(e.length() - e.value().rest_len)*
                (e.node1().position() - e.node2().position())/e.length();
    }
    return f_cumm;
  }
};

/** Damping force function object */
struct DampingForce {
  /** Return the damping force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -d_coeff_*n.value().vel;
  }
  /** Constructor */
  DampingForce(double d_coeff) : d_coeff_(d_coeff) {};
 private:
  double d_coeff_;
};

/** Method to combine force function objects. */
template <typename F1, typename F2, typename F3 = ZeroForce>
struct make_combined_force {
  /** Return the damping force applied to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }
  /** Constructors */
  make_combined_force(F1 f1, F2 f2) : f1_(f1), f2_(f2), f3_(ZeroForce()) {};
  make_combined_force(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {};
 private:
  F1 f1_;
  F2 f2_;
  F3 f3_;
};


//
// GENERALIZED CONSTRAINTS
//

/** Zero constraint function object */
struct ZeroConstraint {
  /** Apply no constraint to a graph @a g. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t; (void) g;
    return;
  }
};

/** Horizontal plane constraint function object */
struct HorzPlaneConstraint {
  /** Constrain a graph @a g to be above a horizontal plane at @a h. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto niter = g.node_begin(); niter != g.node_end(); ++niter) {
      auto n = *niter;
      // If the node is below the plane, set the position to the nearest
      // point on the plane and set the vertical velocity to zero.
      if (n.position().z < h_) {
        n.position().z = h_;
        n.value().vel.z = 0;
      }
    }
    return;
  }
  /** Constructor */
  HorzPlaneConstraint(double h) : h_(h) {};
 private:
  double h_;
};

/** Sphere constraint function object */
struct SphereConstraint {
  /** Constrain a graph @a g to be outside a sphere with radius @a r and 
   *  center @a c. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto niter = g.node_begin(); niter != g.node_end(); ++niter) {
      auto n = *niter;
      // If the node is inside the sphere, set the position to the nearest point 
      // on the surface and set the normal velocity to the surface to zero.
      if (norm_2(n.position() - c_) < r_) {
        Point R = (n.position() - c_)/norm_2(n.position() - c_);
        n.position() = c_ + r_*R;
        n.value().vel -= dot(n.value().vel,R)*R;
      }
    }
    return;
  }
  /** Constructor */
  SphereConstraint(const Point& c, double r) : c_(c), r_(r) {};
 private:
  const Point& c_;
  double r_;
};

/** Fixed constraint function object */
struct FixedConstraint {
  /** Constrain a graph @a g to be fixed at nodes @a n_fixed. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto niter = g_fixed_.node_begin(); niter != g_fixed_.node_end(); ++niter) {
      auto n = *niter;
      // For nodes in the graph of fixed nodes, maintain the same position
      g.node(n.value()).position() = n.position();
    }
    return;
  }
  /** Constructor */
  FixedConstraint(const Graph<int,int>& g_fixed) : g_fixed_(g_fixed) {};
 private:
  const Graph<int,int>& g_fixed_;
};

/** Method to combine constraint function objects. */
template <typename C1, typename C2, typename C3 = ZeroConstraint>
struct make_combined_constraint {
  /** Constrain the graph with the given constriants. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    c1_(g,t); c2_(g,t); c3_(g,t);
    return;
  }
  /** Constructors */
  make_combined_constraint(C1 c1, C2 c2) : c1_(c1), c2_(c2), c3_(ZeroConstraint()) {};
  make_combined_constraint(C1 c1, C2 c2, C3 c3) : c1_(c1), c2_(c2), c3_(c3) {};
 private:
  C1 c1_;
  C2 c2_;
  C3 c3_;
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

    // Diagonal edges
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for the nodes and edges
  for (auto niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
    auto node = *niter;
    node.value().vel = Point(0,0,0);
    node.value().mass = 1./graph.num_nodes();
  }
  for (auto eiter = graph.edge_begin(); eiter != graph.edge_end(); ++eiter) {
    auto edge = *eiter;
    edge.value().rest_len = edge.length();
    edge.value().spr_const = 100;
  }
  // Set fixed nodes in separate graph for FixedConstraint() functor
  std::vector<Point> fixed_nodes;
  fixed_nodes.push_back(Point(0,0,0));
  fixed_nodes.push_back(Point(1,0,0));
  Graph<int,int> g_fixed;
  for (auto niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
    auto node = *niter;
    if (std::find(fixed_nodes.begin(),fixed_nodes.end(),node.position()) 
        != fixed_nodes.end()) {
      g_fixed.add_node(node.position());
      g_fixed.node(g_fixed.num_nodes() - 1).value() = node.index();
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
        //symp_euler_step(graph, t, dt, Problem1Force());
        auto CombinedForce = make_combined_force<
                             GravityForce,MassSpringForce,DampingForce>(
                             GravityForce(),MassSpringForce(),
                             DampingForce(1./graph.num_nodes()));
        auto CombinedConstraints = make_combined_constraint<
                                   HorzPlaneConstraint,SphereConstraint,
                                   FixedConstraint>(
                                   HorzPlaneConstraint(-.75),
                                   SphereConstraint(Point(.5,.5,-.5),.15),
                                   FixedConstraint(g_fixed));
        symp_euler_step(graph, t, dt, CombinedForce, CombinedConstraints);
        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
