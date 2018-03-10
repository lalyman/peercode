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
#include <set>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
//Damping constant
static double cdamp;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store as Edge values */
struct EdgeData {
  double L;
  double K;
  EdgeData() : L(0.1), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
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

    if (!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))){
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if (!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))){
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/** Function object to represent a plane constraint. */
struct PlaneConstraint {
   /** Augment graph g to comply with plane constraint at time t.
   *
   * If a node positions's z coordinate is less than -0.75, we set it to
   * be -0.75 and set its velocity in the z-direction to be 0. */
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    for (auto nfirst = g.node_begin(); nfirst != g.node_end(); ++nfirst) {
      auto n = *nfirst;
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};

/** Function object to represent a sphere constraint. */
struct SphereConstraint {
  /** Augment graph g to comply with sphere constraint at time t.
   *
   * If a node positions's is inside the sphere, change its position
   * to the nearest position on the surface of the sphere and the component
   * of velocity normal to the sphere's surface 0. */
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    Point center (0.5, 0.5, -0.5);
    double radius = 0.15;
    for (auto nfirst = g.node_begin(); nfirst != g.node_end(); ++nfirst) {
      auto n = *nfirst;
      auto dist = n.position() - center;
      if (norm(dist) < radius) {
        auto ri = dist/norm(dist);
        n.position() = center + radius/norm(dist)*dist;
        n.value().vel = n.value().vel - dot(n.value().vel, ri)*ri;
      }
    }
  }
};

/** Function object to represent a destructive sphere constraint. */
struct SphereConstraintRemove {
  /** Augment graph g to comply with sphere constraint at time t.
   *
   * If a node positions's is inside the sphere, remove the node. */
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    Point center (0.5, 0.5, -0.5);
    double radius = 0.15;
    for (auto nfirst = g.node_begin(); nfirst != g.node_end(); ++nfirst) {
      auto n = *nfirst;
      if (norm(n.position() - center) < radius) {
        g.remove_node(n);
      }
    }
  }
};

/** Function object to represent no constraint. */
struct NoConstraint {
  /** Do not act on the graph. */
  template <typename G>
  void operator()(G& g, double t) {
    (void) g;
    (void) t;
  }
};

/** Function object to represent a combined constraint. */
template <typename C1, typename C2>
struct CombinedConstraint {
  C1 c1_;
  C2 c2_;

  /** Apply both constraints to graph g at time t. */
  template <typename G>
  void operator()(G& g, double t){
    c1_(g,t);
    c2_(g,t);
  }
};

/** Create combined constraint from two individual constraints. */
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1_=NoConstraint(), C2 c2_=NoConstraint()) {
  return CombinedConstraint<C1, C2>({c1_, c2_});
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
    (void) t;

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    auto fgravity = (n.value().mass)*Point(0,0,-grav);

    auto fspring = Point(0);

    auto fnet = Point(0);

    auto first = n.edge_begin();
    auto last = n.edge_end();
    for ( ; first != last; ++first) {
      auto inc_edge = *first;
      auto adj_node = inc_edge.node2();
      fspring = fspring - inc_edge.value().K*((n.position()-adj_node.position())/inc_edge.length())*(inc_edge.length() - inc_edge.value().L);
    }
    fnet = fspring + fgravity;
    return fnet;
  }
};

/** Function object to represent Gravity force. */
struct GravityForce {
  /** Return force valued at -mass*gravity for node n at time t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (n.value().mass)*Point(0,0,-grav);
  }
};

/** Function object to represent Mass Spring force. */
struct MassSpringForce {
  /** Return mass spring force aggregated over all nodes incident to node n at time t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    auto fspring = Point(0);

    auto first = n.edge_begin();
    auto last = n.edge_end();
    for ( ; first != last; ++first) {
      auto inc_edge = *first;
      auto adj_node = inc_edge.node2();
      fspring -= inc_edge.value().K*((n.position()-adj_node.position())/inc_edge.length())*(inc_edge.length() - inc_edge.value().L);
    }
    return fspring;
  }
};

/** Function object to represent damping force. */
struct DampingForce {
  /** Return force valued at -c*v for node n at time t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (n.value().vel)*Point(0,-cdamp,0);
  }
};

/** Function object to represent no force. */
struct NoForce {
  /** Return Point(0,0,0) for no force in any direction. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return Point(0,0,0);
  }
};

/** Function object to represent a combined force. */
template <typename F1, typename F2, typename F3>
struct CombinedForce {
  F1 f1_;
  F2 f2_;
  F3 f3_;

  /** Return sum of aggregated forces for node n at time t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n,t) + f3_(n,t);
  }
};

/** Return CombinedForce from input forces. */
template <typename F1=NoForce, typename F2=NoForce, typename F3=NoForce>
CombinedForce<F1, F2, F3> make_combined_force(F1 f1_= NoForce(), F2 f2_ = NoForce(), F3 f3_ = NoForce()) {
  return CombinedForce<F1, F2, F3>({f1_, f2_, f3_});
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
  double mi = 1./graph.num_nodes();
  cdamp = 1./graph.num_nodes();
  double kspring = 100; //spring constant
  auto first = graph.node_begin();
  auto last = graph.node_end();

  for ( ; first != last; ++first) {
    auto n = *first;
    n.value().mass = mi;
    n.value().vel = Point(0);
    auto nfirst = n.edge_begin();
    auto nlast = n.edge_end();
    for ( ; nfirst != nlast; ++nfirst) {
      auto adj_e = *nfirst;
      adj_e.value().L = adj_e.length();
      adj_e.value().K = kspring;
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
  auto CF = make_combined_force(GravityForce(), MassSpringForce());
  auto CC = make_combined_constraint(PlaneConstraint(), SphereConstraintRemove());
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, CF, CC);

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
