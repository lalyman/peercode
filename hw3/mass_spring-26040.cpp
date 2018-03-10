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
#include <math.h>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
// static constexpr double grav = 0.00;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  double c;        //< Node damping constant
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double k;                // force constant
  double equib_length;     // equilibrium length
  EdgeData() : k(0), equib_length(1) {}
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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    // std::cout << "index " << n.index() << " position " << n.position() << " degree " << n.degree() << std::endl;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  g = constraint(g, t);

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

    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) return Point(0, 0, 0);
    Point force = n.value().mass * Point(0.0, 0.0,-grav);
    // iterate over adjacent nodes
    for (auto bonds_iter = n.edge_begin(); bonds_iter != n.edge_end(); ++bonds_iter) {

      // obtain neighbor
      Node neighbor = (*bonds_iter).node2();
      // calculate euclidean distance 
      double dist = (*bonds_iter).length();
      // obtain force constant
      double k = (*bonds_iter).value().k;
      double equib_length = (*bonds_iter).value().equib_length;
      // calculate force
      force += -k * ((n.position() - neighbor.position()) / dist) * (dist - equib_length);

    }
    (void) t; // silence compiler warning
    return force;
  }
};

/** Gravity Force function object for HW2 #3. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force = n.value().mass * Point(0.0, 0.0,-grav);
    (void) t; // silence compiler warning
    return force;
  }
};

/** Mass Spring Force Force function object for HW2 #3. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force = n.value().mass * Point(0.0, 0.0,0.0);
    // iterate over adjacent nodes
    for (auto bonds_iter = n.edge_begin(); bonds_iter != n.edge_end(); ++bonds_iter) {
      
      // obtain neighbor
      Node neighbor = (*bonds_iter).node2();
      // calculate euclidean distance 
      double dist = (*bonds_iter).length();
      // obtain force constant and equilibrium length
      double k = (*bonds_iter).value().k;
      double equib_length = (*bonds_iter).value().equib_length;
      // calculate force
      force += -k * ((n.position() - neighbor.position()) / dist) * (dist - equib_length);

    }
    (void) t; // silence compiler warning
    return force;
  }
};

/** Damping Force function object for HW2 #3. */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force = -n.value().c * n.value().vel;
    (void) t; // silence compiler warning
    return force;
  }  
};

struct PlaneZ {
  template<typename GRAPH>
  GRAPH operator()(GRAPH g, double t) {
    double z = -0.75;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position().z < z) {
        n.position().z = z;
      }
    }  
  (void) t; // silence compiler warning
  return g;
  }
};

struct Sphere {
  template<typename GRAPH>
  GRAPH operator()(GRAPH g, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      const Point& point = (*it).position();
      double dist = sqrt(pow(point.x-center.x, 2) + pow(point.y-center.y, 2) + pow(point.z-center.z, 2));
      if (dist < r) {
        Point R = (point - center) / dist;
        (*it).position()  = R * r + center;
        (*it).value().vel -= std::inner_product((*it).value().vel.begin(), (*it).value().vel.begin(), R.begin(), 0) * R;
      }
    }  
  (void) t; // silence compiler warning
  return g;
  }
};

struct CutSphere {
  template<typename GRAPH>
  GRAPH operator()(GRAPH g, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      const Point& point = (*it).position();
      double dist = sqrt(pow(point.x-center.x, 2) + pow(point.y-center.y, 2) + pow(point.z-center.z, 2));
      if (dist < r) g.remove_node(*it);
    }  
  (void) t; // silence compiler warning
  return g;
  }
};

template <typename C1, typename C2>
class make_combined_constraints {
  public:
    // Constructor
    make_combined_constraints(const C1& c1, const C2& c2)
        : c1_(c1), c2_(c2){}

    template <typename GRAPH>
    GRAPH operator()(GRAPH g, double t) {
      g = c1_(g, t);
      g = c2_(g, t);
    return g;
    }

  private:
    C1 c1_;
    C2 c2_;
};

template <typename F1, typename F2>
class make_combined_force {
  public:
    // Constructor
    make_combined_force(const F1& f1, const F2& f2)
        : f1_(f1), f2_(f2){}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      Point force1 = f1_(n, t);
      Point force2 = f2_(n, t);
    return force1 + force2;
    }

  private:
    F1 f1_;
    F2 f2_;
};

template <typename F1, typename F2, typename F3>
class make_combined_force2 {
  public:
    // Constructor
    make_combined_force2(const F1& f1, const F2& f2, const F3& f3)
        : f1_(f1), f2_(f2), f3_(f3){}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      Point force1 = f1_(n, t);
      Point force2 = f2_(n, t);
      Point force3 = f3_(n, t);
    return force1 + force2 + force3;
    }

  private:
    F1 f1_;
    F2 f2_;
    F3 f3_;
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
    // graph.add_edge(nodes[t[0]], nodes[t[3]]);
    // graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for nodes, if necessary.
  for (auto it2 = graph.node_begin(); it2 != graph.node_end(); ++it2) {
    (*it2).value() = NodeData();
    (*it2).value().mass = 1.0 / graph.num_nodes() ;
    (*it2).value().c = 1.0 / graph.num_nodes() ;
    (*it2).value().vel = Point(0,0,0) ;
  }
  for (auto it3 = graph.edge_begin(); it3 != graph.edge_end(); ++it3) {
    (*it3).value() = EdgeData();
    (*it3).value().k = 50.0 ;
    (*it3).value().equib_length = (*it3).length() ;
  }

  // combine forces and constraints
  make_combined_force2<GravityForce, MassSpringForce, DampingForce> getforce2  = make_combined_force2<GravityForce, MassSpringForce, DampingForce>(GravityForce(), MassSpringForce(), DampingForce());
  make_combined_constraints<PlaneZ, CutSphere> applyconstraint  = make_combined_constraints<PlaneZ, CutSphere>(PlaneZ(), CutSphere());

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  std::cout << sizeof(GraphType::Node) << " " << sizeof(GraphType::Edge) << std::endl;

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

      Point::size_type prev_num_nodes = graph.num_nodes();
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        symp_euler_step(graph, t, dt, getforce2, applyconstraint);

        // Update viewer with nodes' new positions
        if (graph.num_nodes() != prev_num_nodes) {
          viewer.clear();
          node_map.clear();
        }

        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
        prev_num_nodes = graph.num_nodes();
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
