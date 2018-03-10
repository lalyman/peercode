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
  NodeData(const Point& _vel, double _mass) : vel(_vel), mass(_mass) {}
};

struct EdgeData {
  double rest_len;     //< rest length of this edge
  double K;            //< spring constant of this edge
  EdgeData(double _rest_len = 0, double _K = 100) : rest_len(_rest_len), K(_K) {}
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
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

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
   * model that by returning a zero-valued force. 
   */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point force(0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return force;
    }
    force.z = - n.value().mass * grav;

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Edge e = *it;
      double len = e.length();
      EdgeData edge_data = e.value();
      double coef = edge_data.K * (edge_data.rest_len - len) / len;
      force += coef * (n.position() - e.node2().position());
    }

    return force;
  }
};



// Section: Force

  struct GravityForce {
    /** Return the gravity force applying to @a n at time @a t. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) t;
      return Point(0, 0, - n.value().mass * grav);
    }
  };


  struct MassSpringForce {
    /** Return the spring force applying to @a n at time @a t. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) t;
      Point force(0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Edge e = *it;
        double len = e.length();
        EdgeData edge_data = e.value();
        double coef = edge_data.K * (edge_data.rest_len - len) / len;
        force += coef * (n.position() - e.node2().position());
      }

      return force;
    }
  };

  struct DampingForce {
    /** Return the damping force applying to @a n at time @a t. */
    double c;

    DampingForce(double _c = 0.01) : c(_c) {}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) t;
      return - c * n.value().vel;
    }
  };


  struct NoForce {
    /** Return the spring force applying to @a n at time @a t. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) n, (void) t;
      return Point(0);
    }
  };



  template <typename F1, typename F2, typename F3>
  struct Combined_Force {
    F1 force1; 
    F2 force2;
    F3 force3;

    Combined_Force(F1 _force1, F2 _force2, F3 _force3) : force1(_force1), force2(_force2), force3(_force3) {}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      return force1(n, t) + force2(n, t) + force3(n, t);
    }
  };


  template <typename F1=NoForce, typename F2=NoForce, typename F3=NoForce>
  Combined_Force<F1, F2, F3> make_combined_force(F1 force1=F1(), F2 force2=F2(), F3 force3=F3()) {
    return Combined_Force<F1, F2, F3>(force1, force2, force3);
  }



// Section: Constraints

  struct FixedNodeConstraint {
    double dt;

    FixedNodeConstraint(double _dt = 0.001) : dt(_dt) {}

    void operator()(const GraphType& G, double t) {
      (void) t;
      for (auto nit = G.node_begin(); nit != G.node_end(); ++nit) {
        auto n = *nit;
        if (n.position() == Point(0) + n.value().vel * dt) {
          n.position() = Point(0);
          n.value().vel = Point(0);
        } else if (n.position() == Point(1,0,0) + n.value().vel * dt) {
          n.position() = Point(1,0,0);
          n.value().vel = Point(0);
        }
      }
    }
  };


  struct PlaneConstraint {
    double lower_bound_z;

    PlaneConstraint(double _lower_bound_z = -0.75) : lower_bound_z(_lower_bound_z) {}

    void operator()(const GraphType& G, double t) {
      (void) t;
      for (auto nit = G.node_begin(); nit != G.node_end(); ++nit) {
        auto n = *nit;
        if (n.position().z < lower_bound_z) {
          n.position().z = lower_bound_z;
          n.value().vel.z = 0;
        }
      }
    }
  };

  struct SphereConstraint {
    Point c;
    double rSq;

    SphereConstraint(Point _c = Point(0.5, 0.5, -0.5), double _r = 0.15) : c(_c), rSq(_r * _r) {}

    void operator()(const GraphType& G, double t) {
      (void) t;
      for (auto nit = G.node_begin(); nit != G.node_end(); ++nit) {
        auto n = *nit;
        Point R = n.position() - c;
        if (normSq(R) >= rSq) continue;
        R /= norm(R);
        n.position() = c + R * std::sqrt(rSq);
        n.value().vel -= inner_prod(n.value().vel, R) * R;

      }
    }
  };

  struct SphereConstraint_rm {
    Point c;
    double rSq;

    SphereConstraint_rm(Point _c = Point(0.5, 0.5, -0.5), double _r = 0.15) : c(_c), rSq(_r * _r) {}

    void operator()(GraphType& G, double t) {
      (void) t;
      for (auto nit = G.node_begin(); nit != G.node_end(); ) {
        auto n = *nit;
        Point R = n.position() - c;
        if (normSq(R) < rSq) {
          nit = G.remove_node(nit);
        } else {
          ++nit;
        }
      }
    }
  };

  struct NoConstraint {
    void operator()(const GraphType& G, double t) {
      (void) t;
      (void) G;
    }
  };

  template <typename C1, typename C2, typename C3>
  struct Combined_Constraint {
    C1 c1; 
    C2 c2; 
    C3 c3; 

    Combined_Constraint(C1 _c1, C2 _c2, C3 _c3) : c1(_c1), c2(_c2), c3(_c3) {}

    void operator()(GraphType& G, double t) {
      c1(G, t);
      c2(G, t);
      c3(G, t);
    }
  };

  template <typename C1=NoConstraint, typename C2=NoConstraint, typename C3=NoConstraint>
  Combined_Constraint<C1, C2, C3> make_combined_constraint(C1 c1=C1(), C2 c2=C2(), C3 c3=C3()) {
    return Combined_Constraint<C1, C2, C3>(c1, c2, c3);
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
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit) {
    Node n = *nit;
    n.value() = NodeData(Point(0), 1.0 / graph.num_nodes());
    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit) {
      Edge e = *eit;
      e.value() = EdgeData(e.length(), 100);
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
          make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0 / graph.num_nodes())),
          make_combined_constraint(FixedNodeConstraint(dt), PlaneConstraint(-0.75), SphereConstraint_rm()) );

        viewer.clear();
        node_map.clear();
        // Update viewer with nodes' new positions and edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));

        // std::cout << std::endl << std::endl << std::endl << std::endl;
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
