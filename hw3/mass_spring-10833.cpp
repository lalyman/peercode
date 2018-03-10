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

// Mass of Graph
double mi;

//Spring Constant
//double K;

//Damping Coefficient
double c;
//Spring Rest-Length
//double L;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Edge spring constant
  double L;     //< Edge length
  EdgeData() : K(0), L(0) {}
  EdgeData(double k, double l) : K(k), L(l) {}
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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  constraints(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        continue;
    }
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
// struct Problem1Force {
//   /** Return the force applying to @a n at time @a t.
//    *
//    * For HW2 #1, this is a combination of mass-spring force and gravity,
//    * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
//    * model that by returning a zero-valued force. */
//   template <typename NODE>
//   Point operator()(NODE n, double t) {
//     // HW2 #1: YOUR CODE HERE
//     Point total_force = Point(0);
//     Point gravity_force = Point(0);
//     Point spring_force = Point(0);
//     //Contrain the two corners of the cloth by returning a zero force//
//     if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
//         return Point(0,0,0);
//      }
//     gravity_force = (n.value().mass)*Point(0,0,-grav);
//
//     for (auto in_it = n.edge_begin(); in_it != n.edge_end(); ++in_it) {
//         auto edge = *in_it;
//
//         spring_force += -edge.value().K*(
//                         (n.position() - edge.node2().position()) /
//                         (norm(n.position() - edge.node2().position()) )) *
//                         (norm(n.position() - edge.node2().position() )
//                         - edge.value().L);
//     }
//     total_force = spring_force + gravity_force;
//     return total_force;
//     }
// };
struct GravityForce {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //Point total_force = Point(0);
    //Point gravity_force = Point(0);
    //Point spring_force = Point(0);
    //Contrain the two corners of the cloth by returning a zero force//


    (void) t;
    //total_force = gravity_force + spring_force;
    return (n.value().mass)*Point(0,0,-grav);
    }
};
struct MassSpringForce {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //Point total_force = Point(0);
    //Point gravity_force = Point(0);
    Point spring_force = Point(0);
    //Contrain the two corners of the cloth by returning a zero force//
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
    //     return Point(0,0,0);
    // }
    //gravity_force = (n.value().mass)*Point(0,0,-grav);

    for (auto in_it = n.edge_begin(); in_it != n.edge_end(); ++in_it) {
        auto edge = *in_it;
        //std::cout << "edge.value().K " << edge.value().K <<'\n';
        //std::cout << "edge.value().L " << edge.value().L <<'\n';

        spring_force += -edge.value().K*(
                        (n.position() - edge.node2().position()) /
                        (norm(n.position() - edge.node2().position()) )) *
                        (norm(n.position() - edge.node2().position() )
                        - edge.value().L);
    }
    (void) t;
    return spring_force;
    }
};
struct DampingForce {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //Point total_force = Point(0);
    //Point gravity_force = Point(0);
    Point damping_force = Point(0,0,0);
    //Contrain the two corners of the cloth by returning a zero force//
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
    //     return Point(0,0,0);
    // }
    //gravity_force = (n.value().mass)*Point(0,0,-grav);
    damping_force = -c*n.value().vel;
    (void) t;
    return damping_force;
    }
};
struct ZeroForce {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #3: YOUR CODE HERE
    (void) n;
    (void) t;
    return Point(0,0,0);
    }
};

//Define functor for combined force. Make combined force will return this type
// and call its constructor when returning a value (the functor)
template<typename F1, typename F2, typename F3>
struct CombinedForce {
// template<typename F1, typename F2, typename F3>
        CombinedForce(F1 f1, F2 f2, F3 f3)
            : f1_(f1), f2_(f2), f3_(f3) {}
        template<typename NODE>
        Point operator()(NODE n, double t) {
            return f1_(n,t) + f2_(n,t) + f3_(n,t);
            }
        private:
            F1 f1_;
            F2 f2_;
            F3 f3_;
};
template<typename F1= ZeroForce, typename F2= ZeroForce, typename F3= ZeroForce>
CombinedForce<F1,F2,F3> make_combined_force(F1&& f1 = ZeroForce(),
                                            F2&& f2 = ZeroForce(),
                                            F3&& f3 = ZeroForce())
{
    return CombinedForce<F1,F2,F3>(f1,f2,f3);

}

struct ZeroConstraint {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #4, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename G>
  void operator()(G& g, double t) {
    (void) g;
    (void) t;
    return;
    }
};
struct PlaneConstraint {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #4, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename G>
  void operator()(G& g, double t) {
    double z = -0.75;

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        if (n.position().z < z) {
          n.position().z = z;
          n.value().vel.z = 0;
        }

    }
    (void) t;

    }
};

struct SphereConstraint {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #4, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
   template <typename G>
   void operator()(G& g, double t) {
       Point ctr = Point(0.5, 0.5, -0.5);
       double r = 0.15;

       for (auto it = g.node_begin(); it != g.node_end(); ++it) {
           auto n = *it;
           if ((norm(n.position() - ctr)) < r) {
             Point Ri = (n.position() - ctr) / (norm(n.position() - ctr));
             n.position() = ctr + Ri*r;
             n.value().vel = n.value().vel - (dot(n.value().vel,Ri)*Ri);

           }

        }
        (void) t;
        return;
    }
};

struct NewSphereConstraint {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #4, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
   template <typename G>
   void operator()(G& g, double t) {
       Point ctr = Point(0.5, 0.5, -0.5);
       double r = 0.15;

       for (auto it = g.node_begin(); it != g.node_end(); ++it) {
           auto n = *it;
           if ((norm(n.position() - ctr)) < r) {
             g.remove_node(n);
           }

        }
        (void) t;
        return;
    }
};
//Define functor for combined force. Make combined force will return this type
// and call its constructor when returning a value (the functor)
template<typename C1, typename C2>
struct CombinedConstraint {
// template<typename F1, typename F2, typename F3>
        CombinedConstraint(C1 c1, C2 c2)
            : c1_(c1), c2_(c2) {}
        template<typename G>
        void operator()(G& g, double t) {
            c1_(g,t);
            c2_(g,t);
            return;
            }
        private:
            C1 c1_;
            C2 c2_;
};
template<typename C1= ZeroConstraint, typename C2= ZeroConstraint>

CombinedConstraint<C1,C2> make_combined_const(C1&& c1 = ZeroConstraint(),
                                            C2&& c2 = ZeroConstraint())
{
    return CombinedConstraint<C1,C2>(c1,c2);

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
    mi = 1.0/(double)graph.num_nodes();
    c = 1./(double)graph.num_nodes();
    // Zero initial velocity
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      auto n = *it;
      n.value().mass = mi;
      n.value().vel = Point(0);
    }

    for (auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit) {
        auto edge = *eit;
        edge.value().K = 100.;
        edge.value().L = edge.length();
    }
    //Spring Constant
    //K = 100.;

    //Spring Rest-Length
    //L = (*graph.edge_begin()).length();


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
                        make_combined_force(GravityForce(),MassSpringForce()),
                        make_combined_const(PlaneConstraint(),NewSphereConstraint()));

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
