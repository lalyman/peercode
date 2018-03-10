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
  double c;        //< damping coefficient
  NodeData() : vel(0), mass(1), c(0.0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;       //< spring constant
    double L;       //< rest length
    EdgeData() : K(100.), L(1.0) {}
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
 * @tparam G::node_value_type supports graph
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
    if(n.position() != Point(0,0,0) && n.position() != Point(1,0,0))
      n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if(n.position() != Point(0,0,0) && n.position() != Point(1,0,0))
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
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force(0,0,0);
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return force;

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double L = (*it).value().L;
      double K = (*it).value().K;
      double dij = dist_l2(n.position(), (*it).node2().position());
      force += -K * (n.position() - (*it).node2().position()) / dij* (dij - L);
    }
    force += Point(0.,0., -grav*n.value().mass);
    return force;

  }
};


/** Force function object for HW2 #3. */
struct GravityForce {
    /** Return the gravity force applying to @a n at time @a t. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
      return Point(0.,0., -grav*n.value().mass);
    }
};

/** Force function object for HW2 #3. */
struct MassSpringForce {
    /** Return the mass spring interaction force applying to @a n at time @a t.*/
    template <typename NODE>
    Point operator()(NODE n, double t) {
      // HW2 #1: YOUR CODE HERE
      Point force(0,0,0);
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        double L = (*it).value().L;
        double K = (*it).value().K;
        double dij = dist_l2(n.position(), (*it).node2().position());
        force += -K * (n.position() - (*it).node2().position()) / dij* (dij - L);
      }
      return force;

    }
};

/** Force function object for HW2 #3. */
struct DampingForce {
    /** Return the damping applying to @a n at time @a t.*/
    template <typename NODE>
    Point operator()(NODE n, double t) {
      return -n.value().c*n.value().vel;

    }
};



/** Combine force class, for 3 forces. */
template <typename T1, typename T2, typename T3 = void>
class combined_force
{
public:
    combined_force<T1, T2, T3> (T1 f1, T2 f2, T3 f3): f1_(f1), f2_(f2), f3_(f3) {}
    template <typename NODE>
    Point operator() (NODE n, double t) { return f1_(n,t) + f2_(n,t) + f3_(n,t); }
private:
    T1 f1_;
    T2 f2_;
    T3 f3_;
};

/** Combine force class, for 2 forces. */
template <typename T1, typename T2>
class combined_force<T1, T2, void>
{
public:
    combined_force<T1, T2> (T1 f1, T2 f2): f1_(f1), f2_(f2) {}
    template <typename NODE>
    Point operator() (NODE n, double t) { return f1_(n,t) + f2_(n,t); }
private:
    T1 f1_;
    T2 f2_;
};

/**Helper function for constructing combined force. This deduces the type of
 * the forces so the user doesn't have to write it.
 * */
template <typename T1, typename T2>
combined_force<T1,T2> make_combined_force(const T1& f1, const T2& f2) {
  return combined_force<T1,T2>(f1, f2);
}

/**Helper function for constructing combined force. This deduces the type of
 * the forces so the user doesn't have to write it.
 * */
template <typename T1, typename T2, typename T3>
combined_force<T1,T2,T3> make_combined_force(const T1& f1, const T2& f2, const T3& f3) {
  return combined_force<T1,T2,T3>(f1, f2, f3);
}




/** Plane constraint function object for HW2 #4. */
struct PlaneConstraint {
    /** Set the position to the nearest point on the plane.
      * Set the z-component of the Node velocity to zero.
      */
    void operator()(GraphType &graph) {
      double z = -0.75;
      for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        auto n = *it;
        if (n.position()[2] < z) {
          n.position()[2] = z;
          n.value().vel[2] = 0;
        }
      }
    }
};



/** Sphere constraint function object for HW2 #4. */
struct SphereConstraint {
    /** Set the position to the nearest point on the surface of the sphere.
     *  Set the component of the velocity that is normal to the sphere’s surface to zero
     */
    void operator()(GraphType &graph) {
      Point c(0.5, 0.5, -0.5);
      double r = 0.15;
      for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        auto n = *it;

        double d = dist_l2(n.position(), c);

        if (d < r) {
          Point Ri = (n.position() - c) / d;
          n.position() = r * Ri + c;
          double vn = inner_prod(n.value().vel, Ri);
          n.value().vel = n.value().vel - vn * Ri;
        }
      }
    }
};



/** Sphere constraint function object , remove nodes for HW2 #4. */
struct SphereRemoveConstraint {
    /** Remove the node and all of its edges using your remove node function.*/

    void operator()(GraphType &graph) {
      Point c(0.5, 0.5, -0.5);
      double r = 0.15;
      for (auto it = graph.node_begin(); it != graph.node_end();) {
        auto n = *it;
        double d = dist_l2(n.position(), c);
        if (d < r)
          graph.remove_node(n);
        else
            ++it;
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

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }


  // Set initial conditions for your nodes, if necessary.
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().mass = 1. / graph.size();
    (*it).value().c = 1. / graph.size();
     for(auto e_it = (*it).edge_begin(); e_it != (*it).edge_end(); ++e_it)
         (*e_it).value().L = dist_l2((*e_it).node1().position(),(*e_it).node2().position());
  }
    // Set initial conditions for your nodes, if necessary.
  for(auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
      (*it).value().L = dist_l2((*it).node1().position(),(*it).node2().position());
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
          //symp_euler_step(graph, t, dt,Problem1Force());
          symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));

        //uncomment them for Plane Constraint
        PlaneConstraint plane_constraint;
        plane_constraint(graph);

        //uncomment them for Sphere Constraint
//        SphereConstraint sphere_constraint;
//        sphere_constraint(graph);

        //uncomment them for Sphere Remove Constraint
          SphereRemoveConstraint sphere_remove_constraint;
          sphere_remove_constraint(graph);
          viewer.clear();
          node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // If you delete nodes, you need to add edges again
        viewer.add_edges (graph.edge_begin() , graph.edge_end () ,node_map);

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
