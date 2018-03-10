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
  double c;        //< Node damping constant
  NodeData() : vel(0), mass(1), c(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double L;       //< Rest Length
  double K;       //< Spring constant
  EdgeData() : L(1), K(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;




/** Constraint function object for a plane */
struct PlaneConstraint {

  PlaneConstraint(double z) : z(z) {}

  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (dot(n.position(), Point(0,0,1)) < z) {
        n.position() = Point(n.position().x, n.position().y, z);
        n.value().vel -= Point(0, 0, n.value().vel.z);
      }
    }
  }

  private:
      double z;
};


/** Constraint function object for a sphere */
struct SphereConstraint {

  SphereConstraint(Point c, double r) : c(c), r(r) {}

  void operator()(GraphType& g, double t) {
    (void) t;
    // iterate through nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point Ri = n.position() - c;
      if (norm(Ri) < r) {
        n.position() = c + r / norm(Ri) * Ri;
        n.value().vel -= dot(n.value().vel, Ri) / normSq(Ri) * Ri;
      }
    }
  }

  private:
      Point c; 
      double r;
};

/** Constraint function which removes nodes within a sphere */
struct SphereRemoveConstraint {

  SphereRemoveConstraint(Point c, double r) : c(c), r(r) {}

  void operator()(GraphType& g, double t) {
    (void) t;
    // iterate through nodes
    for (auto it = g.node_begin(); it != g.node_end();) {
      auto n = *it;
      Point Ri = n.position() - c;
      if (norm(Ri) < r) {
        g.remove_node(n);
      } else {
        ++it;
      }
    }
  }

  private:
      Point c; 
      double r;
};


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
    n.position() += n.value().vel * dt;
  }

  // Apply Constraints
  PlaneConstraint(-0.75)(g, t);
  SphereConstraint(Point(0.5, 0.5, -0.5), 0.15)(g, t);
  SphereRemoveConstraint(Point(0.5, 0.5, -0.5), 0.15)(g, t);
  // NOTE: The above remove node constraint is not plotting correctly, however, 
  // if you iterate through the nodes and print the their positions, everything
  // is as expected.  I have been unable to find the bug.  Is it possible this 
  // is an unexpected interaction with viewer?  Uncomment the following block to
  // see the positions plotted.  I suggest removing forces.
  for (auto it = g.edge_begin(); it != g.edge_end(); ++it) {
     auto e = *it;
     std::cout << "Index1 : " << g.i2u[e.node1().index()] << " Index 2: " << g.i2u[e.node2().index()] << std::endl;
    //  for (auto t = n.edge_begin(); t != n.edge_end(); ++t) {
    //   auto e = *t;
    //   std::cout << "Node1: " << e.node1().index() << " to Node2: " << e.node2().index() << std::endl;
    // }
  }
  std::cout << "\n" << std::endl;



  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
     auto n = *it;
     // Check if node is anchor point
     if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
      // Update the velocity of the node according to its force
       n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}


/** Force function object for sping */
struct MassSpringForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    // Loop through incident edges to compute spring force
    Point spring = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      // get attributes and values of edge
      Point dir = (*it).direction();
      double dist = (*it).length();
      double L = (*it).value().L;
      double K = (*it).value().K;
      // add force of this edge
      spring += -K * dir * (dist - L) / dist;  
    }
    return spring;
  }

};

/** Force function object for Gravity */
struct GravityForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }

};


/** Force function object for Dampening */
struct DampingForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return - n.value().c * n.value().vel;
  }

};

/** Force function object for null force. */
struct Default {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return Point(0,0,0);
  }

};

/** Force function object for total force. */
template <typename F1, typename F2, typename F3>
struct CombinedForce {

  CombinedForce(F1 f1, F2 f2, F3 f3) : f1(f1), f2(f2), f3(f3) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return f1(n, t) + f2(n, t) + f3(n, t);
  }

  private:
      F1 f1;
      F2 f2;
      F3 f3;
};

/** Function that takes up to three force functors and returns a CombinedForce functor. */
template <typename F1 = Default, typename F2 = Default, typename F3 = Default>
CombinedForce<F1, F2, F3> make_combined_force(F1 f1 = Default(), F2 f2 = Default(), F3 f3 = Default())  {
    return CombinedForce<F1, F2, F3>(f1, f2, f3);
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
#if 1
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for nodes.
  auto N = graph.size();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = 1.0 / N;
    n.value().c = 1.0 / N;
  }
  // Set initial conditions for edges, if necessary.
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    e.value().L = e.length();
    e.value().K = 100;
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

        // Compute force and update nodes
        //symp_euler_step(graph, t, dt, make_combined_force());
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), 
                                                          MassSpringForce(), 
                                                          DampingForce()));

        // clear the viewer's nodes and edges
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
