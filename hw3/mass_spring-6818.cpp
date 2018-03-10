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

// Define the Graph type
using GraphType = Graph<NodeData, double>;
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
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
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    double K = 100.;
    Point massSpringForce(0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      massSpringForce += - K * ((*it).node1().position() - (*it).node2().position()) 
          / (*it).length() * ((*it).length() - (*it).value());
    }
    Point gravityForce(0, 0, - n.value().mass * grav);
    return massSpringForce + gravityForce;
  }
};


/** Gravity Force function object. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0, 0, - n.value().mass * grav);
  }
};

/** Mass Spring Force function object. */
struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    double K = 100.;
    Point massSpringForce(0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      massSpringForce += - K * ((*it).node1().position() - (*it).node2().position()) 
          / (*it).length() * ((*it).length() - (*it).value());
    }
    return massSpringForce;
  }
};


/** Damping Force function object. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return - 1./n.graph_size() * n.value().vel;
  }
};

/** Paired Force function object. */
template <typename F1, typename F2>
struct PairedForce {
  /** Return the Paired force applying to @a n at time @a t.
   */
  F1 f1;
  F2 f2;
  PairedForce(F1 f1, F2 f2) : f1(f1), f2(f2) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t);
  }
};

/** Construct a paired force. */
template <typename F1, typename F2>
PairedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return PairedForce<F1, F2>(f1, f2);
}


/** Tupled Force function object. */
template <typename F1, typename F2, typename F3>
struct TupledForce {
  /** Return the Tupled force applying to @a n at time @a t.
   */
  F1 f1;
  F2 f2;
  F3 f3;
  TupledForce(F1 f1, F2 f2, F3 f3) : f1(f1), f2(f2), f3(f3) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
};

/** Construct a tupled force. */
template <typename F1, typename F2, typename F3>
TupledForce<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return TupledForce<F1, F2, F3>(f1, f2, f3);
}


/** Construct a constant node constraint. */
struct ConstantNodeConstraint {
  template <typename GraphType>
  void operator()(GraphType& graph) {
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  auto n = *it;
      if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
	    n.value().vel = Point(0);
	  }	  
	}
    return;
  }
} constant_node_constraint;


/** Construct a plane constraint. */
struct PlaneConstraint {
  template <typename GraphType>
  void operator()(GraphType& graph) {
    double z = - 0.75;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  auto n = *it;
      if (dot(n.position(), Point(0, 0, 1)) < z) {
	    n.position() += Point(0, 0, z - dot(n.position(), Point(0, 0, 1)));
	    n.value().vel += Point(0, 0, - dot(n.value().vel, Point(0, 0, 1)));
	  }	  
	}
    return;
  }
} plane_constraint;

/** Construct a sphere constraint. */
struct SphereConstraint {
  template <typename GraphType>
  void operator()(GraphType& graph) {
    Point c(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  auto n = *it;
      if (norm(n.position() - c) < r) {
        Point RP = n.position() - c;
        RP = RP / norm(RP);
	    n.position() = c + r * RP;
	    n.value().vel = n.value().vel - dot(n.value().vel, RP) * RP;
	  }	  
	}
    return;
  }
} sphere_constraint;


/** Construct a sphere remove constraint. */
struct SphereRemoveConstraint {
  template <typename GraphType>
  void operator()(GraphType& graph) {
    Point c(0.5, 0.5, -0.5);
    double r = 0.15;
    auto it = graph.node_begin();
    while (it != graph.node_end()) {
	  auto n = *it;
      if (norm(n.position() - c) < r) {
        it = graph.remove_node(it);
	  } else {
	    ++it;
	  }
	}
    return;
  }
} sphere_remove_constraint;




/** Paired Constraint. */
template <typename C1, typename C2>
struct PairedConstraint {
  C1 c1;
  C2 c2;
  PairedConstraint(C1 c1, C2 c2) : c1(c1), c2(c2) {}
  template <typename GraphType>
  void operator()(GraphType& graph) {
    c1(graph);
    c2(graph);
    return;
  }
};

/** Construct a paired constraint. */
template <typename C1, typename C2>
PairedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return PairedConstraint<C1, C2>(c1, c2);
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
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().vel = Point(0);
    (*it).value().mass = 1. / graph.num_nodes();
  }
  
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    (*it).value() = (*it).length();
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
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));
	    make_combined_constraint(sphere_remove_constraint, constant_node_constraint)(graph);
	    
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
