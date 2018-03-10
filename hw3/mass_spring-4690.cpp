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
#include <vector>

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

struct EdgeData {
  double spring_c;
  double rest_length;
  EdgeData() : spring_c(100), rest_length(1) {}
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
 * @tparam C is a function object called as @a constraint(@a g, @a t).
 *           @a constraint does not return anything, but mutates @a g,
 *           passed by reference, to satisfy the constraints where necessary.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    Node n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {

    Node n = *it;

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
  Point operator()(NODE n, double) {
    // HW2 #1: YOUR CODE HERE
    const Point& x = n.position();
    if (x == Point(0,0,0) || x == Point(1,0,0))
      return Point(0,0,0);

    Point spring_force(0,0,0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Edge e = *iter;
      Point diff = x - e.node2().position();
      double diff_norm = norm(diff);

      spring_force -= e.value().spring_c * diff / diff_norm
        * (diff_norm - e.value().rest_length);
    }

    return spring_force + n.value().mass * Point(0,0,-grav);
  }
};

// An abstract struct that provides an interface for making combined force
struct Force {
  virtual Point operator()(Node, double) = 0;
};

struct GravityForce : Force {
  Point operator()(Node n, double) {
    return n.value().mass * Point(0, 0, -grav);
  }
};

struct MassSpringForce : Force {
  Point operator()(Node n, double) {
    Point spring_force(0,0,0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Edge e = *iter;
      Point diff = (n.position() - e.node2().position());
      double diff_norm = norm(diff);
      spring_force -= e.value().spring_c * diff / diff_norm
        * (diff_norm - e.value().rest_length);
    }
    return spring_force;
  }
};

struct DampingForce : Force {
  double damp_c;

  Point operator()(Node n, double) {
    return - damp_c * n.value().vel;
  }

  DampingForce(double c) : damp_c(c) {}
};

struct CombinedForce : Force {
  std::vector<Force*> force_;

  Point operator()(Node n, double t) {
    Point combined(0, 0, 0);
    for (auto& f : force_)
      combined += (*f)(n, t);
    return combined;
  }

  CombinedForce(std::vector<Force*> force) : force_(force) {}
};

/**Combine two forces of different types together into a CombinedForce functor
 * @param[in] f1, f2  a force to be combined
 * @return            a CombinedForce functor that adds up @a f1 and @a f2

 * @tparam F1, F2  a functor struct derived from the abstract Force,
 *                 overriding the pure virtual method operator()(Node, double)

 * @post  the output of @return's operator() is equivalent to that of
 *        calling f1() & f2() separatedly and taking the sum
 */

template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2) {
  std::vector<Force*> force{&f1, &f2};
  return CombinedForce(force);
}

/**Combine three forces of different types together into a CombinedForce functor
 * @param[in] f1, f2, f3  forces to be combined
 * @return                a CombinedForce functor that adds up the three forces
 *
 * @tparam F1, F2, F3  functor structs derived from the abstract Force,
 *                     overriding the pure virtual method operator()(Node, double)
 *
 * @post  the output of @return's operator() is equivalent to that of
 *        calling f1(), f2() & f3() separatedly and taking the sum
 */

template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3) {
  std::vector<Force*> force{&f1, &f2, &f3};
  return CombinedForce(force);
}

struct Constraint {
  virtual void operator()(Graph<NodeData, EdgeData>&, double) = 0;
};

struct ConstantNodeConstraint : Constraint {
  std::vector<Point> anchors_;

  void operator()(Graph<NodeData, EdgeData>& g, double) {
    for (unsigned i = 0; i != anchors_.size(); ++i) {
      Point& x = g.node(i).position();
      if (x != anchors_[i])
        x = anchors_[i];
    }
  }

  ConstantNodeConstraint(std::vector<Point> anchors) : anchors_(anchors) {}
};

struct PlaneConstraint : Constraint {
  double z;

  void operator()(Graph<NodeData, EdgeData>& g, double) {
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      Node n = *iter;
      double distToZ = dot(n.position(), Point(0, 0, 1)) - z;
      if (distToZ < 0) {
        n.position() += Point(0, 0, -distToZ);
        n.value().vel *= Point(1, 1, 0);
      }
    }
  }

  PlaneConstraint(double z_cord) : z(z_cord) {}
};

struct SphereConstraint : Constraint {
  Point c;
  double r;

  void operator()(Graph<NodeData, EdgeData>& g, double) {
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      Node n = *iter;
      double distToC = norm(n.position() - c);
      if (distToC < r) {
        Point R = (n.position() - c) / distToC;
        n.position() = c + R * r;
        n.value().vel -= dot(n.value().vel, R) * R;
      }
    }
  }

  SphereConstraint(Point center, double radius) : c(center), r(radius) {}
};

struct SphereRemovalConstraint : Constraint {
  Point c;
  double r;

  void operator()(Graph<NodeData, EdgeData>& g, double) {
    typename GraphType::node_iterator iter = g.node_begin();
    while (iter != g.node_end()) {
      if (norm((*iter).position() - c) < r)
        g.remove_node(*iter);
      else
        ++iter;
    }
  }

  SphereRemovalConstraint(Point center, double radius)
    : c(center), r(radius) {}
};

struct CombinedConstraint : Constraint {
  std::vector<Constraint*> constraint_;

  void operator()(Graph<NodeData, EdgeData>& g, double t) {
    for (auto& c : constraint_)
      (*c)(g, t);
  }

  CombinedConstraint(std::vector<Constraint*> constraint)
    : constraint_(constraint) {}
};

/**Combine three constraints of different types together
 * into a CombinedConstraint functor
 * @param[in] c1, c2, c3  constraints to be combined
 * @return        void

 * @tparam C1, C2, C3  functor structs derived from the abstract Constraint,
 *                     overriding the pure virtual method operator()
 *                     (Graph<NodeData, EdgeData>, double)

 * @pre  the constraints do not have conflicts that cause infeasibility
 * @post  the output of @return's operator() is equivalent to that of
 *        calling c1(), c2() & c3() separatedly
 */

template <typename C1, typename C2, typename C3>
CombinedConstraint make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  std::vector<Constraint*> constraint{&c1, &c2, &c3};
  return CombinedConstraint(constraint);
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
  for (auto& n : nodes)
    n.value().mass = 1. / graph.num_nodes();

  for (auto iter = graph.edge_begin(); iter != graph.edge_end(); ++iter) {
    Edge e = *iter;
    e.value().rest_length = norm(e.node1().position() - e.node2().position());
  }

  std::vector<Point> anchors {Point(0, 0, 0), Point(1, 0, 0)};
  graph.setAnchors(anchors);

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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt,
          make_combined_force(MassSpringForce(), GravityForce(),
            DampingForce(1./graph.num_nodes())),
          make_combined_constraint(ConstantNodeConstraint(anchors),
            PlaneConstraint(-0.75),
            SphereRemovalConstraint(Point(0.5, 0.5, -0.5), 0.15)));

        // Clear the viewer's nodes and edges
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
