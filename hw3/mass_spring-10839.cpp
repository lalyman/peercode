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
  double K;       //< Edge spring constant
  double L;     //< Edge rest length
  EdgeData() : K(0), L(0) {}
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
 * @tparam G::node_value_type supports NodeData struct.
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
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint on the graph's nodes
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData struct.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint modifies the value of the graph's nodes, or 
 *           removes them.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {

      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }
  constraint(g, t);
  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Problem1Force(double K, double L) : 
    m_K(K), 
    m_L(L) {
    }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point f_spring = Point(0, 0, 0);
    Point x_i = n.position();
    for (auto e=n.edge_begin(); e!=n.edge_end(); ++e) {
      Point x_j = (*e).node2().position();
      f_spring += - m_K * (x_i - x_j) * (norm(x_i - x_j) - m_L) / norm(x_i - x_j);
    }
    return f_spring + n.value().mass * Point(0, 0, -grav);
  }
private:
  double m_K, m_L;
};

/** Force function object for HW2 #2. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #2, this is a combination of mass-spring force and gravity,
   * using the inherent characteristics of an edge,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point f_spring = Point(0, 0, 0);
    Point x_i = n.position();
    for (auto e=n.edge_begin(); e!=n.edge_end(); ++e) {
      Point x_j = (*e).node2().position();
      double K = (*e).value().K;
      double L = (*e).value().L;
      f_spring += - K * (x_i - x_j) * (norm(x_i - x_j) - L) / norm(x_i - x_j);
    }
    return f_spring + n.value().mass * Point(0, 0, -grav);
  }
};

/** Force function object for HW2 #3. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Force function object for HW2 #3. */
struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point f_spring = Point(0, 0, 0);
    Point x_i = n.position();
    for (auto e=n.edge_begin(); e!=n.edge_end(); ++e) {
      Point x_j = (*e).node2().position();
      double K = (*e).value().K;
      double L = (*e).value().L;
      f_spring += - K * (x_i - x_j) * (norm(x_i - x_j) - L) / norm(x_i - x_j);
    }
    return f_spring;
  }
};

/** Force function object for HW2 #3. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t. */
  DampingForce(double c) : 
    m_c(c) {
    }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -m_c * n.value().vel;
  }
private:
  double m_c;
};

/** A combined force function object for HW2 #3. */
template <typename f1, typename f2>
struct CombinedForce {
  /** Return the combination of two forces applying to @a n at time @a t. */
  CombinedForce(f1 force1, f2 force2) :
    m_force1(force1), 
    m_force2(force2) {
    }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return m_force1(n, t) + m_force2(n, t);
  }
private:
  f1 m_force1;
  f2 m_force2;
};

/** Return a combination of two forces
 * @param[in] f1      @a force1
 * @param[in] f2      @a force2
 * @return a CombinedForce object
 *
 * @tparam f1 type of @a force1
 * @tparam f2 type of @a force2
 */
template <typename f1, typename f2>
CombinedForce<f1, f2> make_combined_force(
  f1 force1, f2 force2) {
  CombinedForce<f1, f2> combined_force(force1, force2);
  return combined_force;
}

/** Return a combination of three forces
 * @param[in] f1      @a force1
 * @param[in] f2      @a force2
 * @param[in] f3      @a force3
 * @return a CombinedForce object
 *
 * @tparam f1 type of @a force1
 * @tparam f2 type of @a force2
 * @tparam f3 type of @a force3
 */
template <typename f1, typename f2, typename f3>
CombinedForce<CombinedForce<f1, f2>, f3> make_combined_force(
  f1 force1, f2 force2, f3 force3) {
  CombinedForce<CombinedForce<f1, f2>, f3> combined_force(
    make_combined_force<f1, f2>(force1, force2), force3);
  return combined_force;
}

/** Constraint function object for HW2 #4. */
struct PlaneConstraint {
  /** Return a plane impenetrability constraint applying to every node of @a g at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    double z = -0.75;
    for (auto n=g.node_begin(); n!=g.node_end(); ++n) {
      if ((*n).position().z < z) {
        (*n).position().z = z;
        (*n).value().vel.z = 0;
      }
    }
  return;
  }
};

/** Constraint function object for HW2 #4. */
struct SphereConstraint {
  /** Return a sphere impenetrability constraint applying to every node of @a g at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto n=g.node_begin(); n!=g.node_end(); ++n) {
      Point u = (*n).position();
      if (norm(u - c) < r) {
        Point normal = (u - c) / norm(u - c);
        (*n).position() = r * normal + c;
        (*n).value().vel -= normal * dot((*n).value().vel, normal);
      }
    }
  return;
  }
};

/** A combined constraint function object for HW2 #4. */
template <typename c1, typename c2>
struct CombinedConstraint {
  /** Return the combination of two constraint applying to @a g at time @a t. */
  CombinedConstraint(c1 constraint1, c2 constraint2) :
    m_constraint1(constraint1), 
    m_constraint2(constraint2) {
    }
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    m_constraint1(g, t);
    m_constraint2(g, t);
    return;
  }
private:
  c1 m_constraint1;
  c2 m_constraint2;
};

/** Return a combination of three constraints
 * @param[in] c1      @a constraint1
 * @param[in] c2      @a constraint2
 * @return a CombinedConstraint object
 *
 * @tparam c1 type of @a constraint1
 * @tparam c2 type of @a constraint2
 */
template <typename c1, typename c2>
CombinedConstraint<c1, c2> make_combined_constraint(
  c1 constraint1, c2 constraint2) {
  CombinedConstraint<c1, c2> combined_constraint(constraint1, constraint2);
  return combined_constraint;
}

/** Return a combination of three constraints
 * @param[in] c1      @a constraint1
 * @param[in] c2      @a constraint2
 * @param[in] c3      @a constraint3
 * @return a CombinedConstraint object
 *
 * @tparam c1 type of @a constraint1
 * @tparam c2 type of @a constraint2
 * @tparam c3 type of @a constraint3
 */
template <typename c1, typename c2, typename c3>
CombinedConstraint<CombinedConstraint<c1, c2>, c3> make_combined_constraint(
  c1 constraint1, c2 constraint2, c3 constraint3) {
  CombinedConstraint<CombinedConstraint<c1, c2>, c3> combined_constraint(
    make_combined_constraint<c1, c2>(constraint1, constraint2), constraint3);
  return combined_constraint;
}

/** Constraint function object for HW2 #5. */
struct SphereRemoveConstraint {
  /** Return a sphere impenetrability constraint applying to every node of @a g at 
   * time @a t. Remove the nodes of the graph not respecting the constraint*/
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    auto n=g.node_begin();
    while (n!=g.node_end()) {
      Point u = (*n).position();
      if (norm(u - c) < r) {
        g.remove_node(*n);
      }
      else {
        ++n;
      }
    }
  return;
  }
};

/** Constraint function object for HW2 #5. */
struct PlaneRemoveConstraint {
  /** Return a plane impenetrability constraint applying to every node of @a g at 
   * time @a t. Remove the nodes of the graph not respecting the constraint*/
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    double z = -0.75;
    auto n=g.node_begin();
    while (n!=g.node_end()) {
      if ((*n).position().z < z) {
        g.remove_node(*n);
      }
      else {
        ++n;
      }
    }
  return;
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto i=graph.node_begin(); i!=graph.node_end(); ++i) {
    (*i).value().vel = Point(0, 0, 0);
    (*i).value().mass = 1.0 / double(graph.num_nodes()); 
  }
  //double K = 100.0;
  //double L = (*graph.edge_begin()).length();

  // HW2 #2 YOUR CODE HERE
  // Set initial conditions for your edges, if necessary.
  for (auto e=graph.edge_begin(); e!=graph.edge_end(); ++e) {
    (*e).value().K = 100.0;
    (*e).value().L = (*e).length();
    (*e).twin_value().K = 100.0;
    (*e).twin_value().L = (*e).length();
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        // HW2 #1
        //symp_euler_step(graph, t, dt, Problem1Force(K, L));

        // HW2 #2
        //symp_euler_step(graph, t, dt, Problem2Force());

        // HW2 #3
        //symp_euler_step(graph, t, dt, 
        //  make_combined_force<GravityForce, MassSpringForce, DampingForce>(
        //    GravityForce(), MassSpringForce(), DampingForce(
        //      1.0 / double(graph.num_nodes()))));

        // HW2 #4
        //symp_euler_step(graph, t, dt, 
        //  make_combined_force<GravityForce, MassSpringForce, DampingForce>(
        //    GravityForce(), MassSpringForce(), DampingForce(
        //      1.0 / double(graph.num_nodes()))),
        //  make_combined_constraint<SphereConstraint, PlaneConstraint>(
        //    SphereConstraint(), PlaneConstraint()));

        // HW2 #5
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, MassSpringForce, DampingForce>(
            GravityForce(), MassSpringForce(), DampingForce(
              1.0 / double(graph.num_nodes()))), SphereRemoveConstraint());
        
        // Clear the viewer’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes’ new positions and new edges
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
