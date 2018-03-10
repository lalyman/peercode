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
// Spring constant
double K = 100;
// Rest length - value to be assigned in main()
double L;
// Damping constant - value assigned in main()
double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
// using GraphType = Graph<NodeData>;
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)){
        n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)){
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  constraint(g,t);

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

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }

    /** Gravity force */
    // Point f_grav = n.value().mass*Point(0,0,-1*grav));
    Point f_grav = Point(0,0,-1.0)*n.value().mass*grav;

    /** Spring force. Iterating through neighboring nodes. */
    Point f_spring = Point(0,0,0);
    for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit){
        Edge e = *iit;
        NODE n1 = e.node1();
        NODE n2 = e.node2();
        f_spring += (n1.position()-n2.position()) * (-1.0*K*(e.length() - e.value())/e.length());
    }
    return (f_grav + f_spring);
  }
};

/** Gravity force function object. */
struct GravityForce{
    template <typename NODE>
    Point operator()(NODE n, double t) {
        return Point(0,0,-1.0)*grav*n.value().mass;
    }
};

/** Mass-spring force function object. */
struct MassSpringForce{
    template <typename NODE>
    Point operator()(NODE n, double t) {
        Point f_spring = Point(0,0,0);
        for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit){
            Edge e = *iit;
            NODE neighbor = e.node2();
            f_spring += (n.position()-neighbor.position()) * (-1.0*K*(norm(n.position()-neighbor.position()) - e.value())/norm(n.position()-neighbor.position()));
        }
        return f_spring;
    }
};

/** Damping force function object. */
struct DampingForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        return (-1.0)*c*n.value().vel;
    }
};

/** Default "Dummy" force function object. operator() returns 0-force. */
struct DummyForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        return Point();
    }
};


/** Combining force function object. */
template <typename f1, typename f2, typename f3>
struct CombineForces {
    // Public Constructor.
    CombineForces(const f1& f_g, const f2& f_s, const f3& f_d) : f_g(f_g),
    f_s(f_s), f_d(f_d) {}

    template <typename NODE>
    Point operator()(NODE n, double t){
        return f_g(n,t) + f_s(n,t) + f_d(n,t);
    }

    private:
    f1 f_g;
    f2 f_s;
    f3 f_d;
};

/** Combining force function. Returns CombineForces object. */
template <typename f1 = DummyForce, typename f2 = DummyForce, typename f3 =
DummyForce>
CombineForces<f1, f2, f3> make_combined_force
(const f1& f_g = DummyForce(), const f2& f_s = DummyForce(), const f3& f_d =
DummyForce()) {
  // Takes in variable number of inputs depending on how many forces are
  // relevant to the problem. If any are not added, add 0-force.
  return CombineForces<f1, f2, f3>(f_g, f_s, f_d);
};

/** Planar constraint function object. */
struct PlaneConstraint {
    template <typename G>
    void operator()(G& g, double t){
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit){
            if ((*nit).position().z < -0.75) {
                (*nit).position().z = -0.75;
                (*nit).value().vel.z = 0.0;
            }
        }
        return;
    }
};

/** Spherical constraint function object. */
struct SphereConstraint{
    template <typename G>
    void operator()(G& g, double t){
        Point center = Point(0.5,0.5,-0.5);
        double radius = 0.15;
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit){
            if (norm((*nit).position() - center) < radius) {
                Point unit = ((*nit).position() - center)/norm((*nit).position() - center);
                (*nit).position() = center + radius * unit;
                (*nit).value().vel -= dot((*nit).value().vel, unit)*unit;
            }
        }
        return;
    }
};

/** Spherical rip function object. Nodes which touch sphere are removed from
graph. */
struct SphereConstraint_Rip{
    template <typename G>
    void operator()(G& g, double t){
        Point center = Point(0.5,0.5,-0.5);
        double radius = 0.15;
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit){
            if (norm((*nit).position() - center) < radius) {
                g.remove_node((*nit));
            }
        }
        return;
    }
};

/** Default constraint function object. Has no effect on the graph. */
struct DummyConstraint {
    template <typename G>
    void operator()(G& g, double t){
        return;
    }
};

/** Combining force function object. */
template <typename c1, typename c2>
struct CombineConstraints {
    // Public Constructor
    CombineConstraints(const c1& c_p, const c2& c_s) : c_p(c_p), c_s(c_s) {}

    template <typename G>
    void operator()(G& g, double t){
        c_s(g,t);
        c_p(g,t);
        return;
    }

    private:
    c1 c_p;
    c2 c_s;
};

/** Combining constraint function. Returns CombineConstraints object. */
template <typename c1 = DummyConstraint, typename c2 = DummyConstraint>
CombineConstraints<c1, c2> make_combined_constraints
(const c1& c_p = DummyConstraint(), const c2& c_s = DummyConstraint()) {
  // Takes in variable number of inputs depending on how many constrints are
  // relevant to the problem. If any are not added, treat as null-constraint.
  return CombineConstraints<c1, c2>(c_p, c_s);
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
  graph.clear();

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

  // Getting spring rest length from edge 0 (as all edges should have the same original rest length) -- Problem 1.
  L = (*(graph.edge_begin())).length();

  // Calculating damping constant as given in Problem 2.
  c = 1.0/graph.num_nodes();

  // Setting mass & initial velocities of the nodes
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
      Node n = *nit;
      n.value().vel = Point(0,0,0);
      n.value().mass = 1.0/graph.num_nodes();
  }

  // Storing the original edge lengths as values. These are the rest-lengths of the springs.
  for (auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit){
      Edge e = *eit;
      e.value() = e.length();
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
        // std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());

        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()), make_combined_constraints(PlaneConstraint(),SphereConstraint_Rip()));

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));

        // Clear the viewer’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes’ new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
