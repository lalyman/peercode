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
static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double L;        // rest length
  double K;        // rigidity constant
  EdgeData() : L(0.0), K(100.0) {}
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


template <typename G, typename F>
struct ForceSum2 {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t);
  }
  G& f1;
  F& f2;

  ForceSum2(G& g, F& f): f1(g), f2(f) {}
};

template <typename G, typename F, typename E>
struct ForceSum3 {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t) + f3(n,t);
  }
  G& f1;
  F& f2;
  E& f3;

  ForceSum3(G& g, F& f, E& e): f1(g), f2(f), f3(e) {}
};

template <typename G, typename F>
ForceSum2<G,F> make_combined_force(G f1, F f2) {
  auto a = ForceSum2<G,F>(f1, f2);
  return a;
}

template <typename G, typename F, typename E>
ForceSum3<G,F,E> make_combined_force(G f1, F f2, E f3) {
  auto a = ForceSum3<G,F,E>(f1, f2, f3);
  return a;
}

template<typename C, typename D>
struct CombinedConstraints2 {
  template<typename G>
  void operator()(G& g, double t){
    c1(g, t);
    c2(g, t);
  }

  C c1;
  D c2;

  CombinedConstraints2(C& c, D& d): c1(c), c2(d) {}
};

template<typename C, typename D, typename E>
struct CombinedConstraints3 {
  template<typename G>
  void operator()(G& g, double t){
    c1(g, t);
    c2(g, t);
    c3(g, t);
  }

  C c1;
  D c2;
  E c3;

  CombinedConstraints3(C& c, D& d, E& e): c1(c), c2(d), c3(e) {}
};

// templating make_combined_constraint to take 2, 3, 4, 5, or 6 arguments
template<typename C, typename D>
CombinedConstraints2<C, D> make_combined_constraint(C c, D d) {
  auto a = CombinedConstraints2<C, D>(c, d);
  return a;
}

template<typename C, typename D, typename E>
CombinedConstraints3<C, D, E> make_combined_constraint(C c, D d, E e) {
  auto a = CombinedConstraints3<C, D, E>(c, d, e);
  return a;
}

template<typename C, typename D, typename E, typename F>
CombinedConstraints3<C, D, CombinedConstraints2<E, F>> make_combined_constraint(C c, D d, E e, F f) {
  auto a = CombinedConstraints2<E, F>(e, f);
  auto b = CombinedConstraints3<C, D, CombinedConstraints2<E, F>>(c, d, a);
  return b;
}

template<typename C, typename D, typename E, typename F, typename G>
CombinedConstraints3<C, CombinedConstraints2<D, E>, CombinedConstraints2<F, G>> make_combined_constraint(C c, D d, E e, F f, G g) {
  auto a = CombinedConstraints2<D, E>(d, e);
  auto b = CombinedConstraints2<F, G>(f, g);
  auto z = CombinedConstraints3<C, CombinedConstraints2<D, E>, CombinedConstraints2<F, G>>(c, a, b);
  return z;
}

template<typename D, typename E, typename F, typename G, typename H, typename I>
CombinedConstraints3<CombinedConstraints2<D, E>, CombinedConstraints2<F, G>, CombinedConstraints2<H, I>> make_combined_constraint(D d, E e, F f, G g, H h, I i) {
  auto a = CombinedConstraints2<D, E>(d, e);
  auto b = CombinedConstraints2<F, G>(f, g);
  auto c = CombinedConstraints2<H, I>(h, i);
  auto z = CombinedConstraints3<CombinedConstraints2<D, E>, CombinedConstraints2<F, G>, CombinedConstraints2<H, I>>(a, b, c);
  return z;
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
    // HW2 #1: YOUR CODE HERE
    (void) t;

    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0);

    Point gravity(0, 0, -grav);
    Point mass_spring(0);

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      mass_spring -= (n.position() - (*it).node2().position()) * ((*it).length() - (*it).value().L) * (*it).value().K / (*it).length();
    }
    return mass_spring + n.value().mass * gravity;
  }
};

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point mass_spring(0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      mass_spring -= (n.position() - (*it).node2().position()) * ((*it).length() - (*it).value().L) * (*it).value().K / (*it).length();
    }
    return mass_spring;
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (-c) * n.value().vel;
  }
};

struct ConstantNodeConstraint {
  Point p;
  double dt;

  ConstantNodeConstraint(Point position, double t): p(position), dt(t) {}

  /**
   * Tells whether the Node n violates the constant node Constraint
   * @param  n  node to check
   * @param  dt time increment that moved it
   * @return    boolean, whether it violates the constraint or not
   *
   * Basically we look for nodes that were at the given position but are no more
   */
  template <typename NODE>
  bool violates(NODE n) {
    return n.position() != p && norm(n.position() - n.value().vel * dt - p) < 1e-5;
  }

  /**
   * Fixes the given node, if need be
   * @param dt time increment parameter
   * @param n  node to fix (if it violates the constraint)
   */
  template <typename NODE>
  void fix_node(NODE n) {
    if(violates(n)){
      n.position() = p;
      n.value().vel = Point(0);
    }
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    // for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
    //   fix_node(*it);
    // }
    auto it = graph.node_begin();
    while(it != graph.node_end() && !violates(*it))
      ++it;
    if(it != graph.node_end())
      fix_node(*it);
  }

};

struct PlaneConstraint {
  Point dir_vec;
  double constant;

  PlaneConstraint(const Point vec, const double c): dir_vec(vec), constant(c) {/*std::cout << "constructed " <<dir_vec << " " << constant<< '\n';*/}

  template <typename NODE>
  bool violates(NODE n) {
    return dot(n.position(), dir_vec) < constant;
  }

  template <typename NODE>
  void fix_node(NODE n) {
    if(violates(n)){
      n.position() -= (dot(n.position(), dir_vec) - constant) * dir_vec;
      n.value().vel -= dot(n.value().vel, dir_vec) * dir_vec;
    }
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
      fix_node(*it);
    }
  }

};

struct MovingPlaneConstraint {
  Point dir_vec;
  double constant;
  double amplitude;
  double period;

  MovingPlaneConstraint(const Point vec, const double c, const double a, const double p): dir_vec(vec), constant(c), amplitude(a), period(p) {}

  template <typename NODE>
  bool violates(NODE n, double t) {
    return dot(n.position(), dir_vec) < constant + amplitude * std::sin(t / period);
  }

  template <typename NODE>
  void fix_node(NODE n, double t) {
    if(violates(n, t)){
      n.position() -= (dot(n.position(), dir_vec) - constant - amplitude * std::sin(t / period)) * dir_vec;
      n.value().vel -= dot(n.value().vel, dir_vec) * dir_vec;
    }
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
      fix_node(*it, t);
    }
  }

};

struct DestructivePlaneConstraint {
  Point dir_vec;
  double constant;

  DestructivePlaneConstraint(Point vec, double c): dir_vec(vec), constant(c) {}

  template <typename NODE>
  bool violates(NODE n) {
    // only cuts when too close to the plane
    return std::abs(dot(n.position(), dir_vec) - constant) < 1e-1;
  }

  template <typename NodeIterator, typename G>
  NodeIterator fix_node(NodeIterator& n, G& graph) {
    if(violates(*n)){
      return graph.remove_node(n);
    }
    return ++n;
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    for(auto it = graph.node_begin(); it != graph.node_end();it = fix_node(it, graph)){}
  }
};

struct SphereConstraint {
  Point center;
  double radius;

  SphereConstraint(Point vec, double c): center(vec), radius(c) {}

  template <typename NODE>
  bool violates(NODE n) {
    return norm(n.position() - center) < radius;
  }

  template <typename NODE>
  void fix_node(NODE n) {
    if(violates(n)){
      n.value().vel -= dot(n.value().vel, n.position() - center) * (n.position() - center) / normSq(n.position() - center);
      n.position() += (n.position() - center) * (radius - norm(n.position() - center)) / norm(n.position() - center);
    }
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
      fix_node(*it);
    }
  }
};

struct DestructiveSphereConstraint {
  Point center;
  double radius;

  DestructiveSphereConstraint(const Point vec, const double c): center(vec), radius(c) {}

  template <typename NODE>
  bool violates(const NODE n) {
    return norm(n.position() - center) < radius;
  }

  template <typename NodeIterator, typename G>
  NodeIterator fix_node(NodeIterator& n, G& graph) {
    if(violates(*n)){
      return graph.remove_node(n);
    }
    return ++n;
  }

  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;
    for(auto it = graph.node_begin(); it != graph.node_end();it = fix_node(it, graph)){}
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

  for (auto e = graph.edge_begin(); e != graph.edge_end(); ++e) {
    (*e).value().L = (*e).length();
  }

  double initial_mass(1.0 / graph.num_nodes());
  c = initial_mass;
  for (auto n = graph.node_begin(); n != graph.node_end(); ++n) {
    (*n).value().mass = initial_mass;
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


      auto constraint = make_combined_constraint(
        ConstantNodeConstraint(Point(0), dt),
        ConstantNodeConstraint(Point(1,0,0), dt),
        MovingPlaneConstraint(Point(0,0,1), -0.5, 0.15, 0.1),
        DestructivePlaneConstraint(Point(1,0,0), 0.5),
        // DestructivePlaneConstraint(Point(0,0,1), -0.75)
        // SphereConstraint(Point(0.5,0.5,-0.5), 0.15)
        SphereConstraint(Point(0,0,-0.4), 0.4)
        // DestructiveSphereConstraint(Point(0.5,0.5,-0.5), 0.15)
      );

      auto force = make_combined_force(GravityForce(), MassSpringForce(), DampingForce());

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        // std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());
        // symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));
        symp_euler_step(graph, t, dt,force,constraint);

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
