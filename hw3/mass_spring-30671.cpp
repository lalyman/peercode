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
  NodeData(Point _vel, double _mass) : vel(_vel), mass(_mass) {}
};

//* Custom structure to store spring related data with Edge */
struct EdgeData {
  double rest_length = 1./25.;     //< spring rest length
  double spring_const = 100.;      //< spring hooke's law constant
  EdgeData() {}
  EdgeData(double L, double K) : rest_length(L), spring_const(K) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/*
 * Helper function for debugging purpose. prints the entire graph
 *
 * @param[in]   g     Graph
 * writes to std::cout the complete graph specification.
 * An inspection function that does not modify the graph in anyway.
 */
void print_graph(const GraphType& g) {
  std::cout << g.size() << " nodes and " << g.num_edges() << " edges" << std::endl;
  for (auto np = g.node_begin(); np != g.node_end(); ++np) {
    auto n = *np;
    std::cout << "node index is "<< n.index()
      << ", deg is " << n.degree() << std::endl;
    for (auto ep = n.edge_begin(); ep != n.edge_end(); ++ep) {
      auto v = (*ep).node2();
      std::cout << v.index() << "(" << v.degree() << ") ";
    }
    std::cout << std::endl;
  }
}

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

    if(n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0))
      n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if(n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)){
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/**
 * Constraint is a function object called as operator()(G).
 *
 * @tparam Filter support bool operator()(G::Node) to find nodes that violates the constraint
 * @tparam Transform support void operator()(G::Node) to fix nodes
 */
template<class Filter, class Transform>
struct Constraint {
  Filter filter;
  Transform transform;
  Constraint(Filter _f, Transform _t) : filter(_f), transform(_t) {}
  template<typename G>
  double operator()(G& g, double dt) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (filter(n)) transform(n);
    }
  return dt;
  }
};

/**
 * helper function to make Constraint object from @a f and @a t respectively.
 *
 * @param[in] f support bool operator()(G::Node) to find nodes that violates the constraint
 * @param[in] t support void operator()(G::Node) to fix nodes
 */
template<class F, class T>
Constraint<F, T> make_constraint(F f, T t) {
  return Constraint<F, T>(f, t);
}

/**
 * helper object that is returned by make_composite_constraint()
 *
 * @tparam[in] C1 the first Constraint object to be applied.
 * @tparam[in] C2 the second COnstraint object to be applied.
 */
template<class C1, class C2>
struct CompositeConstraint {
  C1 c1;
  C2 c2;
  CompositeConstraint(C1 _c1, C2 _c2) : c1(_c1), c2(_c2) {}
  template<typename G>
  double operator()(G& g, double dt) {
    c1(g, dt);
    c2(g, dt);
    return dt;
  }
};

/**
 * helper function that composites two or three constraints in sequence.
 *
 * @tparam[in] C1 first Constraint
 * @tparam[in] C2 second Constraint
 * @tparam[in] C3 third COnstraint
 *
 * Note that the constraint are applied in sequence.
 * make_composite_constraint(c1, c2, c3) may not be the same as 
 * make_composite_constraint(c3, c1, c2).
 */
template<class C1, class C2>
CompositeConstraint<C1, C2> make_composite_constraint(C1 c1, C2 c2) {
  return CompositeConstraint<C1, C2>(c1, c2);
};

template<class C1, class C2, class C3>
CompositeConstraint<CompositeConstraint<C1, C2>, C3>  make_composite_constraint(C1 c1, C2 c2, C3 c3) {
  auto c4 = make_composite_constraint(c1, c2);
  return make_composite_constraint(c4, c3);
};

/** Force function object for HW2 #1. */

/**
 * helper object returned by make_combined_force()
 *
 * @tparam[in] F1 an object that supports Point operator()(Graph::Node, double)
 * @tparam[in] F2 another object that supports Point operator()(Graph::Node, double)
 *
 * returns an object that sums the two forces.
 */
template<class F1, class F2>
struct CombinedForce {
  F1 f;
  F2 g;
  CombinedForce() : f(F1()), g(F2()) {}
  CombinedForce(F1 f_, F2 g_) : f(f_), g(g_) {}
  template <typename Node>
  Point operator()(Node n, double t) {
    return f(n, t) + g(n, t);
  }
};

/**
 * helper function to sum two or three forces.
 */
template<class F1, class F2>
CombinedForce<F1, F2> make_combined_force(F1 g, F2 h) {
  return CombinedForce<F1, F2>(g, h);
};

template<class F1, class F2, class F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 a, F2 b, F3 c) {
  auto ab = make_combined_force(a, b);
  return make_combined_force(ab, c);
}

/**
 * forces that simulates the gravity
 * support Point operator()(Node, double).
 */
struct GravityForce {
  template<typename Node>
  Point operator()(Node n, double t) {
    return n.value().mass * Point(0, 0, -grav);
  }
};


/**
 * forces that simulates two point masses connected by a spring.
 * support Point operator()(Node, double).
 */
struct MassSpringForce {
  template<typename Node>
  Point operator()(Node n, double t) {
    Point total_force(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        auto v = (*it).node2();
        double L = (*it).value().rest_length;
        double K = (*it).value().spring_const;
        auto dist = n.position() - v.position();
        total_force -= K * dist * (norm(dist)- L) / norm(dist);
    }
    return total_force;
  }
};

/**
 * forces that simulates damping effect due to friction.
 * support Point operator()(Node, double).
 */
struct DampingForce {
  double c;
  DampingForce(double _c) : c(_c) {}
  template<typename Node>
  Point operator()(Node n, double t) {
    return -c * n.value().vel;
  }
};


struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  //template <typename NODE>
  //Point operator()(NODE n, double t) {
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
        return Point(0, 0, 0);
    Point total_force(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        auto v = (*it).node2();
        double L = (*it).value().rest_length;
        double K = (*it).value().spring_const;
        auto dist = n.position() - v.position();
        total_force -= K * dist * (norm(dist)- L) / norm(dist);
    }
    total_force += n.value().mass * Point(0, 0, -grav);
    return total_force;
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
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  //print_graph(graph);
  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  auto N = graph.size();
  for(auto n = graph.node_begin(); n != graph.node_end(); ++n){
    auto node = *n;
    node.value() = NodeData(Point(0, 0, 0), 1./N);
  }

  for(auto e = graph.edge_begin(); e != graph.edge_end(); ++e) {
    auto edge = *e;
    auto x = edge.node1().position();
    auto y = edge.node2().position();
    edge.value() = EdgeData(norm(x - y), 100.0);
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  //print_graph(graph);
  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  std::cout << "nodes and edge added" << std::endl;
  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;
      //define plane and sphere filter and fix through lambda expression
      auto plane_filter = [](const Node& n) { return n.position().z < -0.75;};
      auto plane_fix = [](Node& n){ n.position().z = -0.75; n.value().vel.z = 0.0;};

      Point circle_center(0.5, 0.5, -0.5);
      double radius = 0.15;
      auto sphere_filter = [circle_center, radius](const Node& n){
        return norm(n.position() - circle_center) < radius;
        };
      auto sphere_fix = [circle_center, radius](Node& n){
        auto u = (n.position() - circle_center)/norm(n.position() - circle_center);
        n.position() = circle_center + u * radius;
        auto &v = n.value().vel;
        v = v - dot(v, u) * u;
        };
      // define constraint and forces
      auto plane_constraint = make_constraint(plane_filter, plane_fix);
      auto sphere_constraint = make_constraint(sphere_filter, sphere_fix);
      auto force = make_combined_force(GravityForce(), MassSpringForce(),
          DampingForce(1./graph.size()));
      auto constraint = make_composite_constraint(plane_constraint, sphere_constraint);
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, Problem1Force());
        //symp_euler_step(graph, t, dt, force);
        constraint(graph, dt);

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
