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

struct EdgeData{
  double K;
  double L;
  EdgeData(): K(0), L(0) {};
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

// NOTE: There are several warnings generated because of not using the paramter
//       t in the forces or the contraints. In fact, we do not need it,
//       although the specification PDF repeatedly mentions it.
/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports node iterators accessible via
 *                            node_begin() and node_end().
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
    n.position() += n.value().vel * dt;
  }
  // Check for constraint violations
  constraint(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
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
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    else{
      return spring_force(n) + gravity_force(n);
    }
  }
  private:
    template <typename NODE>
    Point spring_force(NODE n){
      Point f = Point(0,0,0);
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto e = *it;
        double n = norm(e.node1().position()-e.node2().position());
        f += -e.value().K*(e.node1().position()-e.node2().position())/n*(n - e.value().L);
      }
      return f;
    }

    template <typename NODE>
    Point gravity_force(NODE n){
      return n.value().mass * grav * Point(0,0,-1);
    }
};

// This whole structure to add forces is based on a conversation I had with
// Chris Lazarus, which is based on a discussion he had with Chris Ceka while
// taking the class last year

// Abstract struct from which all out Forces will inherit
struct BaseForce{
  virtual Point operator()(Node n, double t) = 0;
  virtual ~BaseForce() {}
};

// A ProxyForce which ensures inheritance when we are adding forces
template <typename F>
struct ProxyForce : BaseForce{
  F f_;
  ProxyForce(const F& f): f_(f) {}
  Point operator()(Node n, double t){return f_(n,t);}
};

// Model of the force of gravity
struct GravityForce : BaseForce{
  Point operator()(Node n, double t) {
    return n.value().mass * grav * Point(0,0,-1);
  }
};

// Model of the force of the spring
struct MassSpringForce : BaseForce{
  Point operator()(Node n, double t) {
    Point f = Point(0,0,0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto e = *it;
      double n = norm(e.node1().position()-e.node2().position());
      f += -e.value().K*(e.node1().position()-e.node2().position())/n*(n - e.value().L);
    }
    return f;
  }
};

// Model of disipation of energy. Uses global variable c, to maintain the same
// signature on the operator()
struct DampeningForce : BaseForce{
  Point operator()(Node n, double t) {
    return -c*n.value().vel;
  }
};

// Holds a representation of the forces to be applied and applies all of them
// when called on the operator()
struct ForceCollection : BaseForce{
  std::vector<BaseForce*> forces = {};

  ForceCollection() {}

  // Adds a ProxyForce instance, allocated on the heap, to actually own the
  // functor added to the collection
  template <typename F>
  void add_force(F f){
    forces.push_back(new ProxyForce<F>(f));
  }

  // Applies all the operators of the stored forces
  Point operator()(Node n, double t){
    Point p = Point(0,0,0);
    for(auto it = forces.begin(); it != forces.end(); ++it){
      auto f = *it;
      p += (*f)(n,t);
    }
    return p;
  }
};

// Adds two forces to a ForceCollection and returns it
template <typename F1, typename F2>
ForceCollection make_combined_force(F1 f1, F2 f2){
  ForceCollection fc = ForceCollection();
  fc.add_force(f1);
  fc.add_force(f2);
  return fc;
}

// Adds three forces to a ForceCollection and returns it
template <typename F1, typename F2, typename F3>
ForceCollection make_combined_force(F1 f1, F2 f2, F3 f3){
  ForceCollection fc = ForceCollection();
  fc.add_force(f1);
  fc.add_force(f2);
  fc.add_force(f3);
  return fc;
}

// Constraints are based on the same class structure as Forces
// Abstract struct from which all out Forces will inherit
struct BaseConstraint{
  virtual void operator()(GraphType& g, double t) = 0;
  virtual ~BaseConstraint() {}
};

// A ProxyConstraint which ensures inheritance when we are adding contraints
template <typename F>
struct ProxyConstraint : BaseConstraint{
  F f_;
  ProxyConstraint(const F& f): f_(f) {}
  void operator()(GraphType& g, double t){return f_(g,t);}
};

// A constraint to model the two corners of the blanket not moving
struct ConstantNodesConstraint : BaseConstraint{
  void operator()(GraphType& g, double t) {
    // TODO: optimize this
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      double dt = 0.001;
      if (n.position() - n.value().vel * dt == Point(0,0,0)){
        n.position() = Point(0,0,0);
        n.value().vel = Point(0,0,0);
      }
      if (n.position() - n.value().vel * dt == Point(1,0,0)){
        n.position() = Point(1,0,0);
        n.value().vel = Point(0,0,0);
      }
    }
  }
};

// Constraint to simulate the blanket touching a plane
struct PlaneConstraint : BaseConstraint{
  void operator()(GraphType& g, double t) {
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if (n.position()[2] <= -0.75){
        n.position()[2] = -0.75; // TODO: not really, take inner product
        n.value().vel[2] = 0;
      }
    }
  }
};

// Constraint to simulate a sphere the blanket falls onto
struct SphereConstraint : BaseConstraint{
  void operator()(GraphType& g, double t) {
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(norm(n.position() - Point(0.5, 0.5, -0.5)) < 0.15){
        Point R = (n.position() - Point(0.5, 0.5, -0.5))/norm(n.position() - Point(0.5, 0.5, -0.5));
        n.position() =  Point(0.5, 0.5, -0.5) + R * 0.15;
        n.value().vel -= inner_prod(R, n.value().vel)*R;
      }
    }
  }
};

// A constraint to remove every node that touches a sphere centered at
// (0.5, 0.5, -0.5) of radius 0.15
struct SphereRemoveConstraint : BaseConstraint{
  void operator()(GraphType& g, double t) {
    for(auto it = g.node_begin(); it != g.node_end(); ){
      auto n = *it;
      if(norm(n.position() - Point(0.5, 0.5, -0.5)) < 0.15){
        it = g.remove_node(it);
      }else{
        ++it;
      }
    }
  }
};

// Holds a representation of the constraints to be applied and applies all of
// then when called on the operator()
struct ConstraintCollection : BaseConstraint{
  std::vector<BaseConstraint*> constraints = {};

  ConstraintCollection() {}

  // Adds a ProxyConstraint instance, allocated on the heap, to actually own the
  // functor added to the collection
  template <typename C>
  void add_constraint(C c){
    constraints.push_back(new ProxyConstraint<C>(c));
  }

  // Applies all the constraints to the given Graph
  void operator()(GraphType& g, double t){
    for(auto it = constraints.begin(); it != constraints.end(); ++it){
      auto c = *it;
      (*c)(g,t);
    }
  }
};

// Adds two constraints to a ConstraintCollection and returns it
template <typename C1, typename C2>
ConstraintCollection make_combined_constraint(C1 c1, C2 c2){
  ConstraintCollection cc = ConstraintCollection();
  cc.add_constraint(c1);
  cc.add_constraint(c2);
  return cc;
}

// Adds three constraints to a ConstraintCollection and returns it
template <typename C1, typename C2, typename C3>
ConstraintCollection make_combined_constraint(C1 c1, C2 c2, C3 c3){
  ConstraintCollection cc = ConstraintCollection();
  cc.add_constraint(c1);
  cc.add_constraint(c2);
  cc.add_constraint(c3);
  return cc;
}

// TODO: how to predeclare ConstantNodesConstraint to be able to put this code
//       at the top
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  return symp_euler_step(g, t, dt, force, ConstantNodesConstraint());
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
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = 1.0/graph.size();
    // Note that this is preferred over iterating over the edges in the graph
    // since the edge_iterator skips invalid nodes, whereas the
    // incident_iterator will always iterate over all neighbors of a node.
    for(auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it){
      auto e = *e_it;
      e.value().K = 100;
      e.value().L = e.length();
    }
  }

  c = 1.0/graph.num_nodes();

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
        //symp_euler_step(graph, t, dt, Problem1Force());
        ForceCollection force1 = make_combined_force(GravityForce(), MassSpringForce());
        ForceCollection force2 = make_combined_force(GravityForce(), MassSpringForce(), DampeningForce());
        ConstraintCollection c1 = make_combined_constraint(ConstantNodesConstraint(), PlaneConstraint());
        ConstraintCollection c2 = make_combined_constraint(ConstantNodesConstraint(), SphereConstraint());
        ConstraintCollection c3 = make_combined_constraint(ConstantNodesConstraint(), SphereRemoveConstraint());

        //symp_euler_step(graph, t, dt, force2);
        symp_euler_step(graph, t, dt, force2, c3);

        // Clear the viewer, only when running the SphereRemoveConstraint
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // Comment out this line when not using the SphereRemoveConstraint
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
