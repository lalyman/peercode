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

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;     //< Spring constant
  double L;     //< Rest-length
  EdgeData() : K(100), L(0.2) {}
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
template <typename G, typename F, typename C>
  // Compute the t+dt position
  double symp_euler_step(G& g, double t, double dt, F force, C constraint = C()) {
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraint
  constraint(g,t);

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
    (void)t;

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    Point springForce;
    for (auto adj = n.edge_begin(); adj != n.edge_end(); ++adj){
      NODE n2 = (*adj).node2();
      double scalar = -(*adj).value().K / norm(n.position() - n2.position()) * (norm(n.position() - n2.position()) - (*adj).value().L);
      springForce += ((n.position() - n2.position()) * scalar);
    }
    Point result = springForce + n.value().mass*Point(0,0,-grav);
    return result;
  }
};

/* Constant constraint specifying that a given node's position remain fixed */
template <typename NODE>
struct ConstantConstraint{
 public:
  /* Constructor */
  ConstantConstraint(NODE n):node_(n),p_(n.position()){};

  template<typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void)t; (void)g;
    // reset the node's position to its original position
    node_.position() = p_;

    // set the node's velocity to 0
    node_.value().vel = Point(0,0,0);
  }

 private:
   NODE node_;
   Point p_;
};

/* Horizontal plane constraint creating a planar obstacle */
struct PlaneConstraint {
 public:
  /* Constructor */
  PlaneConstraint(double z = -0.75):z_(z){};
  template<typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void)t;
    // iterate over all nodes
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = (*it);
      // if the node's z position is less than z_, reset to z_ and set velocityy to 0
      if(dot(n.position(), Point(0,0,1)) < z_){
        n.position().z = z_;
        n.value().vel.z = 0;
      }
    }
  }
 private:
   double z_;
};

/* Sphere constraint creating a spherical obstacle with center @c and radius @r */
struct SphereConstraint {
 public:
  /* Constructor */
  SphereConstraint(Point c, double r):c_(c), r_(r){};
  template<typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void)t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = (*it);
      // If the point is less than r_ distance from the center of the sphere
      if(norm(n.position() - c_) < r_){
        // reset the point to be on the sphere
        n.position() = nearest_point_on_sphere(n.position());
        // substract the normal component to rom the point's current velocity
        n.value().vel -= normal_component_to_sphere(n.value().vel, n.position());
      }
    }
  }
 private:
  Point c_;
  double r_;
  Point nearest_point_on_sphere(Point p){
    return c_ + (r_ / norm(p - c_)) * (p - c_);
  }
  Point normal_component_to_sphere(Point vel, Point p){
    Point R = (p - c_)/norm(p - c_);
    return dot(vel, R) * R;
  }
};

/* Sphere constraint causing all nodes in contact with the sphere to disappear */
struct SphereDisappearConstraint {
public:
 SphereDisappearConstraint(Point c, double r):c_(c), r_(r){};
 template<typename GRAPH>
 void operator()(GRAPH& g, double t){
   (void)t;
   for(auto it = g.node_begin(); it != g.node_end(); ++it){
     auto n = (*it);
     // if the node's position is less than r_ distance from the center c_, remove the node
     if(norm(n.position() - c_) < r_){
       g.remove_node(n);
     }
   }
 }
private:
 Point c_;
 double r_;
};

/** Combined Constraint functor
* @tparam C1, C2 functors with signature
*   void functor(GRAPH g, double t)
*/
template <typename C1, typename C2>
struct CombinedConstraint {
  public:
    //Constructor
    CombinedConstraint(C1 const1, C2 const2):const1_(const1), const2_(const2){};

    /** Applies both functors C1 and C2 to GRAPH g and double t
    * @tparam GRAPH g
    * @param double t
    *
    */
    template<typename GRAPH>
    void operator()(GRAPH& g, double t){
      const1_(g, t);
      const2_(g, t);
    }
  private:
    C1 const1_;
    C2 const2_;
};

/* Returns a Combined Constraint functor given two input constraints
* @tparam C1, C2 functors with signature
*   void functor(GRAPH g, double t)
*
* @return CombinedConstraint<C1, C2>
*/
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 const1, C2 const2){
  return CombinedConstraint<C1, C2>(const1, const2);
}

/** Force functor returning the gravitational force on a given Node n */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    // F = Point(0,0,m*-g)
    return n.value().mass*Point(0,0,-grav);
  }
};

/** Force functor returning the spring force on a given Node n */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    Point springForce;
    // Iterate over all edges incident to n
    for (auto adj = n.edge_begin(); adj != n.edge_end(); ++adj){
      NODE n2 = (*adj).node2();
      // calculate the scalar portion of the spring force for this adjacent node
      double scalar = -(*adj).value().K / norm(n.position() - n2.position()) * (norm(n.position() - n2.position()) - (*adj).value().L);
      // add to the spring force total
      springForce += ((n.position() - n2.position()) * scalar);
    }
    return springForce;
  }
};

/** Force functor returning the damping force on a given Node n */
struct DampingForce {
  public:
    DampingForce(double c = 0): c_(c){};
    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void)t;
      // F = -c*v
      return -c_*n.value().vel;
    }
  private:
    double c_; // damping constant
};

/** Combined Force functor
* @tparam F1, F2 functors with signature
*   Point functor(NODE n, double t)
*/
template <typename F1, typename F2>
struct CombinedForce {
  public:
    /** Constructor*/
    CombinedForce(F1 force1, F2 force2):force1_(force1), force2_(force2){};

    /** Returns a Point combining the force_1 and force2_ vectors for @n at time @t
    * @tparam NODE n
    * @param double t
    *
    * @return Point representing the combined force
    */
    template <typename NODE>
    Point operator()(NODE n, double t){
      return force1_(n,t) + force2_(n,t);
    }

  private:
    F1 force1_;
    F2 force2_;
};

/* Returns a Combined Force functor given two input forces
* @tparam F1, F2 functors with signature
*   Point functor(NODE n, double t)
*
* @return CombinedForce<F1, F2>
*/
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 force1, F2 force2){
  return CombinedForce<F1, F2>(force1, force2);
}

/* Returns a Combined Force functor given three input forces
* @tparam F1, F2, F3 functors with signature
*   Point functor(NODE n, double t)
*
* @return CombinedForce<F1, CombinedForce<F2, F3>>
*/
template<typename F1, typename F2, typename F3>
CombinedForce<F1, CombinedForce<F2, F3>> make_combined_force(F1 force1, F2 force2, F3 force3){
  return make_combined_force(force1, make_combined_force(force2, force3));
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
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for your nodes and edges
  std::vector<ConstantConstraint<Node>> constConstraints;

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    // set each node's mass to 1/N
    (*it).value().mass = (float)1/(float)graph.num_nodes();
    // v_0 = 0;
    (*it).value().vel = Point(0,0,0);

    // Add constant constraints for the points that we want to fix;
    if((*it).position() == Point(0,0,0) || (*it).position() == Point(1,0,0)){
      constConstraints.push_back(ConstantConstraint<Node>(*it));
    }

    // Rest-length should be equal to the initial length of the edge
    // All spring constants can be initialized to 100
    for (auto adj = (*it).edge_begin(); adj != (*it).edge_end(); ++adj){
      (*adj).value().L = norm((*it).position() - (*adj).node2().position());
      (*adj).value().K = 100;
    }
  }

  // Create a combined force object applying gravity, spring forces, and a damping force with constant 1/N
  auto force = make_combined_force(GravityForce(),MassSpringForce(),DampingForce((float)1/(float)graph.num_nodes()));

  // Create a combined constraint object from the two constant constraints
  auto constant_constraints = make_combined_constraint(constConstraints[0], constConstraints[1]);

  // Create a combined constraint object including the sphere disappear constraint, plane, and constant constraints
  auto constraint = make_combined_constraint(make_combined_constraint(SphereDisappearConstraint(Point(0.5,0.5,-0.5),0.15), PlaneConstraint()), constant_constraints);


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

        symp_euler_step(graph, t, dt, force, constraint);

        // Clear the viewer's nodes and Edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        //if (graph.size() < 100)
          //std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
