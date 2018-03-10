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

struct EdgeData {
  double L;
  double K;
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using NodeType = typename GraphType::node_type;
using EdgeType = typename GraphType::edge_type;


/** Functor returning a Null Force */
struct NullForce {
  /** Return null force. */
  Point operator()(NodeType, double) {
    return Point(0,0,0);
  }
};

/** Functor for a Null Constraint */
struct NullConstraint {
    /** Enforce null constraint. */
    void operator()(GraphType&, double) {return;}
};


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F = NullForce, typename C = NullConstraint>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Enforce constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);

    // Ad-hoc constraint for HW2 #3
    /*if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      n.value().vel = Point(0,0,0); */ // Constrain two corners
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

  Point operator()(NodeType n, double) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0); // Constrain two corners

    Point fSpring = Point(0,0,0);
    Point x_i = n.position();
    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit){
      Point displ = x_i - (*eit).node2().position();
      fSpring += -(*eit).value().K * (norm(displ) - (*eit).value().L)
                              / norm(displ) * displ;
    }
    return fSpring + n.value().mass * Point(0,0,-grav);
  }
};

/** Functor for only Gravity Force */
struct GravityForce {
  /** Return gravitational force. */
  Point operator()(NodeType n, double) {
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Functor for only Mass Spring Force */
struct MassSpringForce {
  /** Return spring force. */
  Point operator()(NodeType n, double) {
    Point fSpring = Point(0,0,0);
    Point x_i = n.position();
    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit){
      Point displ = x_i - (*eit).node2().position();
      fSpring += -(*eit).value().K * (norm(displ) - (*eit).value().L)
                              / norm(displ) * displ;
    }
    return fSpring;
  }
};

/** Functor for only Damping Force */
struct DampingForce {
  double c_;
  DampingForce(const double c) : c_(c) {};
  /** Return damping force. */
  Point operator()(NodeType n, double) {
    return -c_ * n.value().vel;
  }
};

/** Combines multiple force functors */
template <typename F1, typename F2 = NullForce, typename F3 = NullForce>
struct make_combined_force {

  F1 f1_;
  F2 f2_;
  F3 f3_; // Will use NullForce() as default initialisation

  make_combined_force(F1 f1) : f1_(f1) {};
  make_combined_force(F1 f1, F2 f2) : f1_(f1), f2_(f2) {};
  make_combined_force(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {};

  /** Return sum of Forces */
  Point operator()(NodeType n, double t) {
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }
};

/** Find nearest node to point */
NodeType nearest_node(const GraphType& g, const Point& point) {
  return *std::min_element(g.node_begin(), g.node_end(),
           [&](NodeType a, NodeType b){
                return norm(a.position() - point) < norm(b.position() - point);
               });
}

/** Functor enforcing Point Constraint
 * @param[in] holds A vector of Points to hold still
 */
struct PointConstraint {

  std::vector<NodeType> holds_;
  std::vector<Point> points_;

  PointConstraint(const GraphType& g, std::vector<Point> points)
      : points_(points){
    for (auto pit = points.begin(); pit != points.end(); ++pit){
      holds_.push_back(nearest_node(g,*pit));
    }
  };

  /** Enforces Point Constraint
   * @param[in] g Graph to enforce
   * @post For all nodes such that position == holds_[i] for some i,
   *          position remains at holds_[i], and velocity = 0.
   *          I.e. Enforces no-slip internal BC.
   */
  void operator()(GraphType&, double){
    for (unsigned long i = 0; i < holds_.size(); ++i){
      // Enforce all the hold points
      holds_[i].position() = points_[i];
      holds_[i].value().vel = Point(0,0,0);
    }
  }
};

/** Functor enforcing Plane Constraint
 * @param[in] normal The normal vector, i.e. Point defining plane
 * @param[in] offset Normal distance from origin
 */
struct PlaneConstraint {

  Point normal_;
  double offset_;

  PlaneConstraint(Point normal, double offset) : offset_(offset) {
    // Enforce that normal is the unit normal
    double length = norm(normal);
    normal_ = normal/length;
  };

  /** Enforces Plane Constraint
   * @param[in] g Graph to enforce
   * @post For all nodes such that position \cdot normal < offset,
   *          position is at projection onto plane, and velocity.z = 0.
   *          I.e. Enforces no-slip internal BC.
   */
  void operator()(GraphType& g, double){
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit){
      auto n = *nit;

      if (dot(n.position(), normal_) < offset_){ // Crossed the plane!
        // Project back onto plane
        Point position_rel = n.position() - Point(0,0,offset_/normal_.z);
        double depth = dot(position_rel, normal_);
        n.position() -= depth * normal_;
        // Zero normal velocity -- No-slip
        n.value().vel = n.value().vel - dot(n.value().vel, normal_) * normal_;
      }
    }
  }
};

/** Functor enforcing Sphere Constraint
 * @param[in] centre The Point defining sphere centre
 * @param[in] radius The radius of the sphere
 */
struct SphereConstraint {

  Point centre_;
  double radius_;

  SphereConstraint(Point centre, double radius)
      : centre_(centre), radius_(radius) {};

  /** Enforces Sphere Constraint
   * @param[in] g Graph to enforce
   * @post For all nodes such that |position - centre| < radius,
   *          position is at projection onto sphere, and normal velocity = 0.
   *          I.e. Enforces no-slip internal BC.
   */
  void operator()(GraphType& g, double){
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit){
      auto n = *nit;

      if (norm(n.position() -  centre_) < radius_){ // Inside the sphere!
        // Project back onto sphere
        Point position_rel = n.position() - centre_;
        Point R = position_rel / norm(position_rel); // Unit radius vector
        double depth = radius_ - norm(position_rel);
        n.position() += depth * R;
        // Zero normal velocity -- No-slip
        n.value().vel = n.value().vel - dot(n.value().vel, R) * R;
      }
    }
  }
};

/** Functor destroying nodes violating Sphere Constraint
 * @param[in] centre The Point defining sphere centre
 * @param[in] radius The radius of the sphere
 */
struct SphereDestroyConstraint {

  Point centre_;
  double radius_;

  SphereDestroyConstraint(Point centre, double radius)
      : centre_(centre), radius_(radius) {};

  /** Deletes nodes/edges violating Sphere Constraint
   * @param[in] g Graph to enforce
   * @post All nodes such that |position - centre| < radius are deleted.
   */
  void operator()(GraphType& g, double){
    auto nit = g.node_begin();
    while (nit != g.node_end()){
      if (norm((*nit).position() -  centre_) < radius_) // Inside the sphere!
        nit = g.remove_node(nit);
      else
        ++nit;
    }
  }
};

/** Combines multiple constraint functors */
template <typename C1, typename C2 = NullConstraint, typename C3 = NullConstraint>
struct make_combined_constraint {

  C1 c1_;
  C2 c2_;
  C3 c3_; // Will use NullForce() as default initialisation

  make_combined_constraint(C1 c1) : c1_(c1) {};
  make_combined_constraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {};
  make_combined_constraint(C1 c1, C2 c2, C3 c3) : c1_(c1), c2_(c2), c3_(c3) {};

  /** Return sum of Constraints */
  void operator()(GraphType& g, double t) {
    c1_(g,t);
    c2_(g,t);
    c3_(g,t);
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
  std::vector<NodeType> nodes;
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

  // Initialise Mass and Velocity
  for(auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
    (*nit).value().vel = Point(0, 0, 0);
    (*nit).value().mass = 1.0 / graph.num_nodes();
  }

  // Initialise Edge Values
  for(auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit){
    (*eit).value().K = 100.;
    (*eit).value().L = (*eit).length();
  }

  double c_damp = 1.0/graph.num_nodes(); // Damping Coefficient

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

      // HW2 #1
      //symp_euler_step(graph, t, dt, Problem1Force());

      // HW2 #3.0
      /**symp_euler_step(graph, t, dt,
          make_combined_force<GravityForce,MassSpringForce>
                             (GravityForce(),MassSpringForce())); */
      // HW2 #3.1
      /**symp_euler_step(graph, t, dt,
        make_combined_force<GravityForce,MassSpringForce,DampingForce>
                           (GravityForce(),
                            MassSpringForce(),
                            DampingForce(c_damp))); */

      // HW2 #4.0
      /* symp_euler_step(graph, t, dt,
        make_combined_force<GravityForce,MassSpringForce,DampingForce>
                           (GravityForce(),
                            MassSpringForce(),
                            DampingForce(c_damp)),
        make_combined_constraint<PointConstraint,PlaneConstraint>
                           (PointConstraint(graph,{Point(0,0,0),Point(1,0,0)}),
                            PlaneConstraint(Point(0,0,1),-0.75))); */

      // HW2 #4.1
      /* symp_euler_step(graph, t, dt,
        make_combined_force<GravityForce,MassSpringForce,DampingForce>
                           (GravityForce(),
                            MassSpringForce(),
                            DampingForce(c_damp)),
        make_combined_constraint<PointConstraint,SphereConstraint>
                           (PointConstraint(graph,{Point(0,0,0),Point(1,0,0)}),
                            SphereConstraint(Point(0.5,0.5,-0.5),0.15))); */

        // HW2 #5
        symp_euler_step(graph, t, dt,
          make_combined_force<GravityForce,MassSpringForce,DampingForce>
                          (GravityForce(),
                           MassSpringForce(),
                           DampingForce(c_damp)),
          make_combined_constraint<PointConstraint,SphereDestroyConstraint>
                          (PointConstraint(graph,{Point(0,0,0),Point(1,0,0)}),
                           SphereDestroyConstraint(Point(0.5,0.5,-0.5),0.15)));

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
