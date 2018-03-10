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
#include <unordered_map>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Those constants are useful for the first problems
static double L = 0.01;
static double K = 100;
static double c = 1.0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double initial_length;       //< Node velocity
  double K;     //< Node mass
  EdgeData() : initial_length(0), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/////////////////////////////////////////////////////////
// Forces templates

/** Force function */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if ( n . position () == Point (0 ,0 ,0) || n . position () == Point (1 ,0 ,0)){
      return Point (0 ,0 ,0);
    }
    else{
      Point force = n.value().mass * Point(0, 0, -grav);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        // For question 1 only
        // force += K * (1 - L / (*it).length()) * (((*it).node2().position()) - n.position());
        (void) t;
        // Force for question 2
        force += (*it).value().K * (1 - (*it).value().initial_length / (*it).length()) * (((*it).node2().position()) - n.position());
      }
      return force;
    }
  }
};

/** Generic template for a force aggregation.
 * ForceType1 and ForceType2 represent two abstract forces.
 * They must implement the operator()(NODE n, double t).
 * 
 * This structure can take itself as ForceType1 or 2, leading to an easy implementation
 * of make_combined_force. */ 
template <typename ForceType1, typename ForceType2>
struct GenericForce{
  GenericForce(ForceType1 force1, ForceType2 force2):
  force1(force1), force2(force2){
  }
  ForceType1 force1;
  ForceType2 force2;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return force1.operator()(n, t) + force2.operator()(n, t);
}
};

/** Create a generic force combining the effects of two forces. */ 
template < typename Force1, typename Force2 >
GenericForce<Force1, Force2> make_combined_force(Force1 force1, Force2 force2){
  GenericForce<Force1, Force2> force = GenericForce<Force1, Force2>(force1, force2);
  return force;
}
/** Create a generic force combining the effects of three forces. */
template < typename Force1, typename Force2, typename Force3>
GenericForce<Force1, GenericForce<Force2, Force3> > make_combined_force(Force1 force1, Force2 force2, Force3 force3){
  GenericForce<Force2, Force3> force_pair = GenericForce<Force2, Force3>(force2, force3);
  GenericForce<Force1, GenericForce<Force2, Force3> > force = GenericForce<Force1, GenericForce<Force2, Force3> >(force1, force_pair);
  return force;
}

struct GravityForce {
  /** Return the force applying to @a n at time @a t.
   *
   * This is applying gravity to the nodes */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point force = n.value().mass * Point(0, 0, -grav);

    return force;
  }
};

struct MassSpringForce {
  /** Return the force applying to @a n at time @a t.
   *
   * This is the elastic force applied by the edges to the nodes */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point force = Point(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      // force += K * (1 - L / (*it).length()) * (((*it).node2().position()) - n.position());
      force += (*it).value().K * (1 - (*it).value().initial_length / (*it).length()) * (((*it).node2().position()) - n.position());
    }
    return force;
  }

};

struct DampingForce {
  /** Return the force applying to @a n at time @a t.
   *
   * This is friction force due to speed*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point force = -c * n.value().vel;

    return force;
  }
};

///////////////////////////////////////////////////////////////////////
// Constraints Templates

/** Generic template for a constraints aggregation.
 * ConstraintsType1 and ConstraintsType2 represent two abstract constraints.
 * They must implement the operator()(GRAPH& g, double t, double dt).
 * 
 * This structure can take itself as ConstraintsType1 or 2, leading to an easy implementation
 * of make_combined_constr. */
template <typename ConstraintsType1, typename ConstraintsType2>
struct GenericConstraint{
  GenericConstraint(ConstraintsType1 constr1, ConstraintsType2 constr2):
  constr1(constr1), constr2(constr2){}

  ConstraintsType1 constr1;
  ConstraintsType2 constr2;
  template <typename GRAPH>
  void validate(GRAPH& g, double t, double dt) {
    constr1.validate(g, t, dt);
    constr2.validate(g, t, dt);
}
};

/** Create a generic constraint combining the effects of two constraints. */
template < typename Constraint1, typename Constraint2 >
GenericConstraint<Constraint1, Constraint2> make_combined_constr(Constraint1 constr1, Constraint2 constr2){
  GenericConstraint<Constraint1, Constraint2> constr = GenericConstraint<Constraint1, Constraint2>(constr1, constr2);
  return constr;
}

/** Create a generic constraint combining the effects of three forces. */
template < typename Constraint1, typename Constraint2, typename Constraint3>
GenericConstraint<Constraint1,
GenericConstraint<Constraint2, Constraint3> > make_combined_constr(Constraint1 constr1, Constraint2 constr2, Constraint3 constr3){
  GenericConstraint<Constraint2, Constraint3> constr_pair = GenericConstraint<Constraint2, Constraint3>(constr2, constr3);

  GenericConstraint<Constraint1,
                    GenericConstraint<Constraint2, Constraint3> > constr = GenericConstraint<Constraint1, 
                                                                  GenericConstraint<Constraint2, Constraint3> >(constr1, constr_pair);
  return constr;
}

struct FixedCornersConstraint {
  /** Returns a constraint fixing the two corners of the field
   * Checks if an update was performed that corresponded to the corners, and 
   * puts the corners back at their spots.
   * 
   * This requires the use of dt, this is why this argument is shared by the constraint aggregators
   */
  template <typename GRAPH>
  void validate(GRAPH& g, double t, double dt) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if ( (*it).position() == Point(0 ,0 ,0)  + (*it).value().vel * dt ){
        (*it).value().vel = Point(0, 0, 0);
        (*it).position() = Point(0, 0, 0);
      }
      
      else if( (*it).position() == Point(1 ,0 ,0)+ (*it).value().vel * dt ){
        (*it).value().vel = Point(0, 0, 0);
        (*it).position() = Point(1, 0, 0);
    }
    }
  }
};

struct PlaneConstraint {
  /** Creates a constraint on a plane not to cross 
   *  Crossing nodes are displaced to the border of the sphere
  */
  double limit = -0.75;

  template <typename GRAPH>
  void validate(GRAPH& g, double t, double dt) {
    (void) t;
    (void) dt;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if ( (*it).position().z < limit ){
        (*it).position().z = limit;
        (*it).value().vel = Point(0, 0, 0);
    }
    }
  }
};

struct SphereConstraint {
  /** Creates a constraint on a sphere not to cross
   * Crossing nodes are displaced to the border of the sphere
  */

  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  template <typename GRAPH>
  void validate(GRAPH& g, double t, double dt) {
    (void) t;
    (void) dt;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if ( (norm((*it).position() - center) < radius ) ){
        Point Ri = ((*it).position() - center) / norm((*it).position() - center);
        (*it).value().vel -= dot(Ri, (*it).value().vel) * Ri;
        (*it).position() = center + radius * Ri;
    }
    }
  }
};

struct SphereDestrutorConstraint {
  /** Creates a constraint on a sphere not to cross
   * Crossing nodes are destroyed. */

  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  template <typename GRAPH>
  void validate(GRAPH& g, double t, double dt) {
    (void) t;
    (void) dt;
    // The remove_node destroys a node but does not increment
    // Therefore, either we increment or we destroy, to ensure a complete visit of the nodes.
    for (auto it = g.node_begin(); it != g.node_end();){
      if ( (norm((*it).position() - center) < radius ) ){
        g.remove_node(it);
    }
    else{
      ++it;
    }
    }
  }
};
//////////////////////////////////////////////////////////////////////
// Euler updater

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
    if (!( (*it).position() == Point(0 ,0 ,0) || (*it).position() == Point(1 ,0 ,0))){
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
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
  
  // Validate the constraints

  constraint.validate(g, t, dt);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);

  }
  return t + dt;
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

  unsigned graph_size = graph.num_nodes();
  c = 1.0 / graph_size;

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    (*it).value().mass = 1.0 / graph_size;
  }
  for (auto it = graph.edge_begin(); it!= graph.edge_end(); ++it){
    (*it).value().K = 100;
    (*it).value().initial_length = (*it).length();
  }

  L = (*graph.edge_begin()).length();
  K = 100;

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

      // Create the forces and the constraints only once.
      GravityForce myGravity = GravityForce();
      MassSpringForce myMassSpringForce = MassSpringForce();
      DampingForce myDampingForce = DampingForce();

      FixedCornersConstraint myFixedCorners = FixedCornersConstraint();
      // PlaneConstraint myPlane = PlaneConstraint();
      // SphereConstraint mySphere = SphereConstraint();
      SphereDestrutorConstraint myDestructor = SphereDestrutorConstraint();
      
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());
        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce));
        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce));
        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce), myFixedCorners);
        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce),
        //                 make_combined_constr(myFixedCorners, myPlane));
        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce),
        //                 make_combined_constr(myFixedCorners, mySphere));

        // symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce),
        //                 make_combined_constr(myFixedCorners, myPlane, mySphere));

        symp_euler_step(graph, t, dt, make_combined_force(myGravity, myMassSpringForce, myDampingForce),
                        make_combined_constr(myFixedCorners, myDestructor));

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.set_label(t);

        viewer.clear ();
        node_map.clear();
        // Update viewer with nodes ’ new positions and new edges
        viewer.add_nodes(graph.node_begin() ,graph.node_end(), node_map );
        viewer.add_edges(graph.edge_begin() ,graph.edge_end(), node_map );

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
