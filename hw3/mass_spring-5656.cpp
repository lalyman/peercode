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
  double c_damp;        //< Node coefficient of damping
  NodeData() : vel(0), mass(1), c_damp(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Edge spring constant
  double L;       //< Edge spring rest-length
  EdgeData() : K(1), L(1) {}
};

// Define the Graph type
// using GraphType = Graph<NodeData>;
using GraphType = Graph<NodeData,EdgeData>;
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

// Problem 1-3
//template <typename G, typename F>
//double symp_euler_step(G& g, double t, double dt, F force) {

//Problem 4-5
template <typename G, typename F, typename C, typename P>
double symp_euler_step(G& g, double t, double dt, F force, C constraint,P ConstantPosition) {

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    ConstantPosition(n,dt);
    // if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
    //   n.position() += n.value().vel * dt;}
  }
  
  //Problem 4-5
  //Constraints to be applied
  constraint(g,t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    //Do not update if we have the point (0,0,0) or (1,0,0)
    //if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//===============================================================================================
//  FORCES
//===============================================================================================

/** Force function object for HW2 #1. */
/** Struct that defines the total force acting on a node
 *
 */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {

    Point Force_i = Point(0,0,0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Force_i;}
    //iterate through all the edges neighbor of the node @a n
    //calculate the component of spring force due to node j
    for (auto it = n.edge_begin(); it!= n.edge_end(); ++it){
      Edge edge_i = *it;
      Node node_j = edge_i.node2();

      double K_ij = edge_i.value().K;
      Point x_i = n.position();
      Point x_j = node_j.position();

      double L = edge_i.length();
      double L_ij = edge_i.value().L;
      Force_i -= K_ij*(( x_i - x_j )/L)*( L- L_ij);
    }
    double fz_grav = -grav*n.value().mass;

    (void) t;
    return Point(Force_i[0],Force_i[1],Force_i[2]+fz_grav);
  }
};

/** Struct that defines the gravity force
 *  pointing only in the z-direction
 */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //Calculate the z-component of the gravity force
    double fz_grav = -grav*n.value().mass;

    (void) t;
    return Point(0,0,fz_grav);
  }
};

/** Struct that defines the spring force
 *
 */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //initialize the spring force
    Point F_spring_i = Point(0,0,0);
    //iterate through all the edges neighbor of the node @a n
    //calculate the component of spring force due to node j
    for (auto it = n.edge_begin(); it!= n.edge_end(); ++it){
      Edge edge_i = *it;
      Node node_j = edge_i.node2();

      double K_ij = edge_i.value().K;
      Point x_i = n.position();
      Point x_j = node_j.position();

      double L = edge_i.length();
      double L_ij = edge_i.value().L;
      F_spring_i -= K_ij*(( x_i - x_j )/L)*( L- L_ij);
    }

    (void) t;
    return F_spring_i;
  }
};

/** Struct that defines the damping force
 *
 */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //mass_i = 1/N
    double c_damp = n.value().mass;

    (void) t;
    return -c_damp*n.value().vel;
  }
};


/** struct to combine forces by adding their effects.
 *
 */
template <typename f1, typename f2, typename f3>
struct combined_force{
  f1 f1_;
  f2 f2_;
  f3 f3_;
  combined_force(f1 F1,f2 F2, f3 F3):f1_(F1), f2_(F2), f3_(F3){}

  template <typename NODE>
  Point operator()(NODE n, double t){
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
      return f1_(n,t)+f2_(n,t)+f3_(n,t);} 
    else return Point(0,0,0); 
  }
};

/** method to combine three forces by adding their effects.
 *  using the object combined_force.
 *  It takes three forces and call the method operator
 */
template <typename f1, typename f2, typename f3>
combined_force<f1,f2,f3> make_combined_force(f1 F1, f2 F2, f3 F3){
  return combined_force<f1,f2,f3>(F1,F2,F3);
}



/** Struct that defines a zero Force
 *
 */
struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //force with zero intensity
    Point F_zero = Point(0,0,0);

    (void) n;
    (void) t;
    return F_zero;
  }
};

/** method to combine two forces by adding their effects.
 *  using the object combined_force.
 *  It takes two forces and a null force and call the method operator.
 *  method two sum only two forces (third force is zero)
 */
template <typename f1, typename f2>
combined_force<f1,f2,ZeroForce> make_combined_force(f1 F1, f2 F2){
  return combined_force<f1,f2,ZeroForce>(F1,F2,ZeroForce() );
}


/** 
 *  struct that update the position of a node if 
 *  that nodes is not either Point(0,0,0) or Point(1,0,0)
 */
struct ConstantPosition {

  template <typename NODE>
  void operator()(NODE n, double dt){
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
      n.position() += n.value().vel * dt;}
  }
};


//===============================================================================================
//  CONSTRAINTS
//===============================================================================================


/** struct to implement a constraint for a plane
 *  nodes cannot penetrate the the plane z = -0.75
 */
struct PlaneConstraint {
  template <typename GRAPH>

  void operator()(GRAPH& g, double t) {
    //iterate through all the nodes of the graph
    for (auto it = g.node_begin(); it!= g.node_end(); ++it){
      Node node_i = *it;
      // check if node_i violates the constraint of the plane
      if ( node_i.position()[2]< -0.75 ){
        //Set the position to the nearest point on the plane.
        node_i.position()[2] = -0.75;
        //Set the z-component of the Node velocity to zero.
        node_i.value().vel[2] = 0.0;}}

    (void) t;
  }
};

/** struct to implement a constraint for a sphere
 *  nodes cannot penetrate the a sphere with
 *  @a center = (0.5,0.5,-0.5) and @a radius = 0.15
 */
struct SphereConstraint {
  template <typename GRAPH>

  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it!= g.node_end(); ++it){
      Node node_i = *it;
      Point x_i = node_i.position();
      Point dist_from_center = x_i - Point(0.5,0.5,-0.5);
      //check if a node violates this constraint
      if ( norm( dist_from_center ) < 0.15 ){

        Point R_i = ( dist_from_center )/norm( dist_from_center );
        //Set the position to the nearest point on the surface of the sphere
        node_i.position() = Point(0.5,0.5,-0.5) + R_i*0.15;
        // Set the component of the velocity that is normal to the sphere’s surface to zero
        node_i.value().vel -= R_i*inner_prod( R_i , node_i.value().vel );}}

    (void) t;
  }
};


/** struct to implement a constraint for a sphere
 *  nodes that penetrate a sphere are removed
 */
struct SphereConstraint_remove {
  template <typename GRAPH>

  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it!= g.node_end(); ++it){
      Node node_i = *it;
      Point x_i = node_i.position();
      Point dist_from_center = x_i - Point(0.5,0.5,-0.5);
      //check if a node violates this constraint
      if ( norm( dist_from_center ) < 0.15 ){ g.remove_node(node_i);}}

    (void) t;
  }
};


/** struct to combine constraints by adding their effects.
 *
 */
template <typename c1, typename c2>
struct combined_constraint{
  c1 c1_;
  c2 c2_;
  combined_constraint(c1 cons1,c2 cons2):c1_(cons1), c2_(cons2){}

  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    c1_(g,t);
    c2_(g,t);
  }
};


/** method to combine two constraints by adding their effects.
 *  using the object combined_constraint.
 *  It takes two constraints and call the method operator
 */
template <typename c1, typename c2>
combined_constraint<c1,c2> make_combined_constraint(c1 cons1, c2 cons2){
  return combined_constraint<c1,c2>(cons1,cons2);
}

//===============================================================================================
//===============================================================================================

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
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }


  //set initial mass and coefficient of damping for node @a i and
  //set initial spring constant K and spring rest-length L
  int N = graph.num_nodes();
  for (auto it = graph.node_begin(); it!= graph.node_end(); ++it){
      Node node_i = *it;
      node_i.value().mass = 1/double(N);
      node_i.value().c_damp = 1/double(N);
      for (auto jt = node_i.edge_begin(); jt!= node_i.edge_end(); ++jt){
        Edge edge_i = *jt;
        edge_i.value().K = 100;
        edge_i.value().L = edge_i.length();
      }
  }
  //=============================

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

        //Problem 1-2
        // symp_euler_step(graph, t, dt, Problem1Force());

        //Problem 3
        // symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce()) );
        // symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce()) );

        //Problem 4-5
        symp_euler_step(graph, t, dt,
                        make_combined_force(GravityForce(),MassSpringForce(),DampingForce()),
                        make_combined_constraint(PlaneConstraint(), SphereConstraint_remove()),
                        ConstantPosition() );
                                                                                      //SphereConstraint()
                                                                                      //SphereConstraint_remove()

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
