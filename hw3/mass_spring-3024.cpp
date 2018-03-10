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

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;          //< Edge spring constant
  double length;     //< Edge length
  EdgeData() : K(100.0), length(0.0) {}
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
 * @tparam G::node_value_type supports NodeData
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
    n.position() += n.value().vel * dt;
  }

  g = constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    /*
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      n.value().vel = Point(0);
    else
      n.value().vel += force(n, t) * (dt / n.value().mass);
    */
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
    (void) t;
    Point x_i = n.position();

    //keep 2 points static
    if (x_i == Point(0,0,0) || x_i == Point(1,0,0))
      return Point(0,0,0);

    //Spring constant: OK to have a "magic number" hardcoded here since this
    //functor is illustrative and isn't integrated with the larger program
    double K = 100;

    //add up all the spring forces
    Point force = Point(0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){

      auto edge = *it;
      auto adj_node = edge.node2();

      Point x_j = adj_node.position();
      double L = edge.value().length;
      double D = norm_2(x_i - x_j);

      Point adj_force = (-1*K * (x_i - x_j)/D * (D - L));
      force = force + adj_force;
    }
    //add gravity force
    force.z = force.z - (n.value().mass * grav);

    return force;
  }
};

//returns gravity force in the z-direction on a node
struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    Point force = Point(0);
    force.z = force.z - (n.value().mass * grav);
    return force;
  }
};

//returns the sum of all spring forces acting on a node from adjacent nodes
struct SpringForce{

  template <typename NODE>
  Point operator()(NODE n, double t){
    (void)t;
    Point x_i = n.position();
    Point force = Point(0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto edge = *it;
      auto adj_node = edge.node2();
      Point x_j = adj_node.position();
      double K = edge.value().K;
      double L = edge.value().length;
      double D = norm_2(x_i - x_j);
      Point adj_force = (-1*K * (x_i - x_j)/D * (D - L));
      force = force + adj_force;
    }
    return force;
  }
};

//returns a damping force on a node proportional to the velocity of the node
struct DampingForce{
  DampingForce(double num_nodes)
    : num_nodes(num_nodes) {
  }
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return -1 * n.value().vel/num_nodes;
  }
  private:
    double num_nodes;
};

//provides a force that does nothing to a node
//used for for function overloading
struct NoForce{
  NoForce(){}
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t; (void) n;
    return Point(0);
  }
};

//creates a combined force that returns the sum of up to 3 forces
template <typename force1, typename force2, typename force3>
struct combined_force{
  combined_force(force1 f1, force2 f2, force3 f3)
    : f1(f1), f2(f2), f3(f3) {
  }
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
  private:
    force1 f1;
    force2 f2;
    force3 f3;
};

//helper function that handles both overloading and templating for the user
template <typename f1 = NoForce, typename f2 = NoForce, typename f3 = NoForce>
combined_force<f1,f2,f3> make_combined_force(const f1& force1 = NoForce(),
                                             const f2& force2 = NoForce(),
                                             const f3& force3 = NoForce()) {
  return combined_force<f1,f2,f3>(force1, force2, force3);
}

//Constraints

//allows a user to give a vector of nodes that will constrain the node
//to have the same position over all time steps
struct ConstantConstraint{

  //create functor that knows static nodes & their original positions
  ConstantConstraint(std::vector<Node> static_input_nodes){
    for(unsigned i = 0; i < static_input_nodes.size(); i++){
      unsigned index = static_input_nodes[i].index();
      Point position = static_input_nodes[i].position();
      static_nodes.push_back({index, position});
    }
  }

  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    //keep all nodes designated as static at their original positions
    for(unsigned i = 0; i < static_nodes.size(); i++){
      Node n = g.node(static_nodes[i].index);
      n.position() = static_nodes[i].position;
    }
    return g;
  }

  private:
    struct static_node{
      unsigned index;
      Point position;
    };
    std::vector<static_node> static_nodes;
};

//a constraint that prevents nodes from passing through a plane normal to Z-axis
struct ZPlaneConstraint{
  ZPlaneConstraint(double level)
    : level(level) {
  }

  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = *it;
      //if the node passes the plane, project it back onto the plane and
      //remove its velocity in the direction normal to the plane
      if(inner_prod(n.position(), Point(0,0,1)) < level){
        n.position().z = level;
        n.value().vel.z = 0;
      }
    }
    return g;
  }

  private:
    double level;
};

//constraint that prevents a node from passing through a sphere, determined
//by a center and a radius
struct SphereConstraint{
  SphereConstraint(Point center, double radius)
    : center(center), radius(radius) {
  }

  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = *it;
      //if the node passes the sphere, project it back onto the surface
      //and remove its velocity in the direction normal to surface
      if(norm_2(n.position() - center) < radius){
        //find closest point on sphere (i.e, point from center along
        //direction of node position that is length radius away)
        Point direction_to_center = n.position() - center;
        Point unit_direction = direction_to_center/norm_2(direction_to_center);
        n.position() = center + (unit_direction * radius);

        //remove portion of velocity normal to sphere
        Point vel = n.value().vel;
        Point R = (n.position() - center)/norm_2(n.position() - center);
        n.value().vel = vel - inner_prod(vel, R) * R;
      }
    }
    return g;
  }

  private:
    Point center;
    double radius;
};

//constraint that removes nodes from the graph if they pass through the sphere
struct SphereConstraintRemove{
  SphereConstraintRemove(Point center, double radius)
    : center(center), radius(radius) {
  }

  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = *it;
      if(norm_2(n.position() - center) < radius)
        g.remove_node(n);
    }
    return g;
  }

  private:
    Point center;
    double radius;
};

//provides a constraint that does nothing to a node,
//used for for function overloading
struct NoConstr{
  NoConstr(){}
  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    return g;
  }
};

//creates a combined constraint that returns a combination of up to 3
//constraints
template <typename constr1, typename constr2, typename constr3>
struct combined_constraint{
  combined_constraint(constr1 c1, constr2 c2, constr3 c3)
    : c1(c1), c2(c2), c3(c3) {
  }
  template <typename GRAPH>
  GRAPH operator()(GRAPH g, double t){
    (void) t;
    g = c1(g, t);
    g = c2(g, t);
    g = c3(g, t);
    return g;
  }
  private:
    constr1 c1;
    constr2 c2;
    constr3 c3;
};

//helper function that handles both overloading and templating for the user
template <typename c1 = NoConstr, typename c2 = NoConstr, typename c3 = NoConstr>
combined_constraint<c1,c2,c3> make_combined_constraint(const c1& constr1 = NoConstr(),
                                                       const c2& constr2 = NoConstr(),
                                                       const c3& constr3 = NoConstr()) {
  return combined_constraint<c1,c2,c3>(constr1, constr2, constr3);
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

  // Set initial conditions for your nodes, if necessary.
  for(auto it = graph.node_begin(); !(it == graph.node_end()); ++it){
    auto node = *it;
    node.value().vel = Point(0,0,0);
    node.value().mass = 1.0/(float)graph.num_nodes();
  }

  //initial conditions for edges
  for(auto it = graph.node_begin(); !(it == graph.node_end()); ++it){
    auto node = *it;
    for(auto iit = node.edge_begin(); !(iit == node.edge_end()); ++iit){
      auto edge = *iit;
      Point pos1 = edge.node1().position();
      Point pos2 = edge.node2().position();
      edge.value().K = 100;
      edge.value().length = norm_2(pos1 - pos2);
    }
  }

  //Create constraints
  //auto no_constraint = make_combined_constraint();
  ZPlaneConstraint z_c = ZPlaneConstraint(-1.0 * 0.75);
  //SphereConstraint s_c = SphereConstraint(Point(0.5,0.5,-0.5), 0.15);
  SphereConstraintRemove s_c = SphereConstraintRemove(Point(0.5,0.5,-0.5), 0.15);

  //Create forces
  GravityForce g_f = GravityForce();
  SpringForce s_f = SpringForce();
  DampingForce d_f = DampingForce((float) graph.num_nodes());
  //auto combined_force = make_combined_force(g_f, s_f);
  auto combined_force = make_combined_force(g_f, s_f, d_f);

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
      double t_end = 2.5;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        //update constant constraint
        std::vector<Node> static_nodes;
        for(auto it = graph.node_begin(); !(it == graph.node_end()); ++it){
          auto node = *it;
          if(node.position() == Point(0,0,0) || node.position() == Point(1,0,0))
            static_nodes.push_back(node);
        }
        ConstantConstraint c_c = ConstantConstraint(static_nodes);
        auto combined_constraint = make_combined_constraint(c_c, z_c, s_c);

        //symp_euler_step(graph, t, dt, Problem1Force(), no_constraint);
        symp_euler_step(graph, t, dt, combined_force, combined_constraint);

        // Update viewer with nodes' new positions
        // Clear the viewer’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes’ new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        //graph.print_internal_structs();

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
