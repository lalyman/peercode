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
#include <cmath>
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
  double c; //Node drag coefficient
  NodeData() : vel(0), mass(1), c(0) {}
};
/** Custom structure of data to store with Edges */
struct EdgeData {
  double Kij;
  double Lij;
  EdgeData() : Kij(0), Lij(1) {}
};
// Define the Graph type
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
 * @tparam G::node_value_type supports ???????? YOU CHOOSE NodeData
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
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      continue;
   // Fixed(g,t);//make the positions of fixed points same as before
    n.position() += n.value().vel * dt;
   // Fixed(g,t);  
  }
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      continue;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
   //  Fixed(g,t);//
  }

  return t + dt;
}

struct FixedPoint {
  template <typename GraphType>
  void operator()(GraphType& g, double){
    for(auto it = g.node_begin();it!=g.node_end();++it){
      Node n =*it;
      if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        n.value().vel = Point(0,0,0);//setting the velocity to zero
    }
  }
};
/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
 // double L = 0.04; //just to check part1
  template <typename NODE>

  Point operator()(NODE n, double) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    //return Point(0);
    // Check if by adjacent nodes, they mean incident nodes
    if(n.position() == Point(0,0,0) ||
       n.position() == Point(1,0,0) )
      return Point(0,0,0);
    Point spring_force = Point(0,0,0);
    Point gravity_force;
    for(auto it = n.edge_begin();it!=n.edge_end();++it){
      Edge e = *it;
      Node adjacent;
      if(e.node1() == n)
        adjacent = e.node2();
      else
        adjacent = e.node1();
      Point distance = n.position()-adjacent.position();
      spring_force = spring_force + (-e.value().Kij*(distance)*(norm(distance)-e.value().Lij))/norm(distance);
    }
    gravity_force = n.value().mass * Point(0,0,-grav);
    return spring_force + gravity_force;
  }
};
/** This class implements the gravity force. */
struct GravityForce {
  template<typename NODE>
  Point operator()(NODE n, double){
    return  n.value().mass * Point(0,0,-grav);
  }
};

struct MassSpringForce {
  template<typename NODE>
  Point operator()(NODE n, double){
    Point spring_force = Point(0,0,0);
    for(auto it = n.edge_begin();it!=n.edge_end();++it){
      Edge e = *it;
      Node adjacent;
      if(e.node1() == n)
        adjacent = e.node2();
      else
        adjacent = e.node1();
      Point distance = n.position()-adjacent.position();
      spring_force = spring_force + (-e.value().Kij*(distance)*
                     (norm(distance)-e.value().Lij))/norm(distance);
    }
    return spring_force;
  }
};

struct DampingForce {
  template<typename NODE>
  Point operator()(NODE n, double){
    return -n.value().c * n.value().vel;
  }
};
struct PlaneConstraint {
  template<typename GraphType>
  void operator()(GraphType& graph, double){
    for ( auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      Node n = *it;
      if(dot(n.position(),Point(0,0,1)) < -0.75){
        n.position().z = -0.75; //findclosestPoint(graph,n.position());
        n.value().vel.z = 0;// set the z component of velocity to zero
      }
    }
  }
};
struct SphereConstraint {
  template<typename GraphType>
  void operator()(GraphType& graph, double){
    Point centre = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    for ( auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      Node n = *it;
      if(norm(n.position()-centre) < radius){
        Point Ri=( n.position() - centre)/norm(n.position() -centre);
        n.position() = centre + radius*Ri;
        n.value().vel = n.value().vel - dot(n.value().vel,Ri)*Ri;
      }
    }
  }
};
struct SphereConstraint2 {
  template<typename GraphType>
  void operator()(GraphType& graph, double){
    Point centre = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    for ( auto it = graph.node_begin(); it != graph.node_end();) {
      Node n = *it;
      if(norm(n.position()-centre) <= radius){
        graph.remove_node(n);
      }
      else{
        ++it;
      }
    }
  }
};

template<class Functor1, class Functor2, class Functor3> //= DampingForce()>
struct make_combined_force {
  Functor1 f1_;
  Functor2 f2_;
  Functor3 f3_;
  make_combined_force(Functor1 f1, Functor2 f2, Functor3 f3)
                      : f1_(f1), f2_(f2),f3_(f3) { }//Default constructor for f3
  template<typename NODE>
  Point operator() ( NODE n, double t) {
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
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

  // HW2 #1 YOUR CODE HERE
  double massNodes = 1.0/graph.num_nodes();
  double c = 1.0/graph.num_nodes();
  static constexpr double K = 100;
  // Set initial conditions for your nodes, if necessary.
  //std::cout<<"the mass of nodes is"<<massNodes;
  for ( auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    Node n = *it;
    n.value().mass = massNodes; //Setting the mass of each node
    n.value().vel = Point(0,0,0); //Setting the velocity of each node
    n.value().c = c; //Setting the damping coefficient
  }
  for( auto it = graph.edge_begin(); it!= graph.edge_end(); ++it) {
      Edge e = *it;
      e.value().Kij = K;
      e.value().Lij = e.length();
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
      double t_end = 1.0;
      //PlaneConstraint p;
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
       // symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt,Problem1Force(),SphereConstraint2());
 
     //   make_combined_force<GravityForce,MassSpringForce,DampingForce>
      //                      (GravityForce(),MassSpringForce(),DampingForce()),
        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
     /*   if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));*/
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
