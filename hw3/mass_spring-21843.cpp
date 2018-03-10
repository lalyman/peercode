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
  double K; // Constant
  double L; // initial rest length
  EdgeData() : K(0), L(0) {}
};
// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = typename GraphType::size_type;

/** Add a constant constraint to plane.
 *  Skip point (0,0,0) and (1,0,0).
 *  Complexity: O(1)
 */
struct ConstantConstraint {

  // assign value
  ConstantConstraint(size_type i1, size_type i2): i1_(i1),i2_(i2){}
  

  void operator()(GraphType& graph, double t){
    // set its velocity to be 0
    graph.node(i1_).position() = Point(0,0,0);
    graph.node(i1_).value().vel = Point(0,0,0);
    graph.node(i2_).position() = Point(1,0,0);
    graph.node(i2_).value().vel = Point(0,0,0);        
    (void) t;
  }

  size_type i1_;
  size_type i2_;
};

/** Add a constraint to plane
 *  @param[in,out] g  Graph
 *  @param[in]     t  Time

 *  @post If xi*(0, 0, 1) < −0.75, xi's z=0.75, velocity's z=0.
 */
struct PlaneConstraint {
  void operator()(GraphType& graph, double t){
    for (auto i=graph.node_begin();i!=graph.node_end();++i){
      // check if it violates the condition
      if (dot((*i).position(),Point(0,0,1))<-0.75){
        (*i).position().elem[2] = -0.75;
        (*i).value().vel.elem[2] = 0.0;
      }
    }
    (void) t;
  }
};

/** Add a constraint to sphere
 *  @param[in,out] g  Graph
 *  @param[in]     t  Time
 *  @param[in]     c  Center point
 *  @param[in]     r  Radius

 *  @post If |xi − @a c|<@a r, the position is the nearest point on the sphere, 
 *        the component of the velocity that is normal to the sphere’s surface to zero.
 */
struct SphereConstraint {
  Point c = Point(0.5,0.5,-0.5);
  double r = 0.15;
  void operator()(GraphType& graph, double t){
    Point d;
    for (auto i=graph.node_begin();i!=graph.node_end();++i){
      d = (*i).position()-c;
      // check if it violates the condition
      if (norm(d)<r){
        (*i).position() = c+d*r/norm(d);
        (*i).value().vel -= dot((*i).value().vel, (d/norm(d)))*(d/norm(d));
      }
    }
    (void) t;    
  }
};

struct RemoveConstraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& graph, double t){
    auto i = graph.node_begin();
    while(i != graph.node_end()){
      // check if it meets the condition
      if(norm((*i).position() - c) < r){
        i = graph.remove_node(i);
      }
      else{
        ++i;
      }
    }
    (void) t; 
  }  
}; 

/** Combine constraints. */
template<typename C1, typename C2>
struct CombCon{
  C1 c1;
  C2 c2;
  CombCon(C1 con1=C1(), C2 con2 = C2()):c1(con1),c2(con2){}

  void operator()(GraphType& graph, double t){
    //(void) t;
    c1(graph, t);
    c2(graph, t);
  }
};

/** Combine two constraints. */
template<typename C1, typename C2>
CombCon<C1, C2> make_combined_constraints(C1 c1 = C1(), C2 c2 = C2()){
  return CombCon<C1, C2>(c1, c2);
}

/** Combine three constraints. */
template<typename C1, typename C2, typename C3>
CombCon<CombCon<C1, C2>,C3> make_combined_constraints(C1 c1, C2 c2, C3 c3){
  return make_combined_constraints(make_combined_constraints(c1,c2),c3);
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  //auto constraints = make_combined_constraints(ConstantConstraint(),PlaneConstraint(),SphereConstraint());
  //auto constraints = make_combined_constraints(ConstantConstraint(),PlaneConstraint(),RemoveConstraint());
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.position() += n.value().vel * dt;

    /** Problem 2 part */
    // Update the position of the node according to its 
    // if (n.position() != Point(0,0,0) && n.position() != Point(1, 0, 0)){
    //   // x^{n+1} = x^{n} + v^{n} * dt
    //   n.position() += n.value().vel * dt;
    // }
  }

  constraints(g,t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    /** Problem 2 part */
    // // skip point at (0,0,0) and (1,0,0)
    // if(n.position() == Point(0,0,0) || n.position() == Point(1, 0, 0)){
    //   n.value().vel = Point(0, 0, 0);
    // }


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
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)){
      return Point(0,0,0);
    }
    
    Point Spring = Point(0,0,0);
    Point Grav = n.value().mass * Point(0,0,-grav);
    
    for (auto i=n.edge_begin(); i!=n.edge_end();++i){
      /** Problem 1 part */
      //Spring += ((-100)*((n.position() - (*i).node2().position())/
                //(*i).length())*((*i).length()- 0.25));

      // calculate spring force
      Spring += ((-1*(*i).value().K)*((n.position() - (*i).node2().position())/
                (*i).length())*((*i).length()- (*i).value().L));
    }
    (void) t;
    return Spring+Grav;
  }
};

/** Return the gravity force. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Return the spring force. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point Spring = Point(0,0,0);
    for (auto i=n.edge_begin(); i!=n.edge_end();++i){
      // sum up spring
      Spring += ((-1*(*i).value().K)*((n.position() - (*i).node2().position())/
                (*i).length())*((*i).length()- (*i).value().L));
    }
    (void) t;
    return Spring; 
  }
};

/** Return the damping force. */
struct DampingForce {
  double c_;
  DampingForce(double c):c_(c){}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -1*c_*n.value().vel;

  }
};   


/** Struct for combined force, using for the later function. */
template<typename F1, typename F2>
struct CombFor{
  F1 f1;
  F2 f2;
  CombFor(F1 force1 = F1(), F2 force2 = F2()):f1(force1),f2(force2){}

  template <typename NODE>
  Point operator()(NODE n, double t){
    //(void) t;
    return (f1(n, t) + f2(n, t));
  }
};

/** Calculate two forces addition. */
template<typename F1, typename F2>
CombFor<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2 = F2()){
  return CombFor<F1, F2>(f1, f2);
}

/** Calculate three forces addition. */
template<typename F1, typename F2, typename F3>
CombFor<CombFor<F1,F2>,F3> make_combined_force (F1 f1, F2 f2, F3 f3){
  return make_combined_force(make_combined_force(f1, f2),f3);
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

  // initialize mass and vel
  for (auto i = graph.node_begin();i!=graph.node_end();++i){
    (*i).value().mass = 1.0/graph.num_nodes();
    (*i).value().vel = Point(0,0,0);
  }

  //initialize K and L
  for (auto j=graph.edge_begin();j!=graph.edge_end();++j){
    (*j).value().K = 100.0;
    (*j).value().L = (*j).length();
  }


  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  // update with new position
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
      double c = 1.0/graph.num_nodes();


      size_type i1,i2;



      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {

          // check if points need to be skipped
          for (auto i=graph.node_begin();i!=graph.node_end();++i){
            auto n = *i;
            if (n.position()==Point(0,0,0)){
              i1 = n.index();
            }
            if (n.position()==Point(1,0,0)){
              i2 = n.index();
            }        
          }
        
        // cal symp_euler_step
        //symp_euler_step(graph, t, dt, Problem1Force());
        // symp_euler_step(graph, t, dt, 
        //   make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)),
        //   make_combined_constraints(ConstantConstraint(i1,i2),PlaneConstraint(),SphereConstraint()));
        symp_euler_step(graph, t, dt, 
          make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)),
          make_combined_constraints(ConstantConstraint(i1,i2),PlaneConstraint(),RemoveConstraint()));

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
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
