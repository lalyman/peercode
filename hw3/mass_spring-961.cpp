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
/** Custome structure of data to store with Edge */
struct EdgeData{
  double K;
  double L;
  EdgeData(): K(100),L(1){}
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
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F,typename C>
double symp_euler_step(G& g, double t, double dt, F force,C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
 //   if(n.position()!=Point(0,0,0)&&n.position()!=Point(1,0,0))

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
constraint(g,t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

   if(n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
    n.value().vel=Point(0,0,0);}
   else{
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
 
   n.value().vel += force(n, t) * (dt / n.value().mass);
//     if(n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
 //   n.value().vel=Point(0,0,0);}
  }
}
//   constraint(g,t);
 

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
  Point operator()(NODE n, double t){
    // HW2 #1: YOUR CODE HERE
    if(n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
    return Point(0,0,0);
   }
  Point gravity = n.value().mass*Point(0,0,-grav);
  Point spring = Point(0,0,0);
  for(auto it = n.edge_begin();it!= n.edge_end();++it){
  auto e = *it;  
assert(e.node1() == n);
  spring +=-e.value().K*(n.position()-e.node2().position())*(e.length()-e.value().L)/e.length();
 // std::cout<<e.node2().position(); 
 }
 (void) t;
  return  gravity+spring;
   // (void) n; (void) t; (void) grav;    // silence compiler warnings
   // return Point(0);
 } 
};
struct GravityForce{
double grav_;
GravityForce(double grav=9.81):grav_(grav){}
template<typename NODE>
Point operator()(NODE n, double t){
   (void)t;
    return n.value().mass*Point(0,0,-grav_);
  }

};

struct MassSpringForce {
double K_;
MassSpringForce(double K=100):K_(K){}
template<typename NODE>
Point operator()(NODE n, double t){
  (void)t;
 Point spring = Point(0,0,0);
  for(auto it = n.edge_begin();it != n.edge_end();++it){
  auto e = *it;  
  spring +=-e.value().K*(e.node1().position()-e.node2().position())/e.length()*(e.length()-e.value().L);
}
return spring;
}
};

struct DampingForce{
double c_;
DampingForce(double c):c_(c){}
template<typename NODE>
Point operator()(NODE n, double t){
 (void)t; 
 return -n.value().vel*c_;
}
};
//add zero force to realize two force combined
struct Fake{
template<typename NODE>
Point operator()(NODE n, double t){
(void)t;
(void)n;
return Point(0,0,0);
}
};
/** return a combined Force
@Param[in] @a force1 and @a force2 @ force3 are constraints
@return make_combined force combining input forces
*/
template<typename F1,typename F2,typename F3>
struct make_combined_force{
F1 force1_;
F2 force2_;
F3 force3_;
make_combined_force(F1 force1,F2 force2,F3 force3):force1_(force1),force2_(force2),force3_(force3){}
template<typename NODE>
Point operator()(NODE n, double t){
(void)t;
return force1_(n,t)+force2_(n,t)+force3_(n,t);
}
};


struct Plane_Con{
void operator()(GraphType& g,double t){
(void)t;
for(auto it=g.node_begin();it!=g.node_end();++it){
auto n=*it;
if(n.position().z<-0.75){
n.position().z = -0.75;
n.value().vel.z = 0;
   }
  }
 }
};

struct Sphere_Con{
void operator()(GraphType& g, double t){
(void)t;
for(auto it = g.node_begin();it!=g.node_end();++it){
auto n=*it;
auto R = (n.position()-Point(0.5,0.5,-0.5))/norm(n.position()-Point(0.5,0.5,-0.5));
if(norm(n.position()-Point(0.5,0.5,-0.5))<0.15){
n.position() = Point(0.5,0.5,-0.5)+0.15*R;
n.value().vel -= dot(n.value().vel,R)*R;}
    }
  }
};

struct Fix_Con{
void operator()(GraphType& g, double t){
(void)t;
for(auto it = g.node_begin();it!=g.node_end();++it){
auto n =*it;
   if(n.position()==Point(0,0,0))
{   n.position()=Point(0,0,0);
    n.value().vel=Point(0,0,0);}
   if(n.position()==Point(1,0,0)){
   n.position()=Point(1,0,0);
   n.value().vel=Point(0,0,0);
    }
  }
 }
};

struct Remove_Con{
void operator()(GraphType& g, double t){
(void)t;
for(auto it=g.node_begin();it!=g.node_end();++it){
auto n=*it;
if(norm(n.position()-Point(0.5,0.5,-0.5))<0.15){
g.remove_node(n);
    }
   }
  return ;
   }

};

/** Functor to combine two constraints
@Param[in] constraints @a con1 and @a con2
*/
template<typename con1,typename con2>
struct constraints{
con1 c1_;
con2 c2_;
constraints(con1 c1,con2 c2):c1_(c1),c2_(c2){}
void operator()(GraphType& g,double t){
(void)t;
c1_(g,t);
c2_(g,t);
}
};

/** return a combined constraints
@Param[in] @a c1 and @a c2 are constraints
@return constraints combining input constraints
*/

template<typename con1, typename con2>
constraints<con1,con2>make_combined_constraints(con1 c1,con2 c2){return constraints<con1,con2>({c1,c2});
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
  for(auto it = graph.node_begin();it != graph.node_end();++it ){
   auto node = *it;
   node.value().vel = Point(0,0,0);
   node.value().mass = 1.0/graph.num_nodes();
  }

  for(auto it =graph.edge_begin();it != graph.edge_end();++it){
   auto edge = *it;
   edge.value().K = 100;
   edge.value().L = edge.length(); //need a function?
  // std::cout<<edge.length();
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

   Plane_Con Plane;
   Sphere_Con Sphere;
   Fix_Con Fix;
   Remove_Con Remove;
  // auto con = make_combined_constraints(Sphere,Fix);
   auto con = make_combined_constraints(Remove,Fix);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
       // symp_euler_step(graph, t, dt, Problem1Force());
//   symp_euler_step(graph,t,dt,make_combined_force<GravityForce,MassSpringForce,Fake>(GravityForce(),MassSpringForce(),Fake()));

       symp_euler_step(graph,t,dt,make_combined_force<GravityForce,MassSpringForce,DampingForce>(GravityForce(),MassSpringForce(),DampingForce(1.0/graph.num_nodes())),con);
  // symp_euler_step(graph,t,dt,make_combined_force<GravityForce,MassSpringForce,Fake>(GravityForce(),MassSpringForce(),Fake()),Plane);
 

        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(),graph.node_end(),node_map);
        viewer.add_edges(graph.edge_begin(),graph.edge_end(),node_map);

     // Update viewer with nodes' new positions
  //      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
