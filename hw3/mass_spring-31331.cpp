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
#include <iostream>
#include <cmath>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
// Damping constant
static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom stucture of data to store with Edges */
struct EdgeData {
  double K; //< Spring constant K=100
  double L; //< Spring rest-length
  EdgeData() : K(100), L(0) {}
};


// Define the Graph, Node, and Edge type.
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Constraint for the plane.
* If node's position's z value <= -0.75, the node's position will be set
* to the nearest point on the plane and the z-component of the Node's 
* velocity will be set to zero.
*/
struct PlaneConstraint {
  template <typename GRAPH>
  void operator() (const GRAPH& g, const double t) const {
    (void) t;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      Node n = *ni;
      if (inner_prod(n.position(), Point(0,0,1))<-0.75) {
        n.value().vel.z = 0; //Set violated node's velocity.
        n.position().z = -0.75; //Set violated node's position.
      }
    }  
  }
};


/** Constraint for the sphere.
* If node's distance from the center of a sphere is less than the sphere's radius, 
* the node's position will be set to the nearest point on the surface of the sphere 
* and the component of the velocity normal to the surface of the sphere's surface 
* will be set to zero.
*/
struct SphereConstraint {
  template <typename GRAPH>
  void operator() (const GRAPH& g, const double t) const {
    (void) t;
      Point myCenter = Point(0.5,0.5,-0.5); //Center of the sphere
      double myRadius = 0.15; //Radius of sphere
      for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
        Node n = *ni;
        if (norm_2(n.position()- myCenter)<myRadius) {
          Point r = (n.position()-Point(0.5,0.5,-0.5))/norm_2(n.position()-Point(0.5,0.5,-0.5));
          //Set violated node's velocity.
          n.value().vel = n.value().vel - (inner_prod(n.value().vel,r))*r;
          //Set violated node's position.
          n.position() = myCenter + myRadius/norm(n.position() - myCenter) * (n.position()-myCenter);
        }
      }
      
  }
};


struct RemoveSphereConstraint {
  template <typename GRAPH>
  void operator() (GRAPH& g,double t)  {
    (void) t;
      Point myCenter = Point(0.5,0.5,-0.5); //Center of the sphere
      double myRadius = 0.15; //Radius of sphere
      for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
        Node n = *ni;
        if (norm_2(n.position()- myCenter)<myRadius) {
          g.remove_node(n); 
        }
      }
      
  }
};

/** 
Dummy object for make_combined_constraint.
*/
struct noconstraint {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) g; (void) t;
  }
};

/** Combines the imposed constraints. 
*/
template <typename constraint1 = noconstraint, typename constraint2 = noconstraint>
struct make_combined_constraint {

  constraint1 one;
  constraint2 two;
  
  make_combined_constraint(constraint1 one_ = noconstraint(), constraint2 two_ = noconstraint())
  : one(one_), two(two_) {
  }

  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    one(g,t);
    two(g,t);
  }
};


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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Constrain two corners of the cloth.
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)) {
      continue;
    }

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Adds the plane and sphere constraints.
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
    // HW2 #1: YOUR CODE HERE
    (void) t;
    
    // Initialize total force and spring force.
    Point total_force = Point(0,0,0);
    Point spring_force = Point(0,0,0);
    
    Point grav_force = (n.value().mass)*Point(0,0,-grav);
    
    // Constrain two corners of the cloth.
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)) {
      return Point(0,0,0);
    }
    
    // Update spring force by summming over all nodes adjacent to given node.
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      Point x_i = n.position();
      Point x_j = e.node2().position();
      spring_force += -e.value().K*((x_i-x_j)/norm_2(x_i-x_j))*(norm_2(x_i-x_j)-e.value().L);
    }
    total_force = spring_force+grav_force;

    return total_force;

  }
};


/** 
* Return the gravitational force applying to @a n at time @a t.
 */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point grav_force = (n.value().mass)*Point(0,0,-grav);
    return grav_force;

  }
};

/** 
* Return the spring force applying to @a n at time @a t.
 */
struct MassSpringForce {
  template <typename NODE>
    Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force = Point(0,0,0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      Point x_i = n.position();
      Point x_j = e.node2().position();
      spring_force += -e.value().K*((x_i-x_j)/norm_2(x_i-x_j))*(norm_2(x_i-x_j)-e.value().L);
    }
    return spring_force;
  }
};

/** 
* Return the damping force applying to @a n at time @a t.
 */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point damping_force = -(n.value().vel)*c;
    return damping_force;

  }
};

/** 
Dummy object for make_combined_force.
*/
struct nothing {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; (void) n;
    return Point(0);
  }
};

/** Combines the gravitational, spring, and damping forces. 
*/
template <typename force1 = nothing, typename force2 = nothing, typename force3 = nothing>
struct make_combined_force {
  force1 one;
  force2 two;
  force3 three;
  make_combined_force(force1 one = nothing(), force2 two = nothing(), force3 three = nothing())
  : one(one), two(two), three(three) {

  }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return (one(n,t) + two(n,t) + three(n,t));
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
  // Set initial conditions for your nodes, if necessary.
  
  c = 1.0/float(graph.num_nodes());

  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    Edge e = *ei;
    e.value().L = e.length();
  }

  for (auto ai : nodes) {
    ai.value().vel = Point(0);
    ai.value().mass = 1.0/float(graph.num_nodes());
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

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce, DampingForce>(GravityForce(),MassSpringForce(),DampingForce()), 
          make_combined_constraint<PlaneConstraint, RemoveSphereConstraint>(PlaneConstraint(),RemoveSphereConstraint()));

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
