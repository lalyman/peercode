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
  double K;       //< Spring constant
  double L;     //< Rest length
  EdgeData(double ini_K,double ini_L) {
    K = ini_K;
    L = ini_L;
  }
  EdgeData () {} // Default constructor
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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraint to each timestep
  constraint(g);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);

    // Avoid update the time step at these two points
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      n.value().vel = Point(0,0,0);
    }
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

    // Initialize spring force
    Point mass_spring_force = Point(0,0,0);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {

      // Prevent falling to infinity
      if (n.position() == Point(0,0,0) or n.position() == Point(1,0,0)) {
        return Point(0,0,0);
      }

      auto curr_e = *it;

      // Calculating spring fore
      Point node_diff = n.position() - curr_e.node2().position();
      mass_spring_force += -100*node_diff/norm(node_diff)*(norm(node_diff)-curr_e.value().L);
    }

    (void) t;    // silence compiler warnings
    return mass_spring_force + n.value().mass*Point(0,0,-grav);
  }
};

/** Gravity Force object */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;  
    return (n.value().mass*Point(0, 0, -grav));
  }
};

/** Mass Spring Force object */
struct MassSpringForce {
  /** Return the mass spring force applying to @a n at time @a t. 
  *
  * In this object. the value of spring constant @a K and rest length
  * @a L are not constant.
  *
  */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point mass_spring_force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto curr_e = *it;
      double K = curr_e.value().K;
      double L = curr_e.value().L;
      Point node_diff = n.position() - curr_e.node2().position();
      mass_spring_force += -K * (node_diff/norm(node_diff)) * (norm(node_diff)-L);
    }
    (void) t;   // silence compiler warnings
    return mass_spring_force;    
  }
};

/** Damping Force object */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t. */
  double c;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    (void) n;
    return -c * n.value().vel;
  }
};

/** Force addition object */
template <typename force1, typename force2>
struct force_addition {
  /** Return the summation of two forces applying to @a n at time @a t. */
  force1 f1;
  force2 f2;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t);
  }
};

/** Function that is used to construct two force addition
 * @param[in]     f1     First type of force
 * @param[in]     f2     Second type of force
 *
 * @return  The force_addition object that use f1 and f2 as the variables
 */
template <typename force1, typename force2>
force_addition<force1, force2> make_combined_force(force1 f1, force2 f2) {
  return {f1,f2};
}

/** Function that is used to construct two force addition
 * @param[in]     f1     First type of force
 * @param[in]     f2     Second type of force
 * @param[in]     f3     Third type of force
 *
 * @return  The force_addition object that use combined force(f1,f2) and f3 as the variables
 */
template <typename force1, typename force2, typename force3>
force_addition<force_addition<force1,force2>, force3> make_combined_force(force1 f1, force2 f2, force3 f3) {
  // Combine f1 and f2 to one force
  force_addition<force1, force2> f4 = make_combined_force(f1, f2);
  return {f4,f3};
}

/** Plane Constraint Object */
struct plane_con {
  void operator()(GraphType& graph) {
    Point fixed_p = Point(0,0,1);
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      auto curr_node = *it;
      // Check for constraint violation
      if (dot(curr_node.position(),fixed_p) < -0.75) {
        curr_node.value().vel.z = 0.0;
        curr_node.position().z = -0.75;
      }
    }
  }
};

/** Sphere Constraint 1 Object */
struct sphere_con1 {
	Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& graph) {
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      auto curr_node = *it;
      Point diff = curr_node.position() - c;
      Point Ri = diff / norm(diff);
      if (norm(diff) < r) {
        // Component of the velocity that is normal to the sphere's surface to zero
        curr_node.value().vel = curr_node.value().vel - dot(curr_node.value().vel, Ri)*Ri;//dot(curr_node.value().vel,Ri)*Ri;
        // Nearest point on the surface
        curr_node.position() = Ri*r;
      }
    }
  }
};

/** Sphere Constraint 2 Object */
struct sphere_con2 {
	Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& graph) {
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      auto curr_node = *it;
      Point diff = curr_node.position() - c;
      if (norm(diff) < r) {
        graph.remove_node(curr_node);
      }
    }
  }
};

/** Object that combines all the constraints together */
template <typename con1, typename con2, typename con3>
struct con_comb {
  con1 c1;
  con2 c2;
  con3 c3;
  void operator()(GraphType& graph) {
    c1(graph);
    c2(graph);
    c3(graph);
  }
};

/** Function that is used to combine all constraints
 * @param[in]     c1     First type of constraint
 * @param[in]     c2     Second type of constraint
 * @param[in]     c3     Third type of constraint
 *
 * @return  The con_comb object that use c1,c2 and c3 as the variables
 */
template <typename con1, typename con2, typename con3>
con_comb<con1,con2,con3> make_combined_constraint(con1 c1, con2 c2, con3 c3) {
  return {c1,c2,c3};
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
  // Set initial conditions for node nodes
  double c = 1.0 / graph.num_nodes();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto curr_n = *it;
    curr_n.value().vel = Point(0,0,0);
    curr_n.value().mass = 1.0 / graph.num_nodes();
  }

  // Initialize the value for EdgeData
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto curr_edge = *it;
    curr_edge.value() = EdgeData(100.0, curr_edge.length());
  }

  // Initialize all the constraints and combine them into one variable
  plane_con c1;
  sphere_con1 c2;
  sphere_con2 c3;
  auto all_cons = make_combined_constraint(c1, c3, c2);

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
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce{c}));
        // symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), 
        //   MassSpringForce(),DampingForce{c}),c2);
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), 
          MassSpringForce(),DampingForce{c}),all_cons);
        // symp_euler_step(graph, t, dt, Problem1Force(),all_cons);

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
