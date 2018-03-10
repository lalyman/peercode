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

// Damping Constant - will be initialized later;
static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Node velocity
  double L;     //< Node mass
  EdgeData() : K(0), L(0) {
  }
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

  // Constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Skip update steps for (0,0,0) and (1,0,0)
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      // Do nothing
    }
    else {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/** FORCE FUNCTORS 
 ****************************** 
 */

/** Gravity Force function object 
 *
 *  Returns the force of gravity on a node as a Point object
 */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass*Point(0,0,-grav);
  }
};

/** Mass Spring Force function object 
 *
 *  Returns the total spring force acting on a node from adjacent
 *  nodes as a Point object
 */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point f_spring = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto edge = *it;
      Point inc_node_pos = edge.node2().position();
      Point val = -edge.value().K*((n.position() - inc_node_pos)/edge.length())*
                  (edge.length() - edge.value().L);
      f_spring += val;
    }
    return f_spring;
  }
};

/** Damping Force function object 
 *
 *  Returns a damping force (form of friction) on a given node as a 
 *  Point object
 */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -c*n.value().vel;
  }
};

/** Zero Force function object - used as a default in combined force
 *
 *  Returns a zero force Point object
 */
struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    (void) n;
    return Point(0,0,0);
  }
};

/** Combined Force function object 
 *
 *
 *  F1, F2, F3 are the types of the force functors to be combined;
 *  Default type is ZeroForce
 *
 *  Private members are the input force functors
 *
 *  Returns the sum of the forces to be combined as a Point object
 */
template <typename F1 = ZeroForce, 
          typename F2 = ZeroForce, 
          typename F3 = ZeroForce>
struct CombinedForce {

  // Constructor that sets the combined force
  CombinedForce(F1 f1_, F2 f2_, F3 f3_) 
    : f1(f1_), f2(f2_), f3(f3_) {
  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t) + f3(n,t);
  }
  private:
    F1 f1;
    F2 f2;
    F3 f3;
};


/** Make Combined Force function 
 *
 *  Input parameters are the force functors to be combined; default
 *  is the ZeroForce functor
 *
 *  Constructs and initializes a CombineForce functor with these inputs
 *  as members, and returns the resulting CombinedForce functor
 */
template <typename F1 = ZeroForce, 
          typename F2 = ZeroForce,
          typename F3 = ZeroForce>
CombinedForce<F1,F2,F3> make_combined_force(F1 f1 = ZeroForce(),
                                            F2 f2 = ZeroForce(), 
                                            F3 f3 = ZeroForce()) {
  return CombinedForce<F1,F2,F3>(f1,f2,f3);
};


/** CONSTRAINT FUNCTORS 
 **************************
 */

/** Plane constraint 
 *
 *  Functor that prevents nodes from falling below z = -0.75
 *  by changing the position and velocity values of any nodes
 *  that violate the constraint
 */
struct Constraint_Plane {
  template <typename G> 
  void operator()(G& g, double t) {
    (void) t;
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node node = *nit;
      // Set position to nearest point on the plane and velocity
      // z-component to 0 if node violates constraint
      if (node.position().z < -0.75) {
        node.position().z = -0.75;
        node.value().vel.z = 0;
      }
    }
  }
};

/** Sphere constraint 1
 *
 *  Functor that prevents nodes from entering the specified sphere
 *  by changing the position and velocity values of any nodes that
 *  violate the constraints
 */
struct Constraint_Sphere {
  template <typename G>
  void operator()(G& g, double t) {
    // Sphere center and radius
    Point c(0.5,0.5,-0.5);
    double r = 0.15;
    (void) t;
    // If constraint is violated, set position to nearest point on the 
    // surface of the sphere, and component of velocity normal to 
    // sphere's surface equal to 0
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node node = *nit;
      Point dist = node.position() - c;
      Point R = dist/norm(dist);
      if (norm(dist) < r) {
        node.position() = c+r*R;
        node.value().vel = node.value().vel - (dot(node.value().vel,R))*R;
      }
    }

  }
};

/** Sphere constraint 2
 *
 *  Functor that destroys any nodes that come in contact with the
 *  given sphere by removing any nodes that violate the constraint
 */
struct Constraint_Sphere2 {
  template <typename G>
  void operator()(G& g, double t) {
    // Sphere center and radius
    Point c(0.5,0.5,-0.5);
    double r = 0.15;
    (void) t;
    // If constraint is violated, remove the node
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node node = *nit;
      Point dist = node.position() - c;
      if (norm(dist) < r) {
        g.remove_node(node);
      }
    }

  }
};

/** No Constraint - used as default for Combined Constraint
 *
 *  Functor that imposes no constraint on the graph's nodes
 */
struct No_Constraint {
  template <typename G>
  void operator()(G& g, double t) {
    (void) g;
    (void) t;
    // Do nothing
  }
};

/** Combined constraint functor 
 *
 *  Input parameters are the constraint functors to be combined; default
 *  is the No_Constraint functor
 *
 *  Constructs and initializes a Combine_Constraint functor with these inputs
 *  as members, and returns the resulting Combined_Constraint functor
 **/
template <typename C1 = No_Constraint,
          typename C2 = No_Constraint, 
          typename C3 = No_Constraint>
struct Combined_Constraint {
  
  // Constructor that initializes the members of combined constraint
  // to be the constraints we wish to combine
  Combined_Constraint(C1 c1_, C2 c2_, C3 c3_)
    : c1(c1_), c2(c2_), c3(c3_) {
  }
  template <typename G>
  void operator()(G& g, double t) {
    c1(g,t);
    c2(g,t);
    c3(g,t);
  }
  private:
    C1 c1;
    C2 c2;
    C3 c3;
};

/** Make combined constraint function
 *  
 *  Inputs are Constraint functors
 *
 *  Constructs and initializes a Combined_Constraint object with the
 *  given input constraints, and returns the resulting Combined_Constraint
 *  object 
 */

template <typename C1 = No_Constraint, 
          typename C2 = No_Constraint,
          typename C3 = No_Constraint>
Combined_Constraint<C1,C2,C3> make_combined_constraint(C1 c1 = No_Constraint(),
                                                       C2 c2 = No_Constraint(), 
                                                       C3 c3 = No_Constraint()) {
  return Combined_Constraint<C1,C2,C3>(c1,c2,c3);
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
    
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for your nodes, if necessary.

  c = 1/graph.num_nodes(); // Set damping constant

  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    auto node = *ni;

    // Set initial velocity equal to 0 and initial mass equal to
    // 1/number_of_nodes
    node.value().vel = Point(0,0,0);
    node.value().mass = 1.0/(float)graph.num_nodes();

    // Set edge spring constants and rest-lengths using incident edge
    // iterator
    for (auto it = node.edge_begin(); it != node.edge_end(); ++it) {
      auto edge = *it;
      double k = 100;
      double l = edge.length();
      edge.value().K = k; // Set spring constant to 100
      edge.value().L = l; // Set edges to have rest length equal to 
                          // initial length
    }
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
        symp_euler_step(graph, t, dt, 
                        make_combined_force(GravityForce(),
                                            MassSpringForce(),
                                            DampingForce()),
                        make_combined_constraint(Constraint_Sphere2(),
                                                 Constraint_Plane()));

        // Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions and new edges
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
