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
static double c; 

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom strcuture of data to store with Edges */
struct EdgeData {
  double K;      // Spring Constant
  double L;      // Length of edge
  EdgeData(): K(100), L(0){} 
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

//Remove Sphere Constraint removes all nodes that are within a sphere with center 
//(0.5,0.5,-0.5) and radius 0.15. However, since I did not get my remove_edges to 
//work properly, it is commented out

struct RemoveSphereConstraint {
  template <typename GRAPH> 
  void operator()(GRAPH& g, double t) {
    (void) t; 
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15; 
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      Node n = *ni;
      if (norm(n.position() - center) < radius) {
        g.remove_node(n); 
      }
    }
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

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt

    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)) {
      continue; 
    }
    n.position() += n.value().vel * dt;
  }

  // Calls the combined Plane and Sphere constraints which will 
  // manually set the velocities and positions. 
  constraint(g,t);


  /** Call to remove my nodes and edges if it had worked correctly. 
  RemoveSphereConstraint r; 
  r(g,t); */ 

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Sets the gravity forces for each node. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point gravity_force = n.value().mass * Point(0,0,-grav); 
    return (gravity_force);
  }
};

/** Sets the spring forces for each node */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force = Point(0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto def_ei = *ei;
      Point x_i = n.position(); 
      Point x_j = def_ei.node2().position(); 
      spring_force += -def_ei.value().K * (x_i - x_j)/ norm(x_i - x_j) * 
      (norm(x_i - x_j) - def_ei.value().L);
    }
    return(spring_force);
  }
};

/** Sets the damping force for each node. */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point damping_force = -(n.value().vel) * c;  
    return (damping_force);
  }
};

/** Empty force struct to be used in the combined forces struct below  */
struct empty_force {
  template <typename NODE>
  Point operator() (NODE n, double t) {
    (void) t; (void) n; 
    return Point(0); 
  }
}; 

/** This combines the gravity force, spring force, and damping force and allows 
a variation of them to be called at one time.*/
template <typename force1 = empty_force, typename force2 = empty_force, 
typename force3= empty_force>
struct make_combined_force {
  force1 one; 
  force2 two; 
  force3 three; 

  //Intializes them to the empty_force struct incase they are not called. 
  make_combined_force(force1 one_ = empty_force(), force2 two_ = empty_force(), 
    force3 three_ = empty_force()): 
    one(one_),two(two_), three(three_) {}
  template <typename NODE>
  Point operator() (NODE n, double t) {
    (void) t;
    return one(n,t) + two(n,t) + three(n,t); //Add the forces together. 
  }
};

/** Struct from Problem 1 that calculates all the forces together. This does not 
allow for flexibility or generaliziations. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force = Point(0);
    Point gravity_force = n.value().mass * Point(0,0,-grav);  // Calculates the gravity force.
    //Holds the two points constant. 
    if (n.position() == Point(0,0,0)|| n. position() == Point(1,0,0)) {
      return(Point(0,0,0)); 
    }
    //Calculates the spring forces.
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto def_ei = *ei;
      Point x_i = n.position(); 
      Point x_j = def_ei.node2().position(); 
      spring_force += -def_ei.value().K * (x_i - x_j)/ norm(x_i - x_j) * 
      (norm(x_i - x_j) - def_ei.value().L);
    }
    Point total_force = gravity_force + spring_force; //Adds the forces together.
    return(total_force);

  }
};

/** Plane constraint that sets the velocity and position of all nodes that 
are less that -0.75 in their z values.*/
struct PlaneConstraint {
  template <typename GRAPH> 
  void operator()(GRAPH& g, double t){
    (void) t; 
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      Node n = *ni;
      if (n.position().z < -0.75) {
      n.position().z = -0.75;
      n.value().vel.z = 0; 
      } 
    }
  } 
};

/** Sphere Constraint that sets the velocity and position of all nodes
that are within a sphere with center (0.5,0.5,-0.5) and radius 0.15. */
struct SphereConstraint {
  template <typename GRAPH> 
  void operator()(GRAPH& g, double t) {
    (void) t; 
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15; 
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      Node n = *ni;
      if (norm(n.position() - center) < radius) {
        n.position() = center + radius/norm(n.position() - center) * (n.position() - center); 
        Point ri = (n.position() - center)/norm(n.position() - center);
        n.value().vel = n.value().vel - (inner_prod(n.value().vel,ri)) *ri; 
      }
    }
  }
};

/** Empty constraint struct to be used in the combined forces struct below  */
struct empty_constraint {
  template <typename GRAPH>
  void operator() (GRAPH& g, double t) {
    (void) g; 
    (void) t; 
  }
}; 

/** This combines the contraints and allows a variation of them to be called at one time.*/
template <typename constraint1 = empty_constraint, typename constraint2 = empty_constraint>
struct make_combined_constraint {
  constraint1 one; 
  constraint2 two; 

  //Intializes them to the empty_constraint struct incase they are not called. 
  make_combined_constraint(constraint1 one_ = empty_constraint(), 
    constraint2 two_ = empty_constraint()):
  one(one_), two(two_) {}

  template <typename GRAPH>
  void operator()(GRAPH& g, double t) { 
    // Each constraint is called. 
    one(g,t);
    two(g,t);
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

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Initializes the length value for each edge.
  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    Edge def_ei = *ei;
    def_ei.value().L = def_ei.length(); 
  }

  // Initializes the mass and velocity values for each node. 
  for (auto n: nodes) {
    n.value().mass = 1.0/float(graph.num_nodes()); 
    n.value().vel = Point(0);
  }

  // Calculates the damping constant. 
  c = 1.0/float(graph.num_nodes()); 

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
        symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce,DampingForce>(GravityForce(), MassSpringForce(),DampingForce()), 
          make_combined_constraint<PlaneConstraint,RemoveSphereConstraint>(PlaneConstraint(),RemoveSphereConstraint()));

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