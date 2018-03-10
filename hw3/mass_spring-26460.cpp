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
#include <math.h>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Spring constant
static constexpr double spring_constant = 100.0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Constraint function object defined by a plane. */
struct PlaneConstraint {

  /** If node violates contraint of plane, reset node. */
  void operator()(GraphType g, double t) {
    double z_coord = -0.75;

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point n_pos = n.position();
     
      if (dot(n_pos, Point(0,0,1)) < z_coord) {
        Point normal = Point(0, 0, 1); // Normal vector to plane.
        double v = (z_coord - dot(n_pos, normal)) / (dot(normal, normal));
        Point nearest = (n_pos + v) * normal;
        n.position() = nearest;  // Reset position to nearest point on plane.
        n.value().vel.z = 0;  // Reset z--coord of velocity to zero.
      }
    }
  }
};

/** Constraint function object defined by a sphere. */
struct Sphere1Constraint {

  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  /** If node violates contraint of plane, reset node. */
  void operator()(GraphType g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point n_pos = n.position();

      if (norm(n_pos - center) < radius) {
        // Closest point on sphere to a given point p is:
        Point nearest = center + ((radius/(norm(n_pos - center))) *
            (n_pos - center));
        // Reset position to nearest point on plane.
        n.position() = nearest;
        Point R_i = (n_pos - center) / norm(n_pos - center);
        // Reset component of velocity normal to sphere to zero.
        n.value().vel = n.value().vel - ((dot(n_pos, R_i)) * R_i);
      }
    }
  }
};

/** Constraint function object defined by sphere.
 *  Instead of node reset, uses node removal.
 */
struct Sphere2Constraint {

  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  /** If node violates contraint of plane, reset node. */
  void operator()(GraphType g, double t) {
    auto ni = g.node_begin();
    Node n = *ni;

    if (norm(n.position() - center) < radius) {
      Point::size_type result = g.remove_node(n); 
    }
  }
};

/** Generic function object that represents two combined constraints. */
template <class C1, class C2>
struct CombinedConstraints {

  void operator()(GraphType g, double t) {
    constraint1(g, t);
    constraint2(g, t);
  }

  C1 constraint1;
  C2 constraint2;
};

/** Method to combine two constraints. */
template <class C1, class C2>
CombinedConstraints<C1, C2> combine_constraints(C1 c1, C2 c2) {
  CombinedConstraints<C1, C2> combined_constraints;
  combined_constraints.constraint1 = c1; 
  combined_constraints.constraint2 = c2; 
  return combined_constraints;
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // HW2 P3: Skip update for nodes (0,0,0) and (1,0,0).
    if ((n.position() != Point(0,0,0)) && (n.position() != Point(1,0,0))) {
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  // Uncomment to constrain the simulation once edge and node removal fixed..

  //CombinedConstraints<PlaneConstraint, Sphere2Constraint> combined_constraints =
  //  combine_constraints(PlaneConstraint(), Sphere2Constraint());
  //combined_constraints(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Hw2 P3: Skip update for nodes (0,0,0) and (1,0,0).    
    if ((n.position() != Point(0,0,0)) && (n.position() != Point(1,0,0))) {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
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

    // To prevent cloth from falling to infinity, constrain two corners of
    // cloth by returning a zero force.
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    int num_nodes = n.graph()->size();
   
    // Set mass to 1/N where N is # of nodes in graph (constant density).
    n.value().mass = 1.0 / double(num_nodes); 
    
    // Initialize all edges to have rest length equal to initial edge lenth. 
    const double L = 1.0 / (sqrt(num_nodes)-1);

    // Spring force.
    Point f_spring = Point(0, 0, 0); // Initialize to zero vector.

    // Update spring force on node according to values of adjacent nodes.
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      assert(e.node1() == n);
   
      double edge_len = e.length();  
      double K = e.value();  // Edge-specific spring constant K_ij.    
      Point adj_pos = e.node2().position(); // Position of adjacent node.
      Point displacement = (n.position() - adj_pos) / edge_len;
      
      f_spring += (displacement) * (-K * (edge_len - L));
    }

    // Gravity force.
    Point f_gravity = n.value().mass * Point(0, 0, -grav);

    // Return total sum of spring force and gravity force.
    return f_spring + f_gravity;
  }
};

/** Generic function object that represents two combined forces. */
template <class F1, class F2>
struct TwoCombinedF {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return force1(n, t) + force2(n, t);
  }

  F1 force1;
  F2 force2;
};

/** Generic function object that represents three combined forces. */
template <class F1, class F2, class F3>
struct ThreeCombinedF {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return force1(n, t) + force2(n, t) + force3(n, t);
  }

  F1 force1;
  F2 force2;
  F3 force3;
};

/** Method to combine two forces. */
template <class F1, class F2>
TwoCombinedF <F1, F2> make_combined_force(F1 f1, F2 f2) {
  TwoCombinedF <F1, F2> combined_force;
  combined_force.force1 = f1;
  combined_force.force2 = f2;  
  return combined_force;
}

/** Method to combine three forces. */
template <class F1, class F2, class F3>
ThreeCombinedF <F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  ThreeCombinedF <F1, F2, F3> combined_force;
  combined_force.force1 = f1;
  combined_force.force2 = f2;  
  combined_force.force3 = f3;  
  return combined_force;
}

/** Gravity force function object. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Spring force function object. */
struct MassSpringForce {
  /** Return the spring force applying to @a n. */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    int num_nodes = n.graph()->size();
   
    // Set mass to 1/N where N is # of nodes in graph (constant density).
    n.value().mass = 1.0 / double(num_nodes); 
    
    // Initialize all edges to have rest length equal to initial edge lenth. 
    const double L = 1.0 / (sqrt(num_nodes)-1);

    // Spring force.
    Point f_spring = Point(0, 0, 0); // Initialize to zero vector.

    // Update spring force on node according to values of adjacent nodes.
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      assert(e.node1() == n);
   
      double edge_len = e.length();  
      double K = e.value();  // Edge-specific spring constant K_ij.    
      Point adj_pos = e.node2().position(); // Position of adjacent node.
      Point displacement = (n.position() - adj_pos) / edge_len;
      
      f_spring += (displacement) * (-K * (edge_len - L));
    }
    return f_spring;
  }
};

/** Damping force function object. */
struct DampingForce {
  /** Return the damping force applying to @a n. */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    int num_nodes = n.graph()->size();
    double damp_coef = 1.0 / num_nodes;
    Point f_damp = (-1.0 * damp_coef) * n.value().vel;
    return f_damp;
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
    graph.add_edge(nodes[t[0]], nodes[t[1]], spring_constant);
    graph.add_edge(nodes[t[0]], nodes[t[2]], spring_constant);
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]], spring_constant);
    graph.add_edge(nodes[t[2]], nodes[t[3]], spring_constant);
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
        // Use combined force of mass spring, damping, and gravity.
        ThreeCombinedF<GravityForce, MassSpringForce, DampingForce> combinedF =
            make_combined_force(GravityForce(), MassSpringForce(), DampingForce());
         symp_euler_step(graph, t, dt, combinedF);

        // TODO: Uncomment once node and edge removal fixed.
        //viewer.clear();
        //node_map.clear();

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
