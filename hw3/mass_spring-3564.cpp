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
  double length;     //< Edge length
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

    // When ConstantNodeConstraint not enforced, skip the update step for (0, 0, 0) and (1, 0, 0).
    //if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    //  continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // When ConstantNodeConstraint not enforced, skip the update step for (0, 0, 0) and (1, 0, 0).
    //if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    //  continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  // Set constraints
  constraint(g, t);

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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);

    // 1. Spring force
    Point f_spring = Point(0);
    double K = 100; // spring constant
    // Traverse adjacent nodes 
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      double L = (*ei).value().length;
      Node n2 = (*ei).node2();
      double dis = norm(n.position() - n2.position());
      f_spring += (n.position() - n2.position()) / dis * (dis - L);
    }
    f_spring *= -K;

    // 2. Force due to gravity
    Point f_grav = n.value().mass * Point(0, 0, -grav);

    return f_spring + f_grav;
  }
};


/** Function object for gravity. */
struct GravityForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return n.value().mass * Point(0, 0, -grav);
  }
};


/** Function object for mass spring force. */
struct MassSpringForce {

  double K_; // spring constant

  /** MassSpringForce Constructor.
   * @param[in] K Spring Constant.
   */
  MassSpringForce(double K = 100) : K_(K) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {

    Point f_spring = Point(0);

    // Traverse adjacent nodes 
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      double L = (*ei).value().length;
      Node n2 = (*ei).node2();
      double dis = norm(n.position() - n2.position());
      f_spring += (n.position() - n2.position()) / dis * (dis - L);
    }
    f_spring *= -K_;

    return f_spring;
  }
};

/** Function object for damping force */
struct DampingForce {
 
  double c_; // spring constant

  /** DampingForce Constructor.
   * @param[in] c Damping coefficient.
   */
  DampingForce(double c) : c_(c) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return -n.value().vel * c_;
  }
};


/** Function object that returns a combination of forces. 
 * @param[in] Two forces f1 and f2.
 */
template <typename F1, typename F2>
struct CombinedForce {
  F1 f1_;
  F2 f2_;

  /** CombinedForce constructor 
   * @param[in] f1 First force.
   * @param[in] f2 Second force.
   */
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  /** Calculate combined forces. 
   * @param[in] n Node.
   * @param[in] t Time.
   * @return Point object that represents the combination of forces of @a f1 and @a f2.
   */
  Point operator() (Node n, double t) {
    return f1_(n, t) + f2_(n, t);
  }
};


/** Function that returns a combination of two forces.
 * @param[in] f1, f2, f3 Forces that take a node and time as input.
 * @return A CombinedForce object that sum up the two forces.
 */
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1, F2> (f1, f2);
}


/** Force Function that returns a combination of three forces. 
 * @param[in] f1, f2, f3 Forces that take a node and time as input.
 * @return A CombinedForce object that sum up the three forces.
 */
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombinedForce<CombinedForce<F1, F2>, F3> (CombinedForce<F1, F2> (f1, f2), f3);
}


/** Constraint that fixes specific nodes.
 * Complexity: O(1)
 */
struct ConstantNodeConstraint {

  std::vector<GraphType::size_type> indices_; // indices of fixed nodes 

  /** ConstantNodeConstraint Constructor
   * @param[in] v A vector of indices of fixed nodes 
   */
  ConstantNodeConstraint(const std::vector<GraphType::size_type>& v) : indices_(v) {}  

  /** ConstantNodeConstraint Setter
   * @param[in] g Graph.
   * @param[in] t Time.
   * @post The velocity of fixed points is set to 0. 
   */
  void operator() (GraphType& g, double t) {
    for (auto it = indices_.begin(); it != indices_.end(); ++it) {
      auto n = g.node(*it);
      // Fix the node by setting velocity to 0.
      n.value().vel = Point(0, 0, 0); 
    }
  }
};


/** Constraint that defines an impenetrable plane obstacle. 
 * Complexity: O(num_nodes())
 */
struct PlaneConstraint {

  double z_; 

  /** PlaneConstraint Constructor
   * @param[in] z The z-coordinate of the plane. 
   */
  PlaneConstraint(double z) : z_(z) {}  

  /** PlaneConstraint Setter
   * @param[in] g Graph.
   * @param[in] t Time.
   * @post For the node that violates the constraint, 
   *       the z-component of its velocity is set to 0.
   *       its position is set to the nearest point on the plane.  
   */
  void operator() (GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // A node violates this constraint
      if(dot(n.position(), Point(0, 0, 1)) < z_) {
        // Set the position to the nearest point on the plain. 
        n.position().z = z_;
        // Set the z-component of the Node velocity to zero. 
        n.value().vel.z = 0;
      }
    }
  }
};


/** Constraint that defines an impenetrable sphere obstacle. 
 * Complexity: O(num_nodes())
 */
struct SphereConstraint {
  
  Point center_;
  double radius_;

  /** SphereConstraint Constructor
   * @param[in] center The center of the sphere obstacle. 
   * @param[in] radius The radius of the sphere obstacle.
   */
  SphereConstraint(Point center, double radius) : center_(center), radius_(radius) {}  

  /** SphereConstraint Setter
   * @param[in] g Graph.
   * @param[in] t Time.
   * @post For the node that violates the constraint, 
   *       its velocity component normal to the sphere's surface is set to 0.
   *       its position is set to the nearest point on the surface of the sphere.  
   */
  void operator() (GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // A node violates this constraint
      if(norm(n.position() - center_) < radius_) {
        Point R = (n.position() - center_) / norm(n.position() - center_); 
        // Set the position to the nearest point on the surface of the sphere. 
        n.position() = center_ + R * radius_;
        // Set the velocity component normal to the sphere's surface to zero.
        n.value().vel -= dot(n.value().vel, R) * R;
      }
    }
  }
};


/** Constraint that defines a sphere obstacle which can destroy the graph. 
 * Complexity: O(num_nodes())
 */
struct SphereConstraint2 {
  
  Point center_;
  double radius_;

  /** SphereConstraint2 Constructor
   * @param[in] center The center of the sphere obstacle. 
   * @param[in] radius The radius of the sphere obstacle.
   */
  SphereConstraint2(Point center, double radius) : center_(center), radius_(radius) {}  

  /** SphereConstraint2 Setter
   * @param[in] g Graph.
   * @param[in] t Time.
   * @post The node that violates the constraint is removed from the graph. 
   */
  void operator() (GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // A node violates this constraint
      if(norm(n.position() - center_) < radius_) {
        g.remove_node(n);
      }
    }
  }
};


/** Function object that returns a combination of constraints. 
 * @param[in] Two constraints c1 and c2.
 */
template <typename C1, typename C2>
struct CombinedConstraint {
  C1 c1_;
  C2 c2_;

  /** CombinedConstraint constructor 
   * @param[in] c1 First constraint.
   * @param[in] c2 Second constraint.
   */
  CombinedConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}

  /** Calculate combined constraints. 
   * @param[in] g Graph.
   * @param[in] t Time.
   * @post Both constraints have been applied. 
   */
  void operator() (GraphType& g, double t) {
    c1_(g, t);
    c2_(g, t);
  }
};


/** Function that returns a combination of two constraints.
 * @param[in] c1, c2 Constraints that take a graph and time as input.
 * @return A CombinedConstraint object that applies both constraints.
 */
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombinedConstraint<C1, C2> (c1, c2);
}


/** Function that returns a combination of three constraints. 
 * @param[in] c1, c2, c3 Constraints that take a graph and time as input.
 * @return A CombinedConstraint object that applies all three constraints.
 */
template <typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return CombinedConstraint<CombinedConstraint<C1, C2>, C3> (CombinedConstraint<C1, C2> (c1, c2), c3);
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

  std::vector<GraphType::size_type> indices(2);

  // 1. Set mass and get the indices of nodes (0, 0, 0) and (1, 0, 0).
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().mass /= graph.num_nodes();
    if (n.position() == Point(0, 0, 0)) {
      indices[0] = n.index();
    } else if (n.position() == Point(1, 0, 0)){
      indices[1] = n.index();
    }
  }
  
  // 2. Set rest length for all edges to their initial length.
  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    (*ei).value().length = (*ei).length();
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
      //double dt = 0.0005; // for grid3
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        
        // Check whether the index of (0, 0, 0) has changed.
        // Search its new index if it does change. 
        if (indices[0] >= graph.size()) {
          for (auto p1 = graph.node_begin(); p1 != graph.node_end(); ++p1) {
            if ((*p1).position() == Point(0, 0, 0)) {
              indices[0] = (*p1).index();
              break;
            }
          }
        }
    
        // Check whether the index of (1, 0, 0) has changed.
        // Search its new index if it does change. 
        if (indices[1] >= graph.size()) {
          for (auto p2 = graph.node_begin(); p2 != graph.node_end(); ++p2) {
            if ((*p2).position() == Point(1, 0, 0)) {   
              std::cout << (*p2).index() << std::endl;
              indices[1] = (*p2).index();
              break;
            }
          }
        }

        // Set damping coefficient = 1 / N, where N is the number of nodes in the Graph.
        double c = 1.0 / graph.num_nodes();

        //symp_euler_step(graph, t, dt, Problem1Force());
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)), PlaneConstraint(-0.75));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)), SphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)), SphereConstraint2(Point(0.5, 0.5, -0.5), 0.15));
        
        // Apply gravity, mass spring force and damping force.
        // Fix (0, 0, 0) and (1, 0, 0) and apply the sphere obstacle which destroys the graph. 
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)), make_combined_constraint(SphereConstraint2(Point(0.5, 0.5, -0.5), 0.15), ConstantNodeConstraint(indices)));
  
        // Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        /*
        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.set_label(t);
        */

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
