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
  double length;   //< Edge rest length
  double K; //< Edge spring constant
  EdgeData() : length(0), K(0) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;

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
    // HW2 #3: skip steps 
    // if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
    //   n.position() += n.value().vel * dt;
    // }
  }

  // HW2 #4: Find violated constraints after and reset nodes’ 
  // positions and velocities.
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    // HW2 #3: skip steps 
    // if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
    //   n.value().vel += force(n, t) * (dt + n.value().mass);
    // }
  } 

  return t + dt;
}

// FORCES 

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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point xi = n.position();
    Point force = Point(0, 0, - grav * n.value().mass);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Point xj = (*it).node2().position();
      force += - K_ * (xi - xj)  * ((*it).length() - L_) / (*it).length();
    }
    (void) t;
    return force;
  }
  /** Constant spring and rest-length. */
  double K_;
  double L_; 
  
  /** Constructor. */
  Problem1Force(double K, double L) : K_(K), L_(L) {} 
};


/** Force function object for HW2 #2. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #2, each edge may have its own spring constant and rest-length.
   * For HW2 #2, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point xi = n.position();
    Point force = Point(0, 0, - grav * n.value().mass);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Point xj = (*it).node2().position();
      force += - (*it).value().K * (xi - xj) * ((*it).length() - (*it).value().length) / (*it).length();
    }
    (void) t;
    return force;
  }
};


/** Force function object for HW2 #3. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point gravity_force = Point(0, 0, -grav * n.value().mass);
    (void) t;
    return gravity_force;
  }
};


/** Force function object for HW2 #3. */
struct MassSpringForce {
  /** Return the mass spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point xi = n.position();
    Point spring_force = Point(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Point xj = (*it).node2().position();
      spring_force += - (*it).value().K * (xi - xj)  * ((*it).length() - (*it).value().length) / (*it).length();
    }
    (void) t;
    return spring_force;
  }
};


/** Force function object for HW2 #3. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point damping_force = - c_ * n.value().vel;
    (void) t;
    return damping_force;
  }

  /** Construrctor. */
  DampingForce(double c) : c_(c) {}
  
  /** Damping constant. */
  double c_; 
};

// COMBINE FORCES 

/** Struct to combine two forces. */
template<typename F1, typename F2>
struct CombinedForce {
  /** Return the combined force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) { 
    return f1_(n, t) + f2_(n, t);
  }
  F1 f1_; 
  F2 f2_; 
  /** Construrctor. */
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
};


/** Function to combine two forces */
template<typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1, F2>(f1, f2);
}


/** Function to combine three forces */
template<typename F1, typename F2, typename F3>
CombinedForce<F1, CombinedForce<F2, F3>> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombinedForce<F1, CombinedForce<F2, F3>>(f1, CombinedForce<F2, F3>(f2, f3));
}


// HW2 #4: CONSTRAINTS 

/** Constraint function object for HW2 #4. */
struct ConstantConstraint {
  /** Set constraints for one point @a p that
   * that never moves in @a graph at any time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    (void) t;
    graph.node(node_idx_).position() = p_;
  }
  Point p_;
  size_type node_idx_;
  
  /** Construrctor. */
  template <typename GRAPH>
  ConstantConstraint(GRAPH& graph, Point p) : p_(p) {
    for (auto it=graph.node_begin(); it!=graph.node_end(); ++it) {
      Node n = *it;
      if (n.position() == p_) {
        node_idx_ = n.index();
	break;
      }
    } 
  }
};


/** Constraint function object for HW2 #4. */
struct PlaneConstraint {
  /** Find violated constraints after the position update and reset 
   * nodes’ positions and velocities for the plane constraint on the
   * plane z = @a z in @a graph at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    (void) t;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      Node n = *it;
      if (n.position().z < z_) {
	// Set the position to the nearest point on the plane.
        n.position().z = z_;
	// Set the z-component of the Node velocity to zero.
        n.value().vel.z = 0;
      }
    }
  }
  double z_;

  /** Constructor. */
  PlaneConstraint(double z) : z_(z) {} 
};



/** Constraint function object for HW2 #4. */
struct SphereConstraint {
  /** Find violated constraints after the position update and reset 
   * nodes’ positions and velocities for the sphere constraint with 
   * sphere center @a c and radius @a r in @a graph at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    (void) t;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      Node n = *it;
      if (norm(n.position() - c_) < r_) {
	Point R = (n.position() - c_) / norm(n.position() - c_);
	// Set the position to the nearest point on the surface of the sphere.
        n.position() = c_ + r_ * R;
	// Set the component of the velocity that is normal to the sphere’s surface to zero:
        n.value().vel -= dot(R, n.value().vel) * R; 
      }
    }
  }

  /** Sphere center and radius. */
  Point c_;
  double r_;

  /** Constructor. */
  SphereConstraint(Point c, double r) : c_(c), r_(r) {}

};


/** Constraint function object for HW2 #5. */
struct SphereConstraintRemove {
  /** Find violated constraints after the position update and reset 
   * nodes’ positions and velocities for the sphere constraint with 
   * sphere center @a c and radius @a r in @a graph at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    (void) t;
    for (auto it = graph.node_begin(); it!=graph.node_end(); ++it) {
      Node n = *it;
      if (norm(n.position() - c_) < r_) {
        graph.remove_node(n);
      }
    }
  }

  /** Sphere center and radius. */
  Point c_;
  double r_;

  /** Constructor. */
  SphereConstraintRemove(Point c, double r) : c_(c), r_(r) {}

};

/** Struct to combine two constraints. 
 * Assumes that the two constraints don't have any conflict. */
template<typename C1, typename C2>
struct CombinedConstraint {
  /** Return the combined constraint applied to @a graph at time @a t. */
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) { 
    (void) t;
    c1_(graph, t);
    c2_(graph, t);
  }

  C1 c1_; 
  C2 c2_; 
  
  /** Construrctor. */
  CombinedConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
};


/** Function to combine two constraints */
template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombinedConstraint<C1, C2>(c1, c2);
}


/** Function to combine three constraints */
template<typename C1, typename C2, typename C3>
CombinedConstraint<C1, CombinedConstraint<C2, C3>> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return CombinedConstraint<C1, CombinedConstraint<C2, C3>>(c1, CombinedConstraint<C2, C3>(c2, c3));
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
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif    
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) { 
    (*it).value().mass = 1.0 / graph.num_nodes();
    (*it).value().vel = Point(0, 0, 0);
  }

  // HW2 #2 
  // Set initial rest-length spring constant 
  // Can be different from one edge to another 
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    (*it).value().length = (*it).length();
    (*it).value().K = 100.0;
  }
  
  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // FORCES
  
  // HW2 #1 : Spring constant K and rest-length L
  // double K = 100.0;
  // double L = (*graph.edge_begin()).length();
  // auto force = Problem1Force(K, L);
  
  // HW2 #2: Edge-dependant K and L 
  // auto force = Problem2Force();

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
	
        // symp_euler_step(graph, t, dt, Problem1Force());
        // HW2 #1-2-3
	// symp_euler_step(graph, t, dt, force);
        // HW2 #3: Combined Forces
        // auto force = make_combined_force(GravityForce(), MassSpringForce());
        double damping_constant = 1.0 / graph.num_nodes();
        auto force = make_combined_force(GravityForce(), MassSpringForce(), DampingForce(damping_constant));

        // HW2 #4: Add constraint
        Point p1 = Point(0, 0, 0);
        Point p2 = Point(1, 0, 0);
        auto constant_constraint_1 = ConstantConstraint(graph, p1);
        auto constant_constraint_2 = ConstantConstraint(graph, p2);
        auto constant_constraint = make_combined_constraint(constant_constraint_1, constant_constraint_2);
        auto plane_constraint = PlaneConstraint(-0.75); 
        
	// HW2 #5: Constraint
        auto sphere_constraint = SphereConstraintRemove(Point(0.5, 0.5, -0.5), 0.15);
        // auto sphere_constraint = SphereConstraint(Point(0.5, 0.5, -0.5), 0.15);
        auto constraint = make_combined_constraint(constant_constraint, plane_constraint, sphere_constraint);
        symp_euler_step(graph, t, dt, force, constraint);

	// HW2 #5: Update the viewer with nodes's new positions and new edges. 
	viewer.clear();
	node_map.clear();
	viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
	viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        // Update viewer with nodes' new positions
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
