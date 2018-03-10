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
  double spring_constant;     //< Edge spring constant
  double rest_length;         //< Edge rest length
  EdgeData() : spring_constant(100), rest_length(.1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


//
//  CONSTRAINTS
//

/** Constraint function object for constant constraint by fixed point. */
struct ConstantConstraint {
  Point fixed_point_;
  Node n_;

  ConstantConstraint(Point fixed_point, Node n) : fixed_point_(fixed_point), n_(n) {}

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    if (n == n_)
      n.position() = fixed_point_;
  }
};

/** Constraint function object for plane constraint by unit normal to plance and intercept term. */
struct PlaneConstraint {
  Point unit_normal_;
  double intercept_;

  PlaneConstraint(Point unit_normal, double intercept) : unit_normal_(unit_normal), intercept_(intercept) {}

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    double projection = inner_prod(n.position(), unit_normal_);
    if (projection < intercept_)
      n.position() = n.position() + (unit_normal_ * (intercept_ - projection));
  }
};

/** Constraint function object for enforce sphere constraint by center and radius. */
struct SphereConstraint {
  Point center_;
  double radius_;

  SphereConstraint(Point center, double radius) : center_(center), radius_(radius) {}

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    double dist = norm(n.position() - center_);
    if (dist < radius_) {
      Point normal = (n.position() - center_) / dist;
      n.position() += (radius_ - dist) * normal;
      n.value().vel -= inner_prod(n.value().vel, normal) * normal;
    }
  }
};

/** Constraint function object for removal through sphere by center and radius. */
struct SphereConstraintRemove {
  GraphType& g_;
  Point center_;
  double radius_;

  SphereConstraintRemove(GraphType& g, Point center, double radius) : g_(g), center_(center), radius_(radius) {}

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    if (norm(n.position() - center_) < radius_)
      g_.remove_node(n);
  }
};

/** Generic constraint function object to containt all constraint function objects. */
template <typename C1, typename C2>
struct TotalConstraint {
  C1 c1_1_;
  C1 c1_2_;
  C2 c2_;

  TotalConstraint(C1 c1_1, C1 c1_2, C2 c2) : c1_1_(c1_1), c1_2_(c1_2), c2_(c2) {}

  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      c1_1_(n,t);
      c1_2_(n,t);
      c2_(n,t);
    }
  }
};

/** Generic method to combine constraints. */
template <typename C1, typename C2>
TotalConstraint<C1,C2> make_combined_constraints(C1 c1_1, C1 c1_2, C2 c2) {
  TotalConstraint<C1,C2> constraint{c1_1, c1_2, c2};
  return constraint;
}

//
//  SYMPLECTIC EULER STEP
//

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

  // Enforce constraint
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


//
//  FORCES
//

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

    // Corner constraints
    Point point_i = n.position();
    if (point_i == Point(0,0,0) || point_i == Point(1,0,0))
      return Point(0,0,0);

    // Calculate spring force due to each adjacent node
    Point spring_force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double Kij = (*it).value().spring_constant;
      double Lij = (*it).value().rest_length;

      Point point_j = (*it).node2().position();
      double dist_ij = norm(point_i - point_j);
      
      Point spring_force_ij = -Kij * ( (point_i - point_j) * (1. - Lij / dist_ij) );
      spring_force += spring_force_ij;
    }

    // Calculate gravitational force
    Point grav_force = n.value().mass * Point(0,0,-grav);

    return spring_force + grav_force;
  }
};


/** Force function object for gravity. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    // Calculate gravitational force
    Point grav_force = n.value().mass * Point(0,0,-grav);

    return grav_force;
  }
};

/** Force function object for spring force. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    // Calculate spring force due to each adjacent node
    Point point_i = n.position();
    Point spring_force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double Kij = (*it).value().spring_constant;
      double Lij = (*it).value().rest_length;

      Point point_j = (*it).node2().position();
      double dist_ij = norm(point_i - point_j);
      
      Point spring_force_ij = -Kij * ( (point_i - point_j) * (1. - Lij / dist_ij) );
      spring_force += spring_force_ij;
    }

    return spring_force;
  }
};

/** Force function object for damping force. */
struct DampingForce {
  double damping_constant_;

  DampingForce(double damping_constant) : damping_constant_(damping_constant) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    // Calculate damping force
    Point damping_force = -damping_constant_ * n.value().vel;

    return damping_force;
  }
};

/** Force function object for zero force to be used as default. */
struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t;
    return Point(0,0,0);
  }
};

/** Generic force function object to contain all force function objects. */
template <typename F1, typename F2, typename F3 = ZeroForce>
struct TotalForce {
    F1 f1_;
    F2 f2_;
    F3 f3_;

    TotalForce(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      return f1_(n,t) + f2_(n,t) + f3_(n,t);
    }
};

/** Generic method to combine force function objects. */
template <typename F1, typename F2, typename F3 = ZeroForce>
TotalForce<F1,F2,F3> make_combined_forces(F1 f1, F2 f2, F3 f3 = F3()) {
  TotalForce<F1,F2,F3> force{f1, f2, f3};
  return force;
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

  // Set initial conditions for nodes.
  NodeData init_node_data{};
  init_node_data.mass = 1. / ((double) graph.size());
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value() = init_node_data;
  }

  // Set initial conditions for edges.
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    EdgeData init_edge_data{};
    Point pi = (*it).node1().position();
    Point pj = (*it).node2().position();
    init_edge_data.rest_length = norm(pi - pj);
    (*it).value() = init_edge_data;
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  // Add nodes and edges.
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      // Find fixed nodes to constrain.
      Node n1;
      Node n2;
      for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        auto n = *it;
        if (n.position() == Point(0,0,0))
          n1 = n;
        else if (n.position() == Point(1,0,0))
          n2 = n; 
      }

      // Initialize total force function object.
      TotalForce<GravityForce,MassSpringForce,DampingForce> force = make_combined_forces(GravityForce(), MassSpringForce(), DampingForce{1. / ((double) graph.size())});
        
      // Initialize constant constraint function objects.
      ConstantConstraint c1_1{Point(0,0,0), n1};
      ConstantConstraint c1_2{Point(1,0,0), n2};
      
      // Initialize plane constraint function object and accumulate.
      //PlaneConstraint c2{Point(0,0,1), -0.75};
      //TotalConstraint<ConstantConstraint,PlaneConstraint> constraint = make_combined_constraints(c1_1, c1_2, c2);
      
      // Initilize sphere constraint function object and accumulate.
      //SphereConstraint c2{Point(0.5,0.5,-0.5), 0.15};
      //TotalConstraint<ConstantConstraint,SphereConstraint> constraint = make_combined_constraints(c1_1, c1_2, c2);

      // Initialize sphere removal function object and accumulate.
      SphereConstraintRemove c2{graph, Point(0.5,0.5,-0.5), 0.15};
      TotalConstraint<ConstantConstraint,SphereConstraintRemove> constraint = make_combined_constraints(c1_1, c1_2, c2);


      // Run simulation.
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        
        // Run symplectic euler step.
        symp_euler_step(graph, t, dt, force, constraint);

        // Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update the viewer with nodes' new positions and new edges. This causes flickering.
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        
        // Update viewer with nodes' new positions. Use this for nonremoval constraints to remove flickering.
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        
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
