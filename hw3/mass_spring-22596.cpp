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
  double K;
  double L;
  EdgeData() : K(100.0), L(0) {}
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

// HW2: YOUR CODE HERE
// Constraints HW2 Problem 4

struct ZPlaneConstraint {
  double z_;
  ZPlaneConstraint(double z) : z_(z) {}
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = *iter;
      if (dot(n.position(), Point(0, 0, 1)) < z_) {
        n.position()[2] = z_;
        n.value().vel[2] = 0;
      }
    }
  }
};

struct SphereConstraint {
  Point center_;
  double radius_;
  SphereConstraint(Point center, double radius) :
    center_(center), radius_(radius) {}
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = *iter;
      if (norm(n.position() - center_) < radius_) {
        Point vec = n.position() - center_;
        Point unit_vec = vec/norm(vec);
        n.position() =  center_ + unit_vec*radius_;
        n.value().vel = n.value().vel - (dot(n.value().vel, unit_vec))*unit_vec;
      }
    }
  }
};

// Constraint for HW2, Problem 5
struct RemoveSphereConstraint {
  Point center_;
  double radius_;
  RemoveSphereConstraint(Point center, double radius) :
    center_(center), radius_(radius) {}
  template <typename G>
  void operator()(G& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = (*iter);
      if (norm(n.position() - center_) < radius_) {
        unsigned removed_node = g.remove_node(n);
        (void) removed_node;
      }
    }
  }
};

// Combining constraints
template <typename Constraint1, typename Constraint2>
struct CombineConstraints {
  Constraint1 cons1_;
  Constraint2 cons2_;
  CombineConstraints(Constraint1 cons1, Constraint2 cons2):
  cons1_(cons1), cons2_(cons2) {}
  void operator()(GraphType& g, double t) {
    cons1_(g, t);
    cons2_(g, t);
  }
};

template <typename Constraint1, typename Constraint2>
CombineConstraints<Constraint1, Constraint2> make_combined_constraint(
  Constraint1 cons1, Constraint2 cons2) {
    return CombineConstraints<Constraint1, Constraint2> (cons1, cons2);
  }


// TESTS FOR HW2 PROBLEM 4 and 5
ZPlaneConstraint plane_z(-0.75);
SphereConstraint sphere(Point(0.5, 0.5, -0.5), 0.15);
RemoveSphereConstraint rmsphere(Point(0.5, 0.5, -0.5), 0.15);



template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);

    // HW2: YOUR CODE HERE
    // Setting the constant nodes constraints
    // This can be done outside too.
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      n.value().vel = Point(0, 0, 0);
    }
  }

//  auto all_constr = make_combined_constraint(
//                    make_combined_constraint(plane_z, sphere), rmsphere);
//  auto all_constr = make_combined_constraint(plane_z, sphere);
  auto all_constr = make_combined_constraint(rmsphere, plane_z);
  all_constr(g, t);

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
    (void) n; (void) t; (void) grav;    // silence compiler warnings

    Point f_spring  = Point(0, 0, 0);
    Point f_gravity = Point(0, 0, -grav)*n.value().mass;

    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    else {
      for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
        Point neighbour = (*iter).node2().position();
        Point distance  = n.position() - neighbour;
        f_spring -= (*iter).value().K*
        (norm(distance) - (*iter).value().L)/norm(distance)*distance;
      }
      return f_spring + f_gravity;
    }
  }
};


// Graviational force
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t; (void) grav;
    Point f_gravity = Point(0, 0, -grav)*n.value().mass;
    return f_gravity;
  }
};

// Mass Spring force
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t; (void) grav;
    Point f_spring = Point(0, 0, 0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Point neighbour = (*iter).node2().position();
      Point distance  = n.position() - neighbour;
      f_spring -= (*iter).value().K*
        (norm(distance) - (*iter).value().L)/norm(distance)*distance;
    }
    return f_spring;
  }
};

// Damping force
struct DampingForce {
  double c_;
  DampingForce(double c): c_(c) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t; (void) grav;
    Point f_damping = -c_*n.value().vel;
    return f_damping;
  }
};

// Combining forces
template <typename Force1, typename Force2>
struct CombineForces {
  Force1 force1_;
  Force2 force2_;
  CombineForces(Force1 force1, Force2 force2):
  force1_(force1), force2_(force2) {}
  Point operator()(Node n, double t){
    (void) n; (void) t;
    return force1_(n, t) + force2_(n, t);
  }
};

template <typename Force1, typename Force2>
CombineForces<Force1, Force2> make_combined_force(Force1 force1, Force2 force2) {
  return CombineForces<Force1, Force2>(force1, force2);
}

template <typename Force1, typename Force2, typename Force3>
CombineForces<CombineForces<Force1, Force2>, Force3> make_combined_force(Force1 force1,
  Force2 force2, Force3 force3) {
    return CombineForces<CombineForces<Force1, Force2>, Force3> (
      CombineForces<Force1, Force2>(force1, force2), force3);
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
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
    auto node = *iter;
    node.value().vel = Point(0, 0, 0);
    node.value().mass = 1.0/graph.num_nodes();
  }

  for (auto iter = graph.edge_begin(); iter != graph.edge_end(); ++iter) {
    auto edge = *iter;
    edge.value().K = 100.0;
    edge.value().L = edge.length();
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        DampingForce f_damping(1.0/graph.num_nodes());
        auto total_force = make_combined_force(GravityForce(),
                                               MassSpringForce(),
                                               f_damping);

        symp_euler_step(graph, t, dt, total_force);

        //std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());

        // Update viewer with nodes' new positions
//        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
