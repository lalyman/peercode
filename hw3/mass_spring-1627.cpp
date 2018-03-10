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

/** Custom structure of data to store with Edges 
 * This should not be added to Graph class as this 
 * application specific. 
 */
struct EdgeData {
  double K; // spring constant 
  double L; // initial length
  EdgeData(double spr_const, double len) : K(spr_const), L(len) {

  }
  EdgeData() {

  }
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>; // added EdgeData
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
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
      n.position() += n.value().vel * dt; 
    }
  }
  constraint(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
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
    // HW2 #1: YOUR CODE HERE
    Point f_spr = Point(0, 0, 0); // force spring
    Point f_grav = n.value().mass * Point(0, 0, -grav);// force gravitational
    
    // loop through all the edges incident on n and find force
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
        return Point(0, 0, 0); // zero-valued force
      }
      auto e = *it; 
      // using equation from Page 3 of handout
      Point vec = Point(0, 0, 0);
      if (e.node1() == n)
        vec = n.position() - e.node2().position(); // direction vector 
      else 
        vec = n.position() - e.node1().position();

      double norm_vec = norm(vec); // norm is defined in Point.hpp

      f_spr = f_spr - (e.value().K * (vec / norm_vec) * (norm_vec - e.value().L));
    }
    (void) t;// unused at the moment -- silcence warning
    return f_spr + f_grav;
  }
};


struct GravityForce {

  // implement similar to operator() in Problem1Force
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point f_grav = n.value().mass * Point(0, 0, -grav);// force gravitational
    (void) t;
    return f_grav;
  }

};

struct MassSpringForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    Point f_spr = Point(0, 0, 0); // force spring
    
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {

      auto e = *it; 
      // using equation from Page 3 of handout
      Point vec = Point(0, 0, 0);
      if (e.node1() == n)
        vec = n.position() - e.node2().position(); // direction vector 
      else 
        vec = n.position() - e.node1().position();
      double norm_vec = norm(vec); // norm is defined in Point.hpp

      f_spr = f_spr - (e.value().K * (vec / norm_vec) * (norm_vec - e.value().L));
    }
    (void) t; // unused at the moment
    return f_spr;
  }


};

struct DampingForce {
  double damp_coeff; 

  DampingForce (double c) : damp_coeff(c) {

  } 

  DampingForce () {

  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point f_damp = -n.value().vel * damp_coeff;// ffrom Page 7 of handout
    (void) t;
    return f_damp;
  }

};


template<typename Force1, typename Force2>
struct CombinedForce {
  Force1 f1; 
  Force2 f2;

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t);
  }

};

template<typename Force1, typename Force2>
CombinedForce<Force1, Force2> make_combined_force(Force1 f1, Force2 f2) {
  return {f1, f2};
}

template<typename Force1, typename Force2, typename Force3>
CombinedForce<CombinedForce<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3) {
  return make_combined_force(make_combined_force(f1, f2), f3); 
}


// Representation for Plane Constriant
struct PlaneConstraint {

  void operator()(GraphType& g, double t) {

    double c1 = -0.75;

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (dot(n.position(), Point(0, 0, 1)) < c1) {
        n.position().z = c1; // other coordinates are 0
        n.value().vel.z = 0.0;
      }
    }

    (void) t;
  }
};

// Representation for Spherical Constriant1
struct SphereConstraint1 {
  void operator()(GraphType& g, double t) {

    double r = 0.15; // radius
    Point c = Point(0.5, 0.5, -0.5); // center 

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point vec = n.position() - c;
      double d = norm(vec); // d can be zero

      if (d < 0.000000001) { // handle d is 0 
        // center and point coincide. So set the point to some point on the sphere
        n.position().z = n.position().z + r;
        n.value().vel.z = 0.0;
      } else if (d < r) {
          n.position() = c + (vec * (r / d));
          n.value().vel = n.value().vel - dot(n.value().vel, vec / d) * (vec / d);
      }
    }

    (void) t;
  }

};

// Representation for Spherical Constriant2
struct SphereConstraint2 {
  void operator()(GraphType& g, double t) {

    double r = 0.15; // radius
    Point c = Point(0.5, 0.5, -0.5); // center 


    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point vec = n.position() - c;
      double d = norm(vec); 

      if (d < r) {
        g.remove_node(n);
      }

    }

    (void) t;
  }

};

// Combining constraints as we did with forces 
template<typename C1, typename C2>
struct CombinedConstraints {
  C1 c1;
  C2 c2; 

  void operator()(GraphType& g, double t) {
    c1(g, t);
    c2(g, t);
  }
};


template<typename C1, typename C2> 
CombinedConstraints<C1, C2> make_combined_constraints(C1 c1, C2 c2) {
  return {c1, c2};
}

template<typename C1, typename C2, typename C3>
CombinedConstraints<CombinedConstraints<C1, C2>, C3> make_combined_constraints(C1 c1, C2 c2, C3 c3) {
  return make_combined_constraints(make_combined_constraints(c1, c2), c3);

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
  double spr_const = 100;
  double temp_mass = 1.0 / graph.num_nodes();

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0, 0, 0);
    n.value().mass = temp_mass;
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    e.value().K = spr_const; 
    e.value().L = e.length();
  }

  // Calculate all the forces 
  GravityForce f_g;
  MassSpringForce f_s;
  double damp_coeff = 1.0 / graph.num_nodes();
  DampingForce f_d = DampingForce(damp_coeff);

  auto total_f = make_combined_force(f_g, f_s, f_d);

  // Constraints section 
  PlaneConstraint p1;
  SphereConstraint1 s1;
  SphereConstraint2 s2;

  auto all_constraints = make_combined_constraints(p1, s1, s2);

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
        symp_euler_step(graph, t, dt, total_f, all_constraints);

        // Clear the viewers's nodes and edges
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
