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

// HW2 
struct EdgeData {
  double K = 100;
  double L;

  EdgeData() {}
  EdgeData(double springconst, double len) {
      K = springconst;
      L = len; 
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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue; 
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue; 
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  constraints(g, t); 
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
  if (n. position () == Point (0 ,0 ,0) || n. position () == Point (1 ,0 ,0))    // two endpoints
    return Point (0 ,0 ,0);

  Point f_grav = Point(0, 0, (-1.0) * grav * n.value().mass);    //  gravitational force

  Point f_spring = Point(0, 0, 0);
  for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {    // loop thru adj list 
    Node n1 = (*it).node2();
    auto value = (*it).value();
    //if (n.graph.has_edge(n1, n)) {
    //  value = (*it).reverse().value();
    //}

    double distance = norm_2(n.position()-n1.position());
    if (distance) {
      // f_spring += (-1) * 100 * (distance - 1) * (n.position() - n1.position())/ distance;    // update spring force
      f_spring += (-1) * value.K * (distance - value.L) * (n.position() - n1.position())/ distance; 
    }
  }
  (void) t;

  return f_grav +  f_spring;    // gravitational force + spring force
  }
};

struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0, 0, (-1.0) * grav * n.value().mass); 
  }
};

struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
  Point f_spring = Point(0, 0, 0);
  for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {    // loop thru adj list 
    Node n1 = (*it).node2();
    auto value = (*it).value();

    double distance = norm_2(n.position()-n1.position());
    if (distance) {
      // f_spring += (-1) * 100 * (distance - 1) * (n.position() - n1.position())/ distance;    // update spring force
      f_spring += (-1) * value.K * (distance - value.L) * (n.position() - n1.position())/ distance; 
    }
  }
  (void) t;

  return f_spring;    
  }
};

struct DampingForce{
  // members and constructors 
  double c; 
  DampingForce() {}
  DampingForce(double c1):c(c1) {}

  // damping force 
  template<typename NODE> 
  Point operator()(NODE n, double t) {
    (void) t;
    return -1.0 * c * n.value().vel; 
  }
};

// combine two forces
template <typename F1, typename F2>
struct CombineForce {
  // combine two forces
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; 
    return f1(n, t) + f2(n, t);
  }

  // constructor
  F1 f1;
  F2 f2;
  CombineForce(F1 f11, F2 f22):f1(f11), f2(f22) {}
};

// combine 2 forces 
template <typename F1, typename F2>
CombineForce<F1, F2> make_combined_force(F1 f11, F2 f22) {
  return CombineForce<F1, F2>(f11, f22);  
}

// combine 3 forces 
template <typename F1, typename F2, typename F3>
CombineForce<CombineForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  CombineForce<F1, F2> f_combined = make_combined_force(f1, f2);  
  return make_combined_force(f_combined, f3);
}

// constraints
// constant node constraint 
struct ConstantNodeConstraint {
  ConstantNodeConstraint() {}

  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    for (auto i = graph.node_begin(); i != graph.node_end(); ++i) {
      if ((*i).position() == Point(0, 0, 0) || (*i).position() == Point(1, 0, 0)) 
        (*i).value().vel = Point(0, 0, 0);     // constant node 
    }
    (void) t;
  }
};


// plane constraint 
struct PlaneConstraint {
  // constructor 
  double z; 
  PlaneConstraint(){}
  PlaneConstraint(double value): z(value){}

  // find nodes that violate the constraint, and reset them 
  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    for (auto i = graph.node_begin(); i != graph.node_end(); ++i) {     // loop thru all nodes
      if ((*i).position().elem[2] < z) {                            // constraint violated 
        (*i).position().elem[2] = z;  (*i).value().vel.elem[2] = 0;    // z-component 
      }
    }
    (void) t;
  }
};


//sphere constraint
struct SphereConstraint {
  Point c;
  double r;
  SphereConstraint(){}
  SphereConstraint(Point center, double radius) : c(center), r(radius) {}


  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    for (auto i = graph.node_begin(); i != graph.node_end(); ++i) {
      if (norm_2((*i).position() - c) < r) {
        auto Ri = ((*i).position() - c) / norm_2((*i).position() - c);
        (*i).position() = c + Ri * r;
        (*i).value().vel -= dot((*i).value().vel, Ri) * Ri;
      }
    }
    (void) t;
  }
};


// combine two constraints 
template <typename C1, typename C2>
struct CombineConstraint {
  CombineConstraint(C1 constraint1, C2 constraint2) : c1(constraint1), c2(constraint2) {}
  C1 c1;
  C2 c2; 

  template <typename GRAPH>
  void operator()(GRAPH& graph, double t) {
    c1(graph, t);
    c2(graph, t);
  }
};

// combine 2 constraints 
template <typename C1, typename C2>
CombineConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombineConstraint<C1, C2>(c1, c2);
}

// combine 3 constraints 
template <typename C1, typename C2, typename C3>
CombineConstraint<CombineConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return make_combined_constraint(make_combined_constraint(c1, c2), c3);
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

  int nn = graph.num_nodes(); 
  double mm = 1.0/nn;
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().vel = Point(0, 0, 0); 
    (*it).value().mass = mm;
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    (*it).value().L = (*it).length();
    (*it).value().K = 100.0; 
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
        //std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, 
          make_combined_force(GravityForce(), MassSpringForce()),
          make_combined_constraint(ConstantNodeConstraint(), PlaneConstraint(-0.75), SphereConstraint(Point(0.5, 0.5, -0.5), 0.15))); 

        // Updargte viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
