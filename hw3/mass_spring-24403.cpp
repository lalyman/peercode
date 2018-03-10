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
#include <cmath>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static int num_nodes; //Global variable for damping force

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custrom structure of data to store with Edges */
struct EdgeData {
  double K_;        //Spring constant
  double L_;        //Spring rest length
  EdgeData() : K_(1), L_(1) {} //default
  EdgeData(double K, double L) : K_(K), L_(L) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = typename Point::size_type;


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
//template <typename G, typename F, typename C>
template <typename G, typename F, typename C, typename CC>
//double symp_euler_step(G& g, double t, double dt, F force, C constrain) {
double symp_euler_step(G& g, double t, double dt, F force, C constrain, CC constconstrain) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it; 
    
    //fix corners
    constconstrain(n, t);

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  //apply constraint(s) to graph
  //constrain(g, t);
  constrain(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    //fix corners
    constconstrain(n, t);

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Constraint functor that fixes Point (0,0,0) and Point(1,0,0)
 *  Incompatible with RemoveSphere() functor
 */
/*
struct FixConstCorners {
 private:
  size_type n1_idx_; //uid of node at Point (0,0,0)
  size_type n2_idx_; //uid of node at Point (1,0,0)
 public:
  //constructor
  FixConstCorners(size_type n1_idx, size_type n2_idx) : n1_idx_(n1_idx), n2_idx_(n2_idx) {}
  //call operator
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) { 
    (void) t; //silence compiler
    //reset velocity and position of points at n1_idx and n2_idx to be zero
    g.node(n1_idx_).position() = Point(0,0,0);
    g.node(n1_idx_).value().vel = Point(0,0,0);
    g.node(n2_idx_).position() = Point(1,0,0);
    g.node(n2_idx_).value().vel = Point(0,0,0);

    return;
  }
};
*/

/** Constraint functor that fixes Point (0,0,0) and Point(1,0,0)
 *  Compatible with RemoveSphere() functor, but slower than
 *  FixConstCorners()
 */
struct FixCorners{
  //template <typename GRAPH>
  //void operator()(GRAPH& g, double t) { 
  void operator()(Node n, double t) { 
    (void) t; //silence compiler
    if (n.position() == Point(0,0,0) or n.position() == Point(1,0,0)) {
      n.value().vel = Point(0,0,0);
    }
  }
};


/** Constraint functor for plane z = -0.75 */
struct ConstrainPlane {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) { 
    (void) t; //silence compiler
    //iterate through nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      //check if node violates constraint (z component is < -0.75)
      if (n.position().z < -0.75) {
        //set position to nearest point on plane - i.e. keep x and y, just update x
        n.position().z = -0.75;
        //set z component of Node velocity to zero
        n.value().vel.z = 0.0;
      }
    }
    return; 
  } 
};

/** Constraint functor for sphere with center (0.5,0.5,-0.5) and r = 0.15
 *  When node fails constraint it is placed on surface of sphere 
 */
struct ConstrainSphere {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t; //silence compiler
    //iterate through nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      //check if node violates sphere constraint
      Point del = n.position() - Point(0.5,0.5,-0.5);  
      //check whether the norm of this difference is less than 0.15
      if (norm(del) < 0.15) {
        //set position to nearest node on plane, i.e. c + del scaled up by r 
        double scaling_factor = 0.15/norm(del);  //should be > 1
        n.position() = Point(0.5,0.5,-0.5) + scaling_factor*del; 
        //set component of velocity normal to sphere surface to zero
        Point R = del/norm(del);
        n.value().vel -= dot(n.value().vel, R)*R;
      }
    }
    return;
  }
};

/** Constraint functor for sphere with center (0.5,0.5,-0.5) and r = 0.15
 *  When node fails constraint it is removed 
 */
struct RemoveSphere {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t; //silence compiler
    //iterate through nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      //check if node violates sphere constraint
      Point del = n.position() - Point(0.5,0.5,-0.5);  
      //check whether the norm of this difference is less than 0.15
      if (norm(del) < 0.15) {
        //remove node
        g.remove_node(n);
      }
    }
    return;
  }
};

/** Constraint function object that combines two constraints */
template <typename C1, typename C2>
struct combine_two_constraints {
 private:
  C1 c1_; 
  C2 c2_; 
 public:
  //constructor
  combine_two_constraints(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  //call operator
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;    // silence compiler warnings, force is time independent
    //call two functors, then return
    c1_(g, t);
    c2_(g, t);
    return;
  } 
};


template <typename C1, typename C2>
combine_two_constraints<C1,C2> make_combined_constraint (C1 c1, C2 c2) {
  //return functor
  return {c1, c2};
}

/** Constraint function object that combines three constraints */
/*
template <typename C1, typename C2, typename C3>
struct combine_three_constraints {
 private:
  C1 c1_; 
  C2 c2_; 
  C3 c3_; 
 public:
  //constructor
  combine_three_constraints(C1 c1, C2 c2, C3 c3) : c1_(c1), c2_(c2), c3_(c3) {}
  //call operator
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;    // silence compiler warnings, force is time independent
    //call three functors, then return
    c1_(g, t);
    c2_(g, t);
    c3_(g, t);
    return;
  } 
};

template <typename C1, typename C2, typename C3>
combine_three_constraints<C1,C2,C3> make_combined_constraint (C1 c1, C2 c2, C3 c3) {
  //return functor
  return {c1, c2, c3};
}
*/

/** Force function object that combines forces, returning a Point
 *   containing x, y, z coordinates of force
 */
template <typename F1, typename F2>
struct combine_two {
 private:
  F1 f1_; 
  F2 f2_; 
 public:
  //constructor
  combine_two(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
  //call operator
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings, force is time independent
    //combine results from f1_ and f2_ 
    return f1_(n,t) + f2_(n,t);
  } 
};


template <typename F1, typename F2, typename F3>
struct combine_three {
 private:
  F1 f1_; 
  F2 f2_; 
  F3 f3_; 
 public:
  //constructor
  combine_three(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {}
  //call operator
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings, force is time independent
    //combine results from f1_ and f2_ and f3_ 
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  } 

};

template <typename F1, typename F2>
combine_two<F1,F2> make_combined_force (F1 f1, F2 f2) {
  //return functor combine_two
  return {f1, f2};
}

template <typename F1, typename F2, typename F3>
combine_three<F1,F2,F3> make_combined_force (F1 f1, F2 f2, F3 f3) {
  //return functor combine_three
  return {f1, f2, f3};
}

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //silence compiler
    //return point with gravity in z direction
    return Point(0,0, n.value().mass*(-grav));
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //silence compiler
    Point s_force = Point(0,0,0); //define initial spring force
    //add contributions from each edge incident to n
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      //assume first node is spawning node
      Node n1 = e.node1();      
      Node n2 = e.node2();
      //check assumption
      if (e.node2() == n) {
	n2 = e.node1(); //adjacent node
	n1 = e.node2(); //spawning node
      }
      //update spring force with contribution from e
      Point del = n1.position() - n2.position();
      s_force += ((-e.value().K_)*del/e.length())*(e.length() - e.value().L_);      
    }
    return s_force;
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //silence compiler
    //return point of -c*v where c = 1/N and v is velocity of point n
    //unsigned N = n.graph_->num_nodes(); //won't be able to access -- graph_ is private
    //double c = -1.0/double(N);
    double c = -1.0/double(num_nodes); //access global variable
    //double c = -.00001;
    return c*n.value().vel;

  }

};


/** Force function object for HW2 #1. */
struct Problem1Force {
 
 //private:
  //double K_;
  //double L_;

 //public:
  /*Constructor to pass in K and L*/
  //Problem1Force (double K, double L) : K_(K), L_(L) {}

  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;    // silence compiler warnings, force is time independent

    //constrain two corners of cloth
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }

    //define gravity force
    double z_grav = n.value().mass*(-grav);
    //Point g_force = Point(0,0,n.value().mass*(-grav));

    //define spring force
    double x_spring = 0.0; 
    double y_spring = 0.0;
    double z_spring = 0.0;
    //Point s_force = Point(0,0,0);
    //call length on each edge incident to n
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      //assume first node is spawning node
      Node n1 = e.node1();      
      Node n2 = e.node2();
      //check assumption
      if (e.node2() == n) {
	n2 = e.node1(); //adjacent node
	n1 = e.node2(); //spawning node
      }

      double delx = n1.position().x - n2.position().x;
      double dely = n1.position().y - n2.position().y;
      double delz = n1.position().z - n2.position().z;

      x_spring += ((-e.value().K_)*delx/e.length())*(e.length() - e.value().L_);      
      y_spring += ((-e.value().K_)*dely/e.length())*(e.length() - e.value().L_);      
      z_spring += ((-e.value().K_)*delz/e.length())*(e.length() - e.value().L_);      
    }

    return Point(x_spring, y_spring, z_spring + z_grav);
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
  //initial positions are already set from reading in file
  //set initial velocities by calling NodeData constructor on each node in graph
  //set initial masses by updating value().mass for each node in graph

  int N = graph.num_nodes();  
  size_type P000;
  size_type P100;

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value() = NodeData(); //constructor sets all velocities to 0
    n.value().mass = 1.0/double(N);// update mass for each node to be 1/N

    //Find uid's of Point(0,0,0) and Point (1,0,0) - use for FixCorners
    if (n.position() == Point(0,0,0)) {
      P000 = n.index();
    } 
    if (n.position() == Point(1,0,0)) {
      P100 = n.index();
    } 
  }

  //Set spring constant K = 100
  double K = 100.0;
  //Set rest length L to be the initial length of the edges
  //Edge first_edge = *graph.edge_begin(); //choose arbitrary edge
  //double L = first_edge.length(); 

  //#2: populate EdgeData 
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    //e.value() = EdgeData(K, e.length()); //constructor sets up values of K and L
    e.set_value(EdgeData(K, e.length()));  // ensures both entries in adj_ are populated
  }


  // Update global variable num_nodes
  num_nodes = graph.num_nodes();


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
      //double dt = 0.0002;
      double dt = 0.0005;
      double t_start = 0;
      //double t_end = 0.2;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {


        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()), 
 	//   make_combined_constraint(ConstrainPlane(), ConstrainSphere()), FixCorners());

        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()), 
 	                make_combined_constraint(RemoveSphere(), ConstrainPlane()), FixCorners());

        //for RemoveSphere - clear viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map); //for RemoveSphere
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
