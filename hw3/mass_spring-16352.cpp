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
  //HW2
  NodeData(Point vel_, double mass_) : vel(vel_), mass(mass_) {}
};
//HW2
struct EdgeData {
  double K;//spring constant
  double L;//rest-length
  EdgeData() : K(0.0), L(1.0) {}
  EdgeData(double K_, double L_) : K(K_), L(L_) {}
};

// Define the Graph type
//using GraphType = Graph<NodeData>;
using GraphType = Graph<NodeData,EdgeData>;
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt

  //constraint on two points
  if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {}
  else {
    n.position() += n.value().vel * dt;
    }
  }
 
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {}
    else {
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }
  return t + dt;
}

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt

  //constraint on two points
  if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {}
  else {
    n.position() += n.value().vel * dt;
    }
  }
 
  //Constraints HW2
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    constraint(n,t);
  }
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {}
    else {
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
  Problem1Force() : K(1.0), L(1.0) {}
  Problem1Force(double K_, double L_) : K(K_), L(L_) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point (0,0,0);
    //return Point(0);
    //update L value

    Point force (0);
    Point distance (0);
    //Calculate spring force
    for( auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
        //HW2-2
        K=(*ei).value().K;
        L=(*ei).value().L;
        //
        Node other_node = (*ei).node1();
        if ( n == (*ei).node1())
            other_node = (*ei).node2();  // Get other node of the incident edge
        distance = n.position()-other_node.position();
        force += -K*distance*(norm(distance)-L)/norm(distance);
    }
    //adding gravitational force
    force += n.value().mass*( Point (0,0,-grav) );
    return force;
  }
private:
    double K;//spring constant
    double L;//spring rest lenght
};

struct NullForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t;     // silence compiler warnings
    return Point (0);
  }
};
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    return n.value().mass*( Point (0,0,-grav) );
  }
};
struct MassSpringForce {
  MassSpringForce() : K(1.0), L(1.0) {}
  MassSpringForce(double K_, double L_) : K(K_), L(L_) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    Point force (0);
    Point distance (0);
    //Calculate spring force
    for( auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
        K=(*ei).value().K;
        L=(*ei).value().L;
        Node other_node = (*ei).node1();
        if ( n == (*ei).node1())
            other_node = (*ei).node2();  // Get other node of the incident edge
        distance = n.position()-other_node.position();
        force += -K*distance*(norm(distance)-L)/norm(distance);
    }
    return force;
  }
private:
    double K;//spring constant
    double L;//spring rest lenght
};
struct DampingForce {
  DampingForce() : c(1.0) {}
  DampingForce(double c_) : c(c_) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point (0,0,0);
    return -c*n.value().vel;
  }
private:
    double c;//damping constant
};

    
//Function object for combined force
template <class F1, class F2, class F3>
struct CombinedForce {
    CombinedForce() {}
    CombinedForce(F1 f1_, F2 f2_, F3 f3_) : f1(f1_), f2(f2_), f3(f3_) {}

    Point operator()(Node n, double t) {
        return f1(n,t) + f2(n,t) + f3(n,t);
    }
private :
    F1 f1;
    F2 f2;
    F3 f3;
};

template <class F1, class F2, class F3>
CombinedForce<F1,F2,F3> make_combined_force(const F1& f1_, const F2& f2_, const F3& f3_)
{
    return CombinedForce<F1,F2,F3>(f1_, f2_, f3_);
}

template <class F1, class F2>
CombinedForce<F1,F2,NullForce> make_combined_force(const F1& f1_, const F2& f2_)
{
    return CombinedForce<F1,F2,NullForce>(f1_, f2_, NullForce());
}

struct ZPlaneConstraint {
  ZPlaneConstraint() : z(-0.75) {}
  ZPlaneConstraint(double z_) : z(z_) {}
  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) n; (void) t;     // silence compiler warnings
    if ( inner_prod(n.position(),Point(0,0,1)) < z ) {
        n.position().z = z;
        n.value().vel.z = 0.0;
    }
  }
private:
  double z;
};

struct SphereConstraint {
    SphereConstraint() : center(Point(0.5,0.5,-0.5)), radius(0.15) {}
    SphereConstraint(Point center_, double radius_) : center(center_), radius(radius_) {}
  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) n; (void) t;     // silence compiler warnings
    if ( norm(n.position()-center) < radius ) {
        Point unitC2P = (n.position() - center) / norm(n.position() - center);
        n.position() = radius * unitC2P + center;
        n.value().vel = n.value().vel - inner_prod(n.value().vel,unitC2P)*unitC2P;
    }
  }
private:
    Point center;
    double radius; 
};

//Function object for combined constraint 
template <class C1, class C2>
struct CombinedConstraint {
    CombinedConstraint() {}
    CombinedConstraint(C1 c1_, C2 c2_) : c1(c1_), c2(c2_) {}

    void operator()(Node n, double t) {
        c1(n,t); c2(n,t);
    }
private :
    C1 c1;
    C2 c2;
};

template <class C1, class C2>
CombinedConstraint<C1,C2> make_combined_constraint(const C1& c1_, const C2& c2_)
{
    return CombinedConstraint<C1,C2>(c1_, c2_);
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

  NodeData tempNode { Point (0), 1.0/graph.num_nodes() };
  for ( auto ni = graph.node_begin(); ni != graph.node_end(); ++ni)
  {
      (*ni).value() = tempNode;
  }
  double L=(*(graph.edge_begin())).length();

  //HW2-2
  for ( auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei)
  {
      (*ei).value() = EdgeData { 100.0, (*ei).length() };
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
      //double t_end = 5.0;
      double t_end = 3.0;

      Problem1Force problem1Force {100.0,L};//HW2-1
      DampingForce dampingForce {1.0/graph.num_nodes()};//HW2-3
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force());
        //symp_euler_step(graph, t, dt, problem1Force);
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), dampingForce));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), dampingForce), ZPlaneConstraint());
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), dampingForce),
                make_combined_constraint(SphereConstraint(),ZPlaneConstraint()));

        // Update viewer with nodes' new positions
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
