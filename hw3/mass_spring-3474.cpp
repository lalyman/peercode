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

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  double damping;
  NodeData() : vel(0), mass(1), damping(0.01) {}
};

struct EdgeData {
  double K;       //< Edge spring constant
  double length;     //< edge L
  EdgeData() : K(0), length(0) {}
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
double symp_euler_step(G& g, double t, double dt, F force, C constrain){
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  constrain(g, t);
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
    (void) t;
    if(n.position () == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point();
    }
    auto f_spring = Point(0, 0, 0);
    auto pni = n.position();
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      f_spring += -(*it).value().K*(pni - (*it).node2().position())*
        (norm_2(pni - (*it).node2().position()) - (*it).value().length)/
        norm_2(pni - (*it).node2().position());
    }

    f_spring += n.value().mass * Point(0, 0, -grav);
    return f_spring;
  }
};

//Fix two points given
struct ConstConstraint {
  template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
    (void)t;
		for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = (*it);
      if(n.position () == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
        n.value().vel = Point(0, 0, 0);
      }
    }
	}
};

//constraint on a given plane
struct PlaneConstraint {
  PlaneConstraint(double z): z_(z) {}
  template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
    (void)t;
		for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = (*it);
      if(inner_prod(n.position(), Point(0,0,1)) < z_){
          n.position().z = -0.75;
          n.value().vel.z = 0;
        }
    }
	}
  double z_;
};

//constraint on a given sphere
struct SphereConstraint {
	SphereConstraint(Point c, double r): center_(c), radius_(r) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = (*it);
      if(norm_2(n.position() - center_) < radius_){
          Point rescaled = n.position() - center_;
          rescaled = rescaled/norm_2(rescaled);
          n.position() = center_ + rescaled * radius_;
          n.value().vel = n.value().vel - rescaled*dot(n.value().vel, rescaled);
      }
    }
	}
	Point center_;
	double radius_;
};

//constraint when removing
struct SphereConstraint_remove {
	SphereConstraint_remove(Point c, double r): center_(c), radius_(r) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
    (void)t;
    for(unsigned i = 0; i < g.num_nodes();){
      if(norm_2(g.node(i).position() - center_) < radius_) {
        i = g.remove_node(g.node(i));
      }
      else{
        i++;
      }
    }
	}
	Point center_;
	double radius_;
};

//functor combines two constraints
template<typename C1, typename C2>
struct CombinedConstraint {
	CombinedConstraint(C1 c1, C2 c2):c1_(c1), c2_(c2) {}
	template<typename GRAPH>
	void operator()(GRAPH& g, double t) {
		c1_(g, t);
		c2_(g, t);
		(void) t;
	}
	C1 c1_;
	C2 c2_;
};

//Function that returns the combined constraint of two constraints
template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
	return CombinedConstraint<C1, C2>(c1, c2);
}

// Function that returns the combined constraint of three constraints
template<typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
	return CombinedConstraint<CombinedConstraint<C1, C2>, C3>(make_combined_constraint(c1, c2), c3);
}

struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return n.value().mass * Point(0, 0, -grav);
  }
};


struct DampingForce{
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return -n.value().damping*n.value().vel;
  }
};


struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    if(n.position () == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point();
    }
    auto f_spring = Point(0, 0, 0);
    auto pni = n.position();
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      f_spring += -(*it).value().K*(pni - (*it).node2().position())*
        (norm_2(pni - (*it).node2().position()) - (*it).value().length)/
        norm_2(pni - (*it).node2().position());
    }
    return f_spring;
  }
};

//functor that combine two forces
template <typename F1, typename F2>
struct combine_forces{
  F1 force1;
  F2 force2;
  combine_forces(F1 force1_, F2 force2_) : force1(force1_), force2(force2_){
  }
  template <typename NODE>
  Point operator()(NODE n, double t){
    return force1(n, t) + force2(n, t);
  }
};

// Function that returns the combined force of two forces
template<typename F1, typename F2>
combine_forces<F1, F2> make_combined_forces(F1 f1, F2 f2) {
	return combine_forces<F1, F2>(f1, f2);
}

// Function that returns the combined force of three forces
template<typename F1, typename F2, typename F3>
combine_forces<combine_forces<F1, F2>, F3> make_combined_forces(F1 f1, F2 f2, F3 f3) {
	return combine_forces<combine_forces<F1, F2>, F3>(make_combined_forces(f1, f2), f3);
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


  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    (*it).value().mass = 1.0/graph.size();
    (*it).value().vel = Point(0,0,0);
    (*it).value().damping = 1.0/graph.size();
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
    (*it).value().K = 100;
    (*it).value().length = (*it).length();
    graph.get_counterpart((*it)).value().K = 100;
    graph.get_counterpart((*it)).value().length = (*it).length();

  }

  // Print out the stats

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
      //double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        auto myforce = make_combined_forces(GravityForce(), MassSpringForce(), DampingForce());
        auto myconstraint = make_combined_constraint(ConstConstraint(),PlaneConstraint(-0.75),
          SphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        //auto removeconstraint = make_combined_constraint(ConstConstraint(),
        //  SphereConstraint_remove(Point(0.5, 0.5, -0.5), 0.15));
        symp_euler_step(graph, t, dt, myforce, myconstraint);

        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        //constrain2(graph, Point(0.5, 0.5, -0.5), t, 0.15);
        // Update viewer with nodes' new positions
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
  //interrupt_sim_thread = true;
  //sim_thread.join();

  return 0;
}
