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
  double len;   //< Rest length
  double K;     //< spring constant
  EdgeData() : len(1.), K(100.) {}
};

// Define the Graph type
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // hard-code constraints for HW2 #3
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
    //   continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // apply constraint
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // hard-code constraints for HW2 #3
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
    //   continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
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
  Point operator()(NODE n, double) {

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    double K = 100.;
    double len = 0.01; // use same rest length for all edges

    Point f_spring(0.);
    // Accumulate the spring force using the incident iterator
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto dx = n.position() - (*it).node2().position();
      auto norm_dx = norm(dx);
      f_spring += -K*dx/norm_dx*(norm_dx - len); 
    }

    Point f_grav(0.,0.,-n.value().mass*grav);

    return f_spring + f_grav;
  }
};



/** Force function object for HW2 #2. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * Here we use the EdgeData struct of each edge for K and rest length.
   * Otherwise the same with HW2 #1 */
  template <typename NODE>
  Point operator()(NODE n, double) {

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    Point f_spring(0.);
    // Accumulate the spring force using the incident iterator
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      auto dx = n.position() - e.node2().position();
      auto norm_dx = norm(dx);
      f_spring += -e.value().K*dx/norm_dx*(norm_dx - e.value().len);
    }

    Point f_grav(0.,0.,-n.value().mass*grav);

    return f_spring + f_grav;
  }
};


/** Gravity force function object for HW2 #3. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double) {
    return Point(0.,0.,-n.value().mass*grav);
  }
};

/** Spring force function object for HW2 #3. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double) {
    Point f_spring(0.);
    // Accumulate the spring force using the incident iterator
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      auto dx = n.position() - e.node2().position();
      auto norm_dx = norm(dx);
      f_spring += -e.value().K*dx/norm_dx*(norm_dx - e.value().len);
    }
    return f_spring;
  }
};

/** Damping force function object for HW2 #3.
 * @a c the damping constant which equals to graph size. */
struct DampingForce {
  double c;
  template <typename G>
  DampingForce(const G& g) : c(1./g.size()) {}

  template <typename NODE>
  Point operator()(NODE n, double) {
    return -c*n.value().vel;
  }
};


/** function object to combine forces for HW2 #3 */
template <typename F1, typename F2>
struct CombinedForce {
  F1 f1_;
  F2 f2_;
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
  
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
};

/** Return the sum of two forces
 * @param[in] f1,f2 forces to be combined
 * @return the combined force @a f1 + @a f2
 *
 * @tparam F1,F2 function objects called as e.g. @a f1(n, t),
 *               where n is a node of the graph and t is the current time.
 *               @a f1 and @a f2 must return a Point representing
 *               the force vector on Node n at time @a t.
 *
 * It uses the function object CombinedForce and deduces
 * the template parameters of CombinedForce from the types of arguments.
 */
template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1,F2>(f1, f2);
}

/** Return the sum of three forces
 * @param[in] f1,f2,f3 forces to be combined
 * @return the combined force @a f1 + @a f2 + @a f3
 *
 * @tparam F1,F2,F3 function objects called as e.g. @a f1(n, t),
 *                  where n is a node of the graph and t is the current time.
 *                  @a f1, @a f2, and @a f3 must return a Point representing
 *                  the force vector on Node n at time @a t.
 *
 * It uses the function object CombinedForce recursively and deduces
 * the template parameters of CombinedForce from the types of arguments.
 */
template <typename F1, typename F2, typename F3>
CombinedForce<F1,CombinedForce<F2,F3> > make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombinedForce<F1,CombinedForce<F2,F3> >(f1, make_combined_force(f2,f3));
}

/** Constant node constraint for HW2 #4. 
 * Two nodes @a n1 and @a n2 which are initially at 
 * positions @a p1_ and @a p2_, respectively, will stay fixed.
 * We assume that such @a n1 and @a n2 exist in the graph. */
template <typename NODE>
struct ConstantConstraint {
  Point p1_, p2_;
  NODE n1, n2;

  template <typename G>
  ConstantConstraint(G& g, const Point& p1, const Point& p2) : p1_(p1), p2_(p2) {
    int found = 0;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position() == p1_) { n1 = n; ++found; } 
      if (n.position() == p2_) { n2 = n; ++found; } 
      if (found == 2) break;
    }
  }

  template <typename G>
  void operator()(G&, double) {
    n1.position() = p1_;
    n1.value().vel = Point(0,0,0);
    n2.position() = p2_;
    n2.value().vel = Point(0,0,0);
  }
};

/** z-Plane constraint for HW2 #4. 
 * @a z the z-coordiante of the z-Plane we use for the constraint. */
struct ZPlaneConstraint {
  double z;
  ZPlaneConstraint(double z_) : z(z_) {}

  template <typename G>
  void operator()(G& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (dot(n.position(), Point(0,0,1)) < z) {
        n.position().z = z; // the nearest point to the plane.
        n.value().vel.z = 0.; // z-component of velocity set to be zero.
      }
    }
  }
};

/** sphere constraints for HW2 #4. 
 * @a c the center of the sphere.
 * @a r the radius of the sphere. */
struct SphereConstraint {
  Point c;
  double r;
  SphereConstraint(const Point& c_, double r_) : c(c_), r(r_) {}

  template <typename G>
  void operator()(G& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      auto x_c = n.position() - c; // (x-c)
      auto normx_c = norm(x_c); // |x-c|
      if (normx_c < r) {
        auto R_x = x_c/normx_c;
        // the nearest point on the surface of the sphere.
        n.position() = R_x*r + c; 
        // velocity normal to the sphere set to be zero.
        n.value().vel -= dot(n.value().vel, R_x)*R_x;
      }
    }
  }
};

/** The constraint which removes nodes hitting a sphere, for HW2 #5. 
 * @a c the center of the sphere.
 * @a r the radius of the sphere. */
struct RemoveSphereConstraint {
  Point c;
  double r;
  RemoveSphereConstraint(const Point& c_, double r_) : c(c_), r(r_) {}

  template <typename G>
  void operator()(G& g, double) {
    auto it = g.node_begin();
    while (it != g.node_end()) {
      auto n = *it;
      if (norm(n.position() - c) < r) it = g.remove_node(it);
      else ++it;
      // std::cout << (*it).index() << std::endl;
    }
  }
};


/** function object to combine constraints for HW2 #4 */
template <typename C1, typename C2>
struct CombinedConstraints {
  C1 c1_;
  C2 c2_;
  CombinedConstraints(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  
  template <typename G>
  void operator()(G& g, double t) {
    c1_(g,t);
    c2_(g,t);
  }
};

/** Return the combined constraints
 * @param[in] c1,c2 constraints to be combined
 * @return the combined constraint which enforces both @a c1 and @a c2.
 *
 * @tparam C1,C2 function objects called as e.g. @a c1(g, t),
 *               where g is a graph and t is the current time.
 *
 * @pre @a c1 and @a c2 are compatible constraints. 
 *
 * It uses the function object CombinedConstraints and deduces
 * the template parameters of CombinedConstraints from the types of arguments.
 */
template <typename C1, typename C2>
CombinedConstraints<C1,C2> make_combined_constraints(C1 c1, C2 c2) {
  return CombinedConstraints<C1,C2>(c1, c2);
}

/** Return the combined constraints
 * @param[in] c1,c2,c3 constraints to be combined
 * @return the combined constraint which enforces all of @a c1, @a c2, and @a c3.
 *
 * @tparam C1,C2,C3 function objects called as e.g. @a c1(g, t),
 *                  where g is a graph and t is the current time.
 *
 * @pre @a c1, @a c2, and @a c3 are compatible constraints. 
 *
 * It uses the function object CombinedConstraints recursively and deduces
 * the template parameters of CombinedConstraints from the types of arguments.
 */
template <typename C1, typename C2, typename C3>
CombinedConstraints<C1,CombinedConstraints<C2,C3> > 
    make_combined_constraints(C1 c1, C2 c2, C3 c3) {
        return CombinedConstraints<C1,CombinedConstraints<C2,C3> >
               (c1, make_combined_constraints(c2,c3));
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

  // HW2 #1 mass initialization
  // Note v^{0} are initilized to be 0 when nodes are created.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().mass /= graph.num_nodes();
  }

  // HW2 #2 rest length initialization
  // Note spring constants are already initialized to be 100.
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    e.value().len  = norm(e.node1().position() - e.node2().position());
  }

  // HW2 #3 Forces
  auto f_combined = make_combined_force(GravityForce(),MassSpringForce(),
                                            DampingForce(graph));

  // HW2 #4 Constraints
  // ConstantConstraint<Node> c_constant(graph, Point(0,0,0), Point(1,0,0));
  // ZPlaneConstraint c_plane(-0.75);
  // SphereConstraint c_sphere(Point(0.5,0.5,-0.5), 0.15);
  // auto c_combined = make_combined_constraints(c_constant,c_plane,c_sphere);

  // HW2 #5 Constraints
  ConstantConstraint<Node> c_constant(graph, Point(0,0,0), Point(1,0,0));
  ZPlaneConstraint c_plane(-0.75);
  RemoveSphereConstraint c_rsphere(Point(0.5,0.5,-0.5), 0.15);
  auto c_combined = make_combined_constraints(c_constant,c_plane,c_rsphere);

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
        //std::cout << "t = " << t << std::endl;
        // symp_euler_step(graph, t, dt, Problem1Force()); // for HW2 #1
        // symp_euler_step(graph, t, dt, Problem2Force()); // for HW2 #2
        // symp_euler_step(graph, t, dt, f_combined);  // for HW2 #3
        // symp_euler_step(graph, t, dt, f_combined, c_combined);  // for HW2 #4
        symp_euler_step(graph, t, dt, f_combined, c_combined);  // for HW2 #5

        // Update viewer with nodes' new positions
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);

        // re-draw for HW2 #5
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
