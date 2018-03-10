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
static constexpr double K = 100; // spring constant
double damp_const = 0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/* Euclidian distance between x_i and x_j */
double euc_dist(Point& p1, Point& p2) {
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

/* Constaint call set up using virtual operator */
class Constraint {
public:
  virtual void operator()(GraphType &g, double t) {(void) g; (void) t;}
};

/* Returning a combination of all the constraints */
class ComboConstraint:public Constraint {
public:
  ComboConstraint(Constraint *c1, Constraint *c2, Constraint *c3 = nullptr):c1_(c1),
    c2_(c2), c3_(c3) {}
  void operator()(GraphType &g, double t) {
    (*c1_)(g,t);
    (*c2_)(g,t);
    if (c3_ != nullptr) (*c3_)(g,t);
  }
private:
  Constraint *c1_, *c2_, *c3_;
};

/* Corner constraint to keep point(0,0,0) and point(1,0,0) unchanged */
class FixedCornersConstraint:public Constraint {
public:
  void operator()(GraphType &g, double t) {
    (void) t;
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node n = *nit;
      if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        n.value().vel = Point(0,0,0);
      }
    }
  }
};

/* Floor contraint to keep graph from sliding past specific z position */
class ZConstraint:public Constraint {
public:
  ZConstraint(double z):floor_(z) {};
  void operator()(GraphType &g, double t) {
    (void) t;
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node n = *nit;
      if (n.position().z < floor_) {
        n.position().z = floor_;
        n.value().vel.z = 0; // z component of node velocity
      }
    }
  }
private:
  double floor_;
};

/* Graph will wrap around a spherical ball constraint */
class SphereConstraint:public Constraint {
public:
  SphereConstraint(Point center, double radius):center_(center), radius_(radius) {}
  void operator()(GraphType &g, double t) {
    (void) t;
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node n = *nit;
      // if node violates constraint
      if (euc_dist(n.position(), center_) < radius_) {
        n.position() = (n.position() - center_) * (radius_ / euc_dist(n.position(),
          center_)) + center_;
        Point R_i = (n.position() - center_) / euc_dist(n.position(), center_);
        n.value().vel -= dot(n.value().vel,R_i) * R_i;
      }
    }
  }
private:
  Point center_;
  double radius_;
};

/* Remove nodes and edges from center of graph */
class RemoveConstraint:public Constraint {
public:
  RemoveConstraint(Point center, double radius):center_(center), radius_(radius) {}
  void operator()(GraphType &g, double t) {
    (void) t;
    for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
      Node n = *nit;
      if (euc_dist(n.position(), center_) < radius_) {
        g.remove_node(n);
      }
    }
  }
private:
  Point center_;
  double radius_;
};

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
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

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
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE

    (void) n; (void) t; (void) grav;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    Point grav_force = n.value().mass * Point(0,0,-1*grav);

    Point spring_force = Point(0,0,0);
    Point x_i = n.position();
    for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
      Edge e = *iit;
      Node n2 = e.node2();
      Point x_j = n2.position();
      double dist = euc_dist(x_i, x_j);
      double L = e.value();
      Point f = -1 * K * (x_i - x_j) / dist * (dist - L);
      spring_force += f;
    }
    return spring_force + grav_force;
  }
};

/* Return the combination of 2 forces */
template <typename Force1, typename Force2>
struct ComboForce {
  ComboForce(Force1 f1, Force2 f2):f1_(f1), f2_(f2) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
  private:
    Force1 f1_; Force2 f2_;
};

/* Return combination of two forces by calling ComboForce struct */
template <typename Force1, typename Force2>
ComboForce<Force1, Force2> make_combined_force(Force1 f1, Force2 f2) {
  return ComboForce<Force1, Force2>(f1,f2);
}

/* Return combination of 3 forces */
template <typename Force1, typename Force2, typename Force3>
ComboForce<ComboForce<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3) {
  return ComboForce<ComboForce<Force1, Force2>, Force3>(ComboForce<ComboForce<Force1,
    Force2>, Force3>(ComboForce<Force1, Force2>(f1,f2), f3));
}

/* Calculation of the gravity force */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    return n.value().mass * grav * Point(0,0,-1);
  }
};

/* Calculation of mass spring force */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    Point x_i = n.position();
    Point result;
    for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
      Point x_j = (*iit).node2().position();
      double dist_i_j = euc_dist(x_i, x_j);
      double L = (*iit).value();
      result += (-1 * K) * ((x_i-x_j) / dist_i_j) * (dist_i_j - L);
    }
    return result;
  }
};

/* Calculation of damping force */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    return -1 * damp_const * n.value().vel;
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

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  /* Initialize damp constant */
  damp_const = 1.0 / graph.num_edges();

  /* Initialize node masses */
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit)
    (*nit).value().mass = 1.0 / graph.num_nodes();

  /* Initialize edge values (val = distance b/w n1 and n2)*/
  for (auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit)
    (*eit).value() = euc_dist((*eit).node1().position(), (*eit).node2().position());

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  /* Add node and edges to graph */
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.0006;
      double t_start = 0;
      double t_end = 5.0;

      /* Call to the 4 different contraints and combination of them */
      SphereConstraint sc(Point(0.5,0.5,-0.5),0.15);
      FixedCornersConstraint fc;
      ZConstraint zc(-0.75);
      RemoveConstraint rc(Point(0.5,0.5,-0.5),0.15);
      //ComboConstraint cc = ComboConstraint(&zc, &sc, &fc);
      ComboConstraint cc = ComboConstraint(&zc, &rc, &fc);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),
          MassSpringForce(), DampingForce()));
        cc(graph, t); //correct for voilated constraints once symp_e_step is finished

        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        viewer.set_label(t);

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
