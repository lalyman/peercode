/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <chrono>
#include <cmath>
#include <fstream>
#include <memory>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Global for HW2 #1
// static double L = 0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double L; // Spring rest length.
  double K; // Spring constant
  EdgeData() : L(1.0), K(100.0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/// Change a graph's nodes according to a step of the symplectic Euler
///   method with the given node force.
/// @param[in,out] g      Graph
/// @param[in]     t      The current time (useful for time-dependent forces)
/// @param[in]     dt     The time step
/// @param[in]     force  Function object defining the force per node
/// @return the next time step (usually @a t + @a dt)
//
/// @tparam G::node_value_type supports ???????? YOU CHOOSE
/// @tparam F is a function object called as @a force(n, @a t),
///         where n is a node of the graph and @a t is the current time.
///         @a force must return a Point representing the force vector on
///         Node n at time @a t.
/// @tparam C is a function object called as @a constraint(@a g, @a t),
///         where @a g is the graph and @a t is the current time. @a
///         constraint will modify node values when called.
///
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraint.
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


///** Force function object for HW2 #1. */
//struct Problem1Force {
//  /** Return the force applying to @a n at time @a t.
//   *
//   * For HW2 #1, this is a combination of mass-spring force and gravity,
//   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
//   * model that by returning a zero-valued force. */
//  template <typename NODE>
//  Point operator()(NODE n, double t) {
//
//    // (void) t;
//    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
//      return Point(0, 0, 0);
//    }
//
//    Point force = n.value().mass * Point(0, 0, -grav);
//    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
//      double edgelength = it->length();
//      force += (
//        -it->value().K * (n.position() - it->endpoint(n).position())
//           * (edgelength - it->value().L) / edgelength);
//    }
//    return force;
//  }
//};


// Container class that can collect either Forces or Constraints and
// apply them in succession to a node/graph.
template <class T>
class PhysicalSum {

  public:

    // Since forces (constraints) can be of different types, but are
    // all derived from the same abstract base class, our container
    // must store pointers to the forces (constraints).
    void add_object(std::shared_ptr<T> t) {
      objects.push_back(t);
    }

    // Applies each force (constraint) in succession. The return type
    // is Point, but this is irrelevant for the constraints.
    template <class Subject>
    Point operator()(Subject &s, double t) {

      Point total = Point(0, 0, 0);
      for (auto it = objects.begin(); it < objects.end(); ++it) {
        total += (**it)(s, t);
      }
      return total;
    }

  private:

    // Vector of pointers to Forces/Constraints.
    std::vector<std::shared_ptr<T>> objects;
};


// Abstract base class for representing forces.
class Force {
  
  public:
    
    virtual Point operator()(Node &n, double t) = 0;
};


// Applies a downward gravitational force.
class GravityForce : public Force {
  
  public:

    Point operator()(Node &n, double t) {
      (void) t;
      return n.value().mass * Point(0, 0, -grav);
    }
};


// Edges connecting nodes to their neighbors are treated as ideal
// springs, modeled by Hooke's law. Each spring can have its own spring
// constant (set in the edge's value).
class MassSpringForce : public Force {

  public:

    Point operator()(Node &n, double t) {

      (void) t;
      Point force = Point(0, 0, 0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        double edgelength = it->length();
        force += (
          -it->value().K * (n.position() - it->endpoint(n).position())
             * (edgelength - it->value().L) / edgelength);
      }
      return force;
    }
};


// Damping force, proportional to velocity.
class DampingForce : public Force {

  public:

    DampingForce(double c) : C(c) {
    }

    Point operator()(Node &n, double t) {
      (void) t;
      return -C * n.value().vel;
    }

  private:

    double C;
};


// Abstract base class for constraints.
class Constraint {

  public:

    // operator() must return a Point to be compatible with the
    // PhysicalSum container.
    virtual Point operator()(GraphType &g, double t) = 0;
};

// Infinite plane that prevents the fabric from passing through.
class PlaneConstraint : public Constraint {

  public:

    // Plane computations are completely general. That is, by passing
    // in a normal vector and an offset, we can create arbitrary planes.
    PlaneConstraint(Point normal, double offset)
      : N(normal/norm(normal)), c0(offset/norm(normal)) {
    }

    // When the constraint is applied, a Node violating the constraint
    // is adjusted to the closest point on the plane and the component
    // of its velocity normal to the plane is set to zero.
    Point operator()(GraphType &g, double t) {

      (void) t;
      for (auto it = g.node_begin(); it < g.node_end(); ++it) {
        if (dot(N, it->position()) < c0) {
          it->position() -= (dot(it->position(), N) - c0) * N;
          it->value().vel -= dot(it->value().vel, N) * N;
        }
      }
      return Point(0, 0, 0);
    }

  private:

    Point N;
    double c0;
};

// Solid sphere that can handle collisions with the fabric.
class SphereConstraint : public Constraint {

  public:

    SphereConstraint(Point center, double radius) : C(center), r(radius) {
    }

    Point operator()(GraphType &g, double t) {

      (void) t;
      for (auto it = g.node_begin(); it < g.node_end(); ++it) {
        Point Ri = it->position() - C;
        if (norm(Ri) < r) {
          Point Ri_hat = Ri / norm(Ri);
          it->position() += (r - norm(Ri)) * Ri_hat;
          it->value().vel -= dot(it->value().vel, Ri_hat) * Ri_hat;
        }
      }
      return Point(0, 0, 0);
    }
  
  protected:

    Point C;
    double r;
};

/// Solid sphere that causes the fabric to tear upon collision.
class PowerfulSphere : public SphereConstraint {

  public:

    using SphereConstraint::SphereConstraint;

    Point operator()(GraphType &g, double t) {

      (void) t;
      for (auto it = g.node_begin(); it < g.node_end(); ) {
        Point Ri = it->position() - C;
        if (norm(Ri) < r) {
          g.remove_node(*it++);
        } else {
          ++it;
        }
      }
      return Point(0, 0, 0);
    }
};

// Constraint that fixes an arbitrary number of nodes in place. 
class PinConstraint : public Constraint {

  public:

    // Takes in a time-step in order to be able to reset position.
    PinConstraint(double dt) : dt(dt) {
    }

    void add_pin(Node &n) {
      pin_nodes.push_back(n);
    }

    // Constraint is applied in O(|pin_nodes|) (constant time for each
    // node in pin_nodes).
    Point operator()(GraphType &g, double t) {

      (void) t, (void) g;
      for (auto it = pin_nodes.begin(); it < pin_nodes.end(); ++it) {
        // Must reset the position before setting the velocity to 0.
        it->position() -= it->value().vel * dt;
        it->value().vel = Point(0, 0, 0);
      }
      return Point(0, 0, 0);
    }

  private:

    double dt;
    std::vector<Node> pin_nodes;
};

/// Combines two forces or constraints.
/// @param[out] combined Container object that evaluates both
///   forces/constraints when called.
/// @param[in] f The first force/constraint.
/// @param[in] g The second force/constraint.
///
/// @tparam T Either Force or Constraint.
/// @tparam F1 Type of the first force/constraint. Must be a derived
///   class of T.
/// @tparam F2 Type of the second force/constraint. Must be a derived
///   class of T.
///
template <class T, class F1, class F2>
void combine(PhysicalSum<T> &combined, F1 f, F2 g) {
  
  combined.add_object(std::shared_ptr<T>(new F1(f)));
  combined.add_object(std::shared_ptr<T>(new F2(g)));
}

/// Combines three forces or constraints.
/// @param[out] combined Container object that evaluates all three
///   forces/constraints when called.
/// @param[in] f The first force/constraint.
/// @param[in] g The second force/constraint.
/// @param[in] h The third force/constraint.
///
/// @tparam T Either Force or Constraint.
/// @tparam F1 Type of the first force/constraint. Must be a derived
///   class of T.
/// @tparam F2 Type of the second force/constraint. Must be a derived
///   class of T.
/// @tparam F3 Type of the third force/constraint. Must be a derived
///   class of T.
///
template <class T, class F1, class F2, class F3>
void combine(PhysicalSum<T> &combined, F1 f, F2 g, F3 h) {
  
  combined.add_object(std::shared_ptr<T>(new F1(f)));
  combined.add_object(std::shared_ptr<T>(new F2(g)));
  combined.add_object(std::shared_ptr<T>(new F3(h)));
}


/// Find the node with the minimum euclidean distance to a point.
/// @param g  The graph of nodes to search.
/// @param point  The point to use as the query.
/// @return An iterator to the node of @a g with the minimun Eucliean
///           distance to @a point.
///           graph.node_end() if graph.num_nodes() == 0.
/// 
/// @post For all i, 0 <= i < graph.num_nodes(),
///          norm(point - *result) <= norm(point - g.node(i).position())
/// 
Node &nearest_node(GraphType &g, const Point &point) {

  auto it = std::min_element(
    g.node_begin(), g.node_end(),
    // Using lambda function instead of functor to keep things
    // succinct. Check my implementation of a full functor in
    // subgraph.cpp to see that I also know how to do that.
    [point](const Node &u, const Node &v) {
      return norm(u.position() - point) < norm(v.position() - point); } );
  return *it;
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

  // Initialize mass.
  for (auto it = graph.node_begin(); it < graph.node_end(); ++it) {
    it->value().mass = 1.0 / graph.num_nodes();
  }

  // Initialize spring rest length.
  for (auto it = graph.edge_begin(); it < graph.edge_end(); ++it) {
    it->value().L = it->length();
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
      double dt = 0.001/2;
      double t_start = 0;
      double t_end = 5.0;

      // Create a constraint the pins the corners.
      PinConstraint pin(dt);
      pin.add_pin(nearest_node(graph, Point(0, 0, 0)));
      pin.add_pin(nearest_node(graph, Point(1, 0, 0)));

      // Combine relevant forces and desired constraints.
      PhysicalSum<Force> combined_force;
      PhysicalSum<Constraint> combined_constraint;
      combine(
        combined_force, GravityForce(), MassSpringForce(),
        DampingForce(1./graph.num_nodes()));
      combine(
        combined_constraint, pin, PlaneConstraint(Point(0, 0, 1), -1),
        //SphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        PowerfulSphere(Point(0.5, 0.5, -0.5), 0.15));

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(
          graph, t, dt, combined_force, combined_constraint);

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
