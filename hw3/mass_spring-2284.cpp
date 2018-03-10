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
  Point vel;            //< Node velocity, default value of 0
  double mass;          //< Node mass, default value of 1
  double dampCoeff;     //< Node damping coefficient, default value of 0
  NodeData() : vel(0), mass(1), dampCoeff(0) {}
};

/** Custom structure of data to store with Edges. */
struct EdgeData {
  double restL;         //< Edge rest length, default value of 1
  double sprCoeff;      //< Edge spring coefficient, default value of 100
  EdgeData() : restL(1), sprCoeff(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = typename GraphType::size_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint    Function object defining the constraints
 *                              applied to the graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type stores Point vel, double mass, and double
 *                            dampCoeff, each of which can be set for each node
 *                            individually if necessary.
 * @tparam G::edge_value_type stores double restL and double sprCoeff, each of
 *                            each of which can be set for each edge
 *                            individually if necessary.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t), where g is
 *           the graph itself and @a t is the current time. The method scans
 *           for and corrects node positions and/or velocities.
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

  // Apply constraints as necessary
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {

    auto n = *it;
    // Update the position of the node according to its velocity
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//
// CONSTRAINTS
//

/** Parent class for constraint objects.
 */
class GeneralizedConstraint {
  public:
  /** Virtual method that is to be overwritten in child classes
   *
   * Enforces class design for users who wish to create their own constraint
   * objects, which can be written as children who inherit from this class.
   * Relies on GraphType, as specified at the top fo this file. Default virtual
   * method provides no constraint, meaning that no change is to be expected.
   * This is useful in applications below. Additionally, requiring users to
   * inherit publicly from this class ensures that operator() will be
   * implemented as a valid constraint (trivial constraint, in this case).
   */
  virtual void operator()(GraphType& g, double t) {
    (void) g; (void) t;         // Quiet compiler warnings
    return;
  }
};

/** Fixed point constraints for arbitrary vector of Points.
 *
 * Inherits from parent class GeneralizedConstraint. Overrides operator() and
 * stores information about which points are to fixeid.
 */
class PointsConstraint : public GeneralizedConstraint {
  public:
  /**  Overload of operator(), to return constraint at fixed points. 
   *
   * Fixes Nodes at certain Points, as specified at the moment functor is
   * constructed.
   *
   * @param[in/out]     g       Graph for which Nodes are to be checked
   * @param[in]         t       Time at which check is occurring
   *
   * @tparam    GRAPH   Requirements on Graph are given in operator() doc
   *
   * Runs in O(fixedPoints_.size()), which is assumed to be small
   * */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t; (void) g;         // Quiet compiler warning
    for (unsigned ii = 0; ii != fixedPoints_.size(); ++ii) {
      fixedNodes_[ii].position() = fixedPoints_[ii];
      fixedNodes_[ii].value().vel = Point(0,0,0);
    }
  }

  /** Constructor, for building list of indices of fixed points.
   *
   * @param[in] g       Graph on which constraint will act
   * @param[in] fix     Vector of Points that will be held fixed
   *
   * @pre Nodes that are fixed as constraints are not deleted at a later time
   *
   * Disclaimer: I use the established equality test for Points when testing
   * whether or not a Node's position is equal to a point we would like to fix.
   * The equality operator for Points relies on calling operator== on doubles,
   * which we know to be a terrible idea. Potential workarounds: using
   * nearest node to Point (implemented in previous homework) or changing the
   * equality test for Point objects, so that they consider some small ball of
   * radius epsilon.
   *
   * Nodes with constraints applied to them are moved toward the front of the
   * list of nodes. This is done in case PointsConstraint is called
   * simultaneously with some other constraint that utilizes node/edge removal
   * functions. Since we assume that fixed Nodes are not deleted, and since
   * we know that remove_node(n) may call a swap on node n and node LAST,
   * placing fixed nodes at the very front ensures that the invariant can be
   * satisifed, (See listing of member variables, below).
   *
   * @tparam    GRAPH   Graph must provide iterator to nodes, nodes must
   * provide position() method and value() method. Value must contained field
   * vel that can be set to value Point.
   *
   * Runs in O(fix.size() * g.size()), but only runs once when functor is
   * initialized. For small number of fix.size(), this should not be too bad.
   * This allows for reduced complexity of operator(), which is called in every
   * time step.
   */
  template <typename GRAPH>
  PointsConstraint(GRAPH& g, std::vector<Point> fix) {
    size_type count {0};
    for (auto pt = fix.begin(); pt != fix.end(); ++pt) {
      for (auto nt = g.node_begin(); nt != g.node_end(); ++nt) {
        if ((*pt) == (*nt).position()) {
          // Move toward the front, then add to our internal lists
          g.swap_nodes(g.node(count),(*nt));
          fixedNodes_.push_back(g.node(count));
          fixedPoints_.push_back(*pt);
          ++count;
        }
      }
    }
  }

  private:
  // Storage container for fixed nodes and fixed positions
  // Invariant: These two vectors are of same size and indices correspond to
  // each other at all times. Only established matches go in here.
  std::vector<Point> fixedPoints_;
  std::vector<Node> fixedNodes_;
};

/** Floor constraint for floor at arbitrary position.
 *
 * All points below the floor instantaneously have their z-position changed to
 * zero. Thus, strange effects may be encountered if floor is set too high
 * relative to rest of the graph.
 *
 * Inherits from parent class GeneralizedConstraint. Overrides operator() and
 * stores information about where the floor is located.
 */
class FloorConstraint : public GeneralizedConstraint {
  public:
  /** Constructor for constraint functor.
   *
   * @param[in] floor   The z-level at which the floor is to be set.
   */
  FloorConstraint(double floor) : floor_(floor) {}

  /** Overloaded operator(), applies floor constraint.
   *
   * If position is below floor, z-value of position is changed to z-value of
   * floor, and z-direction of velocity is set to zero.
   *
   * @param[in/out]     g       Graph on which constraints are to be applied
   * @param[in]         t       Time at which constraint is to be applied
   *
   * @tparam    GRAPH   graph must offer node iterators, and nodes most offer
   * position() and value(). Value() must have field vel that is a Point.
   *
   * Method is O(n), since we need to scan through all the nodes.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;           // Quiet compiler warning
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position()[2] < floor_) {
        (*it).position()[2] = floor_;
        (*it).value().vel[2] = 0;
      }
    }
  }

  private:
    /** Store z-position of the floor. */
  double floor_;
};

/** Sphere constraint for sphere with arbitrary location and radius.
 *
 * All points determined to be within sphere are placed along the sphere at the
 * point closest to their current location, and their radial velocity (relative
 * to the sphere center) is set to zero. Thus, strange behavior may be
 * encountered if the sphere is made to intersect with the graph at the initial
 * time.
 *
 * Inherits from parent class GeneralizedConstraint. Overrides operator() and
 * stores information about where the sphere is located.
 */
class SphereConstraint : public GeneralizedConstraint {
  public:
  /** Constructor for constraint functor.
   *
   * @param[in] pos     The position of the sphere center.
   * @param[in] rad     The radius of the sphere
   */
  SphereConstraint(Point pos, double rad) : pos_(pos), rad_(rad) {}

  /** Overloaded operator(), applies sphere constraint.
   *
   * If position is within sphere, it is changed to the location closest to the
   * point on the sphere. Velocity in the radial direction is set to zero.
   *
   * @param[in/out]     g       Graph on which constraints will be applied
   * @param[in]         t       Time at which constraint is applied
   *
   * @tparam    GRAPH   Graph must offer node iterators. Nodes must offer
   * position() and value(). Value must have field vel that is of type Point.
   *
   * Method is O(n), since we need to scan through all the nodes.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;           // Quiet compiler warning
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      Point diff = (*it).position() - pos_;
      if (norm_2(diff) < rad_) {
        // Determine new position
        Point& nodePos = (*it).position();
        nodePos = pos_ + rad_ * (diff/norm_2(diff));

        // Determine new velocity
        Point& nodeVel = (*it).value().vel;
        nodeVel -= (dot(nodeVel, diff/norm_2(diff))) * (diff/norm_2(diff));
      }
    }
  }

  private:
    /** Store sphere position and radius */
  Point pos_;
  double rad_;
};


/** Sphere constraint for sphere with arbitrary location and radius. In this
 * case, all nodes that come into contact with the sphere are removed from the
 * graph.
 *
 * Inherits from parent class GeneralizedConstraint. Overrides operator() and
 * stores information about where the sphere is located.
 */
class DeleteSphere : public GeneralizedConstraint {
  public:
  /** Constructor for constraint functor.
   *
   * @param[in] pos     The position of the sphere center.
   * @param[in] rad     The radius of the sphere
   */
  DeleteSphere(Point pos, double rad) : pos_(pos), rad_(rad) {}

  /** Overloaded operator(), applies sphere constraint.
   *
   * If position is within sphere, it is deleted from the graph.
   *
   * @param[in/out]     g       Graph on which constraints will be applied
   * @param[in]         t       Time at which constraint is applied
   *
   * @tparam    GRAPH   Graph must offer node iterators and have remove_node()
   * capability. Nodes must offer position() that returns Point.
   * 
   * Method is O(r * d), where r is the number of nodes that must be removed
   * in each iteration, and d is the degree of nodes.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;           // Quiet compiler warning
    for (auto it = g.node_begin(); it != g.node_end();) {
      Point diff = (*it).position() - pos_;
      if (norm_2(diff) < rad_) {
        g.remove_node(it);              // O(d)
      } else {++it;}
    }
  }

  private:
    /** Store sphere position and radius */
  Point pos_;
  double rad_;
};

/** Constraint object for summing GeneralizedConstraint children.
 *
 * Inherits from class GeneralizedConstraint. Overrides operator(), and keeps
 * track of individual Constraint functors as well. Currently written to
 * combine three constraints, where 1, 2, or 3 of the inputs could be instances
 * of GeneralizedConstraint(), which is the trivial constraint.
 *
 * Templated on the three types that are passed into the constructor
 * @tparam      C1      First type of constraint to be summed
 * @tparam      C2      Second type of constraint to be summmed
 * @tparam      C3      Third type of constraint to be summed
 */
template <typename C1, typename C2, typename C3>
class CombinedConstraint : public GeneralizedConstraint {
  public:
  /** Overload of operator(), applies combined constraints to graph.
   *
   * @param[in/out]     g       Graph for which constraints are applied.
   * @param[in]         t       Time at which constraint is to be calcualated.
   *
   * @tparam    GRAPH   Templated to take in a generic Graph, as long as the
   *                    graph meets all the requirements of consitutent
   *                    constraint objects.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g, t); c2_(g, t); c3_(g,t);
    return;
  };

  /** Constructor for two or three constraint arguments.
   *
   * @param[in] con1    First constraint object to be summed.
   * @param[in] con2    Second constraint object to be summed.
   * @param[in] con3    Third constraint object to be summed. Defaults to
   *                    GeneralizedForce(), which effectively does not change
   *                    the sum of the previous two inputs.
   */
  CombinedConstraint(C1 con1, C2 con2, C3 con3 = GeneralizedConstraint()) :
    c1_(con1),  c2_(con2), c3_(con3) {}

  private:
  /** Storage for individual constraint components. */
  C1 c1_; C2 c2_; C3 c3_;
};

/** Helper function to predict types when constructing a combined constraint.
 *
 * @tparam C1
 * @tparam C2
 * @tparam C3
 * Currently written for up to three constraints. The third argument is
 * optional and defaults to GeneralizedConstraint(), which effectively does
 * nothing since it is the trivial constraint.
 */
template <typename C1, typename C2, typename C3 = GeneralizedConstraint>
CombinedConstraint<C1, C2, C3> make_combined_constraint(C1 con1,
                                  C2 con2, C3 con3 = GeneralizedConstraint()){
  return CombinedConstraint<C1, C2, C3>(con1, con2, con3);
}

//
//FORCES
//

/** Parent class for force objects.
 *
 * Enforces class design for users who wish to create their own force objects,
 * which can be written as children who inherit from this class. Relies on type
 * Node, as specified at the top of this file. Default virtual method provides
 * zero force, which is useful in the applications below. Additionally,
 * requiring users to inherit publicly from this class ensures that operator()
 * will be implemented to return a valid Force (Point(0,0,0) in this case).
 */
class GeneralizedForce {
  public:
  /** Virtual method that is to be overridden in child classes
   *
   * @param[in]         n       Node for which force is to be calculated.
   * @param[in]         t       Time at which force is to be calcualted.
   * @param[out]        Point   Force vector applied to node @a n
   */
  virtual Point operator()(Node n, double t) {
    (void) t; (void) n;         // Quiet compiler warnings
    return Point {0,0,0};
  }
};


/** Force function object for gravity forces.
 *
 * Inherits from parent class GeneralizedForce. Overrides operator().
 */
class GravityForce : public GeneralizedForce {
  public:
  /** Overload of operator(), returns force of gravity for node n.
   *
   * @param[in]         n       Node for which force is to be calculated.
   * @param[in]         t       Time at which force is to be calcualted.
   * @param[out]        Point   Force vector applied to node @a n
   *
   * grav is defined as a constant at the top of this file.
   *
   * @tparam NODE Templated to take in a generic node, as long as that node has
   * member function value() that returns a node_value with member mass of type
   * double.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;           // Quiet compiler warning
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Force function object for spring forces.
 *
 * Inherits from class GeneralizedForce. Overrides operator().
 */
class MassSpringForce : public GeneralizedForce {
  public:
  /** Overload of operator(), returns force of springy edges for node n.
   *
   * @param[in]         n       Node for which force is to be calculated.
   * @param[in]         t       Time at which force is to be calcualted.
   * @param[out]        Point   Force vector applied to node @a n
   *
   * @tparam NODE
   * Templated to take in a generic node, as long as that node offers iterators
   * for adjacent edges. Additionally, edges are required to store double for
   * rest length and double for spring constant inside the edge_value.
   * Furthermore, nodes are required to have member function position() that
   * returns Point.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;           // Quiet compiler warning
    Point force {0,0,0};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto aEdge = *it;
      double K = aEdge.value().sprCoeff;
      double L = aEdge.value().restL;

      auto aNode = aEdge.node2();
      double l2Norm = norm_2(n.position() - aNode.position());
      force += -K*((n.position() - aNode.position())/l2Norm)*(l2Norm - L);
    }

    return force;
  }
};

/** Force function object for damping forces.
 *
 * Inherits from class GeneralizedForce. Overrides operator().
 */
class DampingForce : public GeneralizedForce {
  public:
  /** Overload of operator(), returns damping force for node n.
   *
   * @param[in]         n       Node for which force is to be calculated.
   * @param[in]         t       Time at which force is to be calcualted.
   * @param[out]        Point   Force vector applied to node @a n
   *
   * @tparam NODE
   * Templated to take in a generic node, as long as that node has member
   * function value() that returns a node_value with member dampCoeff of type
   * double.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;           // Quiet compiler warning
    double C = n.value().dampCoeff;
    return (-C) * n.value().vel;
  }
};

/** Force object for summing GeneralizedForce children.
 *
 * Inherits from class GeneralizedForce. Overrides operator(), and keeps track
 * of individual Force functors as well. Currently written to combine three
 * forces, where 1, 2, or 3 of the inputs could be instances of
 * GeneralizedForce(), which adds 0 to the Force Sum. Effectively, this does
 * nothing.
 *
 * Templated on the three types that are passed into the constructor
 * @tparam      F1      First type of Force to be summed
 * @tparam      F2      Second type of Force to be summmed
 * @tparam      F3      third type of Force to be summed
 */
template <typename F1, typename F2, typename F3>
class CombinedForce : public GeneralizedForce {
  public:
  /** Overload of operator(), returns combined force for node n.
   *
   * @param[in]         n       Node for which force is to be calculated.
   * @param[in]         t       Time at which force is to be calcualted.
   * @param[out]        Point   Force vector applied to node @a n
   *
   * @tparam NODE
   * Templated to take in a generic node, as long as that node meets all the
   * requirements of consitutent force objects.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t) + f3_(n,t);
  };

  /** Constructor for two or three force arguments.
   *
   * @param[in] force1  First force object to be summed.
   * @param[in] force2  Second force object to be summed.
   * @param[in] force3  Third force object to be summed. Defaults to
   *                    GeneralizedForce(), which effectively does not change
   *                    the sum of the previous two inputs.
   */
  CombinedForce(F1 force1, F2 force2, F3 force3 = GeneralizedForce()) :
    f1_(force1),  f2_(force2), f3_(force3) {}

  private:
  /** Storage for individual force components. */
  F1 f1_; F2 f2_; F3 f3_;
};

/** Helper function to predict types when constructing a combined force.
 *
 * @tparam      F2
 * @tparam      F2
 * @tparam      F3
 * Currently written for up to three forces. The third argument is optional and
 * defaults to GeneralizedForce(), which effectively does not contribute to the
 * sum of the first two forces.
 */
template <typename F1, typename F2, typename F3 = GeneralizedForce>
CombinedForce<F1, F2, F3> make_combined_force(F1 force1,
                                              F2 force2,
                                              F3 force3 = GeneralizedForce()){
  return CombinedForce<F1, F2, F3>(force1, force2, force3);
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
    // Computation of constants
    (void) t;           // Quiet compiler warning
    double m = n.value().mass;

    // Constrain two points to be fixed in space by returning zero force
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }

    // Calculate force caused by adjacent nodes
    Point force {0,0,0};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto aEdge = *it;
      double K = aEdge.value().sprCoeff;
      double L = aEdge.value().restL;

      auto aNode = aEdge.node2();
      double l2Norm = norm_2(n.position() - aNode.position());
      force += -K*((n.position() - aNode.position())/l2Norm)*(l2Norm - L);
    }

    // Add force of gravity and return
    force += (m * Point(0,0,-grav));
    return force;
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

  // Initialize node mass and damping coefficient values
  // Velocity defaults to 0, no need to set that manually
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    // Change mass for grid3 case, to avoid changing time step.
    (*it).value().mass = graph.num_nodes() < 10000 ?
                                1.0/graph.num_nodes() : 3.0/graph.num_nodes();
    (*it).value().dampCoeff = 1.0/graph.num_nodes();
  }

  // Initialize edge length initial values
  // Spring constant defaults to K = 100, no need to set that manually
  for (auto it = graph.edge_begin(); it !=graph.edge_end(); ++it) {
    auto e = (*it);
    e.value().restL = norm_2(e.node1().position() - e.node2().position());
  }


  // Set up forces ahead of time, for speed.
  auto forces = make_combined_force(GravityForce(),
                                    MassSpringForce(),
                                    DampingForce());

  // Set up constraints ahead of time, for speed.
  // IF YOU'RE NOT DELETING NODES/EDGES, THEN CHANGE SOME VIEWER OPTIONS BELOW
  auto corners = PointsConstraint(graph,
                  std::vector<Point> {Point(0,0,0), Point(1,0,0)});

  auto constraints = make_combined_constraint(corners,
                        FloorConstraint(-0.75),
                        DeleteSphere(Point(0.5,0.5,-0.5),0.15));

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
        symp_euler_step(graph, t, dt, forces, constraints);

        // Clear nodes and edges (COMMENT OUT IF NOT DELETING NODES/EDGES)
        // NOT SUPER CRUCIAL, BUT JUST LOOKS BETTER ON MY MACHINE
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        
        // Update viewer with edges
        // DON'T DO THIS IF NOT DELETING NODES/EDGES
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
