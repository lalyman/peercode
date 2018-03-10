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

//
// TYPE AND CONSTANT DEFINITIONS
//

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Default elastic constant
static constexpr double defaultK = 100.;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double rest_length;  //< Resting length
  double K;            //< Elastic constant
  EdgeData() : rest_length(0.), K(0.) {}
};

// Define the Graph type
using size_type = unsigned; 
using GraphType = Graph<NodeData, EdgeData>;
using NodeType  = typename GraphType::node_type;
using EdgeType  = typename GraphType::edge_type;
using NodeIter  = typename GraphType::node_iterator;
using EdgeIter  = typename GraphType::edge_iterator;
using IncidentIter  = typename GraphType::incident_iterator;

//
// EULER INTEGRATION
//

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraints
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports velocity and mass. G::edge_value supports
 *         resting length and elastic constant.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
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

  // Apply the constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//
// FORCE FUNCTORS
//

/** Null force function object.
  * @brief This force represents the zero force vector. 
  */
struct NullForce {
  /** Return the zero force applying to @a n at time @a t.
   *
   * Used to model points that do not move, or default
   * force applied when we only want a pair of forces instead of three.
   * @return A force == (0,0,0)
   */
  template <typename NODE>  
  Point operator()(NODE, double) {
    return Point (0,0,0);
  }
};

/** Gravity force function object.
  * @brief This force represents the vertical gravity force vector. 
  */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t.
   *
   * @return A force == m * (0,0,-g) */
  Point operator()(NodeType n, double) {
    return n.value().mass * Point (0,0,-grav);
  }
};

/** Elastic spring force function object.
  * @brief This force represents the effects of a system of springs applied
  *        to a single node. 
  */
struct MassSpringForce {
  /** Return the elastic force applied to @a n at time @a t.
   *
   * @return Elastic force applied to @a n by the system of springs. */
  Point operator()(NodeType n, double) {
    Point force{0,0,0};

    // Compute the force due to each edge (spring) incident on the node
    for(IncidentIter it = n.edge_begin(); it != n.edge_end(); ++it) {
      double restLength = (*it).value().rest_length;
      double elasticK   = (*it).value().K;
      force -= elasticK * 
              ((n.position() - (*it).node2().position()) / (*it).length())  * 
              ((*it).length() - restLength);
      } 

    return force;
  }
};

/** Drag damping force function object.
  * @brief This force represents the effects of drag, as a damping force.
  */
struct DampingForce {

  private:
  /** Damping constant c. */
  const double c_;

  public:
  /** Damping force constructor. */
  DampingForce(const double c) : c_(c) {}

  /** Return the damping force applied to @a n at time @a t.
   *
   * @return A force == -c * n.vel */
  Point operator()(NodeType n, double) {
    return -c_ * n.value().vel;
  }
};

/** Combined force function object.
 * @brief This class represents a combination of forces acting on nodes.
 *
 * @tparam FRC1 Type of the first force function to be combined.
 * @tparam FRC2 Type of the second force function to be combined.
 * @tparam FRC3 Type of the third force function to be combined. Defaults to
                the zero force.
 */
template <typename FRC1, typename FRC2, typename FRC3 = NullForce> 
struct CombinedForce {

  private:
  /* Forces to combine */
  FRC1 f1_;
  FRC2 f2_;
  FRC3 f3_;

  public:
  /** Combined force constructor */
  CombinedForce(FRC1 f1, FRC2 f2, FRC3 f3 = FRC3()) : 
      f1_(f1), f2_(f2), f3_(f3) {}

  /** Return the combination of forces applied to @a n at time @a t.
   *
   * @return A force == f1 + f2 +f3 
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t) + f3_(n,t); 
  }
};

/** Create a combined force from two or three force functors.
 * @param[in]     f1  First force to be combined.
 * @param[in]     f2  Second force to be combined
 * @param[in]     f3  Third force to be combined. This is optional.
 * @return the CombinedForce representing the addition of the three individual
 *         forces
 *
 * @tparam FRCX are the types of forces to be combined. See CombinedForce for
 *         more details.
 */
template <typename FRC1, typename FRC2, typename FRC3 = NullForce> 
const CombinedForce<FRC1,FRC2,FRC3> make_combined_force(FRC1 f1, 
                                                        FRC2 f2, 
                                                        FRC3 f3 = FRC3() ) {
  return CombinedForce<FRC1, FRC2, FRC3>(f1,f2,f3);
}

//
// CONSTRAINT FUNCTORS
//

/** Null constraint function object.
  * @brief This constraint represents the null constraint (no restrictions).
  *        It is used as a default constraint for the combined constraints. 
  */
class NullConstraint {

  public:
  template <typename GRAPH>
  void operator()(GRAPH&, double) {
  }
};

/** Fixed point constraint function object.
  * @brief This constraint represents the nodes fixed in one position. 
  */
class FixedPointConstraint {
  private:
  /** Container used to store the nodes and locations to set as constant. */
  std::vector<std::pair<NodeType, Point>> fixnodes_;

  public:
  /** Constructor */
  FixedPointConstraint() {}

  /** Add a node to set to a constant position constraint.
   * @param[in] n Node to set to a constant position.
   * @param[in] p Point where the node should be fixed.
   */
  void addPoint(const NodeType n, const Point p) {
    fixnodes_.push_back(std::make_pair(n,p));
  }

  /** Enforce constant position constraint.
   * @param[in] g Graph for which this constraint needs to be enforced.
   * @param[in] t Current time.
   */
  template <typename GRAPH>
  void operator()(GRAPH&, double) {
    for(auto it = fixnodes_.begin(); it != fixnodes_.end(); ++it) {
      (*it).first.position() = (*it).second;
    }
  }
};

/** Planed constraint function object.
  * @brief This constraint represents a plane aligned with one of the 3D axis
  *        which nodes cannot go through.
  */
class PlaneConstraint {

  private:
    /** Normal of the plane. It should be (1,0,0), (0,1,0) or (0,0,1). */
    Point  position_;
    /** Location of the node along the axis.*/
    double magn_;

  public:
  /** PlaneConstraint constructor */
  PlaneConstraint(Point pos, double mag) : position_(pos), magn_(mag) {}

  /** Enforce plane constraint.
   * @param[in] g Graph for which this constraint needs to be enforced.
   * @param[in] t Current time.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double) {
    for(NodeIter it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      
      // Compute the distance from the node to the plane
      double distance  = dot(n.position(), position_);

      // If the node were to go through
      if(distance < magn_) {
        // Reset the position (on the plane) and zero the velocity normal to 
        // plane.
        n.position()  += position_ * (magn_ - distance);
        n.value().vel *= (Point(1.,1.,1.) - position_);
      }
    }
  }
};

/** Solid sphere constraint function object.
  * @brief This constraint represents a solid sphere which nodes cannot 
  *        go through.
  */
class SolidSphereConstraint {

  private:
    /** Center of the sphere */
    Point  center_;
    /** Radius of the sphere */
    double radius_;

  public:
  /** SolidSphereConstraint constructor */
  SolidSphereConstraint(Point cen, double rad) : center_(cen), radius_(rad) {}

  /** Enforce the solid sphere constraint.
   * @param[in] g Graph for which this constraint needs to be enforced.
   * @param[in] t Current time.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double) {
    for(NodeIter it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      
      // Compute the distance to the center
      double distance  = norm(n.position() - center_);

      // Compute the normalized radial vector
      Point  R = (n.position() - center_) / distance;

      if(distance < radius_) {
        n.position()   = center_ + R * radius_;
        n.value().vel -= dot(n.value().vel, R) * R; 
      }
    }
  }
};

/** Hole sphere constraint function object.
  * @brief This constraint represents a sphere that removes the nodes 
  *        (and its incident edges) that try to go through it.
  */
class HoleSphereConstraint {

  private:
    /** Center of the sphere */
    Point  center_;
    /** Radius of the sphere */
    double radius_;

  public:
  /** HoleSphereConstraint constructor */
  HoleSphereConstraint(Point cen, double rad) : center_(cen), radius_(rad) {}

  /** Enforce the hole sphere constraint.
   * @param[in] g Graph for which this constraint needs to be enforced.
   * @param[in] t Current time.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double) {
    for(NodeIter it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;

      // Compute the distance to the center
      double distance  = norm(n.position() - center_);

      if(distance < radius_) {
        g.remove_node(n);
      }
    }
  }
};

/** Combined constraint function object.
 * @brief This class represents a combination of constraints acting on nodes. 
 * The constraints are applied in the order they are provided in
 * the constructor.
 *
 * @tparam CNSTR1 Type of the first constraint function to be combined.
 * @tparam CNSTR2 Type of the second constraint function to be combined.
 * @tparam CNSTR3 Type of the third constraint function to be combined. 
 *                Defaults to the null constraint.
 */
template <typename CNSTR1, typename CNSTR2, typename CNSTR3 = NullConstraint> 
struct CombinedConstraint {

  private:
  /* Constraints to combine */
  CNSTR1 c1_;
  CNSTR2 c2_;
  CNSTR3 c3_;

  public:
  /** Combined constraint constructor */
  CombinedConstraint(CNSTR1 c1, CNSTR2 c2, CNSTR3 c3 = CNSTR3()) : 
      c1_(c1), c2_(c2), c3_(c3) {}

  /** Apply the combination of constraints to @a g at time @a t.
   *
   * @param[in] g Graph to constraint.
   * @param[in] t Current time.
   * @pre The combined constraints are compatible.
   * @post The Graph is consistent with the constraints provided,
   *       in the order they are executed.
   */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g,t);
    c2_(g,t);
    c3_(g,t); 
  }
};

/** Create a combined constraint from two or three constraint functors.
 * The constraints will be applied in the order provided to this function.
 * @param[in]     f1  First constraint to be combined.
 * @param[in]     f2  Second constraint to be combined
 * @param[in]     f3  Third constraint to be combined. This is optional.
 * @return the CombinedConstraint representing the application of the three 
 * individual constraints.
 *
 * @tparam CNSTRX are the types of constraints to be combined. See
 *         CombinedConstraint for more details.
 */
template <typename CNSTR1, typename CNSTR2, typename CNSTR3 = NullConstraint> 
const CombinedConstraint<CNSTR1,CNSTR2,CNSTR3> 
      make_combined_constraint(CNSTR1 c1, CNSTR2 c2, CNSTR3 c3 = CNSTR3() ) {
  return CombinedConstraint<CNSTR1, CNSTR2, CNSTR3>(c1,c2,c3);
}

//
// NEAREST NODE
//

/** @class MinDistance
 * @brief Functor to compare nodes based on distance to a point.
 *
 * The functor sets a strict ordering based on Node distance
 * to a certain point set on construction.
 */
class MinDistance {
  private:
    /** Point against to which measure euclidean distance */
    const Point point_;
  public:
  /** Construct a MinDistance with a point
   * @param[in] p  Point to measure Node distance
   * @post This MinDistance will consider point @a p as the 
   *       reference to measure Node distance.
   */ 
  MinDistance(const Point p) : point_(p) {}

  /** Compare Node distance to point.
   * This operator compares the euclidan distance of two 
   * Node objects to a point P.
   * @param[in] a  Node to compare.
   * @param[in] b  Node to compare.
   * @pre The Node objects have a position.
   *
   * @return True if @a a is closer to P than @a b, false otherwise.
   */ 
  bool operator() (const NodeType& a, const NodeType& b) {
    // Compute the Nodes' Euclidean distances to the point
    double distanceA = norm(a.position()-point_);
    double distanceB = norm(b.position()-point_);

    return distanceA < distanceB;
  }
}; 

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */
NodeIter nearest_node(const GraphType& g, const Point& point)
{
  // Build our comparator functor
  MinDistance comp(point);

  return std::min_element(g.node_begin(), g.node_end(), comp);
}

//
// MAIN
//

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

  // Set initial conditions for your nodes, if necessary.
  double initial_mass = 1. / graph.size();

  for(NodeIter it = graph.node_begin(); it != graph.node_end(); ++it) {
   (*it).value().vel  = Point{0.,0.,0.};
   (*it).value().mass = initial_mass;
  }

  for(EdgeIter it = graph.edge_begin(); it != graph.edge_end(); ++it) {
   (*it).value().rest_length  = (*it).length();
   (*it).value().K            = defaultK;
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
      double dt = (initial_mass < 1./9999.) ? 0.0005 : 0.001;
      double t_start = 0;
      double t_end = 5.0;
    
      auto forces = make_combined_force(GravityForce(),
                                        MassSpringForce(),
                                        DampingForce(1./graph.size()));

      FixedPointConstraint fp;
      auto anchor1 = nearest_node(graph, Point(0.,0.,0.));
      auto anchor2 = nearest_node(graph, Point(1.,0.,0.));
      fp.addPoint((*anchor1), (*anchor1).position());
      fp.addPoint((*anchor2), (*anchor2).position());

      PlaneConstraint plane(Point(0.,0.,1.), -0.75);
      //SolidSphereConstraint sphere(Point(0.5,0.5,-0.5), 0.15);
      HoleSphereConstraint sphere(Point(0.5,0.5,-0.5), 0.15);

      auto constraints = make_combined_constraint(fp, plane, sphere);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, forces, constraints);

        // Clear  the  viewer ’s nodes  and  edges
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
