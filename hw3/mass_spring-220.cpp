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
static constexpr double spring_const = 100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double len;
  double K;
  EdgeData() : len(1), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using NodeIterator = typename GraphType::node_iterator;

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
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

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

  /**
  * @brief Computes the force
  *
  * @param[in] n The node compute the force for
  * @return The Point force vector
  *
  * @tparam NODE The type GraphType::Node
  *
  * @pre n is a valid NODE object
  * @post result is the force vector computed from the node
  **/
  template <typename NODE>
  Point operator()(NODE n, double) {
    // HW2 #1: YOUR CODE HERE
    Point x_0 = n.position();

    if (x_0 ==  Point (0,0,0) || x_0 ==  Point (1,0,0))
      return Point(0,0,0);

    //Compute the spring force
    Point spring_force = Point(0,0,0);
    for(auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      auto incidentEdge = (*ii);
      Point numerator = x_0 - incidentEdge.node2().position();
      double mag = incidentEdge.length();
      spring_force = spring_force + (-incidentEdge.value().K  *
                            numerator * (mag - incidentEdge.value().len) / mag);
    }
    //Compuute the gravity force
    Point gravity = n.value().mass * Point(0, 0, -grav);
    Point total_force = spring_force + gravity;

    //return total force
    return total_force;
  }
};

/**
 * @brief Compute the Gravity force
 * @note This is a functor
 **/
struct GravityForce {

  /**
  * @brief Computes the gravitational force
  *
  * @param[in] n The node compute the gravitational force for
  * @return The Point gravitational force vector
  *
  * @tparam NODE The type GraphType::Node
  *
  * @pre n is a valid NODE object
  * @post result is the gravitational force vector computed from the node
  **/
  template <typename NODE>
  Point operator()(NODE n, double) {
    //Same code as in Problem1Force
    Point gravity = n.value().mass * Point(0, 0, -grav);
    return gravity;
  }
};

/**
 * @brief Compute the MassSpring force
 * @note This is a functor
 **/
struct MassSpringForce {

  /**
  * @brief Computes the mass spring force
  *
  * @param[in] n The node compute the mass spring force for
  * @return The Point mass spring force vector
  *
  * @tparam NODE The type GraphType::Node
  *
  * @pre n is a valid NODE object
  * @post result is the mass spring force vector computed from the node
  **/
  template <typename NODE>
  Point operator()(NODE n, double) {
    //Same code as in Problem1Force
    Point spring_force = Point(0,0,0);

    for(auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      auto incidentEdge = (*ii);
      Point numerator = n.position() - incidentEdge.node2().position();
      double mag = incidentEdge.length();
      spring_force = spring_force + (-incidentEdge.value().K  * numerator *
                                        (mag - incidentEdge.value().len) / mag);
    }
    return spring_force;
  }
};

/**
 * @brief Compute the Damping force
 * @note This is a functor
 **/
struct DampingForce {
  //Private since it is global across all the edges and nodes
  private:
    double damp;

  public:
    /**
    * @brief Constructor for the functor
    *
    * @param[in] d damping constant
    * @return A DampingForce object
    * @pre none
    * @post result is DampingForce object with the damp = d
    **/
    DampingForce(double d) : damp(d) {}

    /**
    * @brief Computes the damping force
    *
    * @param[in] n The node compute the damping force for
    * @return The Point damping force vector
    *
    * @tparam NODE The type GraphType::Node
    *
    * @pre n is a valid NODE object
    * @post result is the damping force vector computed from the node
    **/
    template <typename NODE>
    Point operator()(NODE n, double) {
      return -damp*n.value().vel;
    }
};

/**
 * @brief Dummy force for the combine function
 * @note This is a functor
 **/
struct ZeroForce {
  /**
  * @brief Computes the zero force
  *
  * @return The zero force
  *
  * @tparam NODE The type GraphType::Node
  *
  * @post result is the zero force
  **/
  template <typename NODE>
  Point operator()(NODE, double) {
    return Point(0, 0, 0);
  }
};

/**
 * @brief Compute the Combined force. Default third argument as zero force
 * @note This is a functor
 **/
template <typename T1, typename T2, typename T3 = ZeroForce>
struct CombinedForce {

  //Store the forces
  private:
    T1 f1;
    T2 f2;
    T3 f3;

  public:
    /**
    * @brief Constructor for the functor
    *
    * @param[in] ff1 functors
    * @param[in] ff2 functors
    * @param[in] ff3 functors
    *
    * @return A CombinedForce object
    * @pre none
    * @post result is CombinedForce object that is initialized with params
    **/
    CombinedForce(T1 ff1, T2 ff2, T3 ff3 = ZeroForce()) : f1(ff1), f2(ff2),
                                                          f3(ff3) {}

    /**
    * @brief Computes the combined force
    *
    * @param[in] n The node compute the damping force for
    * @param[in] t time step
    * @return The combined force
    *
    * @tparam NODE The type GraphType::Node
    *
    * @pre n is a valid NODE object
    * @post result is the combined force
    **/
    template <typename NODE>
    Point operator()(NODE n, double t) {
      return f1(n, t) + f2(n, t) + f3(n, t);
    }
};

/**
* @brief Computes the combined force
*
* @param[in] f1 Functor force
* @param[in] f2 Functor force
* @param[in] f3 Functor force

* @return The CombinedForce object
*
* @tparam T1 GravityForce, MassSpringForce, ZeroForce, or DampingForce
* @tparam T2 GravityForce, MassSpringForce, ZeroForce, or DampingForce
* @tparam T3 GravityForce, MassSpringForce, ZeroForce, or DampingForce
*
* @pre f1, f2, f3 are all valid functor objects.
* @post result is the CombinedForce object with the functors initialized
**/
template <typename T1, typename T2, typename T3 = ZeroForce>
CombinedForce<T1, T2, T3> make_combined_force(T1 f1, T2 f2,
                                              T3 f3 = ZeroForce()) {
  return CombinedForce<T1, T2, T3> (f1, f2, f3);
}

/**
 * @brief Compute the plane constraint
 * @note This is a functor
 **/
struct PlaneConstraint {

  /**
  * @brief Computes the plane constraint z > -0.75
  *
  * @param[in] g The graph to constrain
  *
  * @tparam G The type GraphType
  *
  * @pre g is a valid GraphType object
  * @post Nodes that are constrained should have the appropriate velocities
  *       and positions modified
  **/
  template <typename G>
  void operator()(G& g, double) {
    //Iterate over all the nodes. Those that cross the plane z = -0.75
    //are then mapped back to the plane and have the normal component zeroed
    //out.
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto node = *ni;
      if(dot(node.position(),Point(0,0,1)) < -0.75) {
        node.position() = node.position()*Point(1,1,0) + Point(0,0,-0.75);
        node.value().vel = node.value().vel*Point(1,1,0);
      }
    }
  }
};

/**
 * @brief Compute the point constraint
 * @note This is a functor
 **/
struct PointConstraint {

  //Store the points to constrain
  private:
    std::vector<Node> nodesConstrain;
    std::vector<Point> pointConstrain;

  public:
    /**
    * @brief Constructor for the functor
    *
    * @return A PointConstraint object
    *
    * @pre none
    * @post result is PointConstraint object
    **/
    PointConstraint() {}

    /**
    * @brief Adds the points that we want to constrain
    *
    * @param[in] n Node to constrain
    * @param[in] p point to constrain at
    *
    * @pre none
    * @post result is points and nodes are added to the respective vector
    **/
    void addPointsToConstrain(Node n, Point p) {
      nodesConstrain.push_back(n);
      pointConstrain.push_back(p);
    }

    /**
    * @brief Update the position so that the points are constrained.
    *
    * @param[in] g The graph to constrain
    *
    * @tparam G The type GraphType
    *
    * @pre g is a valid GraphType object
    * @post Nodes that are constrained should have the appropriate velocities
    *       and positions modified
    **/
    template <typename G>
    void operator()(G&, double) {
      for(unsigned i = 0; i < nodesConstrain.size(); i++) {
        nodesConstrain[i].position() = pointConstrain[i];
    }
  }
};

/**
 * @brief Compute the sphere constraint
 * @note This is a functor
 **/
struct SphereConstraint {

  /**
  * @brief Computes the sphere constraint
  *
  * @param[in] g The graph to constrain
  *
  * @tparam G The type GraphType
  *
  * @pre g is a valid GraphType object
  * @post Nodes that are constrained should have the appropriate velocities
  *       and positions modified
  **/
  template <typename G>
  void operator()(G& g, double) {
    //Center and radius of the sphere
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;

    //Check all the nodes
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto node = *ni;
      double mag = norm(node.position() - center);
      //If the node is within the ball
      if(mag < radius) {
        //compute the distance vector
        Point difference = node.position() - center;
        //Compute the position on the ball and update it
        node.position() = center + radius*difference/mag;
        //Compute the normal vector to the ball
        Point R = difference/norm(difference);
        //Zero out the velocity with respect to the normal vector.
        node.value().vel = node.value().vel - dot(node.value().vel, R) * R;
      }
    }
  }
};

/**
 * @brief Compute the sphere constraint and destroy the nodes if they touch
 *        the sphere.
 * @note This is a functor
 **/
struct SphereDestroyConstraint {

  /**
  * @brief Computes the sphere constraint and ddestroy nodes that touch it
  *
  * @param[in] g The graph to modify
  *
  * @tparam G The type GraphType
  *
  * @pre g is a valid GraphType object
  * @post Nodes that touch the sphere should be removed from the graph.
  **/
  template <typename G>
  void operator()(G& g, double) {
    //same as above except remove the nodes now.
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;

    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto node = *ni;
      double mag = norm(node.position() - center);
      if(mag < radius) {
        g.remove_node(node);
      }
    }
  }
};

/**
 * @brief Dummy constraint like with ZeroForce
 * @note This is a functor
 **/
struct ZeroConstraint {

  /**
  * @brief Computes the zero constraint
  *
  * @post Does nothing
  **/
  template <typename G>
  void operator()(G&, double) {
    return;
  }
};

/**
 * @brief Compute the Combined copnstraint. Default third argument is zero
 *        constraint
 * @note This is a functor
 **/
template <typename C1, typename C2, typename C3 = ZeroConstraint>
struct CombinedConstraint {

  //Store the constraints
  private:
    C1 c1;
    C2 c2;
    C3 c3;

  public:
    /**
    * @brief Constructor for the functor
    *
    * @param[in] cc1 functors
    * @param[in] cc2 functors
    * @param[in] cc3 functors
    *
    * @return A CombinedConstraint object
    * @pre none
    * @post result is CombinedConstraint object that is initialized with params
    **/
    CombinedConstraint(C1 cc1, C2 cc2, C3 cc3 = ZeroConstraint()) : c1(cc1),
                                                c2(cc2), c3(cc3) {}

    /**
    * @brief Computes the combined constraint
    *
    * @param[in] g The graph we are constraining
    * @param[in] t time step
    * @return The combined constraint
    *
    * @tparam G The type GraphType
    *
    * @pre g is a valid GraphType object
    * @post result is the combined constraints acting on the graph
    **/
    template <typename G>
    void operator()(G& g, double t) {
      //Activate all the constraints
      c1(g, t);
      c2(g, t);
      c3(g, t);
      return;
    }
};

/**
* @brief Computes the combined constraint
*
* @param[in] c1 Functor constraint
* @param[in] c2 Functor constraint
* @param[in] c3 Functor constraint

* @return The CombinedForce object
*
* @tparam C1 PlaneConstraint, ZeroConstraint, PointConstraint,
              SphereConstraint or SphereDestroyConstraint
* @tparam C2 PlaneConstraint, ZeroConstraint, PointConstraint,
              SphereConstraint or SphereDestroyConstraint
* @tparam C3 PlaneConstraint, ZeroConstraint, PointConstraint,
              SphereConstraint or SphereDestroyConstraint
*
* @pre c1, c2, c3 are all valid functor objects.
* @post result is the CombinedConstraint object with the functors initialized
**/
template <typename C1, typename C2, typename C3 = ZeroConstraint>
CombinedConstraint<C1, C2, C3> make_combined_constraint(C1 c1, C2 c2,
                                                    C3 c3 = ZeroConstraint()) {
  return CombinedConstraint<C1, C2, C3> (c1, c2, c3);
}

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
NodeIterator nearest_node(const GraphType& g, const Point& point) {
  // HW1 #3: YOUR CODE HERE
  NodeIterator min_ni = g.node_begin();
  double min = norm(point - (*g.node_begin()).position());
  for(auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
    auto node = *ni;
    double distance = norm(point - node.position());
    if(distance < min) {
      min = distance;
      min_ni = ni;
    }
  }

  return min_ni;
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
    //
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  //Compute all the node masses and velocities and store them
  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    auto node = *ni;
    node.value().mass = 1/(float)graph.num_nodes();
    node.value().vel = Point(0,0,0);
  }

  //Compute all the edge lengths and the spring constants and store them.
  for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    auto edge = *ei;
    edge.value().len = edge.length();
    edge.value().K = spring_const;
  }

  //Compute the global damping constant
  double damping_const = 1/(float)graph.num_nodes();

  //Initialize the point constraint
  auto pc = PointConstraint();
  auto ni1 = nearest_node(graph, Point(0,0,0));
  auto ni2 = nearest_node(graph, Point(1,0,0));
  pc.addPointsToConstrain(*ni1, Point(0,0,0));
  pc.addPointsToConstrain(*ni2, Point(1,0,0));

  //Create a combined force and combined constraint functors.
  auto F = make_combined_force(GravityForce(), MassSpringForce(),
                               DampingForce(damping_const));
  auto C = make_combined_constraint(pc, SphereDestroyConstraint(),
                                    PlaneConstraint());

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

        //Perform the operations
        symp_euler_step(graph, t, dt, F, C);

        // Update viewer with nodes' new positions
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
