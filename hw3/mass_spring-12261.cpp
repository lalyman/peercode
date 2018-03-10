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

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  //HW2_CJ
  double num_nodes; //< num nodes in graph
  NodeData() : vel(0), mass(1) {}
};

//HW2_CJ
/** Custorm structure of data to store with Edges */
struct EdgeData {
  double K_;  //< Spring constant
  double L_;  //< Rest length
  EdgeData() : K_(1), L_(1) {}
  EdgeData(double K, double L) : K_(K), L_(L) {}
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    //HW2_CJ
    /*
    //skip update step for node (0,0,0) and node (1,0,0)
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      //skip update step
    }
    else {
    */
      // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    //}
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    //HW2_CJ
    /*
    //skip update step for node (0,0,0) and node (1,0,0)
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      //skip update step
    }
    else {
    */
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    //}
  }

  return t + dt;
}

//HW2_CJ

//
// FORCES
//

/*** Gravity force function object ***/
struct GravityForce {
  /** Returns the gravity force applying to @a n at time @a t. **/
  public:
    //Calculate the gravity force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      //constrain two corners of the cloth by returning zero force
      //if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      //  return Point(0,0,0);
      //}
      //force due to gravity
      Point f_grav = Point(0,0, n.value().mass*(-grav));
      
      //silence compiler warning 
      (void) t;
      return f_grav; 
    }
};

/*** Spring force function object ***/
struct MassSpringForce {
  /** Returns the spring force applying to @a n at time @a t. **/
  public:
    //Calculate spring force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      //constrain two corners of the cloth by returning zero force
      //if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      //  return Point(0,0,0);
      //}
      //spring force
      Point f_spring = Point(0,0,0);
      //use incident iterator to iterate through adjacent node
      for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
        auto edge = *ei;
        auto node1 = edge.node1();
        auto node2 = edge.node2();
        //make n equal to node1
        if (node2 == n) {
          node2 = edge.node1();
          node1 = edge.node2();
        }
        auto n1_pos = node1.position();
        auto n2_pos = node2.position();

        f_spring += (-edge.value().K_)*(n1_pos - n2_pos)*(edge.length() - edge.value().L_)/edge.length();
      }
      //silence warning for time
      (void) t;
      return f_spring; 
    }

};

/*** Damping force function object ***/
struct DampingForce {
  /** Returns the damping force applying to @a n at time @a t. **/
  public:
    //Calculate spring force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      Point f_damping = Point(0,0,0);
      double c = 1.0/double(n.value().num_nodes);
      f_damping = -c*(n.value().vel);
      (void) t;
      return f_damping;    
    }
    
};

//HW2_CJ

/** @brief Functor to return the sum of two forces
*   @param[in] f1  templated force functor
*   @param[in] f2  templated force functor
*   @pre NODE n is within graph
*   return sum of two forces
*   
* @tparam F1, F2 are function objects called as @a force(n, @a t),
*                where n is a node of the graph and @a t is the current time.
*                @a force must return a Point representing the force vector on
*                Node n at time @a t.
*/
template <typename F1, typename F2>
struct combine_2forces {
  public:
    //Constructor - two forces
    combine_2forces(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

    //Calculate total force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      return f1_(n,t) + f2_(n,t);
    }

  private:
    F1 f1_;
    F2 f2_;
};


/** @brief Functor to return the sum of three forces
*   @param[in] f1  templated force functor
*   @param[in] f2  templated force functor
*   @param[in] f3  templated force functor
*   @pre NODE n is within graph
*   return sum of three forces
*   
* @tparam F1,F2,F3 are function objects called as @a force(n, @a t),
*                where n is a node of the graph and @a t is the current time.
*                @a force must return a Point representing the force vector on
*                Node n at time @a t.
*/
template <typename F1, typename F2, typename F3>
struct combine_3forces {
  public:
    //Constructor - three forces
    combine_3forces(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {}

    //Calculate total force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      return f1_(n,t) + f2_(n,t) + f3_(n,t);
    }

  private:
    F1 f1_;
    F2 f2_;
    F3 f3_; 
};

/*** Function to call combine_2forces() struct and return sum of 2 forces ***/
template <typename F1, typename F2>
combine_2forces<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return {f1, f2};
}

/*** Function to call combine_3forces() struct and return sum of 3 forces ***/
template <typename F1, typename F2, typename F3>
combine_3forces<F1,F2,F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return {f1, f2, f3};
}

//
// CONSTRAINTS
//

/*** Constraint for desired fixed nodes ***/
struct ConstNodeConstraint {
  public:
    template <typename G>
    void operator()(G& g, double t) {
      //search for nodes in graph that violate the constraint
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        //constant node constraints
        if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
          //reset velocity so these points do not move
          n.value().vel = Point(0,0,0);
        }
      }
      (void) t; //silence compiler warning
    }

};

/*** Constraint for desired z plane ***/
struct PlaneConstraint {
  public:
    template <typename G>
    void operator()(G& g, double t) {
      double z_bound = -0.75;
      //search for nodes in graph that violate the constraint
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        //plane constraints
        if (n.position().z < z_bound) {
          //set the position to the nearest point on the plane
          n.position() = Point(n.position().x, n.position().y, z_bound);
          //set z-component of Node velocity to zero
          n.value().vel.z = 0.0;
        }
      }
      (void) t; //silence compiler warning
    }
};

/*** Constraint for a sphere ***/
struct SphereConstraint {
  public:
    template<typename G>
    void operator()(G& g, double t) {
      //sphere properties
      double radius = 0.15;
      Point center = Point(0.5,0.5,-0.5);
      //search for nodes in graph that violate the constraint
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        Point R = (n.position() - center)/(norm_2(n.position() - center));
        Point diff_center = n.position() - center;
        auto norm_diff_center = norm_2(diff_center);
        //sphere constraint
        if (norm_diff_center < radius) {
          //set the position to the nearest point on the sphere surface
          n.position() = center + (diff_center)*radius/(norm_diff_center);
          //set component of velocity normal to sphere's surface to zero
          n.value().vel = n.value().vel - ((n.value().vel*R)*R);      
        } 
      (void) t; //silence compiler warning
      }
    }
};

/*** Constraint to delete a sphere ***/
struct SphereDelete {
  public:
    template<typename G>
    void operator()(G& g, double t) {
      //sphere properties
      double radius = 0.15;
      Point center = Point(0.5,0.5,-0.5);
      //search for nodes in graph that violate the constraint
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        Point diff_center = n.position() - center;
        auto norm_diff_center = norm_2(diff_center);
        //sphere constraint
        if (norm_diff_center < radius) {
          //remove node and all of its edges
          g.remove_node(n);
        } 
      (void) t; //silence compiler warning
      }
    }
};

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  //HW2_CJ
  public:
    //Constructor 
    //Problem1Force(double K, double L) : K_(K), L_(L) {}

    //Calculate total force
    template <typename NODE>
    Point operator()(NODE n, double t) {
      // HW2 #1: YOUR CODE HERE
      //constrain two corners of the cloth by returning zero force
      if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
      }
      //force due to gravity
      Point f_grav = Point(0,0, n.value().mass*(-grav));

      //spring force
      Point f_spring = Point(0,0,0);
      Point f_total = Point(0,0,0);
      //use incident iterator to iterate through adjacent node
      for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
        auto edge = *ei;
        auto node1 = edge.node1();
        auto node2 = edge.node2();
        //make n equal to node1
        if (node2 == n) {
          node2 = edge.node1();
          node1 = edge.node2();
        }
        auto n1_pos = node1.position();
        auto n2_pos = node2.position();

        f_spring += (-edge.value().K_)*(n1_pos - n2_pos)*(edge.length() - edge.value().L_)/edge.length();
      }

      f_total = f_grav + f_spring;

      //silence warning for time
      (void) t;

      return f_total;
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

  //HW2_CJ
  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double K = 100;  //spring constant
  unsigned int N = graph.num_nodes();  //number of nodes in graph
  double m = 1.0/double(N);  //mass (constant density)

  //initialize Node data members  
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value() = NodeData();
    n.value().mass = m; 
    n.value().num_nodes = N;
  }
 
  //initialize Edge data members
  for (auto et = graph.edge_begin(); et != graph.edge_end(); ++et) {
    auto e = *et;
    double L = e.length();
    //e.value() = EdgeData(K, L);
    e.set_edge_value(EdgeData(K, L));
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
      //double dt = 0.001;
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      //HW2_CJ
      /*** Declare constraint functors ***/
      ConstNodeConstraint ConstNode = ConstNodeConstraint();
      PlaneConstraint PlaneConst = PlaneConstraint();
      SphereConstraint SphereConst = SphereConstraint();
      SphereDelete SphereDel = SphereDelete();
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        /*Problem 3*/
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce()));
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce())); 
        /*Problem 2*/
        //symp_euler_step(graph, t, dt, Problem1Force());
        /*Problem 1*/
        //symp_euler_step(graph, t, dt, Problem1Force(K,L));

        /*** Update Nodes and Edges by re-drawing graph ***/
        //clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();
        //update the viewer with nodes' new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

 
        /*** APPLY CONSTRAINTS TO NODES ***/
        //iterate through nodes to see if constraints are satisfied 
        ConstNode(graph,t);
        PlaneConst(graph,t);
        //SphereConst(graph,t);
        SphereDel(graph,t);

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
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
