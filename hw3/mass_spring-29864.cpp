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
  double cval;     //< Damping factor
  NodeData() : vel(0), mass(1), cval(0) {}
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K; // Edge spring constant
  double L_rest; // Edge rest length
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
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraint
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Manually set points (0,0,0) and (1,0,0) stationary
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        n.value().vel = Point(0,0,0);
    } else {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
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
public:
  template <typename NODE>
  Point operator()(NODE n, double) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);

    // gravitational force
    Point f_grv = Point(0,0,-grav*n.value().mass); 

    // spring force
    Point f_spr = Point(0,0,0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){

        auto current_edge = *i;
        auto n2 = current_edge.node2();

        auto diff = n.position() - n2.position();
        auto diff_norm = norm(diff);
        f_spr -= current_edge.value().K/diff_norm
                *(diff_norm-current_edge.value().L_rest)*diff;
    } 
    return f_grv+f_spr;
  }
};
        
/** Return gravitational force applied to @a n 
*/
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double){
    return Point(0,0,-grav*n.value().mass); 
  }    
};

/** Return damping force applied to @a n 
*/
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double){
    return (-n.value().cval*n.value().vel); 
  }    
};

/** Return spring force applied to @a n 
*/
struct MassSpringForce {
public:
  template <typename NODE>
  Point operator()(NODE n, double) {
    Point f_spr = Point(0,0,0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){

        auto current_edge = *i;
        auto n2 = current_edge.node2();

        auto diff = n.position() - n2.position();
        auto diff_norm = norm(diff);
        f_spr -= current_edge.value().K/diff_norm
                *(diff_norm-current_edge.value().L_rest)*diff;
    } 
    return f_spr;
  }
};


/** Combine two functors */
template <typename F1, typename F2>
struct add_two_functor{
    private:
        F1 f1_;
        F2 f2_;
    public:
        add_two_functor(F1 f1, F2 f2):f1_(f1),f2_(f2){}
        template<typename NODE>
        Point operator()(NODE& n, double t){
            return f1_(n,t) + f2_(n,t);
        }
};

template <typename F1, typename F2>
add_two_functor<F1, F2> make_combined_force(F1 f1, F2 f2){
    return add_two_functor<F1,F2>(f1,f2);
}

/** Combine three functors */
template <typename F1, typename F2, typename F3>
struct add_three_functor{
    private:
        F1 f1_;
        F2 f2_;
        F3 f3_;
    public:
        add_three_functor(F1 f1, F2 f2, F3 f3):f1_(f1),f2_(f2),f3_(f3){}
        template<typename NODE>
        Point operator()(NODE& n, double t){
            return f1_(n,t) + f2_(n,t) + f3_(n,t);
        }

};

template <typename F1, typename F2, typename F3>
add_three_functor<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3){
    return add_three_functor<F1,F2,F3>(f1,f2,f3);
}



/** HW2 Part 4
 *  The plane constraint enforces all nodes and edges in @a graph
 *  to stay above z = -0.75 plane.
 */
struct c_plane{
    Point operator()(GraphType& graph, double){
        for (auto i = graph.node_begin(); i != graph.node_end(); ++i){
          if ((*i).position().z < -0.75){
              (*i).position().z = -0.75;
              (*i).value().vel.z = 0.0;
          }
        }
    return Point(0,0,0);
    }
};

/** HW2 Part 4
 *  The sphere constraint enforces all nodes and edges in @a graph
 *  to stay on and around a sphere centered at (0.5,0.5,-0.5)
 *  with a radius of 0.15.
 */
struct c_sphere{
  Point operator()(GraphType& graph, double){
    double r = 0.15;
    Point c = Point(0.5,0.5,-0.5);
    for (auto i = graph.node_begin(); i != graph.node_end(); ++i){
        Point R = (*i).position() - c;
        if (norm(R) < r){
            R = R/norm(R);
            (*i).position() = c + r/norm((*i).position()-c)*((*i).position()-c);
            (*i).value().vel -= dot((*i).value().vel,R)*R;
        }
    }
    return Point(0,0,0);
  }
};

/** HW2 Part 5
 *  The plane constraint removes all nodes and edges in @a graph
 *  that came into contact with a sphere centered at (0.5,0.5,-0.5)
 *  with a radius of 0.15.
 */
struct c_circle{
  Point operator()(GraphType& graph, double){
    double r = 0.15;
    Point c = Point(0.5,0.5,-0.5);
    for (auto i = graph.node_begin(); i != graph.node_end(); ++i){
        Point R = (*i).position() - c;
        if (norm(R) <= r){
            graph.remove_node(*i);
        }
    }
    return Point(0,0,0);
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

  // Set initial conditions for nodes
  double m_i = 1.0/graph.num_nodes();
  for (auto j = graph.node_begin(); j != graph.node_end(); ++j){
    (*j).value().cval = m_i;
    (*j).value().mass = m_i;
    (*j).value().vel = Point(0,0,0);
    for (auto i = (*j).edge_begin(); i != (*j).edge_end(); ++i){
      (*i).value().K = 100.0;
      (*i).value().L_rest = (*i).length();
    }
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
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, 
            make_combined_force(GravityForce(), MassSpringForce(),DampingForce()),
            make_combined_force(c_circle(),c_plane()));

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        //if (graph.size() < 100)
        //  std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
