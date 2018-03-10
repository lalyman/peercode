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
#include <ctime>

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
  NodeData(Point v, double m) : vel(v), mass(m) {}
};


struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100), L(1.39435) {}
  EdgeData(double k, double l) : K(k), L(l) {}
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
 * 'this' argument has type 'const Point', but method is not marked const          where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    while (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      ++it;
      n = *it;
    }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    constraints(n);
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    while (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      ++it;
      n = *it;
    }
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double) {
     return Point(0,0,-grav*n.value().mass);    
  }
};

struct MassSpringForce {

  template <typename NODE>
  Point operator()(NODE n, double) {
    Point f {0,0,0};
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      NODE n2 = (*it).node2();
      if(n2 == n){
        n2 = (*it).node1();
      }

      Point diff = n.position() - n2.position();
      double K = (*it).value().K;
      double L = (*it).value().L;
      f = f - K*((diff)/norm(diff))*(norm(diff) - L);
    } 
    return f;
  }
};

struct DampingForce {

  DampingForce() = default;
  DampingForce(double c) : c_(c) {}

  template <typename NODE>
  Point operator()(NODE n, double) {
    return -c_*n.value().vel;  
  }

  private:
    double c_ = 0;
};


template <typename F1, typename F2>
struct CombinedForce {

  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  template<typename NODE>
  Point operator()(NODE n, double t) {
      return f1_(n,t) + f2_(n,t); 
  }

  private:
    F1 f1_;
    F2 f2_;
};

template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1,F2>(f1,f2);
}

template <typename F1, typename F2, typename F3>
CombinedForce<F1,CombinedForce<F2,F3>> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombinedForce<F1,CombinedForce<F2,F3>>(f1,make_combined_force(f2,f3));
}

struct PlaneConstraint {
   template <typename NODE>
   bool violate(NODE n) {
     return n.position().z < -.75;
   }

   template <typename NODE>
   void reset(NODE n) {
     n.position().z = -.75;
     n.value().vel.z = 0; 
   }

   template <typename NODE>
   void operator()(NODE n) {
     if(violate(n)){
       reset(n);   
     }
   }
};

struct SphereConstraint {
   template <typename NODE>
   bool violate(NODE n) {
     return norm(n.position() - Point(.5, .5, -.5)) < .15;
   }

   template <typename NODE>
   void reset(NODE n) {
     Point r_i = n.position() - Point(.5, .5, -.5);
     r_i = r_i/norm(r_i);

     n.position() = .15*r_i + Point(.5, .5, -.5);
     n.value().vel = n.value().vel - dot(n.value().vel, r_i)*r_i;
   }

   template <typename NODE>
   void operator()(NODE n) {
     if(violate(n)){
       reset(n);   
     }
   }

};


struct RemovingSphereConstraint {
   template <typename NODE>
   bool violate(NODE n) {
     return norm(n.position() - Point(.5, .5, -.5)) < .15;
   }

   template <typename NODE>
   void del(NODE n) {
     graph_->remove_node(n);
   }

   template <typename NODE>
   void operator()(NODE n) {
     if(violate(n)){
       del(n);   
     }
   }
  
   RemovingSphereConstraint(GraphType* graph) : graph_(graph) {}
  
  private:
    GraphType* graph_;

};

template <typename C1, typename C2>
struct CombinedConstraint {

  CombinedConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}

  template<typename NODE>
  NODE operator()(NODE n) {
    c1_(n);
    c2_(n);
    return n;
  }

  private:
    C1 c1_;
    C2 c2_;
};

template <typename C1, typename C2>
CombinedConstraint<C1,C2> make_combined_constraint(C1 c1, C2 c2){
  return CombinedConstraint<C1,C2>(c1,c2);
}

template <typename C1, typename C2, typename C3>
CombinedConstraint<C1,CombinedConstraint<C2, C3>> make_combined_constraint(C1 c1, C2 c2, C3 c3){
  return CombinedConstraint<C1,CombinedConstraint<C2, C3>>(c1,CombinedConstraint<C2,C3>(c2, c3));
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
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    return Point(0,0,-grav*n.value().mass) + spring_force(n) ;
  }

  template <typename NODE>
  Point spring_force(NODE n) {
    Point f {0,0,0};
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      NODE n2 = (*it).node2();
      if(n2 == n) {
        n2 = (*it).node1();
      }

      Point diff = n.position() - n2.position();
      double K = (*it).value().K;
      double L = (*it).value().L;
      f = f - K*((diff)/norm(diff))*(norm(diff) - L);
    } 
    return f;
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
  
  // Set initial conditions
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  NodeData initial_data {Point(0,0,0), 1.0/nodes.size()};
  for (auto n : nodes) {
    n.value() = initial_data;
  }
  
  for(auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    EdgeData initial_data {100.0, (*it).length()};
    (*it).value() = initial_data;
  }


  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);


  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();



  auto combined_constraint = make_combined_constraint(PlaneConstraint(), RemovingSphereConstraint(&graph));
  auto f = make_combined_force(MassSpringForce(), GravityForce(), DampingForce(1.0/nodes.size()));
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
        symp_euler_step(graph, t, dt, f, combined_constraint);

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

