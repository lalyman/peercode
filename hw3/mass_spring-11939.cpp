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

#include <iostream>
#include <stdexcept>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  double c ;        // Damping constant
  NodeData() : vel(0), mass(.1), c(0.1) {}
};

struct EdgeData {
  double K;       //< Spring Constant
  double L;     //< Rest Length
  EdgeData() : K(100), L(0.0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using IncidentIterator = typename GraphType::IncidentIterator;
using EdgeIterator = typename GraphType::EdgeIterator;

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
template <typename G, typename F,typename C >
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


//------------------------------- BEGIN FORCES --------------------------------------------//

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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);
    Point xi {n.position()};
    Point total_force {n.value().mass*Point(0,0,-grav)}; //initialize with gravity
    for (IncidentIterator j=n.edge_begin(); j!=n.edge_end(); ++j)
    {
        Point xj {(*j).node2().position()};
        double dist {norm(xi-xj)};
        double K  {(*j).value().K};
        double L {(*j).value().L};
        total_force += -K*(xi-xj)/dist*(dist-L);
    }
    (void)t;  
    return total_force;
  }
};

/* Spring force */
struct MassSpringForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        Point xi {n.position()};
        Point total_force {Point(0,0,0)};
        for (IncidentIterator j=n.edge_begin(); j!=n.edge_end(); ++j)
        {
            Point xj {(*j).node2().position()};
            double dist {norm(xi-xj)};
            double K  {(*j).value().K};
            double L {(*j).value().L};
            total_force += -K*(xi-xj)/dist*(dist-L);
        }    
        return total_force;
    }
};

/* Force of gravity*/
struct GravityForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        return n.value().mass*Point(0,0,-grav);
    }
};

/* Damping force, strength is defined by a parameter in node value struct: C*/
struct DampingForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        return -n.value().vel*n.value().c*Point(1,1,1);
    }
};

//Make a CombinedForce struct that will function same as
//F1 or F2 individually, akin to std::pair
template <typename F1, typename F2>
struct CombinedForce {
    F1 m_force1;
    F2 m_force2;
    template <typename NODE>
    Point operator()(NODE n, double t) {
        return m_force1(n,t) + m_force2(n,t);
    }
};

//Takes in 2 force structs
//Return an instance of CombinedForce 
template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1& f1,F2& f2) {
    CombinedForce<F1,F2> combined_force;
    combined_force.m_force1 = f1;
    combined_force.m_force2 = f2;
    return combined_force;
}

//Take in 3 forces, return CombinedForce
template <typename F1, typename F2,typename F3>
CombinedForce< CombinedForce<F1,F2>,F3 > make_combined_force(F1& f1, F2& f2, F3& f3) {
    CombinedForce< CombinedForce<F1,F2>,F3 > combined_force;
    combined_force.m_force1.m_force1 = f1;
    combined_force.m_force1.m_force2 = f2;
    combined_force.m_force2 = f3;
    return combined_force;
}


//------------------------------- BEGIN CONSTRAINTS --------------------------------------------//

/*  Constrains 2 nodes passed in to have fixed points (0,0,0) 
 *  and (1,0,0).  The user should pass the correct nodes 
 *  to be constrained so that a node search does not need to be 
 *  performed on each constraint call.
 *  */
struct ConstantNodeConstraint {
    std::vector<Node> m_nodes;
    ConstantNodeConstraint() {}
    ConstantNodeConstraint(std::vector<Node> nodes) {m_nodes=nodes;}
    template <typename GRAPH>
    void operator()(GRAPH g, double t) {
        (void)t;(void)g;
        //std::cout<<"Constraining..."<<std::endl;
        m_nodes[0].position() = Point(0,0,0);
        m_nodes[1].position() = Point(1,0,0);
        return;
    }
};

/* Constraint that places a plane at z=-0.75. 
 * Nodes positions are fixed to be the nearest point on the plane surface
 * and velocities normal to the plane are set to zero. */
struct PlaneConstraint {
    template <typename GRAPH>
    void operator()(GRAPH& g, double t) {
        (void)t;
        for (auto nit=g.node_begin(); nit!=g.node_end();++nit)
        {
            if ((*nit).position()[2]<-0.75)
            {   
                (*nit).position()[2] = (-0.75);
                (*nit).value().vel[2]=0;
            }
        }
        return;
    }
};

/* Constraint that places a sphere at (0.5,0.5.-0.5) with radius 0.15 
 * Nodes positions are fixed to be the nearest point on sphere surface
 * and velocities normal to the sphere surface are set to zero. */
struct SphereConstraint {
    Point center;
    double radius;
    SphereConstraint(): center(Point(0.5,0.5,-0.5)), radius(0.15) {}
    template <typename GRAPH>
    void operator()(GRAPH& g, double t) {
        (void)t;
        double dist;
        Point unit_dist_vector;
        for (auto nit=g.node_begin(); nit!=g.node_end();++nit)
        {
            dist = norm((*nit).position()-center);
            if (dist < radius)
            {
                unit_dist_vector = ((*nit).position()-center)/dist;
                (*nit).position() = unit_dist_vector*radius+center;
                (*nit).value().vel-= dot(unit_dist_vector,(*nit).value().vel)*unit_dist_vector ;
            }
        }
        return;
    }
};

/* Constraint that deletes all nodes and adjacent edges
 * that touch a sphere at (0.5,0.5.-0.5) with radius 0.15 */
struct SphereDeleteConstraint {
    Point center;
    double radius;
    SphereDeleteConstraint(): center(Point(0.5,0.5,-0.5)), radius(0.15) {}
    template <typename GRAPH>
    void operator()(GRAPH& g, double t) {
        (void)t;
        double dist;
        Point unit_dist_vector;
        for (auto nit=g.node_begin(); nit!=g.node_end();)
        {
            dist = norm((*nit).position()-center);
            if (dist < radius)
            {   
                nit = g.remove_node(nit);
            }
            else { ++nit; }
        }
        return;
    }
};

//Make a CombinedConst struct that will function same as
//F1 or F2 individually, akin to std::pair
template <typename C1, typename C2>
struct CombinedConst {
    C1 m_const1;
    C2 m_const2;
    template <typename GRAPH>
    void operator()(GRAPH& g, double t) {
        m_const1(g,t);
        m_const2(g,t);
        return;
    }
};

//Takes in 2 force structs
//Return an instance of CombinedConst
template <typename C1, typename C2>
CombinedConst<C1,C2> make_combined_constraint(C1& c1,C2& c2) {
    CombinedConst<C1,C2> combined_const;
    combined_const.m_const1 = c1;
    combined_const.m_const2 = c2;
    return combined_const;
}

//Take in 3 forces, return CombinedConst
template <typename C1, typename C2,typename C3>
CombinedConst< CombinedConst<C1,C2>,C3> make_combined_constraint(C1& c1,C2& c2,C3& c3) {
    CombinedConst< CombinedConst<C1,C2>,C3 > combined_const;
    combined_const.m_const1.m_const1 = c1;
    combined_const.m_const1.m_const2 = c2;
    combined_const.m_const2 = c3;
    return combined_const;
}


//------------------------------- BEGIN MAIN --------------------------------------------//
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double graph_size = double(graph.size());
  for(auto it=graph.node_begin();it!=graph.node_end();++it)
  {
      (*it).value().vel = Point(0,0,0);
      (*it).value().mass = 1.0/graph_size;
      (*it).value().c = 1.0/graph_size;
  }

  for(auto eit=graph.edge_begin();eit!=graph.edge_end();++eit)
  {
      (*eit).value().K = 100;
      (*eit).value().L = (*eit).length();

  }
  
  //graph.remove_edge(graph.edge_begin());
  //std::cout << "A"<<std::endl;
  std::vector<Node> fixed_nodes(2);
  for (auto nit =graph.node_begin(); nit!=graph.node_end(); ++nit)
  {

    if ( (*nit).position() == Point(0,0,0) ) 
        {
            fixed_nodes[0] = (*nit);} 
    else if ( (*nit).position() == Point(1,0,0) )
        {fixed_nodes[1] = (*nit);} 
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
     
      auto g_force = GravityForce();
      auto ms_force = MassSpringForce();
      auto d_force = DampingForce();
      auto c_force = make_combined_force(g_force,ms_force,d_force);

      auto cn_const = ConstantNodeConstraint( fixed_nodes );
      auto p_const = PlaneConstraint();
      auto s_const = SphereConstraint();
      auto sd_const = SphereDeleteConstraint();
      auto c_const = make_combined_constraint(cn_const,sd_const,p_const);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {


        //symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, c_force, c_const );
        
        //Add following lines if including node removal
        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
