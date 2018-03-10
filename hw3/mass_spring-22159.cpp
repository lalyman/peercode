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
#include <cmath>
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
  NodeData() : vel(0), mass(1.) {}
};
/** Custom structure of data to store with Nodes */
struct EdgeData {
    double K;       // edge elasticity 
    double len;     // edge length
    EdgeData() : K (100.), len (0.01) {}
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
  // check constraints
  // constraint(g,t);
  // Compute the t+dt position
  //std::cout << "t = " << t << std::endl;
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    assert(not std::isnan(norm_2(n.value().vel)));
    n.position() += n.value().vel * dt;
  }
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    assert(not std::isnan(norm_2(force(n,t))));
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  // check constraint again after velocity update
  constraint(g,t);

  return t + dt;
}

/*  gravity force. node mass must be correctlly specified before
    called */
struct GravityForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return n.value().mass*Point(0,0,-grav);
    }
}; 

/* spring force */
struct MassSpringForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        Point summ {0}, temp {0};
        for (auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
            assert(not std::isnan(norm_2(n.position())));
            assert(not std::isnan(norm_2((*it).node2().position())));
            /*std::cout << "n1 = " << n.index() <<\
                "  n2 = " << (*it).node2().index() <<\
                "  length = " << (*it).length() << std::endl;*/
            temp = (n.position()-(*it).node2().position())*\
            (1-(*it).value().len/(*it).length());
            summ += temp*(-(*it).value().K);
        }
        return summ; 
    }
}; 

/*  damping force. constant c_ is allowed to be customized
    when calling proper constructor. */
struct DampingForce {
    private:
        double c_ = 0.1; // damping constant 
    public:
    DampingForce(double c=0.1):c_(c) {}
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return -c_*n.value().vel;
    }
};

/*  zero force. it is used as as a placeholder. */
struct ZeroForce {
// use as default force 
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) n; (void) t;
        return Point(0,0,0);
    }
};

/*  combined force class. take in functors and calculate force accordingly
    for any nodes*/
template <typename F1=ZeroForce, typename F2=ZeroForce, typename F3=ZeroForce>
struct make_combined_force {
    private: 
        F1 force1_; 
        F2 force2_; 
        F3 force3_;
    public: 
        make_combined_force(F1 force1=ZeroForce(),\
                            F2 force2=ZeroForce(),\
                            F3 force3=ZeroForce())\
                        :force1_(force1), force2_(force2), force3_(force3){}
        template <typename NODE>
        Point operator()(NODE n,double t) {
            return force1_(n,t)+force2_(n,t)+force3_(n,t); 
        }
};

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
    // const double K = 100; // spring constant
    // const double L = 0.01; // rest length of all edges 
    (void) t;
    // fix (0,0,0) and (1,0,0) points 
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);

    Point summ {0}, temp {0};
    for (auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
        temp = (n.position()-(*it).node2().position())*\
            (1-(*it).value().len/(*it).length());
        summ += temp*(-(*it).value().K); 
    }
    summ += n.value().mass*Point(0,0,-grav);
    
    return summ;
  }
};

/* define constraints */
/*  plane constrant */ 
struct PlaneConstraint {
    private:
        double z_ = -0.75; 
    public: 
        PlaneConstraint(){}
        PlaneConstraint(double z): z_(z) {}
        void operator()(GraphType& graph,double t) {
            (void) t;
            for (auto it=graph.node_begin(); it!=graph.node_end(); ++it) {
                if (dot((*it).position(),Point(0,0,1))<z_) {
                    (*it).position()[2] = z_; 
                    (*it).value().vel[2] = 0.;
                }
            }
        }
};

/*  sphere constraint */
struct SphereConstraint {
    private: 
        Point c_ {0.5,0.5,-0.5};
        double r_ {0.15}; 
    public:
        SphereConstraint(){}
        SphereConstraint(const Point& c,double r): c_(c),r_(r) {}
        void operator()(GraphType& graph,double t) {
           (void) t; 
           for (auto it=graph.node_begin(); it!=graph.node_end(); ++it) {
                auto n = *it; 
                Point delta = n.position()-c_; 
                double rn = norm_2(delta);
                if (rn<r_) {
                    n.position() = (r_/rn*delta+c_);
                    n.value().vel -= \
                        dot(n.value().vel,delta)*delta/normSq(delta); 
                }
           }
        }
};

/*  sphere constraint. will remove nodes that hit the sphere. */
struct SphereRmConstraint {
    private:
        Point c_ {0.5,0.5,-0.5};
        double r_ {0.15};
    public:
        SphereRmConstraint(){}
        SphereRmConstraint(const Point& c,double r): c_(c),r_(r) {}
        void operator()(GraphType& graph,double t) {
           (void) t;
           Point delta {};
           for (auto it=graph.node_begin(); it!=graph.node_end(); ++it) {
                auto n = *it;
                delta = n.position()-c_;
                double rn = norm_2(delta);
                /*std::cout << n.position() << "  rn-r_ = " \
                        << rn-r_ << "  node# = " << graph.num_nodes() \
                        << std::endl;*/
				if (rn<r_) {
                    /*n.position() = (r_/rn*delta+c_);
                    n.value().vel -= \
                        dot(n.value().vel,delta)*delta/normSq(delta);*/
                    /*std::cout << n.position() << "  rn-r_ = " \
                        << rn-r_ << "  node# = " << graph.num_nodes() \
                        << std::endl;*/
                    graph.remove_node(n);
                }
           }
        }
};

/*  fixed node constraint. fix node position and set velocity to 0. */
struct ConstantNodeConstraint {
    private:
        std::vector<Point> fix_n_ {Point(0,0,0),Point(1,0,0)}; 
    public: 
        ConstantNodeConstraint(){}
        ConstantNodeConstraint(const std::vector<Point>& fix_n):\
            fix_n_(fix_n) {}
        void operator()(GraphType& graph,double t) {
            (void) t;
            double tol = ((*graph.edge_begin()).length()+\
                (*graph.edge_end()).length())/2000.;
            for (auto it=graph.node_begin(); it!=graph.node_end(); ++it) {
                auto n = *it; 
                for (unsigned i=0; i<fix_n_.size(); i++) {
                    if (norm_2(n.position()-fix_n_[i])<tol) {
                        n.position() = fix_n_[i];
                        n.value().vel = Point(0);
                    }
                }
            }
        }
}; 

/*  no constrait. used like a placeholder like ZeroForce. */
struct NoConstraint {
// use as default constraint 
    void operator()(GraphType graph,double t) 
    {
        (void) graph; (void) t;
    }
}; 

// combine constraints 
template
<typename C1=NoConstraint,typename C2=NoConstraint,typename C3=NoConstraint>
struct make_combined_constraint {
    private:
        C1 constraint1_; 
        C2 constraint2_;
        C3 constraint3_; 
    public: 
        make_combined_constraint(   C1 constraint1=NoConstraint(),\
                                    C2 constraint2=NoConstraint(),\
                                    C3 constraint3=NoConstraint())
        :   constraint1_(constraint1),\
            constraint2_(constraint2),\
            constraint3_(constraint3) {}
        void operator()(GraphType& graph, double t) {
            constraint1_(graph,t);
            constraint2_(graph,t);
            constraint3_(graph,t);
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
    // set mass using num of nodes and set initial length for edges
    for (auto nit=graph.node_begin(); nit!=graph.node_end(); ++nit) {
        (*nit).value().mass = 1./graph.num_nodes();
        for (auto eit=(*nit).edge_begin(); eit!=(*nit).edge_end(); ++eit) {
            (*eit).value().len = (*eit).length();
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        make_combined_force<GravityForce,MassSpringForce,DampingForce> \
            force ((GravityForce()),(MassSpringForce()),\
            (DampingForce(1./graph.num_nodes())));
        make_combined_constraint
            <PlaneConstraint,SphereRmConstraint,ConstantNodeConstraint>
            constraint 
            ((PlaneConstraint()),SphereRmConstraint(),ConstantNodeConstraint());
        symp_euler_step(graph, t, dt, force, constraint);
        
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
