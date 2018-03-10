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
double K;
double c;
double L;
std::map<unsigned,Point> anchors; 

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Struct to store the data for edges. These are currently its length and spring
// constant. In this file the spring constant is fixed for all edges though.
struct EdgeData{
    double len;
    double K;
    EdgeData() :len(0), K(0){}
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
 * @tparam C is a function object called as @a constraint(g,@a t),
 *           where g is the grap the operation takes place on ant @a t is the
 *           current time. It acts on the graph which is passed by reference and
 *           therefor does not return anything.
 */
template <typename G, typename F,typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position

  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  //Check if the constraints are violated.
  constraint(g,t);

  // Compute the t+dt velocity
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
   * model that by returning a zero-valued force. 
   *
   * The forces are vectors represented by Point objects.
   * Gravitational force is the mass of the node times the gravitational
   * constant in the z direction.
   * The spring force is calculated by iterating over all the adjacent nodes,
   * calculating the spring force for each and summing the result.
   * */
  template <typename NODE>
  Point operator()(NODE n, double) {
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)){
        return Point(0,0,0);
    }
    Point Fgrav=Point(0.0,0.0,-1*grav*n.value().mass);
    Point Fspring=Point(0.0,0.0,0.0);
    auto adj_iter = n.edge_begin();
    while(adj_iter!=n.edge_end()){
        Point dist=n.position()-(*adj_iter).node2().position();
        double Enorm = norm(dist);
        Fspring+=-K*(dist/Enorm)*(Enorm-L);
        ++adj_iter;
    }
    return Fgrav+Fspring;
  }
};

/* The force on a node due to gravity.
 * This function returns a Point with a 
 * z-component equal to the mass of the node times the gravitational constant.
 */
 
struct GravityForce{
    template <typename NODE>
    Point operator()(NODE& n,double ){
        return Point(0.0,0.0,-1*grav*n.value().mass);
    } 
};

/* The force on a node due to the spring forces of all adjacent nodes.
 * The force is calculated by iterating over all incident edges, calculating
 * the spring force for each and returning the sum of all the combined forces.
 * The rest length and spring constant for each edge are stored in edge.value().
 */
struct MassSpringForce{
    template <typename NODE>
    Point operator()(NODE& n,double ){
        Point Fspring=Point(0.0,0.0,0.0);
        for(auto adj_iter = n.edge_begin();adj_iter!=n.edge_end();++adj_iter){
            Edge e = *adj_iter;
            Point dist=n.position()-e.node2().position();
            double Enorm = norm(dist);
            Fspring+=-(e.value().K)*(dist/Enorm)*(Enorm-e.value().len);
        }
        return Fspring;
    }
};
/*The damping force on a node
 *This force is anti-parallel to the velocity of the node.
 *The return value is the negative of the velocity times the damping constant.
 */
struct DampingForce{
    template <typename NODE>
    Point operator()(NODE& n,double ){
       return -1*c*n.value().vel;
    }
};
/* This force simulates a impact on the center of the sheet at t=1.
 * For this to work it assumes that the sheet has not fallen yet (i.e. gravity is turned of.)
 * The force works best if all four corners are fixed.
 *
 */
struct ImpactForce{
    ImpactForce(){
        center_=Point(0.5,0.5,0.0);
        radius_=0.05;
    }
    template<typename NODE>
    Point operator()(NODE& n,double t){
        Point line = n.position()-center_;
        if(t>=0.999&&t<1.001&&norm(line)<=radius_){
            return Point(0.0,0.0,-3.0); 
        }
        else{
            return Point(0.0,0.0,0.0);
        }
    }
    private:
        Point center_;
        double radius_;
};
/*A struct that allows two forces to be applied simulatneously on a node.
 *The constructor takes in two function objects that both must take in a node and a double
 *as argument and return a Point.
 *The defined operator() applies both forces on the node and returns the sum of the two results.
 *
 */
template <typename F1,typename F2>
struct combined_2force{
    combined_2force(F1 f1, F2 f2){
        f1_= f1;
        f2_= f2;
    };
    template <typename NODE>
   Point operator()(NODE& n,double t){
        return(f1_(n,t)+f2_(n,t));
    }
    private:
        F1 f1_;
        F2 f2_;
};
/*A struct that allows three forces to be applied simulatneously on a node.
 *The constructor takes in three function objects that all must take in a node and a double
 *as argument and return a Point.
 *The defined operator() applies all forces on the node and returns the sum of the three results.
 *
 */
template <typename F1,typename F2,typename F3=int>
struct combined_3force{
    combined_3force(F1 f1, F2 f2,F3 f3){
        f1_=f1;
        f2_=f2;
        f3_=f3;
    };
    combined_3force(F1 f1, F2 f2){
        f1_= f1;
        f2_= f2;
    };
    template <typename NODE>
    Point operator()(NODE& n,double t){
        return(f1_(n,t)+f2_(n,t)+f3_(n,t));
    }
    private:
        F1 f1_;
        F2 f2_;
        F3 f3_;
};

/*The two function definitions below allow for the construnction of the combined
 * forces structs above without having to template the function call. This makes
 * the implementation more accessible.
 * NOTE the function call can only accept two or three function objects.
 *
 */
template <typename F1,typename F2>
combined_2force<F1,F2> make_combined_force(const F1& f1,const F2& f2){
    return(combined_2force<F1,F2>{f1,f2});
}

template <typename F1,typename F2,typename F3>
combined_3force<F1,F2,F3> make_combined_force(const F1& f1,const F2& f2,const F3& f3){
    return(combined_3force<F1,F2,F3>{f1,f2,f3});
}


/* This struct simulates a flat surface at z=-0.75 which the nodes cannot pass
 * through.
 * If a node passes the surface is position gets updated to lie on the table and
 * the z-component of the velocity is set to zero.
 */
struct ConstraintTable{
    template<typename G>
    void operator()(G& graph,double ){
        for(auto fit= graph.node_begin();fit!=graph.node_end();++fit){
            auto n = *fit;
            if (n.position().z<-0.75){
                n.position().z=-0.75;
                n.value().vel.z=0;
            } 
        }
    }
};

/* This struct fixes the certain points in the graph. Which points are fixed and
 * to which position are stored in the map anchors where the keys are the
 * indices of nodes to be fixed and the value is the psotion they are fixed to.
 * anchors is a global variable that is initializd once the graph has been
 * constructed and fixed afterwards.  
 */
struct FixedPoints{
    template<typename G>
    void operator()(G& graph,double ){
        for(auto fit= graph.node_begin();fit!=graph.node_end();++fit){
            auto n = *fit;
            if(anchors.count(n.index())!=0){
                n.position()=anchors[n.index()];
                n.value().vel=Point(0,0,0);
            }
        }
    }
};

/* This struct is very similar to FixedPoints. However if t>0 we no longer 
 * iteratate through the nodes to enforce the conbstraints. This means that the
 * anchors are released at t=1 and no longer hold up the sheet.
 */
struct ReleaseSheet{
    template<typename G>
    void operator()(G& graph,double t){
        if(t>1.0){
            return;
        }
        for(auto fit= graph.node_begin();fit!=graph.node_end();++fit){
            auto n = *fit;
            if(anchors.count(n.index())!=0){
                n.position()=anchors[n.index()];
                n.value().vel=Point(0,0,0);
            }
        }
    }
};
/* This struct simulates a sphere which the nodes cannot pass through.
 * If a node passes through the sphere we find the projection of the node on the
 * sphere and update the position of the node to that. find the velocity of the
 * node perpendicular to the sphere and subtract that from the velocity.
 * The dimensions of the sphere are fixed during initialization.
 */
struct ConstraintSphere1{
    ConstraintSphere1(){
        center_=Point(0.5,0.5,-0.5);
        radius_=0.15;
    }
    template<typename G>
    void operator()(G& graph,double ){
        for(auto fit= graph.node_begin();fit!=graph.node_end();++fit){
            auto n = *fit;
            Point line = n.position()-center_;
            if(norm(line)<=radius_){
                Point scaled_line=line/norm(line);
                Point scale = radius_*scaled_line;
                n.position()=scale+center_;
                n.value().vel-=dot(n.value().vel,(scaled_line))*(scaled_line);
            }
        } 
    }
    private:
    Point center_;
    double radius_;

};

/* This struct simulates a sphere which removes nodes when they touch it.
 * We use the same check as the previous constraint to see if the node is passed
 * the sphere. If the node is inside the sphere we remove it.
 * Else we increment the node iterator.
 * Since removing a node pops the back of the vector and the anchors are defined
 * by index we need to make sure that we do not remove an anchor without
 * updating @a anchors. When the index of the last element in the vector
 * is in @a anchors we add a new pair to anchors with the index the last element
 * is moved to and the position we want fixed..
 * */
struct ConstraintSphere2{
    ConstraintSphere2(){
        center_=Point(0.5,0.5,-0.5);
        radius_=0.15;
    }
    template<typename G>
    void operator()(G& graph,double ){
        for(auto fit= graph.node_begin();fit!=graph.node_end();){
            auto n = *fit;
            Point line = n.position()-center_;
            if(norm(line)<=radius_){
                if (anchors.count(graph.num_nodes()-1)!=0){
                    anchors[n.index()]=anchors[graph.num_nodes()-1];
                    anchors.erase(graph.num_nodes()-1);
                }
                graph.remove_node(n);
            }            
            else{
                ++fit;
            }
        }
    }
    private:
    Point center_;
    double radius_;
};

/*A struct that allows three constriants to be applied simulatneously on a node.
 *The constructor takes in three function objects that all must take in a node and a double
 *as argument and return a Point.
 *The defined operator() applies the constraints on the graph consecutively.
 *
 */
template <typename C1,typename C2,typename C3>
struct combined_3const{
    combined_3const(C1 c1, C2 c2, C3 c3){
        c1_= c1;
        c2_= c2;
        c3_= c3;
    };
    template <typename G>
   void operator()(G& g,double t){
       c1_(g,t);
       c2_(g,t);
       c3_(g,t);

    }
    private:
        C1 c1_;
        C2 c2_;
        C3 c3_;
};

/*A struct that allows two constriants to be applied simulatneously on a node.
 *The constructor takes in two function objects that all must take in a node and a double
 *as argument and return a Point.
 *The defined operator() applies the constraints on the graph consecutively.
 *
 */
template <typename C1,typename C2>
struct combined_2const{
    combined_2const(C1 c1, C2 c2){
        c1_= c1;
        c2_= c2;
    };
    template <typename G>
   void operator()(G& g,double t){
       c1_(g,t);
       c2_(g,t);
    }
    private:
        C1 c1_;
        C2 c2_;
};

/*The two function definitions below allow for the construnction of the combined
 * constraint structs above without having to template the function call. This makes
 * the implementation more accessible.
 * NOTE the function call can only accept two or three constraint objects.
 */
template <typename C1,typename C2>
combined_2const<C1,C2> make_combined_const(const C1& c1,const C2& c2){
    return combined_2const<C1,C2>{c1,c2};
}

template <typename C1,typename C2,typename C3>
combined_3const<C1,C2,C3> make_combined_const(const C1& c1,const C2& c2,const C3& c3){
    return combined_3const<C1,C2,C3>{c1,c2,c3};
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
  // Set the global spring and damping, and rest length constaints.
  K=100.0;
  c=1.0/graph.num_nodes();
  L=norm(graph.edge(0).node1().position()-graph.edge(0).node2().position());
  
  // Set the mass of the nodes.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().mass=(1.0/graph.num_nodes());
  }

  // Set the rest length of the edges to their initial length.
  // Set the spring constant for each node to be K.
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    e.value().len=(norm(e.node1().position()-e.node2().position()));
    e.value().K=K;
  }

  
  // Fix the points that should be kept fixed throughout the simulation.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    if(n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
        anchors[n.index()]=n.position();
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

  auto force3=make_combined_force(GravityForce(),DampingForce(),MassSpringForce());
  auto const3=make_combined_const(ReleaseSheet(),ConstraintTable(),ConstraintSphere2());
  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        symp_euler_step(graph, t, dt,force3,const3);
        viewer.clear();
        node_map.clear();

         // Update viewer with nodes new positions and new edges
         viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map );
         viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map );
        // Update viewer with nodes' new positions
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
