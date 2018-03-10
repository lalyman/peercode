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

double damping_factor = 0.0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double spring_const;      //< spring constant
  double rest_length;     //<  rest length
  EdgeData() : spring_const(0), rest_length(0) {};
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using IncidentIterator = typename GraphType::incident_iterator;
using NodeIterator = typename GraphType::node_iterator;
using EdgeIterator = typename GraphType::edge_iterator;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G supports iteration over the nodes using methods
 *           node_begin() and node_ends()
 *           G::node_value_type has attributes vel (a point) and mass (a double)
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

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    if (n.position()!=Point(0,0,0)&& n.position()!= Point(1,0,0)){
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position()!=Point(0,0,0)&& n.position()!= Point(1,0,0)){
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and constraints.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in,out]     constraint Constraint object defining the constraint to apply
 *                           in the simulation
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G supports iteration over the nodes using methods
 *           node_begin() and node_ends()
 *           G::node_value_type has attributes vel (a point) and mass (a double)
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a constraint object called as @a constraint(@a g, @a t) where @a
 *           g is the graph and @a t is the current time. It must not return any
 *           value, but after it is applied, the graph is still vallid and now
 *           abides by the constraints.
 */
template <typename G, typename F,typename C>
double symp_euler_step(G& g, double t, double dt, F force,C& constraint) {
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
    (void) t;
    if (n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
      return Point(0,0,0);
    }else{
      Point p = Point(0,0,-grav*n.value().mass);
      for(IncidentIterator it = n.edge_begin();it!=n.edge_end();++it){
        p-=(*it).value().spring_const *((*it).node1().position()-(*it).node2().position())*(1-(*it).value().rest_length/(*it).length());
      }
    return p;
  }
  }
};

struct GravityForce{
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point p = Point(0,0,-grav*n.value().mass);
    return p;
  }
};

struct MassSpringForce{
  /** Return the Mass spring force applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point p = Point(0,0,0);
    for(IncidentIterator it = n.edge_begin();it!=n.edge_end();++it){
      p-=(*it).value().spring_const *((*it).node1().position()-(*it).node2().position())*(1-(*it).value().rest_length/(*it).length());
    }
    return p;
  }
};

struct DampingForce{
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point p = Point(-damping_factor*n.value().vel);
    return p;
  }
};

/**
 * Combination of several forces
 * @param force1 First force to apply
 * @param force2 second force to apply
 *
 * Creates an object that applies several forces successively
 * (such structs can be nested to combine more than 2 forces)
 */
template<typename F1,typename F2>
struct CompositeForce2{
  /** Attributes: 2 forces to apply successively */
  F1 f1;
  F2 f2;

  /**
   * Default constructor
   */
  CompositeForce2(F1 force1 ,F2 force2): f1(force1),f2(force2){};
  /**
   * Application of the forces
   * @param[in] n a node of the the graph
   * @param[in] t time of iteration
   * @return a point representing the combined force vector
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point p = Point(0,0,0);
    p+=f1(n,t);
    p+=f2(n,t);
    return p;
  }
};

/** Create a composite force
 *
 * @param[in]     f1  The first force to apply
 * @param[in]     f2  The second force to apply
 *
 * @return f a combined force object
 *
 * @tparam F1, F2 types of forces. Must implement an operator() that takes a node
 *         and time as an input and return the force applied to the node (a Point)
 *
 * @pre The forces objects f1,f2 must be valid and implement operator() that
 *      returns a point
 *
 * @post f implements f(n,t) where f is a Node and t is a double, the time.
 *       with f(n,t) = @a f1(n,t) + @a f2(n,t)
 */
template<typename F1,typename F2>
CompositeForce2<F1,F2> make_combined_force(F1 f1, F2 f2){
  CompositeForce2<F1,F2> f = CompositeForce2<F1,F2>(f1,f2);
  return f;
}

/** Create a composite forc
 *
 * @param[in]     f1  The first force to apply
 * @param[in]     f2  The second force to apply
 * @param[in]     f3  The third force to apply
 * @return f a combined force object
 *
 * @tparam F1, F2, F3 types of forces. Must implement an operator() that takes a node
 *         and time as an input and return the force applied to the node (a Point)
 *
 * @pre The forces objects f1,f2,f3 must be valid and implement operator() that
 *      returns a point
 *
 * @post f implements f(n,t) where f is a Node and t is a double, the time.
 *       with f(n,t) = @a f1(n,t) + @a f2(n,t) + @a f3(n,t)
 */
template<typename F1,typename F2, typename F3>
CompositeForce2<F1,CompositeForce2<F2,F3>> make_combined_force(F1 f1, F2 f2,F3 f3){
  CompositeForce2<F1,CompositeForce2<F2,F3>> f =CompositeForce2<F1,CompositeForce2<F2,F3>>(f1,make_combined_force(f2, f3));
  return f;
}

/**
 * Basic plane constraint, forbids a half of the space defined by a normal vector
 * and an upper bound for the dot product with the normal vector. Violating nodes
 * are put on the surface of the plane
 */
struct PlaneConstraint{
  Point normal_vector; //the vector normmal to the plan
  double upperbound; // the upper bound of the dot product for the forbidden
  //half space
  /** Basic constructor*/
  PlaneConstraint() : normal_vector(0), upperbound(0){ };
  /** Regular constructor
  * @param[in] p the normal direction
  * @param[in] d the upperbound of the forbidden region in the equation
  *
  * @pre The normal direction cannot be null(the cosntraints are non trivial)
  * @post The forbidden zone is the same, the normal vector has norm 1
  *
  * We normalize the direction in the equation to save computations later
  */
  PlaneConstraint(Point p,double d) : normal_vector(p), upperbound(d){
    assert (p!= Point(0,0,0)); //We do not want to consider the trivial case
    upperbound*=norm_2(p);
    normal_vector/=norm_2(p);
  };
  /**
   * Determines if a node violates the constraint at time t
   * @param  n node to check
   * @param  t time of computation
   * @return   true if the node violates the constraint, false if not
   */
  template <typename NODE>
  bool violates(NODE n, double t){
    (void) t;
    return (inner_prod(n.position(),normal_vector)<upperbound);
  }

  /**
   * Fixes a node parameters if it is in the forbidden zone
   * @param[in,out]  n node fix
   * @param  t time of computation
   *
   * @post The position of the node has been projected on the plane
   * @post The velocity of the node has no component normal to the plane
   */
  template <typename NODE>
  void fix(NODE n,double t){
    (void) t;
    n.value().vel-=inner_prod(n.value().vel,normal_vector)*normal_vector;
    n.position()-=(inner_prod(n.position(),normal_vector)-upperbound)*normal_vector;
  }
  /* *Apply the constraint to a graph
  * @param[in,out]  g a graph on which to apply the constraint
  * @param[in]      t the time of iteration
  *
  * @tparam GRAPH the type of the graph
  *
  * @return 1 if a node has been deleted, 0 if not
  *
  * @pre  @a g is a valid graph
  * @post @a g is a valid graph with the same sumber of nodes and for all
  *          node n in @a g,!(violates(new n, @a t)).
  * @post for all node n in @a g, if !(violates(old n, @a t)), old n = new n
  */
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    for(auto ni = g.node_begin();ni!= g.node_end();++ni){
      if(violates(*ni,t)){
        fix(*ni,t);
      }
    }
    return 0;
  }
  /* * Update the indexes stored in the Constraint.
  * Method has no effect here, but it must be implemented
  */
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
  (void) t;
  (void) g;
  }
};


/**
 * Basic plane constraint, forbids a half of the space defined by a normal vector
 * and an upper bound for the dot product with the normal vector. Violating nodes
 * are destroyed
 */
struct PlaneDestroyerConstraint{
  Point normal_vector; //the vector normmal to the plan
  double upperbound; // the upper bound of the dot product for the forbidden
  //half space
  /** Basic constructor*/
  PlaneDestroyerConstraint() : normal_vector(0), upperbound(0) { };
  /** Regular constructor
  * @param[in] p the normal direction
  * @param[in] d the upperbound of the forbidden region in the equation
  *
  * @pre The normal direction cannot be null(the cosntraints are non trivial)
  * @post The forbidden zone is the same, the normal vector has norm 1
  *
  * We normalize the direction in the equation to save computations later
  */
  PlaneDestroyerConstraint(Point p,double d) : normal_vector(p), upperbound(d){
    assert (p!= Point(0,0,0)); //We do not want to consider the trivial case
    upperbound*=norm_2(p);
    normal_vector/=norm_2(p);
  };

  /**
   * Determines if a node violates the constraint at time t
   * @param  n node to check
   * @param  t time of computation
   * @return   true if the node violates the constraint, false if not
   */
  template <typename NODE>
  bool violates(NODE n, double t){
    (void) t;
    return (inner_prod(n.position(),normal_vector)<upperbound);
  }

  /* *Apply the constraint to a graph
  * @param[in,out]  g a graph on which to apply the constraint
  * @param[in]      t the time of iteration
  *
  * @tparam GRAPH the type of the graph
  *
  * @return 1 if a node has been deleted, 0 if not
  *
  * @pre  @a g is a valid graph
  * @post @a g is a valid graph and for all node n in old @a g,
  *             !(violates(n, @a t)), n.position() and n.value() are preserved
  *          or n is not a  node of new @a g
  */
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    int result {0};
    auto ni = g.node_begin();
    while(ni!= g.node_end()){
      if(violates(*ni,t)){
        ni = g.remove_node( ni);
        result=1;
      }else{
        ++ni;
      }
    }
    return result;
  }

  /* * Update the indexes stored in the Constraint.
  * Method has no effect here, but it must be implemented
  */
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
  (void) t;
  (void) g;
  }
};

/**
 * Basic ball constraint, forbids a ball defined by its center and radius
 * Violating nodes are put on the surface of the ball
 */
struct BallConstraint{
  Point center; //center of the ball
  double radius; //radius
  /** Basic constructor*/
  BallConstraint() : center(0), radius(0){};
  BallConstraint( Point p ,double d) : center(p), radius(d){};

  /**
   * Determines if a node violates the constraint at time t
   * @param  n node to check
   * @param  t time of computation
   * @return   true if the node violates the constraint, false if not
   */
  template <typename NODE>
  bool violates(NODE n, double t){
    (void) t;
    return (norm_2(center - n.position())<radius);
  }

  /**
   * Fixes a node parameters if it is in the forbidden zone
   * @param[in,out]  n node fix
   * @param  t time of computation
   *
   * @post The position of the node has been projected on the ball
   * @post The velocity of the node has no component normal to the ball
   */
  template <typename NODE>
  void fix(NODE n,double t){
    (void) t;
    Point ray {n.position()-center};
    ray/=norm_2(ray);
    n.value().vel-=inner_prod(n.value().vel,ray)*ray;
    n.position()-=(norm_2(n.position()-center)-radius)*ray;
  }

  /* *Apply the constraint to a graph
  * @param[in,out]  g a graph on which to apply the constraint
  * @param[in]      t the time of iteration
  *
  * @tparam GRAPH the type of the graph
  *
  * @return 1 if a node has been deleted, 0 if not
  *
  * @pre  @a g is a valid graph
  * @post @a g is a valid graph with the same sumber of nodes and for all
  *          node n in @a g,!(violates(new n, @a t)).
  * @post for all node n in @a g, if !(violates(old n, @a t)), old n = new n
  */
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    for(auto ni = g.node_begin();ni!= g.node_end();++ni){
      if(violates(*ni,t)){
        fix(*ni,t);
      }
    }
    return 0;
  }

  /* * Update the indexes stored in the Constraint.
  * Method has no effect here, but it must be implemented
  */
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
  (void) t;
  (void) g;
  }
};

/**
 * Basic ball constraint, forbids a ball defined by its center and radius
 * Violating nodes are destroyed
 */
struct BallDestroyerConstraint{
  Point center; //center of the ball
  double radius; //radius
  /** Basic constructor*/
  BallDestroyerConstraint() : center(0), radius(0) {};
  BallDestroyerConstraint( Point p ,double d) : center(p), radius(d) {};

  /**
   * Determines if a node violates the constraint at time t
   * @param  n node to check
   * @param  t time of computation
   * @return   true if the node violates the constraint, false if not
   */
  template <typename NODE>
  bool violates(NODE n, double t){
    (void) t;
    return (norm_2(center - n.position())<radius);
  }

  /* *Apply the constraint to a graph
  * @param[in,out]  g a graph on which to apply the constraint
  * @param[in]      t the time of iteration
  *
  * @tparam GRAPH the type of the graph
  *
  * @return 1 if a node has been deleted, 0 if not
  *
  * @pre  @a g is a valid graph
  * @post @a g is a valid graph and for all node n in old @a g,
  *             !(violates(n, @a t)), n.position() and n.value() are preserved
  *          or n is not a  node of new @a g
  */
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    int result {0};
    auto ni = g.node_begin();
    while(ni!= g.node_end()){
      if(violates(*ni,t)){
        ni = g.remove_node(ni);
        result = 1;
      }else{
        ++ni;
      }
    }
    return result;
  }
  /* * Update the indexes stored in the Constraint.
  * Method has no effect here, but it must be implemented
  */
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
  (void) t;
  (void) g;
  }
};


/**Void constraint, has no effect but shows the prototype of a constraint - useful for testing*/
struct VoidConstraint{
  VoidConstraint(){};
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    (void) t;
    (void) g;
    return 0;
  }
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
  (void) t;
  (void) g;
  }
};

/* * Fixes a point at a given position. When possible, does not iterate over the
* nodes to find the elements but rather uses indexes to access it directly
*/
struct FixedPointConstraint{
  Point fixed_coord; //coordinate of the fixed point
  Point::size_type index; //index of the fixed point
  int initialized;

  /* * Basic constructor*/
  FixedPointConstraint() : fixed_coord(0), index(0),initialized(0) {};
  FixedPointConstraint( Point p) : fixed_coord(p), index(0),initialized(0) {};
  template <typename GRAPH>
  int operator()(GRAPH& g,double t){
    (void) t;
    if (initialized){
      g.node(index).position() = Point(fixed_coord);
      g.node(index).value().vel= Point(0,0,0);
    }else{
      for(auto ni = g.node_begin();ni!= g.node_end();++ni){
        if((*ni).position()==fixed_coord){
           index = (*ni).index();
           initialized = 1;
        }
      }
    }
    g.node(index).value().vel= Point(0,0,0);
    return 0;
  }
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
    initialized *=0;
    (*this)(g,t);
  }


};
/* * Creates an aggregate constraint
We will use nested versions of this struct to define combined constraints with a cardinality larger than 2
 */
template<typename C1,typename C2>
struct CompositeConstraint2{
  C1 c1;
  C2 c2;
  //Basic constructor
  CompositeConstraint2(  C1 ct1, C2 ct2) : c1(ct1), c2(ct2){};

  /* *Apply the constraints to a graph
  * @param[in,out]  g a graph on which to apply the constraint
  * @param[in]      t the time of iteration
  *
  * @tparam GRAPH the type of the graph
  *
  * @return a positive int if a node has been deleted, 0 if not
  *
  * @pre  @a g is a valid graph
  * @post @a g repects the constraints c1 and c2
  */
  template <typename GRAPH>
  int operator()(GRAPH& g,double t) {
    int result {0};

    //We apply the constraints
    result+=c1(g,t);
    result+=c2(g,t);
    // If a node has been deleted, then we need to update the indexes stroed in the constraint
    if (result>0){
      std::cout << "Updating the indexes, t = "<<t << '\n';
      c1.update_indexes(g,t);
    }
    return result;
  }
  /** Update the indexes */
  template <typename GRAPH>
  void update_indexes(GRAPH& g,double t){
    c1.update_indexes(g,t);
    c2.update_indexes(g,t);
  }
};

/**
 * Make combined constraints
 * @param c1 a constraint
 * @param c2 a constraint
 *
 * @return a constraint that applies c1, then c2, and then, if needed, reindexes c1
 */
template<typename C1,typename C2>
CompositeConstraint2<C1,C2> make_combined_constraint(C1 c1, C2 c2){
  CompositeConstraint2<C1,C2> c = CompositeConstraint2<C1,C2>(c1,c2);
  return c;
}

/** Same function for 3 constraints */
template<typename C1,typename C2, typename C3>
CompositeConstraint2<C1,CompositeConstraint2<C2,C3>> make_combined_constraint(C1 c1, C2 c2,C3 c3){
  CompositeConstraint2<C1,CompositeConstraint2<C2,C3>> c = CompositeConstraint2<C1,CompositeConstraint2<C2,C3>>(c1,make_combined_constraint(c2, c3));
  return c;
}

/** Same function for 4 constraints */
template<typename C1,typename C2, typename C3, typename C4>
CompositeConstraint2<C1,CompositeConstraint2<C2,CompositeConstraint2<C3,C4>>> make_combined_constraint(C1 c1, C2 c2,C3 c3,C4 c4){
  CompositeConstraint2<C1,CompositeConstraint2<C2,CompositeConstraint2<C3,C4>>>  c = CompositeConstraint2<C1,CompositeConstraint2<C2,CompositeConstraint2<C3,C4>>> (c1,make_combined_constraint(c2,c3,c4));
  return c;
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for(NodeIterator ni = graph.node_begin();ni!=graph.node_end();++ni){
    (*ni).value().mass = double(1)/double(graph.num_nodes());
  }

  for(EdgeIterator ei = graph.edge_begin();ei!=graph.edge_end();++ei){
    (*ei).value().rest_length = (*ei).length();
    (*ei).value().spring_const = double(100);
  }
  damping_factor = 1.0/double(graph.num_nodes());
  //std::cout << spring_rest_length<< '\n';

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
      auto force = make_combined_force(GravityForce(), MassSpringForce(),DampingForce());
      //auto constraint = BallConstraint(Point(0.5,0.5,-0.5),0.15);
      // auto constraint = PlaneConstraint(Point(0,0,1),-.75);
      auto constraint = make_combined_constraint(FixedPointConstraint(Point(1,0,0)),FixedPointConstraint(Point(0,0,0)),BallConstraint(Point(0.5,0.5,-0.5),0.15),PlaneDestroyerConstraint(Point(0,0,1),-.75));
      //auto constraint = FixedPointConstraint(Point(1,0,0));
      unsigned num_nodes {graph.num_nodes()};
      unsigned num_edges {graph.num_edges()};
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {

        symp_euler_step(graph, t, dt,force,constraint );


        // Update viewer with nodes' new positions

        // If the number of nodes has changed, we need to update the viewer completely
        if ((num_nodes!=graph.num_nodes()) || (num_edges!=graph.num_edges())){
          viewer.clear();
          node_map.clear();
          viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
          viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        }else{
          viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        }
        num_nodes = graph.num_nodes();
        num_edges = graph.num_edges();

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
