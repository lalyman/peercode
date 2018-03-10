/* Nicole Schiavone
 * CME 212 - HW 2
 */

/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 *
 * Includes two static constants that represent positions of Nodes that can
 * be fixed by constraints. Asserts that Nodes with these posiitions must be
 * in the graph. Assert was used in this case to "fail with grace" and avoid
 * undefined behavior if the search for the Nodes with these positions failed.
 *
 * Main includes commented out calls that combine different forces and
 * constraints that were used for testing various problems in the assignment
 */

#include <chrono>
#include <fstream>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

//Positions to constrain as fixed
static const Point FixedNode1Pos = Point(0,0,0);
static const Point FixedNode2Pos = Point(1,0,0);

//Spring constant for all edges
static const double givenSpringConst = 100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double springConst;
  double restLen;
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports structs with a member vel to represent
 *           the node velocity and a member mass to represent the node mass
 * @tparam F is a function object called as @a force(@a n, @a t),
 *           where @a n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node @a n at time @a t.
 * @tparam C is a function object called as @a constraint(@a g, @a t),
 *           where @g is the graph and @a t is the current time.
 *           @a constraint has no return type
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

  //Enforce constraints
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
    (void) t; //Quiet compiler warning

    //Check for fixed points first - hard coded for problems #1 and #2
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    //Gravitational force for the node
    Point gravForce = n.value().mass * Point(0,0,-grav);

    //Spring force for the node
    Point springForce = Point(0);
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii)
      {
	double edgeLen = norm(n.position() - (*ii).node2().position());
	double K = (*ii).value().springConst;
	double L = (*ii).value().restLen;
	springForce += -K*((n.position() - (*ii).node2().position())/edgeLen)*
	  (edgeLen - L);
      }
    Point totalForce = gravForce + springForce;
    return totalForce;
  }
};

/** Gravity force function object */
struct GravityForce {
  /** Return the gravitational force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //Quiet compiler warning

    //Gravitational force for the node
    Point gravForce = n.value().mass * Point(0,0,-grav);
    return gravForce;
  }
};

/** Spring force function object */
struct SpringForce {
  /** Return the spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //Quiet compiler warning
 
    //Spring force for the node
    Point springForce = Point(0);
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii)
      {
	double edgeLen = (*ii).length();
	double K = (*ii).value().springConst;
	double L = (*ii).value().restLen;
	springForce += -K*((n.position() - (*ii).node2().position())/edgeLen)*
	  (edgeLen - L);
      }
    return springForce;
  }
};

/** Damping force function object */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   *  Damping constant is a member of the DampingForce struct
   *  and is a required user input at instantiation */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //Quiet compiler warning

    Point dampForce = -c_ * n.value().vel;
    return dampForce;
  }

  double c_; //damping constant
  DampingForce(double c) : c_(c)
  {
  }
};

//A struct (function object) for combining two force structs
template <typename Force1, typename Force2>
struct combineTwoForces {

  /** Returns the combined force from @a f1_ and @a f2_
   *  applying to @a n at time @a t */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //(void) t;
    return f1_(n,t) + f2_(n,t);
  }

  //Constructor
  combineTwoForces(Force1 f1, Force2 f2) : f1_(f1), f2_(f2)
  {
  }

private:
  Force1 f1_;
  Force2 f2_;
};

//A struct (function object) for combining two force structs
template <typename Force1, typename Force2, typename Force3>
struct combineThreeForces {

  /** Returns the combined force from @a f1_, @a f2_, and @a f3_
   *  applying to @a n at time @a t */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //(void) t;
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }

  //Constructor
  combineThreeForces(Force1 f1, Force2 f2, Force3 f3)
    : f1_(f1), f2_(f2), f3_(f3)
  {
  }

private:
  Force1 f1_;
  Force2 f2_;
  Force3 f3_;
};

/*Function that combines @a f1 and @a f2 into a single force
 *Returns a force struct (function object)
 */
template<typename Force1, typename Force2>
combineTwoForces<Force1,Force2> make_combined_force(Force1 f1, Force2 f2)
{
  return combineTwoForces<Force1,Force2>(f1,f2);
}

/*Function that combines @a f1, @a f2, and @a f3 into a single force
 *Returns a force struct (function object)
 */
template<typename Force1, typename Force2, typename Force3>
combineThreeForces<Force1,Force2,Force3>
make_combined_force(Force1 f1, Force2 f2, Force3 f3)
{
  return combineThreeForces<Force1,Force2,Force3>(f1,f2,f3);
}

//Constraint function object to fix the given Node at FixedNode1Pos
struct ConstrainNode1 {

  /**Fixes the node of the constraint in @a g to be at the fixed
   * position with zero velocity at time @a t */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) g;
    (void) t;
    //Reset position and set velocity to zero
    n_.position() = FixedNode1Pos;
    n_.value().vel = Point(0);

    //size_type nodeIdx = n_.index();
    //g.node(nodeIdx).position() = FixedNode1Pos;
    //g.node(nodeIdx).value().vel = Point(0);
  }

  //Constructor
  ConstrainNode1(Node n) : n_(n)
  {
  }

private:
  Node n_;
};

//Constraint function object to fix the given Node at FixedNode2Pos
struct ConstrainNode2 {

  /**Fixes the node of the constraint in @a g to be at the fixed
   * position with zero velocity at time @a t */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) g;
    (void) t;
    //Reset position and set velocity to zero
    n_.position() = FixedNode2Pos;
    n_.value().vel = Point(0);
    //size_type nodeIdx = n_.index();
    //g.node(nodeIdx).position() = FixedNode2Pos;
    //g.node(nodeIdx).value().vel = Point(0);
  }

  //Constructor
  ConstrainNode2(Node n) : n_(n)
  {
  }

private:
  Node n_;
};

/* Constraint function object to prevent Nodes from passing through a specified
 * z plane. The coordinate of the plane and its direction are declared in the
 * the struct. Currently, the coordinate of the plane is @a z_coord = -0.75
 * 
 * A node violates the constraint if dot(@a n.position(),(0,0,1)) < @a z_coord
 * A node that violates the constraint has its position reset to the nearest
 * point on the plane and its z velocity component set to zero
 */ 
struct PlaneConstraint {

  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    
    //For all nodes, check if it violates the plane
    Point zPlane = Point(0,0,1);
    Point zNormal = Point(0,0,-1);
    double z_coord = -0.75;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni)
      {
	auto n = *ni;
	if (dot(n.position(),zPlane) < z_coord)
	  {
	    //Find nearest point to plane
	    Point toPlane = n.position() - Point(0,0,z_coord);
	    double dist = dot(toPlane,zNormal);
	    Point proj = n.position() - dist*zNormal;

	    //Update node to correct violation
	    n.position() = proj;
	    n.value().vel[2] = 0;
	  }
      }
  }
  
};

/* Constraint function object to prevent Nodes from passing through a specified
 * sphere. The center of the sphere and its radius are declared in the
 * the struct. Currently, the center of the sphere is @a c = (0.5,0.5,-0.5) and
 * its radius is @a r = 0.15
 * 
 * A node violates the constraint if norm(@a n.position() - c) < @a r
 * A node that violates the constraint has its position reset to the nearest
 * point on the sphere and velocity component normal to the sphere set to zero
 */ 
struct SphereConstraint {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    
    //For all nodes, check if it violates the sphere
    Point center = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni)
      {
	auto n = *ni;
	if (norm(n.position()-center) < radius)
	  {
	    //Find nearest point on sphere
	    Point toCenter = n.position() - center;
	    double dist = norm(toCenter);
	    Point scaled = (radius/dist)*toCenter;
	    Point proj = scaled+center;

	    //Calculate new velocity
	    Point R = toCenter/dist;
	    double dotprod = dot(n.value().vel,R);
	    Point updatedVel = n.value().vel - dotprod*R;

	    n.position() = proj;
	    n.value().vel = updatedVel;
	  }
      }
  }
};

/* Constraint function object to remove Nodes that contacts a specified
 * sphere. The center of the sphere and its radius are declared in the
 * the struct. Currently, the center of the sphere is @a c = (0.5,0.5,-0.5) and
 * its radius is @a r = 0.15
 * 
 * A node violates the constraint if norm(@a n.position() - c) < @a r
 * A node that violates the constraint is removed from the graph @g. All edges
 * incident to the node are removed as well.
 */ 
struct RemoveSphere {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    //Store all nodes that violate the sphere, then remove them
    Point center = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    std::vector<Node> nodesToRemove;

    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni)
      {
	auto n = *ni;
	if (norm(n.position()-center) < radius)
	  {
	    nodesToRemove.push_back(n);
	  }
      }
    for (size_type iter = 0; iter < nodesToRemove.size(); iter++)
      {
	size_type result = g.remove_node(nodesToRemove[iter]);
	(void) result;
      }
  }
};

//A struct (function object) for combining two constraint structs
template <typename Con1, typename Con2>
struct combineTwoConstraints {

  /** Applies both constraints @a c1_ and @a c2_ to @a g at time @a t */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g,t);
    c2_(g,t);
  }

  //Constructor
  combineTwoConstraints(Con1 c1, Con2 c2) : c1_(c1), c2_(c2)
  {
  }

private:  
  Con1 c1_;
  Con2 c2_;
};

//A struct (function object) for combining three constraint structs
template <typename Con1, typename Con2,  typename Con3>
struct combineThreeConstraints {

    /** Applies constraints @a c1_, @a c2_, and @a c3_ to @a g at time @a t */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g,t);
    c2_(g,t);
    c3_(g,t);
  }

  //Constructor
  combineThreeConstraints(Con1 c1, Con2 c2, Con3 c3) : c1_(c1), c2_(c2), c3_(c3)
  {
  }

private:
  Con1 c1_;
  Con2 c2_;
  Con3 c3_;
};

/*Function that combines @a c1 and @a c2 into a single constraint
 *Returns a constraint struct (function object)
 */
template<typename Con1, typename Con2> 
combineTwoConstraints<Con1,Con2> make_combined_constraint(Con1 c1, Con2 c2)
{
  return combineTwoConstraints<Con1,Con2>(c1,c2);
}

/*Function that combines @a c1, @a c2, and @a c3 into a single force
 *Returns a constraint struct (function object)
 */
template<typename Con1, typename Con2, typename Con3>
combineThreeConstraints<Con1,Con2,Con3>
make_combined_constraint(Con1 c1, Con2 c2, Con3 c3)
{
  return combineThreeConstraints<Con1,Con2,Con3>(c1,c2,c3);
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
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  //For all nodes, set initial velocities to zero and mass to 1/N
  //Also find and store the fixed nodes
  double nodeMass = 1.0/graph.num_nodes();
  Node fixedN1; 
  Node fixedN2;
  bool foundFixedN1 = false;
  bool foundFixedN2 = false;
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni)
    {
      (*ni).value().vel = Point(0);
      (*ni).value().mass = nodeMass;
      
      if ((*ni).position() == FixedNode1Pos)
	{
	  fixedN1 = (*ni);
	  foundFixedN1 = true;
	}
      if ((*ni).position() == FixedNode2Pos)
	{
	  fixedN2 = (*ni);
	  foundFixedN2 = true;
	}
    }
  assert(foundFixedN1 && foundFixedN2);
  
  //For all edges, set the spring constant to the given const
  //and the rest length to the initial length of the edge
  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei)
    {
      (*ei).value().springConst = givenSpringConst;
      (*ei).value().restLen = (*ei).length();
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

      //Define the damping constant for the current graph
      double dampConst = 1.0/graph.num_nodes();
      
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
	//Choose Euler step parameters

	//Problems 1 and 2
	//symp_euler_step(graph, t, dt, Problem1Force()); //#1 and #2

 	//No damping, fixed nodes 
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2)));

	//No damping, fixed nodes, plane constraint
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2),
	//					 PlaneConstraint()));

	//No damping, fixed nodes, sphere constraint
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2),
	//					 SphereConstraint()));

	//No damping, plane and sphere constraint, no fixed nodes
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		make_combined_constraint(PlaneConstraint(),
	//					 SphereConstraint()));
	
	//No damping, plane constraint, remove sphere, no fixed nodes
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		make_combined_constraint(PlaneConstraint(),
	//					 RemoveSphere()));
	
	//No damping, sphere constraint only
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),SpringForce()),
	//		SphereConstraint());
	
	//Damping, fixed nodes
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),
	//				    SpringForce(),
	//				    DampingForce(dampConst)),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2)));

	//Damping, fixed nodes, plane constraint
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),
	//				    SpringForce(),
	//				    DampingForce(dampConst)),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2),
	//					 PlaneConstraint()));
	
	//Damping, fixed nodes, sphere constraint
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),
	//				    SpringForce(),
	//				    DampingForce(dampConst)),
	//		make_combined_constraint(ConstrainNode1(fixedN1),
	//					 ConstrainNode2(fixedN2),
	//					 SphereConstraint()));
	
	//Damping, fixed nodes, remove sphere
	symp_euler_step(graph, t, dt,
			make_combined_force(GravityForce(),
					    SpringForce(),
					    DampingForce(dampConst)),
			make_combined_constraint(ConstrainNode1(fixedN1),
						 ConstrainNode2(fixedN2),
						 RemoveSphere()));
	
	//Damping, plane and sphere constraint, no fixed nodes
	//symp_euler_step(graph, t, dt,
	//		make_combined_force(GravityForce(),
	//				    SpringForce(),
	//				    DampingForce(dampConst)),
	//		make_combined_constraint(PlaneConstraint(),
	//					 SphereConstraint()));


	//Clear the viewer's nodes and edges
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
