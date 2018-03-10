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
  NodeData() : vel(0), mass(1) {}
};

/** Custon structure of data to store with Edges*/
struct EdgeData {
  double K;       //< Node velocity
  double length;     //< Node mass
  EdgeData() : K(100), length(1) {}
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Predicate for distance of a Node to Point */
struct CirclePredicate {
    template <typename NODE>
    bool operator()(const NODE&n) const {
        return norm_2(n.position() - Point (0.5, 0.5, -0.5)) < 0.15;
    }
};

/** Predicate for height of a Node to z = -0.75. */
struct LevelPredicate {
    template <typename NODE>
    bool operator()(const NODE&n) const {
        return dot(n.position(), Point (0,0,1)) < -0.75;
    }
};


/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>>
{
 public:
  // Get all of the iterator traits and make them our own
  using value_type        = typename std::iterator_traits<It>::value_type;
  using pointer           = typename std::iterator_traits<It>::pointer;
  using reference         = typename std::iterator_traits<It>::reference;
  using difference_type   = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return a node this iterator is pointing to */
    value_type operator*() const {
      return  *it_;
    }
    /** increment the iterator to point to next qualified node */
    filter_iterator& operator++(){
      ++it_;
      while(it_!=end_ && !p_(*it_)){
        ++it_;
      }
      return *this;
    }

    /** Test if two filter_iterators are the same
     * the same filter_iterators point to the same node in the same graph
     */
    bool operator==(const filter_iterator& itr) const{
      return (it_ == itr.it_);
    }
 private:
  Pred p_;
  It it_;
  It end_;
};


/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}



/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and graph constraint.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     cons  Function object defining the constraint to the graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports double mass, double vel
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a cons(g), where g is the graph.
 *           @a cons returns iterator of node that violate constraints.
 */

template <typename G, typename F, typename C>
double symp_euler_step_con1(G& g, double t, double dt, F force, C cons) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  auto filter_itr_begin = ++make_filtered(g.node_begin(), g.node_end(), cons);
  auto filter_itr_end = make_filtered(g.node_end(), g.node_end(), cons);
  for (auto it = filter_itr_begin; it != filter_itr_end; ++it){
      (*it).position().z = -0.75;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  for (auto it = filter_itr_begin; it != filter_itr_end; ++it){
      (*it).value().vel.z = 0;
  }

  return t + dt;
}

/** See description of symp_euler_step_con2*/
template <typename G, typename F, typename C>
double symp_euler_step_con2(G& g, double t, double dt, F force, C cons) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  auto filter_itr_begin = ++make_filtered(g.node_begin(), g.node_end(), cons);
  auto filter_itr_end = make_filtered(g.node_end(), g.node_end(), cons);
  for (auto it = filter_itr_begin; it != filter_itr_end; ++it){
      Point r = (*it).position() - Point(0.5, 0.5, -0.5);
      (*it).position() = r/norm_2(r)*0.15 + Point(0.5, 0.5, -0.5);
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  for (auto it = filter_itr_begin; it != filter_itr_end; ++it){
      Point R = ((*it).position() - Point(0.5, 0.5, -0.5))/norm_2((*it).position() - Point(0.5, 0.5, -0.5));
      (*it).value().vel = (*it).value().vel - dot((*it).value().vel, R)*R;
  }

  return t + dt;
}


template <typename G, typename F, typename C>
double symp_euler_step_destroy(G& g, double t, double dt, F force, C cons) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  auto filter_itr_begin = ++make_filtered(g.node_begin(), g.node_end(), cons);
  auto filter_itr_end = make_filtered(g.node_end(), g.node_end(), cons);
  std::vector<Node> n;
  for (auto it = filter_itr_begin; it != filter_itr_end; ++it){
    n.push_back(*it);
  }
  for (auto it = n.begin(); it != n.end(); ++it){
      g.remove_node(*it);
  }

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
    Point force (0,0,0);
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
        return force;
    }else{
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
            double K = (*it).value().K;
            Point dis = (*it).node1().position() - (*it).node2().position();
            force -= K*(dis)/norm_2(dis)*(norm_2(dis) - (*it).value().length);
        }
        Point gravity (0,0,-grav);
        force += n.value().mass * gravity;
        return force;
    }
  }
};

/** A function object that returns the force applying to NODE n at time t */
template <typename F1, typename F2, typename F3 = void>
class combined_force{
public:
  // Constructor
  combined_force(const F1 f1, const F2 f2, const F3 f3)
      : f1_(f1), f2_(f2), f3_(f3) {
  }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
        return Point (0,0,0);
    }else{
      return  f1_(n ,t) + f2_(n, t) + f3_(n, t);
    }
  }
private:
  F1 f1_;
  F2 f2_;
  F3 f3_;
};

template <typename F1, typename F2>
class combined_force<F1, F2, void>{
public:
  // Constructor
  combined_force(const F1 f1, const F2 f2)
      : f1_(f1), f2_(f2){
  }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
        return Point (0,0,0);
    }else{
      return  f1_(n ,t) + f2_(n, t);
    }
  }
private:
  F1 f1_;
  F2 f2_;
};

/** Helper function for constructing combined_force. This deduces the type of
 * the forces so the user doesn't have to write it.
 */
template <typename F1, typename F2, typename F3>
combined_force<F1, F2, F3> make_combined_force(const F1 f1, const F2 f2, const F3 f3) {
  return combined_force<F1, F2, F3>(f1, f2, f3);
}


template <typename F1, typename F2>
combined_force<F1, F2> make_combined_force(const F1 f1, const F2 f2) {
  return combined_force<F1, F2>(f1, f2);
}

/** Force function object for gravity*/
struct GravityForce {
    template <typename NODE>
    Point operator()(NODE n, double t){
      Point gravity (0,0, -grav);
      Point force = n.value().mass * gravity;
      return force;
    }
};

/** Force function object for mass spring. */
struct MassSpringForce {
    template <typename NODE>
    Point operator()(NODE n, double t){
      Point force (0,0,0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
          double K = (*it).value().K;
          Point dis = (*it).node1().position() - (*it).node2().position();
          force -= K*(dis)/norm_2(dis)*(norm_2(dis) - (*it).value().length);
      }
      return force;
    }
};

/** Force function object for damping force. */
struct DampingForce {
    template <typename NODE>
    Point operator()(NODE n, double t){
        Point force = - n.value().mass * n.value().vel;
        return force;
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double mass = 1./graph.num_nodes();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      (*it).value().mass = mass;
  }
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
      (*it).value().K = 100;
      (*it).value().length = norm_2((*it).node1().position() - (*it).node2().position());
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
      double dt = 0.0001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        viewer.clear();
        node_map.clear();


        /** change function call for different constraints. */
        symp_euler_step_destroy(graph, t, dt, Problem1Force(), CirclePredicate());

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
