#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <functional>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // Internal type for nodes
  struct internal_node;
  
  // Internal type for edges
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_vec[uid].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) && (uid == n.uid));
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (*this == n) {
        return false;
      }
      else if (uid == n.uid) {          // Same index, different graph
        std::less<Graph*> less;
        return less(graph_,n.graph_);    // std::less gives total ordering
      }
      else {
        return (uid < n.uid);           // Else, return based on index
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // Unique node id
    size_type uid;

    /** Private Node constructor (gives access to outer graph and 
     * initializes node's id)
     */
    Node(const Graph* graph, size_type i)
        : graph_(const_cast<Graph*>(graph)), uid(i) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_vec.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    internal_node n_internal(position);
    nodes_vec.push_back(n_internal);

    return node(nodes_vec.size() - 1); // Use node function to return most recent node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= nodes_vec.size()) {
      return false;
    }
    return (node(n.index()) == n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   *
   * If node index argument is greater than number of nodes in graph 
   * (minus 1), then an invalid node is returned
   */
  Node node(size_type i) const {
    if (i < num_nodes()) {
      Node n(this,i);
      return n;
    }
    Node n;
    return n;
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(graph_->edges_vec[uid].n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(graph_->edges_vec[uid].n2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((node1() == e.node1() && node2() == e.node2()) ||
          (node1() == e.node2() && node2() == e.node1())) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (*this == e) {
        return false;
      }
      else if (uid == e.uid) {          // Same index, different graph
        std::less<Graph*> less;
        return less(graph_,e.graph_);    // std::less gives total ordering
      }
      else {
        return (uid < e.uid);           // Else, return based on index
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // Unique edge id
    size_type uid;

    /** Private Edge constructor (gives access to outer graph and 
     * initializes edge's id)
     */
    Edge(const Graph* graph, size_type i)
        : graph_(const_cast<Graph*>(graph)), uid(i) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   *
   * If edge index argument is greater than number of edges in graph 
   * (minus 1), then an invalid edge is returned
   */
  Edge edge(size_type i) const {
    if (i < num_edges()) {
      Edge e(this, i);
      return e;
    }
    Edge e;
    return e; 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (has_edge_helper(a, b) == -1) {
      return false;
    }
    return true;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(has_node(a) && has_node(b));
    if (!has_edge(a,b)) {
      internal_edge e(a,b);
      edges_vec.push_back(e);
      // Return the most recent edge added
      return edge(edges_vec.size() - 1);
    }
    else {
      return edge((size_type)has_edge_helper(a,b));
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_vec.clear();
    edges_vec.clear();
  }

 private:
  
  // Internal node type
  struct internal_node {
    Point position;
    // Constructor for internal_node when passed a point object
    internal_node(const Point& p) : position(p) {
    }
  };

  // Internal edge type
  struct internal_edge {
    size_type n1;
    size_type n2;
    // Constructor for internal_edge when passed two node objects
    internal_edge(const Node& node1, const Node& node2)
        : n1(node1.index()), n2(node2.index()) {
    }
  };

  // Internal data for graph's nodes and edges, stored in vectors
  std::vector<internal_node> nodes_vec;
  std::vector<internal_edge> edges_vec;

  //  Helper Functions
  // ******************

  /** Helper function that checks if a graph has an edge and returns
   *  the index of that edge if true and -1 otherwise
   */
  int has_edge_helper(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    for (size_type i = 0; i < edges_vec.size(); i++) {
      if ((edge(i).node1() == a && edge(i).node2() == b) ||
        (edge(i).node1() == b && edge(i).node2() == a)) {
        return (int)i;
      }
    }
    return -1; 
  }
};

#endif // CME212_GRAPH_HPP
