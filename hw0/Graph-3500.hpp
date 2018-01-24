#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
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

  // Predeclare the internal structs.
  struct internal_node;
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
  Graph()
    : internal_nodes_(), internal_edges_(), size_(0), num_edges_(0) {
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

    /** Return this node's position.
     * @pre This Node is valid and belongs to a Graph.
     */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * @pre This Node is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().uid;
    }

    /** Return pointer to the graph to which this Node belongs. */
    Graph* graph() const {
      return graph_;
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /** Check if (1) graph to which n belongs is same as
       * this nodes's graph and (2) the two nodes have the same
       * index/uid.
       */
      Graph* n_graph = n.graph();
      if ((graph_ == n_graph) && (uid_ == n.index())) {
        return true;
      }
      return false;
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
      Graph* n_graph = n.graph();
      if (graph_ < n_graph) {
        // This node's graph is less than n's graph (by pointer comparison).
        return true;
      }

      if (graph_ > n_graph) {
        // This node's graph is larger than n's graph (by pointer comparison).
        return false;
      }

      // Otherwise the graphs are the same (by pointer comparison) so compare
      // the node index/uid.
      return (uid_ < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to Graph to which this node belongs.
    Graph* graph_;

    // This node's unique ID. The unique ID is the same as the node's index.
    // Both correspond to the order in which nodes where added to the graph.

    // Important note: To support node removal in the future, the uid and
    // index will have to be differentiated.
    size_type uid_;

    // Private Constructor.
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    // Helper method to return appropriate internal_node.
    internal_node& fetch() const {
      return graph_->internal_nodes_[uid_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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

    // Create new internal_node object to add.
    internal_node new_node;
    new_node.position = position;
    new_node.uid = size_;

    // Add to list of nodes.
    internal_nodes_.push_back(new_node);

    // Increment size.
    ++size_;

    // Return a Node that points to newly added node.
    return Node(this, size_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type idx = n.index();

    // If the index of the n is larger than this graph's
    // size (# nodes) then n cannot belong to the graph.
    if (size_ < idx) {
      return false;
    }

    // Otherwise, check if n and the node in this graph with
    // the same index as n are equal.
    Node m = node(idx);
    return (n == m);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
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
      return fetch().node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return fetch().node2;
    }

    /** Return index of edge, a number in the range [0, num_edges).
     * @pre This Edge is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().edge_uid;
    }

    /** Return pointer to graph to which this Edge belongs */
    Graph* graph() const {
      return graph_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if (1) graph to which this edge belongs is same as
      // current graph and (2) the two edges have the same uid.

      // Duplicate edges are checked for when attempting to add edges to a
      // graph, so two edges with the same index/uid that have ended up in
      // the same graph must be equal (represent the same undirected edge
      // between two nodes.)
      Graph* e_graph = e.graph();
      if ((graph_ == e_graph) && (edge_uid_ == e.index())) {
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
      Graph* e_graph = e.graph();
      if (graph_ < e_graph) {
        // This edges's graph is less than e's graph (by pointer comparison).
        return true;
      }

      if (graph_ > e_graph) {
        // This edge's graph is larger than e's graph (by pointer comparison).
        return false;
      }

      // Otherwise the graphs are equal so compare the edge index/uid.
      return (edge_uid_ < e.index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph to which this edge belongs.
    Graph* graph_;

    // This edge's unique ID. The unique ID is the same as the edge's
    // index. Both correspond to the order in which edges were added
    // to the graph.

    // Important note: To support edge removal in the future, the uid and
    // index will have to be differentiated.
    size_type edge_uid_;

    // Private constructor.
    Edge(const Graph* graph, size_type edge_uid)
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }

    // Helper method to return appropriate internal_edge.
    internal_edge& fetch() const {
      return graph_->internal_edges_[edge_uid_];
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < internal_edges_.size(); ++i) {
      internal_edge current_edge = internal_edges_[i];
      Node current_edge_node1 = current_edge.node1;
      Node current_edge_node2 = current_edge.node2;

      // Edges are undirected, so to check if edge already exists
      // we compare against (a,b) and (b,a).
      if ((a == current_edge_node1) && (b == current_edge_node2)) {
        return true;
      }

      if ((b == current_edge_node1) && (a == current_edge_node2)) {
        return true;
      }
    }

    return false;
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
    // Add edge only if edge does not already exists.
    if (!(has_edge(a, b))) {
      // Create new internal_edge object to add.
      internal_edge new_edge;
      new_edge.node1 = a;
      new_edge.node2 = b;
      new_edge.edge_uid = num_edges_;

      // Add to list of edges.
      internal_edges_.push_back(new_edge);

      // Increment number of edges.
      num_edges_++;

      // Return an Edge that points to newly added edge.
      return Edge(this, num_edges_ - 1);
    }
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_nodes_.clear();
    internal_edges_.clear();
    size_ = 0;
    num_edges_ = 0;
  }

 private:

  // Graph class's internals helper functions, data members, and so forth.

  struct internal_node {
    // The position of the node.
    Point position;
    // The uid/index of the node (based on order in which nodes were added).
    size_type uid;
  };

  struct internal_edge {
    // The pair of nodes connected by this edge.
    Node node1;
    Node node2;
    // The uid/index of the edge (based on order in which edges were added).
    size_type edge_uid;
  };

  std::vector<internal_node> internal_nodes_;  // Container for nodes in graph.
  std::vector<internal_edge> internal_edges_;  // Container for edges in graph.
  size_type size_;  // Graph size a.k.a. number of nodes in the graph.
  size_type num_edges_;  // Number of edges in the graph.
};

#endif // CME212_GRAPH_HPP
