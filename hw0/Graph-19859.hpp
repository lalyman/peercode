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

  // predeclare internal structs that will hold data for nodes and edges
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
    :  nodes_(), edges_(), num_nodes_(0), num_edges_(0), next_nuid_(0), next_euid_(0)  {
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
      return graph_->nodes_[uid_].position; //-> to access private member
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //return graph_->nodes_[uid_].uid; //an alternate form, useful later?
      return uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check whether uid's are the same
      // Node can access another node's graph_ because Graph is a friend class
      if (this->uid_ == n.index() and this->graph_ == n.graph_) {
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
      // compare uid's, which are ordered by construction
      if (this->uid_ < n.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    ///pointer to graph object (will be used for boolean operators)
    Graph* graph_;
    /// This element's unique identification number
    size_type uid_;

    /**Private Constructor accessed by Graph to construct valid Node objects*/
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
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

    //create a new internal_node - instantiating a struct
    internal_node next_node;

    //set up attributes of new node
    next_node.position = position;
    next_node.uid = next_nuid_;

    nodes_.push_back(next_node);

    ///increment number of nodes and next_nuid_
    ++next_nuid_;
    ++num_nodes_;

    return Node(this, next_nuid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //return true if this is the same as graph_ (both are memory addresses)
    if (n.graph_ == this) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // check that index is within size of graph (i.e. it already exists)
    assert(i < size());
  
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
      return graph_->edges_[uid_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->edges_[uid_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer to Graph container
    Graph* graph_;
    // This element's unique identification number
    size_type uid_;
    /** Private constructor accessed by Graph to construct valid Edge objects */
    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
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
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check that a and b are in the graph
    if (this != a.graph_ or this != b.graph_) {
      return false;
    }
    //loop through edges_ vector, check if a and b match node1 and node2 in each
    // internal edge
    for (size_type i = 0; i < num_edges(); ++i) {
      //if a = node1 and b = node2, return true  -- using overaloaded operators
      if ( edges_[i].node1 == a and edges_[i].node2 == b) {
        return true;
      }
      //if b = node1 and a = node2, return true
      else if (edges_[i].node1 == b and edges_[i].node2 == a) {
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

    //create a new internal_node - instantiate a struct
    internal_edge next_edge;

    //set up attributes of new edge
    next_edge.uid = next_euid_;
    next_edge.node1 = a;
    next_edge.node2 = b;

    //add new edge to vector
    edges_.push_back(next_edge);

    //update counters
    ++next_euid_;
    ++num_edges_;

    return Edge(this, next_euid_ - 1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // use vector clear to remove all data from nodes_ and edges_
    nodes_.clear();
    edges_.clear();
    // reset counters
    num_nodes_ = 0;
    num_edges_ = 0;
    next_nuid_ = 0;
    next_euid_ = 0;
  }

 private:

  // contains data attributes defining a node
  struct internal_node {
    Point position;
    size_type uid;
  };
  // contains data attributes defining an edge
  struct internal_edge {
    Node node1;
    Node node2;
    size_type uid;
  };

  // STL containers for internal_nodes and internal_edges
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

  // attributes for size of graph
  size_type num_nodes_;
  size_type num_edges_;

  //counters for constructor
  size_type next_nuid_; //node counter
  size_type next_euid_; //edge counter

};

#endif // CME212_GRAPH_HPP
