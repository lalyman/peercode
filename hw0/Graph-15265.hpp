#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
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

  /** Construct an empty graph.
   * Uses default constructors for each of the member variables. List
   * initialization is used.
   */
  Graph() : nodes_(), edges_(), edgeMap_() {}

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
     *
     * Here, pointer graph_ is used as a test of validity. Value will be
     * nullptr, unless the Node is initialized by class Graph, and graph_
     * pointer points toward the comprising Graph object.
     */
    Node() : nodeID_(0), graph_(nullptr) {}
    /** Return this node's position. */
    const Point& position() const {
      assert(graph_);   // Check for validity
      return graph_->nodes_[nodeID_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(graph_);   // Check for validity
      return nodeID_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(graph_ && n.getGraph());   // Check for validity of n and this
      (void) n;          // Quiet compiler warning
      return (graph_ == n.getGraph()) && (nodeID_ == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning. Here, it is
     * implemented by comparing ID numbers, which are unique to Node objects
     * in a graph. Comparison is only valid for Nodes that are part of the same
     * graph.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      assert(graph_ && n.getGraph());   // Check for validity of n and this
      assert(graph_ == n.getGraph());   // Check both nodes in same graph
      (void) n;           // Quiet compiler warning
      return nodeID_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Node ID Number
    size_type nodeID_;

    // Pointer to graph
    Graph* graph_;

    // Return pointer to graph; useful in class Graph member functions
    Graph* getGraph() const {
      return graph_;
    }

    // Constructor for use by Graph object
    Node(size_type ID, Graph* graph) : nodeID_(ID), graph_(graph) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
    nodes_.push_back(position);
    Node returnNode(nodes_.size() - 1, this);
    (void) position;      // Quiet compiler warning
    return returnNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    return this == (const Graph*)n.getGraph();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert((i < size()) && i >= 0);
    Node returnNode(i, (Graph*) this);
    (void) i;             // Quiet compiler warning
    return returnNode;
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
    /** Construct an invalid Edge. graph_ pointer is used as a test of
     * validity. Value will be nullptr unless Edge is constructed by class
     * Graph, and graph_ points to the comprising Graph object.
     */
    Edge() : edgeID_(0), graph_(nullptr) {}

    /** Return a node of this Edge */
    Node node1() const {
      assert(graph_);   // Check for validity
      return graph_->node(graph_->edges_[edgeID_][0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(graph_);   // Check for validity
      return graph_->node(graph_->edges_[edgeID_][1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(graph_ && e.getGraph());   // Check for validity
      assert(graph_ == e.getGraph());   // Check for same graph
      auto pair = graph_->edges_[edgeID_];
      auto ePair = e.getGraph()->edges_[e.getID()];
      bool match = pair[0] == ePair[0] && pair[1] == ePair[1];
      bool crossMatch = pair[1] == ePair[0] && pair[0] == ePair[1];
      (void) e;           // Quiet compiler warning
      return match || crossMatch;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning. Here, we compare
     * using edgeID numbers. Comparison is not valid unless edges belong to the
     * same graph.
     */
    bool operator<(const Edge& e) const {
      assert(graph_ && e.getGraph());   // Check for validity
      assert(graph_ == e.getGraph());   // Check for same graph
      (void) e;           // Quiet compiler warning
      return edgeID_ < e.getID();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Variable for edge ID
    size_type edgeID_;

    // Pointer to graph
    Graph* graph_;

    // Return pointer to graph, for use by class Graph
    Graph* getGraph() const {
      return graph_;
    }

    // Return edge ID
    size_type getID() const {
      return edgeID_;
    }

    // Constructor for use by Graph object
    Edge(size_type ID, Graph* graph) : edgeID_(ID), graph_(graph) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1).
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1).
   */
  Edge edge(size_type i) const {
    assert((0 <= i) && (i < num_edges()));
    (void) i;             // Quiet compiler warning
    return Edge(i, (Graph*) this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1) by using an unordered map keyed by unique
   * string identifiers for pairs of nodes (Explained in detail below).
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    std::string search = repString(a, b);
    auto it = edgeMap_.find(search);
    (void) a; (void) b;   // Quiet compiler warning
    return it != edgeMap_.end();
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1) using the has_edge function written above. The
   * details of this scheme are explained there, and further below.
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(!(a == b));  // Check distinctness
    assert(a.getGraph() && b.getGraph());       // Check validity
    assert(a.getGraph() == b.getGraph());       // Check on same graph
    if (has_edge(a, b)) {
      std::string search = repString(a,b);
      return edge(edgeMap_[search]);
    } else {
      std::string search = repString(a, b);
      std::vector<size_type> entry {a.index(), b.index()};
      edges_.push_back(entry);
      edgeMap_[search] = num_edges()-1;
      (void) a, (void) b;   // Quiet compiler warning
      return edge(num_edges()-1);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    edgeMap_.clear();
  }

 private:
  /** Data container for nodes. Each node is represented as a Point object. The
   * ID number of the node is used to index into the vector. This allows for
   * O(1) return of Point based on NodeID. Vector allows us to add Nodes in
   * O(1) amortized.
   */
  std::vector<Point> nodes_;

  /** Data container for edges. Each edge is represented as a vector with two
   * size_type elements. Each size_type corresponds to the ID number of a node.
   * The ID number of the vector is used to index into the (outer) vector. This
   * ID number can be found O(1), as explained below.
   */
  std::vector<std::vector<size_type>> edges_;

  /** Data container for O(1) edge lookups. Unordered_map key is string
   * representation of Edge. This representation is unique, see below.
   * Unordered_map value is ID number for edge. Allows for constant time lookup
   * of edge ID's, based on just two nodes as input. This is crucial for
   * has_edge() method and add_edge() method.
   */
  std::unordered_map<std::string, size_type> edgeMap_;

  /** Return a unique string representation of node pair, for keys. String
   * representation is simply the two node ID's, separated by a '_' character.
   * The ordering is defined to be lower ID followed by higher ID. Since edges
   * are not repeated, this is sufficient to create a unique mapping from edges
   * to string representations. This allows for using an unordered_map
   * container (see above), which allows us to have O(1) lookup of whether or
   * not an edge connects two nodes, in addition to an O(1) lookup of the edge
   * ID number if such an edge exists. Allows for O(1) has_edge() and
   * add_edge().
   */
  std::string repString(const Node& a, const Node& b) const {
    const Node& low = a < b ? a : b;
    const Node& high = a < b ? b : a;
    std::stringstream search;
    search << std::to_string(low.index()) << "_"
      << std::to_string(high.index());
    return search.str();
  }
};

#endif // CME212_GRAPH_HPP
