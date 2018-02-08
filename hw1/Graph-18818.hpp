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
template <typename V>
class Graph {
 private:

  // Predeclare the internal struct
  struct internal_node;
  struct internal_edge;
  // Declare Vector
  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
  class Node : private totally_ordered<Node> {
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
      return graph_->nodes[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Set this node's value. */
    node_value_type& value() {
      return graph_->nodes[uid_].value;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_->nodes[uid_].value;
    }

    /** Return this node's degree in range [0, graph_size]. */
    size_type degree() const {
      return graph_->nodes[uid_].edges.size();
    }

    /** Return a begining iterator over adjacent nodes to this node. */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    /** Return an end iterator over adjacent nodes to this node. */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, uid_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (uid_ == n.index());
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
      if (graph_ == n.graph_) {
        if (uid_ == n.index()) {
          return false;  // same nodes
        } else {
          return uid_ < n.index(); // same graph different index
        }
      } else {
        return graph_ < n.graph_;  // different graphs
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph
    graph_type* graph_;
    // This element's unique identification number
    size_type uid_;
    /** Private Constructor */
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Add internal node
    internal_node n;
    n.point = position;
    n.value = value;
    nodes.push_back(n);
    // Return node
    return Node(this, nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() < nodes.size()) {
      return nodes[n.index()].point == n.position() && 
             nodes[n.index()].value == n.value();
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.node1() == this->node1()) && (e.node2() == this->node2())) || 
             ((e.node1() == this->node2()) && (e.node2() == this->node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return n1_ < e.n1_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph
    graph_type* graph_;
    // First nodes unique identification number
    size_type n1_;
    // Second nodes unique identification number
    size_type n2_;
    /** Private Constructor */
    Edge(const graph_type* graph, size_type n1, size_type n2)
        : graph_(const_cast<graph_type*>(graph)), n1_(n1), n2_(n2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type n1, size_type n2) const {
    return Edge(this, n1, n2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if graph contains node @a, @b
    if (this->has_node(a) && this->has_node(b)) {
      // Iterate through edges of node @a
      for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
        // Get edge.
        auto edge = *ei;
        // Check if node @b is a neighbor of node @a
        if (edge.node2() == b) {
          return true;
        }
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
    if (!this->has_edge(a, b)) {
      // Add internal edge
      internal_edge e;
      e.index1 = a.index();
      e.index2 = b.index();
      edges.push_back(e);
      // Update Node's edges
      nodes[a.index()].edges.push_back(b.index());
      nodes[b.index()].edges.push_back(a.index());
    }
    // Return Edge
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear(); 
    edges.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Return current node of this NodeIterator. */
    Node operator*() const {
      return graph_->node(pos_);
    }

    /** Increment position of this NodeIterator. */
    NodeIterator& operator++() {
      pos_++;
      return *this;
    }

    /** Check if this NodeIterator is equal to @x NodeIterator. 
    *    @param x NodeIterator to compare to.
    */
    bool operator==(const NodeIterator& x) const {
      return graph_ == x.graph_ && pos_ == x.pos_;
    }

   private:
    friend class Graph;
    // Pointer back to the Graph
    graph_type* graph_;
    // Current node's unique identification number
    size_type pos_;
    /** Private Constructor */
    NodeIterator(const graph_type* graph, size_type pos)
        : graph_(const_cast<graph_type*>(graph)), pos_(pos) {
    }
  };

  /** Return a begining iterator over nodes in this graph. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an end iterator over nodes in this graph. */
  node_iterator node_end() const {
    return NodeIterator(this, nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Return current edge of this IncidentIterator. */
    Edge operator*() const {
      size_type edge_index = graph_->nodes[uid_].edges[pos_];
      return graph_->edge(uid_, edge_index);
    }

    /** Increment position of this IncidentIterator. */
    IncidentIterator& operator++() {
      pos_++;
      return *this;
    }

    /** Check if this IncidentIterator is equal to @x IncidentIterator.
    *    @param x IncidentIterator to compare to.
    */
    bool operator==(const IncidentIterator& x) const {
      return graph_ == x.graph_ && uid_ == x.uid_ && pos_ == x.pos_;
    }

   private:
    friend class Graph;
    // Pointer back to the Graph
    graph_type* graph_;
    // Nodes unique identification number
    size_type uid_;
    // Current neighbor
    size_type pos_;
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, size_type uid, size_type pos)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid), pos_(pos) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Return current edge of this EdgeIterator. */
    Edge operator*() const {
      internal_edge e = graph_->edges[pos_];
      return graph_->edge(e.index1, e.index2);
    }

    /** Increment position of this EdgeIterator. */
    EdgeIterator& operator++() {
      pos_++;
      return *this;
    }

    /** Check if this EdgeIterator is equal to @x EdgeIterator.
    *    @param x EdgeIterator to compare to.
    */
    bool operator==(const EdgeIterator& x) const {
      return graph_ == x.graph_ && pos_ == x.pos_;
    }

   private:
    friend class Graph;
    // Pointer back to the Graph
    graph_type* graph_;
    // Current nodes unique identification number
    size_type pos_;
    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type pos)
        : graph_(const_cast<graph_type*>(graph)), pos_(pos) {
    }
  };

  /** Return a begining iterator over edges in this graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return an end iterator over edges in this graph. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges.size());
  }

 private:

  // internal_node attributes
  struct internal_node {
    Point point;   // position of node
    node_value_type value;  // value of node
    std::vector<size_type> edges;  // neighboring nodes index
  };
  // internal_edge attributes
  struct internal_edge {
    size_type index1;   // index node 1
    size_type index2;   // index node 2
  };
};

#endif // CME212_GRAPH_HPP
