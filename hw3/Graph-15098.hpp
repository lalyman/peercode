#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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

  /** Type of node value. */
  using node_value_type = V;

  /** Type of edge value. */
  using edge_value_type = E;

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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      if (uid_ < graph_->nodes_.size())
        return graph_->nodes_[uid_];
      assert(false);
    }

    Point& position() {
      if (uid_ < graph_->nodes_.size())
        return graph_->nodes_[uid_];
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return node value. 
    * @return  A value of node_value_type Type. 
    *
    * @pre     uid_ < graph_->size(). 
    *
    * Complexity: O(1).
    */
    node_value_type& value() {
      if (uid_ < graph_->size())
        return graph_->values_[uid_];
      assert(false);
    }

    /** Return node value. 
    * @return  A value of node_value_type Type. 
    *
    * @pre     uid_ < graph_->size(). 
    *
    * Complexity: O(1).
    */
    const node_value_type& value() const {
      if (uid_ < graph_->size())
        return graph_->values_[uid_];
      assert(false);
    }

    /** Return node degree. 
    * @return  The degree of the node of size_type Type. 
    *
    * @pre     uid_ < graph_->size(). 
    *
    * Complexity: O(1).
    */
    size_type degree() const {
      return graph_->adj_edges_[uid_].size();
    }

    /** Return an iterator for the node's first incident edge. 
    * @return  An iterator of IncidentIterator Type. 
    *
    * @pre     uid_ < graph_->size(). 
    * @post    *result is the first incident edge of the node or 
    *          result refers to the pass-the-end element if the node
    *          has no incident edges. 
    *
    * Complexity: O(1).
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, uid_, size_type(0));
    }

    /** Return an iterator for the past-the-end element of the node's
    *   incident edges. 
    * @return  An iterator of IncidentIterator Type. 
    *
    * @pre     uid_ < graph_->size().
    * @post    result refers to the pass-the-end element. 
    *
    * Complexity: O(1).
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == n.graph_) && (uid_ == n.uid_));
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
      // HW0: YOUR CODE HERE
      // return ((graph_ == n.graph_) && (uid_ < n.uid_));
      return ((graph_ < n.graph_) || ((graph_ == n.graph_) && (uid_ < n.uid_)));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(position);
    values_.push_back(value);
    adj_edges_[num_nodes()-1];
    return Node(this, num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((this == n.graph_) && (n.uid_ < num_nodes()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_nodes())
      return Node(this, i);
    assert(false);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (uid_ < graph_->edges_.size())
        return graph_->node(uid_node_1_);
      assert(false);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (uid_ < graph_->edges_.size())
        return graph_->node(uid_node_2_);
      assert(false);
    }

    double length() const {
      return norm_2(node1().position()-node2().position());
    }

    edge_value_type& value() {
      if (uid_ < graph_->num_edges())
        return graph_->edge_values_[uid_];
      assert(false);
    }

    const edge_value_type& value() const {
      if (uid_ < graph_->num_edges())
        return graph_->edge_values_[uid_];
      assert(false);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((node1() == e.node1()) && (node2() == e.node2())) || ((node1() == e.node2()) && (node2() == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // return ((graph_ == e.graph_) && (uid_ < e.uid_));
      return ((graph_ < e.graph_) || ((graph_ == e.graph_) && (uid_ < e.uid_)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    size_type uid_node_1_;
    size_type uid_node_2_;

    Edge(const Graph* graph, size_type uid, size_type uid_node_1, size_type uid_node_2)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), uid_node_1_(uid_node_1), uid_node_2_(uid_node_2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges())
      return Edge(this, i, edges_[i][0], edges_[i][1]);
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b)) {
      Edge e;
      std::vector<size_type> adj_a = adj_edges_.at(a.index());
      for (size_type i=0; i<adj_a.size(); i++) {
        e = edge(adj_a[i]);
        if (((a == e.node1()) && (b == e.node2())) || ((a == e.node2()) && (b == e.node1())))
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b) && !(a == b)) {
      Edge e;
      std::vector<size_type> adj_a = adj_edges_.at(a.index());
      for (size_type i=0; i<adj_a.size(); i++) {
        e = edge(adj_a[i]);
        if (((a == e.node1()) && (b == e.node2())) || ((a == e.node2()) && (b == e.node1())))
          return Edge(this, adj_a[i], a.index(), b.index()); // double check later
      }
      std::vector<size_type> new_edge;
      new_edge.push_back(a.index());
      new_edge.push_back(b.index());
      edges_.push_back(new_edge);
      edge_values_.push_back(value);
      adj_edges_[a.index()].push_back(num_edges()-1);
      adj_edges_[b.index()].push_back(num_edges()-1);
      return edge(num_edges()-1);
    }
    assert(false);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    values_.clear();
    edges_.clear();
    edge_values_.clear();
    adj_edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return a node that the iterator refers to. 
    * @return  A node. 
    *
    * @pre     uid_ < graph_->size().
    *
    * Complexity: O(1).
    */
    Node operator*() const {
      return graph_->node(uid_);
    }

    /** Return the next iterator for the node's incident edges. 
    * @return  An iterator of IncidentIterator Type. 
    *
    * @pre     uid_ < graph_->size().
    *
    * Complexity: O(1).
    */
    NodeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Return a boolean indicating whether two iterators are equal. 
    * @para[in] ni    A node iterator
    * @return         true if @a ni is the same as this iterator,
    *                 false otherwise.  
    *
    * @post           result is true if @a ni is the same as this iterator,
    *                 otherwise result is false. 
    *
    * Complexity: O(1).
    */
    bool operator==(const NodeIterator& ni) const {
      return ((graph_ == ni.graph_) && (uid_ == ni.uid_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;

    NodeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return an iterator for the first node in the graph. 
  * @return  An iterator of NodeIterator Type. 
  *
  * @post    result refers to first node in the graph. 
  *
  * Complexity: O(1).
  */
  NodeIterator node_begin() const {
    return NodeIterator(this, size_type(0));
  }

  /** Return an iterator for the past-the-end element of the nodes in the graph. 
  * @return  An iterator of NodeIterator Type. 
  *
  * @post    result refers to the pass-the-end element. 
  *
  * Complexity: O(1).
  */
  NodeIterator node_end() const {
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return an edge that the iterator refers to. 
    * @return  An edge. 
    *
    * @pre     uid_ < graph_->size().
    * @pre     uid_adj_idx_ < graph_->node(uid_).degree(). 
    * @post    The first node of result is graph_->node(uid_). 
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      size_type edge_idx = graph_->adj_edges_[uid_][uid_adj_idx_];
      if (graph_->edges_[edge_idx][0] == uid_)
        return Edge(graph_, edge_idx, uid_, graph_->edges_[edge_idx][1]);
      else
        return Edge(graph_, edge_idx, uid_, graph_->edges_[edge_idx][0]);
    }

    /** Return the next iterator for the incident edge iterator. 
    * @return  An iterator of IncidentIterator Type. 
    *
    * @pre     uid_ < graph_->size().
    * @pre     uid_adj_idx_ < graph_->node(uid_).degree(). 
    *
    * Complexity: O(1).
    */
    IncidentIterator& operator++() {
      uid_adj_idx_++;
      return *this;
    }

    /** Return a boolean indicating whether two iterators are equal. 
    * @para[in] iit   An incident edge iterator
    * @return         true if @a iit is the same as this iterator,
    *                 false otherwise. 
    *
    * @post           result is true if @a iit is the same as this iterator,
    *                 otherwise result is false. 
    *
    * Complexity: O(1).
    */
    bool operator==(const IncidentIterator& iit) const {
      return ((graph_ == iit.graph_) && (uid_ == iit.uid_) && (uid_adj_idx_ == iit.uid_adj_idx_));
    }


   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type uid_; // node index
    size_type uid_adj_idx_; // adj edge index

    IncidentIterator(const Graph* graph, size_type uid, size_type uid_adj_idx)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), uid_adj_idx_(uid_adj_idx) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return an edge that the iterator refers to. 
    * @return  An Edge. 
    *
    * @pre     uid_ < graph_->num_edges().
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      return graph_->edge(uid_);
    }

    /** Return the next iterator for the edges in the graph. 
    * @return  An iterator of EdgeIterator Type. 
    *
    * @pre     uid_ < graph_->num_edges().
    *
    * Complexity: O(1).
    */
    EdgeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Return a boolean indicating whether two iterators are equal. 
    * @para[in] eit   An edge iterator
    * @return         true if @a eit is the same as this iterator,
    *                 false otherwise.  
    *
    * @post           result is true if @a eit is the same as this iterator,
    *                 otherwise result is false. 
    *
    * Complexity: O(1).
    */
    bool operator==(const EdgeIterator& eit) const {
      return ((graph_ == eit.graph_) && (uid_ == eit.uid_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;

    EdgeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return an iterator for the first edge in the graph. 
  * @return  An iterator of EdgeIterator Type. 
  *
  * @post    result refers to first edge in the graph. 
  *
  * Complexity: O(1).
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, size_type(0));
  }

  /** Return an iterator for the past-the-end element of the edges in the graph. 
  * @return  An iterator of EdgeIterator Type. 
  *
  * @post    result refers to the pass-the-end element. 
  *
  * Complexity: O(1).
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /** Remove an edge from a graph.
  * @param[in] a,b   Two nodes of the edge.
  * @return          The new index of the edge with old index old num_edges() - 1.
  *
  * @pre    @a a and @a b are two valid nodes in the graph.
  * @post   The edge is removed from the adjacent edges vector for both @a a and @a b.
  * @post   The edge is removed from the vector storing all the edges.
  * @post   The corresponding edge value is removed from the vector storing all the
  *             edge values.
  * @post   If the edge with old index old num_edges() - 1 is not removed, result is
  *             the new index of this edge. Otherwise, result is the new num_edges().
  *
  * Performs O(num_nodes() + num_edges()) operations.
  */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a,b)) {
      Edge e;
      size_type edge_idx;

      // delete adj_edge
      std::vector<size_type> adj_a = adj_edges_.at(a.index());
      for (size_type i=0; i<a.degree(); ++i) {
        e = edge(adj_a[i]);
        if (((a == e.node1()) && (b == e.node2())) || ((a == e.node2()) && (b == e.node1()))) {
          edge_idx = adj_a[i];
          adj_edges_[a.index()][i] = adj_edges_[a.index()].back();
          adj_edges_[a.index()].pop_back();
          break;
        }
      }
      std::vector<size_type> adj_b = adj_edges_.at(b.index());
      for (size_type i=0; i<b.degree(); ++i) {
        e = edge(adj_b[i]);
        if (((a == e.node1()) && (b == e.node2())) || ((a == e.node2()) && (b == e.node1()))) {
          adj_edges_[b.index()][i] = adj_edges_[b.index()].back();
          adj_edges_[b.index()].pop_back();
          break;
        }
      }

      // delete edge and value
      edges_[edge_idx] = edges_.back();
      edges_.pop_back();
      edge_values_[edge_idx] = edge_values_.back();
      edge_values_.pop_back();

      // adjust edge index
      size_type uid_0 = edges_[edge_idx][0];
      size_type uid_1 = edges_[edge_idx][1];
      for (auto it = adj_edges_[uid_0].begin(); it != adj_edges_[uid_0].end(); ++it) {
        if ((*it) == num_edges())
          (*it) = edge_idx;
      }
      for (auto it = adj_edges_[uid_1].begin(); it != adj_edges_[uid_1].end(); ++it) {
        if ((*it) == num_edges())
          (*it) = edge_idx;
      }
      return edge_idx;
    }
    return num_edges();
  }

  /** Remove an edge from a graph.
  * @param[in] e     An edge.
  * @return          The new index of the edge with old index old num_edges() - 1.
  *
  * @pre    @a e is a valid edge in the graph.
  * @post   The edge @a e is removed from the adjacent edges vector for both @a e.node1()
  *             and @a e.node2().
  * @post   The edge @a e is removed from the vector storing all the edges.
  * @post   The corresponding edge value is removed from the vector storing all the
  *             edge values.
  * @post   If the edge with old index old num_edges() - 1 is not removed, result is
  *             the new index of this edge. Otherwise, result is the new num_edges().
  *
  * Performs O(num_nodes() + num_edges()) operations.
  */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from a graph.
  * @param[in] e_it  An edge_iterator.
  * @return          The new edge_iterator referring to the edge with old index
  *                      old num_edges() - 1.
  *
  * @pre    @a e_it is an edge_iterator in range [old edge_begin(), old edge_end()).
  * @post   The edge *(@a e_it) is removed from the adjacent edges vector for both
  *             *(@a e_it).node1() and *(@a e_it).node2().
  * @post   The edge *(@a e_it) is removed from the vector storing all the edges.
  * @post   The corresponding edge value is removed from the vector storing all the
  *             edge values.
  * @post   If the edge with old index old num_edges() - 1 is not removed, result is
  *             the new edge_iterator referring to this edge. Otherwise, result is
  *             new edge_end().
  *
  * Performs O(num_nodes() + num_edges()) operations.
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type e_idx = remove_edge(*e_it);
    if (e_idx < num_edges())
      return EdgeIterator(this, e_idx);
    else
      return edge_end();
  }

  /** Remove a node from a graph.
  * @param[in] n   A node.
  * @return        The new index of the node with old index old num_nodes() - 1.
  *
  * @pre    @a n is a valid node in the graph.
  * @post   All the adjacent edges of @a n are removed.
  * @post   @a n is removed from the vector storing all the nodes.
  * @post   The corresponding node value is removed from the vector storing all the
  *             node values.
  * @post   node(i).index() == i for all i in the range [0, new num_nodes()).
  * @post   node(n.index()) == n for all node n.
  * @post   If the node with old index old num_nodes() - 1 is not removed, result is
  *             the new index of this node. Otherwise, result is the new num_nodes().
  *
  * Performs O(num_nodes()) operations.
  */
  size_type remove_node(const Node& n) {
    size_type node_idx = n.index();
    // remove adj_edge
    while (adj_edges_[node_idx].size() > 0) {
      remove_edge(edge(adj_edges_[node_idx].back()));
    }

    // remove node and value
    nodes_[node_idx] = nodes_.back();
    nodes_.pop_back();
    values_[node_idx] = values_.back();
    values_.pop_back();

    // adjust node index
    size_type node_end = num_nodes();
    adj_edges_[node_idx] = adj_edges_[node_end];
    adj_edges_.erase(node_end);
    if (node_idx < node_end) {
      for (auto it = adj_edges_[node_idx].begin(); it != adj_edges_[node_idx].end(); ++it) {
        if (edges_[(*it)][0] == node_end)
          edges_[(*it)][0] = node_idx;
        else
          edges_[(*it)][1] = node_idx;
      }
    }
    return node_idx;
  }

  /** Remove a node from a graph.
  * @param[in] n_it  A node_iterator.
  * @return          The new node_iterator referring to the node with old index
  *                      old num_nodes() - 1.
  *
  * @pre    @a n_it is a valid node_iterator in range [old node_begin(), old node_end()).
  * @post   All the adjacent edges of (*@a n) are removed.
  * @post   (*@a n) is removed from the vector storing all the nodes.
  * @post   The corresponding node value is removed from the vector storing all the
  *             node values.
  * @post   node(i).index() == i for all i in the range [0, new num_nodes()).
  * @post   node(n.index()) == n for all node n.
  * @post   If the node with old index old num_nodes() - 1 is not removed, result is
  *             the new node_iterator referring to this node. Otherwise, result is
  *             the new node_end().
  *
  * Performs O(num_nodes()) operations.
  */
  node_iterator remove_node(node_iterator n_it) {
    size_type n_idx = remove_node(*n_it);
    if (n_idx < num_nodes())
      return NodeIterator(this, n_idx);
    else
      return node_end();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> nodes_;
  std::vector<node_value_type> values_;
  std::vector<std::vector<size_type>> edges_;
  std::vector<edge_value_type> edge_values_;
  std::map<size_type, std::vector<size_type>> adj_edges_;

};

#endif // CME212_GRAPH_HPP
