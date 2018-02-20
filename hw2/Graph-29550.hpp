#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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

  // Predeclare the internal struct
  struct internal_node;
  // Declare Vector internal data
  std::vector<internal_node> nodes;
  // Declare Vector of indices
  std::vector<unsigned> i2u;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
  using edge_value_type = E;

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

    /** Update this node's position. */
    Point& position() {
      return graph_->nodes[uid_].point;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes[uid_].index;
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
      return incident_iterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (uid_ == (size_type) graph_->i2u[n.index()]);
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
    return i2u.size();
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
    n.index = i2u.size();
    nodes.push_back(n);
    // Add node id
    size_type uid = nodes.size() - 1;
    i2u.push_back(uid);
    // Return node
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() < this->num_nodes()) {
      return this == n.graph_;
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
    return Node(this, i2u[i]);
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

    /** Set this edge's value. */
    edge_value_type& value() {
      return graph_->edge_values[n1_<n2_?n1_:n2_][n1_<n2_?n2_:n1_];
    }

    /** Return this edge's value. */
    const edge_value_type& value() const {
      return graph_->edge_values[n1_<n2_?n1_:n2_][n1_<n2_?n2_:n1_];
    }

    /** Return current edge length */
    double length() const {
      return norm(this->direction());
    }

    /** Return vector direction of edge */
    Point direction() const {
      Node n1 = this->node1();
      Node n2 = this->node2();
      return n1.position() - n2.position();
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, n2_);
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
      // check if edges are equal
      if (*this != e) {
        // check if first node are equal
        if (this->node1() != e.node1()) {
          return this->node1() < e.node1();
        }
        return this->node2() < e.node2();
      }
      return false;
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
    return total_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    edge_iterator it = edge_begin();
    while (i-- != 0) {
      ++it;
    }
    return *it;
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
      // Update Node's edges
      nodes[i2u[a.index()]].edges.push_back(i2u[b.index()]);
      nodes[i2u[b.index()]].edges.push_back(i2u[a.index()]);
      // increment total number of edges
      total_edges += 1;
    }
    // Return Edge
    return Edge(this, i2u[a.index()], i2u[b.index()]);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    i2u.clear(); 
    edge_values.clear();
    total_edges = 0;
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
    return NodeIterator(this, this->num_nodes());
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
      size_type index = graph_->nodes[uid_].edges[pos_];
      //return graph_->edge(index);
      return Edge(graph_, uid_, index);
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
      return Edge(graph_, graph_->i2u[n1_], graph_->nodes[graph_->i2u[n1_]].edges[n2_]);
    }

    /** Increment position of this EdgeIterator. */
    EdgeIterator& operator++() {
      // check that current node still has unchecked edges
      if (n2_ < graph_->nodes[graph_->i2u[n1_]].edges.size() - 1) {
        n2_ += 1;
        // only increment through edges (a,b) with a < b
        if (!(graph_->i2u[n1_] < graph_->nodes[graph_->i2u[n1_]].edges[n2_])) {
          return ++(*this);
        }
        return *this;
      } else {
        // go to next node
        n1_ += 1;
        n2_ = 0;
        return *this;
      }
    }

    /** Check if this EdgeIterator is equal to @x EdgeIterator.
    *    @param x EdgeIterator to compare to.
    */
    bool operator==(const EdgeIterator& x) const {
      return graph_ == x.graph_ && n1_ == x.n1_ && n2_ == x.n2_;
    }

   private:
    friend class Graph;
    // Pointer back to the Graph
    graph_type* graph_;
    // Pointer back to the Graph
    size_type n1_;
    // Current nodes unique identification number
    size_type n2_;
    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type n1, size_type n2)
        : graph_(const_cast<graph_type*>(graph)), n1_(n1), n2_(n2) {
    }
  };

  /** Return a begining iterator over edges in this graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }

  /** Return an end iterator over edges in this graph. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }





  /** Remove a node from the graph and return 1 or return 0 if node graph does not have node.
   * 
   * @param[in] n A node to be removed
   * @result a size_type with 1 implies success and 0 implies node is not in graph
   *
   * @pre if @a n is in graph then 0 <= j < g.num_nodes()
   * @post For all i < g.size() - 1 and not equal to @ n.index():
   *    - old i2u[i] = new i2u[i]
   *    - nodes[i2u[i]].index = i
   * @post For i = @ n.index():
   *    - new i2u[i] = old i2u[g.size() - 1]
   *    - nodes[i2u[i]].index = i
   * @post new g.num_nodes() = old g.num_nodes() - 1
   *
   *
   * Complexity: O(1)
   */
  size_type remove_node(const Node& n) {
    if (this->has_node(n)) {
      // remove incident edges
      auto it = n.edge_begin();
      auto end = n.edge_end();
      while (it != end) {
        this->remove_edge(*it);
        end = n.edge_end();
      }
      // remove node index from i2u
      i2u[n.index()] = i2u.back();
      // update internal nodes
      nodes[i2u[n.index()]].index = n.index();
      // shrink i2u
      i2u.pop_back();
      // successful
      return 1;
    }
    // graph does not have node
    return 0;
  }


  /** Remove a node from the graph at input node iterator.
   * 
   * @param[in] n_it A node iterator at the element of the node to be removed
   * @result an iterator in the range [@it, new end())
   *
   * @pre @a it != end()
   * @post Invalidates all iterators in the range [@a it, old end())
   * @post For all i < g.size() - 1 and not equal to @ n.index():
   *    - old i2u[i] = new i2u[i]
   *    - nodes[i2u[i]].index = i
   * @post For i = @ n.index():
   *    - new i2u[i] = old i2u[g.size() - 1]
   *    - nodes[i2u[i]].index = i
   * @post new g.num_nodes() = old g.num_nodes() - 1
   *
   *
   * Complexity: O(1)
   */
  node_iterator remove_node(node_iterator n_it) {
    this->remove_node(*n_it);
    return n_it;
  }


  /** Remove edge between two nodes in the graph.
   * 
   * @param[in] a First node of edge
   * @param[in] b Second node of edge
   * @result a size_type with 1 implies success and 0 implies edge is not in graph
   *
   * @pre i2u[@a b.index()] in edges of @a
   * @pre i2u[@a a.index()] in edges of @b
   *
   * @post new g.num_edges() = old g.num_edges() - 1
   *
   * @post new a.degree() = old a.degree() - 1
   * @post new b.degree() = old b.degree() - 1
   *
   * Complexity: O(num_nodes) (average is O(1))
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (this->has_edge(a, b)) {
      // remove b from a's edges
      auto& a_edges = nodes[i2u[a.index()]].edges;
      auto b_index = find(a_edges.begin(), a_edges.end(), i2u[b.index()]);
      a_edges.erase(b_index);
      // remove a from b's edges
      auto& b_edges = nodes[i2u[b.index()]].edges;
      auto a_index = find(b_edges.begin(), b_edges.end(), i2u[a.index()]);
      b_edges.erase(a_index);
      // decrement edge counter
      total_edges -= 1;
      // successful
      return 1;
    }
    // graph does not have edge
    return 0;
  }


  /** Remove edge between in the graph.
   * 
   * @param[in] e An edge to be removed
   * @result a size_type with 1 implies success and 0 implies edge is not in graph
   *
   * @post new g.num_edges() = old g.num_edges() - 1
   *
   * @post new a.degree() = old a.degree() - 1
   * @post new b.degree() = old b.degree() - 1
   *
   * Complexity: O(num_nodes) (average is O(1))
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  /** Remove edge at current edge iterator.
   * 
   * @param[in] e_it An edge iterator at the edge to be removed
   * @result an iterator in the range [@it, new end())
   *
   * @pre @a e_it != end()
   * @post Invalidates all iterators in the range [@a it, old end())
   * @post new g.num_edges() = old g.num_edges() - 1
   * @post new a.degree() = old a.degree() - 1
   * @post new b.degree() = old b.degree() - 1
   *
   * Complexity: O(num_nodes) (average is O(1))
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    this->remove_edge(*e_it);
    return e_it;
  }


 private:

  // internal_node attributes
  struct internal_node {
    Point point;   // position of node
    node_value_type value;  // value of node
    size_type index;  // index in i2u
    std::vector<size_type> edges;  // neighboring nodes index
  };

  // internal map of edge values
  std::map<size_type, std::map<size_type, edge_value_type>> edge_values;

  // number of edges
  size_type total_edges = 0;

};

#endif // CME212_GRAPH_HPP
