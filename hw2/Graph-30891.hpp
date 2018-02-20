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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Data structures that store information for the Node and Edge class
  struct internal_node;
  struct internal_edge;

  //?Disable copy and assignment?

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Type of node value. */
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
   * Return type of Graph::Node::index(), Graph::num_nodes(),
   * Graph::num_edges(), and argument type of Graph::node(size_type)
   */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
  }

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
    Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    node_value_type& value() {
      return fetch().value;
    }

    const node_value_type& value() const {
      return fetch().value;
    }

    size_type degree() const {
      return graph_->node_edges_map_[index_].size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, graph_->node_edges_map_[index_].begin());
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, graph_->node_edges_map_[index_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_->has_node(n) && (n.index_ == index_);
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
      return graph_->has_node(n) && (n.index_ < index_);
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date Position information
    // @pre this node is valid (0 <= index_ < graph_->num_nodes())
    internal_node& fetch() const {
      assert(index_ < graph_->num_nodes());
      return graph_->nodes_[index_];
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type index_;

    /** Private Constructor */
    Node(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), index_(index) {
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
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node
    (const Point& position, const node_value_type& value = node_value_type()) {
    // HW1: YOUR CODE HERE
    // Add the new internal node to graph
    nodes_.push_back({position, value});
    // Return a node that associates with the new internal one
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this;
  }

  /** Remove a node from the graph
   * @param  n  a node object
   * @return  index of the node removed
   *
   * @pre  such a node exists in this graph
   * @post  has_node(n) == false
   * @post  new num_nodes() == old num_nodes() - 1
   * @post  incident edges of the removed node are also removed
   * @post  @return in any container of the graph as a node index
   *        now refers to the formerly last node (with the largest index)
   * @post  all other nodes remain the same
   *
   * Two node_iterators are invalidated:
   *   the one that used to point at the last node, and
   *   the one that used to point at the removed node,
   *   which now points to the former's value
   * Same for outstanding node objects
   * A large number of incident_iterators are invalidated
   * Edge_iterators and outstanding edge objects (proportional to
   * degree of the node) are invalidated
   *
   * Complexity: No more than O(num_nodes()), assuming small degrees
   *
   */
  size_type remove_node(const Node& n) {
    assert(has_node(n));

    // Remove all the edges incident to the node
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it)
      remove_edge(*it);

    int i = n.index_;
    // Update affected edges
    node_type last = node(num_nodes() - 1);
    for (incident_iterator it = last.edge_begin(); it != last.edge_end(); ++it) {
      node_edges_map_[(*it).node2().index()][i]
        = node_edges_map_[(*it).node2().index()][num_nodes() - 1];
      node_edges_map_[(*it).node2().index()].erase(num_nodes() - 1);
      internal_edge& e = (*it).fetch();
      if (e.node1 == num_nodes() - 1)
        e.node1 = i;
      else
        e.node2 = i;
    }

    // Remove from map
    node_edges_map_[i] = node_edges_map_[num_nodes() - 1];
    node_edges_map_.erase(num_nodes() - 1);

    // Remove from vector
    nodes_[i] = nodes_.back();
    nodes_.pop_back();

    return i;
  }

  /** Remove a node from the graph
   * @param  n  a node object
   * @return  index of the node removed
   *
   * @pre  such a node exists in this graph
   * @post  has_node(n) == false
   * @post  new num_nodes() == old num_nodes() - 1
   * @post  incident edges of the removed node are also removed
   * @post  @return in any container of the graph as a node index
   *        now refers to the formerly last node (with the largest index)
   * @post  all other nodes remain the same
   *
   * Two node_iterators are invalidated:
   *   the one that used to point at the last node, and
   *   the one passed in, which now points to the former's value
   * Same for outstanding node objects
   * A large number of incident_iterators are invalidated
   * Edge_iterators and outstanding edge objects (proportional to
   * degree of the node) are invalidated
   *
   * Complexity: No more than O(num_nodes()), assuming small degrees
   *
   */
  size_type remove_node(node_iterator n_it) {
    assert(has_node(*n_it));

    // Remove all the edges incident to the node
    for (incident_iterator& it = n_it->edge_begin(); it != n_it->edge_end(); ++it)
      remove_edge(*it);

    int i = n_it->index_;
    // Update affected edges
    node_type last = node(num_nodes() - 1);
    for (incident_iterator it = last.edge_begin(); it != last.edge_end(); ++it) {
      internal_edge& e = (*it).fetch();
      if (e.node1 == num_nodes() - 1)
        e.node1 = i;
      else
        e.node2 = i;
    }

    // Remove from map
    node_edges_map_[i] = node_edges_map_[num_nodes()];
    node_edges_map_.erase(num_nodes());
    
    // Remove from vector
    nodes_[i] = nodes_.back();
    nodes_.pop_back();

    return n_it;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    node_type node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);
    }

    edge_value_type& value() {
      return fetch().value;
    }

    const edge_value_type& value() const {
      return fetch().value;
    }

    double length() const {
      return norm(node2().position() - node1().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      return graph_ == e.graph_ && index_ == e.index_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      return graph_ < e.graph_ || (graph_ == e.graph_ && index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type index_;
    size_type node1_;
    size_type node2_;

    // Private constructor
    Edge(const graph_type* graph, size_type index,
      size_type node1, size_type node2)
      : graph_(const_cast<graph_type*>(graph)), index_(index),
      node1_(node1), node2_(node2) {
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date information
    // @pre this node is valid (0 <= index_ < graph_->num_edges())
    internal_edge& fetch() const {
      if (index_ >= graph_->num_edges())
        std::cout << node1_ << ' ' << node2_;
      assert(index_ < graph_->num_edges());
      return graph_->edges_[index_];
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
    assert(0 <= i && i < num_edges());
    return Edge(this, i, edges_[i].node1, edges_[i].node2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    assert(has_node(a) && has_node(b));
    size_type s = std::min(a.index_, b.index_);
    size_type l = std::max(a.index_, b.index_);
    if (node_edges_map_.count(s) == 1)
      if (node_edges_map_.at(s).count(l) == 1)
        return true;
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
  edge_type add_edge(const node_type& a, const node_type& b, const edge_value_type& value = edge_value_type()) {
    assert(has_node(a) && has_node(b) && (a != b));

    // Check if such edge already exists
    size_type a_i = a.index_;
    size_type b_i = b.index_;
    if (node_edges_map_.count(a_i) == 1)
      if (node_edges_map_[a_i].count(b_i) == 1)
        return Edge(this, node_edges_map_[a_i][b_i], a_i, b_i);

    // If not, add a new one
    node_edges_map_[a_i][b_i] = num_edges();
    node_edges_map_[b_i][a_i] = num_edges();
    edges_.push_back({a_i, b_i, value});

    // Return an edge that associates with the new internal one
    return Edge(this, num_edges(), a_i, b_i);
  }

  /** Remove an edge from the graph
   * @param  n1, n2  two node objects passed by reference
   * @return  index of the removed edge
   *
   * @pre  an edge between @a n1 and @a n2 exists in this graph
   * @post  has_edge(@a n1, @a n2) == false
   * @post  new num_edges() == old num_edges() - 1
   * @post  @return in any container of the graph as an edge index
   *        now refers to the formerly last edge (with the largest index)
   * @post  all other edges remain the same
   *
   * Two edge_iterators are invalidated:
   *   the one that used to point at the last edge, and
   *   the one that used to point at the removed edge,
   *   which now points to the former
   * Same for outstanding edge objects
   * A large number of incident_iterators are invalidated
   * Node objects and iterators are untouched
   *
   * Complexity: No more than O(log(num_nodes() + log(num_edges())),
   * depending on the degrees of the two nodes of the removed edge
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    assert(has_edge(n1, n2));

    int i = node_edges_map_[n1.index_][n2.index_];
    // Remove from map
    node_edges_map_[n1.index_].erase(n2.index_);
    node_edges_map_[n2.index_].erase(n1.index_);
    
    // Remove from vector
    edges_[i] = edges_.back();
    edges_.pop_back();

    // Update indexes
    node_edges_map_[edges_[i].node1][edges_[i].node2] = i;
    node_edges_map_[edges_[i].node2][edges_[i].node1] = i;

    return i;
  }
  /** Remove an edge from the graph
   * @param  e  an edge object passed by reference
   * @return    index of the removed edge
   *
   * @pre  an edge between @a e.node1() and @a e.node2() exists in this graph
   * @post  has_edge(@a e.node1(), @a e.node2()) == false
   * @post  new num_edges() == old num_edges() - 1
   * @post  @return in any container of the graph as an edge index
   *        now refers to the formerly last edge (with the largest index)
   * @post  all other edges remain the same
   *
   * Two edge_iterators are invalidated:
   *   the one that used to point at the last edge, and
   *   the one that used to point at the removed edge,
   *   which now points to the former
   * Same for outstanding edge objects
   * A large number of incident_iterators are invalidated
   * Node objects and iterators are untouched
   *
   * Complexity: No more than O(log(num_nodes() + log(num_edges())),
   * depending on the degrees of the two nodes of the removed edge
   */
  size_type remove_edge(const Edge& e) {
    assert(has_edge(e.node1(), e.node2()));

    int i = e.index_;
    // Remove from map
    node_edges_map_[e.node1_].erase(e.node2_);
    node_edges_map_[e.node2_].erase(e.node1_);
    
    // Remove from vector
    edges_[i] = edges_.back();
    edges_.pop_back();

    // Update indexes
    node_edges_map_[edges_[i].node1][edges_[i].node2] = i;
    node_edges_map_[edges_[i].node2][edges_[i].node1] = i;

    return i;
  }

  /** Remove an edge from the graph
   * @param  e_it  an edge_iterator to the edge we want to remove
   * @return  iterator to the new edge sitting at this memory address
   *
   * @pre  such an edge exists in this graph
   * @post  has_edge(e_it->node1_, e_it->node2_) == false
   * @post  new num_edges() == old num_edges() - 1
   * @post  @return in any container of the graph as an edge index
   *        now refers to the formerly last edge (with the largest index)
   * @post  all other edges remain the same
   *
   * Two edge_iterators are invalidated:
   *   the one that used to point at the last edge, and
   *   the one passed in, which now points to the former's value
   * Same for outstanding edge objects
   * A large number of incident_iterators are invalidated
   * Node objects and iterators are untouched
   *
   * Complexity: No more than O(log(num_nodes() + log(num_edges())),
   * depending on the degrees of the two nodes of the removed edge
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    assert(has_edge(e_it->node1_, e_it->node2_));

    int i = e_it.index_;
    // Remove from map
    node_edges_map_[e_it->node1_].erase(e_it->node2_);
    node_edges_map_[e_it->node2_].erase(e_it->node1_);
    
    // Remove from vector
    edges_[i] = edges_.back();
    edges_.pop_back();

    // Update indexes
    edges_[i].index_ = i;
    node_edges_map_[edges_[i].node1][edges_[i].node2] = i;
    node_edges_map_[edges_[i].node2][edges_[i].node1] = i;

    return e_it; // Still has the old data, but pointing to a new edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    value_type operator*() const {
      return Node(graph_, index_);
    }

    node_iterator& operator++() {
      index_++;
      return *this;
    }

    bool operator==(const NodeIterator& it2) const {
      return graph_ == it2.graph_ && index_ == it2.index_;
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    /** Private Constructor */
    NodeIterator(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    value_type operator*() const {
      size_type j = iter_->second;
      size_type node2_;
      if (node1_ == graph_->edges_[j].node1)
        node2_ = graph_->edges_[j].node2;
      else
        node2_ = graph_->edges_[j].node1;
      return Edge(graph_, j, node1_, node2_);
    }

    incident_iterator& operator++() {
      ++iter_;
      return *this;
    }

    bool operator==(const incident_iterator& it2) const {
      return graph_ == it2.graph_ && node1_ == it2.node1_
        && iter_ == it2.iter_;
    }

    /*
    Edge operator*() const {
      size_type j = graph_->node_edges_map_[node1_][index_];
      size_type node2_;
      if (node1_ == graph_->edges_[j].node1)
        node2_ = graph_->edges_[j].node2;
      else
        node2_ = graph_->edges_[j].node1;
      return Edge(graph_, j, node1_, node2_);
    }

    incident_iterator& operator++() {
      index_++;
      return *this;
    }

    bool operator==(const incident_iterator& it2) const {
      return graph_ == it2.graph_ && node1_ == it2.node1_
        && index_ == it2.index_;
    }
    */

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node1_;
    std::map<size_type, size_type>::iterator iter_;

    /** Private Constructor */
    IncidentIterator(const graph_type* graph, size_type node1,
      std::map<size_type, size_type>::iterator iter)
      : graph_(const_cast<graph_type*>(graph)), node1_(node1), iter_(iter) {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    value_type operator*() const {
      auto& e = graph_->edges_[index_];
      return Edge(graph_, index_, e.node1, e.node2);
    }
    edge_iterator& operator++() {
      ++index_;
      return *this;
    }
    bool operator==(const EdgeIterator& it2) const {
      return graph_ == it2.graph_ && index_ == it2.index_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  void setAnchors(std::vector<Point> anchors) {
    // Locate the anchor points
    std::vector<unsigned> indexes;
    for (Point& a : anchors)
      for (auto it = node_begin(); it != node_end(); ++it)
        if ((*it).position() == a)
          indexes.push_back((*it).index());

    // Switch the anchor indexes to the top of the container in graph
    for (unsigned i = 0; i != indexes.size(); ++i) {
      node_type original = node(indexes[i]);
      node_type target = node(i);

      std::vector<std::vector<size_type>> temp;
      for (auto it = original.edge_begin(); it != original.edge_end(); ++it) {
        internal_edge& e = (*it).fetch();
        if (has_edge((*it).node2(), target)) {
          swap(node_edges_map_[(*it).node2().index()][i],
            node_edges_map_[(*it).node2().index()][indexes[i]]);
        } else {
          temp.push_back({(*it).node2().index(), i, indexes[i]});
        }

        if (e.node1 == indexes[i])
          e.node1 = i;
        else
          e.node2 = i;
      }

      for (auto it = target.edge_begin(); it != target.edge_end(); ++it) {
        internal_edge& e = (*it).fetch();
        if (!has_edge((*it).node2(), original))
          temp.push_back({(*it).node2().index(), indexes[i], i});

        if (e.node1 == i)
          e.node1 = indexes[i];
        else
          e.node2 = indexes[i];
      }

      for (auto& t : temp) {
        node_edges_map_[t[0]][t[1]] = node_edges_map_[t[0]][t[2]];
        node_edges_map_[t[0]].erase(t[2]);
      }

      swap(node_edges_map_[i], node_edges_map_[indexes[i]]);
      swap(nodes_[i], nodes_[indexes[i]]);
    }
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    // Helper function
  template <typename T>
  void swap(T& a, T& b) {
    T tmp {std::move(a)};
    a = std::move(b);
    b = std::move(tmp);
  }

  struct internal_node {
    Point position;
    node_value_type value;
  };

  struct internal_edge {
    size_type node1;
    size_type node2;
    edge_value_type value;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::map<size_type, std::map<size_type, size_type>> node_edges_map_;
};

#endif // CME212_GRAPH_HPP
