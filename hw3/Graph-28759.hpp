#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V, typename E>
class Graph{
 private:
  struct Node_vals;
  struct Edge_vals;

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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Predeclaration of Node Iterator. */
  class NodeIterator;
  /** Synonym for NodeIterator. */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Predeclaration of Incident Iterator. */
  class IncidentIterator;
  /** SYnonym for IncidentIterator. */
  using incident_iterator = IncidentIterator;

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Attributes of Graph class. For the beginning I have two vectors, I can combine them into one struct later.

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

 public:
  /** Construct an empty graph. */
  Graph()
    : coordinate_(), adjacency_(), edge_count_(0), ind2uid_() {
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

  class Node:private totally_ordered<Node> {
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
      // HW0: YOUR CODE HERE
      return graph_->coordinate_[uid_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->coordinate_[uid_].index_;
    }

    /** Obtain node's value.
     *  Return this node's valye type
     */

    node_value_type& value() {
      return graph_->coordinate_[uid_].value_;
    }

    const node_value_type& value() const{
      return graph_->coordinate_[uid_].value_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((uid_ == n.uid_) && (graph_ == n.graph_));
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
      if ((graph_ == n.graph_) && (uid_ < n.uid_)) {
	       return true;
      }
      if (graph_ < n.graph_) {
	       return true;
      }
      return false;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    size_type degree() const {
      return graph_->adjacency_[uid_].size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, this->degree());
    }

    // HW2: YOUR CODE HERE
    Point& position() {
      return graph_->coordinate_[uid_].point_;
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
    return ind2uid_.size();
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
    Node_vals new_node;
    new_node.point_ = position;
    new_node.index_ = ind2uid_.size();
    new_node.value_ = value;
    coordinate_.push_back(new_node);
    adjacency_.resize(adjacency_.size() + 1);
    ind2uid_.push_back(coordinate_.size() - 1);
    return Node(this, coordinate_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((this == n.graph_) && (n.graph_->coordinate_[n.uid_].index_ < num_nodes()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   * Complexity: O(1).
   */

  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, ind2uid_[i]);
  }

  // HW2L YOUR CODE HERE
  /**Remove specified node and incident edges
   * @param[in] n    Node to be removed.
   * @pre has_node(@a n) == 1
   * @post has_node(@a n) == 0
   * @post NodeIterators and EdgeIterators for @a n are invalidated.
   * @post IncidentIterators going through @a n's incident edges are invalidated.
   * @post If (e.node1() == @a n || e.node2() == @a n) then e is removed.
   * @post new size() == old size() - 1
   * @return size_type  Original index of the node removed.
   * Complexity: O(num_nodes())
  */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }

    // Delete edges adjacent to node n
    while (!adjacency_[n.uid_].empty()) {
      size_type del_nod_ind = adjacency_[n.uid_][0].node_id_;
      remove_edge(n, Node(n.graph_, del_nod_ind));
    }

    // Deleting the element from index to uid vector
    auto iter = ind2uid_.begin();
    while (iter != ind2uid_.end()) {
      if (*iter == n.uid_) {
        *iter = ind2uid_.back();
        ind2uid_.pop_back();
        break;
      } else {++iter;}
    }

    // decrementing the indices
    for (unsigned i = coordinate_[n.uid_].index_; i < ind2uid_.size(); ++i) {
      Node_vals& in = coordinate_.at(ind2uid_[i]);
      in.index_ = i;
    }
    return *iter;
  }

  /**Remove specified node iterator
   * @param[in] n_it    NodeIterator to be removed.
   * @pre has_node(@a n) == 1
   * @post has_node(@a n) == 0
   * @post NodeIterators and EdgeIterators for @a n are invalidated.
   * @post IncidentIterators going through @a n's incident edges are invalidated.
   * @post If (e.node1() == @a n || e.node2() == @a n) then e is removed.
   * @post new size() == old size() - 1
   * @return node_iterator pointing to where original node was pre-removal.
   * Complexity: O(num_nodes())
  */

  node_iterator remove_node(node_iterator n_iter) {
    Node n = (*n_iter);
    size_type i = remove_node(n);
    (void) i;
    return n_iter;
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(uid2_);
    }

    edge_value_type& value() {
      Edge e;
      if (uid1_ < uid2_) {
        e = Edge(graph_, uid1_, uid2_);
      } else {
        e = Edge(graph_, uid2_, uid1_);
      }
      size_type n = e.node2().uid_;

      for (auto iter = graph_->adjacency_[n].begin();
            iter != graph_->adjacency_[n].end(); ++iter) {
        if ((*iter).node_id_ == e.node1().uid_)
          return (*iter).value_;
      }
      return graph_->adjacency_[e.uid1_][e.uid2_].value_;
    }

    const edge_value_type& value() const {
      Edge e;
      if (uid1_ < uid2_) {
        e = Edge(graph_, uid1_, uid2_);
      } else {
        e = Edge(graph_, uid2_, uid1_);
      }
      size_type n = e.node2().uid_;

      for (auto iter = graph_->adjacency_[n].begin();
            iter != graph_->adjacency_[n].end(); ++iter) {
        if ((*iter).node_id_ == e.node1().uid_)
          return (*iter).value_;
      }
      return graph_->adjacency_[e.uid1_][e.uid2_].value_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((((node1() == e.node1()) && (node2() == e.node2())) ||
	       ((node1() == e.node2()) && (node2() == e.node1()))) &&
	      (graph_ == e.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      auto edge_1 = std::make_tuple(graph_, std::min(uid1_, uid2_), std::max(uid1_, uid2_));
      auto edge_2 = std::make_tuple(e.graph_, std::min(e.uid1_, e.uid2_), std::max(e.uid1_, e.uid2_));
      return edge_1 < edge_2;
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid1_;
    size_type uid2_;
    Edge(const Graph* graph, size_type uid1, size_type uid2)
      : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_count_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    edge_iterator iter = edge_begin();
    for (size_type j = 0; j < i; ++j)
      ++iter;
    return *iter;
	}

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (auto iter = adjacency_[a.uid_].begin(); iter < adjacency_[a.uid_].end(); ++iter) {
      if ((*iter).node_id_ == b.uid_) {
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
    // HW0: YOUR CODE HERE
    if (!has_edge(a, b)) {
      std::pair <size_type, size_type> new_edge;
      new_edge.first  = a.index();
      new_edge.second = b.index();
      if (new_edge.first >= adjacency_.size()) {
	adjacency_.resize(new_edge.first + 1);
      }
      if (new_edge.second >= adjacency_.size()) {
	adjacency_.resize(new_edge.second + 1);
      }
      adjacency_[a.index()].push_back({new_edge.second, edge_value_type()});
      adjacency_[b.index()].push_back({new_edge.first, edge_value_type()});
      ++edge_count_;
      return Edge(this, new_edge.first, new_edge.second);
    }
    return Edge(this, a.index(), b.index());
  }


  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator> {
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

    Node operator*() const {
//      return Node(graph_, graph_->coordinate_[graph_->ind2uid_[indx_]].index_);
      return graph_->node(indx_);
    }

    NodeIterator& operator++() {
//      if (indx_ < graph_->ind2uid_.size()) {
//	       ++indx_;
//      }
      ++indx_;
      return *this;
    }

    bool operator==(const NodeIterator& n) const {
      return (graph_ == n.graph_ && indx_ == n.indx_);
//      return (graph_ == n.graph_ && )
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type indx_;

    NodeIterator(const Graph* graph, size_type indx):
      graph_(const_cast<Graph*>(graph)), indx_(indx) {
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
    return NodeIterator(this, this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */

  class IncidentIterator:private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge; 			// Element type
    using pointer           = Edge*;			// Pointers to elements
    using reference         = Edge&;			// Reference to elements
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

    Edge operator*() const {
      return Edge(graph_, uid1_, graph_->adjacency_[uid1_][uid2_].node_id_);
    }

    IncidentIterator& operator++() {
      if (uid2_ < graph_->adjacency_[uid1_].size()) {
	       ++uid2_;
      }
      return *this;
    }

    bool operator==(const IncidentIterator& incident_iter) const {
      return ((graph_ == incident_iter.graph_) &&
	      (uid1_  == incident_iter.uid1_)  &&
	      (uid2_  == incident_iter.uid2_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type uid1_;
    size_type uid2_;

    IncidentIterator(const Graph* graph, size_type uid1, size_type uid2):
      graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */

  class EdgeIterator:private totally_ordered<EdgeIterator> {
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

    Edge operator*() const {
      return Edge(graph_, graph_->adjacency_[uid1_][edge_id_].node_id_, uid1_);
    }

    EdgeIterator& operator++() {
      if (edge_id_ < graph_->adjacency_[uid1_].size()) {
	       ++edge_id_;
      }
      get_edge();
      return *this;
    }

    // Function that gets valid edge
    void get_edge() {
      long adj = graph_->adjacency_.size();
      if (uid1_ >= adj) return;
      if (graph_->adjacency_[uid1_].size() <= edge_id_) {
        ++uid1_;
        edge_id_ = 0;
      }

      while ((uid1_ < adj) && (graph_->adjacency_[uid1_].size() == 0)) {
        ++uid1_;
      }

      while ((uid1_ < adj) && (graph_->adjacency_[uid1_][edge_id_].node_id_ < uid1_)) {
        ++edge_id_;
        if (edge_id_ >= graph_->adjacency_[uid1_].size()) {
          ++uid1_;
          edge_id_ = 0;
        }
      }
      while ((uid1_ < adj) && (graph_->adjacency_[uid1_].size() == 0)) {
        ++uid1_;
      }
    }

    bool operator==(const EdgeIterator& e) const {
      return ((graph_ == e.graph_) &&
	      (uid1_  == e.uid1_)  &&
	      (edge_id_  == e.edge_id_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type uid1_;
    size_type edge_id_;

    EdgeIterator(const Graph* graph, size_type uid1, size_type edge_id_):
      graph_(const_cast<Graph*>(graph)), uid1_(uid1), edge_id_(edge_id_) {
        get_edge();
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


  edge_iterator edge_begin() const {
    return edge_iterator(this, 0, 0);
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, adjacency_.size(), 0);
  }


  // HW2: YOUR CODE HERE
  /**Remove an edge from endpoints
  * @param[in] Node @a n1, first node of the edge to be removed
  * @param[in] Node @a n2, second node of the edge to be removed
  * @pre (has_node(n1) && has_node(n2)) == 1
  * @post If (has_edge(n1, n2) == 1) then new num_edge() == old num_edge() - 1
  * @post If (has_edge(n1, n2) == 0) then new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return size_type  Index of the second node removed from the @a edges_ adjacency list
  * Complexity: O(num_nodes() + num_edges())
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    size_type rem_ind = 0;
    if (n1.uid_ >= adjacency_.size() || n2.uid_ >= adjacency_.size()) {
      return 0;
    }

    for (unsigned int i = 0; i < adjacency_[n1.uid_].size(); ++i) {
      if (adjacency_[n1.uid_][i].node_id_ == n2.uid_) {
        adjacency_[n1.uid_].erase(adjacency_[n1.uid_].begin() + i);
        break;
      }
    }

    for (unsigned int i = 0; i < adjacency_[n2.uid_].size(); ++i) {
      if (adjacency_[n2.uid_][i].node_id_ == n1.uid_) {
        adjacency_[n2.uid_].erase(adjacency_[n2.uid_].begin() + i);
        --edge_count_;
        rem_ind = i;
        break;
      }
    }
    return rem_ind;
  }

  /**Remove an edge
  * @param[in] Edge e, edge to be removed
  * @pre @a e is a valid edge
  * @post If (has_edge(@a e.node1(), @a e.node2()) == 1) then new num_edge() == old num_edge()-1
  * @post If (has_edge(@a e.node1(), @a e.node2()) == 0) then new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return size_type  Index of the second node removed from the @a edges_ adjacency list
  * Complexity: O(num_nodes() + num_edges())
  */

  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /**Remove an edge
  * @param[in] EdgeIterator @a e_it to be removed
  * @pre @a e_it is a valid edge iterator
  * @post If (has_edge(n1, n2) == 1) then new num_edge() == old num_edge() - 1
  * @post If (has_edge(n1, n2) == 0) then new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return EdgeIterator  Points to the next valid edge in the @a edges_ adjacency list.
  * Complexity: O(num_nodes() + num_edges())
  */

  edge_iterator remove_edge(edge_iterator e_it) {
    size_type ind = remove_edge(*e_it);
    (void) ind;
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    coordinate_.clear();
    adjacency_.clear();
    edge_count_ = 0;
    ind2uid_.clear();
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

 // Private attributes of Graph Class
 // Node structure contains information about the coordinate, index and type
 // of the node. Edge structure contains nodes that are adjacent to a given
 // node

  struct Node_vals{
    Point point_;
    size_type index_;
    node_value_type value_;
  };

  struct Edge_vals{
    size_type node_id_;
    edge_value_type value_;
  };

  std::vector<Node_vals> coordinate_;
  std::vector<std::vector<Edge_vals> > adjacency_;

  size_type edge_count_;
  std::vector<size_type> ind2uid_;

};

#endif // CME212_GRAPH_HPP




// #ifndef CME212_GRAPH_HPP
// #define CME212_GRAPH_HPP
//
// /** @file Graph.hpp
//  * @brief An undirected graph type
//  */
//
// #include <algorithm>
// #include <cassert>
// #include <utility>
// #include <vector>
//
// #include "CME212/Util.hpp"
// #include "CME212/Point.hpp"
//
//
// /** @class Graph
//  * @brief A template for 3D undirected graphs.
//  *
//  * Users can add and retrieve nodes and edges. Edges are unique (there is at
//  * most one edge between any pair of distinct nodes).
//  */
// template <typename V, typename E>
// class Graph {
//  private:
//
//   /** Predeclaring internal structs */
//   struct internal_edge;
//   struct internal_node;
//
//  public:
//
//   //
//   // PUBLIC TYPE DEFINITIONS
//   //
//   std::vector<internal_node> nodes_;
//   std::vector<std::vector<internal_edge>> adj_;
//   std::vector<unsigned> i2u_; // map between indices (i) and unique identifiers (uids).
//
//   /** Type of this graph. */
//   using graph_type = Graph<V,E>;
//
//   /** Predeclaration of Node type. */
//   class Node;
//   /** Synonym for Node (following STL conventions). */
//   using node_type = Node;
//   /** Specifies the type of the nodal value - provided by user. */
//   using node_value_type = V;
//
//   /** Predeclaration of Edge type. */
//   class Edge;
//   /** Synonym for Edge (following STL conventions). */
//   using edge_type = Edge;
//   /** Specifies the type of the edge value - also provided by user. */
//   using edge_value_type = E;
//
//   /** Type of node iterators, which iterate over all graph nodes. */
//   class NodeIterator;
//   /** Synonym for NodeIterator */
//   using node_iterator = NodeIterator;
//
//   /** Type of edge iterators, which iterate over all graph edges. */
//   class EdgeIterator;
//   /** Synonym for EdgeIterator */
//   using edge_iterator = EdgeIterator;
//
//   /** Type of incident iterators, which iterate incident edges to a node. */
//   class IncidentIterator;
//   /** Synonym for IncidentIterator */
//   using incident_iterator = IncidentIterator;
//
//   /** Type of indexes and sizes.
//       Return type of Graph::Node::index(), Graph::num_nodes(),
//       Graph::num_edges(), and argument type of Graph::node(size_type) */
//   using size_type = unsigned;
//
//   //
//   // CONSTRUCTORS AND DESTRUCTOR
//   //
//
//   /** Construct an empty graph. */
//   Graph() {}
//
//   /** Default destructor */
//   ~Graph() = default;
//
//   //
//   // NODES
//   //
//
//   /** @class Graph::Node
//    * @brief Class representing the graph's nodes.
//    *
//    * Node objects are used to access information about the Graph's nodes.
//    */
//   class Node : private totally_ordered<Node> {
//    public:
//
//     /** Construct an invalid node.
//      *
//      * Valid nodes are obtained from the Graph class, but it
//      * is occasionally useful to declare an @i invalid node, and assign a
//      * valid node to it later. For example:
//      *
//      * @code
//      * Graph::node_type x;
//      * if (...should pick the first node...)
//      *   x = graph.node(0);
//      * else
//      *   x = some other node using a complicated calculation
//      * do_something(x);
//      * @endcode
//      */
//     Node() {
//     }
//
//     /** Return this node's position. */
//     const Point& position() const {
//       return node_graph->nodes_[node_graph->i2u_[node_index]].position;
//     };
//
//     /** Return a mutable reference to this node's position. */
//     Point& position() {
//       return node_graph->nodes_[node_graph->i2u_[node_index]].position;
//     };
//
//     /** Return this node's index, a number in the range [0, graph_size). */
//     size_type index() const {
//         return node_index;
//       // return node_graph->nodes_[node_graph->i2u_[node_index]].index;
//     };
//
//     /** Method which returns the mutable value of the node. */
//     node_value_type& value(){
//         return node_graph->nodes_[node_graph->i2u_[node_index]].value;
//     };
//
//     /** Const method which returns an IMMUTABLE reference to the node value. */
//     const node_value_type& value() const {
//         return node_graph->nodes_[node_graph->i2u_[node_index]].value;
//     };
//
//     /** Returns the number of incident edges to this node. */
//     size_type degree() const{
//         return node_graph->adj_[node_graph->i2u_[node_index]].size();
//     };
//
//     /** Returns pointer to this nodes "first" incident edge. */
//     incident_iterator edge_begin() const {
//         return IncidentIterator(node_graph, node_index, 0);
//     };
//
//     /** Returns pointer one-past this nodes "last" incident edge. */
//     incident_iterator edge_end() const {
//         return IncidentIterator(node_graph, node_index, this->degree());
//     };
//
//     /** Test whether this node and @a n are equal.
//      *
//      * Equal nodes have the same graph and the same index.
//      */
//     bool operator==(const Node& n) const {
//       return node_graph == n.node_graph && node_index == n.index();
//     }
//
//     /** Test whether this node is less than @a n in a global order.
//      *
//      * This ordering function is useful for STL containers such as
//      * std::map<>. It need not have any geometric meaning.
//      *
//      * The node ordering relation must obey trichotomy: For any two nodes x
//      * and y, exactly one of x == y, x < y, and y < x is true.
//      */
//     bool operator<(const Node& n) const {
//       if (node_index < n.index()){
//         return true;
//       } else {
//         return (node_index == n.index() && node_graph < n.node_graph);
//       }
//     }
//
//    private:
//     // Allow Graph to access Node's private member data and functions.
//     friend class Graph;
//
//     Graph* node_graph;
//     size_type node_index;
//
//     /** Private Constructor */
//     Node(const Graph* graph, size_type index)
//         : node_graph(const_cast<Graph*>(graph)), node_index(index) {
//     }
//   };
//
//   /** Return the number of nodes in the graph.
//    *
//    * Complexity: O(1).
//    */
//   size_type size() const {
//     return i2u_.size();
//   }
//
//   /** Synonym for size(). */
//   size_type num_nodes() const {
//     return size();
//   }
//
//   /** Add a node to the graph, returning the added node.
//    * @param[in] position The new node's position
//    * @param[in] node_value The new node's value
//    * @post new num_nodes() == old num_nodes() + 1
//    * @post result_node.index() == old num_nodes()
//    *
//    * Complexity: O(1) amortized operations.
//    */
//   Node add_node(const Point& position, const node_value_type& node_value =
//       node_value_type()) {
//
//     internal_node new_node;
//     new_node.position = position;
//     new_node.value = node_value;
//     new_node.index = i2u_.size();
//
//     i2u_.push_back(nodes_.size());
//     assert(new_node.index == i2u_.size() - 1);
//
//     nodes_.push_back(new_node);
//     std::vector<internal_edge> new_vec;
//     adj_.push_back(new_vec);
//     assert(nodes_.size() == adj_.size());
//
//     return Node(this, new_node.index);
//   }
//
//   /** Determine if a Node belongs to this Graph
//    * @return True if @a n is currently a Node of this Graph
//    *
//    * Complexity: O(1).
//    */
//   bool has_node(const Node& n) const {
//       for (auto nit = this->node_begin(); nit != this->node_end(); ++nit){
//           if ((*nit).position() == n.position()){
//               return (this == n.node_graph && n.index() == (*nit).index());
//           }
//       }
//       return false;
//   }
//
//   /** Return the node with index @a i.
//    * @pre 0 <= @a i < num_nodes()
//    * @post result_node.index() == i
//    *
//    * Complexity: O(1).
//    */
//   Node node(size_type i) const {
//     return Node(this, i);        // Invalid node
//   }
//
//   /**  If @n is a node in graph: Removes @n and all its incident edges
//    * from graph
//    * Else: Returns 0.
//    *
//    * @pre @n is a valid node.
//    * @post new_num_nodes() = old_num_nodes() - 1
//    * @post new_num_edges() = old_num_edges() - n.degree()
//    * @return 0
//    *
//    * Complexity: O(num_edges())
//    */
//   size_type remove_node(const Node& n){
//       if (has_node(n)){
//           assert(nodes_[i2u_[n.index()]].index == n.index());
//
//           // Removing all neighboring edges
//           size_type degree = n.degree();
//           for (size_type i = 0; i < degree; i++){
//               size_type neighbor = adj_[i2u_[n.index()]][0].index;
//               remove_edge(n, this->node(neighbor));
//           }
//           assert(adj_[i2u_[n.index()]].size() == 0);
//
//           // Removing node from i2u_ vector
//           assert(nodes_[i2u_.back()].index == (i2u_.size()-1));
//
//           size_type last_index = nodes_[i2u_.back()].index;
//           nodes_[i2u_.back()].index = n.index();
//           i2u_[n.index()] = i2u_.back();
//           i2u_.pop_back();
//           assert(nodes_[i2u_[n.index()]].index == n.index());
//
//           // Updating the index to neighboring nodes of i2u_.back()
//           for (size_type i = 0; i < adj_[i2u_[n.index()]].size(); i++){
//               size_type neighbor = adj_[i2u_[n.index()]][i].index;
//               for (size_type j = 0; j < adj_[i2u_[neighbor]].size(); j++){
//                   if (adj_[i2u_[neighbor]][j].index == last_index){
//                       adj_[i2u_[neighbor]][j].index = n.index();
//                   }
//               }
//           }
//           return 0;
//       } else {
//           return 0;
//       }
//   };
//
//   /** If edge exists: Removes edge defined by @n1 and @n2 from graph.
//    * Else: Returns 0.
//    *
//    * @pre @n1 and @n2 are valid nodes
//    * @post new_num_edges() = old_num_edges() - 1
//    *
//    * Complexity: O(num_edges)
//    */
//   size_type remove_edge(const Node& n1, const Node& n2){
//       if (has_edge(n1,n2)){
//           // Remove edge from n1 adjacency list.
//           for (size_type i = 0; i < adj_[i2u_[n1.index()]].size(); i++){
//               if (adj_[i2u_[n1.index()]][i].index == n2.index()){
//                   adj_[i2u_[n1.index()]][i] = adj_[i2u_[n1.index()]].back();
//                   adj_[i2u_[n1.index()]].pop_back();
//                   break;
//               }
//           }
//           // Remove edge from n2 adjacency list.
//           for (size_type i = 0; i < adj_[i2u_[n2.index()]].size(); i++){
//               if (adj_[i2u_[n2.index()]][i].index == n1.index()){
//                   adj_[i2u_[n2.index()]][i] = adj_[i2u_[n2.index()]].back();
//                   adj_[i2u_[n2.index()]].pop_back();
//                   break;
//               }
//           }
//           // Decrease total number of edges by 1.
//           edge_count -= 1;
//           return 0;
//       } else {
//           return 0;
//       }
//   };
//
//   /** Calls main remove_edge function (above). Same specification. */
//   size_type remove_edge(const Edge& e){
//       Node n1 = e.node1();
//       Node n2 = e.node2();
//       return remove_edge(n1,n2);
//   };
//
//   /** If edge is in graph: Removes edge at position marked by iterator
//    * @e_it
//    * Else: Returns null iterator.
//    *
//    * @pre @e_it is a valid edge iterator
//    * @post new_num_edges() = old_num_edges() - 1
//    * @return edge_iterator pointing to beginning of adjacency list
//    *
//    * Complexity: O(num_nodes())
//    */
//   edge_iterator remove_edge(edge_iterator eit){
//       Edge e = *eit;
//       if (has_edge(e.node1(),e.node2())){
//           remove_edge(e);
//       // if (has_edge(e.node1(),e.node2())){
//       //     size_type index1 = e.node1().index();
//       //     size_type index2 = e.node2().index();
//       //     for (unsigned i = 0; i < adj_[i2u_[index1]].size(); i++){
//       //         if (adj_[i2u_[index1]][i].index == index2){
//       //             adj_[i2u_[index1]][i] = adj_[i2u_[index1]].back();
//       //             adj_[i2u_[index1]].pop_back();
//       //         }
//       //     }
//       //     for (unsigned i = 0; i < adj_[i2u_[index2]].size(); i++){
//       //         if (adj_[i2u_[index2]][i].index == index1){
//       //             adj_[i2u_[index2]][i] = adj_[i2u_[index2]].back();
//       //             adj_[i2u_[index2]].pop_back();
//       //         }
//       //     }
//       //     edge_count -= 1;
//           return EdgeIterator(this, e.node1().index());
//       } else {
//           return this->edge_end();
//       }
//   };
//
//
//   //
//   // EDGES
//   //
//
//   /** @class Graph::Edge
//    * @brief Class representing the graph's edges.
//    *
//    * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
//    * are considered equal if they connect the same nodes, in either order.
//    */
//   class Edge : private totally_ordered<Edge> {
//    public:
//     /** Construct an invalid Edge. */
//     Edge() {
//     }
//
//     /** Return a node of this Edge */
//     Node node1() const {
//       return Node(edge_graph, node1_index);
//     }
//
//     /** Return the other node of this Edge */
//     Node node2() const {
//       return Node(edge_graph, node2_index);
//     }
//
//     /** Returns the value of this Edge by reference */
//     edge_value_type& value() {
//         // Value is assigned to edge (n1, n2) where n1.index() < n2.index().
//         if (node1_index < node2_index){
//             for (size_type i = 0; i < edge_graph->adj_[edge_graph->i2u_[node1_index]].size(); i++){
//                 if (edge_graph->adj_[edge_graph->i2u_[node1_index]][i].index == node2_index){
//                     return edge_graph->adj_[edge_graph->i2u_[node1_index]][i].value;
//                 }
//             }
//         } else {
//             // If n1.index() > n2.index(), then find the value associated with edge (n2, n1).
//             for (size_type i = 0; i < edge_graph->adj_[edge_graph->i2u_[node2_index]].size(); i++){
//                 if (edge_graph->adj_[edge_graph->i2u_[node2_index]][i].index == node1_index){
//                     return edge_graph->adj_[edge_graph->i2u_[node2_index]][i].value;
//                 }
//             }
//         }
//         assert(false); // Should never get here
//         return edge_graph->adj_[edge_graph->i2u_[node1_index]][0].value;
//     }
//
//     /** Returns an IMMUTABLE value of this Edge by reference
//      *
//      * Same format as above function.
//      */
//     const edge_value_type& value() const{
//         if (node1_index < node2_index){
//             for (size_type i = 0; i < edge_graph->adj_[i2u_[node1_index]].size(); i++){
//                 if (edge_graph->adj_[i2u_[node1_index]][i].index == node2_index){
//                     return edge_graph->adj_[i2u_[node1_index]][i].value;
//                 }
//             }
//         } else {
//             for (size_type i = 0; i < edge_graph->adj_[i2u_[node2_index]].size(); i++){
//                 if (edge_graph->adj_[i2u_[node2_index]][i].index == node1_index){
//                     return edge_graph->adj_[i2u_[node2_index]][i].value;
//                 }
//             }
//         }
//         assert(false); // Should never get here
//         return edge_graph->adj_[i2u_[node1_index]][0].value;
//     }
//
//     /** Returns Euclidean length of the edge */
//     double length() const {
//         Point& pos1 = edge_graph->nodes_[edge_graph->i2u_[node1_index]].position;
//         Point& pos2 = edge_graph->nodes_[edge_graph->i2u_[node2_index]].position;
//         return norm(pos1 - pos2);
//     }
//
//     /** Test whether this edge and @a e are equal.
//      *
//      * Equal edges represent the same undirected edge between two nodes.
//      */
//     bool operator==(const Edge& e) const {
//       Node node1 = Node(edge_graph, node1_index);
//       Node node2 = Node(edge_graph, node2_index);
//       return (node1 == e.node1() && node2 == e.node2());
//     }
//
//     /** Test whether this edge is less than @a e in a global order.
//      *
//      * This ordering function is useful for STL containers such as
//      * std::map<>. It need not have any interpretive meaning.
//      */
//     bool operator<(const Edge& e) const {
//       Node node1 = Node(edge_graph, node1_index);
//       Node node2 = Node(edge_graph, node2_index);
//
//       if (node1 < e.node1()){
//         return true;
//       } else if (node1 == e.node1() && node2 < e.node2()) {
//         return true;
//       } else {
//         return false;
//       }
//     }
//
//    private:
//     // Allow Graph to access Edge's private member data and functions.
//     friend class Graph;
//
//     Graph* edge_graph;
//     size_type node1_index;
//     size_type node2_index;
//
//     /** Private Constructor */
//     Edge(const Graph* graph, const size_type node1, const size_type node2)
//         : edge_graph(const_cast<Graph*>(graph)), node1_index(node1),
//         node2_index(node2) {}
//   };
//
//   /** Return the total number of edges in the graph.
//    *
//    * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
//    */
//   size_type num_edges() const {
//     return edge_count;
//   }
//
//   /** Return the edge with index @a i.
//    * @pre 0 <= @a i < num_edges()
//    *
//    * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
//    */
//   Edge edge(size_type i) const {
//     EdgeIterator eit = this->edge_begin();
//     std::advance(eit,i);
//     return (*eit);
//   }
//
//   /** Test whether two nodes are connected by an edge.
//    * @pre @a a and @a b are valid nodes of this graph
//    * @return True if for some @a i, edge(@a i) connects @a a and @a b.
//    *
//    * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
//    */
//   bool has_edge(const Node& a, const Node& b) const {
//       for (size_type i = 0; i < adj_[i2u_[a.index()]].size(); i++){
//           if (adj_[i2u_[a.index()]][i].index == b.index()){
//               return true;
//           }
//       }
//       return false;
//   }
//
//   /** Add an edge to the graph, or return the current edge if it already exists.
//    * @pre @a a and @a b are distinct valid nodes of this graph
//    * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
//    * @post has_edge(@a a, @a b) == true
//    * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
//    *       Else,                        new num_edges() == old num_edges() + 1.
//    *
//    * Can invalidate edge indexes -- in other words, old edge(@a i) might not
//    * equal new edge(@a i). Must not invalidate outstanding Edge objects.
//    *
//    * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
//    */
//   Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value = edge_value_type()) {
//     if (this->has_edge(a,b)){
//         return Edge(this, a.index(), b.index());
//     } else {
//         internal_edge newedge1, newedge2;
//         newedge1.index = b.index(); newedge2.index = a.index();
//         newedge1.value = edge_value; newedge2.value = edge_value;
//
//         // Adding pair corresponding to newedge to @a's and @b's adjacency list
//         // std::pair<size_type, internal_edge&> a_pair (b.index(), newedge);
//         // assert(a_pair.first == b.index());
//         // assert(a_pair.second.value == newedge.value);
//         // std::pair<size_type, internal_edge&> b_pair (a.index(), newedge);
//         adj_[i2u_[a.index()]].push_back(newedge1);
//         adj_[i2u_[b.index()]].push_back(newedge2);
//
//         // Increasing edge count by 1
//         edge_count += 1;
//
//         return Edge(this, a.index(), b.index());
//         }
//     }
//
//   /** Remove all nodes and edges from this graph.
//    * @post num_nodes() == 0 && num_edges() == 0
//    *
//    * Invalidates all outstanding Node and Edge objects.
//    */
//   void clear() {
//     edge_count = 0;
//     nodes_.clear();
//     i2u_.clear();
//     adj_.clear();
//   }
//
//   //
//   // Node Iterator
//   //
//
//   /** @class Graph::NodeIterator
//    * @brief Iterator class for nodes. A forward iterator.
//    */
//   class NodeIterator {
//    public:
//     // These type definitions let us use STL's iterator_traits.
//     using value_type        = Node;                     // Element type
//     using pointer           = Node*;                    // Pointers to elements
//     using reference         = Node&;                    // Reference to elements
//     using difference_type   = std::ptrdiff_t;           // Signed difference
//     using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
//
//     /** Construct an invalid NodeIterator. */
//     NodeIterator() {
//     }
//
//     /** Dereferencing operation. Dereferenced NodeIterator returns a Node. */
//     value_type operator*() const {
//         return Node(niterator_graph, id);
//     };
//
//     /** Increment operation. Incrementing a NodeIterator should return a
//      * NodeIterator&.
//      */
//     node_iterator& operator++() {
//         ++id;
//         return *this;
//     }
//
//     /** Checks if this node_iterator is equal to input. */
//     bool operator==(const node_iterator& nit) const {
//         return (niterator_graph == nit.niterator_graph && id == nit.id);
//     };
//
//     /** Checks if this node_iterator is NOT equal to input. */
//     bool operator!=(const node_iterator& nit) const {
//         return !(niterator_graph == nit.niterator_graph && id == nit.id);
//     };
//
//    private:
//     friend class Graph;
//
//     Graph* niterator_graph;
//     size_type id;
//
//     /** Private Constructor */
//     NodeIterator(const Graph* graph, size_type index)
//         : niterator_graph(const_cast<Graph*>(graph)), id(index) {}
//   };
//
//   /** Points to the first element in the collection of nodes. */
//   node_iterator node_begin() const{
//       return NodeIterator(this, 0);
//   };
//
//   /** Points to one past the last element in the collection of nodes. */
//   node_iterator node_end() const{
//       return NodeIterator(this, i2u_.size());
//   };
//
//   //
//   // Incident Iterator
//   //
//
//   /** @class Graph::IncidentIterator
//    * @brief Iterator class for edges incident to a node. A forward iterator.
//    */
//   class IncidentIterator {
//    public:
//     // These type definitions let us use STL's iterator_traits.
//     using value_type        = Edge;                     // Element type
//     using pointer           = Edge*;                    // Pointers to elements
//     using reference         = Edge&;                    // Reference to elements
//     using difference_type   = std::ptrdiff_t;           // Signed difference
//     using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
//
//     /** Construct an invalid IncidentIterator. */
//     IncidentIterator() {
//     }
//
//     /** Dereferencing operation. Dereferenced IncidentIterator returns an Edge.
//      * The Edge returned has node_1 = current node, and node_2 =
//      * neighbor.
//      */
//     value_type operator*() const {
//         size_type node2_index = iiterator_graph->adj_[iiterator_graph->i2u_[node_index]][id].index;
//         return Edge(iiterator_graph, node_index, node2_index);
//     };
//
//     /** Incrementing the incident edge. Returns IncidentIterator. */
//     IncidentIterator& operator++() {
//         ++id;
//         return *this;
//     };
//
//     /** Checks if current iterator is equal to input iterator. Returns a
//      * boolean value.
//      */
//     bool operator==(const IncidentIterator& iit) const {
//         return (iiterator_graph == iit.iiterator_graph && node_index == iit.node_index && id == iit.id);
//     };
//
//     /** Checks if this incident_iterator is NOT equal to input. */
//     bool operator!=(const IncidentIterator& iit) const {
//         return !(iiterator_graph == iit.iiterator_graph && node_index == iit.node_index && id == iit.id);
//     };
//
//    private:
//     friend class Graph;
//     Graph* iiterator_graph;
//     size_type node_index;
//     size_type id;
//
//     /** Private Constructor */
//     IncidentIterator(const Graph* graph, size_type node, size_type index)
//         : iiterator_graph(const_cast<Graph*>(graph)), node_index(node), id(index) {}
//   };
//
//   //
//   // Edge Iterator
//   //
//
//   /** @class Graph::EdgeIterator
//    * @brief Iterator class for edges. A forward iterator. */
//   class EdgeIterator {
//    public:
//     // These type definitions let us use STL's iterator_traits.
//     using value_type        = Edge;                     // Element type
//     using pointer           = Edge*;                    // Pointers to elements
//     using reference         = Edge&;                    // Reference to elements
//     using difference_type   = std::ptrdiff_t;           // Signed difference
//     using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
//
//     /** Construct an invalid EdgeIterator. */
//     EdgeIterator() {
//     }
//
//     /** Dereferencing operator returns an Edge. */
//     value_type operator*() const {
//         // size_type node1_index = eiterator_graph->nodes_[eiterator_graph->i2u_[node_id]].index;
//         size_type node2_index = eiterator_graph->adj_[eiterator_graph->i2u_[node_id]][incident_id].index;
//         return Edge(eiterator_graph, node_id, node2_index);
//     };
//
//     /** Increment operator returns an EdgeIterator. */
//     EdgeIterator& operator++() {
//         if (node_id == eiterator_graph->i2u_.size()){
//             return *this;
//         }
//         ++incident_id;
//         this->fix();
//         if (node_id == eiterator_graph->i2u_.size()){
//             return *this;
//         }
//         while(node_id >= eiterator_graph->adj_[eiterator_graph->i2u_[node_id]][incident_id].index){
//             ++incident_id;
//             this->fix();
//             if (node_id >= eiterator_graph->i2u_.size()){
//                 return *this;
//             }
//         }
//         return *this;
//     };
//
//     /** Checks if input EdgeIterator is equal to current iterator. Graphs and ids must be equal.
//      */
//     bool operator==(const EdgeIterator& eit) const {
//         return (eiterator_graph == eit.eiterator_graph && node_id == eit.node_id && incident_id == eit.incident_id);
//     };
//
//     /** Checks if this edge_iterator is NOT equal to input. */
//     bool operator!=(const EdgeIterator& eit) const {
//         return !(eiterator_graph == eit.eiterator_graph && node_id == eit.node_id && incident_id == eit.incident_id);
//     };
//
//     // Helper function
//     void fix() {
//         if (node_id == eiterator_graph->i2u_.size()){
//             return;
//         }
//         while(node_id >= eiterator_graph->adj_[eiterator_graph->i2u_[node_id]][incident_id].index){
//             if(incident_id >= eiterator_graph->adj_[eiterator_graph->i2u_[node_id]].size()){
//                 incident_id = 0;
//                 ++node_id;
//                 // Handles case where we reach the end of i2u_
//                 if (node_id == eiterator_graph->i2u_.size()) {
//                     return;
//                 }
//             } else {
//                 ++incident_id;
//             }
//         }
//
//         // Handles the case where incident-edge increment reaches the end of the incident vector
//         while (incident_id >= eiterator_graph->adj_[eiterator_graph->i2u_[node_id]].size()){
//             incident_id = 0;
//             ++node_id;
//             // Handles case where we reach the end of i2u_
//             if (node_id >= eiterator_graph->i2u_.size()) {
//                 return;
//             }
//             if (eiterator_graph->adj_[eiterator_graph->i2u_[node_id]].size() != 0){
//                 while(node_id >= eiterator_graph->adj_[eiterator_graph->i2u_[node_id]][incident_id].index){
//                     ++incident_id;
//                     if (incident_id == eiterator_graph->adj_[eiterator_graph->i2u_[node_id]].size()){
//                         break;
//                     }
//                 }
//             }
//         }
//         // assert(node_id <= eiterator_graph->i2u_.size());
//         return;
//     };
//
//    private:
//     friend class Graph;
//     Graph* eiterator_graph;
//     size_type node_id;
//     size_type incident_id;
//
//     /** Private Constructor */
//     EdgeIterator(const Graph* graph, size_type node)
//         : eiterator_graph(const_cast<Graph*>(graph)), node_id(node) {
//             incident_id = 0;
//             this->fix();
//         }
//   };
//
//   // Points to the beginning of the adjacency list
//   edge_iterator edge_begin() const {
//       return EdgeIterator(this, 0);
//   };
//
//   // Points to one-past the end of the adjacency list.
//   edge_iterator edge_end() const {
//       return EdgeIterator(this, i2u_.size());
//   };
//
//  private:
//     struct internal_edge {
//         size_type index;
//         edge_value_type value;
//     };
//
//     struct internal_node {
//         Point position;
//         node_value_type value;
//         size_type index;
//     };
//
//     size_type edge_count;
// };
//
// #endif // CME212_GRAPH_HPP
