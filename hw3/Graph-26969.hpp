/*
Some parts originally from source code in lgemar Github.
*/

#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 *  @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <cstdint>

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
  /** Template Node value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Template Edge value type. */
  using edge_value_type = E;

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

  /** Node data container */
   struct InternalNode : private equality_comparable<InternalNode>{
    size_type uid;
    size_type index;
    Point position;
    node_value_type value;

    InternalNode(size_type uid_, size_type index_, Point position_, node_value_type value_) :
      uid(uid_), index(index_), position(position_), value(value_) {}

    bool operator==(size_type i) const {
      return uid == i;
    }

   };

  /** Edge data container */
  struct InternalEdge : private equality_comparable<InternalEdge>{
    size_type n2_id;
    edge_value_type value;
    size_type index;

    InternalEdge(size_type n2_id_, edge_value_type value_, size_type index_) :
      n2_id(n2_id_), value(value_), index(index_) {}

    bool operator==(size_type i) const {
      return n2_id == i;
    }

  };

  typedef std::vector<InternalEdge> adj_list__subtype;
  typedef std::vector<adj_list__subtype> adj_list__type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : num_edges_(0) {}

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

    /** Return this node's position. This is a getter function. */
    const Point& position() const {
      return graph_->nodes_[uid_].position;
    }

    /** Return non const node position. This is a setter function. */
    Point& position() {
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].index;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     * @return true iff the nodes have the same graph and same uid
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (uid_ == n.uid_);
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return (graph_ < n.graph_ || (graph_ == n.graph_ &&  uid_ < n.uid_ ));
    }

    /** Returns this node's degree. */
    size_type degree() const {
      return graph_->adj_list_[uid_].size();
    }

    /** Returns incident iterator to beginning of this node's edges. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /** Returns iterator pointing right past the last edge for this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Return the value of this node. This is a setter function. */
    node_value_type& value() {
      return graph_->nodes_[uid_].value;
    }

    /** Return the value of this node. This is a getter function. */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Node private members
    Graph* graph_;
    size_type uid_; 

    // private constructor for Graph
    Node(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    size_type next_uid = nodes_.size();
    size_type next_idx = size();
    InternalNode internal_node{next_uid, next_idx, position, val};

    // Add data
    nodes_.push_back(internal_node);
    i2u_.push_back(next_uid);
    adj_list_.push_back(adj_list__subtype());

    return Node(this, next_uid);
  }

  /** Determine if a Node belongs to this Graph.
   * @param[in] n   Node to search.
   * @return  True if @a n is currently a Node of this Graph
   *
   * True if same graph and valid index.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @param[in] i   Node index.
   * @return  Constructed Node.
   *
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
  }

  /** Removes a node from the graph
   *
   * @param[in] n Node to be removed
   * @return 1 if old has_node(@a n), else 0
   *
   * @post new size() == old_size() - result
   * @post new num_edges() = old num_edges() - n.degree()
   *
   * complexity: O(degree)
   *
   * Can invalidate iterators
   * if old has_node(@a n), @a n becomes invalid.
   * All other nodes remain valid, because they are referenced
   * by static uids, rather than indices
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n))
      return 0;

    size_type degree = n.degree();
    size_type uid = n.uid_;
    size_type idx = n.index();

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      size_type uid2 = (*it).node2().uid_;
      auto it2 = std::find(adj_list_[uid2].begin(), adj_list_[uid2].end(), uid);
      adj_list_[uid2].erase(it2);
    }

    // update num_edges_
    num_edges_ -= degree;

    // clear adjacency list of node
    adj_list_[uid].clear();

    // copy last node to here
    i2u_[idx] = i2u_.back();
    nodes_[i2u_.back()].index = idx;
    i2u_.pop_back();

    return 1;
  }

  node_iterator remove_node(node_iterator nit) {
    remove_node(*nit);
    return nit;
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
      return graph_->node(n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2_uid_);
    }

    size_type index() const {
      size_type min = std::min(n1_uid_, n2_uid_);
      size_type max = std::max(n1_uid_, n2_uid_);
      auto it = std::find(graph_->adj_list_[min].begin(), graph_->adj_list_[min].end(), max);
      return (*it).index;
    }

    /** Return value of edge. This is a setter function. */
    edge_value_type& value() {
      size_type min = std::min(n1_uid_, n2_uid_);
      size_type max = std::max(n1_uid_, n2_uid_);
      auto it = std::find(graph_->adj_list_[min].begin(), graph_->adj_list_[min].end(), max);
      return (*it).value;
    }

    /** Return the value of this edge. This is a getter function. */
    const edge_value_type& value() const {
      size_type min = std::min(n1_uid_, n2_uid_);
      size_type max = std::max(n1_uid_, n2_uid_);
      auto it = std::find(graph_->adj_list_[min].begin(), graph_->adj_list_[min].end(), max);
      return (*it).value;
    }

    /** Return the length of this edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& e) const {
      bool nodes_same = (node1() == e.node1() && node2() == e.node2())
                        || (node1() == e.node2() && node2() == e.node1());
      bool graphs_same = graph_ == e.graph_;
      return graphs_same && nodes_same;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& e) const {
      return (graph_ < e.graph_) ||
             (graph_ == e.graph_ && n1_uid_ < e.n1_uid_) ||
             (graph_ == e.graph_ && n1_uid_ == e.n1_uid_ && n2_uid_ < e.n2_uid_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Private data members
    Graph* graph_;
    size_type n1_uid_;
    size_type n2_uid_;

    /** private constructor for Graph */
    Edge(const Graph* graph, size_type n1_uid, size_type n2_uid)
        : graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_uid_(n2_uid) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Test whether two nodes are connected by an edge.
   * @param[in] a   Node 1;
   * @param[in] b   Node 2;
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * @pre @a a and @a b are valid nodes of this graph 
   *
   * Complexity: O(deg(a)).
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto it = std::find(adj_list_[a.uid_].begin(), adj_list_[a.uid_].end(), b.uid_);
    return it != adj_list_[a.uid_].end();
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post has_edge(@a b, @a a) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(deg(a))
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // check if node already exists
    if (has_edge(a, b)) {
      return Edge(this, a.uid_, b.uid_);
    }

    size_type old_num_edges = num_edges();

    // Add to adjacency list
    InternalEdge internal_edge_to_b{b.uid_, val, old_num_edges};
    InternalEdge internal_edge_to_a{a.uid_, val, old_num_edges};

    adj_list_[a.uid_].push_back(internal_edge_to_b);
    adj_list_[b.uid_].push_back(internal_edge_to_a);
    ++num_edges_;

    return Edge(this, a.uid_, b.uid_);
  }

  /** Return edge by index. */
  Edge edge(size_type i) const {
      for (auto it = edge_begin(); it != edge_end(); ++it) {
        auto e = (*it);
        if (e.index() == i)
          return e;
      }
      return Edge();
  }

  /** Remove edge given incident nodes.
   * @param[in] a   Node 1;
   * @param[in] b   Node 2;
   * @return true if edge was removed from adj_map_, false otherwise.
   *
   * @pre 0 <= @a a.index() < num_nodes().
   * @pre 0 <= @a b.index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_edges() == old num_edges() - 1.
   * @post new @a a.degree() = old @a a.degree() - 1.
   * @post new @a b.degree() = old @a b.degree() - 1.
   *
   * Remove edge from adjacency map. Update number for edges and degrees of adjacent nodes.
   *
   * Complexity:  O(@a a.degree() + @a b.degree()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b))
      return 0;

    auto it1 = std::find(adj_list_[a.uid_].begin(), adj_list_[a.uid_].end(), b.uid_);
    adj_list_[a.uid_].erase(it1);

    auto it2 = std::find(adj_list_[b.uid_].begin(), adj_list_[b.uid_].end(), a.uid_);
    adj_list_[b.uid_].erase(it2);

    --num_edges_;

    return 1;
  }

  /** Refer to function above. */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    auto next_e_it = ++e_it;
    remove_edge(e.node1(), e.node2());
    return next_e_it;
  }

  incident_iterator remove_edge(incident_iterator e_it) {
    auto e = *e_it;
    auto next_e_it = ++e_it;
    remove_edge(e.node1(), e.node2());
    return next_e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adj_list_.clear();
    i2u_.clear();
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Dereferences iterator. Return constructed Node. */
    Node operator*() const {
      return graph_->node(idx_);
    }

    /** Increments uid and returns this */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Compare by graph and uid_ */
    bool operator==(const NodeIterator& it) const {
      return (graph_ == it.graph_) && (idx_ == it.idx_);
    }

   private:
    friend class Graph;

    /** Private constructor */
    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {}

    // Private data members
    Graph* graph_;
    size_type idx_;
  };


  /** Return iterator to node with uid 0 */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return iterator to node with uid num_nodes() */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   * */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Dereference operator.
       * @return Constructed Edge.
       * @post   (result_edge.node1().index() == node_index_) && (result_edge.node2().index() == adjacent_node_index)
       *
       * Change node indices of edge to ensure node2 is adjacent node.
       */
    Edge operator*() const {
      size_type n1_uid = n1_ind_;
      size_type n2_uid = graph_->adj_list_[n1_ind_][n2_ind_].n2_id;
      return Edge(graph_, n1_uid, n2_uid);
    }

    /** Goes to next forward edge and returns this
    * @pre 0 <= n1_ind < adj_list_.size()
    * @pre 0 <= n2_ind < adj_list_[n1_ind].size()
    *
    * Allow n1 to go one beyond for iterator.end()
    * n2 must still be in bounds
    *
    * @post (0 <= n1_ind < adj_list_.size() && 0 <= n2_ind < adj_list_[n1_ind].size()) or
    *     (n1_ind = adj_list_.size() && n2_ind == 0
    * */
    IncidentIterator& operator++() {
      ++n2_ind_;
      return *this;
    }

    /** Equality compare by graph and indices */
    bool operator==(const IncidentIterator& it) const {
      return (graph_ == it.graph_) && (n1_ind_ == it.n1_ind_) && (n2_ind_ == it.n2_ind_);
    }

    private:
      friend class Graph;
      friend class Node;

      /** private constructor */
      IncidentIterator(const Graph* graph, size_type n1_ind, size_type n2_ind)
        : graph_(const_cast<Graph*>(graph)), n1_ind_(n1_ind), n2_ind_(n2_ind) {
      }

    // Private data members
    Graph* graph_;
    size_type n1_ind_;
    size_type n2_ind_;
  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereference iterator */
    Edge operator*() const {
      size_type n1_uid = n1_ind_;
      size_type n2_uid = graph_->adj_list_[n1_ind_][n2_ind_].n2_id;
      return Edge(graph_, n1_uid, n2_uid);
    }

    /** Goes to next forward edge and returns this
     * @pre 0 <= n1_ind < adj_list_.size()
     * @pre 0 <= n2_ind < adj_list_[n1_ind].size()
     *
     * @post (0 <= n1_ind < adj_list_.size() && 0 <= n2_ind < adj_list_[n1_ind].size()) or
     *    (n1_ind = adj_list_.size() && n2_ind == 0
     * */
    EdgeIterator& operator++() {
      do {
        ++n2_ind_;
        while (n2_ind_ == graph_->adj_list_[n1_ind_].size() && n1_ind_ < graph_->adj_list_.size()) {
          n2_ind_ = 0;
          ++n1_ind_;
        }
      }
      while (n1_ind_ < graph_->adj_list_.size() && n1_ind_ <= graph_->adj_list_[n1_ind_][n2_ind_].n2_id);

      return *this;
    }

    /** Equality compare by graph and indices. */
    bool operator==(const EdgeIterator& it) const {
      return (graph_ == it.graph_) && (n1_ind_ == it.n1_ind_) && (n2_ind_ == it.n2_ind_);
    }

   private:
    friend class Graph;

    EdgeIterator(const Graph* graph, size_type n1_ind, size_type n2_ind)
        : graph_(const_cast<Graph*>(graph)), n1_ind_(n1_ind), n2_ind_(n2_ind) {}

    // Private data members
    Graph* graph_;
    size_type n1_ind_;
    size_type n2_ind_; 
  };

  /** Return iterator with uid 0 0 */
  edge_iterator edge_begin() const {
    size_type n_id = 0;
    while (adj_list_[n_id].size() == 0)
      ++n_id;

    size_type n1_ind_ = 0;
    size_type n2_ind_ = 0;
    while (n1_ind_ < adj_list_.size() && n1_ind_ <= adj_list_[n1_ind_][n2_ind_].n2_id) {
      ++n2_ind_;
      while (n2_ind_ == adj_list_[n1_ind_].size() && n1_ind_ < adj_list_.size()) {
        n2_ind_ = 0;
        ++n1_ind_;
      }
    }

    return EdgeIterator(this, n1_ind_, n2_ind_);
  }

  /** Return iterator with uid num_edges 0 */
  edge_iterator edge_end() const {
    return EdgeIterator(this, adj_list_.size(), 0);
  }


 private:
  friend class EdgeIterator;
  friend class IncidentIterator;

  // Private data members
  std::vector<InternalNode> nodes_;
  std::vector<size_type> i2u_; 
  adj_list__type adj_list_;
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP