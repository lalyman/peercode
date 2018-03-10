/*Take Peer Code 16564*/
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
template <typename V, typename E>
class Graph {
 private:

 public:

  // PUBLIC TYPE DEFINITIONS
  typedef V node_value_type; 
  typedef E edge_value_type;
  
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    nb_edges_ = 0;
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[node_idx_].point_;
    }

    /** Return this node's modifiable position. */
    Point& position() {
      return graph_->nodes_[node_idx_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_idx_;
    }

    /** Return this node's value.
    *
    * Complexity: O(1).
    */
    node_value_type& value() {
      return graph_->nodes_[node_idx_].value_;
    }

    /** Return this node's value without modifying it.
    *
    * Complexity: O(1).
    */
    const node_value_type& value() const {
      return graph_->nodes_[node_idx_].value_;
    }

    /** Return this node's degree.
    *
    * Complexity: O(1).
    */
    size_type degree() const {
      return graph_->adjacency_[node_idx_].size();
    }

    /** Return an incident edge iterator pointing to the first incident edge.
    * @return IncidentIterator pointing to the first incident edge.
    * @post Dereference is invalid if the node container is empty.
    *
    * Complexity: O(1).
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_idx_, 0);
    }

    /** Return an incident edge iterator pointing to the past-the-end incident edge in adjacency.
    * @return IncidentIterator pointing to the last incident edge.
    * @post Dereference is invalid if the node container is empty.
    *
    * Complexity: O(1).
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && node_idx_ == n.index());
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
      if (graph_ == n.graph_){
	      return node_idx_ < n.index();
      } 
      return (graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_; // Pointer back to Graph 
    size_type node_idx_; // node index 
    Node(const Graph* g, size_type i) : graph_(const_cast<Graph*>(g)), node_idx_(i) {}
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
    Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    nodes_.push_back(NodeData(position, node_value));
    // EdgeData neighbors ();
    adjacency_.push_back(std::vector<EdgeData> ());
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ && num_nodes() > n.index());
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


  /** Remove a node from the graph.
   * @param[in] n   The node to remove.
   * @return 	    1 if @a n was removed, 0 otherwise. 
   *
   * @pre           @a n is a valid Node object. 
   * @post If old graph has_node(@a node), new num_nodes() == old num_nodes() - 1
   *                                       new_num_edges() == old_num_edges() - @a n.degree() 
   *       Else,                           new num_nodes() == old num_nodes()
   *                                       new num_edges() == old num_edges()
   * @post has_node(@a n) == false 
   * @post All edges connecting @a n to any other node are removed.
   * @post All nodes with index > @a n.index() have their index decremented 
   *       by 1.
   * @post NodeIterators and EdgeIterators are invalidated.
   * @post Incident iterator associated to @a n is invalidated. 
   *
   * Complexity: O(num_nodes()) assuming the graph is sparse. 
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }
    size_type node_idx = n.index();
    // Remove edges from adjacency. 
    while (adjacency_[node_idx].size() > 0) {
      remove_edge(n, Node(this, adjacency_[node_idx][0].idx_));
    }
    // Decrement node's indexes that were greater than node_idx. 
    for (size_type i=0; i < num_nodes(); ++i) {
      for (size_type j=0; j < adjacency_[i].size(); ++j) {
        if (adjacency_[i][j].idx_ > node_idx) {
	  --adjacency_[i][j].idx_;
	}
      }
    }
    // Remove node from nodes list. 
    nodes_.erase(nodes_.begin() + node_idx);
    // Remove node from adjacency list. 
    adjacency_.erase(adjacency_.begin() + node_idx);
    return 1; 
  }

  
  /** Remove a node from the graph.
   * @param[in] n_it   A valid node iterator poiting to the node to be deleted.
   * @return 	       A node iterator pointin to the node next to the one that
   *                   was deleted. 
   *
   * @pre           @a n_it is a valid node iterator.
   * @post If old graph has_node(@a node), new num_nodes() == old num_nodes() - 1
   *                                       new_num_edges() == old_num_edges() - @a *n.degree() 
   *       Else,                           new num_nodes() == old num_nodes()
   *                                       new num_edges() == old num_edges()
   * @post has_node(@a *n) == false 
   * @post All edges connecting @a *n to any other node are removed.
   * @post All nodes with index > @a *n.index() have their index decremented 
   *       by 1.
   * @post NodeIterators and EdgeIterators are invalidated.
   * @post Incident iterator associated to @a *n is invalidated. 
   *
   * Complexity: O(num_nodes()) assuming the graph is sparse. 
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
    Edge() {}

    /** Return a node of this Edge. */
    Node node1() const {
      return graph_->node(node_idx_1_);
    }

    /** Return the other node of this Edge. */
    Node node2() const {
      return graph_->node(node_idx_2_);
    }

    /** Return the length of this Edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    } 


    
    /** Return the modifiable value of this Edge.
     * @return edge value 
     *
     * Complexity: O(d) with d the largest degree of a node.
     */
    edge_value_type& value() {
      // Swap nodes if node_idx_1 > node_idx_2 to match the logic
      // used with EdgeIterator.  
      Node node_1 = node1();
      Node node_2 = node2();
      if (node_idx_2_ < node_idx_1_) {
        node_1 = node2();
	node_2 = node1();
      }
      for (size_type i = 0; i < node_1.degree(); ++i) {
        if (graph_->adjacency_[node_1.index()][i].idx_ == node_2.index()) {
          return graph_->adjacency_[node_1.index()][i].value_;
	}
      }
    }

    /** Return the const value of this Edge. */
    const edge_value_type& value() const {
      return value();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool cond_1 = (graph_ == e.graph_);
      bool cond_2 = (node_idx_1_ == e.node_idx_1_ && node_idx_2_ == e.node_idx_2_);
      bool cond_3 = (node_idx_1_ == e.node_idx_2_ && node_idx_2_ == e.node_idx_1_);
      return (cond_1 && (cond_2 || cond_3));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      bool cond_1 = (graph_ == e.graph_);
      if (cond_1){
        if (node_idx_1_ == e.node_idx_1_){
		return (node_idx_2_ < e.node_idx_2_);
	}
	return (node_idx_1_ < e.node_idx_1_);
      }
      return (graph_ < e.graph_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_; 
    size_type node_idx_1_;
    size_type node_idx_2_;
    Edge(const Graph* g, size_type i, size_type j) : 
	    graph_(const_cast<Graph*>(g)), node_idx_1_(i), node_idx_2_(j) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return nb_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    EdgeIterator ei = edge_begin();
    while (i-- > 0) {
      ++ei;
    }
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
      if ((*ei).node2() == b) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    size_type idx_1 = a.index();
    size_type idx_2 = b.index();
    if(!has_edge(a, b)){
	adjacency_[idx_1].push_back(EdgeData(idx_2, val));
	adjacency_[idx_2].push_back(EdgeData(idx_1, val));
	nb_edges_++;
    }
    return Edge(this, idx_1, idx_2);
  }


  /** Remove an edge from the graph.
   * @param[in] a   The first node in the edge to be removed.
   * @param[in] b   The second node in the edge to be removed.
   * @return 	    1 if the edge was removed, 0 otherwise. 
   *
   * @pre  @a a and @a b are valid Node objects. 
   * @post If old graph has_edge(@a a, @a b), new num_edges() == old_num_edges() - 1
   *       Else,                              new num_edges() == old num_edges()
   * @post has_edge(@a a, @a b) == false 
   * @post Edges (@a a, @a b) and (@a b, @a @) are removed from the adjacency_. 
   * @post If old graph has_edge(@a a, @a b), @a a and @a b have their degree
   *       decremented by 1. 
   * @post EdgeIterators are invalidated.
   * @post Incident iterators associated with @a a and @a b are invalidated. 
   *
   * Complexity: O(max_degree). 
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) { 
      return 0;
    }
    for (size_type i=0; i < a.degree(); ++i) {
      if (adjacency_[a.index()][i].idx_ == b.index()) {
        adjacency_[a.index()].erase(adjacency_[a.index()].begin() + i);
	break;
      }
    }
    for (size_type i=0; i < b.degree(); ++i) {
      if (adjacency_[b.index()][i].idx_ == a.index()) {
        adjacency_[b.index()].erase(adjacency_[b.index()].begin() + i);
	break;
      }
    }
    --nb_edges_;
    return 1;
  }


  /** Remove an edge from the graph.
   * @param[in] e   The edge to be removed.
   * @return 	    1 if the edge was removed, 0 otherwise. 
   *
   * @pre  @a e is a valid edge. 
   * @post If old graph has_edge(@a e.node1(), @a e.node2()), 
   *                          new num_edges() == old_num_edges() - 1
   *       Else,                                              
   *                          new num_edges() == old num_edges()
   * @post has_edge(@a e.node1(), @a e.node2()) == false 
   * @post Edge @a e is removed from the adjacency_. 
   * @post If old graph has_edge(@a e.node1(), @a e.node2()), nodes in @a e have their 
   *       degree decremented by 1. 
   * @post EdgeIterators are invalidated.
   * @post Incident iterators associated with nodes in @a e are invalidated. 
   *
   * Complexity: O(max_degree). 
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }
  
  /** Remove an edge from the graph.
   * @param[in] e_it  Edge iterator pointin to the edge to be removed.
   * @return 	      1 if the edge was removed, 0 otherwise. 
   *
   * @pre  @a e_it is a valid edge iterator. 
   * @post If old graph has_edge(@a *e_it.node1(), @a *e_it.node2()), 
   *                     new num_edges() == old_num_edges() - 1
   *       Else,                                              
   *                     new num_edges() == old num_edges()
   * @post has_edge(@a *e_it.node1(), @a *e_it.node2()) == false 
   * @post Edge @a *e is removed from the adjacency_. 
   * @post If old graph has_edge(@a *e.node1(), @a *e.node2()), nodes in @a e 
   *       have their degree decremented by 1. 
   * @post EdgeIterators are invalidated.
   * @post Incident iterators associated with nodes in @a *e are invalidated. 
   *
   * Complexity: O(max_degree). 
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_.clear();
    nb_edges_ = 0;
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
    NodeIterator() {}

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const


    /** Return the node that the iterator point to.
    * @pre 0 <= node_idx_ < graph_->size()
    *
    * Complexity: O(1).
    */
    Node operator*() const{
      return graph_->node(node_idx_);
    }

    /** Return the node iterator pointing to the next node.
    *
    * Complexity: O(1).
    */
    node_iterator& operator++(){
      ++node_idx_;
      return *this;
    }

    /** Test whether this node iterator and @a node_iter are equal.
    *
    * Equal node iterators have the same graph and point to the same node index.
    *
    * Complexity: O(1).
    */
    bool operator==(const node_iterator& node_iter) const{
      return ((graph_ == node_iter.graph_) && (node_idx_ == node_iter.node_idx_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_idx_;

    /** Valid NodeIterator constructor. */
    NodeIterator(const Graph* g, size_type i) :
	    graph_(const_cast<Graph*>(g)), node_idx_(i) {
    }
  };

  /** Return a node iterator pointing to the first element of the node sequence of the graph.
  * @return NodeIterator pointing to the first node.
  * @post Dereference is invalid if the node container is empty.
  *
  * Complexity: O(1).
  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return a node iterator point to the past-the-end element of the node sequence of the graph.
  * @return NodeIterator pointing to the past-the-end element in node sequence.
  * @post Dereference is invalid if the node container is empty.
  *
  * Complexity: O(1).
  */
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
    IncidentIterator() {}

    /** Return the edge that the incident iterator points to.
    * @pre 0 <= incident_counter_ < graph_->adjacency_[node_idx_].size()
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      size_type incident_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_].idx_;
      return Edge(graph_, node_idx_, incident_node_idx_);
    }

    /** Return the incident iterator that points to the next incident node.
    *
    * Complexity: O(1).
    */
    incident_iterator& operator++() {
      ++incident_counter_;
      return *this;
    }

    /** Test whether this incident iterator and @a iit are equal.
    *
    * Equal node iterators have the same graph and point to the same node index.
    *
    * Complexity: O(1).
    */
    bool operator==(const incident_iterator& iit) const {
      return ((graph_ == iit.graph_) && (node_idx_ == iit.node_idx_) && (incident_counter_ == iit.incident_counter_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_idx_;
    size_type incident_counter_;

    /** Valid IncidentIterator constructor. */
    IncidentIterator(const Graph* g, size_type i, size_type j) :
	    graph_(const_cast<Graph*>(g)), node_idx_(i), incident_counter_(j) {
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
    EdgeIterator() {}

    /** Return the edge that the edge iterator points to.
    * @pre 0 <= node_idx_ < graph_->num_nodes() and 0 <= incident_counter_ < graph->adjacency_[node_idx_].size()
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      size_type incident_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_].idx_;
      return Edge(graph_, node_idx_, incident_node_idx_);
    }

    /** Return the edge iterator that points to the next edge.
    *
    * Complexity: O(1).
    */
    edge_iterator& operator++() {
      ++incident_counter_;
      get_next_valid_edge();
      return *this;
    }

    /** Test whether this edge iterator and @a eit are equal.
    *
    * Equal edge iterators have the same graph and point to the same edge.
    *
    * Complexity: O(1).
    */
    bool operator==(const edge_iterator& eit) const {
      return ((graph_ == eit.graph_) && (node_idx_ == eit.node_idx_) && (incident_counter_ == eit.incident_counter_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_idx_;
    size_type incident_counter_;

    /** Valid IncidentIterator constructor. */
    EdgeIterator(const Graph* g, size_type i) : graph_(const_cast<Graph*>(g)), node_idx_(i),  incident_counter_(0) {
	    get_next_valid_edge();
    }

    /** Helper to move the edge iterator to the next valid edge.
    * Make sure that we do not count edges twice by only considering
    * edges where the source is smaller than the incident node.
    */
    void get_next_valid_edge() {
      while (node_idx_ < graph_->adjacency_.size()) {
        while (incident_counter_ < graph_->adjacency_[node_idx_].size()) {
          size_type next_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_].idx_;
          if (node_idx_ < next_node_idx_) {
            return;
          }
          ++incident_counter_;
        }
        incident_counter_ = 0;
        ++node_idx_;
      }
      return;
    }
  };

  /** Return an edge iterator pointing to the first edge in the graph.
  * @return EdgeIterator pointing to the first edge.
  * @post Dereference is invalid if the edge container is empty.
  *
  * Complexity: O(1).
  */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return an edge iterator pointing to the past-the-end element in adjacency of the graph.
  * @return EdgeIterator pointing to the past-the-end element.
  * @post Dereference is invalid if the edge container is empty.
  *
  * Complexity: O(1).
  */
  edge_iterator edge_end() const {
    return EdgeIterator(this, adjacency_.size());
  }

 private:
  struct NodeData {
    Point point_;
    node_value_type value_;
    
    /** Constructor. */
    NodeData(Point point, node_value_type value) : point_(point), value_(value) {}
  };
  
  struct EdgeData {
    size_type idx_;
    edge_value_type value_;

    /** Constructor. */
    EdgeData(size_type idx, edge_value_type value) : idx_(idx), value_(value) {}
  };
  
  std::vector<NodeData> nodes_;
  std::vector<std::vector<EdgeData>> adjacency_;
  size_type nb_edges_;
};

#endif // CME212_GRAPH_HPP
