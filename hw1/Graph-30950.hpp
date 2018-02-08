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

  /** Type of the node values. */
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
    // HW0: YOUR CODE HERE
    nb_nodes_ = 0;
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[node_idx_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value.
    *
    * Complexity: O(1).
    */
    node_value_type& value() {
      return graph_->nodes_[node_idx_].second;
    }

    /** Return this node's value without modifying it.
    *
    * Complexity: O(1).
    */
    const node_value_type& value() const {
      return graph_->nodes_[node_idx_].second;
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
      // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_){
	      return node_idx_ < n.index();
      } 
      return (graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; // Pointer back to Graph 
    size_type node_idx_; // node index 
    Node(const Graph* g, size_type i) : graph_(const_cast<Graph*>(g)), node_idx_(i) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nb_nodes_;
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
    // HW0: YOUR CODE HERE
    std::pair <Point, node_value_type> node_data (position, node_value);
    nodes_.push_back(node_data);
    std::vector<size_type> neighbors;
    adjacency_.push_back(neighbors);
    ++nb_nodes_;
    return Node(this, nb_nodes_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph_ && nb_nodes_ > n.index());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(nb_nodes_ > i);
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
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node_idx_1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node_idx_2_);
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
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
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
    // HW0: YOUR CODE HERE
    return nb_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < nb_edges_);
    return Edge(this, edges_[i].node_1, edges_[i].node_2);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type idx_1 = a.index();
    size_type idx_2 = b.index();
    std::vector<size_type> neighbors = adjacency_[idx_1];
    for(size_type i=0 ; i < neighbors.size(); ++i){
      if (neighbors[i] == idx_2){
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
    size_type idx_1 = a.index();
    size_type idx_2 = b.index();
    if(!has_edge(a, b)){
	node_pair new_edge; 
        new_edge.node_1 = idx_1;
        new_edge.node_2 = idx_2;;
        edges_.push_back(new_edge);
	adjacency_[idx_1].push_back(idx_2);
	adjacency_[idx_2].push_back(idx_1);
	nb_edges_++;
    }
    return Edge(this, idx_1, idx_2);
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
    adjacency_.clear();
    nb_nodes_ = 0;
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
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
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_idx_;

    /** Valid NodeIterator constructor. */
    NodeIterator(const Graph* g, size_type i) :
	    graph_(const_cast<Graph*>(g)), node_idx_(i) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return a node iterator pointing to the first element of the node sequence of the graph.
  * @return NodeIterator pointing to the first node.
  * @post Dereference is invalid if the node container is empty.
  *
  * Complexity: O(1).
  */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Return a node iterator point to the past-the-end element of the node sequence of the graph.
  * @return NodeIterator pointing to the past-the-end element in node sequence.
  * @post Dereference is invalid if the node container is empty.
  *
  * Complexity: O(1).
  */
  node_iterator node_end() const{
    return NodeIterator(this, nb_nodes_);
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

    /** Return the edge that the incident iterator points to.
    * @pre 0 <= incident_counter_ < graph_->adjacency_[node_idx_].size()
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      size_type incident_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_];
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
    // HW1 #3: YOUR CODE HERE
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
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return the edge that the edge iterator points to.
    * @pre 0 <= node_idx_ < graph_->num_nodes() and 0 <= incident_counter_ < graph->adjacency_[node_idx_].size()
    *
    * Complexity: O(1).
    */
    Edge operator*() const {
      size_type incident_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_];
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
    // HW1 #5: YOUR CODE HERE
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
          size_type next_node_idx_ = graph_->adjacency_[node_idx_][incident_counter_];
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

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

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
  struct node_pair{
    size_type node_1;
    size_type node_2;
  };
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.
  std::vector<std::pair<Point, node_value_type>> nodes_;
  std::vector<node_pair> edges_;
  std::vector<std::vector<size_type>> adjacency_;
  size_type nb_nodes_;
  size_type nb_edges_;
};

#endif // CME212_GRAPH_HPP
