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

  /** Predeclaration of Node type. */
  class Node;
  class nodeinfo;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  class edgeinfo;
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
      graph_ = NULL;
      ptr_ = std::shared_ptr<nodeinfo>();
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return ptr_->point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return ptr_->index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Returns the value associated with the node
     *
     * @return the alias of the value associated with the node
     * Performs O(1) operations
     */
    node_value_type& value() {
      return ptr_->value;
    }

    /** Constant function that returns the value associated with the node
     *
     * @return the alias of the value associated with the node
     * @post Cannot change the value through the alias
     * Performs O(1) operations
     */
    const node_value_type& value() const {
      return ptr_->value;
    }

    /** Returns the degree of the node
     *
     * @return the degree of the node, which equals the number of
     *   adjacent edges.
     * Performs O(1) operations
     */
    size_type degree() const{
      return graph_->adj_[index()].size();
    }

    /** Returns an incident iterator which is at the `begin` position
     *
     * @return a forward incident iterator to beginning
     * Performs O(1) operations
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, index(), 0);
    }
    
    /** Returns an incident iterator which is at the `end` position
     *
     * @return an incident iterator to end
     * @post edge_begin() <= edge_end()
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return n.graph_==this->graph_ and n.index()==this->index();
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
      if (n.graph_==this->graph_) return n.ptr_<this->ptr_;
      else return n.graph_<this->graph_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    std::shared_ptr<Graph::nodeinfo> ptr_;
    Graph* graph_;
    Node (const graph_type* g, const size_type &index){
      this->graph_ = const_cast<graph_type*>(g);
      this->ptr_ = graph_->nodes_[index];
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return this->size();
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
    size_type index = this->nodes_.size();
    this->nodes_.push_back(std::shared_ptr<nodeinfo>(new nodeinfo(position, index, value)));
    this->adj_.push_back(std::vector<edgeinfo>());
    return Node(this, index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    size_type idx=n.index();
    return this==n.graph_ and idx<this->size() and n.position()==this->node(idx).position();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
      this->n1_ = size_type(-1);
      this->offset_ = size_type(-1);
      this->graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->get_node2_idx());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->graph_ == e.graph_) and
        ((this->n1_ == e.n1_ and offset_ == e.offset_) or
         (this->get_node2_idx() == e.n1_ and this->n1_ == e.get_node2_idx()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this == &e) return false;
      if (this->graph_ != e.graph_) return this->graph_ < e.graph_;
      if (this->n1_ != e.n1_) return this->n1_ < e.n1_;
      return this->offset_ < e.offset_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type n1_,offset_;
    Graph *graph_;

    size_type get_node2_idx() const {
      return (this->graph_->adj_[n1_][offset_]).index;
    }
    Edge(const Graph* g, size_type i1, size_type o) {
      graph_ = const_cast<graph_type*>(g);
      n1_ = i1;
      offset_ = o;
    }
  };


  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type num = 0;
    for (auto it = adj_.begin(); it != adj_.end(); ++it) {
      num += (*it).size();
    }
    return num / 2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    auto it = edge_begin();
    for (size_type c = 0; c < i and it != edge_end(); ++it, ++c) {}
    return (*it);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type a_idx = a.index();
    for (auto it = adj_[a_idx].begin(); it != adj_[a_idx].end(); ++it) {
      if ((*it).index == b.index()) return true;
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
    size_type idx1 = a.index();
    size_type idx2 = b.index();

    for (auto it = adj_[idx1].begin(); it != adj_[idx1].end(); ++it) {
      if (it->index == idx2){
        //return edges_[ adj_[idx1][(*it).offset].edge_index ];
         if (a<b) return Edge(this, idx1, it->offset);
         else return Edge(this, idx2, adj_[idx2][it->offset].offset);
      }
    }

    int offset1 = adj_[idx1].size();
    int offset2 = adj_[idx2].size();

    adj_[idx1].push_back(edgeinfo(idx2, offset2));
    adj_[idx2].push_back(edgeinfo(idx1, offset1));

    // if (a<b) edges_.push_back(Edge(this, idx1, offset1));
    // else edges_.push_back(Edge(this, idx2, offset2));
    // return edges_[-1];
    if (a<b) return Edge(this, idx1, offset1);
    else return Edge(this, idx2, offset2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    adj_.clear();
    // edges_.clear();
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

    /** Get the node to which the iterator points
     *
     * @return an Node object
     * Complexity: O(1)
     */
    Node operator*() const {
      return graph_->node(current_index_);
    }

    /** Move forward the iterator
     *
     * @return an alias of node iterator
     * Complexity: O(1)
     */
    NodeIterator& operator++(){
      current_index_++;
      return (*this);
    }

    /** Compare whether two iterators are equal
     *
     * @param[in] the other node iterator to compare with
     * @return bool, whether two iterators are the same
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator& n_it) const{
      return graph_ == n_it.graph_ and current_index_ == n_it.current_index_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type *graph_;
    size_type current_index_;
    NodeIterator(const graph_type *g, size_type i){
      graph_ = const_cast<graph_type *>(g);
      current_index_ = i;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Returns an node iterator at the `begin` position
   *
   * @return a forward node iterator to the begin
   * Complexity: O(1)
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Returns an node iterator at the `end` position
   *
   * @return an node iterator to the end
   * Complexity: O(1)
   */
  node_iterator node_end() const{
    return NodeIterator(this, size());
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

    /** Get the edge `pointed` by current iterator
     *
     * @return an Edge object
     * Complexity: O(1)
     */
    Edge operator*() const {
      return Edge(graph_, index_, offset_);
    }

    /** Move forward the iterator
     *
     * @return the alias of an incident iterator
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      offset_++;
      return (*this);
    }

    /** Compare whether two iterators are equal
     *
     * @param[in] the iterator to compare with
     * @return bool, whether two iterators are equal
     */
    bool operator==(const IncidentIterator& i_it) const {
      return graph_ == i_it.graph_ and 
             index_ == i_it.index_ and 
             offset_ == i_it.offset_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type *graph_;
    size_type index_;
    size_type offset_;

    IncidentIterator(const graph_type *g, size_type n1, size_type o1){
      graph_ = const_cast<graph_type*>(g);
      index_ = n1;
      offset_ = o1;
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** Get the edge to which the iterator points
     *
     * @return an Edge object
     * Complexity: O(1)
     */
    Edge operator*() const{
      return Edge(graph_, row_, col_);
    }

    /** Move forward the iterator
     *
     * @return an alias of edge iterator
     * Complexity: O(1)
     */
    EdgeIterator& operator++(){
      col_++;
      find_next();
      return (*this);
    }

    /** Compare whether two iterators are equal
     *
     * @param[in] the other edge iterator to compare with
     * @return bool, whether two iterators are the same
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& eit) const{
        return eit.graph_ == graph_ and 
               row_ == eit.row_ and col_ == eit.col_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type row_;
    size_type col_;

    void find_next() {
      if (row_ == graph_->adj_.size()) return;
      for (;row_ < graph_->size(); ++row_) {
        for (;col_ < (graph_->adj_[row_]).size(); ++col_) {
          if (row_ < (graph_->adj_[row_][col_]).index) return;
        }
        col_ = 0;
      }
    }

    EdgeIterator(const Graph* g, size_type r, size_type c) {
      graph_ = const_cast<Graph*>(g);
      row_ = r;
      col_ = c;
      find_next();
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Returns an edge iterator at the `begin` position
   *
   * @return a forward edge iterator to the begin
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }

  /** Returns an edge iterator at the `end` position
   *
   * @return an edge iterator to the end
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, size(), 0);
  }

 
  class nodeinfo {
    Point point;
    size_type index;
    node_value_type value;
    friend class Graph;
    nodeinfo(const Point &p, size_type i, node_value_type v): 
              point(p), index(i), value(v) {}
  };

  class edgeinfo {
    size_type index; // adjacent node index
    size_type offset; // edge number for this node
    friend class Graph;
    edgeinfo(size_type i, size_type o): index(i), offset(o) {}
  };

  private:
    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    std::vector<std::shared_ptr<nodeinfo>> nodes_;
    std::vector<std::vector<edgeinfo>> adj_; //each node has a collection of edges
    // here adj_[i][j] denotes the edge information for the j th edge for node_i
    // if we define this edge as (node_i, node_i1), adj_[i][j] will stores:
    // adj_[i][j].index == i1
    // adj_[i][j].offset == k, where adj_[i1][k] refers to the same edge
    // std::vector<edge_type> edges_;
};

#endif // CME212_GRAPH_HPP
