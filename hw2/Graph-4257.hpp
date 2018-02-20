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
      // ptr_ = std::shared_ptr<nodeinfo>();
      index_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // return ptr_->point;
      return graph_->nodes_[index_]->point;
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      // return ptr_->point;
      return graph_->nodes_[index_]->point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // return ptr_->index;
      return graph_->nodes_[index_]->index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Returns the value associated with the node
     *
     * @return the alias of the value associated with the node
     * Performs O(1) operations
     */
    node_value_type& value() {
      // return ptr_->value;
      return graph_->nodes_[index_]->value;
    }

    /** Constant function that returns the value associated with the node
     *
     * @return the alias of the value associated with the node
     * @post Cannot change the value through the alias
     * Performs O(1) operations
     */
    const node_value_type& value() const {
      // return ptr_->value;
      return graph_->nodes_[index_]->value;
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
      return n.graph_==graph_ and n.index()==index();
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
      if (n.graph_== graph_) return n.index() < index();
      else return n.graph_ < graph_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // std::shared_ptr<Graph::nodeinfo> ptr_;
    size_type index_;
    Graph* graph_;
    Node (const graph_type* g, const size_type &index){
      graph_ = const_cast<graph_type*>(g);
      // ptr_ = graph_->nodes_[index];
      index_ = index;
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
      n1_ = size_type(-1);
      offset_ = size_type(-1);
      graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, get_node2_idx());
    }

    edge_value_type& value() {
      return *((graph_->adj_[n1_][offset_]).value);
    }

    const edge_value_type& value() const {
      return *((graph_->adj_[n1_][offset_]).value);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) and
        ((n1_ == e.n1_ and offset_ == e.offset_) or
         (get_node2_idx() == e.n1_ and n1_ == e.get_node2_idx()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this == &e) return false;
      if (graph_ != e.graph_) return graph_ < e.graph_;
      if (n1_ != e.n1_) return n1_ < e.n1_;
      return offset_ < e.offset_;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& v = edge_value_type()) {
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

    std::shared_ptr<edge_value_type> v_ptr(new edge_value_type(v));

    adj_[idx1].push_back(edgeinfo(idx2, offset2, v_ptr));
    adj_[idx2].push_back(edgeinfo(idx1, offset1, v_ptr));

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

  bool remove_edge_helper(size_type idx1, size_type idx2) {
    auto it = adj_[idx1].begin();
    for (; it != adj_[idx1].end() && (*it).index != idx2; ++it) {};
    if (it == adj_[idx1].end()) return false;
    it = adj_[idx1].erase(it);
    for (; it != adj_[idx1].end(); ++it) {
      adj_[(*it).index][(*it).offset].offset -= 1;
    }
    return true;
  }

  /* Remove an edge from the graph and return a bool indicating whether succeeds.
   *
   * @param[in] n1, n2 Two end nodes of the edge
   * @return if the edge belongs to the graph, return true, else false.
   *
   * @post if succeeds: new num_edges() = old num_edges() - 1
   *       else: new num_edges() = old num_edges()
   * @post edge(e.index()) == e for all e which is not deleted from the old graph
   * @post edge(i).index() == i for all i with 0 <= i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  bool remove_edge(const Node& n1, const Node& n2) {
    size_type idx1 = n1.index();
    size_type idx2 = n2.index();
    return remove_edge_helper(idx1, idx2) and remove_edge_helper(idx2, idx1);
  }

  /* Remove an edge from the graph and return a bool indicating whether succeeds.
   *
   * @param[in] e Edge object to be removed
   * @return if the edge belongs to the graph, return true, else false.
   *
   * @post if succeeds: new num_edges() = old num_edges() - 1
   *       else: new num_edges() = old num_edges()
   * @post edge(e.index()) == e for all e which is not deleted from the old graph
   * @post edge(i).index() == i for all i with 0 <= i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  bool remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /* Remove a edge from the graph and return the iterator of the next edge in the graph.
   *
   * @param[in] e_it Edge iterator to be removed
   * @return if the edge belongs to the graph, return the next edge iterator in the graph.
   *
   * @post if succeeds: new num_edges() = old num_edges() - 1
   *       else: new num_edges() = old num_edges()
   * @post edge(e.index()) == e for all e which is not deleted from the old graph
   * @post edge(i).index() == i for all i with 0 <= i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove a node from the graph, returning a bool indicating whether the operation succeeds.
   * 
   * @param[in] node The node to be removed
   * @return bool if the node truely belongs to the current graph, return true, else false
   *
   * @post if @a node belongs to the old graph:
   *          new num_nodes() == old num_nodes() - 1
   *          new num_edges() == old num_edges() - node.degree()
   *          remove_node(node) == false, i.e. you can not delete the same node twice
   *       else:
   *          there is no change to the graph
   * @post node(i).index() == i for all i with 0 <= i < num_nodes()
   * @post node(n.index()) == n for all valid node n.
   *
   * Complexity: O(num_nodes())
   */
  bool remove_node(const Node& node) {
    if (!has_node(node)) return false;
    size_type idx = node.index();

    nodes_[idx].swap(nodes_.back());
    nodes_[idx]->index = idx;
    nodes_.pop_back();

    for (auto it_e = adj_[idx].begin(); it_e != adj_[idx].end(); ++it_e) {
      size_type idx2 = (*it_e).index;
      size_type offset2 = (*it_e).offset;

      std::iter_swap(adj_[idx2].begin()+offset2, adj_[idx2].end()-1);
      edgeinfo need_update = adj_[idx2][offset2];
      size_type idx3 = need_update.index;
      size_type offset3 = need_update.offset;
      adj_[idx3][offset3].offset = offset2;

      adj_[idx2].pop_back();
    }

    std::iter_swap(adj_.begin()+idx, adj_.end()-1);
    for (auto it_e = adj_[idx].begin(); it_e != adj_[idx].end(); ++it_e) {
      size_type idx2 = (*it_e).index;
      size_type offset2 = (*it_e).offset;

      adj_[idx2][offset2].index = idx;
    }
    adj_.pop_back();
    return true;
  }

  /** Remove a node from the graph and return the iterator of the next node in the graph.
   *
   * @param[in] it Iterator pointing at a node to be removed.
   * @return a node_iterator pointing at the next node in the graph if the remove succeeds.
   *   If the graph is empty after the operation, return end_node iterator.
   *   Otherwise, return the original iterator.
   *
   * @post The input iterator is invalidated
   * Complexity: O(num_nodes())
   */
  node_iterator remove_node(node_iterator it) {
    remove_node(*it);
    return it;
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
    std::shared_ptr<edge_value_type> value;
    friend class Graph;
    edgeinfo(size_type i, size_type o, std::shared_ptr<edge_value_type> v): 
        index(i), offset(o), value(v) {}
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
