#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
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
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type  = V;
  using edge_value_type  = E;

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
  Graph()
    : nodes_(), adj_() {
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
  class Node : private totally_ordered<Node>{
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

    /** Get the position of a node using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an r-value that cannot be modified
     */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().p_;
    }

    // HW2: YOUR CODE HERE
    /** Get the position of a node using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an l-value that can be assigned to
     */
    Point& position(){
      return fetch().p_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    /** Get the value of a node using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an l-value that can be assigned to
     */
    node_value_type& value(){
      return fetch().value_;
    }

    /** Get the value of a node using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an r-value that cannot be modified
     */
    const node_value_type& value() const{
      return fetch().value_;
    }

    /** Get the degree of a given node. That is, the number of incident edges
     *  to it
     */
    size_type degree() const{
      return g_->adj_.at(uid_).size();
    }

    /** Begin iterator function for edges incident to the current node
     *
     * @return An edge_iterator pointing to the first edge incident to this
     *         node.
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(g_, uid_, 0);
    }

    /** End iterator function for edges incident to the current node
     *
     * @return An incident_iterator pointing one past the last edge incident
     *         to this node.
     */
    incident_iterator edge_end() const{
      return IncidentIterator(g_, uid_, g_->adj_.at(uid_).size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return n.index() == uid_ && n.g_ == g_;
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
      return uid_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* g_;
    size_type uid_;

    /** Private constructor to create valid Nodes. Accessible to friends of
     * this class.
     */
    Node(const Graph* graph, size_type uid)
        : g_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    /** Helper method to return the appropriate internal node.
     * This accesses the underlying data structure in which the Graph is stored.
     */
    internal_node& fetch() const {
      return g_->nodes_.at(uid_);
    }

    bool valid() const {
      return uid_ >= 0 && uid_ < g_->nodes_.size() &&
             g_->nodes_[uid_].idx_ < g_->i2u_.size() &&
             g_->i2u_[g_->nodes_[uid_].idx_] == uid_;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u_.size();
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
  Node add_node(const Point& position,
    const node_value_type& value = node_value_type()){
    // HW0: YOUR CODE HERE
    internal_node n = internal_node(i2u_.size(), position, value);
    i2u_.push_back(nodes_.size());
    nodes_.push_back(n);
    std::vector<internal_edge> v {};
    adj_.push_back(v);
    return Node(this, nodes_.size()-1); // A valid Node Object
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return i2u_.at(nodes_.at(n.index()).idx_) == n.index();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i2u_.at(i));
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return g_->node(g_->nodes_[uid1_].idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return g_->node(g_->nodes_[uid2_].idx_);
    }

    // HW2: YOUR CODE HERE
    /** Get the value of an edge using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an l-value that can be assigned to
     */
    edge_value_type& value(){
      return fetch().value_;
    }

    /** Get the value of an edge using the proxy pattern to access the
     *  underlying graph representation
     *
     * @return an r-value that cannot be modified
     */
    const edge_value_type& value() const{
      return fetch().value_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return e.g_ == g_ &&
             e.node1().index() == uid1_ &&
             e.node2().index() == uid2_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return uid1_ < e.node1().index() ||
             (uid1_ == e.node1().index() &&
              uid2_ < e.node2().index());
    }

    /** Calculates the length of the edge by looking at the norm2 node of the
     * difference between its endpoints.
     */
    double length() const{
      return norm(node1().position() - node2().position());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* g_;
    size_type uid1_, uid2_;

    /** Private constructor to create valid Edges. Accessible to friends of
     * this class.
     */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
        : g_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }

    /** Helper method to return the appropriate internal edge.
     * This accesses the underlying data structure in which the Graph is stored.
     *
     * @return The internal_edge coresponding to the current Edge
     *
     * Complexity: O(max_node_degree). This only loops over one sub-vector of
     *             the underlying adjacency list
     */
    internal_edge& fetch() const {
      for(size_type i = 0; i<g_->adj_.at(uid1_).size(); i++){
        if(g_->adj_.at(uid1_).at(i).uid2_ == uid2_){
          return g_->adj_.at(uid1_).at(i);
        }
      }
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type m = 0;
    for(size_type i=0; i<i2u_.size(); i++){
      m += adj_.at(i2u_[i]).size();
    }
    return m/2; // If adding bi-directional edges on undirected graphs: m/2
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type k = i;
    if(i < num_edges()){
      for(size_type j=0; j<adj_.size(); j++){
        if(k < adj_.at(j).size()){
          return Edge(this, adj_.at(j).at(k).uid1_, adj_.at(j).at(k).uid2_);
        }
        k -= adj_.at(j).size();
      }
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for(size_type i=0; i<adj_.at(a.index()).size(); i++){
      if(adj_.at(a.index()).at(i).uid2_ == b.index()){
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
    if(!has_edge(a,b) && !has_edge(b,a)){
      // Edge (a,b)
      internal_edge e1 = internal_edge(a.index(), b.index());
      adj_.at(a.index()).push_back(e1);
      // Edge (b,a)
      internal_edge e2 = internal_edge(b.index(), a.index());
      adj_.at(b.index()).push_back(e2);

      return Edge(this, a.index(), b.index());
    }
    return Edge(); // Invalid Edge
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
    i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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

    /** Dereference operator for the NodeIterator.
     *
     * @return The Node the iterator is currently pointing to.
     */
    Node operator*() const{
      return Node(g_, g_->i2u_.at(uid_));
    }

    /** Increment operator for the NodeIterator
     *
     * @return An iterator pointing to the next Node in the graph
     */
    node_iterator& operator++(){
      uid_++;
      return *this;
    }

    /** Equallity operator for the NodeIterator
     *
     * Formally, two iterators are the same if
     * - they point to the same Graph (g_) and
     * - they point to the same index (uid_) in the nodes_ vector.
     *
     * @param[in] it The iterator against which we wish to compare the current
     *            iterator
     *
     * @return A boolean value indicating if this iterator is equal to @a it
     */
    bool operator==(const NodeIterator& it) const{
      return g_ == it.g_ && uid_ == it.uid_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /** Private constructor for the NodeIterator. Only accessible to its
     * friends.
     */
    NodeIterator(const Graph* graph, size_type uid) :
      g_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    Graph* g_;
    size_type uid_;
  };

  // HW1 #2: YOUR CODE HERE
  /** Begin iterator function for nodes in the graph
   *
   * @return A node_iterator pointing to the first node in the graph
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** End iterator function for the nodes in the graph
   *
   * @return A node_iterator pointing one past the last node in the graph.
   */
  node_iterator node_end() const{
    return NodeIterator(this, i2u_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    /** Dereference operator for the IncidentIterator.
     *
     * @return The Edge the iterator is currently pointing to.
     */
    Edge operator*() const{
      return Edge(g_, uid1_, g_->adj_.at(uid1_).at(idx_).uid2_);
    }

    /** Increment operator for the IncidentIterator
     *
     * @return An iterator pointing to the next incident edge to the current
     *         node in the graph
     */
    incident_iterator& operator++(){
      if(idx_ < g_->adj_.at(uid1_).size()){
        idx_++;
      }
      return *this;
    }

    /** Equallity operator for the IncidentIterator
     *
     * Formally, two iterators are the same if:
     * - they point to the same Graph (g_) and
     * - they point to the same outer index (uid1_) in the adjacency list
     * - they point to the same inner index (idx_) in the adjacency list
     *
     * @param[in] it The iterator against which we wish to compare the current
     *            iterator
     *
     * @return A boolean value indicating if this iterator is equal to @a it
     */
    bool operator==(const incident_iterator& it) const{
      return g_ == it.g_ &&
             uid1_ == it.uid1_ &&
             idx_ == it.idx_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    /** Private constructor for the IncidentIterator. Only accessible to its
     * friends.
     */
    IncidentIterator(const Graph* graph, size_type uid1, size_type idx) :
      g_(const_cast<Graph*>(graph)), uid1_(uid1), idx_(idx) {
    }
    Graph* g_;
    size_type uid1_, idx_;
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
    /** Dereference operator for the EdgeIterator.
     *
     * @return The Edge the iterator is currently pointing to.
     */
    Edge operator*() const{
      size_type uid1_ = g_->i2u_[idx1_];
      if(uid1_< g_->adj_.size()){
        return Edge(g_, g_->adj_.at(uid1_).at(idx_).uid1_,
                        g_->adj_.at(uid1_).at(idx_).uid2_);
      }
      return Edge();
    }

    /** Increment operator for the EdgeIterator
     *
     * @return An iterator pointing to the next edge in the graph.
     */
    EdgeIterator& operator++(){
      size_type uid1_ = g_->i2u_[idx1_];
      if(uid1_ < g_->adj_.size()){
        idx_++;
        if(g_->adj_.at(uid1_).size() == idx_){
          idx_ = 0;
          uid1_ = g_->i2u_[++idx1_];
        }
        fix_idx();
      }
      return *this;
    }

    /** Equallity operator for the EdgeIterator
     *
     * Formally, two iterators are the same if:
     * - they point to the same Graph (g_) and
     * - they point to the same outer index (idx1_) in the adjacency list
     * - they point to the same inner index (idx_) in the adjacency list
     *
     * @param[in] it The iterator against which we wish to compare the current
     *            iterator
     *
     * @return A boolean value indicating if this iterator is equal to @a it
     */
    bool operator==(const EdgeIterator& it) const{
      return (g_ == it.g_ &&
              idx1_ == it.idx1_ &&
              idx_ == it.idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /** Private constructor for the EdgeIterator. Only accessible to its
     * friends.
     */
    EdgeIterator(const Graph* graph, size_type idx1) :
      g_(const_cast<Graph*>(graph)), idx1_(idx1), idx_(0) {fix_idx();}
    Graph* g_;
    size_type idx1_, idx_;

    /** Keeps incrementing the iterator until it points to the next valid edge
     *  in the adjacency list
     */
    // TODO: maybe clean this method a bit?
    void fix_idx(){
      size_type uid1_ = g_->i2u_[idx1_];
      while(uid1_ < g_->adj_.size()){
        if(g_->adj_.at(uid1_).size() == 0){
          idx_ = 0;
          uid1_ = g_->i2u_[++idx1_];
        }
        else{
          if(g_->adj_.at(uid1_).at(idx_).uid2_ < uid1_) break;
          idx_++;
          if(g_->adj_.at(uid1_).size() == idx_){
            idx_ = 0;
            uid1_ = g_->i2u_[++idx1_];
          }
        }

      }
    }
  };

  // HW1 #5: YOUR CODE HERE
  /** Begin iterator function for edges in the graph
   *
   * @return A edge_iterator pointing to the first edge in the graph
   */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0);
  }

  /** End iterator function for the edges in the graph
   *
   * @return An edge_iterator pointing one past the last edge in the graph.
   */
  edge_iterator edge_end() const{
    return edge_iterator(this, adj_.size());
  }

  //
  // REMOVE
  //

  /** Removes an edge from the graph
   *
   * @param[in] node1 A Node object marking one endpoint of the edge to remove
   * @param[in] node2 A Node object marking the other endpoint of the edge to
   *            remove
   *
   */
  size_type remove_edge(const Node& node1, const Node& node2){
    if(has_edge(node1, node2)){
      size_type i = 0;
      while(adj_.at(node1.index()).at(i).uid2_ != node2.index()){
        i++;
      }
      adj_.at(node1.index()).erase(adj_.at(node1.index()).begin()+i);
      i = 0;
      while(adj_.at(node2.index()).at(i).uid2_ != node1.index()){
        i++;
      }
      adj_.at(node2.index()).erase(adj_.at(node2.index()).begin()+i);
    }
    return 0; // TODO: what should I return?
  }

  /** Removes an edge from the graph
   *
   * @param[in] edge An Edge object  to remove
   *
   */
  size_type remove_edge(const Edge& edge){
    return remove_edge(edge.node1(), edge.node2());
  }

  /** Removes an edge at a position specified with an iterator
   *
   * @param[in] e_it An edge_iterator pointing to the edge to remove
   *
   */
  edge_iterator remove_edge(edge_iterator e_it){
    auto edge = *e_it;
    remove_edge(edge);
    e_it.fix_idx();
    return e_it;
  }

  /** Removes an edge at a position specified with an iterator
   *
   * @param[in] e_it An incident_iterator pointing to the edge to remove
   *
   */
  incident_iterator remove_edge(incident_iterator e_it){
    auto edge = *e_it;
    remove_edge(edge);
    return e_it;
  }

  /** Removes a node in the graph
   *
   * @param[in] node A Node object to delete from the graph
   *
   */
  size_type remove_node(const Node& node){
    if(has_node(node)){
      // First remove edges in other places in the adjacency list
      for(auto it = node.edge_begin(); it != node.edge_end(); ){
        it = remove_edge(it);
      }
      // Then remove the node
      size_type idx = nodes_[node.index()].idx_;
      std::swap(i2u_[idx], i2u_[i2u_.size()-1]);
      nodes_[i2u_.at(idx)].idx_ = idx; // Update the index
      i2u_.pop_back();
    }
    return 0; // TODO: what should I return?
  }

  /** Removes a node at a position specified with an iterator
   *
   * @param[in] n_it A node_iterator pointing to the node to remove
   *
   */
  node_iterator remove_node(node_iterator n_it){
    auto node = *n_it;
    remove_node(node);
    return n_it;
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    //size_type uid_;
    size_type idx_;
    Point p_;
    node_value_type value_;
    internal_node(size_type idx, Point point, node_value_type value = {}) :
      idx_(idx), p_(point), value_(value){}
  };

  struct internal_edge {
    size_type uid1_;
    size_type uid2_;
    edge_value_type value_;
    internal_edge(size_type uid1, size_type uid2, edge_value_type value = {}) :
      uid1_(uid1), uid2_(uid2), value_(value){}
  };

  std::vector<internal_node> nodes_;
  std::vector<size_type> i2u_;
  std::vector<std::vector<internal_edge>> adj_;

};

#endif // CME212_GRAPH_HPP
