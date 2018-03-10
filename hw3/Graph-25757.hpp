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
  typedef V node_value_type;
  typedef E edge_value_type;
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
 struct internal_node;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), adj_(), mapv()
    {
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
  class Node: private totally_ordered<Node>{
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
      return graph_->nodes_[proxy_uid_].point_;
    }
    /** Return this node's position (modifiable). */
    Point& position(){
      // HW0: YOUR CODE HERE
      return graph_->nodes_[proxy_uid_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[proxy_uid_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    size_type&idx(){
      return graph_->nodes_[proxy_uid_].idx_;
    }
    /** Return a reference this node's value, which is modifiable. */

    node_value_type& value(){
      return graph_->nodes_[proxy_uid_].v_;
    }
    /** Return this node's value. cannot modify it.*/
    const node_value_type& value() const{
      return graph_->nodes_[proxy_uid_].v_;
    }
    /** Return this node's degree. */
    size_type degree() const{
      return graph_->adj_[index()].size();
    }
    /** Return a beginning incident iterator for the node. */
    IncidentIterator edge_begin() const{
      return Graph<V, E>::IncidentIterator(graph_,index(),0);
    }
    /** Return a ending incident iterator for the node. */
    IncidentIterator edge_end() const{
      return Graph<V, E>::IncidentIterator(graph_,index(),degree());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ( graph_ == n.graph_ && n.position() == position() &&
            n.index() == proxy_uid_ && value() == n.value()){
        return true;
      }else{
        return false;
      }
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
      if(graph_ < n.graph_){
        return true;
      }
      if(graph_ > n.graph_){
        return false;
      }
      if (index() < n.index()){
        return true;
      }else{
        return false;
      }
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type proxy_uid_;

    Node (const Graph* graph, size_type proxy_uid)
        : graph_(const_cast<Graph*>(graph)), proxy_uid_(proxy_uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return mapv.size();
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
  Node add_node (const Point& position,
            const node_value_type& v = node_value_type()){
    // HW0: YOUR CODE HERE
    size_type its_id = nodes_.size();
    mapv.push_back(its_id);
    nodes_.push_back(internal_node{num_nodes()-1, position, v});
    std::vector<size_type> empty;
    adj_.push_back(empty);
    std::vector<internal_edge> empty_edge;
    edges_.push_back(empty_edge);
    return Node(this, its_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.graph_ == this && n.position() == node(n.index()).position() &&
      n.value() == node(n.index()).value()){
      return true;
    }
    return false;
  }
/*
  void printmapv(unsigned n){
      std::cout << mapv[n] << std::endl;
  }
*/
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const{
    // HW0: YOUR CODE HERE
    return Node(this, mapv[i]);
  }

  void check_invariant(){
    for(size_type i = 0; i < mapv.size(); ++i){
      //std::cout << "i is " << i << std::endl;
      //std::cout << "idx is " << node(i).idx() << std::endl;
      assert(node(i).idx() == i);
    }
  }

  /** Remove an node to the graph, return its current index.
   * @pre @a n has to be in this graph
   * @return the current index of the node @a n
   * @post all iterators which points to (n, end) are invalidated
   * @post @a n in mapv is removed
   * @post all edges connected to @a n in adj_ are removed
   *
   * Complexity: No more than O(num_nodes), assuming a sparse graph
   */
  size_type remove_node(const Node & n){
  //  assert(n.index() < num_nodes());

    for(unsigned it = 0; it < adj_[n.index()].size(); ++it){
      for(auto it_ = adj_[Node(this, adj_[n.index()][it]).index()].begin();
      it_ != adj_[Node(this, adj_[n.index()][it]).index()].end(); ++it_){
        if((*it_) == mapv[n.index()]){
          adj_[Node(this, adj_[n.index()][it]).index()].erase(it_);
          break;
        }
      }
    }
    adj_.erase(adj_.begin() + n.index());


    for(unsigned i = n.index() + 1 ; i < num_nodes(); i++){
      node(i).idx() = node(i).idx() - 1;
    }

    mapv.erase(mapv.begin() + n.index());
    return n.index();

  }

  /** Remove an node to the graph, return its updated current iterator.
   * @pre @a n_it has to be a valid node iterator
   * @return the updated current iterator
   * @post all iterators which points to (n, end) are invalidated
   * @post the node where @a n_it in points to in mapv is removed
   * @post all edges connected to he node where @a n_it in points to
   * in adj_ are removed
   * Complexity: No more than O(num_nodes), assuming a sparse graph
   */
  node_iterator remove_node(node_iterator n_it){
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, uid_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, uid_2);
    }


    /** calculate the length of the edge*/
    double length() const{
      return norm_2(node1().position() - node2().position());
    }

    edge_value_type& value(){
      //std::cout << <<std::endl;
      for(size_type i = 0; i < graph_->edges_[uid_1].size(); i++){
        if(graph_->edges_[uid_1][i].node2 == uid_2){
          return graph_->edges_[uid_1][i].v_;
        }
      }
      assert(false);
    }

    const edge_value_type& value() const{
      for(size_type i = 0; i < graph_->edges_[uid_1].size(); i++){
        if(graph_->edges_[uid_1][i].node2 == uid_2){
          return graph_->edges_[uid_1][i].v_;
        }
      }
      assert(false);
/*      for(size_type i = 0; i < graph_->edges_[graph_->mapv[node1().index()]].size(); i++){
        if(graph_->edges_[graph_->mapv[node1().index()]][i].node2 == graph_->mapv[uid_2]){
          return graph_->edges_[graph_->mapv[node1().index()]][i].v_;
        }
      }
      std::cout << "error edge value" << std::endl;
      assert(false);
      */
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_ == e.graph_ && e.node1() == node1() && e.node2() == node2()){
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_ < e.graph_){
        return true;
      }
      if(graph_ > e.graph_){
        return false;
      }
      if(node1().index() < e.node1().index()){
        return true;
      }else if(node1().index() > e.node1().index()){
        return false;
      }else{
        if(node2().index() < e.node2().index()){
          return true;
        }
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_1, uid_2;

    Edge (const Graph* graph, size_type edge_uid_1, size_type edge_uid_2)
    : graph_(const_cast<Graph*>(graph)), uid_1(edge_uid_1), uid_2(edge_uid_2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type count = 0;
    for(auto it = adj_.begin(); it != adj_.end(); ++it){
      count = count + (*it).size();
    }
    return count/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type k) const {
    // HW0: YOUR CODE HERE
    size_type count = 0;
    for(size_type i = 0; i < adj_.size(); i++){
      for(size_type j = 0; j < adj_[i].size(); j++){
        if(node(adj_[i][j]) < node(i)){
          count++;
        }
        if(count == k){
          return Edge(this, mapv[i], adj_[i][j]);
        }
      }
    }
    std::cout << "edge index error" << std::endl;
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
    for(size_type i = 0; i < adj_[a.index()].size(); i++){
      if(adj_[a.index()][i] == mapv[b.index()]){
        //std::cout << a.index() << ", "<< mapv[b.index()] << std::endl;
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
  Edge add_edge(const Node& a, const Node& b,
            const edge_value_type& v = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_edge(a, b)){
      return Edge(this, mapv[a.index()], mapv[b.index()]);
    }
    if (has_edge(b, a)){
      return Edge(this, mapv[b.index()], mapv[a.index()]);
    }
    std::vector<size_type> add {a.index(), b.index()};
    adj_[a.index()].push_back(mapv[b.index()]);
    adj_[b.index()].push_back(mapv[a.index()]);
    edges_[mapv[a.index()]].push_back(internal_edge{mapv[a.index()],mapv[b.index()],v});
    edges_[mapv[b.index()]].push_back(internal_edge{mapv[b.index()],mapv[a.index()],v});
    return Edge(this, mapv[a.index()], mapv[b.index()]);
  }

  /*Get the other directed edge in edges_*/
  Edge get_counterpart(const Edge& e){
    return Edge(this, mapv[e.node2().index()], mapv[e.node1().index()]);
  }


  /** Remove an edge to the graph, return true if sucessfuly removed it.
   * @pre edge connecting @a a, @a b has to point to this graph
   * @return 1 if removed it, 0 if doesn't have this edge.
   * @post all iterators which points to (current, end) are invalidated
   * @post edge connecting @a a, @a b in adj_ is removed
   *
   * Complexity: No more than O(num_edges), assuming a sparse graph
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(has_edge(a, b)){
      for(auto it = adj_[a.index()].begin(); it != adj_[a.index()].end(); ++it){
        if((*it) == mapv[b.index()]){
          adj_[a.index()].erase(it);
          break;
        }
      }
      for(auto it = adj_[b.index()].begin(); it != adj_[b.index()].end(); ++it){
        if((*it) == mapv[a.index()]){
          adj_[b.index()].erase(it);
          break;
        }
      }
      return 1;
    }
    return 0;
  }

  /** Remove an edge to the graph, return true if sucessfuly removed it.
   * @pre edge @a e has to point to this graph
   * @return 1 if removed it, 0 if doesn't have this edge.
   * @post all iterators which points to (iterator to @a e, end) are invalidated
   * @post edge @a e in adj_ is removed
   *
   * Complexity: No more than O(num_edges), assuming a sparse graph
   */

  size_type remove_edge(const Edge &e){
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge to the graph, return the updated current edge iterator.
   * @pre edge that @a e_it points to has to point to this graph
   * @return updated current edge iterator @a e_it.
   * @post all iterators which points to (@a e_it, end) are invalidated
   * @post edge @a e_it points to in adj_ is removed
   *
   * Complexity: No more than O(num_edges), assuming a sparse graph
   */

  edge_iterator remove_edge ( edge_iterator e_it ){
    remove_edge(*e_it);
    return e_it;

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
    adj_.clear();
    mapv.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    NodeIterator() {
    }
    /**return the node we want to work on
    * @pre index < num_nodes()
    */
    value_type operator*() const {
      return graph_->node(index);
    }
    /**incrementing the index so can iterate
    * @pre index < num_nodes()
    */
    NodeIterator & operator ++(){
      index++;
    return *this;
    }
    /**overloading the == operator for comparison*/
    bool operator ==( const NodeIterator & ni) const{

      if(graph_ == ni.graph_ && index == ni.index){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type index;
    NodeIterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index(index){
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**create a beginning node iterator for the graph*/
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }
  /**create a ending node iterator for the graph*/
  node_iterator node_end() const{
    return node_iterator(this, num_nodes());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>{
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

    /**accessing the corresponding edge for the incident iterator
    * @pre iterate_index < node(spawn_index).degree()
    */
    Edge operator *() const{
      return  Graph<V, E>::Edge(graph_, graph_->mapv[spawn_index],
        graph_->adj_[spawn_index][iterate_index]);
    }

    /**incrementing the inside index
    * @pre iterate_index < node(spawn_index).degree()
    */
    incident_iterator & operator ++(){
      iterate_index++;
      return *this;
    }

    /**overloading the == iterator for comparison*/
    bool operator ==( const incident_iterator & iit ) const{
      if (iit.graph_ == graph_ && iit.spawn_index == spawn_index &&
      iit.iterate_index == iterate_index){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type spawn_index;
    size_type iterate_index;
    IncidentIterator(const Graph* graph, size_type s_index, size_type i_index)
        : graph_(const_cast<Graph*>(graph)), spawn_index(s_index),
              iterate_index(i_index){
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**accessing the corresponding edge for the iterator
    * @pre index < num_edges()
    */
    Edge operator*() const{
      return Edge(graph_, graph_->mapv[index1], graph_->adj_[index1][index2]);
    }

    /**incrementing the index for the iterator
    * @pre index < num_edges()
    */
    EdgeIterator& operator++(){
      increment();
      while(!(index1 == graph_->adj_.size() && index2 == 0) && pred()){
        increment();
      }
      return *this;
    }

    bool pred(){
      if(graph_->mapv[index1] > graph_->adj_[index1][index2]){
          return false;
        }
      return true;
    }

    void increment(){
      if(index2 < graph_->adj_[index1].size() - 1){
        index2++;
      }else{
        index1++;
        index2 = 0;
      }
    }

    void fix(){
      if(index1 == graph_->adj_.size() && index2 == 0){

      }else if(pred()){
        increment();
        while(pred() && !(index1 == graph_->adj_.size() && index2 == 0)){
          increment();
        }
      }
    }

    /**overloading the == operator for comparison*/
    bool operator==(const EdgeIterator& ei) const{
      if(graph_ == ei.graph_ && index1 == ei.index1 && index2 == ei.index2){
        return true;
      }
      return false;
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type index1;
    size_type index2;
    EdgeIterator(const Graph* graph, size_type index1, size_type index2)
        : graph_(const_cast<Graph*>(graph)), index1(index1), index2(index2){
          fix();
    }
  };
/*
  void printedge(){
    std::cout << "Printing edges" << std::endl;
    for(size_type i = 0; i < edges_.size(); ++i){
      for(size_type j = 0; j < edges_[i].size(); ++j){
        std::cout << edges_[i][j].node1 << ", " << edges_[i][j].node2 << "  ";
      }
      std::cout << std::endl;
    }
  }


  void printadj(){
    std::cout << "Printing adj list" << std::endl;
    for(size_type i = 0; i < adj_.size(); ++i){
      for(size_type j = 0; j < adj_[i].size(); ++j){
        std::cout << i << ", " << adj_[i][j] << "  ";
      }
      std::cout << std::endl;
    }
  }

  void printmapv(){
    std::cout << "Printing mapping vector" << std::endl;
    for(size_type i = 0; i < mapv.size(); ++i){
        std::cout << i << ", " << mapv[i] << ", " << nodes_[mapv[i]].idx_<< std::endl;
    }
  }
*/
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**create a beginning edge iterator for the graph*/
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  /**create a ending iterator for the graph*/
  edge_iterator edge_end() const{
    return EdgeIterator(this, adj_.size(), 0);
  }




 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
    size_type idx_;
    Point point_;
    node_value_type v_;
  };

  struct internal_edge{
    size_type node1;
    size_type node2;
    edge_value_type v_;
  };

  std::vector<internal_node> nodes_;
  std::vector<std::vector<internal_edge>> edges_;
  std::vector<std::vector<size_type>> adj_;
  std::vector<size_type> mapv;
};


/*
std::cout << "Printing adj list" << std::endl;
for(auto i = 0; i < graph_->adj_.size(); ++i){
  for(auto j = 0; j < graph_->adj_[i].size(); ++j){
    std::cout << i << ", " << graph_->adj_[i][j] << "  ";
  }
  std::cout << std::endl;
}
*/
#endif // CME212_GRAPH_HPP
