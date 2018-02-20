#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>

// template <typename V>
class Graph {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  using edge_value_type = E;
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
  class Node: private totally_ordered<Node> {
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
      graph_ = nullptr;

    }

    /** Return this node's position. */
    const Point& position() const {

      return this->graph_->internal_nodes[this->ind_].p;
    }

    Point& position() {
      return this->graph_->internal_nodes[this->ind_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->ind_;
    }

    // Supply definitions AND SPECIFICATIONS for:
    /**
    * return reference of the value of the node.
    * @pre 0 <= @a ind_ < @a num of nodes in this graph.
    **/
    node_value_type& value(){
      return (graph_->internal_nodes)[ind_].val;

    }
    /**
    * return const reference of the value of the node.
    * @pre 0 <= @a ind_ < @a num of nodes in this graph.
    **/
    const node_value_type& value() const{
      return (graph_->internal_nodes)[ind_].val;
    }
    /**
    * return degree of the node, i.e, number of outgoing edges.
    * @pre 0 <= @a ind_ < @a num of nodes in this graph.
    **/
    size_type degree() const{
      return (graph_->adj_nodes[ind_]).size();
    }
    /**
    * return incident iterator begin();
    **/
    incident_iterator edge_begin() const{
      return incident_iterator(graph_, ind_, size_type(0));

    }
    /**
    * return incident iterator end() 
    **/
    incident_iterator edge_end() const{
      return incident_iterator(graph_, ind_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->index() == n.index()) && (this->graph_ == n.graph_);
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
      if(this->graph_ == this->graph_){
        return (this->index() < n.index());
      } else {
        return this->graph_ < this->graph_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type ind_;
    Node(const Graph* graph, size_type index):
      graph_(const_cast<Graph*>(graph)), ind_(index) {

      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->internal_nodes.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    Internal_Node new_internal_node{position, val};
    size_type new_index = this->size();
    internal_nodes.push_back(new_internal_node);
    this->adj_nodes.push_back({});
    return Node(this, new_index); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_);
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


  /** Remove the given node in the graph and its incident edges.
  * @param[in] n is a valid node in the graph.
  * return the index of the node to be removed.
  * @pre n is in the graph, has_node(n) is true
  * @pre n == this->node(n.index())
  * @post new size() == old size() - 1
  * @post new num_edges() == old num_edges - degree()
  * @post n and last_node(node with index size() - 1) is invalidated
  * @post node iterator of n and last_node is invalidated.
  * @post All outstanding incident iterators containing  last node
  *  and n are invalidated.
  * @post All outstanding edge iterators associasted with n
  * and last_node are invalidated.
  * Complexity((average degree)^2)
  */
  size_type remove_node(const Node& n){
    // replace current internal node with last node
    auto last_node = node(size() - 1);
    size_type n_index = n.index();

    // change the incident edge infos for last nodes
    size_type count = 0;
    for(auto it = last_node.edge_begin(); it != last_node.edge_end(); ++it){
      auto e = *it;
      size_type n2_index = e.node2().index();
      auto temp_val = std::move(e.value());
      // save edge value in valid position
      // we only save edge value once
      if(n_index < n2_index){
        adj_nodes[last_node.index()][count].val = std::move(temp_val);
      }

      // iterate over incident edges of last nodes,
      // update the new position and edge value to last_node

      for(size_type i = 0; i < adj_nodes[n2_index].size(); i++){
        if(adj_nodes[n2_index][i].tail_node == last_node.index()){
          adj_nodes[n2_index][i].tail_node = n_index;
          if(n2_index < n_index){
            adj_nodes[n2_index][i].val = std::move(temp_val);
          }
        }
      }
      count++;
    }

    // erase incident edges.
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto e = *it;
      size_type n2_index = e.node2().index();
      for(size_type i = 0; i < adj_nodes[n2_index].size(); i++){
        if(adj_nodes[n2_index][i].tail_node == n.index()){
          erase(adj_nodes[n2_index], i);
        }
      }
    }
    // erase node n from internal nodes vector
    erase(internal_nodes, n_index);
    // erase corresponding adj_nodes vector.
    erase(adj_nodes, n_index);
    return n_index;

  }




  /** Remove the node n pointed by the given node_iterator
  * in the graph and its incident edges.
  * @param[in] n_it is a valid node iterator in the graph.
  * return the node_iterator of the old last node(now has index of *n_it).
  * @pre n is in the graph, has_node(n) is true
  * @pre n == this->node(n.index())
  * @post new size() == old size() - 1
  * @post n and last_node(node with index size() - 1) is invalidated
  * @post node iterator of n and last_node is invalidated.
  * @post All outstanding incident iterators containing  last node
  *  and n are invalidated.
  * @post All outstanding edge iterators associasted with n
  * and last_node are invalidated.
  * Complexity(degree^2)
  */
  node_iterator remove_node(node_iterator n_it){
    return NodeIterator(this, remove_node(*n_it));
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Return Euclidean distance between two nodes.
    */
    double length() const{
      return norm(node1().position() - node2().position());
    } 

    /** Return value of the edge.
    */

    edge_value_type& value(){
      size_type n1 = node1_;
      size_type n2 = node2_;  
      if(node1_ > node2_){
        std::swap(n1, n2);
      }
      for(size_type i = 0; i < graph_->adj_nodes[n1].size(); i++){
        if(graph_->adj_nodes[n1][i].tail_node == n2){

          return graph_->adj_nodes[n1][i].val;
        }
      }
      assert(false);
    }
    /** Return value of the edge*/
    const edge_value_type& value() const{
      return value();
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_ == e.graph_ && node1_ == e.node1_ && node2_ == e.node2_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // (void) e;           // Quiet compiler warning
      if(graph_ == e.graph_){
        if(node1_ == e.node1_){
          return node2_ < e.node2_;
        } else {
          return node1_ < e.node1_;
        }
      }
      return graph_ < e.graph_;

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_;
    size_type node2_;

    // constructor using node indices 
    Edge(const Graph* graph, size_type n1, size_type n2):
      graph_(const_cast<Graph*>(graph)),node1_(n1), node2_(n2){      
    }  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(edge_begin(), edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    auto it = edge_begin();
    while(i != 0){
      ++it;
      i--;
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
    for(size_type i = 0; i < adj_nodes[a.index()].size(); i++){
      if(b.index() == adj_nodes[a.index()][i].tail_node){
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
    // check if the graph has the edge
    assert(has_node(a) && has_node(b) && !(a == b));
    if(!has_edge(a, b)){
      auto e_ab = incident_edge(b.index(), edge_value_type());
      auto e_ba = incident_edge(a.index(), edge_value_type());
      adj_nodes[a.index()].push_back(e_ab);
      adj_nodes[b.index()].push_back(e_ba);
    } 
    return Edge(this, a.index(), b.index());
  }

  /** Remove the edge between a and b.
  * return 1 if remove edge successfully, 0 otherwise.
  * @param[in] a is a node in the graph.
  * @param[in] b is a node in the graph.
  * @post if has_edge(a, b) is true, remove 
  * edge connecting @a a and @a b. 
  * new num_edges() = old num_edges() - 1;
  * Outstanding edge iterators after add_edge(a, b)
  * are invalidated.
  * Outstanding incident iterators related to 
  * a and b are invalidted.
  * Complexity : O(max(a.degree(), b.degree()).
  */

  size_type remove_edge(const Node& a, const Node& b){
    if(has_edge(a, b)){
      for(size_type i = 0; i < adj_nodes[a.index()].size(); i++){
        if(adj_nodes[a.index()][i].tail_node == b.index()){
          erase(adj_nodes[a.index()], i);
        }
      }
      for(size_type i = 0; i < adj_nodes[b.index()].size(); i++){
        if(adj_nodes[b.index()][i].tail_node == a.index()){
          erase(adj_nodes[b.index()], i);
        }
      }
      return 1;
    } else{
      return 0;
    }
  }


  /** Remove the edge e
  * return 1 if remove edge successfully, 0 otherwise.
  * @param[in] e is an edge to be removed.
  * @post Denote two end nodes as @a a and @a b,
  * if has_edge(a, b) is true, remove 
  * edge e that connects @a a and @a b. 
  * new num_edges() = old num_edges() - 1;
  * Outstanding edge iterators after add_edge(a, b)
  * are invalidated.
  * Outstanding incident iterators related to 
  * a and b are invalidted.
  * Complexity : O(max(a.degree(), b.degree())).
  */
  size_type remove_edge(const Edge& e){
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1, n2);
  }


  /** Remove the edge e pointed by edge iterator e_it
  * return 1 if remove edge successfully, 0 otherwise.
  * @param[in] e_it is an edge iterator, which points to 
  * and edge to be removed.
  * @post Denote two end nodes as @a a and @a b,
  * if has_edge(a, b) is true, remove 
  * edge e that connects @a a and @a b. 
  * new num_edges() = old num_edges() - 1;
  * Outstanding edge iterators after add_edge(a, b)
  * are invalidated.
  * Outstanding incident iterators related to 
  * a and b are invalidted.
  * Complexity : O(max(a.degree(), b.degree())).
  */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return edge_iterator(this, 0, 0);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_nodes.clear();
    adj_nodes.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    // Supply definitions AND SPECIFICATIONS for:
    /*
    * Dereference the itereator, get the node the iterator 
    * points to.
    * @pre 0 <= @a node_ind < size of graph
    */

    Node operator*() const{
      return Node(this->graph_, node_ind);
    }

    /*
    * Increment the node iterator, return 
    * reference to the incremented node itereator.
    * @pre 0 <= @a node_ind < size of graph
    */
    NodeIterator& operator++(){
      node_ind++;
      return *this;
    }
    /*
    * Return boolean of whether two node itereators are equal
    */
    bool operator==(const NodeIterator& node_iter) const{
      return (this->node_ind == node_iter.node_ind) && (this->graph_ && node_iter.graph_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_ind;
    // Node iterator constructor.
    NodeIterator(const Graph* graph, size_type index):
    graph_(const_cast<Graph*>(graph)), node_ind(index){
    }

  };

  // Supply definitions AND SPECIFICATIONS for:
  /*
  * Return a node iterator pointing to the first node of the graph.
  * @pre Graph is non-empty
  */
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }
  /*
  * Iterator for one node pass the last node of the graph.
  */
  node_iterator node_end() const{
    return node_iterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    /*
    * dereference the incident iterator, return an
    * edge such that node1 is the current node
    * and node2 is the outgoing node.
    * @pre 0 <= @a node_ind_ < graph size
    * @post returned edge, denoted as e, satisfies that:
    * e.node1() == Node(this, node_ind_), e.node2() is the node
    * connects with current node.
    */
    Edge operator*() const{
      assert(node_ind_ < this->graph_->size());
      size_type tail_index = graph_->adj_nodes[node_ind_][iter_ind_].tail_node;
      return Edge(this->graph_, node_ind_, tail_index);



    }
    /*
    * Increment the iterator, return the next incident iterator.
    */
    IncidentIterator& operator++(){
      iter_ind_++;
      return *this;
    }
    /*
    * Return true if the @a graph_, @a node ind_, and @a iter_ind_ 
    * of the two iterators are equal.
    * 
    * Return false otherwise.
    */
    bool operator==(const IncidentIterator& ii) const{
      return graph_ == ii.graph_ && node_ind_ == ii.node_ind_\
      && iter_ind_ == ii.iter_ind_;

    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_ind_;
    size_type iter_ind_;
    // Incident iterator constructor.
    IncidentIterator(const Graph* graph, size_type node_index, size_type iter_index):
    graph_(const_cast<Graph*>(graph)), node_ind_(node_index), iter_ind_(iter_index) {
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

    // Supply definitions AND SPECIFICATIONS for:
    /*
    * Dereference the itereator, get the edge the iterator 
    * points to.
    * @pre 0 <= @a edge_ind_ < number of edges in the graph
    */
    Edge operator*() const{
      assert(node_ind_ < this->graph_->size());
      size_type tail_index = graph_->adj_nodes[node_ind_][iter_ind_].tail_node;
      return Edge(this->graph_, node_ind_, tail_index);
    }

    /*
    * Increment the edge iterator, return 
    * reference to the incremented edge itereator.
    * @pre 0 <= @a edge < number of edges in the graph
    */
    EdgeIterator& operator++(){
      iter_ind_++;
      validate();
      return *this;
    }

    /*
    * Return true if the two iterators are equal.
    * 
    * Return false otherwise.
    */    
    bool operator==(const EdgeIterator& ei) const{
      return graph_ == ei.graph_ && node_ind_ == ei.node_ind_\
      && iter_ind_ == ei.iter_ind_;

    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_ind_;
    size_type iter_ind_;

    /** 
    * Validate the edge iterator
    */
    void validate(){
      while(node_ind_ < graph_->size()){
        while(iter_ind_ < graph_->adj_nodes[node_ind_].size()){
          
          if(node_ind_ < graph_->adj_nodes[node_ind_][iter_ind_].tail_node){
            // it is valid
            return;
          }
          iter_ind_++;
        }
        node_ind_++;
        iter_ind_ = 0;
      }
      return;
    }

    // edge_iterator constructor
    EdgeIterator(const Graph* graph, size_type node_index, size_type iter_index):\
    graph_(const_cast<Graph*>(graph)), node_ind_(node_index), iter_ind_(iter_index){
      validate();
    }

  };

  // Supply definitions AND SPECIFICATIONS for:
  /*
  * Return a edge_iterator pointing to the first edge of the graph.
  * @pre @a num_edges() > 0,  if the number of edges is 0, the returned edge iterator is invalid.
  */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0, 0);
  }
  /*
  * Return edge iterator whose index is equal to @a num_edges.
  */
  edge_iterator edge_end() const{
    return edge_iterator(this, num_nodes(), 0);
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  template<typename T>
  void erase(T& U, size_type i){
    U[i] = U.back();
    U.pop_back();
  }

  struct Internal_Node{
    Point p;
    node_value_type val;
    Internal_Node(Point pos, node_value_type value): p(pos), val(value){}
  };


  struct incident_edge{
    size_type tail_node;
    edge_value_type val;
    incident_edge(size_type n, edge_value_type v):tail_node(n), val(v){
    }
  };

  std::vector<Internal_Node> internal_nodes;
  std::vector<std::vector<incident_edge>> adj_nodes;

};

#endif // CME212_GRAPH_HPP
