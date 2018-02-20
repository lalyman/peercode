#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <map>
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
  struct internal_node;
  struct internal_edge;
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  //using node_value_type = V;

  typedef V node_value_type; 
  typedef E edge_value_type;
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // check if this position is in node list
      if (uid_<graph_->nodes_.size()){
        return graph_->nodes_[uid_].point;
      }
      assert(false);
    }

    Point & position () {
      if (uid_<graph_->nodes_.size()){
        return graph_->nodes_[uid_].point;
      }
      assert(false);        
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // check if the number is in the range
      if (uid_<graph_->nodes_.size()){
        return uid_;
      }
      assert(false);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the value of node. */    
    node_value_type& value(){
      if (uid_<graph_->size()){
        return graph_->nodes_[uid_].value;
      }
      assert(false);
    }

    /** Return the value of node. */ 
    const node_value_type& value() const{
      if (uid_<graph_->size()){
        return graph_->nodes_[uid_].value;
      }
      assert(false);      
    }

    /** Return the number of connected edges of the node. */     
    size_type degree() const{
      return graph_->edge_map_[uid_].size();
    }

    /** Return the first incident iterator,
     *  which is the first edge that connected to the given node 
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,uid_,size_type(0));
    }
    
    /** Return the last incident iterator,
     *  which is the last edge that connected to the given node 
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,uid_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((n.graph_==graph_) && (n.index()==uid_));
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
      return ((graph_<n.graph_)||((n.graph_==graph_) && (n.index()>uid_)));
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
      : graph_(const_cast<Graph*>(graph)),uid_(uid){
    }    
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
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {
    // HW0: YOUR CODE HERE
    // add the new node
    size_type n = size();
    internal_node new_node {position,value};
    nodes_.push_back(new_node);
    // add node to edge map
    edge_map_[nodes_.size()-1];
    return Node(this,n);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((n.graph_ == this)&&(n.uid_<size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i<num_nodes()){
        return Node(this,i);
    }
    assert(false);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (uid_<graph_->edges_.size()){
        return graph_->node(uid_n1_);
      }
      assert(false);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (uid_<graph_->edges_.size()){
        return graph_->node(uid_n2_);
      }
      assert(false);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((e.node1()==node1())&&(e.node2()==node2()))||
        ((e.node2()==node1())&&(e.node1()==node2())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((graph_<e.graph_)||((graph_==e.graph_)&&(this->uid_<e.uid_)));
    }

    /** Return the length of the edge. */
    double length() const {
      return norm(node1().position()-node2().position());
    }

    /** Return the value of the edge. */
    edge_value_type& value() {
      if (uid_<graph_->num_edges()){
        return graph_->edges_[uid_].value;
      }
      assert(false);            
    }  

    /** Return the value of the edge, but constant. */
    const edge_value_type& value() const {
      if (uid_<graph_->num_edges()){
        return graph_->edge_[uid_].value;
      }
      assert(false); 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    size_type uid_n1_;
    size_type uid_n2_;
    Edge(const Graph* graph, size_type uid, size_type uid_n1, size_type uid_n2)
      : graph_(const_cast<Graph*>(graph)),uid_(uid),uid_n1_(uid_n1),uid_n2_(uid_n2){
    }   
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i<num_edges()){
      return Edge(this,i,edges_[i].nid1,edges_[i].nid2);
    }
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
    if (has_node(a) && has_node(b)){
      for (auto i =a.edge_begin();i!=a.edge_end();++i){
        if ((*i).node2() == b) {
          return true;
        }
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value= edge_value_type()) {
    // HW0: YOUR CODE HERE
    
    // check if a and b are valid
    if (!(a==b) && has_node(a) && has_node(b)){ 
      // check if edge ab exists    
      // read the index of edge with node a
      for (auto i =a.edge_begin();i!=a.edge_end();++i){
        if ((*i).node2() == b) {
          return (*i);
        }
      }
      
      // add new edge
      internal_edge new_edge {a.index(),b.index(),value};
      edges_.push_back(new_edge);
      edge_map_[a.index()].push_back(num_edges()-1);
      edge_map_[b.index()].push_back(num_edges()-1);
      return edge(num_edges()-1);

    }
    assert(false);
  }


  /** Removes the edge that is connected to the two given node.  
   * 
   * @pre @a a and @b b are two given nodes.
   *
   * @param[in] a,b Nodes that are the endpts of edge to be removed.
   * 
   * @return the last index of edge if edge @a a and @a b does not exist,
   * @       else, return the next edge index
   *
   * @post If the edge connected to @a a and @a b is in the graph, 
   *       new num_edges() == old num_edges() - 1.
   *       The original edge is invalidated from both edge and edge map.
   *       Else, new num_edges() == old num_edges().
   *
   * Complexity: at most O(a.degree() + b.degree() + u1.degree() + u2.degree()). 
   */  
  size_type remove_edge(const Node& a, const Node& b) {

    // if the edge doesn't exist, return number of edges
    if (!has_edge(a,b)) {
      return num_edges();
    }

    // else, delete the adjacent edges from the edge map.
    else {
      std::vector<size_type> v = edge_map_.at(a.index());
      
      Edge e;
      size_type idx;
      // remove edges from node a map list
      for (size_type i=0;i<a.degree();++i){
        e = edge(v[i]);

        // find the edge connected to a and b
        if (((a==e.node1()) && (b==e.node2())) || ((a==e.node2()) && (b==e.node1()))){
          idx = v[i];
          edge_map_[a.index()][i] = edge_map_[a.index()].back();
          edge_map_[a.index()].pop_back(); 
          break; 
        }
      }


      // remove edges from node b map list
      std::vector<size_type> k = edge_map_.at(b.index());
      for (size_type i=0;i<b.degree();++i){
        e = edge(k[i]);  

        // find the edge connected to a and b   
        if (((a==e.node1()) && (b==e.node2())) || ((a==e.node2()) && (b==e.node1()))){
          edge_map_[b.index()][i] = edge_map_[b.index()].back();
          edge_map_[b.index()].pop_back(); 
          break; 
        }
      }

      // remove edges from edge list
      edges_[idx] = edges_.back();
      edges_.pop_back();

      // update index
      size_type u1 = edges_[idx].nid1;
      size_type u2 = edges_[idx].nid2;


      // update node a index
      for (auto i=edge_map_[u1].begin();i!=edge_map_[u1].end();++i){
        if ((*i)==num_edges()) {
          (*i) = idx;
        }
      }

      // update node b index
      for (auto i=edge_map_[u2].begin();i!=edge_map_[u2].end();++i){
        if ((*i)==num_edges()) {
          (*i) = idx;
        }
      }      
      return idx;

    }
  }

  /** Removes a given edge.  
   * 
   * @pre @a e is the given edge
   *
   * @param[in] e The edge to be removed.
   * 
   * @return the last index of edge if edge @a e does not exist,
   * @       else, return the next edge index
   *
  * @post If @a e is in the graph, new num_edges() == old num_edges() - 1.
   *       The original edge is invalidated from both edge and edge map.
   *       Else, new num_edges() == old num_edges().
   *
   * Complexity: at most O(a.degree() + b.degree() + u1.degree() + u2.degree()). 
   */ 
  size_type remove_edge(const Edge& e) {
    // call remove_edge and take in two nodes of e
    return remove_edge(e.node1(),e.node2());
  }


  /** Removes an edge with the given edge iterator.  
   * 
   * @pre @a e_it is the given edge iterator.
   *
   * @param[in] e_it The edge iterator where points to the edge that to be removed.
   * 
   * @return the end of edge iterator if the edge index is not smaller than number of edges,
   * @       else, return the next edge iterator.
   *
   * @post If @a e_it is in the graph, new num_edges() == old num_edges() - 1.
   *       The original edge is invalidated from both edge and edge map.
   *       Else, new num_edges() == old num_edges().
   *
   * Complexity: at most O(a.degree() + b.degree() + u1.degree() + u2.degree()). 
   */   
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type i = remove_edge(*e_it);

    // check  if the edge is at the boundary
    if (i<num_edges()){
      return EdgeIterator(this,i);  
    }
    else {
      return edge_end();
    }
    
  }


  /** Removes a given node and its connected edges.  
   * 
   * @pre @a n is the given node in graph.
   *
   * @param[in] n The node to be removed.
   * 
   * @return return the index of next node after @a n in the old node list.
   *
   * @post If @a n is in the graph, new num_nodes() == old num_nodes() - 1.
   *       new num_edges() == old num_edges() - n.degree()
   *       g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   *       g.node(n.index()) == n.
   *       The original node is invalidated.
   *       The edges that connected to the node are invalidated.
   *       Else, new num_nodes() == old num_nodes(), new num_edges() == old num_edges().   
   *
   * Complexity: at most O(num nodes()). 
   */   
  size_type remove_node(const Node& n){
    assert(has_node(n));
    // remove all adjacent edges
    size_type idx = n.index();
    while (edge_map_[idx].size()>0){
      remove_edge(edge(edge_map_[idx].back()));
    }

    // rempve the node and node value
    nodes_[idx] = nodes_.back();
    nodes_.pop_back();


    // update node index in the map
    size_type num = num_nodes();
    edge_map_[idx] = edge_map_[num];
    edge_map_.erase(num);
    if (idx < num){
        // find the swapped index and swap
        for (auto i=edge_map_[idx].begin();i!=edge_map_[idx].end();++i){
          if (edges_[*i].nid1==num){
            edges_[*i].nid1 = idx;
          }
          else{
            edges_[*i].nid2 =idx;
          }
        }
    }

    return idx;
  } 


  /** Removes a given node and its connected edges given it node iterator.  
   * 
   * @pre @a n_it is the given node iterator.
   *
   * @param[in] n_it The node iterator of the node to be removed.
   * 
   * @return the next node iterator if it is not in the node iterator list,
   *         else, return the end iterator of the node.  
   *
   * @post If @a n_it is in the graph, new num_nodes() == old num_nodes() - 1.
   *       new num_edges() == old num_edges() - n.degree().
   *       g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   *       g.node(n).index()) == n.
   *       The original node is invalidated.
   *       The edges that connected to the node are invalidated.
   *       Else, new num_nodes() == old num_nodes(), new num_edges() == old num_edges().
   *
   * Complexity: at most O(num nodes()). 
   */ 
  node_iterator remove_node(node_iterator n_it){
    size_type i = remove_node(*n_it);

    // if it is not the boundary, return this node index
    if (i<num_nodes()){
      return NodeIterator(this,i);
    }


    // otherwise, return the last node index
    else{
      return node_end();
    }

  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_nodes() == 0.
   *       g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   *       g.node(n).index()) == n.   
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edges_.clear();
    nodes_.clear();
    edge_map_.clear();
  }


  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator>{
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

    /** Return the Node with index uid_. */
    Node operator*() const{
      return graph_->node(uid_);
    }
    
    /** Return the next node iterator. */
    NodeIterator& operator++(){
      uid_+=1;
      return *this;
    }

    /** Test whether this node iterator and @a a are equal.
     *
     * Equal node iterator have the same graph and the same index.
     */
    bool operator==(const NodeIterator& a) const{
      return ((a.graph_==graph_)&&(a.uid_==uid_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;
    NodeIterator(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)),uid_(uid){
    }     
  };

  // HW1 #2: YOUR CODE HERE
    /** Return the first node iterator. */
    node_iterator node_begin() const{
      return (NodeIterator(this,size_type(0)));

    }
    /** Return the last node iterator. */
    node_iterator node_end() const{
      return (NodeIterator(this,size()));
    }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable<IncidentIterator>{
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

    /** Return the specific edge 
     *  @post return the edge which node 1 is this node
     */
    Edge operator*() const{
      size_type edge_idx = graph_->edge_map_[n_uid_][adj_edge_uid_];
      if (graph_->edges_[edge_idx].nid1==n_uid_){
        return Edge(graph_,edge_idx,n_uid_,graph_->edges_[edge_idx].nid2);
      }
      else{
        return Edge(graph_,edge_idx,n_uid_,graph_->edges_[edge_idx].nid1);
      }
    }
    
    /** Return the next incident iterator which is connected to this node. */    
    IncidentIterator& operator++(){
      adj_edge_uid_+=1;
      return *this;
    }

    /** Check if these two are equal
        These two incident iterators are equal if the have same graph, node, and connected edge.
     */
    bool operator==(const IncidentIterator& a) const{
      return ((graph_==a.graph_) && (a.adj_edge_uid_ == adj_edge_uid_) && (a.n_uid_ == n_uid_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type n_uid_;
    size_type adj_edge_uid_;

    IncidentIterator(const Graph* graph, size_type n_uid, size_type adj_edge_uid)
      : graph_(const_cast<Graph*>(graph)),n_uid_(n_uid),adj_edge_uid_(adj_edge_uid){
    } 
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator>{
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
    /** Return the specific edge. */
    Edge operator*() const{
      return graph_->edge(uid_);
    }

    /** Return the next edge iterator. */  
    EdgeIterator& operator++(){
      uid_+=1;
      return *this;
    }

    /** Check if these two are equal
        Two edge iterators are equal if the have same graph and same edge id.
     */
    bool operator==(const EdgeIterator& a) const{
      return ((a.graph_==graph_)&&(a.uid_==uid_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;
    EdgeIterator(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)),uid_(uid){
    }    

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the first edge iterator in the list. */ 
  edge_iterator edge_begin() const{
    return (EdgeIterator(this,size_type(0)));
  }

  /** Return the last edge iterator in the list. */
  edge_iterator edge_end() const{
    return (EdgeIterator(this,num_edges()));
  }



 private:


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  // struct for node
  struct internal_node {
    Point point; // The node's position
    node_value_type value; // The node's value 
  };

  // struct for edge
  struct internal_edge {
    size_type nid1, nid2; // Indices of the edge's vertices
    edge_value_type value; // The edge's value
  };


  std::vector<internal_node> nodes_; // node vector
  std::vector<internal_edge> edges_; // edge vector
  std::map<size_type,std::vector<size_type>> edge_map_; // adjacent edge map


};

#endif // CME212_GRAPH_HPP


