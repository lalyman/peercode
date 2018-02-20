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

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V,typename E>
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

  /** User-specified value, of type node value type */
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** User-specified value, of type edge value type */
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

  /** struct the defines all the information
  *  for a node (including id, position, value, neighbor nodes)
  */
  struct Data_node{
    //id of the node
    size_type id_;
    //position of the node
    Point position_;
    //value of the node
    node_value_type value_;
    //vector with neighbor nodes and edges's value
    std::vector< std::pair< size_type , edge_value_type >> neighbor_;

    Data_node( size_type id,
               Point position,
               node_value_type value,
               std::vector<std::pair<size_type, edge_value_type>> neighbor):
               id_(id), position_(position), value_(value),neighbor_(neighbor ) {}
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : Node_info_(),index_nodei_(){
  }

  /** Default destructor */
  ~Graph() = default;

  // =========================================================================
  // NODES
  // =========================================================================

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
     */

    Node() {
    }

    /** Return this node's position. */
    Point& position() const {
      return graph_->Node_info_[uid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->Node_info_[uid_].id_;
    }

    /** Return a reference to the node's value
    * @post result == reference to the node's value
    */
    node_value_type& value(){
      return graph_->Node_info_[uid_].value_;
    }

    /** Return the node's value.
    * @post result == node's value.
    */
    const node_value_type& value() const{
      return graph_->Node_info_[uid_].value_;
    }

    /** Return the degree of a node.
    *   Degree of a node is the number of edges incident to it.
    * @post result == @a n_degree (Where n_degree is number of edges
    *                              incident to the current node)
             @a n_degree < @a N (total number of nodes in the graph)
    */
    size_type degree() const{
    return graph_->Node_info_[uid_].neighbor_.size();
    }

    /** Return the first incident iterator of the current node
    * @post result = first incident iterator
    *                result <= @a e (@a e == any other valid incident
    *                                iterator for the current node)
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_ ,uid_ ,0);
    }

    /** Return the last incident iterator of the current node
    * @post result = last incident iterator
    *                result >= @a e (@a e == any other valid incident
    *                                iterator for the current node)
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,uid_,this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        return ((n.uid_ == uid_) && (n.graph_ == graph_));
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

        if ( graph_ == n.graph_) {
            //if nodes belongs to the same graph, check the id of the nodes
            return (uid_ < n.uid_);}
        else {
            //check in global sense which node is greater without any disambiguity
            //using the pointer to the graph
        return (graph_ < n.graph_);}
    }

   private:
   // Allow Graph to access Node's private member data and functions.
   friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    //Pointer to the graph
    Graph* graph_;
    //Id of the node
    size_type uid_;
    //Additional constructor
    Node(const Graph* graph, size_type id)
         : graph_(const_cast<Graph*>(graph)), uid_(id) {
    }

  };

  // ==========================================================================================================
  // ==========================================================================================================

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index_nodei_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position : The new node's position
   * @param[in] value    : The new node's value of type node_value_type
   * @post new @a num_nodes() == @a old num_nodes() + 1
   * @post result_node.index() == @a old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {

    // at position @a i initialize  a vector of pairs for all the nodes and edges' value
    // neighbor to this node @a i
    std::vector<std::pair<size_type, edge_value_type>> neighbor_vector;

    // at position @a i, store position, value , neighbor informationfor for this node @a i
    Node_info_.push_back(Data_node( index_nodei_.size(), position,value,neighbor_vector ));
    // store the index of node @a i (unique identifier)
    index_nodei_.push_back(Node_info_.size()-1);

    return Node(this,num_nodes()-1);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

      //check also if they both belogs to the same graph
      //and if the index is within range
      return ( (n.graph_ == this) && (size()-1 >= n.index()) );
  }

  /** Return the node with index @a i.
   * @pre @a 0 <= @a i < @a num_nodes()
   * @post result_node.index() == @a i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, index_nodei_[i]);
  }

  // =========================================================================
  // EDGES
  // =========================================================================

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
      return Node(graph_,idN1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,idN2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.graph_ == graph_){
        if (    (idN1_ == e.idN1_ && idN2_ == e.idN2_)  ||
                (idN1_ == e.idN2_ && idN2_ == e.idN2_))  { return true;}}
      else return false;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if ( graph_ == e.graph_) {
            //if edges belongs to the same graph, check the id of the nodes
            return ( (idN1_ < e.idN1_ ) && (idN2_ < e.idN2_ )  );}
        else {
            //check in global sense which edge is greater without any disambiguity
            //using the pointer to the graph
        return (graph_ < e.graph_);}
    }

    /** Return the euclidean distance between two nodes connecting
     *  the current edge
     */
    double length() const{
      Node node_1 = node1();
      Node node_2 = node2();
      return norm( node_1.position() - node_2.position() );
    }

    /** Return a reference to the edge's value
    * @post result == reference to the edge's value
    */
    edge_value_type& value(){
      unsigned id_node2 = 0;
      //loop thorugh all the neighbor nodes of @a node 1
      for (unsigned it=0; it < graph_->Node_info_[idN1_].neighbor_.size(); ++it){
        //search in this vector the node == @a node 2
        if (graph_->Node_info_[idN1_].neighbor_[it].first == idN2_){
          id_node2 = it;
          break;}}
      return graph_->Node_info_[idN1_].neighbor_[id_node2].second;
    }

    /** Return the edge's value.
    * @post result == edge's value.
    */
    const edge_value_type& value() const{
      unsigned id_node2 = 0;
      //loop thorugh all the neighbor nodes of @a node 1
      for (unsigned it=0; it < graph_->Node_info_[idN1_].neighbor_.size(); ++it){
        //search in this vector the node == @a node 2
        if (graph_->Node_info_[idN1_].neighbor_[it].first == idN2_){
          id_node2 = it;
          break;}}
      return graph_->Node_info_[idN1_].neighbor_[id_node2].second;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Pointer to the graph
    Graph* graph_;
    //Id of node1
    size_type idN1_;
    //Id of node2
    size_type idN2_;
    //Additional constructor
    Edge(const Graph* graph_pointer, size_type idN1, size_type idN2)
        : graph_(const_cast<Graph*>(graph_pointer)), idN1_(idN1), idN2_(idN2){
    }

   friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  // ==========================================================================================================
  // ==========================================================================================================

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(this->edge_begin(),this->edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type count = 0;
    for (auto it = this->edge_begin(); it != this->edge_end(); ++it){
      if (count==i) return *it;
      ++count;}

    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int it=0; it<Node_info_[a.uid_].neighbor_.size(); it++) {
       if (b.uid_ == Node_info_[a.uid_].neighbor_[it].first) return true;}
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

    //Search if we have already inserted this edge and if so, do not add it,
    if (!has_edge(a,b)){
      //populate the vector neighbor_ with neighbor node index and neighbor edge value
      Node_info_[a.uid_].neighbor_.push_back(std::make_pair(b.uid_,value));
      Node_info_[b.uid_].neighbor_.push_back(std::make_pair(a.uid_,value));
    }

    return Edge(this,a.uid_,b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Node_info_.clear();
    index_nodei_.clear();
  }

  // =========================================================================
  // REMOVE
  // =========================================================================

  /** It removes a node and the neighbor edges to this node from the graph
   * @param[in] @a node1  node to erase
   * @return the index deleted associated to node1
   *
   * @pre @a index of node1 and @a node2 are in the graph
   * @pre 0 <= @a index node1 < @a num_nodes() in the graph
   * @post new num_nodes() == old num_nodes()-1
   * @post new num_edges() < old num_edges()
   *       (if node1 has some neighbor)
   *
   * All the data connected to node1 is erased
   * (included all the edges attached to node1)
   * Also node1 is invalidated.
   * Performs at most O (num_nodes() + num_edges() ) operation
   */
  size_type remove_node(const Node& node1){
    // First delete edges in the adjacency matrix

    //delete all the edges that connect node1 to other nodes
    auto it = node1.edge_begin();
    while (it!=node1.edge_end()){
      Edge edge_i = *it;
      Node node2 = edge_i.node2();
      remove_edge(node1,node2);}

    //delete the index of node1 in the vector @a index_nodei_
    unsigned id_remove = Node_info_[node1.uid_].id_;
    index_nodei_.erase(index_nodei_.begin()+ id_remove);
    //invalidate the node1 giving an index of -1
    Node_info_[node1.uid_].id_ = -1;

    //iterate through all nodes of the graph
    for (auto it = this->node_begin(); it!= this->node_end(); ++it){
      Node node_i = *it;
      // node_i with index higher than node1 is decreased by 1
      if (node_i.uid_ > node1.uid_) {Node_info_[node_i.uid_].id_ -= 1;}}

    return id_remove;
  }

  /** It removes a node and the neighbor edges to this node from the graph
   * @param[in] @a n_it  iterator pointing to the node to be erase
   * @return n_it
   *
   * @pre 0 <= @a index n_it < @a num_nodes() in the graph
   * @post new num_nodes() == old num_nodes()-1
   * @post new num_edges() < old num_edges()
   *       (if n_it is pointing to a node with neighbor)
   * @post n_it points to the next node
   *       (in case there is no next node, it points to edge_end() )
   *
   * All the data connected to node1 is erased
   * (included all the edges attached to node1)
   *  Iterators [n_it, last) are invalidated
   * Performs at most O (num_nodes() + num_edges() ) operation
   */
  node_iterator remove_node(node_iterator n_it){
    Node node1 = *n_it;
    remove_node(node1);
    return n_it;
  }

   /** It removes an edge from the graph
   * @param[in] @a node1  node 1 to erase
   * @param[in] @a node2  node 2 to erase
   * @return 0
   *
   * @pre @a node1 and @a node2 are in the graph
   * @post new num_edges() == old num_edges()-1.
   *
   * All the data connected to the edge connecting
   * @a node1 and @a node2 is erased
   * Performs at most O(num_edges()) operation
   */
  size_type remove_edge(const Node& node1, const Node& node2){

    //loop through all the neighbor nodes of node1
    for (unsigned it=0; it < Node_info_[node1.uid_].neighbor_.size(); ++it){
      //check if id of neighbor node is equal to id of node2
      if (Node_info_[node1.uid_].neighbor_[it].first == node2.uid_){
        //remove the pair element at position begin()+it
        Node_info_[node1.uid_].neighbor_.erase( Node_info_[node1.uid_].neighbor_.begin()+it );
        break;}}

    //loop through all the neighbor nodes of node2
    for (unsigned it=0; it < this->Node_info_[node2.uid_].neighbor_.size(); ++it){
      //check if id of neighbor node is equal to id of node1
      if (Node_info_[node2.uid_].neighbor_[it].first == node1.uid_){
        //remove the element at position begin()+it
        Node_info_[node2.uid_].neighbor_.erase(Node_info_[node2.uid_].neighbor_.begin()+it);
        break;}}

    return 0;
  }


  /** It removes an edge from the graph
   * @param[in] @a edge1  edge to erase
   * @return 0
   *
   * @pre @a node1 and @a node2 are connected
   *      by the edge @a edge1
   * @post new num_edges() == old num_edges()-1
   *
   * All the data connected to the edge @a edge1
   * is erased
   * Performs at most O(num_edges()) operation
   */
  size_type remove_edge(const Edge& edge1){

    //edge1 connects node1 and node2
    Node node1 = edge1.node1();
    Node node2 = edge1.node2();
    //use the remove_edge method (which takes two nodes as input)
    return remove_edge(node1,node2);
  }



  /** It removes an edge from the graph
   * @param[in] @a e_it  edge iterator that points to the edge to delete
   * @return @a e_it edge iterator that points to the next edge
   *                 (if edge exists, otherwise points to end())
   *
   * @pre @a e_it point to @a edge1 and
   *      @a node1 and @a node2 are connected by the edge @a edge1
   * @post new num_edges() == old num_edges()-1
   *
   * All the data connected to the edge @a edge1
   * is erased. The iterator points to the next edge if
   * it is valid, otherwise search for the next iterator's vald edge
   * Performs at most O(num_edges()) operation
   */
  edge_iterator remove_edge(edge_iterator e_it){

    //remove the edge that the iterator is pointing to
    Edge edge_i = *e_it;
    remove_edge(edge_i);

    //take the id of the node that the iterator is pointing to
    size_type iter_node_id = e_it.graph_->index_nodei_[e_it.n_idx_];
    if ( e_it.id_ == e_it.graph_->Node_info_[iter_node_id].neighbor_.size() ){
      e_it.id_=0;
      ++e_it.n_idx;
      //loop till the node index is less than total number of nodes
      while( e_it.n_idx_ < e_it.graph_->num_nodes() ){
        size_type node_id = e_it.graph_->index_nodei_[e_it.n_idx_];
        //loop till iterator id is less than the number of neighbor edges to node @a i
        while (e_it.id_ < e_it.graph_->Node_info_[node_id].neighbor_.size()){
          //if we find an id greater than node @a i, return the next edge in the neighbor_
          if (e_it.graph_->Node_info_[node_id].neighbor_[e_it.id_].first > node_id){
             return e_it;}
          ++e_it.id_;}

        e_it.id_= 0;
        ++e_it.n_idx_;}
    }
    //if edges invalid, return end()
    return e_it;
  }

  // ==========================================================================================================
  // ==========================================================================================================

  // =========================================================================
  // NODE ITERATOR
  // =========================================================================

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

    // PROBLEM 2

    /** Return the iterator's value
    * @post result is a node with the same iterator's id
    *       call the method node to return that node
    *       (from this graph and with the same @a id)
    */
    Node operator*() const{
      return graph_->node(id_);
    }

    /** Return the next id (iterator)
    * @post result is a NodeIterator @a id_ incremented of 1
    */
    NodeIterator& operator++(){
      ++id_;
      return *this;
    }

    /** Test if this NodeIterator and @a n are equal
    * If they are the same,
    * these have same graph and @a id_.
    * @post result is true or false (equality or not)
    */
    bool operator==(const NodeIterator& n) const{
      return ((n.graph_ == graph_) && (n.id_ == id_));
    }

   private:
    friend class Graph;
      //Pointer to the graph
      const graph_type* graph_;
      //Id of the iterator
      size_type id_;
      //Additional constructor
      NodeIterator(const Graph* graph, size_type id)
       : graph_(graph), id_(id) {
      }
  };

  // ==========================================================================================================
  // ==========================================================================================================

  /** This function returns the first node iterator
  * @post result == node_iterator
  *        *iter.graph_ == this.
  *        iter.id()    == 0
  */
  node_iterator node_begin() const{
     return NodeIterator(this,0);
  }

  /** This function returns the last node iterator
  * @post result == node_iterator  such that
  *        *iter.graph_ == this
  *        iter.id()    == n (total number of nodes)
  */
  node_iterator node_end() const{
     return NodeIterator(this, this->num_nodes());
  }

  // =========================================================================
  // INCIDENT ITERATOR
  // =========================================================================

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>   {
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

    // PROBLEM 3

    /** Return the value of the iterator
    * @post result is a edge
    *       where edge.graph_ == pointer to the graph
    *             edge.node1 = @a id of the node (node @a n)
    *             edge.node2 = node @a ni connecting to node @a id
    */
    Edge operator*() const{
      return Edge(graph_,                                             //pointer to the graph
                  node_id_,                                           //node n
                  graph_->Node_info_[node_id_].neighbor_[id_].first); //node ni (connected to node n)
    }

    /** Return the next id (iterator)
    * @post result is a IncidentIterator @a id_ incremented of 1
    */
    IncidentIterator& operator++(){
      ++id_;
      return *this;
    }

    /** Test if this IncidentIterator and @a e are equal
    * If they are the same,
    * these have same graph, @a node_id_ and @a id_.
    * @post result is true or false (equality or not)
    */
    bool operator==(const IncidentIterator& e) const{
      return (  (graph_   == e.graph_)     &&
                (node_id_ == e.node_id_)   &&
                (id_   == e.id_) );
    }

   private:
    friend class Graph;
      //Pointer to the graph
      const graph_type* graph_;
      //Id of the node
      size_type node_id_;
      //Id of the iterator
      size_type id_;
      //Additional constructor
      IncidentIterator(const Graph* graph_pointer,size_type node_id, size_type id)
        : graph_(graph_pointer), node_id_(node_id), id_(id){
      }
  };

  // ==========================================================================================================
  // ==========================================================================================================

  // =========================================================================
  // EDGE ITERATOR
  // =========================================================================

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

    // Problem 5

    /** Return the iterator's value
    * @post result is an edge with the same iterator's id
    *       call the method edge to return that edge
            (from this graph and with the same @a id)
    */
    Edge operator*() const{
      return Edge(graph_,
                  graph_->index_nodei_[n_idx_],             //id of node1
                  graph_->Node_info_[graph_->index_nodei_[n_idx_]].neighbor_[id_].first); //id of node2 (neighbor of node1)
    }


   /** Return the next id (iterator)
    * This method return an id of the iterator that point to the edge(i+1)
    * ( the previous id iterator was pointing to edge(i) )
    * @post result is a EdgeIterator @a id_
    */
    EdgeIterator& operator++(){
      ++id_;

      //loop till node index is within range [ 0, index_nodei_.size() )
      while (n_idx_ < graph_->num_nodes()){
        //take the id of the node @a i (unique identifier)
        size_type node_id_i = graph_->index_nodei_[n_idx_];
        //loop till iterator is less than the number of neighbor edges to node @a i
        while (id_ < graph_->Node_info_[node_id_i].neighbor_.size()){
          //if we find an id greater than node @a i, return
          if ( node_id_i < graph_->Node_info_[node_id_i].neighbor_[id_].first  ){
             return *this;}
          ++id_;}
        //start again the first while loop till the condition above is satisifed
        id_= 0;
        ++n_idx_;}

      return *this;
    }

    /** Test if this EdgeIterator and @a e are equal
    * If they are the same,
    * these have same graph and indeces.
    * @post result is true or false (equality or not)
    */
    bool operator==(const EdgeIterator& e) const{
      return ( (id_==e.id_) &&
               (graph_==e.graph_) &&
               (n_idx_==e.n_idx_) );
    }

   private:
    friend class Graph;
      //Pointer to the graph
      const graph_type* graph_;
      //Id of the iterator
      size_type id_;
      //Id of the node
      size_type n_idx_;
      //Additional constructor
      EdgeIterator(const Graph* graph_pointer, size_type id, size_type n_idx)
       : graph_(graph_pointer), id_(id), n_idx_(n_idx){
     }
  };

  // ==========================================================================================================
  // ==========================================================================================================


 /** This function returns the first edge iterator
  * @post result == edge_iterator
  *        *iter.graph_ == this.
  *        iter.id() and iter.n_idx are such that
  *        the iterator is pointing to the first edge in neighbor_
  */
  edge_iterator edge_begin() const{
    auto ei = EdgeIterator(this,0,0);

    //loop till id is within range
    while (ei.n_idx_< ei.graph_->num_nodes()){
      // (unique identifier)
      size_type node_id = ei.graph_->index_nodei_[ei.n_idx_];

      //consider the node that has node_id and loop through the neighbor edges
      while (ei.id_ < ei.graph_->Node_info_[node_id].neighbor_.size()){
        //if id of the node is greater than the node_id, return the first edge in the neighbor_
        if (ei.graph_->Node_info_[node_id].neighbor_[ei.id_].first > node_id){
           return ei;}
        ++ei.id_;}

      ei.id_= 0;
      ++ei.n_idx_;}
    //retrun last() if there is no element in neighbor
    return ei;
  }


 /** This function returns the last edge iterator
  * @post result == edge_iterator
  *        *iter.graph_ == this
  *        iter.id()    == n (total number of nodes)
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this,0,this->num_nodes());
  }


 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** @Node_info_ :
  *  Vector of objects Data_node
  *  At position i, it has all the information relative to node i
  */
  std::vector<Data_node> Node_info_;

 /** @index_nodei_ :
  *  it maps indeces and unique identifiers;
  */
  std::vector<size_type> index_nodei_;

};

#endif // CME212_GRAPH_HPP
