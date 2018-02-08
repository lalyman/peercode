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
  Graph() : Nodes_() , EdgesPosition_(), neighbor_() {
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

    /** Return this node's position. */
    const Point& position() const {
        return graph_->Nodes_[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        return uid_;
    }

    // HW1: PROBLEM 1

    /** Return a reference to the node's value
    * @post result == reference to the node's value
    */
    node_value_type& value(){
      return graph_->Nodes_[uid_].second;
    }

    /** Return the node's value.
    * @post result == node's value.
    */
    const node_value_type& value() const{
      return graph_->Nodes_[uid_].second;
    }

    // HW1: PROBLEM 3

    /** Return the degree of a node.
    *   Degree of a node is the number of edges incident to it.
    * @post result == @a n_degree (Where n_degree is number of edges
    *                              incident to the current node)
             @a n_degree < @a N (total number of nodes in the graph)
    */
    size_type degree() const{
      return graph_->neighbor_[uid_].size();
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
      return Nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {

    // at position @a i, store position and value for this node @a i
    Nodes_.push_back(std::make_pair(position,value));
    // at position @a i initialize  a vector of pairs for all the nodes and edges
    // neighbor to this node @a i
    neighbor_.push_back(std::vector<std::pair<size_type, size_type>> ());

    return Node(this, num_nodes()-1);
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
      // assert(i < num_nodes()-1);
      return Node(this,i);
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
      return graph_->node(idN1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(idN2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.uid_ == uid_) && (e.graph_ == graph_)){
        if (    (idN1_ == e.idN1_ && idN2_ == e.idN2_)  ||
                (idN1_ == e.idN2_ && idN2_ == e.idN2_))  { return true;}}
      else return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if ( graph_ == e.graph_) {
            //if edges belongs to the same graph, check the id of the edges
            return (uid_ < e.uid_);}
        else {
            //check in global sense which edge is greater without any disambiguity
            //using the pointer to the graph
        return (graph_ < e.graph_);}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Pointer to the graph
    Graph* graph_;
    //Id of the edge
    size_type uid_;
    //Id of node1
    size_type idN1_;
    //Id of node2
    size_type idN2_;
    //Additional constructor
    Edge(const Graph* graph_pointer, size_type uid, size_type idN1, size_type idN2)
        : graph_(const_cast<Graph*>(graph_pointer)), uid_(uid), idN1_(idN1), idN2_(idN2){
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
      // return n_edges ;
      return EdgesPosition_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, EdgesPosition_[i].first, EdgesPosition_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i=0; i<neighbor_[a.uid_].size(); i++) {
       if (b.uid_ == neighbor_[a.uid_][i].first)  return true;}
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

    //Search if we have already inserted this edge and if so, do not add it,
    //just return the edge already inserted
    for (unsigned int i=0; i<neighbor_[a.uid_].size(); i++) {
       if (b.uid_ == neighbor_[a.uid_][i].first) {
          return Edge(this, neighbor_[a.uid_][i].second, a.uid_, b.uid_);}}

    //add edges with node a < node b
    if (a<b) {EdgesPosition_.push_back(std::make_pair(a.uid_,b.uid_));}
    else {EdgesPosition_.push_back(std::make_pair(b.uid_,a.uid_));}

    //populate the vector neighbor_ with node and its neighbor
    neighbor_[a.uid_].push_back(std::make_pair(b.uid_,EdgesPosition_.size()-1));
    neighbor_[b.uid_].push_back(std::make_pair(a.uid_,EdgesPosition_.size()-1));
    //return the edge just created
    return Edge(this,EdgesPosition_.size()-1,a.uid_,b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      Nodes_.clear();
      EdgesPosition_.clear();
      neighbor_.clear();
  }

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
    *             edge.uid_   == @a id of the edge connecting node @a id
    *             edge.node1 = @a id of the node (node @a n)
    *             edge.node2 = node @a ni connecting to node @a id
    */
    Edge operator*() const{
      return Edge(graph_,                                   //pointer to the graph
                  graph_->neighbor_[node_id_][id_].second,  //edge connecting node n -- ni
                  node_id_,                                 //node n
                  graph_->neighbor_[node_id_][id_].first);  //node ni (connected to node n)
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
      return graph_->edge(id_);
    }

    /** Return the next id (iterator)
    * @post result is a EdgeIterator @a id_ incremented of 1
    */
    EdgeIterator& operator++(){
      ++id_;
      return *this;
    }

    /** Test if this EdgeIterator and @a e are equal
    * If they are the same,
    * these have same graph and @a id_.
    * @post result is true or false (equality or not)
    */
    bool operator==(const EdgeIterator& e) const{
      return ( (graph_==e.graph_) && (id_==e.id_));
    }

   private:
    friend class Graph;
      //Pointer to the graph
      const graph_type* graph_;
      //Id of the iterator
      size_type id_;
      //Additional constructor
      EdgeIterator(const Graph* graph, size_type id)
        : graph_(graph), id_(id){
      }
  };

  // ==========================================================================================================
  // ==========================================================================================================

  // PROBLEM 5

  /** This function returns the first edge iterator
  * @post result == edge_iterator
  *        *iter.graph_ == this.
  *        iter.id()    == 0
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }

  /** This function returns the last edge iterator
  * @post result == edge_iterator
  *        *iter.graph_ == this
  *        iter.id()    == n (total number of edges)
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->num_edges());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** @Nodes_ :
  *  vector with position and node value type
  *  at index i we have the position of node i and its value
  */
  std::vector<std::pair<Point, node_value_type>> Nodes_;

  /** @EdgesPosition_ :
  *  Vector of pairs of type size_type.
  *  At position i, it has a pairs of nodes for each edge i.
  */
  std::vector< pair<size_type, size_type> > EdgesPosition_;

  /** @neighbor_ :
  *  Vector of vectors which contain pairs of type size_type.
  *  At position i, it has a vector of pairs {node,edge}
  *  with nodes neighbor to node i and
  *  with edges neighbor to node i
  */
  std::vector< std::vector< std::pair< size_type, size_type >>> neighbor_;

};

#endif // CME212_GRAPH_HPP