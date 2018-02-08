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
  using node_value_type = V;
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
    : nodes_(), adj_(), connect()
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
  class Node {
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
      return graph_->nodes_[node_uid_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return a reference this node's value, which is modifiable. */
    node_value_type& value(){
      return graph_->nodes_[node_uid_].v_;
    }
    /** Return this node's value. cannot modify it.*/
    const node_value_type& value() const{
      return graph_->nodes_[node_uid_].v_;
    }
    /** Return this node's degree. */
    size_type degree() const{
      return graph_->adj_[index()].size();
    }
    /** Return a beginning incident iterator for the node. */
    IncidentIterator edge_begin() const{
      return Graph<V>::IncidentIterator(graph_, node_uid_, 0);
    }
    /** Return a ending incident iterator for the node. */
    IncidentIterator edge_end() const{
      return Graph<V>::IncidentIterator(graph_, node_uid_, degree());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ( graph_ == n.graph_ && n.position() == position() &&
            n.index() == node_uid_ && value() == n.value()){
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
      if (node_uid_ < n.index()){
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
    size_type node_uid_;

    Node (const Graph* graph, size_type node_uid)
        : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node (const Point & position,
            const node_value_type & v = node_value_type ()){
    // HW0: YOUR CODE HERE
    nodes_.push_back(internal_node{position, v});
    std::vector<size_type> empty;
    adj_.push_back(empty);

    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.graph_ == this && n.position() == node(n.node_uid_).position() &&
      n.value() == node(n.node_uid_).value()){
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
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
  class Edge {
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
    return connect.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, connect.at(i)[0], connect.at(i)[1]);
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
      if(adj_[a.index()][i] == b.index()){
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
    if (has_edge(a, b)){
      return Edge(this, a.index(), b.index());
    }
    if (has_edge(b, a)){
      return Edge(this, b.index(), a.index());
    }
    std::vector<size_type> add {a.index(), b.index()};
    connect.insert(std::pair<size_type,std::vector<size_type>>(connect.size(),
              add));
    adj_[a.index()].push_back(b.index());
    adj_[b.index()].push_back(a.index());
    return Edge(this, a.index(), b.index());
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
    connect.clear();
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
      return  Graph<V>::Edge(graph_, spawn_index,
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
      return graph_->edge(index);
    }

    /**incrementing the index for the iterator
    * @pre index < num_edges()
    */
    EdgeIterator& operator++(){
      index++;
      return *this;
    }

    /**overloading the == operator for comparison*/
    bool operator==(const EdgeIterator& ei) const{
      if(graph_ == ei.graph_ && index == ei.index){
        return true;
      }
      return false;
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type index;
    EdgeIterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index(index){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**create a beginning edge iterator for the graph*/
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**create a ending iterator for the graph*/
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges());
  }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
    Point point_;
    node_value_type v_;
  };

  std::vector<internal_node> nodes_;

  std::vector<std::vector<size_type>> adj_;
  std::map<size_type,std::vector<size_type>> connect;
};

#endif // CME212_GRAPH_HPP
