#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>


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

  // Declarations of important internal types
  // later in the Graph's definition.

  // Predeclare the internal struct
  struct internal_Ndata;
  struct internal_Edata;

//  friend class NodeIterator; 

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

  using node_value_type = V;

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
  : node_(), edge_(), next_Nindex_(0), next_Eindex_(0), num_nodes_(0), 
    num_edges_(0) { 
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return get_node().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return get_node().index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return get_node().val;
    }
    const node_value_type& value() const{
      return get_node().val;
    }

    size_type edge_by_id ( size_type id) const{
      return get_node().i_edge[id];
    }

    size_type degree() const{
      return get_node().degree;
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(Gnodes_, nindex_, 0);
    }

    incident_iterator edge_end() const{
      return IncidentIterator(Gnodes_, nindex_, degree());
    }

    void set_val(node_value_type x){
      get_node().val = x;
    }

    size_type node_edge_id(size_type i){
      return get_node().i_edge[i];
    }

    internal_Edata& get_edge(size_type i) const {
      return Gnodes_->edge_[i];
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */

    bool operator==(const Node& n) const {
      return get_node().index == n.index();
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
      return get_node().index < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    Graph* Gnodes_;
    // Unique id for a node
    size_type nindex_;
    /** Private Constructor */
    Node(const Graph* Gnodes, size_type index)
    : Gnodes_(const_cast<Graph*>(Gnodes)), nindex_(index) {
    }
    /** Helper method to return the object
      * This loops over the nodes until it finds the node with the
      * correct index.
      */
    internal_Ndata& get_node() const {
      assert(Gnodes_->size() > nindex_);
      return Gnodes_->node_[nindex_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  void set_node_val( size_type i, node_value_type val){
    assert(i < num_nodes());
    node(i).set_val(val);
  }

  /** Return a proxy object for node @a i. */
 /** Return the node with index @a i.
  * @pre 0 <= @a i < num_nodes()
  * @post result_node.index() == i
  *
  * Complexity: O(1).
  */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, i);
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {

    internal_Ndata new_node;

    new_node.position = position;
    new_node.index = next_Nindex_;
    new_node.degree = 0;

    node_.push_back(new_node);
    ++num_nodes_;
    ++next_Nindex_;

    return Node(this, next_Nindex_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.index() < size();
  }

  void set_node_value (Node& n, int i){
    n.val = i;
  }

  size_type node_edge_id(Node& n, size_type i){
    return n.edge_by_id(i);
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
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Gedges_->node(get_edge().n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Gedges_->node(get_edge().n2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use for private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    Graph* Gedges_;
    // This edge's unique identification number
    size_type eindex_;
    // Private Constructor 
    Edge(const Graph* Gedges, size_type index)
    : Gedges_(const_cast<Graph*>(Gedges)), eindex_(index) {
    }

    /*Edge(const Graph* Gedges, size_type index)
    : Gedges_(const_cast<Graph*>(Gedges)), eindex_(index) {
    }*/

    internal_Edata& get_edge() const {
      return Gedges_->edge_[eindex_];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {

    if ( adj_list[a.index()][b.index()] or 
         adj_list[b.index()][a.index()] ) { return true;}
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {

    if (not has_edge(a,b)) {
      internal_Edata new_edge;
      new_edge.n1 = a.index();
      new_edge.n2 = b.index();
      new_edge.index = next_Eindex_;

      // update adjacency list
      adj_list[a.index()][b.index()] = true;
      // add node1 info
      node_[new_edge.n1].degree++;
      node_[new_edge.n1].i_edge.push_back(new_edge.index);
      // add node2 info
      node_[new_edge.n2].degree++;
      node_[new_edge.n2].i_edge.push_back(new_edge.index);

      edge_.push_back(new_edge);
      ++num_edges_;
      ++next_Eindex_;
      return Edge(this, next_Eindex_-1);
    }
    return Edge(this, next_Eindex_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    next_Nindex_ = 0;
    next_Eindex_ = 0;
    num_nodes_ = 0;
    num_edges_ = 0;
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
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Giter_->node(i_); 
    }
    
    NodeIterator& operator++() {
      i_++;
      return *this;
    }

    bool operator==(const NodeIterator& a) const {
      return i_ == a.i_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // Pointer back to the Graph container
    Graph* Giter_; 
    size_type i_;

    /** Private Constructor */
    NodeIterator(const Graph* Giter, size_type index)
    : Giter_(const_cast<Graph*>(Giter)), i_(index) {
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const{
    return NodeIterator(this, size());
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
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
       size_type ii;
       Node temp_n;
       temp_n = G_iter_->node(n_);
       ii = temp_n.node_edge_id(i_);
       return G_iter_->edge(ii);
    }
    IncidentIterator& operator++(){
       i_++;
       return *this;
    }
    bool operator==(const IncidentIterator& a) const{
        return G_iter_==a.G_iter_ and i_ == a.i_;
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE

    Graph* G_iter_;
    size_type n_;
    size_type i_;

    /** Private Constructor */
    IncidentIterator(const Graph* G_iter, size_type n, size_type index)
    : G_iter_(const_cast<Graph*>(G_iter)), n_(n), i_(index) {
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
    Edge operator*() const{
       return G_iter_->edge(i_); 
    }
    EdgeIterator& operator++(){
       ++i_;
       return *this; 
    }
    bool operator==(const EdgeIterator& a) const{
       return G_iter_==a.G_iter_ and i_ == a.i_; 
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    // Pointer back to the Graph container
    Graph* G_iter_;
    size_type i_;

    /** Private Constructor */
    EdgeIterator(const Graph* G_iter, size_type index)
    : G_iter_(const_cast<Graph*>(G_iter)), i_(index) {
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
      return EdgeIterator(this, num_edges());
  }

 private:

   // Graph class's internals:
   // helper functions, data members, and so forth.

   // Internal data for nodes class
   struct internal_Ndata {
     Point position;
     node_value_type val;
     size_type index;
     size_type degree;
     std::vector<size_type> i_edge;
   };

   // Internal data for edge class
   struct internal_Edata {
     size_type n1;
     size_type n2;
     size_type index;
   }; 

   std::vector<internal_Ndata> node_;
   std::vector<internal_Edata> edge_;

   /* adjacency list used to speedup the lookup of an existing edge
    * adj[node1.index][node2.index] = true for existing edge, otherwise 
    * it's false
    */
   std::unordered_map<size_type,std::unordered_map<size_type, bool>> adj_list;

   size_type next_Nindex_;
   size_type next_Eindex_;
   size_type num_nodes_;
   size_type num_edges_;
};

#endif // CME212_GRAPH_HPP
