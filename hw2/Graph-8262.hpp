#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <utility>

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

  // Declarations of important internal types
  // later in the Graph's definition.

  // Predeclare the internal struct
  struct internal_Ndata;
  struct internal_Edata;

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
  Graph() 
  : node_(), edge_(), i2u_(), next_Nindex_(0), next_Eindex_(0), num_nodes_(0), 
    num_edges_(0), num_base_edges_(0) { 
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
    Point& position(){
      return get_node().position;
    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return get_node().index;
    }
    size_type i2u_size() const {
      return i2u_.size();
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return get_node().val;
    }
    const node_value_type& value() const{
      return get_node().val;
    }

    size_type edge_by_id (size_type id) const{
      return get_node().i_edge[id];
    }

    size_type degree() const{
      return get_node().degree;
    }

    size_type adj_size() {
      return get_node().adj_node.size();
    }

    size_type adj_node(size_type id) const{
      return get_node().adj_node[id];
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
      return Gnodes_ == n.Gnodes_ and get_node().index == n.index();
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
      return Gnodes_->node_[nindex_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1)
   * base graph size
   */
  size_type size() const {
    return num_nodes_;
  }

  /** Synonym for size()
   * accounts for removed nodes */
  size_type num_nodes() const {
    return i2u_.size();

  }

  /* This function takes care of removing node @a a and its incident edges
   * @param[in] a node to be removed from the graph
   *
   * @pre the node is part of the graph and has a valid index in i2u_ array
   *
   * @post 1. node a removed from i2u_ array, such that its only 
   *          visible via the base graph
   *       2. nodes are not ordered in i2u_ array i2u_[i] != node(i)
   *       3. i2u_size()(new) = i2u_size()(old) - 1
   *       4. num_nodes() (new) = num_nodes() (old) - 1
   *       5. all node incidents are removed from adj_node array
   *       6. degree of node a is zero
   *       7. num_edges() (new) = num_edges() (old) + (number of edges removed)
   *       8. node_end from node_iterator is affected by updating num_nodes()
   *       9. edge_end from incident_iterator is affected by updated degree of
   *          incident nodes
   *       8. edge_end from indicent_iterator is not affected since it preserves
   *          the full num_base_edges()
   *
   *  Complexity of O(1) assuming degree of incidents is very small
   *                    compared to graph size
   */
  void remove_node (const Node& a){
     Graph* g_ = a.Gnodes_;

     if (not has_node(a)) return;
     auto rm_i = a.index();
     //assert(rm_i < g_->i2u_.size())
     auto base_i = g_->i2u_[rm_i]; // node index in base graph

     /* loop over incidents with node a */
     for (size_type ii = 0; ii < node(rm_i).degree(); ii++) {
       auto i_adj = node(a.index()).adj_node(ii);
       auto n_adj = node_base(i_adj);

       auto n_degree = n_adj.degree();
       /* find and remove incident edges */
       for (size_type in = 0; in < n_degree; in++){
         if (n_adj.adj_node(in) == base_i){
            auto i2u_adj = g_->i2u_[i_adj];
            g_->node_[i2u_adj].adj_node[in] = 
                   std::move(g_->node_[i2u_adj].adj_node.back());
            g_->node_[i2u_adj].i_edge[in] = 
                   std::move(g_->node_[i2u_adj].i_edge.back());
            g_->node_[i2u_adj].adj_node.pop_back();
            g_->node_[i2u_adj].i_edge.pop_back();
            g_->node_[i2u_adj].degree--;
            g_->num_edges_--;
         }
       }
     }
     g_->node_[base_i].degree = 0;
     g_->node_[base_i].adj_node.clear();
     g_->node_[base_i].i_edge.clear();
 
     if (i2u_.size() > 0){
        g_->i2u_[rm_i] = std::move( g_->i2u_.back() );
        g_->i2u_.pop_back();
     }
  }

  node_iterator remove_node(node_iterator it){
    assert(it != node_end());
    auto n = *it;
    remove_node(n);
    return it;
  }

  /* no specifications since it's not being used for remove sphere  */
  bool remove_edge(const Node& a, const Node& b){
    if ( not has_edge(a, b) ) return false;
    Graph* g_ = a.Gnodes_;
    auto a_ = i2u_[a.index()];
    auto b_ = i2u_[b.index()];
    auto a_deg = Node(this, a_).degree();
    auto b_deg = Node(this, b_).degree();

    for (size_type i = 0; i<a_deg; i++){
       if (b_ == Node(this, a_).adj_node(i)) {
          g_->node_[a_].adj_node[i] = 
                       std::move( g_->node_[a_].adj_node.back() );
          g_->node_[a_].adj_node.pop_back();
          g_->node_[a_].degree--;
          break;
       }
    }
    for (size_type i = 0; i<b_deg; i++){
       if (a_ == Node(this, b_).adj_node(i)) {
          g_->node_[b_].adj_node[i] = 
                       std::move( g_->node_[b_].adj_node.back() );
          g_->node_[b_].adj_node.pop_back();
          g_->node_[b_].degree--;
          break;
       }
    }
    num_edges_--;
    return not has_edge(a,b);
  }
  bool remove_edge(Edge& e){
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1, n2);
  }
  edge_iterator remove_edge(edge_iterator it){
    auto e = *it;
    remove_edge(e);
    return it;
  }
  void set_node_val( size_type i, node_value_type val){
    assert(i < num_nodes());
    node_base(i).set_val(val);
  }

  /** Return a proxy object for node @a i. */
 /** Return the node with index @a i.
  * @pre 0 <= @a i < num_nodes()
  * @post result_node.index() == i
  *
  * Complexity: O(1).
  */
  Node node(size_type i) const {
    if (i2u_.size() <= i){
        std::cout << "size i2u " << i2u_.size() << std::endl;
        std::cout << "node index: "<< i << std::endl;
        assert(false);
    }
    return Node(this, i2u_[i]);
  }

  Node node_base(size_type i) const {
    return Node(this, i);
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  //Node add_node(const Point& position, const node_value_type& = node_value_type()) {
  Node add_node(const Point& position, const node_value_type& value) {
    (void) value;
    internal_Ndata new_node;

    new_node.position = position;
    new_node.index = next_Nindex_;
    new_node.degree = 0;
    i2u_.push_back(num_nodes_);

    node_.push_back(new_node);
    ++num_nodes_;
    ++next_Nindex_;

    return Node(this, next_Nindex_-1);
  }

  Node add_node(const Point& position) {
    internal_Ndata new_node;

    new_node.position = position;
    new_node.index = next_Nindex_;
    new_node.degree = 0;
    i2u_.push_back(num_nodes_);
    
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
    return ((n.index() < num_nodes()) and (n.Gnodes_ == this));
  }

  void set_node_value (Node& n, int i){
    n.val = i;
  }

  size_type node_edge_id(Node& n, size_type i){
    auto nn = n.index();
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
      return Gedges_->node_base(get_edge().n1);
    }
    /** Return the other node of this Edge */
    Node node2() const {
      return Gedges_->node_base(get_edge().n2);
    }
    double length() const{
      return norm(node1().position() - node2().position());
    }
    size_type index() const {
      return get_edge().index;
    }
    void set_val(edge_value_type x){
      get_edge().val = x;
    }
    edge_value_type& value() {
      return get_edge().val;
    }
    const edge_value_type& value() const {
      return get_edge().val;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==( const Edge& a ) const {

      if (Gedges_==a.Gedges_) {
         if ((get_edge().n1 == a.node1().index() and 
              get_edge().n2 == a.node2().index()) 
             or 
             (get_edge().n1 == a.node2().index() and
              get_edge().n2 == a.node1().index()))
         return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& a) const {
      size_type e_n1_i, e_n2_i, e_max, e_min;
      size_type ee_n1_i, ee_n2_i, ee_max, ee_min;

      e_n1_i = get_edge().n1; e_n2_i = get_edge().n2;
      ee_n1_i =  a.node1().index(); ee_n2_i =  a.node2().index();

      if (e_n1_i<e_n2_i){ e_max = e_n2_i; e_min = e_n1_i;}
      else { e_max = e_n1_i; e_min = e_n2_i;}

      if (ee_n1_i<ee_n2_i){ ee_max = ee_n2_i; ee_min = ee_n1_i;}
      else { ee_max = ee_n1_i; ee_min = ee_n2_i;}

      if ( (Gedges_ < a.Gedges_) or ((Gedges_ == a.Gedges_) 
           and ((e_min < ee_min) or ((e_min == ee_min) and 
                                     (e_max < ee_max))))) {
        return true; 
      }
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

  size_type num_base_edges() const{
    return num_base_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_base_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {

    auto a_ = Node(this, i2u_[a.index()]);
    auto b_ = Node(this, i2u_[b.index()]);
    auto a_deg = a_.degree();
    auto b_deg = b_.degree();

    for (size_type i=0; i < a_deg ; i++) {
       if (b_.index() == Node(this, a_.index()).adj_node(i)) {
          for (size_type j=0; j < b_deg ; j++) {
             if (a_.index() == Node(this, b_.index()).adj_node(j)) {
                return true; }}}}
    return false;

  }

  /* complexity of O(num_edges) */
  size_type find_edge (Node& n, size_type adj_i) const{
     auto n1_id = n.index(); // base graph index
     auto n2_id = n.adj_node(adj_i); // base graph index

     for (size_type i = 0; i<num_base_edges(); i++){
        if ((n1_id == Edge(this, i).node1().index()  and 
             n2_id == Edge(this, i).node2().index())  or 
            (n1_id == Edge(this, i).node2().index()  and
             n2_id == Edge(this, i).node1().index())){
          return i;
        }
     }
     assert(false);
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

      auto i2u_a = a.index();
      auto i2u_b = b.index();

      assert (i2u_a < size());
      assert (i2u_b < size());

      // add node1 info
      node_[i2u_a].degree++;
      node_[i2u_a].i_edge.push_back(new_edge.index);
      node_[i2u_a].adj_node.push_back(b.index());
      // add node2 info
      node_[i2u_b].degree++;
      node_[i2u_b].i_edge.push_back(new_edge.index);
      node_[i2u_b].adj_node.push_back(a.index());

      edge_.push_back(new_edge);

      ++num_edges_;
      ++next_Eindex_;
      ++num_base_edges_;
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
    num_base_edges_ = 0;
    node_.clear();
    edge_.clear();
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
    return NodeIterator(this, num_nodes());
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
       temp_n = G_iter_->node_base(n_);
       ii = temp_n.node_edge_id(i_);
       /* using find_edge method with complexity O(num_edges) made 
          the simulation extremely slow. I know that STL containers
          should be used effeciently, but this compromised the quality
          that's why i have adj_node and i_edge containers */
       //return  G_iter_->edge(G_iter_->find_edge(temp_n,i_));

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
     return EdgeIterator(this, num_base_edges());
  }

 private:

   // Graph class's internals:
   // helper functions, data members, and so forth.

   // Internal data for nodes class
   struct internal_Ndata {
     Point position;
     node_value_type val {};
     size_type index;
     size_type degree;
     std::vector<size_type> i_edge;
     std::vector<size_type> adj_node;
   };

   // Internal data for edge class
   struct internal_Edata {
     size_type n1;
     size_type n2;
     size_type index;
     edge_value_type val;
   }; 

   std::vector<internal_Ndata> node_;
   std::vector<internal_Edata> edge_;
   std::vector<size_type> i2u_;

   size_type next_Nindex_;
   size_type next_Eindex_;
   size_type num_nodes_;
   size_type num_edges_;
   size_type num_base_edges_;
};

#endif // CME212_GRAPH_HPP
