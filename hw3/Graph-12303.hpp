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
template<typename V, typename E>
class Graph {
 private:
  struct nodeinfo{
    Point p_;
    V v_;
    unsigned idx_;
  };
  struct edgeinfo{
    std::vector<unsigned> edge_vec_;    //store node uid
    E e_;
    unsigned idx_;
  };
  std::vector<nodeinfo> Nodes;
  std::vector<unsigned> i2u_;
  std::vector<edgeinfo> Edges;
  std::vector<unsigned> i2u_e;
  std::vector<std::vector<unsigned>> Adj;  
  unsigned size0;
  unsigned num_edges0;

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
  using graph_type = Graph<V, E>;

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
  typedef V node_value_type;
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    size0 = 0;
    num_edges0 = 0;
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
      uid_n = 0;
      //graph_n = nullptr;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      const Point& result_point = graph_n->Nodes[uid_n].p_;
      return result_point;
      //return Point();
    }

    /** Modifiable node position. */
    Point& position() {
      assert(valid());
      Point& result_point = graph_n->Nodes[uid_n].p_;
      return result_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //assert(valid());
      unsigned result = graph_n->Nodes[uid_n].idx_;
      return result;
      //return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return reference of value of the node  
     */
    node_value_type& value(){
      assert(valid());
      node_value_type& result_v = graph_n->Nodes[uid_n].v_;
      return result_v;
    }

    /** Return a const reference of value of the node
     */ 
    const node_value_type& value() const{
      assert(valid());
      const  node_value_type& result_v = graph_n->Nodes[uid_n].v_;      
      return result_v;
}

    /** Return number of edges incident to the node */
    size_type degree() const{
      assert(valid());
      return graph_n->Adj[uid_n].size();
    }

    /** Return begin of incident iterator */
    incident_iterator edge_begin() const{
      assert(valid());
      return IncidentIterator(graph_n, index(), 0);
    }

    /** Return end of incident iterator of this node */
    incident_iterator edge_end() const{
      assert(valid());
      return IncidentIterator(graph_n, index(), degree());
    } 

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_n == n.graph_n) && (uid_n == n.uid_n));
      //(void) n;          // Quiet compiler warning
      //return false;
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
      //assert(graph_n == n.graph_n);
      return uid_n < n.uid_n;
      //(void) n;           // Quiet compiler warning
      //return false;
    }

    /** Check if this is a valid node */
    bool valid() const {
      return uid_n >= 0 && uid_n < graph_n->Nodes.size() && graph_n->Nodes[uid_n].idx_ < graph_n->i2u_.size() && graph_n->i2u_[graph_n->Nodes[uid_n].idx_] == uid_n;
    }
 
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* graph_n;
    size_type uid_n;

    Node(const graph_type* graph0, size_type uid0): graph_n(const_cast<graph_type*>(graph0)), uid_n(uid0){
    } 
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return number of edges incident to the node with unique id @a id0 */
  size_type node_degree(size_type uid0) const{
    assert(uid0 < Nodes.size());
    return Adj[uid0].size();
  }
 
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
   *
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    size0 += 1;
    Point result_point= const_cast<Point&>(position);
    Nodes.push_back(result_point);
    return Node(this, size0 - 1);
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
  }*/


  /** Add a node to the graph, returning the added node.
   *  @param[in] position: The new node's position
   *  @param[in] value0: The new node's value. Optional.
   *  @post new num_nodes() == old num_nodes()
   *  @post result_node.index() == old num_nodes()
   *  Complexity: O(1)
   */
  Node add_node(const Point& position, const node_value_type& value0 = node_value_type()){
    size0 += 1;
    nodeinfo result_nodeinfo;
    result_nodeinfo.p_ = const_cast<Point&>(position);
    result_nodeinfo.v_ = const_cast<node_value_type&>(value0);
    result_nodeinfo.idx_ = i2u_.size();
    Nodes.push_back(result_nodeinfo);
    unsigned temp_uid = Nodes.size() - 1;
    i2u_.push_back(temp_uid);
    Adj.push_back({});
    return Node(this, temp_uid);
}

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.valid();
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((0 <= i) && (i < size0));
    return Node(this, i2u_[i]);
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
      node1_uid = 0;
      node2_uid = 0;
      edge_uid = 0;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type idx = graph_e->Nodes[node1_uid].idx_;
      return graph_e->node(idx);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type idx = graph_e->Nodes[node2_uid].idx_;
      return graph_e->node(idx);      // Invalid Node
    }
   
    /** Swap node1 and node2 */
    void swap() { 
      int temp = node1_uid;
      node1_uid = node2_uid;
      node2_uid = temp;
    }
 
    /** Get Euclidean length of the edge */
    double length() const{
      assert(valid());
      return norm(node1().position() - node2().position());
    }  

    /** Return edge value */
    edge_value_type& value(){
      assert(valid());
      return graph_e->Edges[edge_uid].e_;
    }

    /** Return edge value, const version */
    const edge_value_type& value() const{
      assert(valid());
      const edge_value_type& result = graph_e->Edges[edge_uid].e_;
      return result;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_e != e.graph_e) return false;
      bool temp = ((node1_uid == e.node1_uid) && (node2_uid == e.node2_uid)) || ((node1_uid == e.node2_uid) && (node2_uid == e.node1_uid)) ;
      return temp && edge_uid == e.edge_uid;
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //return ((node1().index() + node2().index()) < (e.node1().index() + e.node2().index()));
      if(graph_e != e.graph_e) return graph_e < e.graph_e;
      return value() < e.value();
      //(void) e;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_e;
    size_type node1_uid;
    size_type node2_uid;
    size_type edge_uid;
    Edge(const graph_type* graph_e0, size_type node1_uid0, size_type node2_uid0, size_type edge_uid0): graph_e(const_cast<graph_type*>(graph_e0)), node1_uid(node1_uid0), node2_uid(node2_uid0), edge_uid(edge_uid0){
      assert(node1_uid0 != node2_uid0);
    }

    /** Check if this is a valid edge */
    bool valid() const {
      return node1().valid() && node2().valid() && edge_uid >= 0 && edge_uid < graph_e->Edges.size() && graph_e->Edges[edge_uid].idx_ < graph_e->i2u_e.size() && graph_e->i2u_e[graph_e->Edges[edge_uid].idx_] == edge_uid;
    } 
    
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return i2u_e.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    //if(i < num_edges0) return Edge();
    return Edge(this, Edges[i2u_e[i]].edge_vec_[0], Edges[i2u_e[i]].edge_vec_[1], i2u_e[i]);
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_n == b.graph_n);
    if(!(has_node(a) && has_node(b))) return false;
    size_type a0 = a.uid_n;
    size_type b0 = b.uid_n;
    bool test = false;/*
    size_type i = 0;
    while(i < i2u_e.size()){
      size_type i0 = i2u_e[i];
      test = (a0 == Edges[i0].edge_vec_[0] && b0 == Edges[i0].edge_vec_[1]) || (a0 == Edges[i0].edge_vec_[1] && b0 == Edges[i0].edge_vec_[0]);
      if(test){
        return true;
      }
      i++;
    }
    return false;*/
    for(auto i0 : Adj[a0]){
      test = (a0 == Edges[i0].edge_vec_[0] && b0 == Edges[i0].edge_vec_[1]) || (a0 == Edges[i0].edge_vec_[1] && b0 == Edges[i0].edge_vec_[0]);
      if(test) return true;
    }
    return false;
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;
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
  /**Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_n == b.graph_n);
    assert(has_node(a) && has_node(b));
    size_type a0 = a.index();
    size_type b0 = b.index();
    bool test = false;
    size_type i = 0;
    while(i < num_edges0){
      test = ((a0 == Edges[i][0] && b0 == Edges[i][1]) || (a0 == Edges[i][1] && b0 == Edges[i][0]));
      if(test){
        return Edge(this, a.index(), b.index(), i);
      }
      i++;
    }
    num_edges0 += 1;
    std::vector<size_type> temp = {a0, b0};
    Edges.push_back(temp);
    Adj[a0].push_back(num_edges0 - 1);
    Adj[b0].push_back(num_edges0 - 1);
    return Edge(this, a0, b0, num_edges0 - 1);
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }*/

  /** Add a edge 
   *  @param[in] a, b: two nodes of the edge
   *  @param[in] value0: edge value, optional */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value0 = edge_value_type()){
    assert(a.graph_n == b.graph_n);
    assert(has_node(a) && has_node(b));
    size_type a0 = a.uid_n;
    size_type b0 = b.uid_n;
    bool test = false;
    size_type i = 0;
    while(i < num_edges0){ 
      int i0 = i2u_e[i];
      test = ((a0 == Edges[i0].edge_vec_[0] && b0 == Edges[i0].edge_vec_[1]) || (a0 == Edges[i0].edge_vec_[1] && b0 == Edges[i0].edge_vec_[0]));
      if(test){
        return Edge(this, a0, b0, i0);
      }
      i++;
    }
    num_edges0 += 1;
    std::vector<size_type> temp = {a0, b0};
    edgeinfo result_edgeinfo;
    result_edgeinfo.edge_vec_ = temp;
    result_edgeinfo.e_ = const_cast<edge_value_type&>(value0);
    result_edgeinfo.idx_ = i2u_e.size();
    int temp_uid = Edges.size();
    Edges.push_back(result_edgeinfo);
    i2u_e.push_back(temp_uid);
    Adj[a0].push_back(temp_uid);
    Adj[b0].push_back(temp_uid);
    return Edge(this, a0, b0, temp_uid);
  }  

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Nodes.clear();
    i2u_.clear();
    Edges.clear();
    i2u_e.clear();
    Adj.clear();
    size0 = 0;
    num_edges0 = 0; 
    // HW0: YOUR CODE HERE
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>  {
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereference NodeIterator
     *  @post: return the Node pointed by NodeIterator*/
    Node operator*() const{
      //assert(itr_ni < graph_ni->size());
      return graph_ni->node(itr_ni);
    }

    /** Move the iterator forward by 1
     *  @post: return the node_iterator that points to the next node in graph*/ 
    node_iterator& operator++(){
      ++itr_ni;
      return *(this);
    }

    /** Test whether this NodeIterator equals to @a n_itr0.
     *  Equal NodeIterators point to the same Node. */
    bool operator==(const node_iterator& n_itr0) const{
      //return graph_ni->node(itr_ni) == *n_itr0;
      return  graph_ni == n_itr0.graph_ni && itr_ni == n_itr0.itr_ni; 
}



   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_ni;
    size_type itr_ni;    //node index
    /** Construct a NodeIterator. @a itr_ni is initialized as 0 by default. */
    NodeIterator(const graph_type* graph0): graph_ni(const_cast<graph_type*>(graph0)), itr_ni(0){
    }
    NodeIterator(const graph_type* graph0, size_type itr0): graph_ni(const_cast<graph_type*>(graph0)), itr_ni(itr0){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return begin of node_iterator in this graph */
  node_iterator node_begin() const{
    return NodeIterator(this);
  }

  /** Return end of node_iterator in this graph */
  node_iterator node_end() const{
    return NodeIterator(this, size0);
  }
  

  //
  // Incident Iterator
  //

  /** A public interface to use @a Adj 
   *  @a n_id0 is index of the first node
   *  @a i is ith neighbour */
  /** size_type Get_Adj(size_type n_id0, size_type i) const {
    assert(n_id0 < size0);
    size_type uid0 = i2u_[n_id0];    
    assert(i < Adj[uid0].size());
    return Adj[uid0][i];
  }*/

  /** Return edge using @a Adj
   *  @a n_id0 is index of the first node
   *  @a i is ith neighbour */
  Edge Get_Adj(size_type n_id0, size_type i) const {
    assert(n_id0 < size0);
    size_type uid0 = i2u_[n_id0];    
    assert(i < Adj[uid0].size());
    size_type uid1 = Adj[uid0][i];
    return Edge(this, Edges[uid1].edge_vec_[0], Edges[uid1].edge_vec_[1], uid0);
  }

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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    
    /** Dereference IncidentIterator */
    Edge operator*() const{
      assert(node_id < graph_ii->size());
      Node ntemp = graph_ii->node(node_id);
      assert(ntemp.valid());
      assert(itr_ii < ntemp.degree());
      Edge result = graph_ii->edge(graph_ii->Edges[graph_ii->Adj[graph_ii->i2u_[node_id]][itr_ii]].idx_);
      if(result.node1() != ntemp){
        result.swap();
      }
      return result;
    }


    /** Move Incident Iterator forward by one */
    IncidentIterator& operator++() {
      ++itr_ii;
      return *(this);
    }

    /** Return true if this incident iterator equals to @a itr_ii0 */
    bool operator==(const IncidentIterator& itr_ii0) const {
      return (graph_ii == itr_ii0.graph_ii) && (node_id == itr_ii0.node_id) && (itr_ii == itr_ii0.itr_ii);
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_ii;
    size_type node_id;    //node index
    size_type itr_ii;    //index in Adj
    IncidentIterator(const graph_type* graph_ii0, size_type node_id0, size_type itr_ii0): graph_ii(const_cast<graph_type*>(graph_ii0)), node_id(node_id0), itr_ii(itr_ii0){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>  {
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

    /** Dereference Edge Iterator. Return the edge it points to. */
    Edge operator*() const{
      return graph_ei->edge(itr_ei);    
    }

    /** Move the itrator forward by one */
    EdgeIterator& operator++(){
      ++itr_ei;     
      return *(this);
    }

    /** Return true if this == @a e_itr0 */
    bool operator==(const EdgeIterator& e_itr0) const{
      //return graph_ei->edge(itr_ei) == *e_itr0;
      return graph_ei == e_itr0.graph_ei && itr_ei == e_itr0.itr_ei;
    }  

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_ei;
    size_type itr_ei;
    EdgeIterator(const graph_type* graph_ei0, size_type itr_ei0): graph_ei(const_cast<graph_type*>(graph_ei0)), itr_ei(itr_ei0){
    }
    EdgeIterator(const graph_type* graph_ei0): graph_ei(const_cast<graph_type*>(graph_ei0)), itr_ei(0){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /** Return edge iterator that points to begin */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Return edge iterator that points to end */
  edge_iterator edge_end() const{
    return EdgeIterator(this, i2u_e.size());
  }

  /** HW2 remove method for node 
   *  @param[out] return if node is successfully removed
   *  @pre: @a r_node is a valid node already existing in the graph
   *  @post: r_node.valid() is false
   *  @post: new size() == old size() - 1 
   *  @post: for all node that have node.index() < @a r_node.index()
             old node.index() == new node.index()
             for all node that have node.index() > @a r_node.index()
             old node.index() == new node.index() + 1 
   *  @post: all edges incident to this node are also removed
             all elements in Adj that has unique id of @a r_node are erased 
   *  Complexity: O(degree()*log(num_nodes()))*/   
  bool remove_node(const Node& r_node){
    //assert(r_node.valid());
    if(!has_node(r_node)) return false;
    size_type r_idx = r_node.index();
    size_type r_uid = r_node.uid_n;
    //remove related edges
    for(size_type itr0 = 0; itr0 < Adj[r_uid].size(); itr0++){
      size_type e_uid = Adj[r_uid][itr0];
      //Modify Adj[]
      for(size_type id0 = 0; id0 < size0; id0++){
        size_type itr1 = i2u_[id0];
        if(itr1 != r_uid){
          auto found = std::find(Adj[itr1].begin(), Adj[itr1].end(), e_uid); 
          if(found != Adj[itr1].end())
          Adj[itr1].erase(found);
        }
      }
      //Modify i2u_e[]
      i2u_e.erase(std::find(i2u_e.begin(), i2u_e.end(), e_uid));
    }
    //Modify Edges[].idx_
    for(size_type itr_uid = 0; itr_uid < Edges.size(); ++itr_uid){
      auto found = std::find(i2u_e.begin(), i2u_e.end(), itr_uid);
      if(found != i2u_e.end()){
        size_type index0 = std::distance(i2u_e.begin(), found);
        Edges[itr_uid].idx_ = index0;
      }
    }
    
    //Modify num_edges
    num_edges0 = i2u_e.size();

    //Delete content of Adj[r_uid]
    Adj[r_uid] = {};

    //Delete the node
    for(size_type i = r_idx + 1; i < size0; i++){
      Nodes[i2u_[i]].idx_ -= 1; 
    }
    i2u_.erase(i2u_.begin() + r_idx);
    --size0;
    return true;
  }

  /** HW2 remove node iterator 
   *  @param[out]: return the input node iterator @a r_n_it
   *  @pre: @a r_n_it is not out of scope
   *  @post: same as remove_node() 
   *  @Complexity: same as remove_node() */ 
  node_iterator remove_node(node_iterator r_n_it){
    assert(r_n_it.itr_ni < i2u_.size());
    node_type temp = *r_n_it;
    remove_node(temp);
    return r_n_it;
  }

  /** HW2 remove edge by nodes
   *  @param[out]: return whether the removal is successful
   *  @pre: graph has nodes @a r_node1 and @a r_node2
   *  @pre: has_edge(r_node1, r_node2) is true
   *  @post: edge.valid() == false 
   *  @Complexity: O(num_edges()) */

  bool remove_edge(const Node& r_node1, const Node& r_node2){
    if(!(has_node(r_node1) && has_node(r_node2))) return false;
    if(!has_edge(r_node1, r_node2)) return false;
      size_type a0 = r_node1.uid_n;
      size_type b0 = r_node2.uid_n;
      bool test = false;
      size_type i = 0;
      unsigned i0;
      while((!test) && (i < i2u_e.size())){ 
        i0 = i2u_e[i];
        test = ((a0 == Edges[i0].edge_vec_[0] && b0 == Edges[i0].edge_vec_[1]) || (a0 == Edges[i0].edge_vec_[1] && b0 == Edges[i0].edge_vec_[0]));
        i++;
      }
      i--;    //i is index, i0 is uid of edge
      //Modify Adj[]
      Adj[a0].erase(std::find(Adj[a0].begin(), Adj[a0].end(), i0));
      Adj[b0].erase(std::find(Adj[b0].begin(), Adj[b0].end(), i0));
      //Modify i2u_e and Edges[]
      for(size_type j = i + 1; j < num_edges0; j++){
        Edges[i2u_e[j]].idx_ -= 1; 
      }
      i2u_e.erase(std::find(i2u_e.begin(), i2u_e.end(), i0));
      //Modify number of edges
      --num_edges0;
      return true;   
  }

  /** HW2 remove edge by edge 
   *  same as bool remove_edge(const Node& r_node1, const Node& r_node2) */
  bool remove_edge(const Edge& r_edge){
    if(!has_edge(r_edge)) return false;
    size_type uid = r_edge.edge_uid;
    size_type idx = Edges[uid].idx_;
    size_type a0 = r_edge.node1_uid;
    size_type b0 = r_edge.node2_uid;
    //Modify Adj[]
    Adj[a0].erase(std::find(Adj[a0].begin(), Adj[a0].end(), uid));
    Adj[b0].erase(std::find(Adj[b0].begin(), Adj[b0].end(), uid));
    //Modify i2u_e and Edges[]
    for(size_type j = idx + 1; j < num_edges0; j++){
      Edges[i2u_e[j]].idx_ -= 1; 
    }
    i2u_e.erase(i2u_e.begin() + idx);
    //Modify number of edges
    --num_edges0;
    return true;
  }

  /** HW2 remove edge iterator 
   *  @param[out]: return the input edge iterator @a r_e_it
   *  @pre: @a r_e_it is not out of scope
   *  @post: same as remove_edge() 
   *  @Complexity: same as remove_edge() */ 
  edge_iterator remove_edge(edge_iterator r_e_it){
    assert(r_e_it.itr_ei < i2u_e.size()); 
    edge_type temp = *r_e_it;
    remove_edge(temp);
    return r_e_it;
  }
};

#endif // CME212_GRAPH_HPP
