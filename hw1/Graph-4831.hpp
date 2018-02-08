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
template <typename V>
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
      graph_ = nullptr;

    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      return this->graph_->internal_nodes[this->ind_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->ind_;
    }

    // HW1: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
      // (void) n;          // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
      // (void) n;           // Quiet compiler warning
      if(this->graph_ == this->graph_){
        return (this->index() < n.index());
      } else {
        return this->graph_ < this->graph_;
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
    // HW0: YOUR CODE HERE
    // (void) n;            // Quiet compiler warning
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // (void) i;             // Quiet compiler warning

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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // (void) e;           // Quiet compiler warning

      return ((node1() == e.node1()) && (node2() == e.node2())) \
      || ((node1() == e.node2()) && (node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // (void) e;           // Quiet compiler warning
      if(e.graph_ == this->graph_){
        return this->edge_ind_ < e.edge_ind_;
      } else {
        return this->graph_ < this->graph_;
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
    size_type node1_;
    size_type node2_;
    size_type edge_ind_;

    // constructor using node indices and edge_index
    Edge(const Graph* graph, size_type n1, size_type n2, size_type ind):
      graph_(const_cast<Graph*>(graph)),node1_(n1), node2_(n2), edge_ind_(ind) {      
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return internal_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // (void) i;             // Quiet compiler warning
    return Edge(this, internal_edges[i].node1, internal_edges[i].node2, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // (void) a; (void) b;   // Quiet compiler warning
    // make sure the two nodes are in the graph
    if(this->has_node(a) && this->has_node(b)) {
      for(auto node_index: adj_nodes[a.index()]){
        if(b.index() == node_index){
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    // (void) a, (void) b;   // Quiet compiler warning

    // check if the graph has the edge

    assert(has_node(a) && has_node(b) && !(a == b));

    for(auto node_index: adj_nodes[a.index()]){
      if(b.index() == node_index){
        return Edge(this, a.index(), b.index(), -1);
      }
    }

    Internal_Edge new_edge = {a.index(), b.index()};
    adj_nodes[a.index()].push_back(b.index());
    adj_nodes[b.index()].push_back(a.index());
    internal_edges.push_back(new_edge);
    return Edge(this, a.index(), b.index(), num_edges() - 1);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    internal_nodes.clear();
    internal_edges.clear();
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

    // HW1 #2: YOUR CODE HERE
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
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_ind;
    // Node iterator constructor.
    NodeIterator(const Graph* graph, size_type index):
    graph_(const_cast<Graph*>(graph)), node_ind(index){
    }

  };

  // HW1 #2: YOUR CODE HERE
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
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
      return Edge(this->graph_, node_ind_, this->graph_->adj_nodes[node_ind_][iter_ind_], -1);

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
    // HW1 #3: YOUR CODE HERE
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /*
    * Dereference the itereator, get the edge the iterator 
    * points to.
    * @pre 0 <= @a edge_ind_ < number of edges in the graph
    */
    Edge operator*() const{
      return graph_->edge(edge_ind_);
    }

    /*
    * Increment the edge iterator, return 
    * reference to the incremented edge itereator.
    * @pre 0 <= @a edge < number of edges in the graph
    */
    EdgeIterator& operator++(){
      edge_ind_++;
      return *this;
    }

    /*
    * Return true if the @a graph_, @a edge_ind_,  
    * of the two iterators are equal.
    * 
    * Return false otherwise.
    */    bool operator==(const EdgeIterator& ei) const{
      return graph_ == ei.graph_ && edge_ind_ == ei.edge_ind_;

    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_ind_;

    // edge_iterator constructor
    EdgeIterator(const Graph* graph, size_type index):
    graph_(const_cast<Graph*>(graph)), edge_ind_(index){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /*
  * Return a edge_iterator pointing to the first edge of the graph.
  * @pre @a num_edges() > 0,  if the number of edges is 0, the returned edge iterator is invalid.
  */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0);
  }
  /*
  * Return edge iterator whose index is equal to @a num_edges.
  */
  edge_iterator edge_end() const{
    return edge_iterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct Internal_Node{
    Point p;
    node_value_type val;
    Internal_Node(Point pos, node_value_type value): p(pos), val(value){}
  };

  struct Internal_Edge{
    size_type node1;
    size_type node2;
    Internal_Edge(size_type n1, size_type n2): node1(n1), node2(n2){}
  };



  std::vector<Internal_Node> internal_nodes;
  std::vector<Internal_Edge> internal_edges;
  std::vector<std::vector<size_type>> adj_nodes;


};

#endif // CME212_GRAPH_HPP
