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
template <typename V = double, typename E = double>
class Graph {
 private:

  /*** Predeclare internal structs for node and edge ***/
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

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
  Graph() 
    : nodes_(), edges_(), num_nodes_(0), num_edges_(0), next_node_uid_(0), next_edge_uid_(0) {
   /*** Graph constructed with zero size ***/
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
  // using node_value_type = V;
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

    /* Creates an invalid Node object */
    Node() {
      /*** Invalid Node; must be created by calling Graph methods  ***/
    }

    /*** Return this node's position. ***/
    const Point& position() const {
      /*** access this node's position (via its uid) through graph_ ptr ***/
      return graph_->nodes_[node_uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      /*** access this node's index (uid) through graph_ ptr ***/
      return graph_->nodes_[node_uid_].node_uid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Add given value to this node
     *
     */
    node_value_type& value(node_value_type value)
    {
      graph_->nodes_[node_uid_].value = value;
      return graph_->nodes_[node_uid_].value;
    }

    /** Return the value of this node
     *
     */
    const node_value_type& value() const
    {
      return graph_->nodes_[node_uid_].value;
    }

    /** Return the position of this node
     * Access the interal node through the graph pointer
     * and the node_uid_ that are passed to instantiate an
     * object of the Node class. Return the position from
     * the struct
     */
    Point& position()
    {
        return graph_->nodes_[node_uid_].position;
    }

    /** Return the degree of this node
     *
     */
    size_type degree() const
    {
      return graph_->adj_list[node_uid_].size();
    }
    
    /** Return an IncidentIterator beginning at the first
     *  edge for this Node. The order depends on how edges
     *  were added to the graph.
     */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(graph_, node_uid_);
    }

    /** Return an IncidentIterator of the last edge for this
     *  Node. The order depends on how edges were added to the graph
     */
    incident_iterator edge_end() const
    {
      return IncidentIterator(graph_, node_uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /*** Test if this node's ID and graph ptr are equivalent to n's ***/
      return ((n.index() == this->node_uid_) && (n.graph_ == this->graph_));
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
      /*** Compare node uids for ordering purposes ***/
      return (n.index() > this->node_uid_);
    }

   private:
    /*** Private data members and methods for Node ***/
    //Pointer back to the Graph container
    Graph* graph_;
    //This node's unique identification number
    size_type node_uid_;
    // Value for node
    // node_value_type V;
    /* Private constructor for Graph::Node object */
    Node(const Graph* graph, size_type node_uid) 
      : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid) {
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    /*** size_ == num_nodes_ ***/
    return num_nodes_;
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

  /*** Add new internal_node to nodes_ vector of internal_nodes ***/
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type()) {
    //instantiate an internal node to add to nodes_ vector
    internal_node new_node;
    //assign data attributes to new_node
    new_node.position = position;
    new_node.node_uid = next_node_uid_;
    new_node.value = node_val;
    //append new_node to nodes_ vector
    nodes_.push_back(new_node);
    // push back empty pair into the adjacency list
    std::vector<std::pair<size_type, size_type>> empty_list;
    adj_list.push_back(empty_list);
    //update number of nodes and next node ID
    ++num_nodes_;
    ++next_node_uid_;
    return Node(this, next_node_uid_-1);   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /*** same graph if ptr to memory for node is the same as graph_ ptr ***/ 
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    /*** return a proxy object for node @ i ***/
    assert(i < num_nodes());
    return Node(this, i);
  }

  size_type remove_node(const Node&)
  {
    
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
      /*** Invalid Edge; must be created by calling Graph methods ***/
    }

    /** Return a node of this Edge */
    Node node1() const {
      /*** access this edge's node1 member through graph_ ptr ***/
      return graph_->edges_[edge_uid_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      /*** access this edge's node2 member through graph_ ptr ***/
      return graph_->edges_[edge_uid_].node2; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((this->graph_ == e->graph_) && (this->edge_uid_ == e.edge_uid_));
      // (void) e;           // Quiet compiler warning
      // return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->edge_uid_ < e.edge_uid_);
      // (void) e;           // Quiet compiler warning
      // return false;
    }

    /** Return the length of this edge. It is a constant set to 1.0
     *  for problem 1
     */
    double length() const
    {
      return 1.0;
    }

    // Set the value of the edge to the argument passed and return the value
    edge_value_type& value(edge_value_type value)
    {
      graph_->edges_[edge_uid_].value = value;
      return graph_->edges_[edge_uid_].value;
    }

    // Return the value of the edge
    const edge_value_type& value() const
    {
      return graph_->edges_[edge_uid_].value;
    }

   private:
    /*** Private data members and methods for Edge ***/
    // Pointer back to the Graph container
    Graph* graph_;
    //This Edge's unique ID number
    size_type edge_uid_;
    /* private constructor Graph::Edge */
    Edge(const Graph* graph, size_type edge_uid)
      : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    /*** number of edges present in the graph ***/
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /*** return a proxy object for edge @ i ***/
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    /*** Verify that a and b nodes are witin this graph_ ***/
    if (this == a.graph_ && this == b.graph_) {
      //iterate through edges_ vector to see if an edge contains a and b as nodes
      for (size_type i=0; i < num_edges(); ++i) {
        //Case1: check if a == node1 and b == node2
        if (a.node_uid_ == edges_[i].node1.node_uid_ && b.node_uid_ == edges_[i].node2.node_uid_) {
          return true;
        } 
        //Case2: check if a == node2 and b == node1
        if (a.node_uid_ == edges_[i].node2.node_uid_ && b.node_uid_ == edges_[i].node1.node_uid_) { 
          return true;
        } 
      }
    }
    //return false if a and b are not within an Edge in edges_
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
  /*** Add new internal_edge to edges_ vector of internal_edges ***/
  Edge add_edge(const Node& a, const Node& b) {
    //instantiate an internal edge to add to edges_ vector
    internal_edge new_edge;
    //assign data attributes to new_edge
    new_edge.node1 = a;
    new_edge.node2 = b;
    new_edge.edge_uid = next_edge_uid_;
    //append new_edge to edges_ vector
    edges_.push_back(new_edge);
    // Update the adjacency list
    adj_list[a.index()].push_back(std::make_pair(b.index(),next_edge_uid_));
    adj_list[b.index()].push_back(std::make_pair(a.index(),next_edge_uid_));
    //update number of edges and next edge ID
    ++num_edges_;
    ++next_edge_uid_;

    return Edge(this, next_edge_uid_-1);       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    /*** clear nodes_ and edges_ vectors ***/ 
    nodes_.clear();
    edges_.clear();
    //reset number of nodes and edges
    num_nodes_ = 0;
    num_edges_ = 0;
    //reset ID numbers
    next_node_uid_ = 0;
    next_edge_uid_ = 0;
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

    /** Return a Node instance for this node ID
     */
    Node operator*() const
    {
      return Node(graph_, node_iter_uid_);
    }

    /** Increment the pointer to the node in the graph by one
     *  to the next node
     */
    NodeIterator& operator++()
    {
      if (node_iter_uid_ < graph_->size())
        ++node_iter_uid_;
      return *this;
    }

    /** Test whether this NodeIterator and a given one are the same
     * Test by comparing if both of them belong to the same graph
     * as well as if they have the same node ID.
     */
    bool operator==(const NodeIterator& ni) const
    {
      return ((graph_ == ni.graph_) && (node_iter_uid_ == ni.node_iter_uid_));
    }

   private:
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_iter_uid_;
    // Private constructor for NodeIterator
    NodeIterator(const Graph* graph, size_type node_iter_uid=0)
      : graph_(const_cast<Graph*>(graph)), node_iter_uid_(node_iter_uid){
      }
    friend class Graph;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return a NodeIterator beginning at the first node of the
   *  graph. This order depends on how nodes were added to the
   *  graph.
   */
  node_iterator node_begin() const
  {
    return NodeIterator(this, 0);
  }

  /** Return a NodeIterator of the last node of the graph. This
   *  order depends on how nodes were added to the graph.
   */
  node_iterator node_end() const
  {
    return NodeIterator(this, size());
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

    /** Return an Edge instance of this edge incident on this node
     */
    Edge operator*() const
    {
      return Edge(graph_, graph_->adj_list[node_uid_][node_nbr_idx_].second);
    }

    /** Increment the pointer to the edge in the graph by one
     *  to the next edge that is incident on the same node
     */
    IncidentIterator& operator++()
    {
      if (node_nbr_idx_ < graph_->adj_list[node_uid_].size())
        ++node_nbr_idx_;
      return *this;
    }

    /** Test whether this IncidentIterator and a given one are the same
     * Test by comparing if both of them belong to the same graph,
     * they have the same node ID, and have the same neighbour index
     */
    bool operator==(const IncidentIterator& iit) const
    {
      return ((graph_ == iit.graph_) && (iit.node_uid_ == node_uid_) && (iit.node_nbr_idx_ == node_nbr_idx_));
    }

   private:
    Graph* graph_;
    size_type node_uid_;
    size_type node_nbr_idx_;
    // Private constructor for Incident Iterator
    IncidentIterator(const Graph* graph, size_type node_uid, size_type node_nbr_idx=0)
      : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid), node_nbr_idx_(node_nbr_idx){
      }
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

    /** Return an Edge instance of this edge
     */
    Edge operator*() const
    {
      return Edge(graph_, edge_iter_uid_);
    }


    /** Increment the pointer to the edge in the graph by one
     *  to the next edge
     */
    EdgeIterator& operator++()
    {
      if (edge_iter_uid_ < graph_->edges_.size())
        ++edge_iter_uid_;
      return *this;
    }

    /** Test whether this EdgeIterator and a given one are the same
     * Test by comparing if both of them belong to the same graph
     * and have the same egde index
     */
    bool operator==(const EdgeIterator& ei) const
    {
      return ((graph_ == ei.graph_) && (edge_iter_uid_ == this->edge_iter_uid_));
    }

   private:
    Graph* graph_;
    size_type edge_iter_uid_;
    EdgeIterator(const Graph* graph, size_type edge_iter_uid=0)
      : graph_(const_cast<Graph*>(graph)), edge_iter_uid_(edge_iter_uid){
      }
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:


  /** Return an EdgeIterator beginning at the first edge of the
   *  graph. This order depends on how edges were added to the
   *  graph.
   */
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, 0);
  }

  /** Return an EdgeIterator of the last edge of the graph. This
   *  order depends on how edges were added to the graph.
   */
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, edges_.size());
  }

 private:
  /*** Graph class internals for Node and Edge proxys ***/
  struct internal_node {
    Point position; 
    size_type node_uid;
    node_value_type value;
  };

  struct internal_edge {
    Node node1;
    Node node2;
    size_type edge_uid;
    edge_value_type value;
  };

  /*** Graph class vectors of internals ***/
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_list;

  /*** Graph class data members ***/ 
  size_type num_nodes_;
  size_type num_edges_;
  size_type next_node_uid_;
  size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP