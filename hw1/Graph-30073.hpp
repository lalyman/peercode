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
template <typename V = double>
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
  using node_value_type = V;
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

    /** Set and return node's value
	 * @param node_value_val of type node_value_type
	 * @return node_value_val after setting it in the proxy class for this node
	 */
    node_value_type& value(node_value_type node_value_val)
    {
      graph_->nodes_[node_uid_].node_value_ = node_value_val;
      return graph_->nodes_[node_uid_].node_value_;
    }
    /** Return node's value
	 * @return node_value_ from the proxy class for this node
	 */
    const node_value_type& value() const
    {
      return graph_->nodes_[node_uid_].node_value_;
    }
    /** Rerturn node's degree
	 * @return degree computed as the size of this node's adjacency list
	 */
    size_type degree() const
    {
      return graph_->adj_list_[node_uid_].size();
    }
    /* @return pointer to the first edge having this point as one endpoint
	 */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(graph_, node_uid_, 0);
    }
    /*@return pointer to the next after last edge having this point as one endpoint
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
      if ((n.index() == this->node_uid_) && (n.graph_ == this->graph_)) {
        return true;
      }
      else {
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
      /*** Compare node uids for ordering purposes ***/
      if (n.index() > this->node_uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    /*** Private data members and methods for Node ***/
    //Pointer back to the Graph container
    Graph* graph_;
    //This node's unique identification number
    size_type node_uid_;
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    //instantiate an internal node to add to nodes_ vector
    internal_node new_node;
    //assign data attributes to new_node
    new_node.position = position;
    new_node.node_uid = next_node_uid_;
    new_node.node_value_ = node_value_type();
    //append new_node to nodes_ vector
    nodes_.push_back(new_node);
    //update number of nodes and next node ID 
    ++num_nodes_;
    ++next_node_uid_;
    //update adjacency list
    std::vector<std::pair<size_type, size_type>> empty_neighbours;
    adj_list_.push_back(empty_neighbours);
    return Node(this, next_node_uid_-1);   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /*** same graph if ptr to memory for node is the same as graph_ ptr ***/ 
    if (this == n.graph_) {
      return true;
    }
    else {
      return false;
    }
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

  //
  // EDGES
  //

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
      return (graph_ == e.graph_ && edge_uid_ == e.edge_uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return edge_uid_ < e.edge_uid_;
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
    //update number of edges and next edge ID
    ++num_edges_;
    ++next_edge_uid_;
    //update adjacency list
    adj_list_[a.index()].push_back(std::make_pair (b.index(), next_edge_uid_ - 1));
    adj_list_[b.index()].push_back(std::make_pair (a.index(), next_edge_uid_ - 1));
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

    /** @return the current node this iterator is pointing to
	 */
    Node operator*() const
    {
      return Node(graph_, node_uid_val_);
    }
    /** Make the iterator point to the next node
	 * @return iterator to the next node
	 */
    NodeIterator& operator++()
    {
      if(node_uid_val_ < graph_->size())
        node_uid_val_ ++;
      return *(this);
    }
    /** Check if two node iterators are equal
	 * @param ni of type NodeIterator
	 * @return true if ni is the same as this node iterator
	 */
    bool operator==(const NodeIterator& ni) const
    {
      return (graph_ == ni.graph_ && node_uid_val_ == ni.node_uid_val_);
    }

   private:
    Graph* graph_;
    size_type node_uid_val_;
    /* private constructor Graph::Edge */
    NodeIterator(const Graph* graph, size_type node_uid_val = 0)
      : graph_(const_cast<Graph*>(graph)), node_uid_val_(node_uid_val) {
    }
    friend class Graph;
  };

  /** @return node iterator to the first node
	*/
  node_iterator node_begin() const
  {
    return NodeIterator(this, 0);
  }
  /** @return node iterator to the node after the last
	*/
  node_iterator node_end() const
  {
    return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    /** @return the edge pointed to by this incident iterator
	 */
    Edge operator*() const
    {
      return Edge(graph_, graph_->adj_list_[node_uid_val_][node_neighbour_index_].second);
    }
    /** @return next incident iterator to the current one
	 */
    IncidentIterator& operator++()
    {
      if(graph_->adj_list_[node_uid_val_].size() > node_neighbour_index_)
        node_neighbour_index_++;
      return *(this);
    }
    /** Check if two incident iterators are equal
	 * @param iitm of type IncidentIterator
	 * @return true if iitm is the same as this incident iterator
	 */
    bool operator==(const IncidentIterator& iitm) const
    {
      return (graph_ == iitm.graph_ && node_uid_val_ == iitm.node_uid_val_ && node_neighbour_index_ == iitm.node_neighbour_index_);
    }

   private:
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_uid_val_;
    size_type node_neighbour_index_;
    /* private constructor Graph::Edge */
    IncidentIterator(const Graph* graph, size_type node_uid_val, size_type node_neighbour_index = 0)
      : graph_(const_cast<Graph*>(graph)), node_uid_val_(node_uid_val), node_neighbour_index_(node_neighbour_index) {
    }
    friend class Graph;
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

    /** @return edge pointed to by this edge iterator
	 */
    Edge operator*() const
    {
      return Edge(graph_, edge_uid_val_);
    }
    /** @return next edge iterator to the current one
	 */
    EdgeIterator& operator++()
    {
      if(edge_uid_val_ < graph_->num_edges())
        edge_uid_val_++;
      return *(this);
    }
    /** Check if two edge iterators are equal
	 * @param eitm of type EdgeIterator
	 * @return true if eitm is the same as this edge iterator
	 */
    bool operator==(const EdgeIterator& eitm) const
    {
      return (graph_ == eitm.graph_ && edge_uid_val_ == eitm.edge_uid_val_);
    }

   private:
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_uid_val_;
    /* private constructor Graph::Edge */
    EdgeIterator(const Graph* graph, size_type edge_uid_val = 0)
      : graph_(const_cast<Graph*>(graph)), edge_uid_val_(edge_uid_val) {
    }
    friend class Graph;
  };

  /** @return edge iterator to the first edge of this graph
	 */
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, 0);
  }
  /** @return edge iterator to the edge after the last edge of this graph
	 */
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, num_edges());
  }

 private:
  /*** Graph class internals for Node and Edge proxys ***/
  struct internal_node {
    Point position; 
    size_type node_uid;
    node_value_type node_value_;
  };

  struct internal_edge {
    Node node1;
    Node node2;
    size_type edge_uid;
  };

  /*** Graph class vectors of internals ***/
  // Node and Edge proxies
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  // Adjacency list
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_list_;

  /*** Graph class data members ***/ 
  size_type num_nodes_;
  size_type num_edges_;
  size_type next_node_uid_;
  size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP
