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
template <typename V = int, typename E = int>
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
  Graph() {
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

    /*** Return this node's position for modification purposes. ***/
    Point& position() {
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
   */
    node_value_type& value()
    {
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
      return (n.node_uid_ == node_uid_) && (n.graph_ == graph_);
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
      if(graph_ == n.graph_)
        return node_uid_ < n.node_uid_;
      return graph_ < n.graph_;
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
    return node_ids_.size();
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
    new_node.node_uid = node_ids_.size();
    new_node.node_value_ = node_val;
    //append node_ids_
    node_ids_.push_back(nodes_.size());
    //append new_node to nodes_ vector
    nodes_.push_back(new_node);
    //update adjacency list
    std::vector<internal_edge> empty_neighbours;
    adj_list_.push_back(empty_neighbours);
    return Node(this, nodes_.size() - 1);   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_) && (node_ids_[n.index()] == n.node_uid_);
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
    return Node(this, node_ids_[i]);
  }

  /** Remove a node from the graph.
   * @param[in] n A valid node to remove
   * @return Index of the node removed
   *
   * @pre @a n is valid node
   * @pre graph.node(@a n.index()) == @a n
   * @pre has_node(@a n) == true
   * @post @a n is invalidated
   * @post new num_nodes() == old num_nodes() - 1
   * @post has_node(@a n) == false
   * @post All edges adjacent to @a n becomes invalid
   * @post new num_edges() == old num_edges() - n.degree()
   * @post Outstanding NodeIterators and EdgeIterators invalidated
   *
   * Complexity: O(max_degree^2)
   */
  size_type remove_node(const Node& n)
  {
    //Remove incident edges
    for(size_type iit = 0; iit < adj_list_[n.node_uid_].size(); ++iit)
    {
      size_type other_node_uid = adj_list_[n.node_uid_][iit].node_uid;
      size_type idx = 0;
      bool found_edge = false;
      while(!found_edge && idx < adj_list_[other_node_uid].size())
      {
        if(adj_list_[other_node_uid][idx].node_uid == n.node_uid_)
        {
          adj_list_[other_node_uid][idx] = adj_list_[other_node_uid].back();
          adj_list_[other_node_uid].pop_back();
          found_edge = true;
        }
        ++idx;
      }
    }
    adj_list_[n.node_uid_].clear();

    size_type idx = n.index();
    nodes_[node_ids_[node_ids_.size() - 1]].node_uid = idx;
    
    //Remove node from node_ids_
    node_ids_[idx] = node_ids_.back();
    node_ids_.pop_back();
    
    return idx;
  }

  /** Remove a node from the graph.
   * @param[in] n_it A valid node iterator to remove
   * @return Index of the node removed
   *
   * @pre @a n = @a *n_it is valid node
   * @pre graph.node(@a n.index()) == @a n
   * @pre has_node(@a n) == true
   * @post @a n is invalidated
   * @post new num_nodes() == old num_nodes() - 1
   * @post has_node(@a n) == false
   * @post All edges adjacent to @a n becomes invalid
   * @post new num_edges() == old num_edges() - n.degree()
   * @post Outstanding NodeIterators and EdgeIterators invalidated
   *
   * Complexity: O(max_degree^2)
   */
  node_iterator remove_node(node_iterator n_it)
  {
    auto node_curr = *n_it;
    size_type idx = remove_node(node_curr);
    return NodeIterator(this, idx);
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
  using edge_value_type = E;
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      /*** Invalid Edge; must be created by calling Graph methods ***/
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node_1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node_2_uid_);
    }

    /*** return length of edge ***/
    double length() const {
      return norm(node1().position() - node2().position()); 
    }

    /*** return value type of edge ***/
    edge_value_type& value()
    {
      size_type edge_idx;
      if(node_1_uid_ < node_2_uid_)
      {
        for(size_type i = 0; i < graph_->adj_list_[node_1_uid_].size(); ++i)
          if(graph_->adj_list_[node_1_uid_][i].node_uid == node_2_uid_)
            edge_idx = i;
        return graph_->adj_list_[node_1_uid_][edge_idx].edge_value_;
      }
      for(size_type i = 0; i < graph_->adj_list_[node_2_uid_].size(); ++i)
        if(graph_->adj_list_[node_2_uid_][i].node_uid == node_1_uid_)
          edge_idx = i;
      return graph_->adj_list_[node_2_uid_][edge_idx].edge_value_;
    }

    /*** return value type of edge ***/
    const edge_value_type& value() const
    {
      return value();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool graph_check = graph_ == e.graph_;
      bool edge_1_check = ((node_1_uid_ == e.node_1_uid_) && (node_2_uid_ == e.node_2_uid_));
      bool edge_2_check = ((node_1_uid_ == e.node_2_uid_) && (node_2_uid_ == e.node_1_uid_));
      return graph_check && (edge_1_check || edge_2_check);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_ == e.graph_)
      {
        if(node_1_uid_ == e.node_1_uid_)
          return node_2_uid_ < e.node_2_uid_;
        return node_1_uid_ < e.node_1_uid_;
      }
      return graph_ < e.graph_;
    }

   private:
    /*** Private data members and methods for Edge ***/
    // Pointer back to the Graph container
    Graph* graph_;
    //The nodes defining this edge
    size_type node_1_uid_;
    size_type node_2_uid_;
    /* private constructor Graph::Edge */
    Edge(const Graph* graph, size_type node_1_uid, size_type node_2_uid)
      : graph_(const_cast<Graph*>(graph)), node_1_uid_(node_1_uid), node_2_uid_(node_2_uid) {
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
    return std::distance(edge_begin(), edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /*** return a proxy object for edge @ i ***/
    assert(i < num_edges());
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    /*** Verify that a and b nodes are witin this graph_ ***/
    if (this == a.graph_ && this == b.graph_)
      //iterate through adjacency list to check if the nodes pair already exists
      for(size_type i = 0; i < adj_list_[a.node_uid_].size(); ++i)
        if(adj_list_[a.node_uid_][i].node_uid == b.node_uid_)
          return true;
    //return false if a and b are not within an Edge
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_val = edge_value_type()) {
    if(has_edge(a, b))
      return Edge(this, a.node_uid_, b.node_uid_);
    //instantiate internal edges to add to adjacency list
    internal_edge new_edge;
    //assign data attributes to internal edges and update adjacency list
    new_edge.node_uid = a.node_uid_;
    new_edge.edge_value_ = edge_val;
    adj_list_[b.node_uid_].push_back(new_edge);
    new_edge.node_uid = b.node_uid_;
    adj_list_[a.node_uid_].push_back(new_edge);
    return Edge(this, a.node_uid_, b.node_uid_);
  }

  /** Remove an edge from the graph.
   * @param[in] n1 Valid node of edge to remove
   * @param[in] n2 Valid node of edge to remove
   * @return Column index on adj_list_ where edge is removed if
   *         has_edge(@a n1, @a n2). Row index is smaller of node uids
   *         Otherwise, degree of node which has smaller uid.
   *
   * @pre @a n1 and @a n2 valid nodes
   * @post new num_edges() == old num_edges() - 1
   * @post has_edge(@a n1, @a n2) == false
   * @post Outstanding EdgeIterators invalidated.
   * @post Outstanding IncidentIterators associated with both @a n1 and
   *       @a n2 invalidated.
   *
   * Complexity: O(max_degree)
   */
  size_type remove_edge(const Node& n1, const Node& n2)\
  {
    if(!has_edge(n1, n2))
      return n1.node_uid_ < n2.node_uid_ ? n1.degree() : n2.degree();
    //Remove the two edges
    size_type idx = 0;
    bool found_edge = false;
    while(!found_edge && idx < adj_list_[n1.node_uid_].size())
    {
      if(adj_list_[n1.node_uid_][idx].node_uid == n2.node_uid_)
      {
        adj_list_[n1.node_uid_][idx] = adj_list_[n1.node_uid_].back();
        adj_list_[n1.node_uid_].pop_back();
        found_edge = true;
      }
      ++idx;
    }
    size_type idx2 = 0;
    found_edge = false;
    while(!found_edge && idx2 < adj_list_[n2.node_uid_].size())
    {
      if(adj_list_[n2.node_uid_][idx2].node_uid == n1.node_uid_)
      {
        adj_list_[n2.node_uid_][idx2] = adj_list_[n2.node_uid_].back();
        adj_list_[n2.node_uid_].pop_back();
        found_edge = true;
      }
      ++idx2;
    }
    return n1.node_uid_ < n2.node_uid_ ? idx : idx2;
  }

  /** Remove an edge from the graph.
   * @param[in] e Valid edge to remove
   * @return Column index on adj_list_ where edge is removed if
   *         has_edge(@a n1, @a n2). Row index is smaller of node uids
   *         Otherwise, degree of node of edge which has smaller uid.
   *
   * @pre @a e is valid edge
   * @post new num_edges() == old num_edges() - 1
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * @post Outstanding EdgeIterators invalidated.
   * @post Outstanding IncidentIterators associated with both @a n1 and
   *       @a n2 invalidated.
   *
   * Complexity: O(max_degree)
   */
  size_type remove_edge(const Edge& e)
  {
    return remove_edge(e.node1(), e.node2());
  }
  
  /** Remove an edge from the graph.
   * @param[in] e_it Valid edge iterator to remove
   * @return Edge iterator to the removed edge
   *
   * @pre e = @a *e_it is valid edge
   * @post new num_edges() == old num_edges() - 1
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * @post Outstanding EdgeIterators invalidated.
   * @post Outstanding IncidentIterators associated with both @a n1 and
   *       @a n2 invalidated.
   *
   * Complexity: O(max_degree)
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    auto edge_curr = *e_it;
    size_type node_neighbour_index = remove_edge(edge_curr.node1(), edge_curr.node2());
    size_type next_node_idx = edge_curr.node1().node_uid_ < edge_curr.node2().node_uid_ ? edge_curr.node1().node_uid_ : edge_curr.node2().node_uid_;
    return EdgeIterator(this, next_node_idx, node_neighbour_index);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    /*** clear nodes_, adj_list_ and edges_ vectors ***/ 
    nodes_.clear();
    node_ids_.clear();
    adj_list_.clear();
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
      if(node_uid_val_ < graph_->size())
        return Node(graph_, graph_->node_ids_[node_uid_val_]);
      return Node(graph_, graph_->nodes_.size());
    }
    /** Make the iterator point to the next node
	 * @return iterator to the next node
	 */
    NodeIterator& operator++()
    {
      if(node_uid_val_ < graph_->size())
        node_uid_val_++;
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
    return NodeIterator(this, size());
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
      return Edge(graph_, node_uid_val_, graph_->adj_list_[node_uid_val_][node_neighbour_index_].node_uid);
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
      return Edge(graph_, node_uid_val_, graph_->adj_list_[node_uid_val_][node_neighbour_index_].node_uid);
    }
    /** @return next edge iterator to the current one
	 */
    EdgeIterator& operator++()
    {
      if(node_uid_val_ < graph_->nodes_.size())
        node_neighbour_index_++;
      find_valid_edge();
      return *(this);
    }
    /** Check if two edge iterators are equal
	 * @param eitm of type EdgeIterator
	 * @return true if eitm is the same as this edge iterator
	 */
    bool operator==(const EdgeIterator& eitm) const
    {
      return (graph_ == eitm.graph_ && node_uid_val_ == eitm.node_uid_val_ && node_neighbour_index_ == eitm.node_neighbour_index_);
    }

   private:
    Graph* graph_;
    size_type node_uid_val_;
    size_type node_neighbour_index_;
    /* private constructor Graph::Edge */
    EdgeIterator(const Graph* graph, size_type node_uid_val = 0, size_type node_neighbour_index = 0)
      : graph_(const_cast<Graph*>(graph)), node_uid_val_(node_uid_val), node_neighbour_index_(node_neighbour_index) {
        find_valid_edge();
    }
    friend class Graph;
    void find_valid_edge()
    {
      while(node_uid_val_ < graph_-> adj_list_.size())
      {
        while(node_neighbour_index_ < graph_-> adj_list_[node_uid_val_].size())
        {
          if(node_uid_val_ < graph_-> adj_list_[node_uid_val_][node_neighbour_index_].node_uid)
            return;
          ++node_neighbour_index_;
        }
        ++node_uid_val_;
        node_neighbour_index_ = 0;
      }
    }
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
    return EdgeIterator(this, adj_list_.size());
  }

 private:
  /*** Graph class internals for Node and Edge proxys ***/
  struct internal_node {
    Point position; 
    size_type node_uid;
    node_value_type node_value_;
  };

  struct internal_edge {
    size_type node_uid;
    edge_value_type edge_value_;
  };

  /*** Graph class vectors of internals ***/
  // Node proxies
  std::vector<internal_node> nodes_;
  std::vector<size_type> node_ids_;
  // Adjacency list
  std::vector<std::vector<internal_edge>> adj_list_;

};

#endif // CME212_GRAPH_HPP
