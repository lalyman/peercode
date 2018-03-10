#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/* Nicole Schiavone
 * CME 212 -  HW3
 */

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

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
  //Declare struct to hold the point and value of a node
  struct NodeInfo;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node value in this graph */
  using node_value_type = V;

  /** Type of edge value in this graph */
  using edge_value_type = E;

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
  Graph(): nodes_(), index2unique_(), adj_(), edgeValues_(){
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
    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position, can be modified */
    Point& position() {
      assert(valid());
      return this->graph_->nodes_[this->uniqueID_].position_;
    }

    /** Return this node's position, cannot be modified */
    const Point& position() const {
      assert(valid());
      return this->graph_->nodes_[this->uniqueID_].position_;
    }

    /** Return this node's value, can be modified */
    node_value_type& value() {
      assert(valid());
      return this->graph_->nodes_[this->uniqueID_].value_;
    }

    /** Return this node's value, cannot be modified. */
    const node_value_type& value() const {
      assert(valid());
      return this->graph_->nodes_[this->uniqueID_].value_;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return this->graph_->nodes_[this->uniqueID_].idx_;
    }

    /** Return the number of edges incident to this node */
    size_type degree() const {
      assert(valid());
      return this->graph_->adj_[this->uniqueID_].size();
    }

    /** Return an IncidentIterator at the beginning of the list of
     * edges connected to this node
     */
    IncidentIterator edge_begin() const {
      assert(valid());
      return IncidentIterator(this->graph_, this->uniqueID_, 0);
    }

    /** Return an IncidentIterator at the end of the list of
     * edges connected to this node
     */
    IncidentIterator edge_end() const {
      assert(valid());
      return IncidentIterator(this->graph_, this->uniqueID_, this->degree());
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same unique ID.
     */
    bool operator==(const Node& n) const {
      assert(valid() && n.valid());
      return ((this->graph_ == n.graph_) && (this->uniqueID_ == n.uniqueID_));
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
      assert(valid());
      if (this->graph_ == n.graph_)
	{
	  return (this->uniqueID_ < n.uniqueID_);
	}
      else
	{
	  if (this->uniqueID_ == n.uniqueID_)
	    return (this->graph_ < n.graph_);
	  else
	    return (this->uniqueID_ < n.uniqueID_);
	}      
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Node belongs to
    size_type uniqueID_; //The unique identifier of the node

    // Private Constructor
    Node(const Graph* graph, size_type uniqueID)
      : graph_(const_cast<Graph*>(graph)), uniqueID_(uniqueID)
    {
    }

    //Private method to check if a node is valid
    //by checking representation invariants
    bool valid() const {
      return uniqueID_ >= 0
	&& uniqueID_ < graph_->nodes_.size()
	&& graph_->nodes_[uniqueID_].idx_ < graph_->index2unique_.size()
	&& graph_->index2unique_[graph_->nodes_[uniqueID_].idx_] == uniqueID_;
    }

  }; //End of Node subclass
  
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->index2unique_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
		const node_value_type& value = node_value_type()) {
    size_type addedNodeIndex = this->num_nodes();
    size_type addedNodeUniqueID = this->nodes_.size();
    
    NodeInfo n = {position,value,addedNodeIndex};
    this->nodes_.push_back(n); //increases total number of nodes ever
    this->index2unique_.push_back(addedNodeUniqueID); //increases num_nodes

    //Each Node needs to be initialized in nodeConnections for Edge management
    this->adj_.push_back(std::vector<size_type>());

    Node addedNode = Node(this,addedNodeUniqueID);
    assert(addedNode.valid());
    return addedNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    assert(n.valid());
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    //Assert preconditions - can't use node method to add a new node
    assert(i >= 0);
    assert(i < num_nodes());

    // The current index i maps to a unique node identifier based on
    // the index2unique vector
    Node returnNode = Node(this, this->index2unique_[i]);
    assert(returnNode.valid());
    assert(this->nodes_[this->index2unique_[i]].idx_ == i);
    return returnNode;
  }

   /** Removes @a n and returns the number of nodes removed
   * (which is one if remove is successful). Also removes all edges
   * incident to @a n
   *
   * @pre @a n must be a valid node in this graph (checked by assert)
   * @post @a n is invalidated and removed from the graph
   * @post All previous NodeIterators, EdgeIterators, and IncidentIterators
   *       are invalidated
   * @post new num_nodes() = old num_nodes()-1
   * @post new num_edges() = old num_edges - @a n.degree()
   *
   * Complexity: Only removing the node takes O((num_nodes())), assuming the
   * Graph is sparse implying that O(n.degree()) operations are much less
   * significant in complexity
   *
   * Justification for Complexity:
   *  Operations to remove the node:
   *     Erases @n from the list of valid nodes: O(num_nodes()) at worst
   *     Updates non-static indices of all valid nodes: O(num_nodes())
   *  Also has to remove all edges connected to @a n: as meets the provided
   *     specifications, remove_edge takes O(log(num_edges())) so removing
   *     all incident edges is O(n.degree())*O(log(num_edges())) which is
   *     approximately O(log(num_edges())) since the Graph is sparse 
   */
  size_type remove_node(const Node& n) {
    assert(n.valid());
    size_type removeIdx = n.index();
    size_type n_uid = this->index2unique_[n.index()];

    //Remove connected edges
    size_type numEdgesToRemove = n.degree();
    std::vector<Edge> edgesToRemove;
    for (size_type iter = 0; iter < numEdgesToRemove; iter++)
      {
	Edge e = Edge(this,n_uid,this->adj_[n_uid][iter]);
	edgesToRemove.push_back(e);
      }
    for (size_type iter = 0; iter < numEdgesToRemove; iter++)
      {
	size_type result = remove_edge(edgesToRemove[iter]);
	(void) result;
      }
    
    //Remove the node by removing it from index2unique
    auto result =
      this->index2unique_.erase(this->index2unique_.begin()+removeIdx);
    (void) result;

    //Iterate through nodes_ to update NodeInfo struct indices
    for (auto nsi = this->nodes_.begin(); nsi != this->nodes_.end(); ++nsi)
      {
	if ((*nsi).idx_ > removeIdx)
	  (*nsi).idx_--;
      }    
    return 1;
  }

  /** Removes the node pointed to by @a n_it and returns a node_iterator
   * pointing to the new node_end() of this graph. Also removes all edges
   * incident to @a *n_it
   *
   * @pre @a *n_it must be a valid edge in this graph
   * @post @a *n_it is invalidated and removed from the graph
   * @post All previous NodeIterators, EdgeIterators, and IncidentIterators
   *       are invalidated
   * @post new num_nodes() = old num_nodes()-1
   * @post new num_edges() = old num_edges - @a n.degree()
   *
   * Complexity: Only removing the node takes O((num_nodes())), assuming the
   * Graph is sparse implying that O(n.degree()) operations are much less
   * significant in complexity
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type removeIdx = remove_node(*n_it);
    (void) removeIdx; //Quiet compiler warning
    return this->node_end();
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      size_type node1Idx = this->graph_->nodes_[this->node1UniqueID_].idx_;
      return this->graph_->node(node1Idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type node2Idx = this->graph_->nodes_[this->node2UniqueID_].idx_;
      return this->graph_->node(node2Idx);
    }

    /** Return this Edge's value, can be modified */
    edge_value_type& value() {
      size_type n1_uid = this->node1UniqueID_;
      size_type n2_uid = this->node2UniqueID_;
      std::pair<size_type,size_type> edgeNodes =
	this->graph_->createEdgeKey(n1_uid, n2_uid);

      return this->graph_->edgeValues_[edgeNodes];
    }

    /** Return this Edge's value, cannot be modified. */
    const edge_value_type& value() const {
      size_type n1_uid = this->node1UniqueID_;
      size_type n2_uid = this->node2UniqueID_;
      std::pair<size_type,size_type> edgeNodes =
	this->graph_->createEdgeKey(n1_uid, n2_uid);

      return this->graph_->edgeValues_[edgeNodes];
    }

    /** Return the length of this Edge */
    double length() const {
      return norm(this->node1().position() - this->node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //Edges must be in the same graph to be equal
      if (this->graph_ != e.graph_)
	return false;

      size_type e_n1_uid = this->graph_->index2unique_[(e.node1()).index()];
      size_type e_n2_uid = this->graph_->index2unique_[(e.node2()).index()];

      if (this->node1UniqueID_ == e_n1_uid &&
	  this->node2UniqueID_ == e_n2_uid)
	return true;

      if (this->node1UniqueID_ == e_n2_uid &&
	  this->node2UniqueID_ == e_n1_uid)
	return true;

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->graph_ == e.graph_)
	{
	  if (this->node1UniqueID_ == e.node1UniqueID_)
	    return this->node2UniqueID_ < e.node2UniqueID_;
	  else
	    return this->node1UniqueID_ < e.node1UniqueID_;
	}

      if (this->node1UniqueID_ == e.node1UniqueID_ &&
	  this->node2UniqueID_ == e.node2UniqueID_)
	return this->graph_ < e.graph_;
      else
	{
	  if (this->node1UniqueID_ == e.node1UniqueID_)
	    return this->node2UniqueID_ < e.node2UniqueID_;
	  else
	    return this->node1UniqueID_ < e.node1UniqueID_;
	} 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Edge belongs to
    size_type node1UniqueID_; //The unique identifier of Node 1 of the Edge
    size_type node2UniqueID_; //The unique identifier of Node 2 of the Edge

    // Private Constructor
    Edge(const Graph* graph,
	 size_type node1UniqueID, size_type node2UniqueID)
      : graph_(const_cast<Graph*>(graph)),
	node1UniqueID_(node1UniqueID), node2UniqueID_(node2UniqueID)
    {
    }
  };
  //End of Edge subclass

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edgeValues_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    //Assert preconditions - can't use edge method to add a new edge
    assert(i >= 0);
    assert(i < num_edges());

    return *(this->edge_begin()+i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //Assert preconditions - can only test for edges
    //between two nodes in the graph
    assert(this->has_node(a));
    assert(this->has_node(b));

    //Use incident iterator to see if these two nodes are connected
    for (auto ii = a.edge_begin(); ii != a.edge_end(); ++ii)
      {
	if (b == (*ii).node2())
	  return true;
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
  Edge add_edge(const Node& a, const Node& b,
		const edge_value_type& value = edge_value_type()) {
    //Assert preconditions - can only add edges between two nodes in the graph
    //and the nodes must be distinct
    assert(this->has_node(a));
    assert(this->has_node(b));
    assert(!(a == b));
    
    //Get the unique identifiers of the nodes
    size_type a_uid = this->index2unique_[a.index()];
    size_type b_uid = this->index2unique_[b.index()];
    
    //Return the Edge if it is already in the Graph
    for (auto ii = a.edge_begin(); ii != a.edge_end(); ++ii)
      {
	if (b == (*ii).node2())
	  {
	    return Edge(this,a_uid,b_uid);
	  }
      }
    
    //Edge is not in Graph, so create Edge
    Edge addedEdge(this, a_uid, b_uid);
    
    //Add the edge to the adjacency vector
    this->adj_[a_uid].push_back(b_uid);
    this->adj_[b_uid].push_back(a_uid);

    //Add the edge value to the edge value map
    std::pair<size_type,size_type> edgeNodes = createEdgeKey(a_uid, b_uid);
    this->edgeValues_.insert(std::make_pair(edgeNodes,value));
    return addedEdge;
  }

  /** Removes Edge @a e and returns the number of edges removed
   * (which is one if remove is successful)
   *
   * @pre @a e must be a valid edge in this graph
   *
   * @post @a e is invalidated and removed from the graph
   * @post All previous EdgeIterators and IncidentIterators are invalidated
   * @post new num_edges() = old num_edges()-1
   *
   * Complexity: No more than O(log(num_edges())), assuming the Graph is sparse,
   * implying that O(n.degree()) operations are much less significant
   * in complexity
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1();
    Node n2 = e.node2();

    assert(n1.valid());
    assert(n2.valid());
    
    return remove_edge(n1,n2);
  }

  /** Removes the edge connected by Nodes @a a and @a b and returns the number
   * of edges removed (which is one if remove is successful)
   *
   * @pre @a Edge(this,@a a.uniqueID,@a b.uniqueID) must be a valid edge
   * in this graph
   *
   * @post @a Edge(this,@a a.uniqueID,@a b.uniqueID) and
   *       @a Edge(this,@a b.uniqueID,@a a.uniqueID) are invalidated
   *       and removed from the graph
   * @post All previous EdgeIterators and IncidentIterators are invalidated
   * @post new num_edges() = old num_edges()-1
   *
   * Complexity: No more than O(log(num_edges())), assuming the Graph is sparse,
   * implying that O(n.degree()) operations are much less significant
   * in complexity
   *
   * Justification for complexity:
   *   Searches through two node adjacency lists: O(@a a.degree + @a b.degree)
   *   Erases edge from these two lists: O(@a a.degree + @a b.degree) worst case
   *   Erases edge from the map of all values: O(log(num_edges()))
   */
  size_type remove_edge(const Node& a, const Node& b)
  {
    assert(a.valid());
    assert(b.valid());

    size_type a_uid = a.uniqueID_;
    size_type b_uid = b.uniqueID_;

    size_type origNumEdges = this->num_edges();
    std::pair<size_type,size_type> edgeNodes = createEdgeKey(a_uid, b_uid);
    this->edgeValues_.erase(edgeNodes);
    if (this->num_edges() == origNumEdges)
      return 0; //Edge was not in the graph

    //Erase from the Node a adjacency
    size_type numEdgesToSearch_a = a.degree();
    for (size_type iter = 0; iter < numEdgesToSearch_a; iter++)
      {
	if (b_uid == this->adj_[a_uid][iter])
	  {
	    this->adj_[a_uid].erase(this->adj_[a_uid].begin()+iter);
	  }
      }

    //Erase from the Node b adjacency
    size_type numEdgesToSearch_b = b.degree();
    for (size_type iter = 0; iter < numEdgesToSearch_b; iter++)
      {
	if (a_uid == this->adj_[b_uid][iter])
	  {
	    this->adj_[b_uid].erase(this->adj_[b_uid].begin()+iter);
	  }
      }

    return 1;
  }

  /** Removes the edge pointed to by @a e_it and returns an edge_iterator
   * pointing to the new edge_end() of this graph
   *
   * @pre @a *e_it must be a valid edge in this graph
   * @post @a *e_it is invalidated and removed from the graph
   * @post All previous EdgeIterators and IncidentIterators are invalidated
   * @post new num_edges() = old num_edges()-1
   *
   * Complexity: No more than O(log(num_edges())), assuming the Graph is sparse
   * implying that O(n.degree()) operations are much less significant
   * in complexity
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    size_type result = remove_edge(*e_it);
    (void) result;
    return this->edge_end();
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this->nodes_.clear();
    this->index2unique_.clear();
    this->adj_.clear();
    this->edgeValues_.clear();
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

    /** Return the current Node of this NodeIterator */
    Node operator*() const
    {
      Node returnNode = Node(this->graph_,this->graph_->index2unique_[NIidx_]);
      assert(returnNode.valid());
      assert(this->graph_->nodes_[this->graph_->index2unique_[NIidx_]].idx_
	     == NIidx_);
 
      return returnNode;
    }

    /** Increment this NodeIterator to the next Node.
     * Returns updated this NodeIterator
     */
    NodeIterator& operator++()
    {
      this->NIidx_++;
      return *this;
    }

    /** Returns a copy of this NodeIterator that has been incremented @a iter
     * times
     */
    NodeIterator& operator+(size_type iter)
    {
      auto thiscopy = this;
      for (size_type i = 0; i < iter; i++)
	++(*thiscopy);
      return *thiscopy;
    }
    
    /** Test whether this NodeIterator and @a ni are equal.
     *
     * Equal NodeIterators have the same graph and the same index
     * The NodeIterator index indicates the Node to which it is pointing
     */
    bool operator==(const NodeIterator& ni) const
    {
      return ((this->graph_ == ni.graph_) && (this->NIidx_ == ni.NIidx_));
    }

   private:
    friend class Graph;

    Graph* graph_; //The graph of the NodeIterator
    size_type NIidx_; //The iterator index in the list of valid nodes

    // Private Constructor
    NodeIterator(const Graph* graph, size_type NIidx)
      : graph_(const_cast<Graph*>(graph)), NIidx_(NIidx)
    {
    }

  };
  // End of NodeIterator subclass

  /** Returns a NodeIterator pointing to the beginning
   * of this graph's valid nodes list
   */
  NodeIterator node_begin() const
  {
    return NodeIterator(this,0);
  }

  /** Returns a NodeIterator pointing to the end
   * of this graph's valid nodes list
   */
  NodeIterator node_end() const
  {
    return NodeIterator(this,this->num_nodes());
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
    
    /** Return the current Edge of this IncidentIterator */
    Edge operator*() const
    {
      size_type neighborUid = this->graph_->adj_[rootNodeUid_][idx_];
      return Edge(graph_, rootNodeUid_, neighborUid);
    }

    /** Increment this IncidentIterator to the next Edge connected to the node
     * Returns updated this IncidentIterator
     */
    IncidentIterator& operator++()
    {
      this->idx_++;
      return *this;
    }

    /** Returns a copy of this EdgeIterator that has been incremented @a iter
     * times
     */
    IncidentIterator& operator+(size_type iter)
    {
      auto thiscopy = this;
      for (size_type i = 0; i < iter; i++)
	++(*thiscopy);
      return *thiscopy;
    }

    /** Test whether this IncidentIterator and @a ii are equal.
     * 
     * Equal IncidentIterators have the same graph, the same root node, and
     * the same index in the iterations
     */
    bool operator==(const IncidentIterator& ii) const
    {
      if ((this->graph_ == ii.graph_) &&
	  (this->rootNodeUid_ == ii.rootNodeUid_) &&
	  (this->idx_ == ii.idx_))
	return true;

      return false;
    }

   private:
    friend class Graph;

    Graph* graph_; //Pointer to the root Node's graph
    size_type rootNodeUid_; //Unique identifier of the root node
    size_type idx_; //Iterator index in the list of incident edges of root Node

    //Private constuctor
    IncidentIterator(const Graph* graph, size_type rootNodeUid, size_type idx)
      : graph_(const_cast<Graph*>(graph)), rootNodeUid_(rootNodeUid),
	idx_(idx)
    {
    }
  }; // End of IncidentIterator subclass

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

    /** Return the current Edge of this EdgeIterator */
    Edge operator*() const
    {
      size_type neighborUid = this->graph_->adj_[currNodeUid_][currEdgeIdx_];
      return Edge(graph_, currNodeUid_, neighborUid);
    }

    /** Increment this EdgeIterator to the next Edge.
     * Returns updated this EdgeIterator
     */
    EdgeIterator& operator++()
    {
      //Increment EdgeIdx to the next edge we haven't seen yet,
      //which means that the current Node is less than its neighbor
      //If EdgeIdx reaches the end of the current Nodes adjacent node list,
      //increment to the next Node that has connected edges 
      do
	{
	  currEdgeIdx_++;
	  if (currEdgeIdx_ == graph_->adj_[currNodeUid_].size())
	    {
	      currEdgeIdx_ = 0;
	      //Increment the Node, but check for empty Node adj vectors
	      do
		{
		  currNodeUid_++; //unique IDs will always be consecutive
		  if (currNodeUid_ == graph_->adj_.size())
		    return *this;
		}
	      while (graph_->adj_[currNodeUid_].empty());
	    }
	}
      while (currNodeUid_ > graph_->adj_[currNodeUid_][currEdgeIdx_]);
      
      return *this;
    }

    /** Returns a copy of this EdgeIterator that has been incremented @a iter
     * times
     */
    EdgeIterator& operator+(size_type iter)
    {
      auto thiscopy = this;
      for (size_type i = 0; i < iter; i++)
	++(*thiscopy);
      return *thiscopy;
    }
    
    /** Test whether this EdgeIterator and @a ni are equal.
     *
     * Equal EdgeIterators have the same graph, the same current root Node,
     * and the same current Edge index
     */
    bool operator==(const EdgeIterator& ei) const
    {
      return ((this->graph_ == ei.graph_) &&
	      (this->currNodeUid_ == ei.currNodeUid_) &&
	      (this->currEdgeIdx_ == ei.currEdgeIdx_));
    }

   private:
    friend class Graph;

    Graph* graph_; //The graph of the EdgeIterator
    size_type currNodeUid_; //Current root node
    size_type currEdgeIdx_; //Current edge connected to the current root node

    // Private Constructor
    EdgeIterator(const Graph* graph, size_type currNodeUid,
		 size_type currEdgeIdx)
      : graph_(const_cast<Graph*>(graph)), currNodeUid_(currNodeUid),
	currEdgeIdx_(currEdgeIdx)
    {
    }    
  };
  // End of EdgeIterator class

  /** Returns an EdgeIterator pointing to the beginning of the valid edges
   * in this graph
   */
  EdgeIterator edge_begin() const
  {
    //Find first non-empty adjacency vector
    size_type beginNodeUid = 0;
    while (this->adj_[beginNodeUid].empty())
      beginNodeUid++;

    return EdgeIterator(this,beginNodeUid,0);
  }
  
  /** Returns an EdgeIterator pointing to the end of the valid edges
   * in this graph
   */
  EdgeIterator edge_end() const
  {
    return EdgeIterator(this,this->adj_.size(),0);
  }
  
 private:
  //A vector of NodeInfo structs that coorespond to the Nodes in the graph
  std::vector<NodeInfo> nodes_;

  //A vector of unique identifiers for the Nodes
  //Maps indices to unique identifiers
  std::vector<size_type> index2unique_;

  //Adjacency vector that is indexed by unique Node identifiers. Each unique
  //Node identifier has a list of the unique Node identifiers it's connected to
  std::vector<std::vector<size_type>> adj_;

  //Map containing the Edge values, using the pair of the unique Node
  //identifiers that make up the Edge as the key
  std::map<std::pair<size_type,size_type>,E> edgeValues_;

  //Struct to hold the point, value, and a non-static index of a Node
  //The index member allows for Node lookup by index in constant time
  struct NodeInfo {
    Point position_;
    V value_;
    size_type idx_;
  };

  //Helper function to produce the unique key for the edgeValues map given
  //two Nodes unique identifiers
  std::pair<size_type,size_type> createEdgeKey(size_type a_uid, size_type b_uid)
  {
    std::pair<size_type,size_type> edgeNodes;
    if (a_uid < b_uid)
      edgeNodes = std::make_pair(a_uid,b_uid);
    else
      edgeNodes = std::make_pair(b_uid,a_uid);

    return edgeNodes;
  }
};

#endif // CME212_GRAPH_HPP
