#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/* Nicole Schiavone
 * CME 212 -  HW1
 */

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

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
  Graph(): nodes_(), edges_(), nodeConnections_(){
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

    /** Return this node's position. */
    const Point& position() const {
      return this->graph_->nodes_[this->nodeIndex_].position_;
    }

    /** Return this node's value, can be modified */
    node_value_type& value() {
      return this->graph_->nodes_[this->nodeIndex_].value_;
    }

    /** Return this node's value, cannot be modified. */
    const node_value_type& value() const {
      return this->graph_->nodes_[this->nodeIndex_].value_;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->nodeIndex_;
    }

    /** Return the number of edges incident to this node */
    size_type degree() const {
      return this->graph_->nodeConnections_[this->nodeIndex_].size();
    }

    /** Return an IncidentIterator at the beginning of the list of
     * edges connected to this node
     */
    IncidentIterator edge_begin() const {
      return IncidentIterator(this->graph_, this->nodeIndex_, 0);
    }

    /** Return an IncidentIterator at the end of the list of
     * edges connected to this node
     */
    IncidentIterator edge_end() const {
      return IncidentIterator(this->graph_, this->nodeIndex_, this->degree());
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->graph_ == n.graph_) && (this->nodeIndex_ == n.nodeIndex_));
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
      if (this->graph_ == n.graph_)
	{
	  return (this->nodeIndex_ < n.nodeIndex_);
	}
      else
	{
	  if (this->nodeIndex_ == n.nodeIndex_)
	    return (this->graph_ < n.graph_);
	  else
	    return (this->nodeIndex_ < n.nodeIndex_);
	}      
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Node belongs to
    size_type nodeIndex_; //The Node index

    // Private Constructor
    Node(const Graph* graph, size_type nodeIndex)
      : graph_(const_cast<Graph*>(graph)), nodeIndex_(nodeIndex)
    {
    }

  };
  //End of Node subclass
  
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->nodes_.size();
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
    // Create new Node with the current number of nodes as the index
    Node addedNode(this, this->num_nodes());
    
    // Add new NodeInfo struct with position and value
    // Increases number of nodes
    NodeInfo n = {position, value};
    this->nodes_.push_back(n);

    //Each Node needs to be initialized in nodeConnections for Edge management
    this->nodeConnections_.push_back(std::vector<std::vector<size_type> >());
    
    return addedNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
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

    // Use constructor to return a Node with specific index
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return this->graph_->node(this->node1Idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return this->graph_->node(this->node2Idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //Edges must be in the same graph to be equal
      if (this->graph_ != e.graph_)
	return false;

      if ((this->node1Idx_ == e.node1()) && this->node2Idx_ == e.node2())
	return true;

      if ((this->node1Idx_ == e.node2()) && this->node2Idx_ == e.node1())
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
	return this->edgeIndex_ < e.edgeIndex_;
      if (this->edgeIndex_ == e.edgeIndex_)
	return this->graph_ < e.graph_;
      else
	return this->edgeIndex_ < e.edgeIndex_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Edge belongs to
    size_type edgeIndex_; //The Edge index
    size_type node1Idx_; //The index of Node 1 of the Edge
    size_type node2Idx_; //The index of Node 2 of the Edge

    // Private Constructor
    Edge(const Graph* graph, size_type edgeIndex, 
	 size_type node1Idx, size_type node2Idx)
      : graph_(const_cast<Graph*>(graph)), edgeIndex_(edgeIndex),
	node1Idx_(node1Idx), node2Idx_(node2Idx)
    {
    }
  };
  //End of Edge subclass

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edges_.size();
    //return this->edgectr;
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

    //Use constructor to return an Edge with a specific index
    size_type n1Idx = this->edges_[i][0];
    size_type n2Idx = this->edges_[i][1];
    return Edge(this,i,n1Idx,n2Idx);
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

    //Iterate through the node connectivity matrix to see if the two nodes
    //are already connected
    for (size_type iter = 0;
	 iter < this->nodeConnections_[a.index()].size(); iter++)
      {
	if (this->nodeConnections_[a.index()][iter][0] == b.index())
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
  Edge add_edge(const Node& a, const Node& b) {
    //Assert preconditions - can only add edges between two nodes in the graph
    //And the nodes must be distinct
    assert(this->has_node(a));
    assert(this->has_node(b));
    assert(!(a == b));

    //Check if the Edge is in the Graph
    //If it is, find the Edge and return it. If not, create a new Edge
    if (has_edge(a,b))
      {
	size_type foundEdgeIdx;
	for (size_type iter = 0;
	     iter < this->nodeConnections_[a.index()].size(); iter++)
	  {
	    if (this->nodeConnections_[a.index()][iter][0] == b.index())
	      {
		foundEdgeIdx = this->nodeConnections_[a.index()][iter][1];
	      }
	  }
	return Edge(this,foundEdgeIdx,a.index(),b.index());
      }	
    else
      {
	//Create the new Edge
	Edge addedEdge(this, this->num_edges(), a.index(), b.index());

	//Add the edge connection to the node connectivity vector
	std::vector<size_type> temp1 = {b.index(),this->num_edges()};
	this->nodeConnections_[a.index()].push_back(temp1);
	std::vector<size_type> temp2 = {a.index(),this->num_edges()};
	this->nodeConnections_[b.index()].push_back(temp2);
	
	// Add the edge to the vector of Node index pairs
	std::vector<size_type> temp3 = {a.index(), b.index()};
	this->edges_.push_back(temp3);
	
	return addedEdge;
      }
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this->nodes_.clear();
    this->edges_.clear();
    this->nodeConnections_.clear();
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
      return Node(this->graph_,this->idx_);
    }

    /** Increment this NodeIterator to the next Node.
     * Returns updated this NodeIterator
     */
    NodeIterator& operator++()
    {
      this->idx_++;
      return *this;
    }

    /** Test whether this NodeIterator and @a ni are equal.
     *
     * Equal NodeIterators have the same graph and the same index
     * The NodeIterator index indicates the Node to which it is pointing
     */
    bool operator==(const NodeIterator& ni) const
    {
      return ((this->graph_ == ni.graph_) && (this->idx_ == ni.idx_));
    }

   private:
    friend class Graph;

    Graph* graph_; //The graph of the NodeIterator
    size_type idx_; //The iterator index of the list of nodes

    // Private Constructor
    NodeIterator(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx)
    {
    }

  };
  // End of NodeIterator subclass

  /** Returns a NodeIterator pointing to the beginning
   * of this graph's nodes list
   */
  NodeIterator node_begin() const
  {
    return NodeIterator(this,0);
  }

  /** Returns a NodeIterator pointing to the end
   * of this graph's nodes list
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
      std::vector<size_type> edgeInfo =
	this->graph_->nodeConnections_[rootNodeIndex_][idx_];

      return Edge(this->graph_, edgeInfo[1], rootNodeIndex_, edgeInfo[0]);
    }

    /** Increment this IncidentIterator to the next Edge connected to the node
     * Returns updated this IncidentIterator
     */
    IncidentIterator& operator++()
    {
      this->idx_++;
      return *this;
    }

    /** Test whether this IncidentIterator and @a ii are equal.
     * 
     * Equal IncidentIterators have the same graph, the same root node, and
     * the same index in the iterations
     */
    bool operator==(const IncidentIterator& ii) const
    {
      if (this->graph_ == ii.graph_)
	{
	  if (this->rootNodeIndex_ == ii.rootNodeIndex_)
	    {
	      if (this->idx_ == ii.idx_)
		return true;
	    }
	}

      return false;
    }

   private:
    friend class Graph;

    Graph* graph_; //Pointer to the root Node's graph
    size_type rootNodeIndex_; //Index of the root node
    size_type idx_; //Iterator index in the list of incident edges

    //Private constuctor
    IncidentIterator(const Graph* graph, size_type rootNodeIndex, size_type idx)
      : graph_(const_cast<Graph*>(graph)), rootNodeIndex_(rootNodeIndex),
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
      size_type n1Idx = this->graph_->edges_[idx_][0];
      size_type n2Idx = this->graph_->edges_[idx_][1];
      return Edge(this->graph_,this->idx_,n1Idx,n2Idx);
    }

    /** Increment this EdgeIterator to the next Edge.
     * Returns updated this EdgeIterator
     */
    EdgeIterator& operator++()
    {
      this->idx_++;
      return *this;
    }

    /** Test whether this EdgeIterator and @a ni are equal.
     *
     * Equal EdgeIterators have the same graph and the same index
     * The EdgeIterator index indicates the Edge to which it is pointing
     */
    bool operator==(const EdgeIterator& ei) const
    {
      return ((this->graph_ == ei.graph_) && (this->idx_ == ei.idx_));
    }

   private:
    friend class Graph;

    Graph* graph_; //The graph of the EdgeIterator
    size_type idx_; //The iterator index in the list of edges

    // Private Constructor
    EdgeIterator(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx)
    {
    }    
  };
  // End of EdgeIterator class

  /** Returns an EdgeIterator pointing to the beginning
   * of this graph's edges list
   */
  EdgeIterator edge_begin() const
  {
    return EdgeIterator(this,0);
  }
  
  /** Returns an EdgeIterator pointing to the end
   * of this graph's nodes list
   */
  EdgeIterator edge_end() const
  {
    return EdgeIterator(this,this->num_edges());
  }
  
 private:
  //A vector of Points that coorespond to the Nodes in the graph
  std::vector<NodeInfo> nodes_;

  //Vector to hold the edges. The inner vector of size_types contains
  //indices for nodes a and b for each edge
  std::vector<std::vector<size_type> > edges_;
  
  //Vector of all Nodes. Each node has an associated vector of vectors.
  //The inner most vector contains the index of a Node connected to the
  //original Node and the edge that connects them.
  //Node A:<  <NodeB1, EdgeIndex1>, <NodeB2, EdgeIndex2> ... >
  std::vector<std::vector<std::vector<size_type> > > nodeConnections_;

  //Struct to hold the point and value of a node
  struct NodeInfo {
    Point position_;
    V value_;
  };  
};

#endif // CME212_GRAPH_HPP
