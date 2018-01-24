#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** Nicole Schiavone
 * CME 212
 * HW0
 */

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
class Graph {
 private:
  // No internal types used for HW0

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
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graph_->nodes_[this->nodeIndex_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->nodeIndex_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Check for graph matching, then index matching
      if (this->graph_ == n.graph_)
	{
	  if (this->nodeIndex_ == n.nodeIndex_)
	      return true;
	  else
	      return false;
	}
      else
	  return false;
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
      return (this->nodeIndex_ < n.nodeIndex_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Node belongs to
    size_type nodeIndex_; //The integer index of the Node

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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // Create new Node with the current number of nodes as the index
    Node newNode(this, this->num_nodes());
    
    // Add new Node Point to vector of node positions, increases number of nodes
    this->nodes_.push_back(position);

    //Each Node needs to be initialized in nodeConnections for Edge management
    this->nodeConnections_.push_back(std::vector<std::vector<size_type> >());
    
    return newNode;
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
    Node specificNode(this, i);
    return specificNode;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return this->graph_->edges_[this->edgeIndex_][0];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return this->graph_->edges_[this->edgeIndex_][1];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //Edges must be in the same graph to be equal
      if (this->graph_ != e.graph_)
	return false;

      //Check if Nodes match - order does not matter
      if ((this->node1() == e.node1()) &&
	  this->node2() == e.node2())
	return true;

      if ((this->node1() == e.node2()) &&
	  this->node2() == e.node1())
	return true;

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return  (this->edgeIndex_ < e.edgeIndex_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_; //Pointer to the Graph that Edge belongs to
    size_type edgeIndex_; //The integer index of the Node

    // Private Constructor
    Edge(const Graph* graph, size_type edgeIndex)
      : graph_(const_cast<Graph*>(graph)), edgeIndex_(edgeIndex)
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
    Edge specificEdge(this, i);
    return specificEdge;
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
	Edge foundEdge;
	for (size_type iter = 0;
	     iter < this->nodeConnections_[a.index()].size(); iter++)
	  {
	    if (this->nodeConnections_[a.index()][iter][0] == b.index())
	      {
		foundEdge =
		  this->edge(this->nodeConnections_[a.index()][iter][1]);
	      }
	  }
	return foundEdge;
      }	
    else
      {
	//Create the new Edge
	Edge newEdge(this, this->num_edges());

	//Add the edge connection to the node connectivity vector
	std::vector<size_type> temp1 = {b.index(),this->num_edges()};
	this->nodeConnections_[a.index()].push_back(temp1);
	std::vector<size_type> temp2 = {a.index(),this->num_edges()};
	this->nodeConnections_[b.index()].push_back(temp2);
	
	// Add the edge to the vector of Node index pairs
	std::vector<Node> temp3 = {a,b};
	this->edges_.push_back(temp3);

	return newEdge;
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

 private:
  //A vector of Points that coorespond to the Nodes in the graph
  std::vector<Point> nodes_;

  //Vector to hold the edges. The inner vector of Nodes contains
  //Node a and b for each edge
  std::vector<std::vector<Node> > edges_;

  //Vector of all Nodes. Each node has an associated vector of vectors.
  //The inner most vector contains the index of a Node connected to the
  //original Node and the edge that connects them.
  std::vector<std::vector<std::vector<size_type> > > nodeConnections_;
};

#endif // CME212_GRAPH_HPP
