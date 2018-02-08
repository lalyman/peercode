#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. Edges are unique (there is at
 *  most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  /** Use this space for declarations of important internal types you need
      later in the Graph's definition. */
  
  // Predeclaration of struct to hold NodeInfo
  struct NodeInfo;

  // Predeclaration of function to convert node indices to edge index
  unsigned long node2Edge();
  // Predeclaration of struct to hold both node indices of the endpoints of an edge
  struct edgeEndpts;
  // Predeclaration of function to convert edge index to node indices of its endpoints
  edgeEndpts edge2Node();

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // Type of this graph.
  using graph_type = Graph;

  /** Type of indexes and sizes.
   *  Return type of Graph::node2Edge() and Graph::edge2Node(), Graph::Node::index(), 
   *  Graph::num_nodes(), Graph::num_edges(), and argument type of Graph::node(size_type) 
   */
  using size_type = unsigned;
  
  /** Larger size needed for edge_index hash values to prevent overflow and ensure invertible
   *  mapping (see node2Edge() and edge2Node() )*/
  using size_type_long = unsigned long;

  // Predeclaration of Node type.
  class Node;
  /// Synonym for Node (following STL conventions).
  using node_type = Node;
  /// Synonym for Node Value Type, specified since Graph is a template
  using node_value_type = V;

  /// Predeclaration of Edge type.
  class Edge;
  /// Synonym for Edge (following STL conventions).
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


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  // Construct an empty graph
  Graph() : nodeLocations_(), nodeAdjacency_(), edgeHashes_() {
  }

  // Default destructor
  ~Graph() = default;


  //
  // NODES
  //

  /** @class Graph::Node
   *  @brief Class representing the graph's nodes.
   *
   *  Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>{

   public:

    // Construct an invalid/uninitialized node.
    Node() {
    }

    // Return this node's position.
    const Point& position() const {

      /** position is a member of the struct stored in the Graphs's 
       *  nodeLocations_ array indexed by the Node's index. 
       *  Here, ``this'' refers to the Node. */
      return this->graph_->nodeLocations_[this->index_].pos_;
    }

    // Return this node's index, a number in the range [0, graph_size).
    size_type index() const {

      return this->index_;
    }

    // Return this node's value, as a reference (which could be changed).
	node_value_type& value() {

      /** value is a member of the struct stored in the Graphs's 
       *  nodeLocations_ array indexed by the Node's index. 
       *  Here, ``this'' refers to the Node. */
	  return this->graph_->nodeLocations_[this->index_].val_;
	}
    
    // Return this node's value, as a reference (which CANNOT be changed).
    const node_value_type& value() const {

      return this->graph_->nodeLocations_[this->index_].val_;
    }   

    // Return the degree of this node, i.e. how many incident edges it has
    size_type degree() const {

      /** Find the size of the vector of adjacent nodes stored at
       *  node's index in the vector nodeAdjacency_*/
      return this->graph_->nodeAdjacency_[this->index_].size();
    }

    /** Return the IncidentIterator to start iterating over the incident edges */
    IncidentIterator edge_begin() const {

      // Start iterating with the first edge
  	  IncidentIterator startingEdge = IncidentIterator(this->index_, 0, this->graph_);
      return startingEdge;
    }

    /** Return the IncidentIterator to end iterating over the incident edges */
    IncidentIterator edge_end() const{

      // Stop iterating with the last edge
  	  IncidentIterator endingEdge = IncidentIterator(this->index_, this->degree(), this->graph_);
      return endingEdge;
    }

    /** Test whether this node and @a n are equal.
     *  @param[in] @a n Node to test.
     *
     *  Equal nodes have the same graph and the same index. */
    bool operator==(const Node& n) const {

      // Check if node n and this node belong to same graph
      if(this->graph_ != n.graph_) {
        return false;  // if not, they aren't equal
      }
      // Check if node n and this node have the same index
      if(this->index_ != n.index_) {
        return false;  // if not, they aren't equal
      }
      return true;  // otherwise, they are equal
    }

    /** Test whether this node is less than @a n in a global order.
     *  @param[in] @a n Node to test.
     *
     *  The node ordering relation must obey trichotomy: For any two nodes x
     *  and y, exactly one of x == y, x < y, and y < x is true. 
     *
     *  Global ordering of graphs as follows: for graphs a and b created
     *  in that order, global ordering as follows, with (g,n) tuple denoting
     *  (graph letter, node index):
     *  (a,0), (a,1), ..., (a,size(a)-1), (b,0), ..., (b,size(b)-1) */
    bool operator<(const Node& n) const {

      // Check if node n and this node belong to same graph
      if(this->graph_ != n.graph_) {
      	/** for nodes in different graphs, use std::operator< to 
      	 *  return a comparison of the graph pointer values */
        return this->graph_ < n.graph_;
      }else {
      	// for nodes in the same graph, compare indices
	    return this->index_ < n.index_;
      }
    }

   private:

    // Allow Graph to access Node's private member data and functions
    friend class Graph;
  
    // Store the index of this node
    size_type index_;

    // Store a pointer to the graph this node belongs to
    Graph* graph_;

    // Constructor for Node (valid) 
    Node(size_type index, const Graph* graph)
     : index_(index), graph_(const_cast<Graph*>(graph)) {

     }
  };

  /** Return the number of nodes in the graph.
   *  Complexity: O(1). */
  size_type size() const {

    /** nodeLocations_ holds the information for all nodes in graph,
        so size() returns the number of nodes */
    return this->nodeLocations_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {

    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @parma[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {

    /// Construct a new node
    Node newNode = Node(this->num_nodes(), this);

    // Construct a new NodeInfo struct to hold position and value
    NodeInfo newNodeInfo = {position, value};

    // add new NodeInfo struct to nodeLocations_ at index node_index 
    this->nodeLocations_.push_back(newNodeInfo);
    
    /** initialize nodeAdjacency[this->num_nodes] with empty vector for later use
        storing adjacencies of each node */
    this->nodeAdjacency_.push_back(std::vector<size_type>());

    return newNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    // Check if the Graph of node n is the same as this Graph
    if(n.graph_== this) {
      return true;  // if so return true
    }

    return false;  // otherwise return false
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {

    // assert that node index is correct
  	assert( (i >=0) && (i <= this->num_nodes()) );
    
    /** Create new node with index i, return this. This can be used along with
        the Graph's nodeLocations_ array to determine the Point position */
    Node newNode = Node(i,this);

    return newNode; 
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
  class Edge : private totally_ordered<Edge>{

   public:
    
    //Construct an invalid/uninitialized Edge.
    Edge() {
    }

    // Return a node of this Edge
    Node node1() const {

      // Get node with index node1_
      Node newNode = this->graph_->node(this->node1_);

      return newNode; 
    }

    // Return the other node of this Edge
    Node node2() const {
      
      // Get node with index node2_
      Node newNode = this->graph_->node(this->node2_);

      return newNode;
    }

    /** Test whether this edge and @a e are equal.
     *
     *  Equal edges represent the same edge between two nodes, and will have the same
     *  edge index given by node2Edge(node1, node2).
     *
     */
    bool operator==(const Edge& e) const {

      // Check if edge e and this edge belong to same graph
      if(this->graph_ != e.graph_) {
      	
      	/// edges in different graphs aren't equal
        return false;
      } else {

	    // For edges in the same graph, compare hashed indices
        
        // Get index of current edge
        size_type_long thisEdgeIndex = this->graph_->node2Edge(this->node1_,this->node2_);

        // Get index other edge e
        size_type_long eEdgeIndex = e.graph_->node2Edge(e.node1_,e.node2_);
	    
	    return thisEdgeIndex < eEdgeIndex;
      }
    }
    
    /** Test whether this edge is less than @a e in a global order.
	 *
     *  Global ordering of graphs as follows: for graphs a and b created
     *  in that order, global ordering as follows, with (g,e) tuple denoting
     *  (graph letter, edge index):
     *  (a,min_index), ..., (a,max_index), (b,min_index), ..., (b,max_index) */
    bool operator<(const Edge& e) const {

      // Check if edge e and this edge belong to same graph
      if(this->graph_ != e.graph_) {

      	/** for edges in different graphs, use std::operator< to 
      	 *  return a comparison of the graph pointer values */
        return this->graph_ < e.graph_;
      } else {

	    // For edges in the same graph, compare hashed indices

	    // Get index of current edge
        size_type_long thisEdgeIndex = this->graph_->node2Edge(this->node1_,this->node2_);

        // Get index other edge e
        size_type_long eEdgeIndex = e.graph_->node2Edge(e.node1_,e.node2_);
	    
	    return thisEdgeIndex < eEdgeIndex;
      }
    }
    
   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Store the index of node1 of this edge
    size_type node1_;

    // Store the index of node2 of this edge
    size_type node2_;


    // Store a pointer to the graph this edge belongs to
    Graph* graph_;

    // Constructor for Edge (valid)
    Edge(size_type node1, size_type node2, const Graph* graph)
     : node1_(node1), node2_(node2), graph_(const_cast<Graph*>(graph)) {
     }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {

    // edgeHashes_ stores the hashed indices of all the edges in the graph
    return this->edgeHashes_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {

  	// assert edge index is valid
  	assert( (i >= 0) && (i < this->num_edges()) );

  	// Get node index endpoints from the hashed index for this global edge index
  	edgeEndpts ee = this->edge2Node(edgeHashes_[i]); 

    // Create new edge with these endpoints and return it.
    Edge newEdge = Edge(ee.a, ee.b, this);

    return newEdge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

  	// assert that a and b are valid nodes of the graph
  	assert( this->has_node(a) && this->has_node(b) );

    // get node index of both endpoint nodes
    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // Search through nodeAdjacency_ vector of tuples to see if b and a are neightbors
    for(unsigned i=0; i < this->nodeAdjacency_[a_ind].size(); i++) {
      size_type adj_node = this->nodeAdjacency_[a_ind][i];

      // if the adjcacent node index is the same as b's index, then there is an edge
      if(adj_node == b_ind){
        return true;
      }
    }
    return false;  // no edge
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   *  @pre @a a and @a b are distinct valid nodes of this graph
   *  @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *  @post has_edge(@a a, @a b) == true
   *  @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   *  Can invalidate edge indices -- in other words, old edge(@a i) might not
   *  equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {

  	// assert that a and b are DISTINCT valid nodes of the graph
  	assert( this->has_node(a) && this->has_node(b) && !(a == b) );

    // get node index of both endpoints
    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // check if edge already in graph
    if(this->has_edge(a, b)) {
      
      // construct edge with node1 = a, node2 = b
      Edge existingEdge = Edge(a_ind, b_ind, this);
      return existingEdge;
    }
    // if edge is not in the graph
    else {

      // update adjacency information of node a
      this->nodeAdjacency_[a_ind].push_back(b_ind);

      // update adjacency information of node b
      this->nodeAdjacency_[b_ind].push_back(a_ind);

      // Make new edge
      Edge newEdge = Edge(a_ind, b_ind, this);

      // Store hashed index of this edge
      this->edgeHashes_.push_back(this->node2Edge(a_ind, b_ind));

      return newEdge;
    }
  }

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Use vector clear to remove all nodes and edges from the graph
    this->nodeLocations_.clear();
    this->nodeAdjacency_.clear();
    this->edgeHashes_.clear();
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

    /** Return the Node currently in NodeIterator */
    Node operator*() const {

      // Get the node
      Node currNode = this->graph_->node(this->nodeCount_); 
 	  return currNode;
    }

    /** Update current NodeIterator object to the next node */
    NodeIterator& operator++() {

      // Increment current nodeCount of this NodeIterator
      this->nodeCount_++;

      /* Return the (dereferenced) NodeIterator object, 
       * @a not a pointer to the object */
      return *this;
    }

    /** Test if this NodeIterator and @a ni  are equal
     *  @param @a ni is the NodeIterator to compare to current NodeIterator
     *
     *  Equal NodeIterators have the same graph and same nodeCount */
    bool operator==(const NodeIterator& ni) const {

      // Check if NodeIterator ni and this NodeIterator belong to same graph
      if(this->graph_ != ni.graph_) {
        return false;  // if not, they aren't equal
      }
      // Check if NodeIterator ni and this NodeIterator have the same count
      if(this->nodeCount_ != ni.nodeCount_) {
        return false;  // if not, they aren't equal
      }
      return true;  // otherwise, they are equal
    }

   private:

    friend class Graph;
    
    // Store the iterator value of this NodeIterator, i.e. count of which node
    size_type nodeCount_;

    // Store a pointer to the graph this NodeIterator belongs to
    Graph* graph_;

    // Constructor for NodeIterator (valid) 
    NodeIterator(size_type nodeCount, const Graph* graph)
     : nodeCount_(nodeCount), graph_(const_cast<Graph*>(graph)) {
    }
  };

  /** Return the NodeIterator to start iterating over the nodes */
  NodeIterator node_begin() const {
  	
  	// Start iterating with the first node
  	NodeIterator startingNode = NodeIterator(0, this);
    return startingNode;
  }

  /** Return the NodeIterator to end iterating over the nodes */
  NodeIterator node_end() const {
    
    // Stop iterating after the last node
  	NodeIterator endingNode = NodeIterator(this->num_nodes(), this);
    return endingNode;
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

    /** Return the Edge currently in IncidentIterator */
    Edge operator*() const {

    	// Get the index of the current node adjacent to incidentNode_
    	size_type adjNode = this->graph_->nodeAdjacency_[this->incidentNode_][this->edgeCount_];
    	
    	// Get the edge with these two endpoints, with incidentNode_ = node1 and adjNode = node2 of this edge
    	Edge currEdge = Edge(incidentNode_, adjNode, this->graph_);

    	return currEdge;
    }

    /** Update current IncidentIterator object to the next edge count*/    
    IncidentIterator& operator++() {

      // Increment current edge count of this IncidentIterator
      this->edgeCount_++;

      /* Return the (dereferenced) NodeIterator object, 
       * @a not a pointer to the object */
      return *this;
    }
    
    /** Test if this IncidentIterator and @a ii  are equal
     *  @param @a ii is the IncidentIterator to compare to current IncidentIterator
     *
     *  Equal IncidentIterators have the same graph, incidentNode_ and edgeCount_ */
    bool operator==(const IncidentIterator& ii) const {

      // Check if ii and this IncidentIterator belong to same graph
      if(this->graph_ != ii.graph_) {
      	return false;  // if not, they aren't equal
      }
      // Check if ii and this IncidentIterator have the same incidentNode_ index
      if(this->incidentNode_ != ii.incidentNode_) {
      	return false;  // if not, they aren't equal
      }
      // Check if ii and this IncidentIterator have the same edgeCount_
      if(this->edgeCount_ != ii.edgeCount_) {
      	return false;  // if not, they aren't equal
      }
      return true; // otherwise, they are equal
    }

   private:

    friend class Graph;
    
    // Store the index of the node over whose incident edges we are iterating
    size_type incidentNode_;

    // Store the iterator value of this IncidentIterator, i.e. incident edge counter
    size_type edgeCount_;

    // Store a pointer to the graph this IncidentIterator belongs to
    Graph* graph_;

    // Constructor for IncidentIterator (valid) 
    IncidentIterator(size_type incidentNode, size_type edgeCount, const Graph* graph)
     : incidentNode_(incidentNode),edgeCount_(edgeCount), graph_(const_cast<Graph*>(graph)) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    /** Return the Edge currently in EdgeIterator */
    Edge operator*() const {

      // Get the edge
      Edge currEdge = this->graph_->edge(this->edgeCount_); 
 	  return currEdge;
    }

    /** Update current EdgeIterator object to the next edge */
    EdgeIterator& operator++() {

      // Increment current edgeCount of this EdgeIterator
      this->edgeCount_++;

      /* Return the (dereferenced) EdgeIterator object, 
       * @a not a pointer to the object */
      return *this;
    }
    
    /** Test if this EdgeIterator and @a ei  are equal
     *  @param @a ei is the EdgeIterator to compare to current EdgeIterator
     *
     *  Equal EdgeIterators have the same graph and same edgeCount */
    bool operator==(const EdgeIterator& ei) const {

      // Check if EdgeIterator ei and this EdgeIterator belong to same graph
      if(this->graph_ != ei.graph_) {
        return false;  // if not, they aren't equal
      }
      // Check if EdgeIterator ei and this EdgeIterator have the same count
      if(this->edgeCount_ != ei.edgeCount_) {
        return false;  // if not, they aren't equal
      }
      return true;  // otherwise, they are equal
    }

   private:
    friend class Graph;

    /** Store the iterator value of this EdgeIterator, i.e. what edge count 
     *  use as (index into edge hash table to get edge information) */
    size_type edgeCount_;

    // Store a pointer to the graph this EdgeIterator belongs to
    Graph* graph_;

    // Constructor for EdgeIterator (valid) 
    EdgeIterator(size_type edgeCount, const Graph* graph)
     : edgeCount_(edgeCount), graph_(const_cast<Graph*>(graph)) {
    }


  };

  /** Return the EdgeIterator to start iterating over the edges */
  EdgeIterator edge_begin() const {
  	
  	// Start iterating with the first edge
  	EdgeIterator startingEdge = EdgeIterator(0, this);
  	return startingEdge;
  }
  
  /** Return the EdgeIterator to stop iterating over the edges */
  EdgeIterator edge_end() const {
    
    // Stop iterating after the last edge
  	EdgeIterator endingEdge = EdgeIterator(this->num_edges(), this);
    return endingEdge;
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Vector to store Node Point locations, with the index in the vector 
   *  equal to the Node's index **/
  std::vector<NodeInfo> nodeLocations_;

  /** Vector to store the adjacent nodes of each node, ordered by node index. 
   *
   *  Note that if (a,b) is an edge, then b is listed in the adjacency list
   *  of a AND a is listed in the adjacency list of b */
  std::vector< std::vector<size_type> > nodeAdjacency_;

  /** Store the edge index hashes for graph (with this index from node2Edge) */ 
  std::vector<size_type_long> edgeHashes_;

  /** Struct for storing the various members of our Nodes */
  struct NodeInfo {

 	// location of node in space
  	Point pos_;

  	// value of the node
  	node_value_type val_;
  };

  /** Function to turn pair of edge endpoints into a @a unique edge index @a e
   *  This is the standard Cantor Pairing Function (see https://en.wikipedia.org/wiki/Pairing_function)
   *  
   *  @pre (@a a, @a b) are valid indices of nodes of an edge in the graph
   *
   *  @return @a e is the index of the edge connecting nodes @a a and @a b. This is an unsigned long to 
   *		  prevent overflow, as @a e could be larger than an unsigned int.
   *
   *  @post node2Edge(@a a, @a b) = node2Edge(@a b, @a a), as in HW 1  we only have undirected edges. 
   *        If this is not the desired behavior, the marked line below can be commented out, as in
   *		general the pairing function will NOT return node2Edge(a,b) = node2Edge(b,a). 
   *
   *  Complexity: O(1)
   */
  size_type_long node2Edge(size_type a, size_type b) {
  	
  	// Check Node objects with indices a and b
  	Node aNode = node(a);
  	Node bNode = node(b);

  	// Check that we have a valid edge
  	assert(this->has_edge(aNode, bNode));

  	/** THESE LINES CAN COMMENTED OUT TO ACCOMODATE DIRECTED EDGES, i.e. so that
  	*   node2Edge(a, b) != node2Edge(b, a) */
  	// To enforce node2Edge(a, b) = node2Edge(b, a), let k1 = min(a,b) and k2 = max(a,b)
  	size_type k1 = fmin(a,b);
  	size_type k2 = fmax(a,b);

  	// Cantor pairing function
  	size_type_long e = 0.5*(k1 + k2)*(k1 + k2 + 1) + k2;
    return e;
  }

  /** Struct for storing the Node indices of the endpoints of an Edge */
  struct edgeEndpts{
  	size_type a;
  	size_type b;
  };

  /** Function to turn an edge index @e into the edge endpoints (@a a, @a b).
   *  This is the standard Cantor Pairing Function (see https://en.wikipedia.org/wiki/Pairing_function)
   *  
   *  @pre @a e is a valid index of an edge in the graph
   * 
   *  @return @ s is a struct containing {@ a a and @a b}, which are the node indices of the 
   *          endpoint nodes of edge @a e
   *
   *  @post edge2Node(@a e) returns (@a a, @a b) with @a a < @a b, as in HW 1 we can
   *        only have undirected edges, so all edge indices are created with this convention.
   *        Just keep in mind that we prescribed node2Edge(@a a, @a b) = node2Edge(@a b, @a a) 
   *
   *  Complexity: O(1)
   */
  edgeEndpts edge2Node(size_type_long e) const {

  	// Cantor pairing function inversion intermediary variables
  	size_type_long x = floor( 0.5*(sqrt(8*e + 1) - 1) );
  	size_type_long y = 0.5*(x*x + x);

  	// Cantor pairing function inversion return values
  	size_type b = e - y;
  	size_type a = x - b;

  	// Wrap endpoint indices in a struct to return
  	edgeEndpts ee = {a, b};
  	return ee;
  }

};

#endif // CME212_GRAPH_HPP