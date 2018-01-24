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
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. Edges are unique (there is at
 *  most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  /** Use this space for declarations of important internal types you need
      later in the Graph's definition. */

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // Type of this graph.
  using graph_type = Graph;

  // Predeclaration of Node type.
  class Node;
  /// Synonym for Node (following STL conventions).
  using node_type = Node;

  /// Predeclaration of Edge type.
  class Edge;
  /// Synonym for Edge (following STL conventions).
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  // Construct an empty graph.
  Graph() : nodeLocations_(), nodeAdjacency_(), edgeEndpoints_() {
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
  class Node {

   public:

    // Construct an invalid/uninitialized node.
    Node() {
    }

    // Return this node's position.
    const Point& position() const {
      /** position is stored in the Graphs's nodeLocations_ array indexed by 
        * the Node's index. Here, ``this'' refers to the Node. */
      return this->graph_->nodeLocations_[this->index_];
    }

    // Return this node's index, a number in the range [0, graph_size).
    size_type index() const {
      return this->index_;
    }

    /** Test whether this node and @a n are equal.
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
     *  The node ordering relation must obey trichotomy: For any two nodes x
     *  and y, exactly one of x == y, x < y, and y < x is true. */
    bool operator<(const Node& n) const {
      if(this->index_ < n.index_) {
        return true;
      }
      return false;
    }

   private:

    // Allow Graph to access Node's private member data and functions
    friend class Graph;
  
    // Store the index of this node
    size_type index_;

    // Store a pointer to the graph this node belongs to
    const Graph* graph_;

    // Constructor for Node (valid) 
    Node(size_type index, const Graph* graph)
     : index_(index), graph_(const_cast<Graph*>(graph)) {

     }
  };

  /** Return the number of nodes in the graph.
   *  Complexity: O(1). */
  size_type size() const {
    /** nodeLocations_ holds the positions of all nodes in graph,
        so size() returns the number of nodes */
    return this->nodeLocations_.size();
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
    /// Construct a new node
    Node newNode = Node(this->num_nodes(), this);

    // add new node position to nodeLocations_ at index node_index 
    this->nodeLocations_.push_back(position);
    
    /** initialize nodeAdjacency[this->num_nodes] with empty vectors for later use
        storing adjacencies of each node */
    this->nodeAdjacency_.push_back(std::vector<std::vector<size_type>>());

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
  	assert( (i >=0) && (i < num_nodes()) );
    
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
  class Edge{

   public:
    
    //Construct an invalid/uninitialized Edge.
    Edge() {
    }

    // Return a node of this Edge
    Node node1() const {

      // retrieve the first node's index from the Graph's edgeEndpoint vector of Edge endpoint tuples
      size_type node_ind = this->graph_->edgeEndpoints_[this->index_][0];

      // Get node with index node_ind
      Node newNode = this->graph_->node(node_ind);

      return newNode; 
    }

    // Return the other node of this Edge
    Node node2() const {
      
      // retrieve the second node's index from the Graph's edgeEndpoint vector of Edge endpoint tuples
      size_type node_ind = this->graph_->edgeEndpoints_[this->index_][1];

      // Get node with index node_ind
      Node newNode = this->graph_->node(node_ind);

      return newNode;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      // Get endpoint nodes of current edge
      size_type thisNode0 = this->graph_->edgeEndpoints_[this->index_][0];
      size_type thisNode1 = this->graph_->edgeEndpoints_[this->index_][1];

      // Get endpoint nodes of other edge e
      size_type otherNode0 = e.graph_->edgeEndpoints_[e.index_][0];
      size_type otherNode1 = e.graph_->edgeEndpoints_[e.index_][1];

      /** if the first endpoint of current edge is not equal to either endpoint of
          other edge, then these can't be the same edge */
      if( (thisNode0 != otherNode0) && (thisNode0 != otherNode1)) {
        return false;
      } 
  
      /** if the first endpoint of current edge is equal to either endpoint of
          other edge... */
      else if(thisNode0 == otherNode0) {
        // ...then the remaining endpoints of each edge need to be equal as well
        if(thisNode1 != otherNode1){

          return false;
        }
      } 

      /** if the first endpoint of current edge is equal to either endpoint of
          other edge... */
      else {
        // ...then the remaining endpoints of each edge need to be equal as well
        if(thisNode1 != otherNode0){
          
          return false;
        }
      }

      return true;
    }

    /** Test whether this edge is less than @a e in a global order. */
    bool operator<(const Edge& e) const {
      if(this->index_ < e.index_) {
        return true;
      }
      return false;
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    

    // Store the index of this edge
    size_type index_;

    // Store a pointer to the graph this edge belongs to
    const Graph* graph_;

    // Constructor for Edge (valid)
    Edge(size_type index, const Graph* graph)
     : index_(index), graph_(const_cast<Graph*>(graph)) {
     }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    /** edgeEndpoints holds the indices of endpoint nodes for all
        edges in graph, so size() returns the number of edges */
    return this->edgeEndpoints_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  	// assert edge index is valid
  	assert( (i >= 0) && (i < num_edges()) );

    /** Create new edge with index i, return this. This can be used along with
        the Graph's edgeEndpoints_ array to determine the endpoint nodes */
    Edge newEdge = Edge(i, this);
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
  	assert( has_node(a) && has_node(b) );

    // get node index of both endpoint nodes
    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // Search through nodeAdjacency_ vector of tuples to see if b and a are neightbors
    for(unsigned i=0; i < this->nodeAdjacency_[a_ind].size(); i++) {
      size_type adj_node = this->nodeAdjacency_[a_ind][i][0];

      // if the adjcacent node index is the same as b's index, then there is an edge
      if(adj_node == b_ind){
        return true;
      }
    }
    return false;  // no edge
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
  	// assert that a and b are DISTINCT valid nodes of the graph
  	assert( has_node(a) && has_node(b) && !(a == b) );

    // get node index of both endpoints
    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // check if edge already in graph
    if(has_edge(a, b)) {
      
      size_type edge_ind;
      // Find edge index for (a,b) endpoints
      for(unsigned i=0;i < this->nodeAdjacency_[a_ind].size();i++) {
        size_type adj_node = this->nodeAdjacency_[a_ind][i][0];
        
        // save edge index once its found
        if(adj_node == b_ind){
          edge_ind = this->nodeAdjacency_[a_ind][i][1];      
        }
      }
      
      // construct node with index edge_ind
      Edge existingEdge = edge(edge_ind);
      return existingEdge;
    }
    // if edge is not in the graph
    else {
      // edge_ind is the next edge to be added
      size_type edge_ind = this->num_edges();

      // update adjacency information of node a
      std::vector<unsigned> temp_b = {b_ind, edge_ind};
      this->nodeAdjacency_[a_ind].push_back(temp_b);

      // update adjacency information of node b
      std::vector<unsigned> temp_a = {a_ind, edge_ind};
      this->nodeAdjacency_[b_ind].push_back(temp_a);

      // Make new edge
      Edge newEdge = Edge(edge_ind, this);

      // add edge endpoints for new edge
      std::vector<unsigned> temp_ab = {a_ind, b_ind};
      this->edgeEndpoints_.push_back(temp_ab);

      return newEdge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Use vector clear to remove all nodes and edges from the graph
    this->nodeLocations_.clear();
    this->nodeAdjacency_.clear();
    this->edgeEndpoints_.clear();
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Vector to store Node Point locations, with the index in the vector 
      equal to the Node's index **/
  std::vector<Point> nodeLocations_;

  /** Vector to store the adjacent nodes of each node, along with the edge
      that connects the two nodes **/
  std::vector<std::vector<std::vector<unsigned>>> nodeAdjacency_;

  /** Vector to store the indices of the two endpoints of each edge **/
  std::vector<std::vector<unsigned>> edgeEndpoints_;

};

#endif // CME212_GRAPH_HPP
