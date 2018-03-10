#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. Edges are unique (there is at
 *  most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:

  /** Use this space for declarations of important internal types you need
      later in the Graph's definition. */
  
  // Predeclaration of struct to hold NodeInfo
  struct NodeInfo;
  
  public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // Type of this graph.
  using graph_type = Graph;

  // Type of indexes and sizes.
  using size_type = unsigned;

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
  /// Synonym for Edge Value Type, specified since Graph is a template
  using edge_value_type = E;

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
  Graph() : graph_nodes_(), adj_(), edge_values_(), active_nodes_() {
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
   *
   *  Note: the following methods are only defined for active nodes (i.e. not removed) from the graph
   */
  class Node : private totally_ordered<Node>{

   public:

    // Construct an invalid/uninitialized node.
    Node() {
    }

    /** Return this node's position, as a reference (which could be changed)
     *  position is a member of the struct stored in the Graphs's 
     *  graph_nodes_ vector indexed by the Node's identifier node_id_. */
	  Point& position() {
      return graph_->graph_nodes_[node_id_].pos_;
    }

    // Return this node's position, as a reference (which CANNOT be changed).
    const Point& position() const {
      return position();
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     *  This non-static value is the position of this node_id in the Graph's active_nodes_ vector. */
    size_type index() const {
      // find the node index from the node identifier
      auto it = std::find(graph_->active_nodes_.begin(), graph_->active_nodes_.end(), node_id_);

      return std::distance(graph_->active_nodes_.begin(), it);
    }

    /** Return this node's value, as a reference (which could be changed).
     *  value is a member of the struct stored in the Graphs's 
     *  graph_nodes_ vector indexed by the Node's identifier node_id_. */
	  node_value_type& value() {
	    return graph_->graph_nodes_[node_id_].val_;
	  }
    
    // Return this node's value, as a reference (which CANNOT be changed).
    const node_value_type& value() const {
      return value();
    }   

    /** Return the degree of this node, i.e. how many incident edges it has
     *  from the size of its vector of adjacent nodes. */
    size_type degree() const {
      return graph_->adj_[node_id_].size();
    }

    /** Return the IncidentIterator to start iterating over the incident edges */
    IncidentIterator edge_begin() const {
      // Start iterating with the first edge incident to this node
  	  return IncidentIterator(graph_, node_id_);
    }

    /** Return the IncidentIterator to end iterating over the incident edges */
    IncidentIterator edge_end() const{
      // Stop iterating with the last edge incident to this node
  	  return IncidentIterator(graph_, node_id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *  @param[in] @a n Node to test.
     *
     *  Equal nodes have the same graph and the same identifier. */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) && (node_id_ == n.node_id_) );
    }

    /** Test whether this node is less than @a n in a global order.
     *  @param[in] @a n Node to test.
     *
     *  The node ordering relation must obey trichotomy: For any two nodes x
     *  and y, exactly one of x == y, x < y, and y < x is true. 
     *
     *  Global ordering of graphs as follows: for graphs a and b created
     *  in that order, global ordering as follows, with (g,n) tuple denoting
     *  (graph letter, node identifier):
     *  (a,0), (a,1), ..., (a,size(a)-1), (b,0), ..., (b,size(b)-1) */
    bool operator<(const Node& n) const {

      if(graph_ != n.graph_) {
      	/** for nodes in different graphs, use std::operator< to 
      	 *  return a comparison of the graph pointer values */
        return graph_ < n.graph_;
      } else {
      	// for nodes in the same graph, compare indices
	      return node_id_ < n.node_id_;
      }
    }

   private:

    // Allow Graph to access Node's private member data and functions
    friend class Graph;

    // Store a pointer to the graph this node belongs to
    Graph* graph_;

    // Store the index of this node
    size_type node_id_;

    /** Constructor for Node (valid) 
     *  @param[in] index  Index of node to construct
     *  @param[in] graph  Graph pointer for graph of node to construct */
    Node(const Graph* graph, size_type node_id)
     : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {
     }
  };

  /** Return the number of active nodes in the graph.
   *  Complexity: O(1). */
  size_type size() const {
    return active_nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // get new unique node identifier
    size_type node_id = graph_nodes_.size();

    // mark this node as active
    active_nodes_.push_back(node_id);

    // Construct a new NodeInfo struct to hold position, value and unique node identifier
    NodeInfo newNodeInfo = NodeInfo(node_id, position, value);

    // add new NodeInfo struct to graph_nodes_ at position node.index() 
    graph_nodes_.push_back(newNodeInfo);

    // add this node to the adj_ vector (to later store adj node ids)
    adj_.push_back(std::vector<size_type>());


    return Node(this, node_id);
  }

  
  /** Remove a node (and all incident edges) from the graph, returning the number of edges deleted.
   * @param[in] n   The node to remove
   * @return the number of edge deleted (should equal n.degree()) 
   *
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edgees() - degree(n)
   * @post for nodes with old index j > index i = n.index(), new index j = old index j -1
   *
   * Complexity: O((num_nodes()) operations to simply remove the node
   *       Removing the adjacent edges is O(n.degree()*ln(num_edges())) ~ O(ln(num_edges())) 
   *       since n.degree() is assumed to be much smaller than num_edges()
   */
  size_type remove_node(const Node& n) {
    // keep count of removed edges
    size_type edges_deleted_count = 0; 

    // to remove a node, first remove all adjacent edges (in reverse order)
    for (auto ii = n.edge_end(); ii != n.edge_begin(); ) {
      // decrement in the loop to ensure that we remove n.edge_begin()
      --ii;

      // remove the edge, update counter if successful
      Edge e = *ii;
      size_type ret = remove_edge(e);
      if (ret)
        ++edges_deleted_count;
    }

    // then remove node by marking it inactive (delete from active_nodes)
    // find position of node in active_nodes
    auto n_it = std::find(active_nodes_.begin(), active_nodes_.end(), n.node_id_);
    // remove node from the list
    active_nodes_.erase(n_it);

    return edges_deleted_count;
  }


  /** Remove a node (and all incident edges) from the graph, returning the number of edges deleted.
   * @param[in] n_it   Iterator to the node to remove
   * @return iterator to the new node at original index of n_it, 
   *         unless original index = num_nodes() -1, then return the iterator to the NEW ending index  
   *
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edgees() - degree(n)
   * @post for nodes with old index j > index i = n.index(), new index j = old index j -1
   * @post all iterators to values after the removed node will be invalidated
   *
   * Complexity: O((num_nodes()) operations to simply remove the node
   *       Removing the adjacent edges is O(n.degree()*ln(num_edges())) ~ O(ln(num_edges())) 
   *       since n.degree() is assumed to be much smaller than num_edges()
   */
  node_iterator remove_node(node_iterator n_it) {
    // turn this node iterator into a node
    Node n = (*n_it);

    /** return this iterator (to same node index, which will be a new node) if n_it doesn't point
     *  to the last node, otherwise return iterator to new last index*/
    if (n_it.nodeit_idx_ == num_nodes() -1)
      --n_it;

    // remove this node
    remove_node(n);

    return n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_== this);
  }

  /** Return the node with index @a i.
   *  @pre 0 <= @a i < num_nodes()
   *  @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Return invalid node if invalid index provided
  	if ( (i < 0) || (i >= num_nodes()) )
      return Node();
    //assert((i>=0) && (i < num_nodes()));

    // Otherwise, create new node with the proper node identifier (found in active_nodes_ vector) and return this.     
    size_type node_id = active_nodes_[i];
    return Node(this, node_id);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   *  Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   *  are considered equal if they connect the same nodes, in either order.
   *
   *  Note: the following methods should only be called for valid edges.
   */
  class Edge : private totally_ordered<Edge>{

   public:
    
    // Construct an invalid/uninitialized Edge.
    Edge() {
    }

    // Return a node of this Edge
    Node node1() const {
      // Return node with identifier edge_node1_id_
      return (Node(graph_, edge_node1_id_));
    }

    // Return the other node of this Edge
    Node node2() const {
      // Return node with identifier edge_node2_id_
      return (Node(graph_, edge_node2_id_));
    }

    /** Return this edge's value, as a reference (which could be changed).
     *  @return  value of this edge, which is a value of the map stored in the Graphs's edge_values_, 
     *     with key determined by the endpoints of the edges (ordered (min, max) node id).
     *
     *  Complexity: O(num_edges)  
     */
    edge_value_type& value() {
      // search for this edge in edge_values_ by (endpt1, endpt2) key..return iterator to this edge
      size_type min_id = std::min(edge_node1_id_, edge_node2_id_);
      size_type max_id = std::max(edge_node1_id_, edge_node2_id_);

      auto e_it = graph_->edge_values_.find(std::make_pair(min_id, max_id));

      // error if we can't find the edge - should not be attemping to get value of an invalid edge
      assert(e_it != graph_->edge_values_.end());

      // otherwise, return found edge
      return (e_it->second);
    }
    
    /** Return this edge's value, as a reference (which CANNOT be changed)
     *  @return  value of this edge, which is a value of the map stored in the Graphs's edge_values_, 
     *     with key determined by the endpoints of the edges (ordered (min, max) node id).
     *
     *  Complexity: O(num_edges)  
     */
    const edge_value_type& value() const {
      return value();
    }

    // Return the length of this edge (the euclidean distance between endpoints)
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     *  Equal edges represent the same edge between two nodes, and will have the same
     *  edge endpoint nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if edge e and this edge belong to same graph
      return ( (graph_ == e.graph_) && 
               ( ((e.edge_node1_id_ == edge_node1_id_)  && (e.edge_node2_id_ == edge_node2_id_)) 
               || ((e.edge_node1_id_ == edge_node2_id_)  && (e.edge_node2_id_ == edge_node1_id_)) ));
    }
    
    /** Test whether this edge is less than @a e in a global order.
	   *
     *  Global ordering  doens't have an interpretive meaning.
     *  It is defined as follows: for graphs a and b created
     *  in that order, global ordering as follows, with (g,i) tuple denoting
     *  (graph pointer value, node1 id):
     *  (a,min_id), ..., (a,max_id), (b,min_id), ..., (b,max_id) 
     */
    bool operator<(const Edge& e) const {
      // Check if edge e and this edge belong to same graph
      if(graph_ != e.graph_) {
      	/** for edges in different graphs, use std::operator< to 
      	 *  return a comparison of the graph pointer values */
        return graph_ < e.graph_;
      } 
      else {
	      // For edges in the same graph, first compare node1 ids
        // if they're the same, return comparison of node2 ids
	      if (edge_node1_id_ == e.edge_node1_id_)
          return edge_node2_id_ < e.edge_node2_id_;
        //otherwise, return comparision of node1 ids
        else
          return edge_node1_id_ < e.edge_node1_id_;
      }
    }
    
   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Store a pointer to the graph this edge belongs to
    Graph* graph_;

    // Store the node identifier of node1 of this edge
    size_type edge_node1_id_;

    // Store the node identifier of node2 of this edge
    size_type edge_node2_id_;

    /** Constructor for Edge (valid)
     *  @param[in] graph  Pointer to Graph for this edge 
     *  @param[in] edge_node1_id  Node identifier of first endpoint of edge
     *  @param[in] edge_node2_id  Node identifier of second endpoint of edge
     */
    Edge(const Graph* graph, size_type edge_node1_id, size_type edge_node2_id)
     : graph_(const_cast<Graph*>(graph)), edge_node1_id_(edge_node1_id), edge_node2_id_(edge_node2_id) {
     }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_values_.size();
  }

  /** Return the edge with index @a i.
   *  @pre 0 <= @a i < num_edges()
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  	// Return invalid edge if invalid index is provided
  	if ( (i < 0) || (i >= num_edges()) )
      return Edge();
    //assert((i >= 0) && (i < num_edges()));

  	// Return iterator to ith edge
    auto e_it = edge_begin() + i;
    return (*e_it);
  }

  /** Test whether two nodes are connected by an edge.
   *  @pre @a a and @a b are valid nodes of this graph
   *  @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
  	// return false if a and b are not valid nodes of the graph
  	if( !(has_node(a)) || !(has_node(b)) )
      return false;

    // Search through node a's adjacent nodes to see if b and a are neighbors
    for(auto ii = a.edge_begin(); ii != a.edge_end(); ++ii) {
      Node adj_node = (*ii).node2();

      // if the adjacent node is equal to node b, then there is an edge
      if(adj_node == b)
        return true;
    }
    return false;  // no edge
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   *  @pre @a a and @a b are distinct valid nodes of this graph
   *  @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *  @post has_edge(@a a, @a b) == true
   *  @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *        Else,                        new num_edges() == old num_edges() + 1.
   *
   *  Can invalidate edge indices -- in other words, old edge(@a i) might not
   *  equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // get node identifiers of both endpoints
    size_type a_id = active_nodes_[a.index()];
    size_type b_id = active_nodes_[b.index()];

    // check if edge already in graph by iterating over incident edges of node a
    for (auto ii = a.edge_begin(); ii != a.edge_end(); ++ii) {
      if (b == (*ii).node2()) {
        // return edge with node1 = a, node2 = b
        return Edge(this, a_id, b_id);
      }
    }

    // If we don't already have the edge, make a new one!
    // update adjacency information of node a
    adj_[a_id].push_back(b_id);

    // update adjacency information of node b
    adj_[b_id].push_back(a_id);

    // update edge values
    size_type min_id = std::min(a_id, b_id);
    size_type max_id = std::max(a_id, b_id);
    
    edge_values_.insert(std::make_pair(std::make_pair(min_id, max_id), value));

    // return new edge
    return Edge(this, a_id, b_id);
  }

  /** Remove an edge from the graph, returning the number of edges deleted.
   *  @param[in] n1   Endpoint 1 of the edge to remove
   *  @param[in] n2   Endpoint 2 of the edge to remove
   *  @return the number of edges deleted (should equal 1 unless theres an error) 
   *          If !has_edge(n1, n2), then return 0     
   *
   *  @pre @a n1 and @a n2 must be valid nodes in the graph.
   *  @post new num_edges() == old num_edgees() - 1
   *
   *  Complexity: O(ln(num_edges()) operations
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // keep count of removed edges
    size_type edge_deleted_count = 0; 

    // remove edge from edge_values_
    size_type min_id = std::min(n1.node_id_, n2.node_id_);
    size_type max_id = std::max(n1.node_id_, n2.node_id_);
    
    size_type original_size = edge_values_.size();
    edge_values_.erase(std::make_pair(min_id, max_id));
    
    // return 0 if the remove is unsuccessful
    if (edge_values_.size() == original_size)
      return edge_deleted_count;

    // update number of edge deleted if the remove is successful
    edge_deleted_count++;

    // then remove n2 from adj_[n1]...
    auto n1_it = std::find(adj_[n1.node_id_].begin(), adj_[n1.node_id_].end(), n2.node_id_);
    // remove node from the list
    adj_[n1.node_id_].erase(n1_it);

    // ...and remove n1 from adj_[n2]
    auto n2_it = std::find(adj_[n2.node_id_].begin(), adj_[n2.node_id_].end(), n1.node_id_);
    // remove node from the list
    adj_[n2.node_id_].erase(n2_it);

    return edge_deleted_count;
  }

  /** Remove an edge from the graph, returning the number of edges deleted.
   *  @param[in] e   Edge to remove
   *  @return the number of edges deleted (should equal 1 unless theres an error) 
   *          If !has_edge(n1, n2), then return 0     
   *
   *  @pre @a e must be a valid edge in the graph.
   *  @post new num_edges() == old num_edgees() - 1
   *
   *  Complexity: O(ln(num_edges()) operations
   */
  size_type remove_edge(Edge& e) {
    Node node1 = Node(this, e.edge_node1_id_);
    Node node2 = Node(this, e.edge_node2_id_);
    return remove_edge(node1, node2);
  }

  /** Remove an edge from the graph, returning the number of edges deleted.
   *  @param[in] e_it   Iterator for edge to remove
   *  @return the number of edges deleted (should equal 1 unless theres an error) 
   *          If !has_edge(n1, n2), then return 0     
   *
   *  @pre @a e_it must be an iterator to a valid edge in the graph.
   *  @post new num_edges() == old num_edgees() - 1
   *  @post invalidates all iterators to edges after the removed edge
   *
   *  Complexity: O(ln(num_edges()) operations
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // dereference edge_iterator
    Edge e = (*e_it); 

     /** return this iterator (to same edge index) if e_it doesn't point to the last edge, 
     *  otherwise return iterator to new last index*/
    if ((e_it.curr_root_idx_ == num_nodes() -1) && (e_it.leaf_count_ == node(e_it.curr_root_idx).degree()-1))
      --e_it;

    // remove this edge
    remove_edge(e);

    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Use vector clear to remove all nodes and edges from the graph
    adj_.clear();
    edge_values_.clear();
    active_nodes_.clear();
    graph_nodes_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   *  @brief Iterator class for nodes. A forward iterator. 
   * 
   *  Iterates over all ACTIVE (i.e. not deleted) nodes of the graph
   */
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

    /** Return the Node currently in NodeIterator, i.e this node with index nodeit_idx_ */
    Node operator*() const {
      return graph_->node(nodeit_idx_); 
    }

    /** Update current NodeIterator object to the next node */
    NodeIterator& operator++() {
      // Increment current node index of this NodeIterator
      nodeit_idx_++;

      // Return the (dereferenced) NodeIterator object
      return *this;
    }

    /** Test if this NodeIterator and @a ni  are equal
     *  @param @a ni is the NodeIterator to compare to current NodeIterator
     *
     *  Equal NodeIterators have the same graph and same node index 
     */
    bool operator==(const NodeIterator& ni) const {
      return ( (graph_ == ni.graph_) && (nodeit_idx_ == ni.nodeit_idx_) );
    }

    /** Return a NodeIterator that is offset past current iterator
     *  @param @a offset  is the amount to increment the nodeit_idx_ of the current iterator
     *
     *  @pre offset must be a valid offset (i.e. NodeIterator.nodeit_idx_ + offset < NodeIterator.node_end())
     *  @post initial NodeIterator is unchanged: a copy with an updated nodeit_idx_ is returned
     */
    NodeIterator& operator+(const size_type offset) const {
      auto ni_out = this;
      for (size_type i = 0; i <= offset; ++i)
        ++(*ni_out);
      return *ni_out;;
    }

    /** Update current NodeIterator object to the previous edge */    
    NodeIterator& operator--() {
      /// Increment current node index of this NodeIterator
      --nodeit_idx_;

      // Return the (dereferenced) NodeIterator object, 
      return *this;
    }

   private:

    friend class Graph;
    
    // Store a pointer to the graph this NodeIterator belongs to
    Graph* graph_;

    // Store the iterator value of this NodeIterator, i.e. count of which node
    size_type nodeit_idx_;

    /** Constructor for NodeIterator (valid) 
     *  @param[in] nodeit_idx  node index the iterator starts with (default to 0)
     *  @param[in] graph       Pointer to the graph of this iterator 
     */
    NodeIterator(const Graph* graph, const size_type nodeit_idx = 0)
     : graph_(const_cast<Graph*>(graph)), nodeit_idx_(nodeit_idx) {
    }
  };

  /** Return the NodeIterator to start iterating over the nodes */
  NodeIterator node_begin() const {
  	// Start iterating with the first (zeroeth) node
  	return NodeIterator(this);
  }

  /** Return the NodeIterator to end iterating over the nodes */
  NodeIterator node_end() const {    
    // Stop iterating after the last node
  	return NodeIterator(this, num_nodes());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   *  @brief Iterator class for edges incident to a node. A forward iterator. 
   *
   *  Note: should only be called for active root nodes. Iterates over all edges incident to this root node.
   */
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
    	// Get the node identifier of the current adjacent node (to node with node id root_id_)
    	size_type adj_id = graph_->adj_[root_id_][leaf_count_];

    	// Return the edge with these two endpoints, with root node = node1 and adj node = node2 of this edge
    	return Edge(graph_, root_id_, adj_id);
    }

    /** Update current IncidentIterator object to the next edge */    
    IncidentIterator& operator++() {
      // Increment current leaf count of this IncidentIterator
      ++leaf_count_;

      // Return the (dereferenced) NodeIterator object
      return *this;
    }
    
    /** Test if this IncidentIterator and @a ii  are equal
     *  @param @a ii is the IncidentIterator to compare to current IncidentIterator
     *
     *  Equal IncidentIterators have the same graph, root_id_ and leaf_count_ 
     */
    bool operator==(const IncidentIterator& ii) const {
      return ( (graph_ == ii.graph_) &&
               (root_id_ == ii.root_id_) &&
               (leaf_count_ == ii.leaf_count_) );
    }

    /** Return an IncidentIterator that is offset past current iterator
     *  @param @a offset  is the amount to increment the leaf_count_ of the current iterator
     *
     *  @pre offset must be a valid offset (i.e. IncidentIterator.leaf_count_ + offset < IncidentIterator.edge_end())
     *  @post initial IncidentIterator is unchanged: a copy with an updated leaf_count_ is returned
     */
    IncidentIterator& operator+(const size_type offset) {
      auto ii_out = this;
      for (size_type i = 0; i <= offset; ++i)
        ++(*ii_out);
      return *ii_out;
    }

    /** Update current IncidentIterator object to the previous edge */    
    IncidentIterator& operator--() {
      // Decrement current leaf count of this IncidentIterator
      --leaf_count_;

      // Return the (dereferenced) NodeIterator object, @a not a pointer to the object
      return *this;
    }

   private:

    friend class Graph;
    
    // Store a pointer to the graph this IncidentIterator belongs to
    Graph* graph_;

    // Store the node identifier of the node over whose incident edges we are iterating
    size_type root_id_;

    // Store the iterator value of this IncidentIterator, i.e. incident edge counter for node root
    size_type leaf_count_;

    /** Constructor for IncidentIterator (valid)
     *  @param[in] root_id  index of root node 
     *  @param[in] leaf_count     position of the iterator (default = 0)
     *  @param[in] graph         pointer to graph for incident iterator
     */
    IncidentIterator(const Graph* graph, const size_type root_id, const size_type leaf_count = 0)
     : graph_(const_cast<Graph*>(graph)), root_id_(root_id),leaf_count_(leaf_count) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   *  @brief Iterator class for edges. A forward iterator. 
   */
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
      // get current root node identifier
      size_type curr_root_id = graph_->active_nodes_[curr_root_idx_];
      // get current leaf node identifier
      size_type leaf_id = graph_->adj_[curr_root_id][leaf_count_];
      
      // return this edge
      return Edge(graph_, curr_root_id, leaf_id);
    }

    /** Update current EdgeIterator object to the next edge 
     *  @return next EdgeIterator object
     *    If we can't iterate further, return EdgeIterator object with:
     *          curr_root_id_ = num_nodes and leaf_count_ = 0 <=> this.edge_end()
     */
    EdgeIterator& operator++() {
      // increment to next leaf node
      leaf_count_++;

      // iterate over at most all possible root indices (i.e. all active nodes in graph)
      for( ; curr_root_idx_ != graph_->num_nodes(); ++curr_root_idx_) {
        // get current root node identifier
        size_type curr_root_id = graph_->active_nodes_[curr_root_idx_];

        // iterate over all incident edges of current root node
        size_type root_degree = graph_->node(curr_root_idx_).degree();
        for ( ; leaf_count_ != root_degree; ++leaf_count_) {
          // get current leaf node identifier
          size_type leaf_id = graph_->adj_[curr_root_id][leaf_count_];
          
          // don't double count edges
          if (curr_root_id < leaf_id)
            // Return the (dereferenced) EdgeIterator object
            return *this;
        }

        // iterated through all incident edges of curr root node: move to next root node and reset leaf count
        leaf_count_ = 0;
      }

      // Return the (dereferenced) EdgeIterator object (will be this.edge_end() if we get there)
      return *this;
    }

    /** Update current EdgeIterator object to the previous edge 
     *  @return previous EdgeIterator object
     *    If we can't iterate further, return EdgeIterator object with:
     *          curr_root_id_ = 0 and leaf_count_ = 0 <=> this.edge_begin()
     */  
    EdgeIterator& operator--() {
      // decrement leaf count
      --leaf_count_;

      // iterate over at most all possible root indices (i.e. all active nodes in graph)
      for( ; curr_root_idx_ >= 0; --curr_root_idx_) {
        // get current root node identifier
        size_type curr_root_id = graph_->active_nodes_[curr_root_idx_];

        // iterate over all incident edges of current root node
        for ( ; leaf_count_ >= 0; --leaf_count_) {
          // get current leaf node identifier
          size_type leaf_id = graph_->adj_[curr_root_id][leaf_count_];
          
          // don't double count edges
          if (curr_root_id < leaf_id)
            // Return the (dereferenced) EdgeIterator object
            return *this;
        }
   
        // iterated through all incident edges of curr root node: move to next root node and reset leaf count
        leaf_count_ = graph_->node(curr_root_idx_).degree() - 1;
      }

      // Return the (dereferenced) EdgeIterator object (will be this.edge_begin() if we get here)
      leaf_count_ = 0;
      curr_root_idx_ = 0;
      return *this;
    }
    
    /** Test if this EdgeIterator and @a ei are equal
     *  @param @a ei is the EdgeIterator to compare to current EdgeIterator
     *
     *  Equal EdgeIterators have the same graph, curr_root_idx and leaf_count 
     */
    bool operator==(const EdgeIterator& ei) const {
      return ( (graph_ == ei.graph_) 
                && (curr_root_idx_ == ei.curr_root_idx_) 
                && (leaf_count_ == ei.leaf_count_) );
    }

    /** Return an EdgeIterator that is offset past current iterator
     *  @param @a offset  is the amount to increment the (curr_root_ind_, leaf_count_) of the current iterator
     *
     *  @pre offset must be a valid offset (i.e. EdgeIterator + offset < EdgeIterator.edge_end())
     *  @post initial EdgeIterator is unchanged: a copy with an updated (curr_root_ind_, leaf_count_) is returned
     */
    EdgeIterator& operator+(const size_type offset) {
      auto ei_out = this;
      for (size_type i = 0; i <= offset; ++i)
        ++(*ei_out);
      return *ei_out;
    }

   private:

    friend class Graph;

    // Store a pointer to the graph this EdgeIterator belongs to
    Graph* graph_;

    // Store current root node index (so we only iterate over ACTIVE root nodes)
    size_type curr_root_idx_;

    // Store the iterator value of this EdgeIterator, i.e. what edge count for the current root node
    size_type leaf_count_;

    /** Constructor for EdgeIterator (valid) 
     *  @param[in] curr_root_node  index current root node (default = 0)
     *  @param[in] leaf_count  position of iterator (default = 0)
     *  @param[in] graph      pointer to graph of iterator
     */
    EdgeIterator(const Graph* graph, const size_type curr_root_idx = 0, const size_type leaf_count = 0)
     : graph_(const_cast<Graph*>(graph)), curr_root_idx_(curr_root_idx), leaf_count_(leaf_count) {
    }


  };

  /** Return the EdgeIterator to start iterating over the edges */
  EdgeIterator edge_begin() const {
  	// Start iterating with the first edge
  	return EdgeIterator(this);
  }
  
  /** Return the EdgeIterator to stop iterating over the edges */
  EdgeIterator edge_end() const {
    // Stop iterating after the last edge
  	return EdgeIterator(this, num_nodes());
  }

 private:
  
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Struct for storing the various members of our Nodes */
  struct NodeInfo {

    // Unique identifier for this node
    size_type node_id_;

    // location of node in space
    Point pos_;

    // value of the node
    node_value_type val_;

    /** Constructor for NodeInfo
     *  @param[in] node_id  Unique node identifier of this node
     *  @param[in] pos      Position of the node in space
     *  @param[in] val      Value of the node
     */
    NodeInfo(const size_type node_id, const Point& pos, node_value_type val) 
      : node_id_(node_id), pos_(pos), val_(val) {
    }
  };

  /** Vector to store Node-specific information
   *  Note: this vector includes all nodes ever added to graph, even if they are later made inactive (i.e. removed) */
  std::vector<NodeInfo> graph_nodes_;

  /** Vector to store the node adjacency information  
   *  Note: this vector retains an empty slot for all nodes ever added to graph,
   *        even if they are later made inactive (i.e. removed) */ 
  std::vector<std::vector<size_type>> adj_;

  /** Map to store edge values, indexed by node pair of endpoint node identifiers  
   *  Note: when a node or edge is removed, the corresponding edges ARE removed from this vector */ 
  std::map<std::pair<size_type, size_type>, edge_value_type> edge_values_;

  /** Vector to track the nodes currently "active" (i.e. not removed) in the graph, by holding the unique identifiers of the
   *  active nodes (i.e. the position of their node information in graph_nodes_ and adjacency information in adj_) */
  std::vector<size_type> active_nodes_;

};

#endif // CME212_GRAPH_HPP