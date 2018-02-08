#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
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

  /** Synonym for Value */
  using node_value_type = V;

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
  class Node : private totally_ordered<Node>{
   public:
    
    Node() {
    }

    const Point& position() const {
    	if(nid_ < graph_->node_.size())
			return graph_->node_[nid_];

		assert(false);
    }

    size_type index() const {
    	if( nid_ < graph_->node_.size() )
	    	return nid_; 

		assert(false);
    }

    // HW1: YOUR CODE HERE
    /** 
    *  Get the user-specified value of each node
	*/
    node_value_type& value() {
    	if( nid_ < graph_->node_.size() )
	    	return graph_->nodevalue_[nid_]; 

		assert(false);
    }

    const node_value_type& value() const{
    	if( nid_ < graph_->node_.size() )
	    	return graph_->nodevalue_[nid_];

		assert(false);
    }

    void set_value(node_value_type vlu) {
    	val_ = vlu;
    }

    size_type degree() const{
    	return graph_->connect.at(nid_).size();
    }

    incident_iterator edge_begin() const{
    	return IncidentIterator(graph_, 0, nid_);
    }

    incident_iterator edge_end() const{
    	return IncidentIterator(graph_, degree(), nid_);
    }

    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
    	if(nid_ == n.nid_ and graph_ == n.graph_) 
    		return true;
    	else return false;
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
    	if((nid_ < n.nid_ and graph_ == n.graph_) || graph_<n.graph_) 
    		return true;
    	else return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_;
    size_type nid_;
    node_value_type val_;

    // private constructor
    Node(const Graph* graph, size_type nid)
	        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
	}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
		return node_.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type ()) {
    // HW0: YOUR CODE HERE
  	node_.push_back(position);
  	nodevalue_.push_back(val);
  	//add a key in connect mapping and link it to an empty vector
  	std::vector<size_type> emp;
  	connect[node_.size()-1] = emp;
	return Node(this, node_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    // if the node is not in this graph
    if(n.graph_!= this) 
		return false;
	// otherwise,
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //assert(i < size());
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
    	return Node(graph_, aidx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
    	return Node(graph_, bidx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	    if(graph_ == e.graph_ and eid_ == e.eid_) 
	    	return true;
	    else return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	if((eid_ < e.eid_ and graph_ == e.graph_ ) || graph_<e.graph_) 
    		return true;
    	else return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions. 
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_;
    size_type eid_;
    size_type aidx_;
    size_type bidx_;

    // private constructor
    Edge(const Graph* graph, size_type eid, size_type aidx, size_type bidx)
	        : graph_(const_cast<Graph*>(graph)), eid_(eid), 
	        aidx_(aidx), bidx_(bidx) {
	    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
	if(i >= edge_.size()) {
        return Edge();
    }
    return Edge(this, i, edge_[i][0], edge_[i][1]) ;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (auto e : connect.at(a.index())){
    	if (edge_[e][0]==b.index() or edge_[e][1]==b.index())
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
  	if (a == b) 
  		return Edge();
  	
  	if(has_edge(a, b)){
		for (auto e : connect.at(a.index())) {
			if (this->edge(e).node1() == a and this->edge(e).node2() == b) 
    		// if (edge_[e][0]==b.index() or edge_[e][1]==b.index())
    			return edge(e);
    		if (this->edge(e).node2() == a and this->edge(e).node1() == b)
    			return Edge(this,edge(e).eid_,a.index(),b.index());
		}
  	}

	std::vector<size_type> temp;
	temp.push_back(a.index());
	temp.push_back(b.index());
	edge_.push_back(temp);

  	connect[a.index()].push_back(edge_.size()-1);
  	connect[b.index()].push_back(edge_.size()-1);

  	return Edge(this,edge_.size()-1,a.index(),b.index());
	
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
  	node_.clear();
  	edge_.clear();
  	connect.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    /** 
  	* Dereference the node iterator
  	* @return The node with current node_id
    */
    Node operator*() const{
    	if (iid_ < graph_->size())
	    	return graph_->node(iid_);

	    assert(false);
    }

    /** 
	  * Increment the current node_id and returns the node iterator
    */
    NodeIterator& operator++() {
    	iid_++;
    	return *this;
    }

    /** 
  	* Check if two edge iterators are equal by comparing their edge_id 
  	* and if they are in the same graph.
    */
    bool operator==(const NodeIterator& i) const {
    	return (iid_ == i.iid_ and graph_ == i.graph_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type iid_;

    // private constructor
    NodeIterator(const Graph* graph, size_type iid)
	        : graph_(const_cast<Graph*>(graph)), iid_(iid) {
	}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** 
  * Returns a node iterator with node id 0, the begin of a node iterator
  */
  node_iterator node_begin() const{
  		return NodeIterator(this, 0);
  }
  /** 
  * Returns the end of node iterator
  */
  node_iterator node_end() const{
  		return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. 
   */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** 
	* Dereference the incident iterator
	* @return The edge with one end connected to node with node_id nidx
    */
    Edge operator*() const {
    	if (incid_ < graph_->node(nidx_).degree()) {
			size_type temp_eid = graph_->connect.at(nidx_)[incid_];
			return graph_->edge(temp_eid);
    	}

    	assert(false);
    }

    /** 
	* Increment the current edge_id 
	* @return The incident iterator representing the incident edge
    */
    IncidentIterator& operator++() {
    	incid_++;
    	return *this;
    }

    /** 
	* Check if two incident iterators are equal by comparing their host node id, 
	* positions in the mapping vector and if they are in the same graph
    */
    bool operator==(const IncidentIterator& i) const {
    	return (incid_ == i.incid_ and graph_ == i.graph_ 
    		and nidx_ == i.nidx_); 
    }

   private:
    friend class Graph; 
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type incid_; // index/position in the vector of edge_id 
    size_type nidx_;

    // private constructor
    IncidentIterator(const Graph* graph, size_type incid, size_type nidx)
		: graph_(const_cast<Graph*>(graph)), incid_(incid), nidx_(nidx) {
	}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()  {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

	/** 
	* Dereference the edge iterator
	* @return The edge with edge_id eiid_
    */
    Edge operator*() const {
    	if (eiid_<graph_->num_edges())
	    	return graph_->edge(eiid_);

	    assert(false);    
	}

	/** 
	* Increment the current edge_id 
	* @return The edge iterator 
    */
    EdgeIterator& operator++() {
    	eiid_++;
    	return *this;
    }

    /** 
	* Check if two edge iterators are equal by comparing their node_id, 
	* and if they are in the same graph
    */
    bool operator==(const EdgeIterator& e) const {
    	return (graph_ == e.graph_ and eiid_ == e.eiid_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type eiid_; 

    // private constructor
    EdgeIterator(const Graph* graph, size_type eiid)
     : graph_(const_cast<Graph*>(graph)), eiid_(eiid) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** 
  * Returns an edge iterator with edge id 0, the begin of the edge iterator
  */
  edge_iterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }
  /** 
  * Returns the end of edge iterator
  */
  edge_iterator edge_end() const {
  	return EdgeIterator(this,num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

 	std::vector<Point> node_;
 	std::vector<node_value_type> nodevalue_;
 	std::vector<std::vector<size_type>> edge_;
	// Mapping from one node (node_id) to the a set of 
	// edges (edge_id) that has one end on this node.
 	std::map<size_type,std::vector<size_type>> connect;
};

#endif // CME212_GRAPH_HPP

