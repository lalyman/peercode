#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
/** Graph as a template */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  //struct internal_element;
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  // Predeclare the internal struct
  struct NodeInfo;
  struct EdgeInfo;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node values. */
  using node_value_type = V;

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
    //: sizenode_(0), sizeedge_(0) {
    // HW0: YOUR CODE HERE
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
  /** Derive operator != from a class's operator == 
   * Derive operator>, operator<=, and operator>=.
   */
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_ -> node_[uid_].P_;
	  //return graph_->mynodes[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }
	/** Return this node's value */
	node_value_type& value(){
		return const_cast<V&>(graph_ -> node_[uid_].V_);
	}
	/** Return this node's value */
	const node_value_type& value() const {
		return graph_ -> node_[uid_].V_;
	}

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return the vector that contains uid of adjacent nodes 
     * This is usefule for Problem 3.*/
    std::vector<size_type> adjacent() const {
        return graph_ -> adj_[uid_];
    }
    /** Return the number of edges incident to the current node*/
    size_type degree() const{
        return graph_ -> adj_[uid_].size();
    }
    /** Return incident iterator pointing the current node index*/
    incident_iterator edge_begin() const {
        return IncidentIterator(this, uid_, 0);
    }
    /** Return incident iterator pointing the last incident  node index*/
    incident_iterator edge_end() const {
        return IncidentIterator(this, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
//	  if(this->graph_->node_[this->uid_].P_ == n.graph_->node_[n.uid_].P_ and
//		 this->graph_->node_[this->uid_].V_ == n.graph_->node_[n.uid_].P_)
	  if(this->graph_ == n.graph_ and this->uid_ == n.uid_) 
	  //if(this->uid_ == n.uid_) 
		  return true;
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
	  if (this->graph_ == n.graph_ and this->uid_ < n.uid_)
		  return true;
	  else
		  return false;
    }

   private:
    // HW0: YOUR CODE HERE
	Graph* graph_;
	size_type uid_;
	Node(const Graph* graph, size_type uid)
		: graph_(const_cast<Graph*>(graph)), uid_(uid) {
	}
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	// Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
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
    return node_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * @post adj_.size() should increase by 1
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
	node_.push_back({position, val});
	adj_.push_back(std::vector<size_type>());
	//return Node(this, node_.size() - 1);
	return node(num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
	if( this == n.graph_ and n.uid_ <= node_.size()-1)
		return true;
	else
		return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if( i < num_nodes() ) 
    	return Node(this, i);
	return Node();
  }

  /** Return the node with index j th adjacent node on i th node.
   */
  Node node(size_type i, size_type j) const {
    if( i < num_nodes() and j < node(i).degree() ) 
    	return node(adj_[i][j]);
    else
	    return Node();
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
	  //return Node(graph_, graph_ -> edge_[uid_].uid1_);
	  return graph_ -> node( graph_ -> edge_[uid_].uid1_ );
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
	  //return Node(graph_, graph_ -> edge_[uid_].uid2_);
	  return graph_ -> node( graph_ -> edge_[uid_].uid2_ );
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	  if ( graph_ == e.graph_ and uid_ == e.uid_ and
		  ( (node1() == e.node1() and node2() == e.node2()) or
		    (node2() == e.node1() and node1() == e.node2()) ) )
		  return true;
	  else
		  return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
	  if ( graph_ == e.graph and uid_ < e.uid_ )	 
		 return true;
	  else
	  	 return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
	graph_type* graph_;
	size_type uid_;
	Edge(const graph_type* graph, size_type uid)
		: graph_(const_cast<graph_type*>(graph)), uid_(uid) {
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
    if ( i <= edge_.size()-1 )
		return Edge(this, i);  
	else
		return Edge();
  }

  /** Return the edge that comprise node a and node b.
   */
  Edge edge(const node_type& a, const node_type& b) const {
      for(size_type i = 0; i < num_edges(); ++i) {
		  if ( (edge(i).node1() == a and edge(i).node2() == b) or
			   (edge(i).node1() == b and edge(i).node1() == b) )
		     return Edge(this, i);
      }
      return Edge();  // Invalid Edge
  }
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
	if (has_node(a) and has_node(b)){
		for (size_type i = 0; i < num_edges(); ++i){
			if ( (edge(i).node1() == a and edge(i).node2() == b) or
			     (edge(i).node1() == b and edge(i).node1() == a) )
				return true;
		}
    return false;
	}
	else
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
	if ( has_edge(a, b) or a == b )
		return Edge(); // Invalid Edge
	else {
		edge_.push_back({a.index(), b.index()});
		adj_[a.index()].push_back(b.index());
		adj_[b.index()].push_back(a.index());
		return edge(num_edges() - 1);
	}
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
	node_.clear();
	edge_.clear();
	adj_.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return the Node that the iterator is pointing */
    Node operator*() const { 
        return graph_ -> node(uid_);
	}
    /** Increase the iterator to the next node index */
    NodeIterator& operator++() {
        ++uid_;
        return *this;
    }
    /** Test whether this iterator and @a node_it are equal.
     */
    bool operator==(const NodeIterator& node_it) const {
        if ( graph_ == node_it.graph_ and uid_ == node_it.uid_ )
            return true;
        else
            return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
	graph_type* graph_;
	size_type uid_;
	NodeIterator(const graph_type* graph, const size_type uid)
	  : graph_(const_cast<graph_type*>(graph)), uid_(uid){
	  }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the node iterator pointing the first node index*/ 
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }
  /** Return the node iterator pointing the last node index*/ 
  node_iterator node_end() const {
      return NodeIterator(this, num_nodes());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** */
    Edge operator*() const {
        return graph_ -> edge(graph_ -> node(uid1_),graph_ -> node(uid1_, uid2_));
    }
    /** Increase the adjacent node index */
    IncidentIterator& operator++() {
      ++uid2_;
      return *this; 
    }
    /** Test whether the two incident iterator pointing the same edge.*/
    bool operator==(const IncidentIterator& inci_it) const {
        if ( graph_ == inci_it.graph_ and uid1_ == inci_it and uid2_ == inci_it.uid2_ )
            return true;
        else
            return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
	graph_type* graph_;
	size_type uid1_; // node index where we want to iterate on
	size_type uid2_; // neiboring nodes' index
	IncidentIterator(const graph_type* graph, const size_type uid1, const size_type uid2)
	  : graph_(const_cast<graph_type*>(graph)), uid1_(uid1), uid2_(uid2) {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return the edge that edge iterator is pointing */
    Edge operator*() const {
        return graph_ -> edge(uid_);
    }
    /** Return the iterator after increasing it's index by 1 */
    EdgeIterator& operator++() {
        ++uid_;
        return *this;
    }
    /** Test wether the two iterators are equal */
    bool operator==(const EdgeIterator& edge_it) const {
        if ( graph_ == edge_it.graph_ and uid_ == edge_it.uid_ )
            return true;
        else
            return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	graph_type* graph_;
	size_type uid_;
	EdgeIterator(const graph_type* graph, const size_type uid)
	  : graph_(const_cast<graph_type*>(graph)), uid_(uid){
	  }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the edge iterator pointing the first element of edge container */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }
  /** Return the edge iterator pointing the last element of edge container */
    edge_iterator edge_end() const {
        return EdgeIterator(this, num_edges());
    }

 private:
  // HW0: YOUR CODE HERE
  struct NodeInfo {
	  Point P_;
	  node_value_type V_;
  };

  struct EdgeInfo {
	  size_type uid1_, uid2_;
  };

  std::vector<NodeInfo> node_;
  std::vector<EdgeInfo> edge_;
  std::vector<std::vector<size_type> > adj_;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
