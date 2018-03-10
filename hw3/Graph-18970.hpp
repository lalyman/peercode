#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * Using peer code 23158 as starting point
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
template <typename V, typename E>
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
  class Node ;
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

    Point& position() {
    	return graph_->nodes_[node_uid_].position;
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

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    node_value_type& value(node_value_type val) { // argument added 
    	graph_->nodes_[node_uid_].value = val;
    	return graph_->nodes_[node_uid_].value;
    }

    const node_value_type& value() const {
    	return graph_->nodes_[node_uid_].value;
    }

    node_value_type& value() {
    	return graph_->nodes_[node_uid_].value;
    }

    size_type degree() const {
    	return graph_->adj[node_uid_].size();
    }
    
    incident_iterator edge_begin() const {
    	return IncidentIterator(graph_, node_uid_, 0);
    }
    
    incident_iterator edge_end() const {
    	return IncidentIterator(graph_, node_uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /*** Test if this node's ID and graph ptr are equivalent to n's ***/
      if ((n.node_uid_ == this->node_uid_) && (n.graph_ == this->graph_)) {
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
      if (n.node_uid_ > this->node_uid_) {
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
    return int_id.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    size_type i = nodes_.size();
    size_type j = int_id.size();
    internal_node n1(j, position, val);
    nodes_.push_back(n1);
    int_id.push_back(i);
    std::vector<internal_edge> adj_edges;
    adj.push_back(adj_edges);
    return Node(this, i);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /*** same graph if ptr to memory for node is the same as graph_ ptr ***/ 
    if (this == n.graph_ && int_id[n.index()] == n.node_uid_) {
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
    // assert(i < num_nodes());
    return Node(this, int_id[i]);
  }


  // additional node methods for HW2 // 

  size_type remove_node(const Node& n) {
  	for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
  		Edge e = *it;
  		remove_tail_node(e.node2(), n);
  	}
  	adj[n.node_uid_].clear();

  	size_type k = n.index();
  	for (size_type i = k; i < int_id.size(); i++) {
  		--nodes_[int_id[i]].node_uid;
  	}

  	int_id.erase(int_id.begin() + k);

  	return k;
  }

  node_iterator remove_node(node_iterator n_it) {
  	Node n = *n_it;
  	return NodeIterator(this, remove_node(n_it));
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
      /*** access this edge's node1 member through graph_ ptr ***/
      return Node(graph_, n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      /*** access this edge's node2 member through graph_ ptr ***/
      return Node(graph_, n2_uid_);
    }

    // adding functions for HW2

    edge_value_type& value() {
    	size_type id;
    	if (n1_uid_ > n2_uid_) {
    		id = graph_->match_tail_node(Node(graph_, n2_uid_), Node(graph_, n1_uid_));
    		return graph_->adj[n2_uid_][id].value;
    	} else {
    		id = graph_->match_tail_node(Node(graph_, n1_uid_), Node(graph_, n2_uid_));
    		return graph_->adj[n1_uid_][id].value;
    	}
    }

    const edge_value_type& value() const {
    	return value(); 
    }

    double length() {
    	return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool isgraphsame = (graph_ == e.graph_);
      bool c1 = (n1_uid_ == e.n1_uid_) && (n2_uid_ == e.n2_uid_);
      bool c2 = (n1_uid_ == e.n2_uid_) && (n2_uid_ == e.n1_uid_);
      return isgraphsame && (c1 || c2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ && n1_uid_ == e.n1_uid_) 
      	return n2_uid_ < e.n2_uid_;
      if (graph_ == e.graph_ && n1_uid_ != e.n1_uid_) 
      	return n1_uid_ < e.n1_uid_;

      return graph_ < e.graph_;

    }

   private:
    /*** Private data members and methods for Edge ***/
    // Pointer back to the Graph container
    Graph* graph_;
   
    size_type n1_uid_;
    size_type n2_uid_;
    /* private constructor Graph::Edge */
    Edge(const Graph* g, size_type id1, size_type id2)
      : graph_(const_cast<Graph*>(g)), n1_uid_(id1), n2_uid_(id2) {
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
    // more efficient than storing a variable for it
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
    size_type i = match_tail_node(a, b);
    return i != a.degree();
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
    
    size_type id1 = a.node_uid_;
    size_type id2 = b.node_uid_;

    // add if edge already not there
    if (!has_edge(a, b)) {
    	internal_edge e1(id2, edge_value_type());
    	internal_edge e2(id1, edge_value_type());
    	adj[id1].push_back(e1);
    	adj[id2].push_back(e2);
    }
       
    return Edge(this, id1, id2);   
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    /*** clear nodes_ and edges_ vectors ***/ 
    nodes_.clear();
    adj.clear();
    int_id.clear();
  }


  // additional methods for HW2

  size_type remove_edge (const Node& n1, const Node& n2) { 
  	if (!has_edge(n1, n2)) {
  		if (n1.node_uid_ < n2.node_uid_)
  			return n1.degree();
  		else 
  			return n2.degree();
  	}

  	size_type i = remove_tail_node(n1, n2);
  	size_type j = remove_tail_node(n2, n1);

  	if (n1.node_uid_ < n2.node_uid_)
  		return i;
  	else 
  		return j;

  }

  size_type remove_edge ( const Edge& e) {
  	return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge (edge_iterator e_it) {
  	Edge e = *e_it;
  	size_type col = remove_edge(e.node1().node_uid_, e.node2().node_uid_);
  	return EdgeIterator(this, e.node1().node_uid_, col);
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
    Node operator*() const {
    	return graph_->node(itr_ind);
    }

    NodeIterator& operator++() {
    	if (itr_ind < graph_->size()) {
    		itr_ind = itr_ind + 1;
    	}
    	return *this;
    }

    bool operator==(const NodeIterator& nd) const { 
    	// graph and index should be same 
    	return ((graph_ == nd.graph_) && (itr_ind == nd.itr_ind));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type itr_ind;
    NodeIterator(const Graph* g, size_type i = 0) : 
    	graph_(const_cast<Graph*>(g)), itr_ind(i) {

    	}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
  	return NodeIterator(this, 0); // node_iterator is synonym for NodeIterator
  }

  node_iterator node_end() const {
  	return NodeIterator(this, this->size());
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
    Edge operator*() const {
    	size_type id1 = row_;
    	size_type id2 = graph_->adj[row_][col_].edge_uid;
    	return Edge(graph_, id1, id2);
    }


    IncidentIterator& operator++() {
    	++col_;
    	return *this;
    }
  
    bool operator==(const IncidentIterator& nd) const {
    	return ((graph_ == nd.graph_) && (row_ == nd.row_) && 
    		(col_ == nd.col_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type row_;
    size_type col_; 

    IncidentIterator(const Graph* g, size_type ro, size_type co) : 
    	graph_(const_cast<Graph*>(g)), row_(ro), col_(co) {

    }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
    	size_type id1 = row_;
    	size_type id2 = graph_->adj[row_][col_].edge_uid;
    	return Edge(graph_, id1, id2);
    }

    EdgeIterator& operator++() {
    	++col_;
    	shift();
    	return *this;
    }
    
    bool operator==(const EdgeIterator& nd) const { 
    	return ((graph_ == nd.graph_) && (row_ == nd.row_) && 
    		(col_ == nd.col_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type row_;
    size_type col_;

    EdgeIterator(const Graph* g, size_type ro) : 
    	graph_(const_cast<Graph*>(g)), row_(ro), col_(0) {
    		shift();
    	}

    EdgeIterator(const Graph* g, size_type ro, size_type co) : 
    	graph_(const_cast<Graph*>(g)), row_(ro), col_(co) {
    		shift();
    	}

    // helper function to shift the edge iterator 
    void shift() {
    	// chk if row and col index are valid
    	while(row_ < graph_->adj.size()) {
    		while(col_ < graph_->adj[row_].size()) {
    			size_type id1 = row_;
    			size_type id2 = graph_->adj[row_][col_].edge_uid;
    			if (id1 < id2) 
    				return; // no operation required
    			++col_;
    		}
    		++row_; // shift to next row 
    		col_ = 0; // reset 
    	}
    	return;

    }


  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
  	return EdgeIterator(this, adj.size());
  }

 private:

 	// additional internal helper methods 

  template<typename T>
  void erase(std::vector <T> &vec, unsigned int i) {
    vec[i] = vec.back();
    vec.pop_back();
  }

 	size_type match_tail_node(const Node& h, const Node& t) const {
 		size_type id = t.node_uid_;
 		auto match = [&id](Edge e) {
 			return e.node2().node_uid_== id;
 		}; 
 		IncidentIterator it = std::find_if(h.edge_begin(), h.edge_end(), match);
 		// return the distance from the head to it
 		return std::distance(h.edge_begin(), it);

 	}

 	size_type remove_tail_node(const Node& h, const Node& t) {
 		size_type id = match_tail_node(h, t);
 		if (id != h.degree()) {
 			erase(adj[h.node_uid_], id);
 		}
 		return id;
 	}


  /*** Graph class internals for Node and Edge proxys ***/
  // The following section was significantly changed to make the code leight weight
  struct internal_node {
  	size_type node_uid;
    Point position;
    node_value_type value;

    internal_node(size_type i, const Point& p, const node_value_type& v) :
    node_uid(i), position(p), value(v) {

    }
  };

  struct internal_edge {
    // Removing storing two nodes to make it light weight
    // Node node1;
    // Node node2;
    size_type edge_uid;
    edge_value_type value;

    internal_edge(size_type id, const edge_value_type& v) : edge_uid(id), value(v) {

    }

  };

  /*** Graph class vectors of internals ***/
  std::vector<internal_node> nodes_;
  std::vector<std::vector<internal_edge>> adj;
  std::vector<size_type> int_id;
  // std::vector<internal_edge> edges_;
  // added for HW1
  // std::vector<std::vector<std::pair<size_type, size_type>>> adjacency_list_; 

  /*** Graph class data members ***/ 
  // size_type num_nodes_;
  // size_type num_edges_;
  // size_type next_node_uid_;
  // size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP
