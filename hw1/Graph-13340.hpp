#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>  //TODO: #1
class Graph {
 private:

  // predeclare internal structs that will hold data for nodes and edges
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

  /** Type of node value, a user-specified value */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    :  nodes_(), edges_(), adj_(), next_nuid_(0), next_euid_(0) {
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
      return graph_->nodes_[uid_].position; //-> to access private member
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //return graph_->nodes_[uid_].uid; //an alternate form, useful later?
      return uid_;
    }

    /* Return this node's value */
    node_value_type& value() {
      return graph_->nodes_[uid_].val;
    }

    /* Return this node's value */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].val;
    }

    /* Returns number of edges incident to this Node */
    size_type degree() const {
      return graph_->adj_[uid_].size();
    }

    /* Returns the first edge iterator for an IncidentIterator */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }
    /* Returns the last edge iterator for an IncidentIterator */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check whether uid's are the same
      // Node can access another node's graph_ because Graph is a friend class
      return (this->uid_ == n.index() and this->graph_ == n.graph_);
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
      // check whether nodes are in the same graph	
      // compare uid's, which are ordered by construction
      return (this->uid_ < n.index() and this->graph_ == n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    ///pointer to graph object (will be used for boolean operators)
    Graph* graph_;
    /// This element's unique identification number
    size_type uid_;

    /**Private Constructor accessed by Graph to construct valid Node objects*/
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //returns size of nodes_ vector
    return nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {

    //create a new internal_node - instantiating a struct
    internal_node next_node;

    //set up attributes of new node
    next_node.position = position;
    next_node.val = value;

    nodes_.push_back(next_node);

    //push a placeholder vector in adjacency matrix
    adj_.emplace_back(std::vector<std::tuple<size_type, size_type>>());

    ///increment node counter
    ++next_nuid_;

    return Node(this, next_nuid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //return true if this is the same as graph_ (both are memory addresses)
    return (n.graph_ == this and n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // check that index is within size of graph (i.e. it already exists)
    assert(i < size());
  
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
      //return node at location of first node
      return graph_->node(n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return node at this location of second node
      return graph_->node(n2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // check whether uid and graph are the same
      // Edge can access another Edge's member attributes b/c Graph is friend
      if (graph_ == e.graph_) {
        return ((e.node1() == node1() and e.node2() == node2())
                or (e.node1() == node2() and e.node2() == node1()));
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // check that the edges are in same graph, and compare uid's
      if (graph_ == e.graph_) {
        return (node1() < e.node1() or (node1() == e.node1() and node2() < e.node2()));
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer to Graph container
    Graph* graph_;
    // This element's unique identification number
    size_type uid_;
    // Unique identification number of one node
    size_type n1_uid_;
    // Unique identification number of the other node
    size_type n2_uid_;
    /** Private constructor accessed by Graph to construct valid Edge objects */
    Edge(const Graph* graph, size_type uid, size_type n1_uid, size_type n2_uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid), n1_uid_(n1_uid), n2_uid_(n2_uid)  {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    //returns size of edges_ vector
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    internal_edge new_edge = edges_[i]; //access ith edge in vector
    //define n1 and n2 of edge
    size_type n1_idx = new_edge.node1_uid; 
    size_type n2_idx = new_edge.node2_uid; 
    //pass to constructor 
    return Edge(this, i, n1_idx, n2_idx);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check that a and b are in the graph
    if (this != a.graph_ or this != b.graph_) {
      return false;
    }
    //check internal vector of adjacency matrix to see if b is in a's vector
    for (size_type i = 0; i < adj_[a.index()].size(); i++) {
      // check first element of tuple (i.e. uid for node)
      if (std::get<0>(adj_[a.index()][i]) == b.index()) { 
        return true;
      }
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

    //check whether edge already exists
    if (has_edge(a, b)) {

      //loop through a's internal vector in adj_, return edge of matching node
      for (size_type i = 0; i < adj_[a.index()].size(); i++) {
        if (std::get<0>(adj_[a.index()][i]) == b.index()) { 
          return edge(std::get<1>(adj_[a.index()][i]));    
        }
      } 
    }

    //create a new internal_node - instantiate a struct
    internal_edge next_edge;
    //set up attributes of new edge
    next_edge.node1_uid = a.index();
    next_edge.node2_uid = b.index();
    //add new edge to vector
    edges_.push_back(next_edge);
    //update edge counter
    ++next_euid_;

    //add entries to adjacency matrix
    adj_[a.index()].push_back(std::make_tuple(b.index(), next_euid_ - 1));
    adj_[b.index()].push_back(std::make_tuple(a.index(), next_euid_ - 1));

    //return a new edge object
    return Edge(this, next_euid_ - 1, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // use vector clear to remove all data from nodes_, edges_, and adj_
    nodes_.clear();
    edges_.clear();
    adj_.clear();
    // reset counters
    next_nuid_ = 0;
    next_euid_ = 0;
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
    /* Dereferences the node iterator to return a Node object*/
    Node operator*() const {
      return graph_->node(index_);
    }
    /* Increments the node iterator*/
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }
    /* Checks whether this NodeIterator equals x */ 
    bool operator==(const NodeIterator& x) const {
      return (this->graph_ == x.graph_ and this->index_ == x.index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // pointer to graph object
    Graph* graph_;
    // counter for node iterator class
    size_type index_;
    /* Private constructor for node_iterator class */
    NodeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /* Returns the first node iterator  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /* Returns the last node iterator */
  node_iterator node_end() const {
    return NodeIterator(this, this->size());
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

    // HW1 #3: YOUR CODE HERE
    /* Dereferences the incident iterator to return an Edge object */ 
    Edge operator*() const {
      size_type edge_index = std::get<1>(graph_->adj_[node1_uid_][index_]);
      return graph_->edge(edge_index);
    }
    /* Increments the incident iterator */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }
    /* Check whether this IncidentIterator equals x */
    bool operator==(const IncidentIterator& x) const {
      return (this->graph_ == x.graph_ and this->index_ == x.index_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // pointer to graph object
    Graph* graph_;
    // uid of node that spawned the incident iterator
    size_type node1_uid_;
    // counter for incident iterator class
    size_type index_;
    /* Private constructor for incident_iterator class */
    IncidentIterator(const Graph* graph, size_type node1_uid, size_type index)
      : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), index_(index) {
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
    /* Dereferences the edge iterator to return an Edge object */
    Edge operator*() const {
      return graph_->edge(index_);
    }
    /* Increments the edge iterator */
    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }
    /* Checks whether this EdgeIterator is the same as x */
    bool operator==(const EdgeIterator& x) const {
      return (this->graph_ == x.graph_ and this->index_ == x.index_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    //pointer to graph object
    Graph* graph_;
    //counter for edge iterator class
    size_type index_;
    /* Private constructor for EdgeIterator class */
    EdgeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  /* Returns the first edge iterator */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /* Returns the last edge iterator */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // contains data attributes defining a node
  struct internal_node {
    Point position;
    node_value_type val;
  };
  // contains data attributes defining an edge:
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
  };

  // STL containers for internal_nodes and internal_edges
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<std::vector<std::tuple<size_type, size_type>>> adj_;

  //counters for constructor
  size_type next_nuid_; //node counter
  size_type next_euid_; //edge counter

};

#endif // CME212_GRAPH_HPP
