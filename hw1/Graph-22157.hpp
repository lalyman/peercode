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
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  // Predeclare the internal structs.
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Type of the node value. */
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
  Graph()
    : internal_nodes_(), internal_edges_(), adj_() {
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

    /** Return this node's position.
     * @pre This Node is valid and belongs to a Graph.
     */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * @pre This Node is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().uid;
    }
    
    /** Return the number of edges incident to this node.
     * @pre This Node is valid and belongs to a Graph.
     */
    size_type degree() const {
      return graph_->adj_[uid_].size(); 
    }
    
    /** Return the beginning incident iterator for this node.
     * This incident iterator corresponds to the first edge
     * in the list of edges incident to this node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(0, uid_, graph_); 
    }

    /** Return the end incident iterator for this node..
     * This incident iterator corresponds to the last edge
     * in the list of edges incident to this node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(degree(), uid_, graph_);
    }

    /** Return the value of this node.
     * Note that this is returned by reference and can be
     * modified.
     */
    node_value_type& value() {
      return fetch().value;
    }

    const node_value_type& value() const {
      return fetch().value;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /** Check if (1) graph to which n belongs is same as
       * this nodes's graph and (2) the two nodes have the same
       * index/uid.
       */
      if ((graph_ == n.graph_) && (uid_ == n.index())) {
        return true;
      }
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
      Graph* n_graph = n.graph_;
      if (graph_ < n_graph) {
        // This node's graph is less than n's graph (by pointer comparison).
        return true;
      }

      if (graph_ > n_graph) {
        // This node's graph is larger than n's graph (by pointer comparison).
        return false;
      }

      // Otherwise the graphs are the same (by pointer comparison) so compare
      // the node index/uid.
      return (uid_ < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to Graph to which this node belongs.
    Graph* graph_;

    // This node's unique ID. The unique ID is the same as the node's index.
    // Both correspond to the order in which nodes where added to the graph.

    // Important note: To support node removal in the future, the uid and
    // index will have to be differentiated.
    size_type uid_;

    // Private Constructor.
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    // Helper method to return appropriate internal_node.
    internal_node& fetch() const {
      return graph_->internal_nodes_[uid_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return internal_nodes_.size();
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
    node_value_type default_value = {};
    return add_node(position, default_value);
  }
 
  Node add_node(const Point& position, const node_value_type& value) {

    size_type current_size = internal_nodes_.size();

    // Create new internal_node object to add.
    internal_node new_node;
    new_node.position = position;
    new_node.value = value;
    new_node.uid = current_size;

    // Add to list of nodes.
    internal_nodes_.push_back(new_node);

    // Return a Node that points to newly added node.
    return Node(this, current_size);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type idx = n.index();

    // If the index of n is larger than this graph's
    // size (# nodes) then n cannot belong to the graph.
    if (internal_nodes_.size() < idx) {
      return false;
    }

    // Otherwise, check if n and the node in this graph with
    // the same index as n are equal.
    Node m = node(idx);
    return (n == m);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
      size_type node1_uid = fetch().node1_uid;
      return graph_->node(node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type node2_uid = fetch().node2_uid;
      return graph_->node(node2_uid);
    }

    /** Return index of edge, a number in the range [0, num_edges).
     * @pre This Edge is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().edge_uid;
    }

    /** Return pointer to graph to which this Edge belongs */
    Graph* graph() const {
      return graph_;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if (1) graph to which this edge belongs is same as
      // current graph and (2) the two edges have the same uid.

      // Duplicate edges are checked for when attempting to add edges to a
      // graph, so two edges with the same index/uid that have ended up in
      // the same graph must be equal (represent the same undirected edge
      // between two nodes.)
      Graph* e_graph = e.graph();
      if ((graph_ == e_graph) && (edge_uid_ == e.index())) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      Graph* e_graph = e.graph();
      if (graph_ < e_graph) {
        // This edges's graph is less than e's graph (by pointer comparison).
        return true;
      }

      if (graph_ > e_graph) {
        // This edge's graph is larger than e's graph (by pointer comparison).
        return false;
      }

      // Otherwise the graphs are equal so compare the edge index/uid.
      return (edge_uid_ < e.index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph to which this edge belongs.
    Graph* graph_;

    // This edge's unique ID. The unique ID is the same as the edge's
    // index. Both correspond to the order in which edges were added
    // to the graph.

    // Important note: To support edge removal in the future, the uid and
    // index will have to be differentiated.
    size_type edge_uid_;

    // Private constructor.
    Edge(const Graph* graph, size_type edge_uid)
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }
    
    /** Swap node1 and node2 for this edge. The graph is undirected and
     *  so switching which node is considered node1 with that for node2
     *  does not change how we can use the edge. However, if a particular
     *  ordering of the nodes in an edge needs to be enforced in the
     *  future, this method should be reworked.
     */
    void swap_nodes() const {
      internal_edge& internal_edge = graph_->internal_edges_[edge_uid_];
      size_type n1_uid = internal_edge.node1_uid;
      size_type n2_uid = internal_edge.node2_uid;

      internal_edge.node1_uid = n2_uid;
      internal_edge.node2_uid = n1_uid; 
    }

    // Helper method to return appropriate internal_edge.
    internal_edge& fetch() const {
      return graph_->internal_edges_[edge_uid_];
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return internal_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < internal_edges_.size(); ++i) {
      internal_edge current_edge = internal_edges_[i];
      size_type current_edge_node1_uid = current_edge.node1_uid;
      size_type current_edge_node2_uid = current_edge.node2_uid;

      // Edges are undirected, so to check if edge already exists
      // we compare against (a,b) and (b,a).
      if ((a.index() == current_edge_node1_uid) &&
		(b.index() == current_edge_node2_uid)) {
        return true;
      }

      if ((b.index() == current_edge_node1_uid) &&
		(a.index() == current_edge_node2_uid)) {
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
    // Add edge only if edge does not already exists.
    if (!(has_edge(a, b))) {
      size_type current_num_edges = internal_edges_.size();	

      // Create new internal_edge object to add.
      internal_edge new_edge;
      size_type node1_uid = a.index();
      size_type node2_uid = b.index();

      new_edge.node1_uid = node1_uid;
      new_edge.node2_uid = node2_uid;
      new_edge.edge_uid = current_num_edges;

      // Add to list of edges.
      internal_edges_.push_back(new_edge);

      // Add edge to { node_uid : [adjacent edge uids] } adjacency map.
      // First node.
      if (adj_.find(node1_uid) == adj_.end()) {
        adj_.insert(std::pair<size_type,
		std::vector<size_type>>(node1_uid, {current_num_edges}));
      } else {
        adj_[node1_uid].push_back(current_num_edges);  
      }
      // Second node.
      if (adj_.find(node2_uid) == adj_.end()) {
        adj_.insert(std::pair<size_type,
		std::vector<size_type>>(node2_uid, {current_num_edges}));
      } else {
        adj_[node2_uid].push_back(current_num_edges);  
      }

      // Return an Edge that points to newly added edge.
      return Edge(this, current_num_edges);
    }
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_nodes_.clear();
    internal_edges_.clear();
    adj_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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
    
    /** Construct a valid NodeIterator. */
    NodeIterator(size_type i, Graph* g) : current_{i}, graph_{g} {}

    // OPERATORS.

    Node operator*() const {
      return Node(graph_, current_); 
    }

    node_iterator& operator++() {
      ++current_;
      return *this;
    }

    bool operator==(const node_iterator& rhs) const {
      return (graph_ == rhs.graph_) && (current_ == rhs.current_);
    }
    
    bool operator!=(const node_iterator& rhs) const {
      return (graph_ != rhs.graph_) || (current_ != rhs.current_);
    }
    
   private:
    friend class Graph;
    size_type current_;
    Graph* graph_;
  };

  Graph* this_graph = this;

  /** Return the beginning node iterator.
   *  This node iterator corresponds to the first node
   *  added to the graph.
   */
  node_iterator node_begin() const {
    return NodeIterator(0, this_graph);
  }

  /** Return the end node iterator.
   *  This node iterator corresponds to the last node
   *  added to the graph.
   */
  node_iterator node_end() const {
    return NodeIterator(internal_nodes_.size(), this_graph);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}
    
    /** Construct a valid IncidentIterator. */
    IncidentIterator(size_type i, size_type n, Graph* g) :
      current_{i}, spawn_node_uid_{n}, graph_{g}, adj_edges_{g->adj_[n]} {}

    // OPERATORS.

    Edge operator*() const {
      size_type edge_uid = adj_edges_[current_];
      Edge edge = Edge(graph_, edge_uid);
      
      // Make sure spawn node always returned by node1(). 
      if (edge.node2().index() == spawn_node_uid_) {
        edge.swap_nodes();
      }      

      return edge;
    }
    
    IncidentIterator& operator++() {
      ++current_;
      return *this;
    }
    
    bool operator==(const IncidentIterator& rhs) const {
      return (graph_ == rhs.graph_) && (current_ == rhs.current_)
	&& (spawn_node_uid_ == rhs.spawn_node_uid_);
    }
    
    bool operator!=(const IncidentIterator& rhs) const {
      return (graph_ != rhs.graph_) || (current_ != rhs.current_)
	|| (spawn_node_uid_ != rhs.spawn_node_uid_);
    }

   private:
    friend class Graph;
      size_type current_;
      size_type spawn_node_uid_;
      Graph* graph_;
      std::vector<size_type> adj_edges_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
  
    /** Construct a valid EdgeIterator. */
    EdgeIterator(size_type i, Graph* g) : current_{i}, graph_{g} {}

    // OPERATORS.

    Edge operator*() const {
      return Edge(graph_, current_);
    }

    EdgeIterator& operator++() {
      ++current_;
      return *this;
    }
    
    bool operator==(const EdgeIterator& rhs) const {
      return (graph_ == rhs.graph_) && (current_ == rhs.current_);
    }
    
    bool operator!=(const EdgeIterator& rhs) const {
      return (graph_ != rhs.graph_) || (current_ != rhs.current_);
    }

   private:
    friend class Graph;
    size_type current_;
    Graph* graph_;
  };

  /** Return the beginning edge iterator for this graph.
   *  This edge iterator corresponds to the first edge
   *  added to the graph.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this_graph);
  }
  
  /** Return the end edge iterator for this graph.
   *  This edge iterator corresponds to the first edge
   *  added to the graph.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(internal_edges_.size(), this_graph);
  }

 private:

  // Graph class's internals helper functions, data members, and so forth.

  struct internal_node {
    // The position of the node.
    Point position;
    // The value of the node.
    node_value_type value;
    // The uid/index of the node (based on order in which nodes were added).
    size_type uid;
  };

  struct internal_edge {
    // The indices/uids of the pair of nodes connected by this edge.
    size_type node1_uid;
    size_type node2_uid;
    // The uid/index of the edge (based on order in which edges were added).
    size_type edge_uid;
  };

  std::vector<internal_node> internal_nodes_;  // Container for nodes in graph.
  std::vector<internal_edge> internal_edges_;  // Container for edges in graph.

  // Container for map of node uids to a list of uids of incident edges.
  std::map<size_type, std::vector<size_type>> adj_;
};

#endif // CME212_GRAPH_HPP
