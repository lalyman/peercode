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
template <typename V, typename E>
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
  using graph_type = Graph<V, E>;

  /** Type of the node value. */
  using node_value_type = V;
  
  /** Type of the edge value. */
  using edge_value_type = E;

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
    : internal_nodes_(), internal_edges_(), adj_(), i2u_nodes_(), i2u_edges_() {
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
    
    /** Return this node's position.
     * @pre This Node is valid and belongs to a Graph.
     */
    Point& position() {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * @pre This Node is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().idx;
    }
    
    /** Return the number of edges incident to this node.
     * @pre This Node is valid and belongs to a Graph.
     */
    size_type degree() const {
      return graph_->adj_[uid_].size(); 
    }
    
    /** Return pointer to graph to which this Node belongs */
    Graph* graph() const {
      return graph_;
    }
    
    /** Return node uid. */
    size_type uid() const {
      return uid_;
    }
    
    /** Return the beginning incident iterator for this node.
     * This incident iterator corresponds to the first edge
     * in the list of edges incident to this node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(0, Node(graph_, uid_), graph_); 
    }

    /** Return the end incident iterator for this node..
     * This incident iterator corresponds to the last edge
     * in the list of edges incident to this node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(degree(), Node(graph_, uid_), graph_);
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
       * uid.
       */
      if ((graph_ == n.graph_) && (uid_ == n.uid())) {
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
      return (uid_ < n.uid());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to Graph to which this node belongs.
    Graph* graph_;

    // This node's unique ID.
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
    return i2u_nodes_.size();
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
    
    size_type current_internal_size = internal_nodes_.size();
    size_type current_registered_size = i2u_nodes_.size();

    // Create new internal_node object to add.
    internal_node new_node;
    new_node.position = position;
    new_node.value = value;
    new_node.uid = current_internal_size;
    new_node.idx = current_registered_size;
    
    // Add to internal list of nodes.
    internal_nodes_.push_back(new_node);

    // Add to index to uid vector.
    i2u_nodes_.push_back(current_internal_size); 

    // Add entry adjacency list to be filled with this nodes neighbors.
    adj_.push_back({});

    // Return a Node that points to newly added node.
    return Node(this, current_internal_size);
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
    if (i2u_nodes_.size() < idx) {
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
    size_type uid = i2u_nodes_[i];
    return Node(this, uid);
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
      //return graph_->node(node1_uid);
      return Node(graph_, node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type node2_uid = fetch().node2_uid;
      //return graph_->node(node2_uid);
      return Node(graph_, node2_uid);
    }

    /** Return index of edge, a number in the range [0, num_edges).
     * @pre This Edge is valid and belongs to a Graph.
     */
    size_type index() const {
      return fetch().idx;
    }
    
    /** Return the value of this node. */
    const edge_value_type& value() const {
      return fetch().value;
    }

    /** Return pointer to graph to which this Edge belongs */
    Graph* graph() const {
      return graph_;
    }
   
    size_type uid() const {
      return edge_uid_;
    }
    
    /** Return length of edge (euclidean distance between the two end nodes) */
    double length() const {
      Node node1 = this->node1();
      Node node2 = this->node2();
      return norm(node1.position() - node2.position());
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
      if ((graph_ == e_graph) && (edge_uid_ == e.uid())) {
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
      return (edge_uid_ < e.uid());
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
    return i2u_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type edge_uid = i2u_edges_[i];
    return Edge(this, edge_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return Uid of edge if for some @a i, edge(@a i) connects @a a and @a b.
   *         If not, returns -1;
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  int has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < i2u_edges_.size(); ++i) {
      size_type edge_uid = i2u_edges_[i];
      internal_edge current_edge = internal_edges_[edge_uid];
      size_type current_edge_node1_uid = current_edge.node1_uid;
      size_type current_edge_node2_uid = current_edge.node2_uid;

      // Edges are undirected, so to check if edge already exists
      // we compare against (a,b) and (b,a).
      if ((a.uid() == current_edge_node1_uid) &&
		(b.uid() == current_edge_node2_uid)) {
        return edge_uid;
      }

      if ((b.uid() == current_edge_node1_uid) &&
		(a.uid() == current_edge_node2_uid)) {
        return edge_uid;
      }
    }

    return -1;
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
    edge_value_type default_value = {};
    return add_edge(a, b, default_value);
  }

  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value) {
    // Add edge only if edge does not already exists.
    int has_edge_uid = has_edge(a,b);
    if (has_edge_uid == -1) {

      size_type current_internal_edges = internal_edges_.size();	
      size_type current_registered_edges = i2u_edges_.size();	

      // Create new internal_edge object to add.
      internal_edge new_edge;
      size_type node1_uid = a.index();
      size_type node2_uid = b.index();

      new_edge.node1_uid = node1_uid;
      new_edge.node2_uid = node2_uid;
      new_edge.edge_uid = current_internal_edges;
      new_edge.value = value;
      new_edge.idx = current_registered_edges;

      // Add to list of edges.
      internal_edges_.push_back(new_edge);

      // Add to edge uid to index vector.
      i2u_edges_.push_back(current_internal_edges);

      // Update adjacency list.
      adj_[node1_uid].push_back(node2_uid);
      adj_[node2_uid].push_back(node1_uid);

      // Return an Edge that points to newly added edge.
      return Edge(this, current_internal_edges);
    }
    return Edge(this, has_edge_uid);
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
    i2u_nodes_.clear();
    i2u_edges_.clear();
  }

  //
  // Node and Edge Removal
  // TODO: Fix implementation of remove_node() and remove_edge().
  //
 
  /** Remove a node (and all its incident edges) from the graph.
   * @param[in] n_it    NodeIterator corresponding to node to remove.
   * @return the next valid node iterator preceeding removed node.
   *
   * @post The idx value for internal_nodes beyond *@n_it decreases by 1.
   * and the idx value for the internal node corresponding to *@n_it
   * is i2u_nodes_.size() to indicate it is not "registered" in the index. 
   * @post All edges incident to @n are removed from adjacency list
   * and idx values of internal_nodes are updated according
   * (see remove_edge()).
   */ 
  node_iterator remove_node(node_iterator n_it) {
    Node node = *n_it;
    size_type result  = remove_node(node);
    return n_it++;
  }

  /** Remove a node (and all its incident edges) from the graph.
   * @param[in] n    Node to remove.
   * @return the idx of the next valid node preceeding removed node.
   *
   * @post The idx value for internal_nodes beyond @n decreases by 1.
   * and the idx value for the internal node corresponding to @n
   * is i2u_nodes_.size() to indicate it is not "registered" in the index. 
   * @post All edges incident to @n are removed from adjacency list
   * and idx values of internal_nodes are updated according
   * (see remove_edge()).
   */ 
  size_type remove_node(const Node& n) {

    // For each node v in n's adj list, remove the (n, v)/(v, n) edge. 
    std::vector<size_type> adj_nodes = adj_[n.uid()];
    for (auto it = adj_nodes.begin(); it != adj_nodes.end(); ++it) {
      Node adj_node = node(n.index());

      // TODO: Fix this. Investigate bug which may cause this call to
      // remove_edge() to trigger a seg fault.
      //int result = remove_edge(n, adj_node);  
    }

    // Remove node from node index to uid vector.
    size_type index = n.index();
    i2u_nodes_.erase(i2u_nodes_.begin() + index);
    
    // Decrement the idx values of internal node objects following n.
    for (auto it = internal_nodes_.begin() + index;
        it != internal_nodes_.end(); ++it) {
      if (it == internal_nodes_.begin() + index) {
        internal_node& int_node = *it;
        int_node.idx = i2u_nodes_.size();
      } else {
        internal_node& int_node = *it;
        int_node.idx = int_node.idx - 1;
      }
    }

    return index + 1;
  }
  
  /** Remove an edge from the graph.
   * @param[in] e_it    Iterator corresponding to edge to remove.
   * @return the next valid EdgeIterator preceeding that of the removed node..
   *
   * @post The idx value for internal_edges beyond *@e_it decreases by 1.
   * and the idx value for the internal edges corresponding to *@e_it
   * is i2u_edges_.size() to indicate it is not "registered" in the index. 
   */ 
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge edge = *e_it;
    remove_edge(edge);
    return e_it++;
  }

  /** Remove an edge from the graph.
   * @param[in] a    Node belonging to edge to remove.
   * @param[in] b    Other node belonging to edge to remove.

   * @return the idx of the next valid edge preceeding that of the removed edge.
   *
   * @post The idx value for internal_edges beyond edge with @a and @b as
   * endpoints decreases by 1 and the idx value for the internal edges
   * corresponding to the edge with endpoints @a and @b is i2u_edges_.size()
   * to indicate it is not "registered" in the index. 
   */ 
  size_type remove_edge(const Node& a, const Node& b) {
    size_type edge_uid = has_edge(a, b);
    size_type result = remove_edge(Edge(this, edge_uid));
    return result;
  }
  
  /** Remove an edge from the graph.
   * @param[in] e    Edge to remove.

   * @return the idx of the next valid edge preceeding that of the removed edge.
   *
   * @post The idx value for internal_edges beyond @e decreases by 1.
   * and the idx value for the internal edges corresponding to @e
   * is i2u_edges_.size() to indicate it is not "registered" in the index. 
   */ 
  size_type remove_edge(const Edge& e) {

    // Remove edge from edge index to uid vector.
    size_type index = e.index();
    i2u_edges_.erase(i2u_edges_.begin() + index);

    // Decrement the idx values of internal edge objects following e.
    for (auto it = internal_edges_.begin() + index;
            it != internal_edges_.end(); ++it) {
      if (it == internal_edges_.begin() + index) {
        internal_edge& int_edge = *it;
        int_edge.idx = i2u_edges_.size();
      } else {
        internal_edge& int_edge = *it;
        int_edge.idx = int_edge.idx - 1;
      }
    }     
    
    // Remove this edge as it's represented in adjacency list.
    size_type n1_uid = e.node1().uid();
    size_type n2_uid = e.node2().uid();
    std::vector<size_type>& n1_adj = adj_[n1_uid];
    std::vector<size_type>& n2_adj = adj_[n2_uid];

    auto position1 = std::find(n1_adj.begin(), n1_adj.end(), n2_uid);
    if (position1 != n1_adj.end()) {
       n1_adj.erase(position1);
    }
    
    auto position2 = std::find(n2_adj.begin(), n2_adj.end(), n1_uid);
    if (position2 != n2_adj.end()) {
      n2_adj.erase(position2);
    }

    return index + 1; 
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
      size_type current_uid = graph_->i2u_nodes_[current_];
      return Node(graph_, current_uid);
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
    return NodeIterator(i2u_nodes_.size(), this_graph);
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
    IncidentIterator(size_type i, Node n, Graph* g) :
      current_{i}, spawn_node_{n}, graph_{g}, adj_nodes_{g->adj_[n.uid()]} {}

    // OPERATORS.

    Edge operator*() const {
      size_type adj_node_uid = adj_nodes_[current_];
      Node adj_node = Node(graph_, adj_node_uid);
      int has_edge_uid = graph_-> has_edge(spawn_node_, adj_node);
      Edge edge = Edge(graph_, has_edge_uid);

      // Make sure spawn node always returned by node1(). 
      if (edge.node2() == spawn_node_) {
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
	&& (spawn_node_ == rhs.spawn_node_);
    }
    
    bool operator!=(const IncidentIterator& rhs) const {
      return (graph_ != rhs.graph_) || (current_ != rhs.current_)
	|| (spawn_node_ != rhs.spawn_node_);
    }

   private:
    friend class Graph;
      size_type current_;
      Node spawn_node_;
      Graph* graph_;
      std::vector<size_type> adj_nodes_;
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
      size_type current_uid = graph_->i2u_edges_[current_];
      return Edge(graph_, current_uid);
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
    return EdgeIterator(i2u_edges_.size(), this_graph);
  }

 private:

  // Graph class's internals helper functions, data members, and so forth.

  struct internal_node {
    // The position of the node.
    Point position;
    // The value of the node.
    node_value_type value;
    // The uid of the node (based on order in which nodes were added).
    size_type uid;
    // The index of the node.
    size_type idx;
  };

  struct internal_edge {
    // The uids of the pair of nodes connected by this edge.
    size_type node1_uid;
    size_type node2_uid;
    // The uid of the edge (based on order in which edges were added).
    size_type edge_uid;
    // The value of the edge.
    edge_value_type value;
    // The index of the edge.
    size_type idx;
  };

  std::vector<internal_node> internal_nodes_;  // Container for nodes in graph.
  std::vector<internal_edge> internal_edges_;  // Container for edges in graph.

  // Container for node uid and corresponding list of uids of incident nodes.
  std::vector<std::vector<size_type>> adj_;

  std::vector<size_type> i2u_nodes_;  // Container for *registered* nodes.
  std::vector<size_type> i2u_edges_;  // Container for *registered* edges.
};

#endif // CME212_GRAPH_HPP
