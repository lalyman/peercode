#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <functional>
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

  // Internal type for nodes
  struct internal_node;
  
  // Internal type for node adjacency objects
  struct adj_obj;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value. */
  using node_value_type = V;

  /** Type of edge value. */
  using edge_value_type = E;

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
      return graph_->nodes_vec[graph_->i2u[uid]].position;
    }

    /** Setter function - sets the position of the node */
    Point& position() {
      return graph_->nodes_vec[graph_->i2u[uid]].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid;
    }


    /** Return node value */

    /** Setter function - sets the value of the node */
    node_value_type& value() {
      return graph_->nodes_vec[graph_->i2u[uid]].val;
    }

    /** Getter function - retrieves the value of the node */
    const node_value_type& value() const {
      return graph_->nodes_vec[graph_->i2u[uid]].val;
    }

    /** Returns the number of incident edges for a given node */
    size_type degree() const {
      return graph_->nodes_vec[graph_->i2u[uid]].adj.size();
    }

    /** Returns the starting iterator needed to iterate over the
     *  the incident edges of a given node
     */
    incident_iterator edge_begin() const {
      IncidentIterator begin(this, graph_, 0);
      return begin;
    }

    /** Returns the ending iterator indicating when to stop iterating
     *  the incident edges of a given node
     */
    incident_iterator edge_end() const {
      IncidentIterator end(this, graph_, degree());
      return end;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) && (uid == n.uid));
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
      if (*this == n) {
        return false;
      }
      else if (uid == n.uid) {          // Same index, different graph
        std::less<Graph*> less;
        return less(graph_,n.graph_);    // std::less gives total ordering
      }
      else {
        return (uid < n.uid);           // Else, return based on index
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // Unique node id
    size_type uid;

    /** Private Node constructor (gives access to outer graph and 
     *  initializes node's id)
     */
    Node(const Graph* graph, size_type i)
        : graph_(const_cast<Graph*>(graph)), uid(i) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& val = node_value_type()) {
    size_type node_id = nodes_vec.size();
    i2u.push_back(node_id);
    internal_node n_internal(i2u.size()-1, position, val);
    nodes_vec.push_back(n_internal);

    return node(i2u.size() - 1); // Use node function to return most recent node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= i2u.size()) {
      return false;
    }
    return (node(n.index()) == n); // Equality operator for nodes checks
                                   // index equality and graph pointer equality
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   *
   * If node index argument is greater than number of nodes in graph 
   * (minus 1), then an invalid node is returned
   */
  Node node(size_type i) const {
    Node n(this,i);
    return n;
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
      return graph_->node(n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2);
    }

    /** Setter method for edge value */
    edge_value_type& value() {
      
      // Find location of n2 in adjacency vector
      size_type i;
      for (i = 0; i < graph_->nodes_vec[graph_->i2u[n1]].adj.size(); i++) {
        if (graph_->i2u[n2] == graph_->nodes_vec[graph_->i2u[n1]].adj[i].id) {
          break;
        }
      }
      return graph_->nodes_vec[graph_->i2u[n1]].adj[i].val;
    }

    /** Getter method for edge value */
    const edge_value_type& value() const {

      // Find location of n2 in adjacency vector
      size_type i;
      for (i = 0; i < graph_->nodes_vec[graph_->i2u[n1]].adj.size(); i++) {
        if (graph_->i2u[n2] == graph_->nodes_vec[graph_->i2u[n1]].adj[i].id) {
          break;
        }
      }
      return graph_->nodes_vec[graph_->i2u[n1]].adj[i].val;
    }

    /** Return the Euclidean distance between this Edge's two nodes */
    double length() const {
      Point n1_pos = node1().position();
      Point n2_pos = node2().position();
      return (norm(n1_pos - n2_pos));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((node1() == e.node1() && node2() == e.node2()) ||
          (node1() == e.node2() && node2() == e.node1())) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * The hierarchy of checks is:
     * 1. Check if edges are equal based on == operator
     * 2. Check if edges connect same node indices but are from different graphs
     * 3. Check if smallest node indices of the two edges are equal
     * 4. Return true if edge 1's smallest node index is less than edge 2's
     */
    bool operator<(const Edge& e) const {
      if (*this == e) {
        return false;
      }
      else if (smaller(n1,n2) == smaller(e.n1,e.n2) &&
               larger(n1,n2) == larger(e.n1,e.n2)) {  // Connects same node indices, 
        std::less<Graph*> less;                      // but from different graph
        return less(graph_,e.graph_);    // std::less gives total ordering
      }
      else if (smaller(n1, n2) == smaller(e.n1, e.n2)) {
        return (larger(n1, n2) < larger(e.n1, e.n2));
      }
      return (smaller(n1,n2) < smaller(e.n1, e.n2));
    }
      


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the outer graph
    Graph* graph_;

    // Index of edge's node 1
    size_type n1;

    // Index of edge's node 2
    size_type n2;

    /** Helper function that returns the smaller of two distinct
     *  size_type elements.
     */
    size_type smaller(size_type a, size_type b) const {
      if (a < b) {
        return a;
      }
      return b;
    }

    /** Helper function that returns the larger of two distinct
     *  size_type elements.
     */
    size_type larger(size_type a, size_type b) const {
      if (a < b) {
        return b;
      }
      return a;
    }

    /** Private Edge constructor that gives access to outer graph and
     * initializes the indices of the nodes it connects
     */
    Edge(const Graph* g, size_type ind_1, size_type ind_2)
          : graph_(const_cast<Graph*>(g)), n1(ind_1), n2(ind_2) {

    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type result = 0;
    for (size_type i = 0; i < i2u.size(); i++) {
      result += node(i).degree();
    }
    return result/2; // Correct for double counting of edges
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   *
   * If edge index argument is greater than number of edges in graph 
   * (minus 1), then an invalid edge is returned
   */

  Edge edge(size_type i) const {

    if (i > num_edges()-1) {
      Edge e;
      return e; // Edge index argument greater than num_edges, return invalid edge
    }
    edge_iterator iter = edge_begin();
    // Advance beginning edge iterator i times
    for (size_type k = 0; k < i; ++k) {
      ++iter;
    }
    return *iter;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
      if ((*it).node2() == b) {
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
    if (!has_edge(a,b)) {
      nodes_vec[i2u[a.index()]].adj.push_back({i2u[b.index()],edge_value_type()});
      nodes_vec[i2u[b.index()]].adj.push_back({i2u[a.index()],edge_value_type()});
    }
    Edge e(this, a.index(), b.index());
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_vec.clear();
    i2u.clear();
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


    /** Dereference operator - returns the node that the given node
     *  iterator refers to
     */
    Node operator*() const {
      return graph_->node(id);
    }
    
    /** Increment operator - returns a NodeIterator that has been advanced
     *  one position
     */
    NodeIterator& operator++() {
      id++;
      return *this;
    }

    /** Equality operator - returns true if NodeIterator objects
     *  belong to same graph and refer to the same node
     */
    bool operator==(const NodeIterator& iter) const {
      return ((this->graph_ == iter.graph_) && (this->id == iter.id));
    }

   private:
    friend class Graph;

    /** Pointer to outer graph */
    Graph* graph_;

    /** Index of node that iterator is currently on */
    size_type id;
    
    /** Private NodeIterator constructor (gives access to outer graph and 
     *  keeps track of current index)
     */
    NodeIterator(const Graph* g, size_type i)
        : graph_(const_cast<Graph*>(g)), id(i) {
    }
  };


  /** Returns the starting iterator needed to iterate over the
   *  the nodes of a graph
   */
  node_iterator node_begin() const {
    NodeIterator begin(this, 0);
    return begin;
  }
  
  /** Returns the ending iterator indicating when to stop iterating over
   *  the nodes of a graph
   */
  node_iterator node_end() const {
    NodeIterator end(this, num_nodes());
    return end;
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

    /** Dereference operator - returns the edge that the given IncidentIterator
     *  refers to
     */
    Edge operator*() const {
      return Edge(graph_, node_->index(), 
                  graph_->nodes_vec[graph_->nodes_vec[graph_->i2u[node_->index()]].adj[id].id].index);
    }
    
    /** Increment operator - returns an IncidentIterator that has been advanced
     *  by one position
     */
    IncidentIterator& operator++() {
      id++;
      return *this;
    }
    
    /** Equality operator - returns true if the two IncidentIterators are
     *  generated by the same node and refer to the same edge
     *
     *  Input is a second incident iterator
     */
    bool operator==(const IncidentIterator& inc_iter) const {
      return ((node_ == inc_iter.node_) && (id == inc_iter.id));
    }

   private:
    friend class Graph;

    /** Pointer to node that generated the iterator. */
    Node* node_;

    /** Pointer to outer graph. */
    Graph* graph_;

    /** Index of edge that iterator is currently on */
    size_type id;
    
    /** Private IncidentIterator constructor (gives access to outer graph and 
     *  node, and keeps track of current index)
     */
    IncidentIterator(const Node* n, const Graph* g, size_type i)
        : node_(const_cast<Node*>(n)), graph_(const_cast<Graph*>(g)), id(i) {
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
    
    /** Dereference operator - returns the edge that the given EdgeIterator
     *  refers to
     */
    Edge operator*() const {
      Edge e(graph_, node_id, 
             graph_->nodes_vec[graph_->nodes_vec[graph_->i2u[node_id]].adj[adj_id].id].index);
      return e;
    }
    
    /** Increment operator - returns an EdgeIterator that has been advanced
     *  by one position.
     *
     *  Skips duplicate edges - always returns and EdgeIterator that refers
     *  to the next unique edge, except when it reaches the end iterator
     */
    EdgeIterator& operator++() {
      internal_node current_node = graph_->nodes_vec[graph_->i2u[node_id]];

      while (node_id != graph_->num_nodes()) {
        // Iterate through incident edges of current node, if it has any
        if ((int)adj_id < (int)current_node.adj.size()-1) {
          ++adj_id;
          // Use min indexing to determine whether or not to skip a duplicate edge.
          if (node_id < 
              graph_->nodes_vec[current_node.adj[adj_id].id].index) {
            return *this;
          }
        }

        // No more incident edges for current node, move to next node
        else {
          ++node_id;
          adj_id = 0;  // Reset index for adjacent nodes
          if (node_id == graph_->num_nodes()) {
            break;
          }
          current_node = graph_->nodes_vec[graph_->i2u[node_id]];
          if (current_node.adj.size() != 0) {
            if (node_id < 
                graph_->nodes_vec[current_node.adj[adj_id].id].index) {
                return *this;
            }
          }
        }
      }
      return *this; // At end iterator
    }

    /** Equality operator - returns true if the two EdgeIterators belong to
     *  the same graph and refer to the same edge
     *
     *  Input is a second edge iterator
     */

    bool operator==(const EdgeIterator& iter_2) const {
      return (graph_ == iter_2.graph_ &&
              node_id == iter_2.node_id &&
              adj_id == iter_2.adj_id);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    // Pointer to outer graph
    Graph* graph_;

    // Internal node index
    size_type node_id;

    // Index for node adjacent to given internal node
    size_type adj_id;

    /** Private EdgeIterator constructor (gives access to outer graph and 
     *  keeps track of current index)
     *
     *  Increments internal data to ensure the constructed iterator begins
     *  at a valid edge
     */
    EdgeIterator(const Graph* g, size_type n_id)
        : graph_(const_cast<Graph*>(g)), node_id(n_id) {

      adj_id = 0;
      if (node_id == graph_->num_nodes()) {
        // Do nothing, at end iterator
      }
      else {  // Advance iterator to first valid edge
        internal_node current_node = graph_->nodes_vec[graph_->i2u[node_id]];
        while (current_node.adj.size() == 0) {
          ++node_id;
          if (node_id == graph_->num_nodes()) {
            break;
          }
          current_node = graph_->nodes_vec[graph_->i2u[node_id]];
        }
      }
    }
  };

  
  /** Returns the starting iterator needed to iterate over the
   *  the edges of a graph
   */
  edge_iterator edge_begin() const {
    EdgeIterator begin(this, 0);
    return begin;
  }

  /** Returns the ending iterator indicating when to stop iterating over
   *  the edges of a graph
   */
  edge_iterator edge_end() const {
    EdgeIterator end(this, i2u.size());
    return end;
  }

  /** REMOVAL METHODS
   *********************
   */

  /** Removes the input node, and returns the index of the node that
   *  was removed (now refers to a different node)
   *
   *  Input is a node object; function will check if graph has node
   *  and remove it and any edges it is connected to if it is.
   *  Otherwise, nothing will be done.
   *
   *  Returns the input nodes's index
   */
  size_type remove_node(const Node& n) {
    // Delete node from adjacent nodes' adjacency vectors using incident
    // iterator
    if (has_node(n)) {
      for (auto nit = n.edge_begin(); nit != n.edge_end(); ++nit) {
      Node n2 = (*nit).node2();
      remove_edge_helper(n2,n);
      }
      // Delete node itself using copy-pop_back strategy
      i2u[n.index()] = i2u[i2u.size()-1]; // Overwrite id we wish to remove
      i2u.pop_back();
      nodes_vec[i2u[n.index()]].index = n.index(); // Set the internal node id to
                                                 // node index (location in i2u)
    }
    return n.index();
  }

  /** Removes the node referred to by the input node iterator
   *
   *  Input is a node iterator; function will check if graph has node
   *  referred to by the node iterator and remove it and any edges 
   *  it is connected to if it is. Otherwise, nothing will be done.
   *
   *  Returns a Node Iterator pointing to the new node that has
   *  replaced the position of the removed node
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type new_ind = remove_node(*n_it);
    return NodeIterator(this, new_ind);
  }

  /** Removes the edge specified by the two input nodess
   *
   *  Input are two nodes; function will check if graph has an
   *  edge whose endpoints are the two nodes, and delete
   *  the edge if it does. Otherwise, nothing is done.
   *
   *  Returns a 1 if an edge was removed and 0 otherwise
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1,n2)) {
      remove_edge_helper(n1,n2);
      remove_edge_helper(n2,n1);
      return 1;
    }
    return 0;
  }

  /** Removes the input edge
   *
   *  Input is an edge; function will check if graph has the
   *  given edge and delete the edge if it does. Otherwise, 
   *  nothing is done.
   *
   *  Returns 0 after completion
   */
  size_type remove_edge(const Edge& edge) {
    remove_edge(edge.node1(), edge.node2());
    return 0; 
  }

  /** Removes the edge referred to by the given input edge iterator
   *
   *  Input is an edge iterator; function will check if graph
   *  has the edge specified by the input edge iterator and delete
   *  the edge if it does. Otherwise, nothing is done..
   *
   *  Returns an edge iterator object that points to the beginning
   *  of the adjacent edge vector for node_1 of the dereferenced
   *  input edge iterator
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return EdgeIterator(this, *e_it.node1().index());
  }

 private:
  
  /** Helper Functions 
   ***********************
   */

  /** Remove Edge helper - removes n2 from n1's adjacency vector using a
   *  copy-pop_back strategy
   *
   *  Inputs are two nodes, n1 and n2. Removal is directed one-way;
   *  n1 will not be removed from n2's adjacency vector
   *
   *  No return type (void)
   */
  void remove_edge_helper(const Node& n1, const Node& n2) {
    size_type i;
      for (i = 0; i < n1.degree(); i++) {
        if (i2u[n2.index()] == nodes_vec[i2u[n1.index()]].adj[i].id) {
          break;
        }
      }
      // Overwrite adjacent object associated with n2 with last object
      // in adjacency vector
      nodes_vec[i2u[n1.index()]].adj[i] = nodes_vec[i2u[n1.index()]].adj[n1.degree()-1];
      nodes_vec[i2u[n1.index()]].adj.pop_back();
  }

  // Adjacency object - contains the id of the adjacent node and a value
  // associated with the edge connecting those nodes
  struct adj_obj {
    size_type id;
    edge_value_type val;
  };

  // Internal node type
  struct internal_node {
    size_type index; // Index in i2u vector
    Point position;
    node_value_type val;
    std::vector<adj_obj> adj; // Nodes that the internal node is connected to
    
    // Constructor for internal_node when passed a point object
    internal_node(size_type id, const Point& p, const node_value_type& v) 
        : index(id), position(p), val(v) {
    }
  };


  // Internal data for graph's nodes including an adjacency list, stored in vectors
  std::vector<internal_node> nodes_vec;

  /** Vector that contains a node mapping for the internal nodes
   *
   *  Invariant: For all 0 <= i < num_nodes:
   *             nodes_vec[i2u[i]].index = i
   *
   *  Note: Using i2u allows for addition and removal of nodes to be
   *  handled entirely with a copy-pop_back strategy on a vector of
   *  size_types rather than a vector of expensive objects
   */
  std::vector<size_type> i2u;
};

#endif // CME212_GRAPH_HPP
