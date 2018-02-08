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
template <typename V>
class Graph {
 private:

  // Internal type for nodes
  struct internal_node;
  
  // Internal type for edges
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value. */
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
      return graph_->nodes_vec[uid].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid;
    }

    // HW1: YOUR CODE HERE

    /** Return node value */

    /** Setter function - sets the value of the node */
    node_value_type& value() {
      return graph_->nodes_vec[uid].val;
    }

    /** Getter function - retrieves the value of the node */
    const node_value_type& value() const {
      return graph_->nodes_vec[uid].val;
    }

    /** Returns the number of incident edges for a given node */
    size_type degree() const {
      return graph_->nodes_vec[uid].adj.size();
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
    return nodes_vec.size();
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
    internal_node n_internal(position, val);
    nodes_vec.push_back(n_internal);

    return node(nodes_vec.size() - 1); // Use node function to return most recent node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= nodes_vec.size()) {
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
    size_type smaller(size_type& a, size_type& b) const {
      if (a < b) {
        return a;
      }
      return b;
    }

    /** Helper function that returns the larger of two distinct
     *  size_type elements.
     */
    size_type larger(size_type& a, size_type& b) const {
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
    for (size_type i = 0; i < nodes_vec.size(); i++) {
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
    for (size_type i = 0; i < a.degree(); ++i) {
      if (nodes_vec[a.index()].adj[i] == b.index()) {
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
      nodes_vec[a.index()].adj.push_back(b.index());
      nodes_vec[b.index()].adj.push_back(a.index());
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

    /** Equality operator - returns true if NodeIterator objects have
     *  belong to same graph and refer to the same node
     */
    bool operator==(const NodeIterator& iter) const {
      return ((this->graph_ == iter.graph_) && (this->id == iter.id));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

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

  // HW1 #2: YOUR CODE HERE

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
    NodeIterator end(this, nodes_vec.size());
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

    // HW1 #3: YOUR CODE HERE

    /** Dereference operator - returns the edge that the given IncidentIterator
     *  refers to
     */
    Edge operator*() const {
      return Edge(graph_, node_->index(), 
                  graph_->nodes_vec[node_->index()].adj[id]);
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
     */
    bool operator==(const IncidentIterator& inc_iter) const {
      return ((node_ == inc_iter.node_) && (id == inc_iter.id));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

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

    // HW1 #5: YOUR CODE HERE
    
    /** Dereference operator - returns the edge that the given EdgeIterator
     *  refers to
     */
    Edge operator*() const {
      Edge e(graph_, node_id, graph_->nodes_vec[node_id].adj[adj_id]);
      return e;
    }
    
    /** Increment operator - returns an EdgeIterator that has been advanced
     *  by one position.
     *
     *  Skips duplicate edges - always returns and EdgeIterator that refers
     *  to the next unique edge, except when it reaches the end iterator
     */
    EdgeIterator& operator++() {
      internal_node current_node = graph_->nodes_vec[node_id];

      while (node_id != graph_->nodes_vec.size()) {
        // Iterate through incident edges of current node, if it has any
        if (adj_id < current_node.adj.size()-1) {
          ++adj_id;
          // Use min indexing to determine whether or not to skip a duplicate edge.
          if (node_id < current_node.adj[adj_id]) {
            return *this;
          }
        }

        // No more incident edges for current node, move to next node
        else {
          ++node_id;
          adj_id = 0;  // Reset index for adjacent nodes
          if (node_id == graph_->nodes_vec.size()) {
            break;
          }
          current_node = graph_->nodes_vec[node_id];
          if (node_id < current_node.adj[adj_id]) {
            return *this;
          }
        }
      }
      return *this; // At end iterator
    }

    /** Equality operator - returns true if the two EdgeIterators belong to
     *  the same graph and refer to the same edge
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
      if (node_id == graph_->nodes_vec.size()) {
        // Do nothing, at end iterator
      }
      else {  // Advance iterator to first valid edge
        internal_node current_node = graph_->nodes_vec[node_id];
        while (current_node.adj.size() == 0) {
          ++node_id;
          if (node_id == graph_->nodes_vec.size()) {
            break;
          }
          current_node = graph_->nodes_vec[node_id];
        }
      }
    }
  };

  // HW1 #5: YOUR CODE HERE
  
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
    EdgeIterator end(this, nodes_vec.size());
    return end;
  }

 private:
  
  // Internal node type
  struct internal_node {
    Point position;
    node_value_type val;
    std::vector<size_type> adj; //Nodes that the internal node is connected to
    
    // Constructor for internal_node when passed a point object
    internal_node(const Point& p, const node_value_type& v) 
        : position(p), val(v) {
    }
  };


  // Internal data for graph's nodes including an adjacency list, stored in vectors
  std::vector<internal_node> nodes_vec;
};

#endif // CME212_GRAPH_HPP
