#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
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
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of the node values. */
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /**
    * @brief Return a reference to this node's value
    *
    * @return a reference to an object of type node_value_type
    */
    node_value_type& value() {
      return graph_->nodes_[index_].value;
    }

    /**
    * @brief Return a reference to this node's value
    *
    * @return a reference to a constant object of type node_value_type
    * @post this will not be modified
    */
    const node_value_type& value() const {
      return const_cast<node_value_type& >(graph_->nodes_[index_].value);
    }

    /**
    * @brief Return this node's degree (number of incident edges)
    *
    * @return the number of edges incident to the current node in its graph
    * @post this will not be modified
    */
    size_type degree() const {
      return graph_->edges_[index_].size();
    }

    /**
    * @brief Return an iterator for the position belonging to the first edge incident to this node
    *
    * @return an iterator of type incident_iterator
    * @post this will not be modified
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }

    /**
    * @brief Return an iterator for the position one past the final edge incident to this node
    *
    * @return an iterator of type incident_iterator
    * @post this will not be modified
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph_ == graph_ && n.index() == index();
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
      if (n == *this)
        return false;
      if (n.graph_ == graph_)
        return (index_ < n.index());
      return (this < &n);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // pointer back to the graph container
    Graph* graph_;
    // This node's index
    size_type index_;

    /** Constructor used by the graph class to construct valid nodes. */
    Node(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {};

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
   * @post result_node.value() == value
   * @post result_node.position() == position
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type curr_size = size();
    internal_node n = {position, value};
    nodes_.push_back(n);
    edges_.push_back(std::vector<internal_edge>{});
    return Node(this, curr_size);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= size())
      return false;
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i >= size()){
      return Node();
    };
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_index_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_index_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (node1() == e.node1() && node2() == e.node2())
        return true;
      if (node1() == e.node2() && node2() == e.node1())
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (e == *this)
        return false;
      if (e.graph_ == graph_ && e.node1_index_ != node1_index_)
        return (node1_index_ < e.node1_index_);
      if (e.graph_ == graph_ && e.node2_index_ != node2_index_)
        return (node2_index_ < e.node2_index_);
      return (this < &e);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_index_;
    size_type node2_index_;

    /** Constructor used by the graph class to construct valid edges. */
    Edge(const Graph* graph, size_type node1_index, size_type node2_index)
    : graph_(const_cast<Graph*>(graph)), node1_index_(node1_index), node2_index_(node2_index) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i >= num_edges()){
      return Edge();
    };
    for (size_type j = 0; j < size(); ++ j){
      for (size_type k = 0; k < edges_[j].size(); ++ k){
        if (edges_[j][k].edge_index == i)
            return Edge(this, j, edges_[j][k].node2_index);
      }
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (!has_node(a) || !has_node(b)){
      return false;
    }
    for (size_type i = 0; i < edges_[a.index()].size(); ++ i){
      if (edges_[a.index()][i].node2_index == b.index())
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
    if (!has_node(a) || !has_node(b) || a == b){
      return Edge();
    }
    if (!has_edge(a,b)){
      edges_[a.index()].push_back(internal_edge{b.index(), num_edges_});
      edges_[b.index()].push_back(internal_edge{a.index(), num_edges_});
      num_edges_++;
    }
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_ = std::vector<internal_node>();
    edges_ = std::vector<std::vector<internal_edge>>();
    num_edges_ = 0;
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

    /** @brief Return the Node at the iterator's current position
    *
    * @return a Node
    * @post this will not be modified
    */
    Node operator*() const {
      return graph_->node(index_);
    }

    /** @brief Return a reference to a NodeIterator at one past the iterator's current position
    *
    * @return a reference to a NodeIterator
    */
    NodeIterator& operator++() {
      index_++;
      return *this;
    }

    /** @brief Test whether this node iterator and @a n are equal
    *
    * @return a boolean value; true only if both iterators belong to the same graph and have the same node position
    * @post this will not be modified
    */
    bool operator==(const NodeIterator& n) const {
      return n.graph_ == graph_ && n.index_ == index_;
    }

   private:
    friend class Graph;

    NodeIterator(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {};

    Graph* graph_;
    size_type index_;

  };

  /** @brief Return an interator for the position of the node at index 0
  *
  * @return an interator of type node_iterator
  * @post this will not be modified
  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** @brief Return an interator for one past the final node in this graph (index num_nodes())
  *
  * @return an interator of type node_iterator
  * @post this will not be modified
  */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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

    /** @brief Return the Edge at the iterator's current position
    *
    * @return an Edge
    * @post this will not be modified
    */
    Edge operator*() const {
      return Edge(graph_, node1_index_, graph_->edges_[node1_index_][node2_index_].node2_index);
    }

    /** @brief Return a reference to an IncidentIterator at one past the iterator's current position
    *
    * @return a reference to an IncidentIterator
    */
    IncidentIterator& operator++(){
      node2_index_++;
      return *this;
    }

    /** @brief Test whether this incident iterator and @a i are equal
    *
    * @return a boolean value; true only if both iterators belong to the same graph and have the same position
    * @post this will not be modified
    */
    bool operator==(const IncidentIterator& i) const {
      return i.graph_ == graph_ && i.node2_index_ == node2_index_ && i.node1_index_ == node1_index_;
    }

   private:
    friend class Graph;

    IncidentIterator(const Graph* graph, size_type node1_index, size_type node2_index)
    : graph_(const_cast<Graph*>(graph)), node1_index_(node1_index), node2_index_(node2_index) {};

    Graph* graph_;
    size_type node1_index_;
    size_type node2_index_;

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
    EdgeIterator() {
    }

    /** @brief Return the Edge at the iterator's current position
    *
    * @return an Edge
    * @post this will not be modified
    */
    Edge operator*() const {
      return *edge_;
    }

    /** @brief Return an EdgeIterator for the position one past this iterator's position
    *
    * @return a reference to an EdgeIterator
    * @post if result != edge_end(), (*result).node1().index() < (*result).node2().index()
    */
    EdgeIterator& operator++() {
      ++edge_;
      while (node_ != graph_-> node_end()) {
        while (edge_ != (*node_).edge_end()){
          // only return the current edge if node1 < node2, otherwise, keep searching.
          // This is to avoid returning duplicate edges
          if ((*edge_).node1().index() < (*edge_).node2().index()){
            return *this;
          }
          ++edge_;
        }
        // Increment the node iterator once we get to the end of the current node's incident edges
        // As long as the node iterator is still pointing to a valid node, reset the incident iterator and continue
        ++node_;
        if (node_ != graph_->node_end()){
          edge_ = (*node_).edge_begin();
        }
      }
      return *this;
    }

    /** @brief Test whether this EdgeIterator and @a ei are equal
    *
    * @return a boolean value; true only if both edge iterators belong to the same graph and have the same position
    * @post this will not be modified
    */
    bool operator==(const EdgeIterator& ei) const {
      return graph_ == ei.graph_ && edge_== ei.edge_ && node_ == ei.node_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    NodeIterator node_;
    IncidentIterator edge_;

    EdgeIterator(const Graph* graph, NodeIterator node, IncidentIterator edge )
    : graph_(const_cast<Graph*>(graph)), node_(node), edge_(edge) {};

  };

  /** @brief Return an iterator for the edge at position 0 in the graph
  *
  * @return an edge_iterator
  * @post this will not be modified
  */
  edge_iterator edge_begin() const{
    if(num_edges() == 0){
      return edge_end();
    }
    return EdgeIterator(this, node_begin(), node(0).edge_begin());
  }

  /** @brief Return an iterator for the position one past the last edge in the graph
  *
  * @return an edge_iterator
  * @post this will not be modified
  */
  edge_iterator edge_end() const{
    if(num_edges() == 0){
      return EdgeIterator(this, node_end(), IncidentIterator(this, 0, 0));
    }
    return EdgeIterator(this, node_end(), node(num_nodes()-1).edge_end());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
    Point position;
    node_value_type value;
  };

  struct internal_edge{
    size_type node2_index;
    size_type edge_index;
  };

  // Internal representation of nodes
  std::vector<internal_node> nodes_;

  // Internal representation of edges.
  // This is an adjacency list representation
  // edges[i][j] = k means that there is an edge between nodes_[i] & nodes_[k]
  std::vector<std::vector<internal_edge>> edges_;

  // Quick reference to the number of edges in the graph
  size_type num_edges_ = 0;

};

#endif // CME212_GRAPH_HPP
