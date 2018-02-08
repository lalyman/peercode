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

  struct internal_node;
  struct internal_edge;

 //private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Predeclaration of iterator types. */
  class node_iterator;
  class edge_iterator;
  class incident_iterator;

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    /** Return this node's degree. */
    size_type degree() const {
      return fetch().incident_edge_index.size();
    }

    /** Return a reference to this node's value*/
    node_value_type& value() {
      // HW1: YOUR CODE HERE
      return fetch().value;
    }

    /** Return this node's value*/
    const node_value_type& value() const{
      //: HW1 YOUR CODE HERE
      return fetch().value;
    }

    /** Return this node's graph pointer*/
    Graph* graph() const {
      // HW0: YOUR CODE HERE
      return graph_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (index() == n.index()) && (graph() == n.graph());
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
      // HW0: YOUR CODE HERE
      return (index() < n.index());
    }

    incident_iterator edge_begin() const {
      return incident_iterator(this, 0);
    }

    incident_iterator edge_end() const {
      return incident_iterator(this, degree());
    }

   private:
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;     // Pointer back to the parent Graph
    size_type index_;  // this node's index, a number in the range [0, graph_size)
    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
    /** Helper method to return the appropriate element.*/
    internal_node& fetch() const {
      return graph_->nodes_[index_];
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //friend class incident_iterator;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
    // HW0: YOUR CODE HERE
    size_type new_index {num_nodes()};
    nodes_.push_back({position, value, new_index});
    return Node(this, new_index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (Node(this, n.index()) == n)
      return true;

    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return node1_; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1()) && (node2() == e.node2())) ||
             ((node1() == e.node2()) && (node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (node1() < e.node1()) || (node2() < e.node2());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Edge(const Node& n1, const Node& n2)
        : node1_(n1), node2_(n2) {
    }
    Edge(const Graph* graph, size_type n1_index, size_type n2_index)
        : node1_(Node(graph, n1_index)), node2_(Node(graph, n2_index)) {
    }

    Node node1_;
    Node node2_;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edges_[i].node1_index, edges_[i].node2_index);
  }

  Edge edgerev(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edges_[i].node2_index, edges_[i].node1_index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type i = 0; i < num_edges(); i++) {
      Edge test_edge(this, edges_[i].node1_index, edges_[i].node2_index);
      if (((test_edge.node1() == a) && (test_edge.node2() == b)) || 
          ((test_edge.node1() == b) && (test_edge.node2() == a)))
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
    // HW0: YOUR CODE HERE
    edges_.push_back({a.index(), b.index()});
    a.fetch().incident_edge_index.push_back(num_edges()-1);
    b.fetch().incident_edge_index.push_back(num_edges()-1);
    return Edge(a, b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

  /** @class Graph::node_iterator
   * @brief Class for iterating over the graph's nodes.
  */
  class node_iterator : private equality_comparable<node_iterator> {

   public:
    node_iterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

    Node operator*() const {
      return graph_->node(index_);
    }

    node_iterator& operator++() {
      index_++;
      return *this;
    }

    bool operator==(const node_iterator& iter) const {
      return (iter.graph_ == graph_) && (iter.index_ == index_);
    }

   private:

    Graph* graph_;     // Pointer back to the parent Graph
    size_type index_;  // this iterator's current node index, a number in the range [0, graph_size)

  };

  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  node_iterator node_end() const {
    return node_iterator(this, size());
  }

  /** @class Graph::incident_iterator
   * @brief Class for iterating over the graph's nodes.
  */
  class incident_iterator : private equality_comparable<incident_iterator> {
   public:
    incident_iterator(const Node* node, size_type index)
        : node_(const_cast<Node*>(node)), index_(index) {
    }

    Edge operator*() const {
      size_type edge_index {node_->fetch().incident_edge_index[index_]};
      if (node_->graph()->edge(edge_index).node1() == *node_)
        return node_->graph()->edge(edge_index);
      else
        return node_->graph()->edgerev(edge_index);
    }

    incident_iterator& operator++() {
      index_++;
      return *this;
    }

    bool operator==(const incident_iterator& iit) const {
      return (*(iit.node_) == *node_) && (iit.index_ == index_);
    }

   private:

    Node* node_;     // Pointer back to the parent Graph
    size_type index_;  // this iterator's current node index, a number in the range [0, graph_size)

  };

  /** @class Graph::edge_iterator
   * @brief Class for iterating over the graph's edges.
  */
  class edge_iterator : private equality_comparable<edge_iterator> {
   public:
    edge_iterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

    Edge operator*() const {
      return graph_->edge(index_);
    }

    edge_iterator& operator++() {
      index_++;
      return *this;
    }

    bool operator==(const edge_iterator& iter) const {
      return (iter.graph_ == graph_) && (iter.index_ == index_);
    }

   private:

    Graph* graph_;     // Pointer back to the parent Graph
    size_type index_;  // this iterator's current node index, a number in the range [0, graph_size)

  };

  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }


 private:

  // HW0: YOUR CODE HERE
  // Internal type for graph nodes
  struct internal_node {
    Point position;         // The position of a node
    node_value_type value;  // The value of a node
    size_type index;        // The id of a node
    std::vector<size_type> incident_edge_index;   // List of incident edge indices
  };
  // Internal type for graph edges
  struct internal_edge {
    size_type node1_index;
    size_type node2_index;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
