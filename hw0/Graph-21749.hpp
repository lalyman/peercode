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
class Graph {
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
  class Node {
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
      assert(index_ < graph_->size());
      return graph_->nodes_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
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

    /** Constructor used by the graph class to construct valid edges. */
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    size_type curr_size = size();
    nodes_.push_back(position);
    edges_.push_back(std::vector<size_type>{});
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
  class Edge {
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
    size_type result = 0;
    for (size_type i = 0; i < size(); ++ i)
      result += edges_[i].size();
    return result;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    size_type curr_index = 0;
    for (size_type j = 0; j < size(); ++ j){
      if (i >= curr_index + edges_[j].size())
        curr_index += edges_[j].size();
      else {
        for (size_type k = 0; k < edges_[j].size(); ++ k){
          if (curr_index == i)
            return Edge(this, j, edges_[j][k]);
          curr_index += 1;
        }
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
    assert(has_node(a));
    assert(has_node(b));
    size_type min_index = std::min(a.index(),b.index());
    size_type max_index = std::max(a.index(),b.index());
    for (size_type i = 0; i < edges_[min_index].size(); ++ i){
      if (edges_[min_index][i] == max_index)
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
    assert(has_node(a));
    assert(has_node(b));
    size_type min_index = std::min(a.index(),b.index());
    size_type max_index = std::max(a.index(),b.index());
    if (!has_edge(a,b)){
      edges_[min_index].push_back(max_index);
    }
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_ = std::vector<Point>();
    edges_ = std::vector<std::vector<size_type>>();
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal representation of nodes
  std::vector<Point> nodes_;

  // Internal representation of edges.
  // This is an adjacency list representation
  // edges[i][j] = k means that there is an edge between nodes_[i] & nodes_[k]
  std::vector<std::vector<size_type>> edges_;

};

#endif // CME212_GRAPH_HPP
