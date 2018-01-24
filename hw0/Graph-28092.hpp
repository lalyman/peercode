#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
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

  /** Type of Node and Edge indexes and sizes. */
  using size_type = unsigned;

  //** Type of data structure maintaining edge data, a nested unordered_map */
  using edge_map_type = std::unordered_map<size_type,
      std::unordered_map<size_type, size_type>>;

  //
  // PRIVATE INTERNAL DATA
  //

 private:
  /** Internal structure for holding Node data */
  std::vector<Point> nodes_;

  /** Internal structure for holding Edge data */
  edge_map_type edges_;

  /** Maintain number of edges */
  size_type num_edges_;

 public:
  //
  // CONSTRUCTOR AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_(), edges_(), num_edges_(0)
  {
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
    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (&graph_ == &n.graph_ && uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any geometric meaning.
     *
     * The node ordering relation obeys trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return (&graph_ == &n.graph_ && uid_ < n.uid_) || (&graph_ < &n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer to graph that contains node
    Graph* graph_;

    // Index of node in graph
    size_type uid_;

    // Private Constructor
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type) nodes_.size();
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
    nodes_.push_back(position);
    return Node{this, size()-1};
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < size() && nodes_[n.index()] == n.position());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node{this, i};
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1_ == e.node1_ && node2_ == e.node2_) ||
          (node1_ == e.node2_ && node2_ == e.node1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any interpretive meaning.
     *
     * The edges are ordered first by comparing node1_, then node2_.
     */
    bool operator<(const Edge& e) const {
      return (node1_ < e.node1_) || (node1_ == e.node1_ && node2_ < e.node2_);
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer to graph that contains edge
    Graph* graph_;

    // Nodes of the edge
    size_type node1_;
    size_type node2_;

    // Private constructor
    Edge(const Graph* graph, const Node& a, const Node& b) :
      graph_(const_cast<Graph*>(graph)), node1_(a.index()), node2_(b.index()) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  Edge edge(size_type i) const {
    size_type node1_index = 0;    // index of first node
    size_type node2_index = 0;    // index of second node
    size_type cumul_sum = 0;      // find termination for first node search
    size_type prev_cumul_sum = 0; // cumulative index to get internal index

    // Find first node whose map contains index i
    while(cumul_sum <= i){
      if(edges_.find(node1_index) != edges_.end()){
        prev_cumul_sum = cumul_sum;
        cumul_sum += edges_.at(node1_index).size();
      }
      node1_index += 1;
    }

    node1_index -= 1;

    // Find second node with internal index i-prev_cumul_sum
    for(auto& n2 : edges_.at(node1_index)){
      if(n2.second == i-prev_cumul_sum){
        node2_index = n2.first;
        break;
      }
    }

    return Edge{this, node(node1_index), node(node2_index)};
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool a_b_edge;      // true iff edge (a, b) is in edges_
    bool b_a_edge;      // true iff edge (b, a) is in edges_
    
    a_b_edge = edges_.find(a.index()) != edges_.end()
              && edges_.at(a.index()).find(b.index())
                 != edges_.at(a.index()).end();

    b_a_edge = edges_.find(b.index()) != edges_.end()
               && edges_.at(b.index()).find(a.index())
                  != edges_.at(b.index()).end();

    return a_b_edge || b_a_edge;
  }

  /** Add an edge to the graph, or return the current edge if it already exists
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Does not invalidate outstanding Edge objects.
   *
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b) {

    // Edge (a, b) does not exist, add to Graph
    if(!has_edge(a, b) && !has_edge(b,a)){
      edges_[a.index()][b.index()] = (size_type) edges_[a.index()].size();
      num_edges_ += 1;
    }

    return Edge{this, a, b};
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    num_edges_ = 0;
  }

};

#endif // CME212_GRAPH_HPP

