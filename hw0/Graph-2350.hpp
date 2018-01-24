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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Data structures that store information for the Node and Edge class
  struct internal_node;
  struct internal_edge;

  //?Disable copy and assignment?

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
   * Return type of Graph::Node::index(), Graph::num_nodes(),
   * Graph::num_edges(), and argument type of Graph::node(size_type)
   */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
  }

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
    Node() {
      // HW0: YOUR CODE HERE
    }

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

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_->has_node(n) && (n.index_ == index_);
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
      return graph_->has_node(n) && (n.index_ < index_);
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date Position information
    // @pre this node is valid (0 <= index_ < graph_->num_nodes())
    internal_node& fetch() const {
      assert(index_ < graph_->num_nodes());
      return graph_->nodes_[index_];
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* graph_;
    size_type index_;

    /** Private Constructor */
    Node(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
    }
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
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    // Add the new node object
    node_objs_.push_back(Node(this, num_nodes()));
    // Add the new internal node to graph
    nodes_.push_back({position});
    // Return a node that associates with the new internal one
    return node_objs_.back(); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
    return node_objs_[i];
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    node_type node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node_objs_[fetch().node1];
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node_objs_[fetch().node2];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      return node1() == e.node1() && node2() == e.node2();
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      return node1() < e.node1();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type index_;

    // Private constructor
    Edge(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date information
    // @pre this node is valid (0 <= index_ < graph_->num_edges())
    internal_edge& fetch() const {
      assert(index_ < graph_->num_edges());
      return graph_->edges_[index_];
    }

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
    assert(0 <= i && i < num_edges());
    return edge_objs_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    // HW0: YOUR CODE HERE
    // In the double-layered nodes_edge_map_, the smaller index
    // of the two nodes for each edge is used as key in first layer
    size_type s = std::min(a.index_, b.index_);
    size_type l = std::max(a.index_, b.index_);
    if (nodes_edge_map_.count(s) == 1)
      if (nodes_edge_map_.at(s).count(l) == 1)
        return true;
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
  edge_type add_edge(const node_type& a, const node_type& b) {
    // HW0: YOUR CODE HERE
    assert(this->has_node(a) && this->has_node(b));
    // Check if such edge already exists
    size_type s = std::min(a.index_, b.index_);
    size_type l = std::max(a.index_, b.index_);
    if (has_edge(a, b))
      return edge_objs_[nodes_edge_map_[s][l]];
    // Add the new edge to nodes_edge_map_
    nodes_edge_map_[s][l] = num_edges();
    // Add the new edge object
    edge_objs_.push_back(Edge(this, num_edges()));
    // Add the new internal edge to graph
    edges_.push_back({s, l});

    // Return an edge that associates with the new internal one
    return edge_objs_.back();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    Point position;
  };

  struct internal_edge {
    const size_type node1;
    const size_type node2;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<Node> node_objs_;
  std::vector<Edge> edge_objs_;
  std::map<size_type, std::map<size_type, size_type>> nodes_edge_map_;
};

#endif // CME212_GRAPH_HPP
