#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
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
    Node(): uid_(0), graph_(nullptr)  {}
    Node(const graph_type* graph, size_type uid): uid_(uid), graph_(const_cast<graph_type *> (graph)){}



    /** Return this node's position. */
    const Point& position() const {
      if(this->index() < graph_->num_nodes())
        return graph_->nodes_[this->index()];
      assert(false);

    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return size_type(uid_);
    }

    graph_type* graph() const {
      return graph_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->graph_ ==  n.graph() && this->uid_ == n.index());
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
       return (this->uid_ < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
      size_type uid_;
      graph_type * graph_;

    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
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
    num_nodes_ += 1;
    nodes_.push_back(position);
    return Node(this,num_nodes_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i >= num_edges()) {
      std::cout << "In func node: @pre 0 <= @a i < num_nodes()" <<std::endl;
      return Node();
      }
    else
      return Node(this,i);
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
    Edge(): graph_(nullptr), uid_(0), n1_(0), n2_(0){
    }
    /** Construct a valid Edge. */
    Edge(const graph_type* graph, size_type n1, size_type n2, size_type uid):
            graph_(const_cast<graph_type *>(graph)),uid_(uid), n1_(n1), n2_(n2){}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, n2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.node1()==this->node1() && e.node2()==this->node2()) ||
              (e.node1()==this->node2() && e.node2()==this->node1()));
    }


    /** Return this edge's index, a number in the range [0, graph_edge_size). */
    size_type index() const {
      return size_type(uid_);
    }


      /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->index() < e.index());
    }

   private:
      graph_type * graph_;
      size_type uid_;
      size_type n1_, n2_;
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
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
    if(i >= num_edges()) {
      std::cout << "In func edge: @pre 0 <= @a i < num_edges()" <<std::endl;
      return Edge();
    }
    else{
      const std::pair<size_type ,size_type > &p = edges_[i];
      return Edge(this, p.first, p.second, i);
    }

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(a.index() >= num_nodes() || b.index() >= num_nodes()){
      std::cout <<"In func has_edge: @pre @a a and @a b are valid nodes of this graph" <<std::endl;
      return false;
    }

    for(std::vector<std::pair<size_type , size_type >>::const_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      if((it->first == a.index() && it->second == b.index())||
         (it->first == b.index() && it->second == a.index()))
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
    if(!has_edge(a,b)){
      num_edges_ += 1;
      edges_.push_back(std::make_pair(a.index(),b.index()));
      return Edge(this, a.index(), b.index(), num_edges_ - 1);

    }
    else{
      return Edge();
    }
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
    num_edges_ = 0;
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<Point> nodes_;
  std::vector<std::pair<size_type , size_type >> edges_;
  size_type num_nodes_;
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
