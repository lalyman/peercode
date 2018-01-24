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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

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
  Graph()
   : next_node_uid_(0), next_edge_uid_(0) {
  };

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
    Node()
      : graph_(nullptr){
    }

    /** Return this node's position. */
    const Point& position() const {
      size_type node_id = this->uid_;
      return graph_->nodes_[node_id].location;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(this->graph_ == n.graph_)
        if(this->uid_ == n.uid_)
          return true;
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
      //We dont check if this is part of the same graph since it needs to work
      //on a global order, we simply check uid_s
      if(this->uid_ < n.uid_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    //Pointer back to the graph container
    Graph* graph_;
    //This element's unique identification number
    size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

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

    internal_node_element node_to_add = {position, next_node_uid_};
    nodes_.push_back(node_to_add);
    ++next_node_uid_;
    return node(next_node_uid_ - 1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(node(n.uid_) == n)
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
    assert(i < nodes_.size());
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
    Edge()
      : graph_(nullptr){
    }

    /** Return a node of this Edge */
    Node node1() const {
      size_type edge_uid = this->uid_;
      internal_edge_element edge_like_object = graph_->edges_[edge_uid];
      size_type node_uid = edge_like_object.node1_uid;
      Node node_to_return = graph_->node(node_uid);
      return node_to_return;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type edge_uid = this->uid_;
      internal_edge_element edge_like_object = graph_->edges_[edge_uid];
      size_type node_uid = edge_like_object.node2_uid;
      Node node_to_return = graph_->node(node_uid);
      return node_to_return;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      //if both edges empty, true
      if(this->graph_ == nullptr && e.graph_ == nullptr)
        return true;

      //if either edge empty (but not both) throw false
      if(this->graph_ == nullptr && e.graph_ != nullptr)
        return false;
      if(this->graph_ != nullptr && e.graph_ == nullptr)
        return false;

      //else go ahead and access the edges, checking node equality
      if(this->node1().uid_ == e.node1().uid_
      && this->node2().uid_ == e.node2().uid_)
        return true;
      if(this->node1().uid_ == e.node2().uid_
      && this->node2().uid_ == e.node1().uid_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //returns true if either node1 is less than other node1, OR
      //returns true if node1s are tied but node2 is less
      if(this->node1().uid_ < e.node1().uid_)
        return true;
      if(this->node1().uid_ == e.node1().uid_
      && this->node2().uid_ <  e.node2().uid_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    //Pointer back to the graph container
    Graph* graph_;
    //This element's unique identification number
    size_type uid_;

    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < edges_.size());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(edge_finder(a,b) == Edge()){
      return false;
    } else {
      return true;
    }
  }

  Edge edge_finder(const Node& a, const Node& b) const {
    size_type node1_uid = a.uid_;
    size_type node2_uid = b.uid_;

    //Exhaustive search for matching edge
    for(size_type i = 0; i < edges_.size(); ++i){
      size_type comparison_node1_uid = edge(i).node1().uid_;
      size_type comparison_node2_uid = edge(i).node2().uid_;

      if(node1_uid == comparison_node1_uid && node2_uid == comparison_node2_uid)
        return edge(i);
      if(node1_uid == comparison_node2_uid && node2_uid == comparison_node1_uid)
        return edge(i);
    }

    return Edge();
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
    //don't need to check node validity, that is pre-condition assumed true

    //find the edge if it already exists
    Edge edge_if_exists = edge_finder(a,b);

    //add it if it doesn't exit already
    if(edge_if_exists == Edge()){

      internal_edge_element edge_to_add = {a.uid_, b.uid_, next_edge_uid_};
      edges_.push_back(edge_to_add);
      ++next_edge_uid_;

      //if(next_edge_uid_%1000 == 0)
      //  std::cout << next_edge_uid_ << std::endl;

      return edge(next_edge_uid_ - 1);
    //otherwise, return the edge
    }else{
      return edge_if_exists;
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    nodes_.clear();
    next_node_uid_ = 0;
    next_edge_uid_ = 0;
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //internal node & edge like structures, will store an array of them
  struct internal_node_element{
    Point location; //The point held by a node
    size_type uid; //The unique identification of an element
  };

  struct internal_edge_element{
    size_type node1_uid; //first node of edge
    size_type node2_uid; //second node of edge
    size_type uid; //The unique identification of an element
  };

  std::vector<internal_node_element> nodes_;
  std::vector<internal_edge_element> edges_;
  size_type next_node_uid_;
  size_type next_edge_uid_;

};

#endif // CME212_GRAPH_HPP
