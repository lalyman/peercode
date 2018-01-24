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

  std::vector<Point> Nodes;
  std::vector<std::vector<unsigned>> Edges;
  unsigned size0;
  unsigned num_edges0;

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
    size0 = 0;
    num_edges0 = 0;
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
    Node() {
      uid_n = 0;
      //graph_n = nullptr;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      const Point& result_point = graph_n->Nodes[uid_n];
      return result_point;
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_n;
      //return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_n == n.graph_n) && (uid_n == n.index()));
      //(void) n;          // Quiet compiler warning
      //return false;
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
      assert(graph_n == n.graph_n);
      return uid_n < n.index();
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type uid_n;
    graph_type* graph_n;

    Node(const graph_type* graph0, size_type uid0): graph_n(const_cast<graph_type*>(graph0)), uid_n(uid0){
    } 
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
    // HW0: YOUR CODE HERE
    return size0;
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
    // HW0: YOUR CODE HERE
    size0 += 1;
    Point result_point= const_cast<Point&>(position);
    Nodes.push_back(result_point);
    return Node(this, size0 - 1);
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.index() < size0);
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((0 <= i) && (0 < size0));
    return Node(this, i);
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
      node1_id = 0;
      node2_id = 0;
      edge_id = 0;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_e->node(node1_id);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_e->node(node2_id);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1()) && (node2() == e.node2())) || ((node1() == e.node2()) && (node2() == e.node1())) ;
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((node1().index() + node2().index()) < (e.node1().index() + e.node2().index()));
      //(void) e;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_e;
    size_type node1_id;
    size_type node2_id;
    size_type edge_id;
    Edge(const graph_type* graph_e0, size_type node1_id0, size_type node2_id0, size_type edge_id0): graph_e(const_cast<graph_type*>(graph_e0)), node1_id(node1_id0), node2_id(node2_id0), edge_id(edge_id0){
      assert(node1_id0 < graph_e->size0 && node2_id0 < graph_e->size0 && node1_id0 != node2_id0);
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges0);
    return Edge(this, Edges[i][0], Edges[i][1], i);
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_n == b.graph_n);
    size_type a0 = a.index();
    size_type b0 = b.index();
    bool test = false;
    size_type i = 0;
    while(i < size0){
      test = ((a0 == Edges[i][0] && b0 == Edges[i][1]) || (a0 == Edges[i][1] && b0 == Edges[i][0]));
      if(test){
        return true;
      }
      i++;
    }
    return false;
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;
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
    assert(a.graph_n == b.graph_n);
    size_type a0 = a.index();
    size_type b0 = b.index();
    bool test = false;
    size_type i = 0;
    while(i < num_edges0){
      test = ((a0 == Edges[i][0] && b0 == Edges[i][1]) || (a0 == Edges[i][1] && b0 == Edges[i][0]));
      if(test){
        return Edge(this, a.index(), b.index(), i);
      }
      i++;
    }
    num_edges0 += 1;
    std::vector<size_type> temp = {a0, b0};
    Edges.push_back(temp);
    return Edge(this, a.index(), b.index(), num_edges0 - 1);
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Nodes.clear();
    Edges.clear();
    size0 = 0;
    num_edges0 = 0; 
    // HW0: YOUR CODE HERE
  }

 //private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
