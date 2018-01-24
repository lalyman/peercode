#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // access nodes vector containing Point data
      return (*graph_).nodes[index_];
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
      // sufficient to check node_index since we don't have remove node method
      return (index_ == n.index()) && (graph_ == n.graph_);
      
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
      return (index_ < n.index()) && (graph_ == n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    const graph_type* graph_; // Pointer back to the Graph container - 8 bytes in 64 bit
    size_type index_; // This node's index - 4 bytes as unsigned integer

    // Private constructor for friend class Graph
    Node(const graph_type* graph, size_type index)
      : graph_(graph), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // return size of nodes vector
    return nodes.size();
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
    // push back position into nodes vector
    // return construction using new node size
    nodes.push_back(position);
    return Node(this, nodes.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // check if valid index
    if  (n.index() >= nodes.size())
      return false;
    return (node(n.index()) == n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // return Node construction with index
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // access node1 data of edge vector
      return (*graph_).edges[index_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      //access node2 data of edge vector
      return (*graph_).edges[index_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // sufficient to check edge_index since we don't have remove edge method
      return (index_ == e.index_) && (graph_ == e.graph_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (graph_ == e.graph_) && (index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    const graph_type* graph_; // Pointer back to the Graph container - 8 bytes
    size_type index_; // This edge's index - 4 bytes

    // Private constructor for friend class Graph
    Edge(const graph_type* graph, size_type index)
      : graph_(graph), index_(index) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // return size of edges vector
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // return Edge construction given index
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    /** check if adjacency map contains node a
     *  if not, given edge is not in graph
     *  if yes, return if node b is in adjacency map of node a
     */

    auto adj_nodesA = adj_map.find(a.index());
    if (adj_nodesA == adj_map.end()) {
      return false;
    }
    else {
      auto edgeB = (adj_nodesA->second).find(b.index());
      return (edgeB != (adj_nodesA->second).end());
    }
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
    /** Check if adjacency map contains node a
      *   if yes, check if node b is in adjacency map of node a
      *       - if yes, return edge index contained in map
              - if no, insert node b and edge index in adjacency map of node a
        else,
            insert node a in adjacency map as key and adjacency map of node a as value, initialized to contain edge to b 
            check if adjacency map contains node b
              - if yes, insert node a and edge index in adjacency map of b
              - if no, insert node b in adjacency map as key and adjacency map of node b as value, initialized to to contain edge to a
    */

    size_type size_edges = edges.size();
    size_type indexA = a.index();
    size_type indexB = b.index();

    // pointers to check if adjacency map contains nodes a, b
    auto adj_nodesA = adj_map.find(indexA);
    auto adj_nodesB = adj_map.find(indexB);

    // conditions on existence in adjacency maps, with appropriate insertions
    if (adj_nodesA != adj_map.end()) {
      auto edgeB = (adj_nodesA->second).find(indexB);
      if (edgeB != (adj_nodesA->second).end()) {
        return edge(edgeB->second);
      } else {
        (adj_nodesA->second).insert({indexB, size_edges});
      }
    } else {
      std::unordered_map<size_type,size_type> mapA = {{indexB, size_edges}};
      adj_map.insert({indexA, mapA}); 

      if (adj_nodesB != adj_map.end()) {
        (adj_nodesB->second).insert({indexA, size_edges});
      } else {
        std::unordered_map<size_type,size_type> mapB = {{indexA, size_edges}};
        adj_map.insert({indexB, mapB});
      }
    }

    // add edge element to edges vector, to keep track of edge indices
    internal_edge edge = {a, b};
    edges.push_back(edge);

    return Edge(this, size_edges);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    edges.clear();
    adj_map.clear();
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // internal type for edges
  struct internal_edge {
    Node node1;
    Node node2;
  };

  std::vector<Point> nodes; // vector of Nodes, index in vector is node_index
  std::vector<internal_edge> edges; // vector of Edges, index in vector is edge_index
  std::unordered_map<size_type,std::unordered_map<size_type,size_type>> adj_map; // adjacency map with following structure - {node_index: {adjacent_node_index: edge_index}}
};

#endif // CME212_GRAPH_HPP
