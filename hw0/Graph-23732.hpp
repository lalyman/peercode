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

  // Predeclaring the interal structures. 
  struct internal_node; 
  struct interal_edge;

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
      : nodes(), edges() {  
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
    }

    /** Return this node's position. */
    const Point& position() const {
     return (graph->nodes[this-> node_id_]).position; 
    } 

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
     return (graph->nodes[this-> node_id_]).node_id; 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph == this -> graph && n.node_id_ == this -> node_id_);
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
      return (n.graph == this -> graph && this -> node_id_ < n.node_id_);
    }

  private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Allow Graph to access Node's private member data and functions.
    friend class Graph; 
    Graph* graph;
    size_type node_id_;
    // Private Node constructor.
    Node(const Graph* g, size_type id)
      : graph(const_cast<Graph*>(g)), node_id_(id) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
    internal_node new_node;
    new_node.position = position; 
    new_node.node_id = size();
    nodes.push_back(new_node);
    return Node(this,new_node.node_id); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return(n.graph == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i<size());           
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return (graph->edges[this-> edge_id_]).first_node; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return (graph->edges[this-> edge_id_]).second_node; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.edge_id_ == this -> edge_id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this -> edge_id_ < e.edge_id_);
    }

   private:
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type edge_id_;
    size_type first_node_id;
    size_type second_node_id;
    // Private Edge constructor.
    Edge(const Graph* g, size_type id, size_type first_node, size_type second_node)
      : graph(const_cast<Graph*>(g)), edge_id_(id), first_node_id(first_node), 
      second_node_id(second_node) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<edges.size());           
    return Edge(this, i, edges[i].first_node.index(), edges[i].second_node.index());
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(a.graph != this || b.graph != this) {
      return false;
    } else {
       for (unsigned i = 0; i < edges.size(); i++) {
          if ((edges[i].first_node == a && edges[i].second_node == b) ||
           (edges[i].first_node == b && edges[i].second_node == a))  {return true;}
       }
       return false;
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
   // Checking if a and b are distinct.  
   if (a == b) {
    return Edge();
   } 
   // Checking if the current edge already exists.
   else if (has_edge(a,b)) {
    for (unsigned i = 0; i < edges.size(); i++) {
     if ((edges[i].first_node == a && edges[i].second_node == b) ||
     (edges[i].first_node == b && edges[i].second_node == a)) {
      return Edge(this, i, a.index(), b.index());
     } 
    }
   }
   // If distinct and does not exist, then create a new node and add it.
   else {
      internal_edge new_edge;
      new_edge.first_node = a; 
      new_edge.second_node = b;
      new_edge.edge_id = edges.size();
      edges.push_back(new_edge);
      return Edge(this,new_edge.edge_id, new_edge.first_node.index(),
      new_edge.second_node.index());
    } return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this -> nodes.clear();
    this -> edges.clear();
  }

 private:
  // Internal type node and edge.
  struct internal_node {
    Point position; // Position of the node.
    size_type node_id; // The identification for the node. Each is unique.
  };

  struct internal_edge {
  Node first_node; // First node of the edge.
  Node second_node;  //Second node of the edge.
  size_type edge_id; // The identification for the edge. Each is unique.
  };

  // Initailizing the vectors of the interal structs. 
  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;
  
};

#endif // CME212_GRAPH_HPP
