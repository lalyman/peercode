#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

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
  class SecretNode;
  class SecretEdge;
  std::vector<SecretNode> Nodes;
  std::vector<SecretEdge> Edges;

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
  	//vectors already initialized to empty vectors
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
      g_ = nullptr;
      ind_ = size_type(-1);
    }

    /** Return this node's position. */
    const Point& position() const {
      return g_->Nodes[ind_].get_position();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return ind_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      (void) n;          // Quiet compiler warning
      return (this->ind_ == n.index()) && (this->g_ == n.g_);
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
      (void) n;           // Quiet compiler warning
      if (*this == n) {
      	//same node so not less than
      	return false;
      }
      else if (this->ind_ == n.index()) {
      	//same index, different graph
      	std::less<Graph*> graph_ptr_less;
      	return graph_ptr_less(this->g_, n.g_);
      }
      else{
      	return (this->ind_ < n.index());
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* g_;
    size_type ind_;

    Node(const Graph* g, size_type ind) {
      g_ = const_cast<Graph*>(g);
      ind_ = ind;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return Nodes.size();
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
    (void) position;      // Quiet compiler warning
    size_type i = Nodes.size();
    Nodes.push_back(SecretNode(position));
    return Node(this, i);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    if (n.index() < num_nodes()) {
        return (this->node(n.index()) == n);
    }
    else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    (void) i;             // Quiet compiler warning
    return Node(this,i);        // Invalid node
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
      g_ = nullptr;
      eid_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, g_->Edges[eid_].get_node1());
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, g_->Edges[eid_].get_node2());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      Node n1 = this->node1();
      Node n2 = this->node2();

      return (((n1 == e.node1()) && (n2 == e.node2())) 
      	|| ((n1 == e.node2()) && (n2 == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if(*this == e) {
      	//same edge
      	return false;
      }
      else if (this->eid_ == e.eid_) {
      	//same edge id, different graphs
      	std::less<Graph*> graph_ptr_less;
      	return graph_ptr_less(this->g_, e.g_);
      }
      else{
      	//different edge ids
      	return (this->eid_ < e.eid_);
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* g_;
    size_type eid_;

    Edge(const Graph* g, size_type eid) {
      g_ = const_cast<Graph*>(g);
      eid_ = eid;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return Edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    (void) i;             // Quiet compiler warning
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    (void) a; (void) b;   // Quiet compiler warning
    for (SecretEdge e : Edges) {
      if(((e.get_node1() == a.index()) && (e.get_node2() == b.index())) || 
         ((e.get_node2() == a.index()) && (e.get_node1() == b.index()))) {
        return true;
      }
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
    (void) a, (void) b;   // Quiet compiler warning
    for (size_type i = 0; i < Edges.size(); ++i) {
      SecretEdge e = Edges[i];
      if(((e.get_node1() == a.index()) && (e.get_node2() == b.index())) || 
         ((e.get_node2() == a.index()) && (e.get_node1() == b.index()))) {
        return Edge(this, i);
      }
    }
    size_type eid = Edges.size();
    Edges.push_back(SecretEdge(a.index(), b.index()));
    return Edge(this, eid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Nodes.clear();
    Edges.clear();
  }

 private:
  //Class to represent Nodes internally
  class SecretNode {
    private:
        Point p_;
    public:
      const Point& get_position() const {
        return p_;
      }
      SecretNode(const Point& position) {
        p_ = position;
      }
  };
  //Class to represent Edges internally
  class SecretEdge {
    private:
        size_type n1_;
        size_type n2_;
    public:
      size_type get_node1() const {
        return n1_;
      }
      size_type get_node2() const {
        return n2_;
      }
      SecretEdge(const size_type a, const size_type b) {
        n1_ = a;
        n2_ = b;
      }
  };
};

#endif // CME212_GRAPH_HPP
