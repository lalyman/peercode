#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
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
  Graph() : n_nodes(0), n_edges(0), next_nid(0), next_eid(0),
            nodes(new std::unordered_map<size_type, Point>()),
            edges(new std::unordered_map<
                    size_type,
                    std::pair<size_type, size_type>>()),
            edges_index(new std::map<
                          std::pair<size_type, size_type>,
                          size_type>()){
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
    Node(){
      // HW0: DONE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: DONE
      // std::cout << "fetching position of node " << index() << '\n';
      return graph_->get_position(index());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: DONE
      return nid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: DONE
      return index() == n.index() && graph_ == n.graph_;
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
      // HW0: DONE
      return index() < n.index();
    }

   private:
    size_type nid_;
    graph_type* graph_;

    Node(const size_type id, const graph_type* graph) :
        nid_(id), graph_(const_cast<graph_type*>(graph)) {
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /* HW0: DONE
    Use this space to declare private data members and methods for Node
    that will not be visible to users, but may be useful within Graph.
    i.e. Graph needs a way to construct valid Node objects */
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: DONE
    return n_nodes;
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
    // HW0: DONE
    std::pair<size_type, Point> newpair(next_nid, position);
    nodes->insert(newpair);
    ++n_nodes;
    return Node(next_nid++, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: DONE
    return nodes->find(n.index()) != nodes->end();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: DONE
    assert(i < size());
    return Node(i, this);
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
    Edge(){
      // HW0: DONE
    }

    size_type index() const {
      return eid_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: DONE
      return Node(nid_a_, graph_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: DONE
      return Node(nid_b_, graph_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return nid_a_ == e.nid_a_ && nid_b_ == e.nid_b_ && graph_ == e.graph_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return index() < e.index();
    }

   private:

     size_type eid_, nid_a_, nid_b_;
     graph_type* graph_;

     Edge(const size_type id, const size_type id1, const size_type id2, const graph_type* graph) :
      eid_(id), graph_(const_cast<graph_type*>(graph)){
        if (id1 < id2) {
          nid_a_ = id1;
          nid_b_ = id2;
        }
        else {
          nid_a_ = id2;
          nid_b_ = id1;
        }
      }

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    /* HW0: DONE
    * Use this space to declare private data members and methods for Edge
    * that will not be visible to users, but may be useful within Graph.
    * i.e. Graph needs a way to construct valid Edge objects
    */
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: DONE
    return n_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    std::pair<size_type, size_type> pair_i(edges->at(i));
    return Edge(i, pair_i.first, pair_i.second, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: DONE
    std::pair<size_type, size_type> pair(a.index(), b.index());
    return edges_index->find(pair) != edges_index->end();
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
    // HW0: DONE
    size_type index;
    if (has_edge(a, b)) {
      std::pair<size_type, size_type> pair(a.index(), b.index());
      return Edge(edges_index->at(pair), a.index(), b.index(), this);
    }
    else {
      index = next_eid;
      std::pair<size_type, std::pair<size_type, size_type>> newpair(
        next_eid,
        std::pair<size_type, size_type>(a.index(), b.index()));
      std::pair<std::pair<size_type, size_type>, size_type> newpair_index1(
        std::pair<size_type, size_type>(a.index(), b.index()),
        next_eid
      ), newpair_index2(
        std::pair<size_type, size_type>(b.index(), a.index()),
        next_eid
      );
      edges->insert(newpair);
      edges_index->insert(newpair_index1); edges_index->insert(newpair_index2);
      ++n_edges;
      return Edge(next_eid++, a.index(), b.index(), this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes->clear(); edges->clear(); edges_index->clear();
    next_nid = 0;  next_eid = 0; n_nodes = 0; n_edges = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  size_type n_nodes, n_edges;
  size_type next_nid, next_eid;
  std::unordered_map<size_type, Point>* nodes;
  std::unordered_map<size_type, std::pair<size_type, size_type>>* edges;
  std::map<std::pair<size_type, size_type>, size_type>* edges_index;

  Point& get_position(const size_type id){
    return nodes->at(id);
  }

};

#endif // CME212_GRAPH_HPP
