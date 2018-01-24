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

using namespace std;

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
  }

  /** Default destructor */
  ~Graph() = default;

  // =========================================================================
  // NODES
  // =========================================================================

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
        return graph_->Positions_[uid_];
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
        return ((n.uid_ == uid_) && (n.graph_ == graph_));
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

        if ( graph_ == n.graph_) {
            //if nodes belongs to the same graph, check the id of the nodes
            return (uid_ < n.uid_);}
        else {
            //check in global sense which node is greater without any disambiguity
            //using the pointer to the graph
        return (graph_ < n.graph_);}
    }

   private:
   // Allow Graph to access Node's private member data and functions.

   //Pointer to the graph
   Graph* graph_;
   //Id of the node
   size_type uid_;
   //Additional constructor
   Node(const Graph* graph, size_type id)
        : graph_(const_cast<Graph*>(graph)), uid_(id) {
   }

   friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  // ==========================================================================================================
  // ==========================================================================================================

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
      return Positions_.size();
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
      Positions_.push_back(position);
      return Node(this, num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

      //check also if they both belogs to the same graph
      //and if the index is within range
      return ( (n.graph_ == this) && (size()-1 >= n.index()) );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      assert(i < num_nodes());
      return Node(this,i);
  }

  // =========================================================================
  // EDGES
  // =========================================================================

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
        return  graph_->node(graph_->EdgesPosition_[uid_].first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return  graph_->node(graph_->EdgesPosition_[uid_].second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        //check that both belongs to the same graph and also the indeces
        return ((e.uid_ == uid_) && (e.graph_ == graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if ( graph_ == e.graph_) {
            //if edges belongs to the same graph, check the id of the edges
            return (uid_ < e.uid_);}
        else {
            //check in global sense which edge is greater without any disambiguity
            //using the pointer to the graph
        return (graph_ < e.graph_);}
    }

   private:
   // Allow Graph to access Edge's private member data and functions.
   //Pointer to the graph
   Graph* graph_;
   //Id of the edge
   size_type uid_;
   //Additional constructor
   Edge(const Graph* graph, size_type id)
        : graph_(const_cast<Graph*>(graph)), uid_(id) {
   }

   friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  // ========================================================================
  // ========================================================================

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      // return n_edges ;
      return EdgesPosition_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      assert(i < num_edges());
      return Edge(this,i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

      pair<size_type, size_type> p1;
      if (a.index() > b.index()) { p1 = std::make_pair(a.index(), b.index());}
      else { p1 = std::make_pair(b.index(), a.index());}

      if ((std::find(EdgesPosition_.begin(), EdgesPosition_.end(), p1) != EdgesPosition_.end())) {
          return true;}
      else return false;
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

      if (a==b) { std::cerr << "Node a and Node b have to be distinct \n";}

      if (has_edge(a, b)) {
          return Edge(this,whichIndex( a, b));
      }else{
          if  (a.index()> b.index() ) {
              EdgesPosition_.push_back(std::make_pair(a.index(), b.index()));
          }else{EdgesPosition_.push_back(std::make_pair(b.index(), a.index()));}
              return Edge(this,num_edges()-1);
      }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      Positions_.clear();
      EdgesPosition_.clear();     
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> Positions_;
  std::vector< pair<size_type, size_type> > EdgesPosition_;

  size_type whichIndex(const Node& a, const Node& b) const {
      pair<size_type, size_type> p1;
      if (a.index() > b.index()) { p1 = std::make_pair(a.index(), b.index());}
      else { p1 = std::make_pair(b.index(), a.index());}

      auto it = std::find(EdgesPosition_.begin(), EdgesPosition_.end(), p1) ;
      auto idx = std::distance(EdgesPosition_.begin(), it);
      return idx;
  }
};

#endif // CME212_GRAPH_HPP
