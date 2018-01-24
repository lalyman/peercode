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

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  std::vector<Point>* nodes_;

  std::vector<std::pair<size_type, size_type>>* edges_list_;
  std::map<std::pair<size_type, size_type>,size_type>* rev_edges_list_;
public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    nodes_ = new std::vector<Point>();
    edges_list_ = new std::vector<std::pair<size_type, size_type>>();
    rev_edges_list_ = new std::map<std::pair<size_type, size_type>,size_type>;
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
      graph_ = nullptr;
      uid_ = size_type (-1);
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      Point& point = (*(graph_->nodes_))[uid_];
      return point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return size_type(uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool result =(graph_== n.graph_);
      result = result && (uid_ == n.uid_);
      return result;
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
      const Point point1 = position();
      const Point point2 = n.position();
      return (point1.x<point2.x)||(point1.x==point2.x && point1.y<point2.y)||(point1.x==point2.x&&point1.y==point2.y&&point1.z<point2.z);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    const Graph*  graph_;

    size_type  uid_;

    Node(const Graph* graph, size_type& uid) {
      graph_ = graph;
      uid_ = uid;
      // HW0: YOUR CODE HERE
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
    return nodes_->size();
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
    size_type old_size = size();
    nodes_->push_back(position);
    Node new_node = Node(this,old_size);

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_==this)&&(n.uid_ < this->size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    const Node node = Node(this,i);
    return node;        // Invalid node
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
      graph_ = nullptr;
      uid_ = size_type(-1);
      id1_ = size_type(-1);
      id2_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(id1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(id2_)  ;    // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_==e.graph_ && ((id1_==e.id1_&&id2_==e.id2_)||(id2_==e.id1_&&id1_==e.id2_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    const Graph * graph_;
    size_type uid_;
    size_type id1_;
    size_type id2_;

    Edge(const Graph* graph, const size_type& uid, const size_type& id1, const size_type& id2) {
      graph_ = graph;
      uid_ = uid;
      id1_ = id1;
      id2_ = id2;
      // HW0: YOUR CODE HERE
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
    return edges_list_->size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    std::pair<size_type, size_type> edge_pair = (*edges_list_)[i];
    const Edge edge = Edge(this,i,edge_pair.first,edge_pair.second);
    return edge;        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    return rev_edges_list_->count(std::pair<size_type, size_type>(a.uid_,b.uid_));
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
    std::map<std::pair<size_type, size_type>,size_type>::iterator found_edge=  rev_edges_list_->find(std::pair<size_type, size_type>(a.uid_,b.uid_));
    if (found_edge== rev_edges_list_->end()){
     const std::pair<size_type, size_type> pair1  = std::pair<size_type, size_type>(a.uid_,b.uid_);
     const std::pair<size_type, size_type> pair2 = std::pair<size_type, size_type>(b.uid_,a.uid_);

     const size_type old_sz = num_edges();
     edges_list_->push_back(pair1);
     rev_edges_list_->insert(std::pair<std::pair<size_type, size_type>,size_type>(pair1,old_sz));
     rev_edges_list_->insert(std::pair<std::pair<size_type, size_type>,size_type>(pair2,old_sz));
     return edge(old_sz);
    }else{
     return edge(found_edge->second);
    }
           // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //destroy items
    nodes_->clear();
    edges_list_->clear();
    rev_edges_list_->clear();
    // HW0: YOUR CODE HERE
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
