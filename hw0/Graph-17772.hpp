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
  struct internal_node;
  struct internal_edge;
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
  Graph()
    // HW0: YOUR CODE HERE
    : nodes_(), edges_(), size_(0), next_uid_(0), nedges_(0), next_euid_(0) {
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
      // NO CODE ADDED - CREATE INVALID NODE
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE (changed return line)
      return fetch().pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE (changed return line)
      return fetch().uid;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (this->graph_ == n.graph_ &&
	      this->uid_ == n.uid_);
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
      return (this->uid_ < n.fetch().uid);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; //Pointer to graph container
    size_type uid_; //node's identification number
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    // Function to return node information
    internal_node& fetch() const {
      for (size_type i = 0; i < graph_->size(); ++i)
	if (graph_->nodes_[i].uid == uid_)
	  return graph_->nodes_[i];
      assert(false);
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
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
    internal_node* new_nodes = new internal_node[size_ + 1];
    // Copy the current nodes  into a new array
    for (size_type i = 0; i < size_; ++i)
      new_nodes[i] = nodes_[i];
    //Set the point and uid for the new element
    new_nodes[size_].pt = position;
    new_nodes[size_].uid = next_uid_;
    //Delete the old nodes and reassign values
    delete[] nodes_;
    nodes_ = new_nodes;
    ++size_;
    ++next_uid_;
    return Node(this, next_uid_-1); // Points to new element
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size());
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
      // HW0: YOUR CODE HERE (done, nothing added)
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return fetch().nd1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch().nd2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((this->node1() == e.node1()) && (this->node2() == e.node1())) ||
	      ((this->node2() == e.node2()) && (this->node1() == e.node2())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (uid_ < e.fetch().uid);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_; // pointer back to the graph
    size_type uid_; // edge's identification number
    //Constructor
    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    //fetching function to return information about the edge
    internal_edge& fetch() const {
      for (size_type i=0; i<graph_->num_edges(); ++i)
	if (graph_->edges_[i].uid == uid_)
	  return graph_->edges_[i];
      assert(false);
    }
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE --> changed the return statement
    return nedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());          
    return Edge(this, i);        // Return Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Loop through edges, if the right one is found break the loop
    bool edge_exists = false;
    for (size_type i = 0; i<num_edges(); ++i) {
      if (((edges_[i].nd1 == a) && (edges_[i].nd2 == b)) ||
	  ((edges_[i].nd2 == a) && (edges_[i].nd1 == b))) {
	edge_exists = true;
	break;
      }
    }
    
    return edge_exists;
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
    size_type edgenum = 0;
    if (has_edge(a,b) == false) {
      internal_edge* new_edges = new internal_edge[nedges_ + 1];
      // Copy edges into a new array
      for (size_type i = 0; i < nedges_; ++i)
	new_edges[i] = edges_[i];
      // Set the nodes and uid for the new edge
      new_edges[nedges_].nd1 = a;
      new_edges[nedges_].nd2 = b;
      new_edges[nedges_].uid = next_euid_;
      // Delete the old edges and reassign values
      delete[] edges_;
      edges_ = new_edges;
      ++nedges_;
      ++next_euid_;
      edgenum = next_euid_-1;
    } else {
      for (size_type i = 0; i<num_edges(); ++i) {
	if (((edges_[i].nd1 == a) && (edges_[i].nd2 == b)) ||
	    ((edges_[i].nd2 == a) && (edges_[i].nd1 == b))) {
	  edgenum = edges_[i].uid;
	  break;
	}
      }
    }
    // Returns an edge that points to the new element
    return Edge(this, edgenum);       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    delete[] edges_;
    delete[] nodes_;
    size_ = 0; next_uid_ = 0; nedges_ = 0; next_euid_ = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // Internal type for set node

  // internal node struct contains a point and uid
  struct internal_node {
    Point pt;
    size_type uid;
  };

  // internal edge struct contains two nodes and a uid
  struct internal_edge {
    Node  nd1;
    Node  nd2;
    size_type uid;
  };

  internal_node* nodes_;
  internal_edge* edges_;
  size_type size_;
  size_type next_uid_;
  size_type nedges_;
  size_type next_euid_;

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
