#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
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

  /*struct internal_node {
    Point position;
    unsigned int nid;
  };

  struct internal_edge {
    unsigned int nid1, nid2;
    unsigned int eid;
  };

  internal_node* nodes_;
  internal_edge* edges_;
  unsigned int nodes_size_;
  unsigned int edges_size_;
  unsigned int next_nid_;
  unsigned int next_eid_;

  unsigned int nodes_[0][2];
  unsigned int edges_[0][2];*/

  std::vector<Point> nodes_pos_;
  std::vector<unsigned int> nodes_ids_;
  std::vector<std::vector<unsigned int> > edges_nds_;
  std::vector<unsigned int> edges_ids_;
  unsigned int next_nid_;
  unsigned int next_eid_;

  /*std::vector<Point> nodes_pos_;
  std::map<unsigned int, unsigned int> nodes_map_;
  std::vector<std::vector<unsigned int> > edges_nds_;
  std::map<unsigned int, unsigned int> edges_map_;
  unsigned int next_nid_;
  unsigned int next_eid_;*/

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

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
    next_nid_ = 0;
    next_eid_ = 0;
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
      //unsigned int id = index();
      /*Point pos = graph_->nodes_pos_[index()];
      std::cout << "Node " << nid_ << " " << index() << " : [ "
                << pos.x << " "
                << pos.y << " "
                << pos.z << " ]" << std::endl;*/
      return graph_->nodes_pos_[index()];
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      /*std::vector<unsigned>::iterator it;
      it = std::find(graph_->nodes_ids_.begin(), graph_->nodes_ids_.begin(), nid_);
      unsigned int id = std::distance(graph_->nodes_ids_.begin(), it);
      return id;*/
      for (unsigned i=0; i<graph_->nodes_ids_.size(); ++i) {
        if (graph_->nodes_ids_[i] == nid_) {
          return i;
        }
      }
      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ && nid_ == n.nid_) {
        return true;
      }
      //(void) n;          // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
      if (nid_ < n.nid_) {
        return true;
      }
      return false;
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    unsigned int nid_;
    Node(const Graph* graph, unsigned int nid)
        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_ids_.size();
    //return 0;
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
    nodes_pos_.push_back(position);
    nodes_ids_.push_back(next_nid_);
    Node new_node_ = Node(this, next_nid_);
    ++next_nid_;
    /*std::cout << "node " << next_nid_-1 
              << "   " << position.x 
              << "  " << position.y
              << "  " << position.z 
              << "   " << new_node_.graph_ << "  " << new_node_.nid_ << std::endl;*/
    return new_node_;
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
    if (this == n.graph_) {
      for (unsigned i=0; i<nodes_ids_.size(); ++i) {
        if (nodes_ids_[i] == n.nid_) {
          return true;
        }
      }
    }
    return false;
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
    Node new_node = Node(this, nodes_ids_[i]);
    return new_node;
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      unsigned nid = graph_->edges_nds_[edge_index()][0];
      Node new_node = Node(graph_, nid);
      return new_node;
      //return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      unsigned nid = graph_->edges_nds_[edge_index()][1];
      Node new_node = Node(graph_, nid);
      return new_node;
      //return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (eid_ == e.eid_) {
        return true;
      }
      return false;
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (eid_ < e.eid_) {
        return true;
      }
      return false;
      //(void) e;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    unsigned int eid_;

    Edge(const Graph* graph, unsigned eid)
        : graph_(const_cast<Graph*>(graph)), eid_(eid) {
    }
    
    size_type edge_index() const {
      /*std::vector<unsigned>::iterator it;
      it = std::find(graph_->edges_ids_.begin(), graph_->edges_ids_.begin(), eid_);
      unsigned int id = std::distance(graph_->edges_ids_.begin(), it);
      return id;*/
      for (unsigned i=0; i<graph_->edges_ids_.size(); ++i) {
        if (graph_->edges_ids_[i] == eid_) {
          return i;
        }
      }
      return size_type(-1);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_ids_.size();
    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    Edge new_edge = Edge(this, edges_ids_[i]);
    return new_edge;
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
    for (unsigned i=0; i<edges_ids_.size(); ++i) {
      unsigned nid0 = edges_nds_[i][0];
      unsigned nid1 = edges_nds_[i][1];
      if ((nid0==a.nid_ && nid1==b.nid_) || (nid0==b.nid_ && nid1==a.nid_)) {
        return true;
      }
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
    for (unsigned i=0; i<edges_ids_.size(); ++i) {
      unsigned nid0 = edges_nds_[i][0];
      unsigned nid1 = edges_nds_[i][1];
      if ((nid0==a.nid_ && nid1==b.nid_) || (nid0==b.nid_ && nid1==a.nid_)) {
        return edge(i);
      }
    }
    std::vector<unsigned int> nids = {a.nid_, b.nid_};
    edges_nds_.push_back(nids);
    edges_ids_.push_back(next_eid_);
    Edge new_edge = Edge(this, next_eid_);
    ++next_eid_;
    return new_edge;
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_pos_.clear();
    nodes_ids_.clear();
    edges_nds_.clear();
    edges_ids_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
