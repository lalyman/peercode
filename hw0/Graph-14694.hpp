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
  using size_type = unsigned;
  /** Type of this graph. */
  using graph_type = Graph;
  /** Predeclaration of Node and Edge type. */
  /** Synonym for Node and Edge (following STL conventions). */
  class Node;
  using node_type = Node;
  class Edge;
  using edge_type = Edge;

private: 
  std::vector<Point> node_elements_;
  std::vector<Edge> edge_elements_;
  size_type node_size_;
  size_type edge_size_;

public:

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */

  // CONSTRUCTORS AND DESTRUCTOR
  /** Construct an empty graph. */
  Graph(): node_elements_(), edge_elements_(), node_size_(0), edge_size_(0){ }

  /** Default destructor */
  ~Graph() = default;


  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class reXpresenting the graph's nodes.
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
    Node() { }


    /** Return this node's position. */
    const Point& position() const { 
        return graph_address_->node_elements_[index_]; 
    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return index_; }
    const Graph* graph_address() const {return graph_address_;}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.index() == index_ and n.graph_address() == graph_address_);
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
        if ( n.graph_address() != graph_address_){ 
            return false; 
        } else if (n.index() > index_){
            return true;
        } else {
            return false;
        }

    }
   private:
    friend class Graph;
    Graph* graph_address_;
    size_type index_;
    Node(const Graph* graph_address, const size_type index): 
        graph_address_(const_cast<Graph*>(graph_address)), index_(index) {}
  };

  /** Return the number of nodes in the graph.
   */
  size_type size() const {
    return node_size_;
  }
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
    node_elements_.emplace_back(position);
    node_size_++;
    return Node(this, node_size_-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ( this == n.graph_address() );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < node_size_){
        return Node(this, i); 
    } else {
        return Node();
    }
  }

public:


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
    Edge() {  }

    /** Return a node of this Edge */
    Node node1() const {
        return *node1_address_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return *node2_address_;
    }


    size_type index_edge() const {
        return index_edge_;
    } 

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if ( e.node1() == *node1_address_ and e.node2() == *node2_address_){
            return true;
        } else if ( e.node2() == *node1_address_ and e.node1() == *node2_address_){
            return true;
        } else {
            return false;
        }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if (index_edge_ < e.index_edge()) {
            return true;
        } else {
            return false;
        }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Node* node1_address_;
    Node* node2_address_;
    size_type index_edge_;

    Edge(const Node* node1_address, const Node* node2_address, const size_type index_edge): 
        node1_address_(const_cast<Node*>(node1_address)), 
        node2_address_(const_cast<Node*>(node2_address)), index_edge_(index_edge) {}

  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size_;  
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
        return edge_elements_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
        for (size_type i = 0; i < edge_size_; i++){
            Edge current_edge = edge_elements_[i];
            if (current_edge.node1() == a and current_edge.node2() == b){
                return true;
            } else if (current_edge.node2() == a and current_edge.node1() == b) {
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

    if (a.graph_address() != b.graph_address()){
        return Edge();
    } else {
        if(has_edge(a,b) == true){
            return Edge(&a, &b, edge_size_);
        } else {
            Edge new_edge = Edge(&a, &b, edge_size_);
            edge_elements_.emplace_back(new_edge);
            edge_size_++;
            return new_edge;
        }
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_elements_ = {}; 
    node_size_ = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
