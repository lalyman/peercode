#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
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
      // check if this position is in node list
      if (uid_<graph_->nodes_.size()){
        return graph_->nodes_[uid_];
      }
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // check if the number is in the range
      if (uid_<graph_->nodes_.size()){
        return uid_;
      }
      assert(false);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((n.graph_==graph_) && (n.index()==uid_));
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
      return ((n.graph_==graph_) && (n.index()>uid_));

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)),uid_(uid){
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
    // add the new node
    nodes_.push_back(position);
    // add node to edge map
    edge_map_[nodes_.size()-1];
    return Node(this,size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ((n.graph_ == this)&&(n.uid_<size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i<num_nodes()){
        return Node(this,i);
	}
	assert(false);
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
      if (uid_<graph_->edges_.size()){
        return graph_->node(graph_->edges_[uid_][0]);
      }
      assert(false);
      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (uid_<graph_->edges_.size()){
        return graph_->node(graph_->edges_[uid_][1]);
      }
      assert(false);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((e.node1()==node1())&&(e.node2()==node2()))||
        ((e.node2()==node1())&&(e.node1()==node2())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->uid_<e.uid_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)),uid_(uid){
    }    


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i<num_edges()){
      return Edge(this,i);
    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (has_node(a) && has_node(b)){
      std::vector<size_type> temp = edge_map_.at(a.index());
      for (size_type i=0; i<temp.size();i++){
        if (((edge(temp[i]).node1()==a) &&(edge(temp[i]).node2()==b)) ||
          ((edge(temp[i]).node1()==b) &&(edge(temp[i]).node2()==a))){
          return true;
        }
      }
    }
    else{
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
    // check if a and b are valid
    if (!(a==b) && has_node(a) && has_node(b)){ 
      // check if edge ab exists    
      if (has_edge(a,b)){
        // read the index of edge with node a
        std::vector<size_type> temp = edge_map_.at(a.index());
        for (size_type i=0; i<temp.size();i++){
          if (((edge(temp[i]).node1()==a) &&(edge(temp[i]).node2()==b)) ||
             ((edge(temp[i]).node1()==b) &&(edge(temp[i]).node2()==a))){
            return edge(temp[i]);
          }
        }      
      }
      
      // add new edge
      else {
        std::vector<size_type> lst_edge;
        lst_edge.push_back(a.index());
        lst_edge.push_back(b.index());
        edges_.push_back(lst_edge);
        edge_map_[a.index()].push_back(num_edges()-1);
        edge_map_[b.index()].push_back(num_edges()-1);
        return edge(num_edges()-1);
      }
    }
    assert(false);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    nodes_.clear();
    edge_map_.clear();
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> nodes_;
  std::vector<std::vector<size_type>> edges_;
  std::map<size_type,std::vector<size_type>> edge_map_;
};

#endif // CME212_GRAPH_HPP
