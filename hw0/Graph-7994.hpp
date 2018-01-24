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
    size_ = 0;
    num_edges_ = 0;
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
     /* Forward the request to graph_, which looks up the
      * appropriate internal node's position
      */
      return graph_->fetch_node(id_).p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.get_graph() == graph_) && (n.get_id() == id_);
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
      return id_ < n.get_id();
    }

    size_type get_id() const {
      return id_;
   }

   const graph_type* get_graph() const{
     return graph_;
   }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A node is characterized by its graph and its index
    graph_type* graph_; 
    size_type id_;

    Node(const graph_type* g, size_type id) : graph_(const_cast<graph_type*>(g)),
                                              id_(id) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
    // Create a new internal node, add it 
    // to the graph's vector of nodes, and
    // return a proxy associated with this node

    InternalNode n = {position};
    internal_nodes.push_back(n);
    size_ += 1;
    return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.get_graph() == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i && i < size_);
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
      // Forwards the request to graph_
      // which looks up the appropriate internal
      // edge and returns its first node n1_
      return graph_->fetch_edge(index_).n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->fetch_edge(index_).n2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // First check that the two edges are from the same graph
      if(e.get_graph() != graph_) {
        return false;
      }
      else {
        Node n1_ = this->node1();
        Node n2_ = this->node2();
        
        //Edges are undirected, so check both orientations
        bool ans = (n1_ == e.node1() && n2_ == e.node2());
        ans = ans || (n1_ == e.node2() && n2_ == e.node1());
        return ans; 
      }
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return index_ < e.get_index();
    }

    size_type get_index() const {
      return index_;
    }

    const graph_type* get_graph() const{
      return graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // As with Nodes, Edges are characterized by their graph
    // and their index

    graph_type* graph_;
    size_type index_;
    
    Edge(const graph_type* g_, size_type index) : 
        graph_(const_cast<graph_type*>(g_)),
        index_(index) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
      return num_edges_;
    }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    Edge edge(size_type i) const {
      return Edge(this, i);
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    bool has_edge(const Node& a, const Node& b) const {
      // Iterate over all the internal edges, checking each to see
      // if any has endpoints @a a and @a b
      for(auto e : internal_edges){
        if(e.are_endpoints(a, b)){
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
     // Similar to has_edge, iterate over internal_edges
     // returning a proxy associated with the InternalEdge with
     // endpoints a and b if such an InternalEdge exists
     for(size_type i = 0; i < num_edges_; ++i){
        if(internal_edges.at(i).are_endpoints(a,b)){
          return Edge(this, i);
        }
      }
     // If we don't find an edge with these endpoints,
     // we need to add it, increment num_edges_
     // and return a proxy associated with it
  
     internal_edges.push_back(InternalEdge(a,b));
     ++num_edges_;
     return Edge(this, num_edges_ - 1);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    void clear() {
      internal_edges.clear();
      internal_nodes.clear();
      size_ = 0;
      num_edges_ = 0;
    }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
   size_type size_;
   size_type num_edges_;

   struct InternalNode{
     Point p;
     InternalNode(Point pos) : p(pos) {}
   };

   struct InternalEdge{
     node_type n1_;
     node_type n2_;

     InternalEdge(Node n1, Node n2) : n1_(n1), n2_(n2) {}

     // Helper function to check if n1 and n2 
     // are the endpoints of this internal edge.
     bool are_endpoints(const Node& n1, const Node& n2) {
       return((n1_ == n1 && n2_ == n2)) || (n1_ == n2 && n2_ == n1);
     }
   };
   // We store the internal objects in vectors for simplicity.
   // For the edges, an unordered map keyed by the endpoints
   // would be an alternative that would allow has_edge and
   // add_edge to run in O(1).

   std::vector<InternalNode> internal_nodes;
   std::vector<InternalEdge> internal_edges;

   const InternalNode& fetch_node(size_type id) const{
     assert (0 <= id && id < size_);
     return internal_nodes.at(id);
   }

   const InternalEdge& fetch_edge(size_type id) const{
     assert (0 <= id && id < num_edges_);
     return internal_edges.at(id);
   }


};


#endif // CME212_GRAPH_HPP
