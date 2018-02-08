#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Type of Nodes */
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //** Type of data structure maintaining edge data, a nested unordered_map */
  using adjacency_map_type = std::unordered_map<size_type,
      std::unordered_set<size_type>>;

  //
  // PRIVATE INTERNAL DATA
  //

 private:
  /** Internal structure for holding Node data */
  struct internal_node {
    Point position;
    node_value_type value;
  };

  /** Internal structure for holding Edge data */
  struct internal_edge {
    size_type node1;
    size_type node2;
  };

  /** Vector of internal_node structs */
  std::vector<internal_node> nodes_;

  //** Vector of edge data allowing fast lookup by index */
  std::vector<internal_edge> edges_;

  /** Adjacency map that holds edges with each node appearing as a key */
  adjacency_map_type adjacency_map_;

  /** Maintain number of edges */
  size_type num_edges_;

 public:
  //
  // CONSTRUCTOR AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_(), edges_(), adjacency_map_(), num_edges_(0)
  {
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
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return (graph_->nodes_[uid_]).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Value of node */
    node_value_type& value(){
      return (graph_->nodes_[uid_]).value;
    }

    /** Value of node as const */
    const node_value_type& value() const{
      return (graph_->nodes_[uid_]).value;
    }

    /**
     * Return the number of nodes that are connected to this node via an edge.
     */
    size_type degree() const{
      if(graph_->adjacency_map_.find(uid_) != graph_->adjacency_map_.end()){
        return (size_type) (graph_->adjacency_map_[uid_]).size();
      }else{
        return 0;
      }
    }

    /**
     * Return an IncidentIterator to an edge with node1 equal to this node
     * @post Returned IncidentIterator iit satisfies
     *  (*iit).node1() == (*this).index() and
     *  and has_edge((*iit).node1(), (*iit).node2())
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,
                              uid_,
                              (graph_->adjacency_map_[uid_]).begin());
    }

    /**
     * Return an IncidentIterator to one past the end of nodes incident
     * to this node
     * @post Returned IncidentIterator iit cannot be dereferenced
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,
                              uid_,
                              (graph_->adjacency_map_[uid_]).end());
    }

    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any geometric meaning.
     *
     * The node ordering relation obeys trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return (graph_ == n.graph_ && uid_ < n.uid_) || (graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph<V>;

    // Pointer to graph that contains node
    Graph* graph_;

    // Index of node in graph
    size_type uid_;

    // Private Constructor
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type) nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    nodes_.push_back({position, value});
    return Node{this, size()-1};
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node{this, i};
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1_ == e.node1_ && node2_ == e.node2_) ||
          (node1_ == e.node2_ && node2_ == e.node1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any interpretive meaning.
     *
     * The edges are ordered first by comparing node1_, then node2_.
     */
    bool operator<(const Edge& e) const {
      return (node1_ < e.node1_) || (node1_ == e.node1_ && node2_ < e.node2_);
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph<V>;

    // Pointer to graph that contains edge
    Graph* graph_;

    // Nodes of the edge
    size_type node1_;
    size_type node2_;

    // Private constructor
    Edge(const Graph* graph, const size_type a, const size_type b) :
      graph_(const_cast<Graph*>(graph)), node1_(a), node2_(b) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  Edge edge(size_type i) const {
    return Edge{this, edges_[i].node1, edges_[i].node2};
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    return (adjacency_map_.find(a.index()) != adjacency_map_.end()
              && adjacency_map_.at(a.index()).find(b.index())
                 != adjacency_map_.at(a.index()).end());
  }

  /** Add an edge to the graph, or return the current edge if it already exists
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b) {

    // Edge (a, b) does not exist, add to Graph
    if(!has_edge(a, b)){
      adjacency_map_[a.index()].insert(b.index());
      adjacency_map_[b.index()].insert(a.index());
      edges_.push_back({a.index(), b.index()});
      num_edges_ += 1;
    }

    return Edge{this, a.index(), b.index()};
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    adjacency_map_.clear();
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator :private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator () {
    }

  /** Return node pointed to by NodeIterator
   * @pre *this != node.end()
   */
  Node operator*() const{
    return graph_->node(idx_);
  }

  /** Increment NodeIterator
   * @pre *this != node.end()
   */
  node_iterator& operator++(){
    idx_++;
    return *this;
  }

  /**
   * Return true if this NodeIterator points to the same
   * node as NodeIterator nit and returns false otherwise
   * @param nit Second NodeIterator for comparison
   * @return true if **this == *nit, false otherwise
   */
  bool operator==(const node_iterator& nit) const{
    return (nit.idx_ == idx_);
  }

 private:
  friend class Graph<V>;

  // Pointer to graph that contains iterator
  Graph* graph_;

  // Index of node that is currently iterated in graph
  size_type idx_;

  // Private Constructor
  NodeIterator(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
  }
};

/**
 * Return an iterator to the first node
 * @post Returned NodeIterator nit satisfies has_node(*nit)
 */
node_iterator node_begin() const{
  return NodeIterator(this, 0);
}

/**
 * Return an iterator to one past the end of all nodes
 * @post Returned NodeIterator nit cannot be dereferenced
 */
node_iterator node_end() const{
  return NodeIterator(this, num_nodes());
}

//
// Incident Iterator
//

/** @class Graph::IncidentIterator
 * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator :private equality_comparable<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /**
     * Return Edge object where node1 is the node being iterated
     * @pre *this != N.edge_end() where N is the root node
     * @return Edge e where e.node1().index() == N and has_edge(e.node1(), e.node2())
     */
    Edge operator*() const{
      return Edge{graph_, node_, *iter_};
    }

    /**
     * Increment this IncidentIterator
     * @pre *this != n.edge_end() where @a n is the root node
     */
    IncidentIterator& operator++(){
      iter_++;
      return *this;
    }

    /**
     * Return true if this iterator and iit point to the same edge
     * @param iit second IncidentIterator for comparison
     * @return true if **this == *iit and false otherwise
     */
    bool operator==(const IncidentIterator& iit) const{
      return node_ == iit.node_ && iter_ == iit.iter_;
    }

   private:
    friend class Graph<V>;

    // Pointer to graph that contains iterator
    Graph* graph_;

    // Origin node
    size_type node_;

    // Pointer to element in adjacency_map_[node_]
    std::unordered_set<size_type>::iterator iter_;

    // Private Constructor
    IncidentIterator(const Graph* graph,
                     size_type node,
                     std::unordered_set<size_type>::iterator iter)
        : graph_(const_cast<Graph*>(graph)),
          node_(node),
          iter_(iter)
    {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator :private equality_comparable<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /**
     * Return Edge that is referenced by iterator
     * @pre *this != edge_end()
     * @return Edge @a e where has_edge(e.node1(), e.node2())
     */
    Edge operator*() const{
      return graph_->edge(idx_);
    }

    /**
     * Increment this EdgeIterator
     * @pre *this != edge_end()
     * @return Reference to EdgeIterator eit where
     *  has_edge((*eit).node1(), (*eit).node2()) or eit == edge_end()
     */
    EdgeIterator& operator++(){
      ++idx_;
      return *this;
    }

    /**
     * Compares two EdgeIterators
     * @param eit second EdgeIterator reference for comparison
     * @return true if **this == *eit, false otherwise
     */
    bool operator==(const EdgeIterator& eit) const{
      return idx_ == eit.idx_;
    }

   private:
    friend class Graph<V>;

    // Pointer to graph that contains iterator
    Graph* graph_;

    // Index of edge that is currently iterated in graph
    size_type idx_;

    // Private Constructor
    EdgeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }

  };

  /**
    * Return an iterator to the first edge
    * @post Returned EdgeIterator eit satisfies
    *   has_edge((*eit).node1(), (*eit).node2())
    */
  edge_iterator edge_begin() const{
    return EdgeIterator{this, 0};
  }

  /**
   * Return an iterator to one past the end of all edges
   * @post Returned EdgeIterator eit cannot be dereferenced
   */
  edge_iterator edge_end() const {
    return EdgeIterator{this, num_edges()};
  }

};

#endif // CME212_GRAPH_HPP

