#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Type of Nodes */
  using node_value_type = V;
  
  /** Type of Edges */
  using edge_value_type = E;

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

  //** Type of data structure maintaining adjacency of nodes,
  // a nested unordered_map */
  using adjacency_map_type = std::unordered_map<size_type,
      std::unordered_map<size_type, size_type>>;

  //
  // PRIVATE INTERNAL DATA
  //

 private:
  /** Internal structure for holding Node data */
  struct internal_node {
    size_type uid;
    Point position;
    node_value_type value;
  };

  /** Internal structure for holding Edge data */
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    edge_value_type value;
  };

  /** Vector of internal_node structs */
  std::vector<internal_node> nodes_;

  /** Map from node uid's to indexes */
  std::unordered_map<size_type, size_type> node_map_;

  /** Next uid for nodes */
  size_type next_uid_ = 0;

  /** Vector of internal_edge structs */
  std::vector<internal_edge> edges_;

  /** Adjacency map that holds edges with each node uid appearing as a key */
  adjacency_map_type adjacency_map_;

 public:
  //
  // CONSTRUCTOR AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_(), edges_(), adjacency_map_()
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

    /** Return this node's position by reference. */
    Point& position(){
      return (graph_->nodes_[graph_->node_map_[uid_]]).position;
    }

    /** Return this node's position. */
    const Point& position() const {
      return (graph_->nodes_[graph_->node_map_[uid_]]).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->node_map_[uid_];
    }

    /** Value of node */
    node_value_type& value(){
      return (graph_->nodes_[graph_->node_map_[uid_]]).value;
    }

    /** Value of node as const */
    const node_value_type& value() const{
      return (graph_->nodes_[graph_->node_map_[uid_]]).value;
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
     *  (*iit).node1() == *this and
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
     * Equal nodes have the same graph and the same uid.
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
    friend class Graph<V,E>;

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
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type())
  {
    nodes_.push_back({next_uid_, position, value});
    node_map_[next_uid_] = (size_type) nodes_.size()-1;
    return Node{this, next_uid_++};
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (node_map_.find(n.uid_) != node_map_.end());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node{this, nodes_[i].uid};
  }

  /**
   * Remove Node @a n and all edges that contain Node @a n from this graph
   * @param n Node for removal
   * @return 0 if @a n is not a node in this graph
   *    1 if @a n is a node in this graph and removed successfully
   * @post Node @a n is invalidated
   * @post All edges containing Node @a n are invalidated
   * @post All node iterators and edge iterators are invalidated
   * @post Incident iterators with @a n as the root or where @a n
   *    is an adjacent node are invalidated
   * @post Indexes of nodes and edges may not be the same after removal
   * @post old num_edges() == new num_edges() + old n.degree() if return is 1
   *    old num_edges() == new num_edges() if return is 0
   * @post old num_nodes() == new num_nodes() + 1 if return is 1
   *    old num_nodes() == new num_nodes() if return is 0
   *
   * Complexity: O(n.degree())
   */
  size_type remove_node(const Node& n){
    if(!has_node(n))
      return 0;

    // Remove Edges incident to Node n
    while(n.degree() > 0){
      remove_edge(*(n.edge_begin()));
    }

    size_type idx = n.index();
    nodes_[idx] = nodes_[nodes_.size()-1];
    node_map_[nodes_[idx].uid] = idx;
    nodes_.pop_back();
    node_map_.erase(n.uid_);

    return 1;
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
      return Node{graph_, node1_uid_};
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node{graph_, node2_uid_};
    }

    /** Value of edge */
    edge_value_type& value(){
      size_type idx = graph_->adjacency_map_[node1_uid_][node2_uid_];
      return (graph_->edges_[idx]).value;
    }

    /** Value of edge as const */
    const edge_value_type& value() const{
      size_type idx = graph_->adjacency_map_[node1_uid_][node2_uid_];
      return (graph_->edges_[idx]).value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes
     * in the same graph.
     */
    bool operator==(const Edge& e) const {
      return graph_ == e.graph_ &&
          ((node1_uid_ == e.node1_uid_ && node2_uid_ == e.node2_uid_) ||
          (node1_uid_ == e.node2_uid_ && node2_uid_ == e.node1_uid_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any interpretive meaning.
     *
     * The edges are ordered first by comparing node1_, then node2_.
     */
    bool operator<(const Edge& e) const {
      return graph_ < e.graph_ ||
          (graph_ == e.graph_ && ((node1_uid_ < e.node1_uid_) ||
              (node1_uid_ == e.node1_uid_ && node2_uid_ < e.node2_uid_)));
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph<V,E>;

    // Pointer to graph that contains edge
    Graph* graph_;

    // Nodes of the edge
    size_type node1_uid_;
    size_type node2_uid_;

    // Private constructor
    Edge(const Graph* graph, const size_type a, const size_type b) :
      graph_(const_cast<Graph*>(graph)), node1_uid_(a), node2_uid_(b) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return (size_type) edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_nodes())
   */
  Edge edge(size_type i) const {
    return Edge{this, edges_[i].node1_uid, edges_[i].node2_uid};
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    return (adjacency_map_.find(a.uid_) != adjacency_map_.end()
              && adjacency_map_.at(a.uid_).find(b.uid_)
                 != adjacency_map_.at(a.uid_).end());
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
      adjacency_map_[a.uid_][b.uid_] = num_edges();
      adjacency_map_[b.uid_][a.uid_] = num_edges();
      edges_.push_back({a.uid_, b.uid_,0});
    }

    return Edge{this, a.uid_, b.uid_};
  }

  /**
   * Remove Edge with nodes n1, n2 from this graph
   * @param n1 First node that defines the edge for removal
   * @param n2 Second node that defines the edge for removal
   * @return 0 if (@a n1, @a n2) is not an edge in this graph
   *    1 if (@a n1, @a n2) is an edge in this graph and removed successfully
   * @post Edge defined by (@a n1, @a n2) is invalidated
   * @post All edge iterators are invalidated
   * @post Incident iterators with @a n1 or @a n2 as the root are invalidated
   * @post Indexes of edges may not be the same after removal
   * @post old num_edges() == new num_edges() + 1 if return is 1
   *  old num_edges() == new num_edges() if return is 0
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& n1, const Node& n2){
    if(!has_node(n1) || !has_node(n2) || !has_edge(n1, n2))
      return 0;

    size_type e_idx = adjacency_map_[n1.uid_][n2.uid_];

    edges_[e_idx] = edges_[num_edges()-1];
    adjacency_map_[edges_[e_idx].node1_uid][edges_[e_idx].node2_uid] = e_idx;
    adjacency_map_[edges_[e_idx].node2_uid][edges_[e_idx].node1_uid] = e_idx;

    edges_.pop_back();
    adjacency_map_[n1.uid_].erase(n2.uid_);
    adjacency_map_[n2.uid_].erase(n1.uid_);

    return 1;
  }

  /**
   * Remove Edge @a e from this graph
   * @param e Edge for removal from this graph
   * @return 0 if @a e is not an edge in this graph
   *    1 if @a e is an edge in this graph and removed successfully
   * @post Edge @a e is invalidated
   * @post All edge iterators are invalidated
   * @post Incident iterators with e.node1() or e.node2() as
   *    the root are invalidated
   * @post Indexes of edges may not be the same after removal
   * @post old num_edges() == new num_edges() + 1 if return is 1
   *  old num_edges() == new num_edges() if return is 0
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge& e){
    if(!has_node(e.node1()) ||
        !has_node(e.node2()) ||
        !has_edge(e.node1(), e.node2()))
      return 0;

    size_type n1_uid = e.node1_uid_;
    size_type n2_uid = e.node2_uid_;
    size_type e_idx = adjacency_map_[n1_uid][n2_uid];

    edges_[e_idx] = edges_[num_edges()-1];
    adjacency_map_[edges_[e_idx].node1_uid][edges_[e_idx].node2_uid] = e_idx;
    adjacency_map_[edges_[e_idx].node2_uid][edges_[e_idx].node1_uid] = e_idx;

    edges_.pop_back();
    adjacency_map_[n1_uid].erase(n2_uid);
    adjacency_map_[n2_uid].erase(n1_uid);

    return 1;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    node_map_.clear();
    edges_.clear();
    adjacency_map_.clear();
    next_uid_ = 0;
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
    bool operator==(const node_iterator& nit) const {
      return (nit.idx_ == idx_);
    }

  private:
    friend class Graph<V,E>;

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
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  /**
   * Remove Node *(@a nit) and all edges that contain Node *(@a nit) from this graph
   * @param nit NodeIterator pointing to node for removal
   * @return NodeIterator pointing to valid node or node_end()
   * @post Node *(@a nit) is invalidated
   * @post All edges containing Node *(@a nit) are invalidated
   * @post All node iterators and edge iterators are invalidated
   * @post Incident iterators with @a n as the root or where @a n
   *    is an adjacent node are invalidated
   * @post Indexes of nodes and edges may not be the same after removal
   * @post old num_edges() == new num_edges() + old (*nit).degree()
   * @post old num_nodes() == new num_nodes() + 1
   *
   * Complexity: O((*nit).degree())
   */
  node_iterator remove_node(node_iterator nit){
    Node n = *nit;
    remove_node(n);
    return nit;
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
      return Edge{graph_, node_, iter_->first};
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
    friend class Graph<V,E>;

    // Pointer to graph that contains iterator
    Graph* graph_;

    // Origin node
    size_type node_;

    // Pointer to element in adjacency_map_[node_]
    std::unordered_map<size_type, size_type>::iterator iter_;

    // Private Constructor
    IncidentIterator(const Graph* graph,
                     size_type node,
                     std::unordered_map<size_type, size_type>::iterator iter)
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
    friend class Graph<V,E>;

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

  /**
   * Remove Edge *(@a eit) from this graph
   * @param e Edge for removal from this graph
   * @return EdgeIterator pointing to a valid edge in this graph or edge_end()
   * @post Edge *(@a eit) is invalidated
   * @post All edge iterators are invalidated
   * @post Incident iterators with (*(@a eit)).node1() or(*(@a eit)).node2()
   *    as the root are invalidated
   * @post Indexes of edges may not be the same after removal
   * @post old num_edges() == new num_edges() + 1
   *
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator eit){
    Edge e = *eit;
    remove_edge(e);
    return eit;
  }

};

#endif // CME212_GRAPH_HPP

