#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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
template <typename V, typename E>
class Graph {
 private:

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
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_value_type = E;

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

 private:

  // Graph class's internals: helper functions, data members, and so forth.

  // Struct containing info for node
  struct nodeinfo{
   Point p_;
   std::vector<size_type> adj_; // Holds edge id for all neighbours
   node_value_type v_;
   size_type idx_;            // Non-unique index into i2u_
  };

  // Struct containing edge node pairs and edge values
  struct edgeinfo{
    size_type nid1_;
    size_type nid2_;
    edge_value_type v_;
  };


  /** Nodes: Vector of struct nodeinfo, indexed by _unique_ NID */
  std::vector<nodeinfo> nodes_;

  /** Edges: Vector of struct edgeinfo, indexed by Edge Index. */
  std::vector<edgeinfo> edges_;

  /** Map from arbitrary index to unique NID.
      Deletion should modify lightweight i2u_ rather than nodes_. */
  std::vector<size_type> i2u_;


 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return g_->nodes_[nid_].p_;
    }

    /** Return modifiable reference to node's position. */
    Point& position() {
      return g_->nodes_[nid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return g_->nodes_[nid_].idx_;
    }

    /** Return this node's value, e.g. for modification. */
    node_value_type& value() {
      return g_->nodes_[nid_].v_;
    }

    /** Return this node's value as const. */
    const node_value_type& value() const {
      return g_->nodes_[nid_].v_;
    }

    /** Return list of edges connected to this node as const. */
    const std::vector<size_type>& adjacency() const {
      return g_->nodes_[nid_].adj_;
    }

    /** Number of connected nodes, i.e. with an edge shared with this node. */
    size_type degree() const {
      // Adjacency list stores all neighbours' adjinfo in adjacency_[nid_]
      return adjacency().size();
    }

    /** Iterator to first neighbour edge */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, g_->nodes_[nid_].idx_, 0);
    }

    /** Iterator to last neighbouring edge */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, g_->nodes_[nid_].idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return g_ == n.g_ && nid_ == n.nid_;
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
      return nid_ < n.nid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Graph needs a way to construct valid Node objects
    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : g_(const_cast<Graph*>(graph)), nid_(graph->i2u_[index]) {
    }

    // Pointer to graph
    Graph* g_;
    // Unique NID of the proxy node
    size_type nid_;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortised operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    size_type idx_new = num_nodes();
    size_type nid_new = nodes_.size();
    nodes_.push_back({position,{},value,idx_new});  // Add to node list.
    i2u_.push_back(nid_new);                        // Add to active node list.
    return node(idx_new); // Return new last proxy to node.
  }

  /** Removes node a and all connected edges from this graph.
   *  @param[in] a Node to be removed.
   *  @return New node index of the element following the deleted node.
   *  @post has_node(a) == false
   *  @post new num_nodes() == old num_nodes() - 1
   *  @post new num_edges() == old num_edges() - a.degree()
   *
   *  Complexity:   O(num_nodes()) amortised operations.
   *  Invalidates:  Nodes with nid_ == a.nid_
   *                Edges with node1() == a || node2() == a
   *                NodeIterators pointing to an invalidated node
   *                EdgeIterators pointing behind e_min,
   *                      where e_min is the smallest invalidated edge.
   *                IncidentIterators of _a_ and any neighbouring node
   */
  size_type remove_node(const Node& a){

    while (a.degree() != 0){
      remove_edge(a.edge_begin()); // Remove all connected edges
    }

    size_type idx_delete = a.index();
    i2u_.erase(i2u_.begin() + idx_delete); // Remove element from i2u_

    for (size_type i = idx_delete; i < num_nodes(); ++i)
      nodes_[i2u_[i]].idx_--; // Decrement all idx_ if idx_ > idx_delete

    return idx_delete; // Return idx_ of new element in the next position
  }

  /** Overloaded function to remove node using NodeIterator.
   *  @param[in] nit NodeIterator whose node is to be removed.
   *  @return New node iterator for the element following the deleted element.
   */
  node_iterator remove_node(node_iterator nit){
    remove_node(*nit);
    return nit; // NodeIterator uses index which remains unchanged.
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *           Must not be remnant deleted node.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.g_  &&  nodes_[i2u_[n.idx_]].idx_ == n.idx_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return 1st node of this Edge */
    Node node1() const {
      // Ordering of node1 and node2 is coordinated in IncidentIterator
      return g_->node(g_->nodes_[g_->edges_[eidx_].nid1_].idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return g_->node(g_->nodes_[g_->edges_[eidx_].nid2_].idx_);
    }

    /** Euclidean distance between nodes of the edge. **/
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return this edge's value, e.g. for modification. */
    edge_value_type& value() {
      return g_->edges_[eidx_].v_;
    }

    /** Return this edge's value as const. */
    const edge_value_type& value() const {
      return g_->edges_[eidx_].v_;
    }

    /** Return the edge index, a number in [0,num_edges) */
    size_type index() const {
      return eidx_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Ensure same graph, and same ID
      return g_ == e.g_ && eidx_ == e.eidx_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (g_ == e.g_)
        return eidx_ < e.eidx_; // Return sensible comparitor within graph
      else
        return g_ < e.g_; // Silly to compare, but lets just use address
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Graph needs a way to construct valid Edge objects
    /** Private Constructor */
    Edge(const Graph* graph, size_type index)
        : g_(const_cast<Graph*>(graph)), eidx_(index) {
    }

    // Pointer to graph
    Graph* g_;
    // Index of the proxy Edge
    size_type eidx_;

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
    return Edge(this,i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre _a_ and _b_ are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Index into edge_ using a, and search for b
    for(size_type i = 0; i < a.degree(); i++){
      if(edges_[a.adjacency()[i]].nid1_ == b.nid_
          || edges_[a.adjacency()[i]].nid2_ == b.nid_){
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
    size_type eidx_new = num_edges();

    // Index into edge_ using a, and search for b
    for(size_type i = 0; i < a.degree(); i++){
      if(edges_[a.adjacency()[i]].nid1_ == b.nid_
          || edges_[a.adjacency()[i]].nid2_ == b.nid_){
        return Edge(this,a.adjacency()[i]);
      }
    }

    // Update edges list
    edges_.push_back({a.nid_,b.nid_,{}});

    // Update _both_ neighbour adjacencies
    nodes_[a.index()].adj_.push_back(eidx_new);
    nodes_[b.index()].adj_.push_back(eidx_new);
    return Edge(this,eidx_new);
  }

  /** Removes Edge e from this graph.
   *  @param[in] e Edge to be removed.
   *  @return New edge index of the element following the deleted edge.
   *  @post has_edge(e.node1(),e.node2()) == false
   *  @post new num_edges() == old num_edges() - 1
   *
   *  Complexity:   Complexity: O(num_edges())
   *  Invalidates:  Edges with eidx_ == e.eidx_
   *                EdgeIterators and IncidentIterators pointing behind _e_
   */
  size_type remove_edge(const Edge& e){

    size_type eidx_delete = e.index();
    edges_.erase(edges_.begin() + eidx_delete); // Remove element from edges_

    // Remove edge from adj_, and decrement the others if eidx_ > eidx_delete
    for (size_type ni = 0; ni < i2u_.size(); ++ni){ // Over active nodes
      nodeinfo& n = nodes_[i2u_[ni]];
      for (size_type ai = 0; ai < n.adj_.size(); ++ai){ // Over adj
        if (n.adj_[ai] == eidx_delete){
          n.adj_.erase(n.adj_.begin() + ai);
          ai--; // Will have to check this one again
        }
        else if (n.adj_[ai] > eidx_delete){
          n.adj_[ai]--;
        }
      }
    }
    return eidx_delete; // Return eidx_ of new edge in the next position
  }

  /** Overloaded function to remove edge connected by Nodes a & b.
   *  @param[in] a,b Nodes whose edge is to be removed.
   *  @return New edge index of the element following the deleted edge.
   */
  size_type remove_edge(const Node& a, const Node& b){
    Edge e = add_edge(a,b); // Get the current edge
    return remove_edge(e);
  }

  /** Overloaded function to remove edge pointed to by EdgeIterator
   *  @param[in] eit EdgeIterator to the edge to be removed.
   *  @return New EdgeIterator to the element following the deleted edge.
   */
  edge_iterator remove_edge(edge_iterator eit){
    remove_edge(*eit);
    return eit; // EdgeIterator uses index which remains unchanged.
  }

  /** Overloaded function to remove edge pointed to by IncidentIterator
   *  @param[in] eit IncidentIterator to the edge to be removed.
   *  @return New IncidentIterator to the element following the deleted edge.
   */
  incident_iterator remove_edge(incident_iterator eit){
    remove_edge(*eit);
    return eit; // IncidentIterator uses index which remains unchanged.
  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Public Constructor */
    NodeIterator(const Graph* graph, size_type index)
        : g_(const_cast<Graph*>(graph)), current_idx_(index) {
    }

    /** Retrieve current node of the iterator. */
    Node operator*() const {
      return g_->node(current_idx_);
    }

    /** Iterate to next active node. */
    NodeIterator& operator++() {
      current_idx_++;
      return *this;
    }

    /** Test whether this iterator and _ni_ are equal.
     *
     * Equal iterators have the same graph and point to same node index.
     */
    bool operator==(const NodeIterator& nit) const {
      // Ensure same graph, and same current index
      return g_ == nit.g_ && current_idx_ == nit.current_idx_;
    }

   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // Index of current node this iterator points to
    size_type current_idx_;
  };

  /** Iterator to first node. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Iterator to last node. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Public Constructor */
    IncidentIterator(const Graph* graph, size_type centre, size_type index)
        : g_(const_cast<Graph*>(graph)), centre_idx_(centre),
          neighbour_nidx_(index) {
    }

    /** Retrieve current edge of the incident iterator.
     * @post Returned edge has nodes ordered, {centre_nid_,neighbour_nid_}. */
    Edge operator*() const {
      size_type current_eid = g_->nodes_[g_->i2u_[centre_idx_]].adj_[neighbour_nidx_];
      if(g_->edges_[current_eid].nid2_ == g_->i2u_[centre_idx_]) {
        // Flip order around in edges_
        size_type temp = g_->edges_[current_eid].nid2_;
        g_->edges_[current_eid].nid2_ = g_->edges_[current_eid].nid1_;
        g_->edges_[current_eid].nid1_ = temp;
      }
      return g_->edge(current_eid);
    }

    /** Iterate to next edge connected to centre node. */
    IncidentIterator& operator++() {
      neighbour_nidx_++;
      return *this;
    }

    /** Test whether this iterator and _ii_ are equal.
     *
     * Equal iterators have the same graph & centre,
     * and point to the same edge.
     */
    bool operator==(const IncidentIterator& iit) const {
      // Ensure same graph, centre, and same current ID
      return   g_ == iit.g_
            && centre_idx_ == iit.centre_idx_
            && neighbour_nidx_ == iit.neighbour_nidx_;
    }


   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // Index of centre node
    size_type centre_idx_;
    // Index of current neighbour node within adj_, \el [0,centre.degree)
    size_type neighbour_nidx_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Public Constructor */
    EdgeIterator(const Graph* graph, size_type index)
        : g_(const_cast<Graph*>(graph)), current_eidx_(index) {
    }

    /** Retrieve the current edge of the iterator. */
    Edge operator*() const {
      return g_->edge(current_eidx_);
    }

    /** Iterate to next non-redundant edge in the graph. */
    EdgeIterator& operator++() {
      current_eidx_++;
      return *this;
    }

    /** Test whether this iterator and _ei_ are equal.
     *
     * Equal iterators have the same graph & point to the same edge.
     */
    bool operator==(const EdgeIterator& ei) const {
      // Ensure same graph and same current ID
      return g_ == ei.g_ && current_eidx_ == ei.current_eidx_;
    }


   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // ID of current edge this iterator points to
    size_type current_eidx_;
  };


  /** Iterator to first edge. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Iterator to last edge. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

};

#endif // CME212_GRAPH_HPP
