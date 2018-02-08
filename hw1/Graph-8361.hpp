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
template <typename V>
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

  /** Type of node/edge ID. */
  using id_type = unsigned;

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

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nid_;
    }

    /** Return this node's value, e.g. for modification. */
    node_value_type& value() {
        return g_->nodes_[nid_].v_;
    }

    /** Return this node's value as const. */
    const node_value_type& value() const {
        return g_->nodes_[nid_].v_;
    }

    /** Number of connected nodes, i.e. with an edge shared with this node. */
    size_type degree() const {
      // Adjacency list stores all neighbours' edgeinfo in adjacency_[nid_]
      return g_->adjacency_[nid_].size();
    }

    /** Iterator to first neighbour edge */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, nid_, 0);
    }

    /** Iterator to last neighbouring edge */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, nid_, degree()-1);
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
    Node(const Graph* graph, id_type index)
        : g_(const_cast<Graph*>(graph)), nid_(index) {
    }

    // Pointer to graph
    Graph* g_;
    // Index of the proxy Node relating to the vector in the superclass
    id_type nid_;

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
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    size_type nid_new = num_nodes();
    nodes_.push_back({position,value}); // Add to node to list.
    adjacency_.push_back({}); // Add new empty vector for neighbours.
    return node(nid_new); // Return new last proxy to node.
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.g_;
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

    /** Return (the smallest) node of this Edge */
    Node node1() const {
      // g_->edges_[eid_] is the tuple of node IDs comprising this edge
      return g_->node(g_->edges_[eid_][0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // g_->edges_[eid_] is the tuple of node IDs comprising this edge
      return g_->node(g_->edges_[eid_][1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Ensure same graph, and same ID
      return g_ == e.g_ && eid_ == e.eid_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return eid_ < e.eid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Graph needs a way to construct valid Edge objects
    /** Private Constructor */
    Edge(const Graph* graph, size_type index)
        : g_(const_cast<Graph*>(graph)), eid_(index) {
    }

    // Pointer to graph
    Graph* g_;
    // ID of the proxy Edge (Use edges_ vector to retrieve Node IDs)
    id_type eid_;

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
      if(adjacency_[a.nid_][i].neighbour_nid_ == b.nid_){
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
    size_type eid_new = num_edges();

    // Index into edge_ using a, and search for b
    for(size_type i = 0; i < a.degree(); i++){
      if(adjacency_[a.nid_][i].neighbour_nid_ == b.nid_){
        return Edge(this,adjacency_[a.nid_][i].neighbour_eid_);
      }
    }

    // Update edges list
    edges_.push_back({{a.nid_,b.nid_}});

    // Update _both_ neighbour adjacencies
    adjacency_[a.nid_].push_back({b.nid_, eid_new});
    adjacency_[b.nid_].push_back({a.nid_, eid_new});
    return Edge(this,eid_new);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_.clear();
    edges_.clear();
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
        : g_(const_cast<Graph*>(graph)), current_nid_(index) {
    }

    /** Retrieve current node of the iterator. */
    Node operator*() const {
      return g_->node(current_nid_);
    }

    /** Iterate to next node. */
    NodeIterator& operator++() {
      current_nid_++;
      return *this;
    }

    /** Test whether this iterator and _ni_ are equal.
     *
     * Equal iterators have the same graph and point to same node id.
     */
    bool operator==(const NodeIterator& nit) const {
      // Ensure same graph, and same current ID
      return g_ == nit.g_ && current_nid_ == nit.current_nid_;
    }

   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // ID of current node this iterator points to
    id_type current_nid_;
  };

  /** Iterator to first node. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Iterator to last node. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes() - 1);
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
    IncidentIterator(const Graph* graph, id_type centre, id_type index)
        : g_(const_cast<Graph*>(graph)), centre_nid_(centre),
          neighbour_id_(index) {
    }

    /** Retrieve current edge of the incident iterator.
     * @post Returned edge has nodes ordered, {centre_nid_,neighbour_id_}. */
    Edge operator*() const {
      id_type neighbour_eid =
                  g_->adjacency_[centre_nid_][neighbour_id_].neighbour_eid_;
      if(g_->edges_[neighbour_eid][1] == centre_nid_) {
        // Flip order around in edges_
        id_type temp = g_->edges_[neighbour_eid][1];
        g_->edges_[neighbour_eid][1] = g_->edges_[neighbour_eid][0];
        g_->edges_[neighbour_eid][0] = temp;
      }
      return g_->edge(neighbour_eid);
    }

    /** Iterate to next edge connected to centre node. */
    IncidentIterator& operator++() {
      neighbour_id_++;
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
            && centre_nid_ == iit.centre_nid_
            && neighbour_id_ == iit.neighbour_id_;
    }


   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // ID of centre node
    id_type centre_nid_;
    // ID of neighbour number, \el [0,centre_nid_.degree() )
    id_type neighbour_id_;
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
    EdgeIterator(const Graph* graph, id_type index)
        : g_(const_cast<Graph*>(graph)), current_eid_(index) {
    }

    /** Retrieve the current edge of the iterator. */
    Edge operator*() const {
      return g_->edge(current_eid_);
    }

    /** Iterate to next non-redundant edge in the graph. */
    EdgeIterator& operator++() {
      current_eid_++;
      return *this;
    }

    /** Test whether this iterator and _ei_ are equal.
     *
     * Equal iterators have the same graph & point to the same edge.
     */
    bool operator==(const EdgeIterator& ei) const {
      // Ensure same graph and same current ID
      return g_ == ei.g_ && current_eid_ == ei.current_eid_;
    }


   private:
    friend class Graph;

    // Pointer to graph
    Graph* g_;
    // ID of current edge this iterator points to
    id_type current_eid_;
  };


  /** Iterator to first edge. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Iterator to last edge. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges() - 1);
  }


 private:
  // Graph class's internals: helper functions, data members, and so forth.

  // Struct containing info for node
  struct nodeinfo{
    Point p_;
    node_value_type v_;
  };

  // Struct containing infor for edge
  struct edgeinfo{
    id_type neighbour_nid_;
    id_type neighbour_eid_;
  };

  /** Nodes: Vector of struct nodeinfo
        where the Node ID is just the index into the vector **/
  std::vector<nodeinfo> nodes_;

  /** Adjacency List: Vector of vectors, where edges_[i] is a vector of all
        adjacent nodes to node i. **/
  std::vector<std::vector<edgeinfo>> adjacency_;
  // Edges: Vector holding the indexed pair of node IDs comprising the edge.
  std::vector<std::array<id_type,2>> edges_;

};

#endif // CME212_GRAPH_HPP
