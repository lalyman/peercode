#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <list>
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
  // Predeclare the internal structs
  struct internal_node;

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
  /** Synonym for template type parameter V. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for template type parameter E. */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes() {
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
      // The default constructor should create an invalid node, so we leave it empty.
    }

    /** Return this node's position. */
    const Point& position() const { return fetch_node().point; }

    Point& position() { return fetch_node().point; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return fetch_node().idx; }

    /** Return this node's value. */
    node_value_type& value() { return fetch_node().nv; }

    /** Return this node's value as const. */
    const node_value_type& value() const { return fetch_node().nv; }

    /** Return the degree of this node */
    size_type degree() const { return fetch_node().deg; }

    /** Return the iterator which points to the first incident edge */
    incident_iterator edge_begin() const {
      // the incident edges are implemented in two lists, edges1 and edges2
      // so if edges1 is empty, beginning iterator should be edges2.begin()
      if (fetch_node().edges1.empty())
        return incident_iterator(graph_, nid_, fetch_node().edges2.begin()); 
      else
        return incident_iterator(graph_, nid_, fetch_node().edges1.begin()); 
    }

    /** Return the iterator which points to one past the last incident edge */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, nid_, fetch_node().edges2.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return nid_ == n.nid_&& graph_ == n.graph_;
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
      return nid_ < n.nid_ || (nid_ == n.nid_ && graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type nid_;

    /** Private Constructor */
    Node(const Graph* graph, size_type nid)
        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }

    internal_node& fetch_node() const {
      assert(nid_ < graph_->nodes.size());
      return graph_->nodes[nid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2n.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] nv (optional) The new node's value. 
   *               If not provided, use the default value is used.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nv = node_value_type()) {
    internal_node new_node {position, nv, num_nodes(), {}, {}, {}, 0};
    i2n.push_back(nodes.size());
    nodes.push_back(new_node);
    return Node(this,nodes.size()-1);
  }

  /** Remove a node from the graph
   * @param[in] n The node to be removed
   * @return 1 if @a n is removed, 0 if the graph does not have @a n.
   * @post new num_nodes() == old num_nodes() - 1
   * @post All edges incident to @a n are removed.
   * @post References to @a n is invalidated.
   * @post The following iterators are invalidated, including the end() iterator:
   *       1) node_iterators at or after @a n
   *       2) edge_iterators at any of the removed edges
   *       3) incident_iterators at any of the removed edges
   *
   * Complexity: O(num_nodes()), assuming the max degree is of O(1), 
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) return 0;

    // remove edges incident to this node.
    for (auto n1 : nodes[n.nid_].edges1) {
      nodes[n1].edges2.remove(n.nid_);
      nodes[n1].edges2_v.erase(n.nid_);
      --(nodes[n1].deg);
    }
    for (auto n2 : nodes[n.nid_].edges2) {
      nodes[n2].edges1.remove(n.nid_);
      --(nodes[n2].deg);
    }
    nodes[n.nid_].edges1.clear();
    nodes[n.nid_].edges2.clear();
    nodes[n.nid_].edges2_v.clear();
    nodes[n.nid_].deg = 0;

    // remove the node by erasing it in i2n
    size_type idx = n.index();
    i2n.erase(i2n.begin() + idx);
    for (size_type i = idx; i != i2n.size(); ++i) 
      nodes[i2n[i]].idx = i;
    return 1;
  }

  /** Remove a node from the graph
   * @param[in] n_it node_iterator to the node to be removed.
   * @return node_iterator following the removed node. 
          If @a n_it refers to the last node, the end() iterator is returned.
   * @pre @a n_it must be dereferenceable. Thus @a n_it != end() iterator.
   * @post new num_nodes() == old num_nodes() - 1
   * @post All edges incident to the removed node are removed.
   * @post The following iterators are invalidated, including the end() iterator:
   *       1) node_iterators at or after @a n
   *       2) edge_iterators at any of the removed edges
   *       3) incident_iterators at any of the removed edges
   *
   * Complexity: O(num_nodes()), assuming the max degree is of O(1), 
   */
  node_iterator remove_node(node_iterator n_it) {
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_ && n.nid_ < nodes.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i2n[i]);
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, nid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, nid2_);
    }

    /** Return length of this edge from the current node positions */
    double length() const {
      return norm(graph_->nodes[nid1_].point - graph_->nodes[nid2_].point);
    }

    /** Return this edge's value. */
    edge_value_type& value() {
      size_type x = std::min(nid1_, nid2_);
      size_type y = std::max(nid1_, nid2_);
      return graph_->nodes[x].edges2_v[y];
    }

    /** Return this edge's value. */
    const edge_value_type& value() const {
      size_type x = std::min(nid1_, nid2_);
      size_type y = std::max(nid1_, nid2_);
      return graph_->nodes[x].edges2_v[y];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) 
              && ((nid1_==e.nid1_ && nid2_==e.nid2_) || (nid1_==e.nid2_ && nid2_==e.nid1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type small = std::min(nid1_, nid2_);
      size_type large = std::max(nid1_, nid2_);
      size_type esmall = std::min(e.nid1_, e.nid2_);
      size_type elarge = std::max(e.nid1_, e.nid2_);
      return (small < esmall) || ((small == esmall) && (large < elarge))
                || ((small == esmall) && (large == elarge) && graph_ < e.graph_); 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type nid1_, nid2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type nid1, size_type nid2)
        : graph_(const_cast<Graph*>(graph)), nid1_(nid1), nid2_(nid2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type count = 0;
    for (size_type idx = 0; idx < num_nodes(); ++idx) 
      count += nodes[i2n[idx]].deg;
    assert(count%2 == 0); // we've counted each edge twice
    return count/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *(std::next(edge_begin(), i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    const Node& x = std::min(a,b);
    const Node& y = std::max(a,b);
    auto& x_edges = nodes[x.nid_].edges2;
    auto it = std::find(x_edges.begin(), x_edges.end(), y.nid_);
    if (it == x_edges.end()) 
      return false;
    else
      return true;
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
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& ev = edge_value_type()) {
    if (!has_edge(a,b)) {
      const Node& x = std::min(a,b);
      const Node& y = std::max(a,b);
      nodes[x.nid_].edges2.push_back(y.nid_);
      ++(nodes[x.nid_].deg);
      nodes[x.nid_].edges2_v[y.nid_] = ev;
      nodes[y.nid_].edges1.push_back(x.nid_);
      ++(nodes[y.nid_].deg);
    }
    return Edge(this, a.nid_, b.nid_);
  }

  /** Remove an edge from the graph
   * @param[in] n1,n2 The nodes of the edge to be removed.
   * @return 1 if the edge is removed, 0 if the graph does not have the edge.
   * @post new num_edges() == old num_edges() - 1
   * @post References to the edge is invalidated.
   * @post edge_iterators and incident_iterators at the removed edge are
   *       invalidated, including the end() iterators.
   *
   * Complexity: O(max degree)
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1,n2)) return 0;

    const Node& x = std::min(n1,n2);
    const Node& y = std::max(n1,n2);
    nodes[x.nid_].edges2.remove(y.nid_);
    nodes[x.nid_].edges2_v.erase(y.nid_); // remove the edge value too.
    --(nodes[x.nid_].deg);
    nodes[y.nid_].edges1.remove(x.nid_);
    --(nodes[y.nid_].deg);
    return 1;
  }

  /** Remove an edge from the graph
   * @param[in] e The edge to be removed.
   * @return 1 if @a e is removed, 0 if the graph does not have @a e.
   * @post new num_edges() == old num_edges() - 1
   * @post References to the edge is invalidated.
   * @post edge_iterators and incident_iterators at @a e are
   *       invalidated, including the end() iterators.
   *
   * Complexity: O(max degree)
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph
   * @param[in] e_it edge_iterator to the edge to be removed.
   * @return edge_iterator following the removed edge. 
          If @a e_it refers to the last edge, the end() iterator is returned.
   * @pre @a e_it must be dereferenceable. Thus @a e_it != edge_end() iterator.
   * @post new num_edges() == old num_edges() - 1
   * @post References to the edge is invalidated.
   * @post edge_iterators and incident_iterators at the removed edge are
   *       invalidated, including the end() iterators.
   *
   * Complexity: O(max degree)
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    ++e_it;
    remove_edge(e);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    i2n.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Return the current node */
    Node operator*() const { return Node(graph_, graph_->i2n[idx_]); }
    
    /** Return the NodeIterator which points to the next node */
    NodeIterator& operator++() { 
      ++idx_;
      return *this;
    }

    /** Test whether this NodeIterator and @a x are equal
    *
    * Equal NodeIterators belong to the same graph and points to the same node. */
    bool operator==(const NodeIterator& x) const { 
      return graph_ == x.graph_ && idx_ == x.idx_;
    }


   private:
    friend class Graph;

    Graph* graph_;
    size_type idx_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type idx) 
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

  /** Return the iterator which points to the first node */
  node_iterator node_begin() const { return node_iterator(this, 0); }

  /** Return the iterator which points to one past the last node */
  node_iterator node_end() const { return node_iterator(this, num_nodes()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>  {
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

    /** Return the current edge. 
    * @post result_edge.node1() == the node that spawns this IncidentIterator.
    * @post result_edge.node2() == the adjacent node.
    */
    Edge operator*() const { return Edge(graph_, nid_, *it_); }

    /** Return the IncidentIterator which points to the next incident edge. */
    IncidentIterator& operator++() {
      ++it_;
      // incident edges are implemented in two list, edges1 and edges2.
      // so we connect edges1.end() to edges2.begin()
      if (it_ == graph_->nodes[nid_].edges1.end()) 
        it_ = graph_->nodes[nid_].edges2.begin();
      return *this;
    }

    /** Test whether this IncidentIterator and @a x are equal
    *
    * Equal IncidentIterators belong to the same graph and points to the same edge. */
    bool operator==(const IncidentIterator& x) const {
      return graph_ == x.graph_ && nid_ == x.nid_ && it_ == x.it_;
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type nid_;
    std::list<size_type>::iterator it_;

    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type nid, std::list<size_type>::iterator it) 
        : graph_(const_cast<Graph*>(graph)), nid_(nid), it_(it) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    /** Return the current edge. */
    Edge operator*() const { return Edge(graph_, graph_->i2n[idx_], *it_); }

    /** Return the EdgeIterator which points to the next edge. */
    EdgeIterator& operator++() {
      ++it_;
      // Only edges2 is used so that we visit each edge exactly once.
      // Connect node(i).edges2.end() to node(i+1).edges2.begin()
      while (true) {
        if (it_ == graph_->nodes[graph_->i2n[idx_]].edges2.end()) {
          ++idx_;
          // if at one past the last node, break.
          if (idx_ == graph_->num_nodes()) break;
          // if not, connect to the next node.
          it_ = graph_->nodes[graph_->i2n[idx_]].edges2.begin();
        } else break;
      }
      return *this;
    }

    /** Test whether this EdgeIterator and @a x are equal
    *
    * Equal EdgeIterators belong to the same graph and points to the same edge. */
    bool operator==(const EdgeIterator& x) const {
      if (graph_ == x.graph_ && idx_ == x.idx_) {
        if (idx_ == graph_->num_nodes()) return true; // for the end iterator
        else return it_ == x.it_;
      } else return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type idx_;
    std::list<size_type>::iterator it_; // for the end iterator, leave it_ undefined.

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type idx, std::list<size_type>::iterator it) 
        : graph_(const_cast<Graph*>(graph)), idx_(idx), it_(it) {}
    /** Private Constructor for the begin and end iterator */
    EdgeIterator(const Graph* graph, size_type idx) 
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {
          // For end iterator, leave it_ undefined.          
          // For others, find the index >= idx of the first node with an edge in edges2
          while (idx_ != graph_->size()) {
            if (graph_->nodes[graph_->i2n[idx_]].edges2.empty()) ++idx_;
            else {
              it_ = graph_->nodes[graph_->i2n[idx_]].edges2.begin();
              break;
            }
          }
        }
  };

  /** Return the iterator which points to the first edge */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  /** Return the iterator which points to one past the last edge */
  edge_iterator edge_end() const {
    return edge_iterator(this, num_nodes());
  }

 private:
  // Internal type for nodes
  struct internal_node {
    Point point;
    node_value_type nv;
    size_type idx; // Index of this node, which can be different from nid of this node.

    // adjacency list for edges.
    std::list<size_type> edges1; // node id lower than this node
    std::list<size_type> edges2; // node id bigger than this node
    std::map<size_type, edge_value_type> edges2_v; // values of edges in edges2
    size_type deg;    // degree. should satisfy deg == edges1.size() + edges1.size()
  };

  std::vector<internal_node> nodes;
  // helper member to map from index to node id, i.e. i2n[idx] = nid 
  std::vector<size_type> i2n;

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

};

#endif // CME212_GRAPH_HPP
