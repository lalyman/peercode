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
template <typename V>
class Graph {
 private:

 public:

  /** Type of this graph. */
  using graph_type = Graph;
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
  using uid_type = unsigned int;

  struct nodeinfo {
    Point p_;
    node_value_type v_;
  };

  /** Construct an empty graph. */
  Graph() {

  }

  /** Default destructor */
  ~Graph() = default;

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

    }

    /* a proxy for some Point data in the node_points list */
    Node(size_type uid, const Graph *g) {
      g_ = g;
      uid_ = uid;
    }

    /** Return this node's position. */
    const Point& position() const {
      return g_->node_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /* Obtain the value associated with the node uid_ */
    node_value_type& value() {
      return const_cast<Graph*>(g_)->node_[uid_].v_;
    }

    /* Obtain the value associated with the node uid_ */
    const node_value_type& value() const {
      return g_->node_[uid_].v_;
    }

    /* Number of nodes connected to uid_ */
    size_type degree() const {
      return g_->adj_[uid_].size();
    }

    /* Starting edge to iterate from will be beginning edge 0 */
    incident_iterator edge_begin() const {
      return IncidentIterator(uid_, 0, g_);
    }

    /* The final edge to iterate till using total number of
     * nodes connected to teh nonexistent end edge
     */
    incident_iterator edge_end() const {
      return IncidentIterator(uid_, degree(), g_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return uid_ == n.uid_ && g_ == n.g_;
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
      return uid_ < n.uid_ && g_ == n.g_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    const Graph *g_;
    uid_type uid_;
    std::vector<uid_type> neighbors() {return g_->adj_[uid_];}
 };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_.size();
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
    node_.push_back(nodeinfo{position, value});
    adj_.push_back(std::vector<uid_type>());
    return Node(num_nodes() - 1, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.index() < num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    return Node(i, this);
  }

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(const Node& uid1, const Node& uid2, const Graph *g) {
      uid1_ = uid1.index();
      uid2_ = uid2.index();
      g_ = g;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(uid1_, g_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(uid2_, g_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(g_ != e.g_) return false;
      return (uid1_ == e.uid1_ && uid2_ == e.uid2_)
        || (uid1_ == e.uid2_ && uid2_ == e.uid1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(g_ != e.g_) return g_ < e.g_;
      if (uid1_ == e.uid1_) {
        return uid2_ < e.uid1_;
      } else {
        return uid1_ < e.uid1_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const Graph* g_;
    uid_type uid1_;
    uid_type uid2_;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type sum = 0;
    for(auto it = adj_.begin(); it != adj_.end(); ++it)
      sum += (*it).size();
    return sum/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0 && i < num_edges());
    EdgeIterator eit = edge_begin();
    for(size_type idx = 0; idx < i; ++idx) {
      ++eit;
    }
    return *eit;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    uid_type uid_small = std::min(a.index(), b.index());
    uid_type uid_large = std::max(a.index(), b.index());
    auto adj_small = adj_[uid_small];
    return std::find(adj_small.begin(), adj_small.end(), uid_large) != adj_small.end();
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
    if (has_edge(a,b)) return Edge(a,b, this);
    adj_[a.index()].push_back(b.index());
    adj_[b.index()].push_back(a.index());
    Edge e(a,b, this);
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_.clear();
    for (auto it = adj_.begin(); it != adj_.end(); ++it)
      (*it).clear();
    adj_.clear();
  }

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
    NodeIterator(uid_type uid, const Graph* g) {
      uid_ = uid;
      g_ = g;
    }

    /* Return an actual node */
    Node operator*() const {
      return Node(uid_, g_);
    }

    /** Increment uid_ and return a pointer to itself that we are
     *  deferencing from.
     */
    NodeIterator& operator++() {
      ++uid_;
      return *this;
    }

    /** Check if two node iterators are the same. They are the
     *  same if they have the same parent graph and same index.
     */
    bool operator==(const NodeIterator& other) const {
      return uid_ == other.uid_ && g_ == other.g_;
    }

   private:
    friend class Graph;
    uid_type uid_;
    const Graph* g_;
  };

  /* The first starting point from where to begin iteration */
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }

  /* Last node where iteration will end */
  node_iterator node_end() const {
    return NodeIterator(num_nodes(), this);
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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

    IncidentIterator(uid_type uid1, size_type neighbor_idx, const Graph* g) {
      uid1_ = uid1;
      neighbor_idx_ = neighbor_idx; // retreive a particular neighbor of uid1
      g_ = g;
    }

    /* Return edge connecting node1 and node2 */
    Edge operator*() const {
      uid_type uid2 = g_->adj_[uid1_][neighbor_idx_];
      return Edge(Node(uid1_, g_), Node(uid2, g_), g_);
    }

    /* Increment the neighbors of node uid1 one by one */
    IncidentIterator& operator++() {
      ++neighbor_idx_;
      return *this;
    }

    /* Check if two edges are neighbors, within the graph and connected */
    bool operator==(const IncidentIterator& other) const {
      return neighbor_idx_ == other.neighbor_idx_ && g_ == other.g_ && uid1_ == other.uid1_;
    }

   private:
    friend class Graph;
    uid_type uid1_;
    size_type neighbor_idx_;
    const Graph* g_;
  };

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

    EdgeIterator(uid_type uid1, size_type neighbor_idx, const Graph* g) {
      uid1_ = uid1;
      neighbor_idx_ = neighbor_idx;
      g_ = g;
    }

    /* Obtain edge connected by 2 nodes within the graph */
    Edge operator*() const {
      Node n1 = Node(uid1_, g_);
      Node n2 = Node(g_->adj_[uid1_][neighbor_idx_], g_);
      return Edge(n1, n2, g_);
    }

    /* Increment neighbor_idx_ by one and if the last neighbor of uid_
     * is reached, move to the right in next node of all nodes. Then
     * get to the top of the next node neightbor list.
     */
    EdgeIterator& operator++() {
      ++neighbor_idx_;
      if(neighbor_idx_ == g_->adj_[uid1_].size()) {
        ++uid1_;
        neighbor_idx_ = 0;
      }
      return *this;
    }

    bool operator==(const EdgeIterator& other) const {
      return uid1_ == other.uid1_ && neighbor_idx_ == other.neighbor_idx_ &&  g_ == other.g_;
    }

   private:
    friend class Graph;
    uid_type uid1_;
    size_type neighbor_idx_;
    const Graph* g_;
  };

  /* The first edge in the graph to start iterating from */
  edge_iterator edge_begin() const {
    EdgeIterator it(0,0,this);
    return it;
  }

  /* The last edge in the graph to end iteration */
  edge_iterator edge_end() const {
    EdgeIterator it(num_nodes(),0,this);
    return it;
  }

 private:

  /* a node containing the point in a struct */
  std::vector<nodeinfo> node_;
  /*  Vector of vectors where each node in outer vector represents node_i.
   *  The vector at that location contains the uid's of all the nodes it
   *  connects to via edges
   */
  std::vector<std::vector<uid_type>> adj_;
};

#endif // CME212_GRAPH_HPP
