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

  // Predeclare the internal structs
  struct internal_node;
  struct internal_edge;
  
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_(), edges_() {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch_node().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value, read-only. */
    const node_value_type& value() const {
      return fetch_node().value;
    }

    /** Return this node's value and possibly modify it. */
    node_value_type& value() {
      return fetch_node().value;
    }

    /** Return the number of edges incident to this node. */
    size_type degree() const {
      return graph_->adj_.at(index_).size();
    }

    /** Return an iterator pointing to the first incident edge to this node. */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }

    /** Return an iterator referring to the pass-the-end position in this node's IncidentEdge list. */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->graph_ == n.graph_ && this->index_ == n.index_;
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
      // first compare graph indices, then compare node indices
      return (this->graph_ < n.graph_) || (this->graph_ == n.graph_ && this->index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph
    Graph* graph_;
    // This node's index number
    size_type index_;
    /** Private Constructor */
    Node(const Graph* graph, size_type index) 
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
    /** Helper method to return the appropriate node. */
    internal_node& fetch_node() const {
      return graph_->nodes_[index_];
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type nid = num_nodes();
    internal_node new_node {position, value};
    nodes_.push_back(new_node);
    return Node(this, nid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(nid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(nid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->nid1_ == e.nid1_ && this->nid2_ == e.nid2_) || (this->nid1_ == e.nid2_ && this->nid2 == e.nid1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // first compare graph indices, then compare edge indices
      return (this->graph_ < e.graph) || (this->graph_ == e.graph_ && this->index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Allow IncidentIterator to access Edge's private constructor.
    friend class IncidentIterator;
    // Pointer back to the Graph
    Graph* graph_;
    // This edge's index number
    size_type index_;
    // Index number of one vertex
    size_type nid1_;
    // Index number of the other vertex
    size_type nid2_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type index, size_type nid1, size_type nid2)
        : graph_(const_cast<Graph*>(graph)), index_(index), nid1_(nid1), nid2_(nid2) {
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
    return Edge(this, i, edges_[i].nid1, edges_[i].nid2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type nid1 = a.index();
    size_type nid2 = b.index();
    auto pos = adj_.find(nid1);
    if (pos != adj_.end()) {
      auto it = adj_.at(nid1).begin();
      for (; it != adj_.at(nid1).end(); ++it) {
        if (*it == nid2) {
          return true;
        }
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
    size_type nid1 = a.index();
    size_type nid2 = b.index();
    if (has_edge(a, b)){
      return Edge(this, 0, nid1, nid2); // invalidate edge index
    } else {
      size_type index = num_edges();
      internal_edge new_edge {nid1, nid2};
      edges_.push_back(new_edge);
      adj_[nid1].push_back(nid2); // update adjacency list: n1 -> n2
      adj_[nid2].push_back(nid1); // update adjacency list: n2 -> n1
      return Edge(this, index, nid1, nid2);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear(); 
    edges_.clear();
    adj_.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator that returns the Node at current position. */
    Node operator*() const {
      return graph_->node(count_);
    }

    /** Increment operator that advances position in the Graph's Node list. */
    NodeIterator& operator++() {
      count_ += 1;
      return *this;
    }

    /** Equality operator that tests whether two NodeIterators are equal. 
     * 
     * Equal NodeIterators belong to the same Graph and have the same current position.
     */
    bool operator==(const NodeIterator& x) const {
      return graph_ == x.graph_ && count_ == x.count_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Pointer back to the Graph
    Graph* graph_;
    // Current position in the Graph's Node list
    size_type count_;
    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type count) 
        : graph_(const_cast<Graph*>(graph)), count_(count) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator pointing to the first Node in the Graph. */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  /** Return an iterator referring to the pass-the-end position in the Graph's Node list. */
  NodeIterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator that returns the Edge at current position. */
    Edge operator*() const {
      size_type nid1 = index_;
      size_type nid2 = graph_->adj_.at(index_).at(count_);
      return Edge(graph_, 0, nid1, nid2); // invalidate edge index 
    }

    /** Increment operator that advances position in the Node's IncidentEdge list. */
    IncidentIterator& operator++() {
      count_ += 1;
      return *this;
    }

    /** Equality operator that tests whether two IncidentIterators are equal. 
     * 
     * Equal IncidentIterators belong to the same Graph and have the same Node index and current position.
     */
    bool operator==(const IncidentIterator& x) const {
      return graph_ == x.graph_ && index_ == x.index_ && count_ == x.count_;
    }

   private:
    // Allow Graph to access IncidentIterator's private member data and functions. 
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Allow Node to access IncidentIterator's private constructor.
    friend class Node;
    // Pointer back to the Graph
    Graph* graph_;
    // Index of the Node
    size_type index_;
    // Current position in the Node's IncidentEdge list
    size_type count_; 
    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type index, size_type count) 
        : graph_(const_cast<Graph*>(graph)), index_(index), count_(count) {
    }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator that returns the Edge at current position. */
    Edge operator*() const {
      return graph_->edge(count_);
    }

    /** Increment operator that advances position in the Graph's Edge list. */
    EdgeIterator& operator++() {
      count_ += 1;
      return *this;
    }

    /** Equality operator that tests whether two EdgeIterators are equal. 
     * 
     * Equal EdgeIterators belong to the same Graph and have the same current position.
     */
    bool operator==(const EdgeIterator& x) const {
      return graph_ == x.graph_ && count_ == x.count_; 
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer back to the Graph
    Graph* graph_;
    // Current position in the Graph's Edge list
    size_type count_; 
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type count) 
        : graph_(const_cast<Graph*>(graph)), count_(count) {
    }
    
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an iterator pointing to the first Edge in the Graph. */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /** Return an iterator referring to the pass-the-end position in the Graph's Edge list. */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    Point point; // The node's position
    node_value_type value; // The node's value 
  };

  struct internal_edge {
    size_type nid1, nid2; // Indices of the edge's vertices
  };
  
  // Node container
  std::vector<internal_node> nodes_;
  // Edge container
  std::vector<internal_edge> edges_;
  // Adjacency list
  std::map<size_type, std::vector<size_type>> adj_;

};

#endif // CME212_GRAPH_HPP
