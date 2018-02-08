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

  // predeclaration of internal structs
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Predeclaration of Node value type */
  using node_value_type = V;
  
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
    
    // attribute initialization
    Node() : graph_(nullptr), index_(-1) {
    }

    /** Return this node's position. */
    const Point& position() const {
      if (graph_ == nullptr) throw std::runtime_error("graph_ is nullptr");
      return fetch_node().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return type of this node's value and allow it to be re-assigned */
    node_value_type& value() {
      return fetch_node().value;
    }

    /** Return this node's value (without allowing it to be re-assigned). */
    const node_value_type& value() const {
      return fetch_node().value;
    }

    /** Return the number of edges pertaining to a node */
    size_type degree() const {
      return graph_ -> node_bonds_[index_].size();
    }

    /** begin iterator for incident nodes */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }

    /** end iterator for incident nodes */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && index_ == n.index_;
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
      if (graph_ != n.graph_) return graph_ < n.graph_; 
      return index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // pointer to original Graph object
    Graph* graph_;
    // index of Node in original Graph 
    size_type index_;
    // Private verbatim constructor
    Node(const Graph* graph, size_type index) : 
        graph_(const_cast<Graph*> (graph)), index_(index) {}

    // call to internal nodes for proxy structure
    internal_node& fetch_node() const {
        return graph_ -> nodes_[index_];
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
  // Node add_node(const Point& position) {
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // current number of nodes
    size_type prev_num_nodes = num_nodes();
    // create new internal node instance
    internal_node new_internal_node {position, prev_num_nodes, node_value_type()};
    // append new internal node to nodes_ vector
    nodes_.push_back(new_internal_node);
    // obtain instance of new node 
    Node new_node = Node(this, new_internal_node.index);
    // create empty vector of neighbor (nodes with shared edge) indices
    std::vector<size_type> neighbor_inds;
    // push to neighborlist vector
    node_bonds_.push_back(neighbor_inds);
    
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_ && n.index() < nodes_.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < nodes_.size()) return Node(this, i);
    else return Node();
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
    Edge(){
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) return false;
      if (index1_ == e.index1_ && index2_ == e.index2_) return true;
      if (index1_ == e.index2_ && index2_ == e.index1_) return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
       size_type e_ind1 = e.node1().index();
       size_type e_ind2 = e.node2().index();
       return (index1_ < e_ind1) && (index1_ < e_ind2) && (index2_ < e_ind1) && (index2_ < e_ind2);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type index1_;
    size_type index2_;

    // Private verbatim constructor
    Edge(const Graph* graph, size_type index1, size_type index2) : 
        graph_(const_cast<Graph*>(graph)),
        index1_(index1), index2_(index2) {}
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
    if (i < edges_.size()) return Edge(this, edges_[i].index1, edges_[i].index2);
    else return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    if (std::find(node_bonds_[a.index()].begin(), node_bonds_[a.index()].end(), b.index()) != node_bonds_[a.index()].end())
        return true;
    else
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
    bool val_has_edge = has_edge(a, b);
    
    if (!val_has_edge) {
        internal_edge new_edge {a.index(), b.index()};
        edges_.push_back(new_edge);
        node_bonds_[a.index()].push_back(b.index());
        node_bonds_[b.index()].push_back(a.index());
    }
    Edge nedge = Edge(this, a.index(), b.index());
    // map_edges_.insert(std::make_pair(std::make_pair(a.index(), b.index()), nedge));

    return nedge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    node_bonds_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    /** De-reference operator for node iterator- returns Node object at current index*/
    Node operator*() const {
      return Node(graph_, index_);
    }

    /** Increments iterator */ 
    node_iterator& operator++() {
      ++index_ ;
      return *this;
    }

    /** Equivalency Operator for NodeIterator */ 
    bool operator==(const NodeIterator& iter1) const {
      if (iter1.graph_ == graph_ && iter1.index_ == index_)
        return true;
      else
        return false;
    }

   private:
    friend class Graph;
    // pointer to original Graph object
    Graph* graph_;
    // index of Node in original Graph 
    size_type index_;
    // Private verbatim constructor
    NodeIterator(const Graph* graph, size_type index) : 
      graph_(const_cast<Graph*> (graph)), index_(index) {}
  };

  /** begin iterator for nodes */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** end iterator for nodes */
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
    IncidentIterator() {
    }

    /** De-reference operator for incident node iterator- returns Edgeobject at current index*/
    Edge operator*() const {
      size_type ind = graph_ -> node_bonds_[n_index_][index_];
      return Edge(graph_, n_index_, ind);
    }
   
    /** Increments iterator */ 
    IncidentIterator& operator++() {
      ++index_; 
      return *this;
    }

    /** Equivalency Operator for NodeIterator */ 
    bool operator==(const IncidentIterator& iter1) const {
      if (iter1.n_index_ == n_index_ && iter1.index_ == index_)
        return true;
      else
        return false;
    }

   private:
    friend class Graph;
    // pointer to original graph object
    Graph* graph_;
    // index of node
    size_type n_index_;
    // index of iterator
    size_type index_;
    // private verbatim constructor
    IncidentIterator(const Graph* graph, size_type n_index, size_type index) : 
        graph_(const_cast<Graph*> (graph)), n_index_(n_index), index_(index) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    /** Incident edge iterator- returns Edgeobject at current index*/
    Edge operator*() const {
      size_type ind = graph_ -> node_bonds_[n1_index_][n2_index_];
      return Edge(graph_, n1_index_, ind);
    }

    /** Increments iterator */ 
    EdgeIterator& operator++() {
      // update iterator index
      ++index_;
      // avoid accessing invalid memory addresses
      if (index_ < graph_->num_edges()) {
        // update neighbor list index
        ++n2_index_;
        // move to new node if current neighbor list exhausted
        if (n2_index_ == Node(graph_,n1_index_).degree()) {
          n2_index_=0; 
          ++n1_index_;
        }
        // retrieve second node index
        size_type ind = graph_ -> node_bonds_[n1_index_][n2_index_];

        // skip already visited edges
        while (ind < n1_index_) {
          ++n2_index_;
          // move to new node if current neighbor list exhausted
          if (n2_index_ == Node(graph_,n1_index_).degree()) {
            n2_index_=0; 
            ++n1_index_;
          }
          // update second node index
          ind = graph_ -> node_bonds_[n1_index_][n2_index_];
        }
      }
      return *this;
    }

    /** Equivalency Operator for NodeIterator */ 
    bool operator==(const EdgeIterator& iter1) const {
      if (iter1.graph_ == graph_ && iter1.index_ == index_)
        return true;
      else
        return false;
    }

   private:
    friend class Graph;
    // pointer to original graph object
    Graph* graph_;
    // index of edge
    size_type index_;
    // index of first pertinent node
    size_type n1_index_;
    // index of second pertinent node index in bond list of first node
    size_type n2_index_;
    // private verbatim constructor
    EdgeIterator(const Graph* graph, size_type index) : graph_(const_cast<Graph*> (graph)), index_(index), n1_index_(0), n2_index_(0) {}
  };

  /** begin iterator for edges */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** end iterator for edges */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // proxy pattern implementation
  struct internal_node {
    Point point;
    size_type index;
    node_value_type value;
  };

  struct internal_edge {
    size_type index1, index2;
  };
 
  // Is this ok in terms of memory requirements? (ok to make new copy of Edge?) 
  // std::map<std::pair<int, int>, Edge> map_edges_; 
  std::vector<std::vector<size_type>> node_bonds_;
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
