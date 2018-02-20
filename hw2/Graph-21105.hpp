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
template <typename V, typename E>
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

  /** Predeclaration of Edge value type */
  using edge_value_type = E;
  
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

  /** Remove a node.
   *
   * @param[in] n The node to be removed
   * @result[in] Removes a node from the graph, and returns the
   *    index of the removed node. Also removes all edges incident to
   *    to the removed node and references to it in other nodes
   *
   * @pre  @a n.index() < num_nodes()
   *
   * @post @a new num_nodes() == old num_nodes()-1
   * @post @a new num_edges() == old num_edges()-n.degree()
   * @post @a old (node.index() == num_nodes-1) == new (node.index() == n.index())
   *    the node formerly containing the last node index, now contains the index of the removed node
   *
   * @post @a old edge.index != new edge.index
   *    edges are re-indexed post edge removal
   *    see remove edge documentation for specifics
   *
   * @post @a old node.fetch_node.edge_inds_ != node.fetch_node.edge_inds_
   *    nodes incident to the re-indexed node have their edges re-indexed 
   *    see remove edge documentation for specifics
   *
   * @post @a new i2u.size() == old i2u.size()-1
   *    unique identifier to invalidated internal node is removed 
   * @post @a new i2e.size() == old i2e.size()-n.degree()
   *    unique identifiers to invalidated internal edges are removed 
   * 
   * Complexity: O(n.degree()).
   *
   */
  size_type remove_node(const Node& n) {
    // determine node index
    size_type index = n.index();

    // remove affiliated edges
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it) remove_edge((*it));

    // swap index to end of i2u
    std::swap(i2u_[index], i2u_[num_nodes()-1]);
    // assign new index value to swapped node identifier
    nodes_[i2u_[index]].index = index;
    // remove node reference from i2u
    i2u_.pop_back();

    // correct edge index pairing to reflect new index assignment
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        if ((*it).node1().index() == num_nodes()) {
          edges_[i2e_[(*it).index()]].index1 = index; 
        }
        else if ((*it).node2().index() == num_nodes()) {
          edges_[i2e_[(*it).index()]].index2 = index; 
        }
    }

    // return index of removed node
    return  index;
  }
    
  /** Remove a node.
   *
   * @param[in] n_it The iterator point at the node to be removed
   * @result[in] Removes a node from the graph, and returns the
   *    index of the removed node. Also removes all edges incident to
   *    to the removed node and references to it in other nodes
   *
   * @pre  @a n.index() < num_nodes()
   *
   * @post @a new num_nodes() == old num_nodes()-1
   * @post @a new num_edges() == old num_edges()-n.degree()
   * @post @a old (node.index() == num_nodes-1) == new (node.index() == n.index())
   *    the node formerly containing the last node index, now contains the index of the removed node
   *
   * @post @a old edge.index != new edge.index
   *    edges are re-indexed post edge removal
   *    see remove edge documentation for specifics
   *
   * @post @a old node.fetch_node.edge_inds_ != node.fetch_node.edge_inds_
   *    nodes incident to the re-indexed node have their edges re-indexed 
   *    see remove edge documentation for specifics
   *
   * @post @a new i2u.size() == old i2u.size()-1
   *    unique identifier to invalidated internal node is removed 
   * @post @a new i2e.size() == old i2e.size()-n.degree()
   *    unique identifiers to invalidated internal edges are removed 
   * 
   * Complexity: O(n.degree()).
   *
   */
  size_type remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return (*n_it).index(); 
  }

  /** Remove an edge.
   *
   * @param[in] e The edge to be removed
   * @result[in] Removes an edge from the graph, and returns the index of the removed edge
   *
   * @pre  @a e.index() < num_edges()
   *
   * @post @a new num_edges() == old num_edges()-1
   * @post @a old (edge.index() == num_edges-1) == new (edge.index() == e.index())
   *    the edge formerly containing the last edge index, now contains the index of the removed edge
   * @post @a new i2e.size() == old i2e.size()-n.degree()
   *    unique identifiers to invalidated internal edges are removed 
   * @post @a old node.fetch_node.edge_inds_ != node.fetch_node.edge_inds_
   *    nodes incident to the edge with the last index have their edge index modified
   *    to reflect the change
   * 
   * 
   * Complexity: O(n1.degree() + n2.degree()).
   *
   */
  size_type remove_edge(const Edge& e) {
    // determine edge index
    size_type index = e.index();

    // obtain indices of nodes in last edge index
    size_type ind1 = edges_[i2e_[num_edges()-1]].index1;
    size_type ind2 = edges_[i2e_[num_edges()-1]].index2;

    // reassign the new edge index to the nodes
    auto it1 = std::find(nodes_[i2u_[ind1]].edge_inds_.begin(), nodes_[i2u_[ind1]].edge_inds_.end(), num_edges()-1);
    auto it2 = std::find(nodes_[i2u_[ind2]].edge_inds_.begin(), nodes_[i2u_[ind2]].edge_inds_.end(), num_edges()-1);
    (*it1) = index;
    (*it2) = index;

    // swap index to end of i2e
    std::swap(i2e_[index], i2e_[num_edges()-1]);
    // assign new index value to swapped index
    edges_[i2e_[index]].index = index;
    // remove edge reference from i2e
    i2e_.pop_back();

    return index;
  }

  /** Remove an edge.
   *
   * @param[in] n1 First node incident to edge to be removed
   * @param[in] n2 Second node incident to edge to be removed
   * @result[in] Removes an edge from the graph, and returns true if edge existed and was removed
   *
   * @pre  @a e.index() < num_edges()
   *
   * @post @a new num_edges() == old num_edges()-1
   * @post @a old (edge.index() == num_edges-1) == new (edge.index() == e.index())
   *    the edge formerly containing the last edge index, now contains the index of the removed edge
   * @post @a new i2e.size() == old i2e.size()-n.degree()
   *    unique identifiers to invalidated internal edges are removed 
   * @post @a old node.fetch_node.edge_inds_ != node.fetch_node.edge_inds_
   *    nodes incident to the edge with the last index have their edge index modified
   *    to reflect the change
   * 
   * Complexity: O(num_edges()).
   *
   */
  bool remove_edge(const Node& n1, const Node& n2) {
    size_type b_index = n1.index();
    size_type a_index = n2.index();
    auto has_it = std::find_if(edge_begin(), edge_end(), [a_index, b_index](decltype(*edge_begin()) e) {
      return ((a_index == e.node1().index() && b_index == e.node2().index()) || (a_index == e.node2().index() && b_index == e.node1().index()));
    });
    if (has_it != edge_end()) {
      remove_edge(*has_it); 
      return true;
    }
    return false;
  }

  /** Remove an edge.
   *
   * @param[in] e_it The iterator pointing at the edge to be removed
   * @result[in] Removes an edge from the graph, and returns the index of the removed edge
   *
   * @pre  @a e.index() < num_edges()
   *
   * @post @a new num_edges() == old num_edges()-1
   * @post @a old (edge.index() == num_edges-1) == new (edge.index() == e.index())
   *    the edge formerly containing the last edge index, now contains the index of the removed edge
   * @post @a new i2e.size() == old i2e.size()-n.degree()
   *    unique identifiers to invalidated internal edges are removed 
   * @post @a old node.fetch_node.edge_inds_ != node.fetch_node.edge_inds_
   *    nodes incident to the edge with the last index have their edge index modified
   *    to reflect the change
   * 
   * Complexity: O(n1.degree() + n2.degree()).
   *
   */
  size_type remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return (*e_it).index();
  }


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

    /** Return this node's position. */
    Point& position() {
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
      std::cout << fetch_node().value << std::endl;
      return fetch_node().value;
    }

    /** Return the number of edges pertaining to a node */
    size_type degree() const {
      return fetch_node().edge_inds_.size();
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
        return graph_ -> nodes_[graph_ -> i2u_[index_]];
    }

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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  // Node add_node(const Point& position) {
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // current number of nodes
    size_type prev_num_nodes = num_nodes();
    // create empty array of edge neighbors
    std::vector<size_type> edge_indices;
    // create new internal node instance
    internal_node new_internal_node {position, prev_num_nodes, node_value, edge_indices};
    // update i2u index to accept new node
    i2u_.push_back(nodes_.size());
    // append new internal node to nodes_ vector
    nodes_.push_back(new_internal_node);
    // obtain instance of new node 
    Node new_node = Node(this, new_internal_node.index);
    
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
    if (i < i2u_.size()) return Node(this, i);
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

    size_type index() const {
        return index_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index2_);
    }

    /** Return type of this edge's value and allow it to be re-assigned */
    edge_value_type& value() {
      return fetch_edge().value;
    }

    /** Return this edge's value (without allowing it to be re-assigned). */
    const edge_value_type& value() const {
      return fetch_edge().value;
    }

    /** Return length of the edge */ 
    double length() const{
      const Point& point = node1().position();
      const Point& point1 = node2().position();
      double dist = sqrt(pow(point.x-point1.x, 2) + pow(point.y-point1.y, 2) + pow(point.z-point1.z, 2));
      return dist;
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
    size_type index_;
    size_type index1_;
    size_type index2_;


    // Private verbatim constructor
    Edge(const Graph* graph, size_type index, size_type index1, size_type index2) : 
        graph_(const_cast<Graph*>(graph)),
        index_(index), index1_(index1), index2_(index2) {}

    // call to internal nodes for proxy structure
    internal_edge& fetch_edge() const {
        return graph_ -> edges_[graph_ -> i2e_[index_]];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2e_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < i2e_.size()) return Edge(this, edges_[i2e_[i]].index, edges_[i2e_[i]].index1, edges_[i2e_[i]].index2); 
    else return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   
    size_type b_index = b.index();
    size_type a_index = a.index();
    auto has_it = std::find_if(edge_begin(), edge_end(), [a_index, b_index](decltype(*edge_begin()) e) {
        return ((a_index == e.node1().index() && b_index == e.node2().index()) || (a_index == e.node2().index() && b_index == e.node1().index()));
    });
    if (has_it != edge_end())
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

  Edge add_edge(const Node& a, const Node& b, const edge_value_type& = edge_value_type()) {
    bool val_has_edge = has_edge(a, b);

    // current number of edges
    size_type prev_num_edges = num_edges();
      
    if (!val_has_edge) {
      internal_edge new_edge {prev_num_edges, a.index(), b.index(), edge_value_type()};
      i2e_.push_back(edges_.size());
      edges_.push_back(new_edge);
      a.fetch_node().edge_inds_.push_back(prev_num_edges);
      b.fetch_node().edge_inds_.push_back(prev_num_edges);
    }
    else {
      for (incident_iterator it = a.edge_begin(); it != a.edge_end(); ++it) {
        if ((*it).node2().index() == b.index()) {
          prev_num_edges = (*it).index();
          break;
        }
      }
    }
    Edge nedge = Edge(this, prev_num_edges, a.index(), b.index());

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
    i2u_.clear();
    i2e_.clear();
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
      size_type ind = graph_ -> nodes_[graph_ -> i2u_[n_index_]].edge_inds_[index_];
      size_type n1_index = graph_ -> edges_[graph_ -> i2e_[ind]].index1;
      size_type n2_index = graph_ -> edges_[graph_ -> i2e_[ind]].index2;
      if (n_index_ == n2_index) return Edge(graph_, ind, n2_index, n1_index );
      return Edge(graph_, ind, n1_index, n2_index ); 
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

    /** De-reference operator for edge iterator- returns Edge object at current index*/
    Edge operator*() const {
      size_type n1_index = graph_ -> edges_[graph_ -> i2e_[index_]].index1;
      size_type n2_index = graph_ -> edges_[graph_ -> i2e_[index_]].index2;
      return Edge(graph_, index_, n1_index, n2_index);
    }

    /** Increments iterator */ 
    edge_iterator& operator++() {
      ++index_ ;
      return *this;
    }

    /** Equivalency Operator for EdgeIterator */ 
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
    std::vector<size_type> edge_inds_;
  };

  struct internal_edge {
    size_type index, index1, index2; 
    edge_value_type value;
  };

  // vector of unique node identifiers
  std::vector<size_type> i2u_;  
  // vector of unique edge identifiers
  std::vector<size_type> i2e_; 

  // nodes
  std::vector<internal_node> nodes_;
  // edges
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
