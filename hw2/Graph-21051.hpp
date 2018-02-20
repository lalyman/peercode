#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 *  @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <cstdint>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. Edges are unique (there is at
 *  most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
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
  /** Template Node value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Template Edge value type. */
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

  /** Node struct data container. */
  struct InternalNode {
    Point position;
    node_value_type value;
    size_type degree;
    size_type index;
  };

  /* Edge struct data container. */
  struct InternalEdge {
    edge_value_type value;
    size_type index;
    size_type node1_index;
    size_type node2_index;
  };


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : num_nodes_(0), num_edges_(0) {}

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

    /** Getter for this node's position. */
    const Point& position() const { return graph_->nodes_[uid_].position; }

    /** Setter for this node's position. */
    Point& position() { return graph_->nodes_[uid_].position; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { InternalNode n = graph_->nodes_[uid_]; return n.index; }

    /** Return a reference to this node's value. This is a setter function for node value. */
    node_value_type& value() { InternalNode& node = graph_->nodes_[uid_]; return node.value; }

    /** Return a const reference to this node's value. This is a getter function for node value.  */
    const node_value_type& value() const { InternalNode& node = graph_->nodes_[uid_]; return node.value; }

    /** Return this node's degree. */
    size_type degree() const { nodes_[uid_].degree; }

    /** Return IncidentIterator pointing to the first element in the unordered map container of incident edges to this node. */
    incident_iterator edge_begin() const { return IncidentIterator(graph_, uid_, graph_->adj_map_[uid_].begin()); }

    /** Return IncidentIterator pointing to the past-the-end element in the unordered map container of incident edges to this node.  */
    incident_iterator edge_end() const { return IncidentIterator(graph_, uid_, graph_->adj_map_[uid_].end()); }

    /** Test whether this node and @a n are equal. Equal nodes have the same graph and the same index. */
    bool operator==(const Node& n) const { return (uid_ == n.uid_) && (graph_ == n.graph_); }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     * Compare graph_ pointers as casted integer types. If pointers match, compare indices.
     * That is, less than if one of the two conditions holds
     *    - this.graph_ < n.graph_
     *    - this.graph_ == n.graph_ && this.uid_ < n.uid_
     */
    bool operator<(const Node& n) const {
      auto g = reinterpret_cast<std::uintptr_t>(graph_);
      auto n_g = reinterpret_cast<std::uintptr_t>(n.graph_);
      return (g < n_g) || (graph_ == n.graph_ && uid_ < n.uid_);
    }


   private:
    friend class Graph;         // Allow Graph to access Node's private member data and functions.
    
    // Node private members
    graph_type* graph_;         // Pointer back to the Graph container - 8 bytes in 64 bit.
    size_type uid_;             // This Node's uid_a                   - 4 bytes as unsigned integer.

    // Private constructor accessible by friend class Graph.
    Node(const graph_type* graph, size_type uid) : graph_(const_cast<graph_type*>(graph)), uid_(uid) {}   
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return num_nodes_; }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position  The new node's position.
   * @param[in] value     The new node's value.
   * @return  A Node proxy to input node data.
   *
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.position() == @a position
   * @post result_node.value() == @a value
   *  
   * Create an InternalNode object with given data and add to nodes_ unordered map
   *
   * Complexity: O(1).
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    InternalNode internal_node{position, value, 0, num_nodes_};
    i2u_.push_back(nodes_.size());
    nodes_.push_back(internal_node);
    ++num_nodes_;
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph.
   * @param[in] n   Node to search.
   * @return  True if @a n is currently a Node of this Graph
   *
   * True if same graph and valid index.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const { return (n.graph_ == this && n.index() < num_nodes_); }

  /** Return the node with index @a i.
   * @param[in] i   Node index.
   * @return  Constructed Node.
   *
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { 
    size_type uid = i2u_[i];
    return Node(this, uid);
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

    /** Return a node of this Edge */
    Node node1() const { return Node(graph_, node1_uid_); }

    /** Return the other node of this Edge */
    Node node2() const { return Node(graph_, node2_uid_); }

    /** Return a reference to this edge's value. This is a setter function for edge value. */
    edge_value_type& value() { InternalEdge& edge = graph_->adj_map_[node1_uid_][node2_uid_]; return edge.value; }

    /** Return a const reference to this edge's value. This is a getter function for edge value.  */
    const edge_value_type& value() const { InternalEdge& edge = graph_->adj_map_[node1_uid_][node2_uid_]; return edge.value; }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Equal if they contain the same graph and same index.
     */
    bool operator==(const Edge& e) const { return (graph_ == e.graph_ && node1_uid_ == e.node1_uid_ && node2_uid_ == e.node2_uid_); }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * Compare graph_ pointers as casted integer type. If pointers match, compare indices.
     * That is, less than if one of the two conditions holds
     *    - this.graph_ < n.graph_
     *    - this.graph_ == n.graph_ && this.index_ < n.index_
     */
    bool operator<(const Edge& e) const {
      auto g = reinterpret_cast<std::uintptr_t>(graph_);
      auto n_g = reinterpret_cast<std::uintptr_t>(e.graph_);
      return (g < n_g) ||
             (graph_ == e.graph_ && node1_uid_ < e.node1_uid_) ||
             (graph_ == e.graph_ && node1_uid_ == e.node1_uid_ && node2_uid_ < e.node2_uid_);
    }


   private:
    friend class Graph;       // Allow Graph to access Edge's private member data and functions.
    
    // Edge private members
    graph_type* graph_;       // Pointer back to the Graph container - 8 bytes in 64 bit.
    size_type node1_uid_;     // Node 1 uid                          - 4 bytes in 64 bit.
    size_type node2_uid_;     // Node 2 uid                          - 4 bytes in 64 bit.

    // Private constructor accessible by friend class Graph.
    Edge(const graph_type* graph, size_type node1_uid, size_type node2_uid)
      : graph_(const_cast<graph_type*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const { return num_edges_; }

  /** Return the edge with index @a i.
   * @param[in] i   Edge index.
   * @return  Constructed Edge.
   *
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges()).
   */
  Edge edge(size_type i) const {
    for (auto it = edge_begin(); it != edge_end(); ++it) {
      auto edge = *it;
      auto node1_uid = i2u_[edge.node1().index()];
      auto node2_uid = i2u_[edge.node2().index()];
      auto map = adj_map_.at(node1_uid);
      if (map.at(node2_uid).index == i)
        return Edge(this, node1_uid, node2_uid);
    }
    return Edge();
}

  /** Test whether two nodes are connected by an edge.
   * @param[in] a   Node 1;
   * @param[in] b   Node 2;
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * @pre @a a and @a b are valid nodes of this graph 
   *
   * Complexity: O(1) average, O(num_edges()) worst case.
   */
  bool has_edge(const Node& a, const Node& b) const {
    /** Check if adjacency map contains node @a a
          if not, given edge is not in graph
          if yes, return if node @a b is in adjacency map of node @a a
    */
    size_type uid_a = i2u_[a.index()];
    size_type uid_b = i2u_[b.index()];
    auto adj_nodes_to_a = adj_map_.find(uid_a);
    if (adj_nodes_to_a == adj_map_.end()) {
      return false;
    }
    else {
      auto edge_to_b = (adj_nodes_to_a->second).find(uid_b);
      return (edge_to_b != (adj_nodes_to_a->second).end());
    }
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @param[in] a   Node 1;
   * @param[in] b   Node 2;
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(1) average, O(num_edges()) worst case.
   */
  Edge add_edge(const Node& a, const Node& b) {
    /* Check if adjacency map contains node @a a
         if yes, check if node @a b is in adjacency map of node @a a
             - if yes, return edge index contained in map
              - if no, insert node @a b and edge index in adjacency map of node @a a
        else,
            insert node @a a in adjacency map as key and adjacency map of node @a a as value, initialized to contain edge to @a b 
            check if adjacency map contains node @a b
              - if yes, insert node @a a and edge index in adjacency map of @a b
              - if no, insert node @a b in adjacency map as key and adjacency map of node @a b as value, initialized to to contain edge to @a a
    */

    // get node uid's
    size_type uid_a = i2u_[a.index()];
    size_type uid_b = i2u_[b.index()];

    // if edge already in adjacency map, return edge
    if (has_edge(Node(this, uid_a), Node(this, uid_b))) 
      return Edge(this, uid_a, uid_b);

    // form InternalEdge element
    InternalEdge edge{edge_value_type(), num_edges_, uid_a, uid_b};

    // pointers to check if adjacency map contains nodes a, b
    auto adj_nodes_to_a = adj_map_.find(uid_a);
    auto adj_nodes_to_b = adj_map_.find(uid_b);

    // Insert edge in adjacency map of a
    if (adj_nodes_to_a != adj_map_.end()) {
      (adj_nodes_to_a->second).insert({uid_b, edge});
    } else {
      std::unordered_map<size_type,InternalEdge> adj_map_a = {{uid_b, edge}};
      adj_map_.insert({uid_a, adj_map_a});
    }

    // Insert edge in adjacency map of b
    if (adj_nodes_to_b != adj_map_.end()) {
        (adj_nodes_to_b->second).insert({uid_a, edge});
    } else {
        std::unordered_map<size_type,InternalEdge> adj_map_b = {{uid_a, edge}};
        adj_map_.insert({uid_b, adj_map_b});
    }

    // Increase number of edges
    ++num_edges_;

    // Increase degree of node a and node b
    nodes_[uid_a].degree++;
    nodes_[uid_b].degree++;

    return Edge(this, uid_a, uid_b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_nodes_ = 0;
    num_edges_ = 0;
    nodes_.clear();
    i2u_.clear();
    adj_map_.clear();
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
    NodeIterator() { }

    /** Dereference operator. Return constructed Node.  */
    Node operator*() const { return Node(graph_, *iter_); }
    
    /** Pre-increment operator. */
    NodeIterator& operator++() { ++iter_; return *this; }

    /** Equality operator. */
    bool operator==(const NodeIterator& iit) const { return iter_ == iit.iter_; }

   private:
    friend class Graph;         // Allow Graph to access NodeIterator's private member data and functions. 

    // NodeIterator private data members.
    graph_type* graph_;                                         // Pointer to graph.  
    typename std::vector<size_type>::const_iterator iter_;      // Iterator for i2u_ vector.

    // Private constructor accessible by friend class Graph.
    NodeIterator(const graph_type* graph, typename std::vector<size_type>::const_iterator iter) : graph_(const_cast<graph_type*>(graph)), iter_(iter) {}
  };

  /** Return NodeIterator pointing to the first element in i2u_ vector. */
  node_iterator node_begin() const { return NodeIterator(this, i2u_.begin()); }
  
  /** Return NodeIterator pointing to the past-the-end element in i2u_ vector. */
  node_iterator node_end() const { return NodeIterator(this, i2u_.end()); }


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
    IncidentIterator() { }

    /** Dereference operator.
     * @return Constructed Edge.
     * @post   (result_edge.node1().index() == node_index_) && (result_edge.node2().index() == adjacent_node_index)
     *
     * Change node indices of edge to ensure node2 is adjacent node.
     */
    Edge operator*() const { return Edge(graph_, node_uid_, iter_->first); }
    
    /** Pre-increment operator. */
    IncidentIterator& operator++() { ++iter_; return *this; }

    /** Equality operator. */
    bool operator==(const IncidentIterator& iit) const { return iter_ == iit.iter_; }

   private:
    friend class Graph;         // Allow Graph to access IncidentIterator's private member data and functions.

    // IncidentIterator private data members
    graph_type* graph_;                                                       // Pointer to graph.   
    size_type node_uid_;                                                      // uid of outgoing node.
    typename std::unordered_map<size_type,InternalEdge>::iterator iter_;      // Iterator for incident edges unordered map.

    // Private constructor accessible by friend class Graph.
    IncidentIterator(const graph_type* graph, size_type node_uid, typename std::unordered_map<size_type,InternalEdge>::iterator iter)
     : graph_(const_cast<graph_type*>(graph)), node_uid_(node_uid), iter_(iter) {}

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
    EdgeIterator() { }

    /** Dereference operator. */
    Edge operator*() const { return Edge(graph_, outer_iter_->first, inner_iter_->first); }

    /** Pre-increment operator. */
    EdgeIterator& operator++() { 
      ++inner_iter_;
      if (inner_iter_ == (outer_iter_->second).end()) {
        ++outer_iter_;
        if (outer_iter_ != graph_->adj_map_.end())
          inner_iter_ = (outer_iter_->second).begin();
      }
      return *this;
    }

    /** Equality operator. */
    bool operator==(const EdgeIterator& iit) const {
      if (outer_iter_ == (graph_->adj_map_).end())
        return (iit.outer_iter_ == (graph_->adj_map_).end());
      return (inner_iter_ == iit.inner_iter_);
    }

   private:
    friend class Graph;         // Allow Graph to access EdgeIterator's private member data and functions.

    graph_type* graph_;                                                                                                 // Pointer to graph.
    typename std::unordered_map<size_type,std::unordered_map<size_type,InternalEdge>>::const_iterator outer_iter_;      // Outer iterator through nodes in adjacency map.
    typename std::unordered_map<size_type,InternalEdge>::const_iterator inner_iter_;                                    // Inner iterator through adjacent nodes for given node in outer iterator.

    // Private constructor accessible by friend class Graph
    EdgeIterator(const graph_type* graph, typename std::unordered_map<size_type,std::unordered_map<size_type,InternalEdge>>::const_iterator outer_iter, typename std::unordered_map<size_type,InternalEdge>::const_iterator inner_iter)
      : graph_(const_cast<graph_type*>(graph)), outer_iter_(outer_iter), inner_iter_(inner_iter) { }
  };

  /** Return EdgeIterator pointing to the first element in adjacency map of internal edges. */
  edge_iterator edge_begin() const { return EdgeIterator(this, adj_map_.begin(), (adj_map_.begin()->second).begin()); }
  
  /** Return EdgeIterator pointing to the past-the-end element in adjaceny map. */
  edge_iterator edge_end() const { return EdgeIterator(this, adj_map_.end(), (adj_map_.begin()->second).end()); }


  //
  //  REMOVE NODES & EDGES
  //

  /** Remove node given Node reference.
   * @param[in] n   Node;
   * @return uid of node that replaced @a n in iud_ vector.
   *
   * @pre 0 <= @a n.index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_nodes() == old num_nodes() - 1.
   * @post new num_edges() == old num_edges() - n.degree().
   * @post new iud_[@a n.index()] == old iud_[back].
   *
   * Remove entry at @a n.index() from iud_ vector and replace with entry at end of vector.
   *
   * Complexity: O(@a n.degree())
   */
  size_type remove_node(const Node& n) {
    auto it = n.edge_begin();
    while (it != n.edge_end())
      it = remove_edge(it);

    size_type n_idx = n.index();
    size_type replacement_uid = i2u_.back();
    i2u_[n_idx] = replacement_uid;
    i2u_.pop_back();
    nodes_[replacement_uid].index = n_idx;
    --num_nodes_;
    
    return replacement_uid;
  }

  /** Remove node given node iterator.
   * @param[in] n_it   NodeIterator;
   * @return NodeIterator of node that replaced that of @a n_it in in iud_ vector.
   *
   * @pre 0 <= @a (*n_it).index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_nodes() == old num_nodes() - 1.
   * @post new num_edges() == old num_edges() - n.degree().
   * @post new iud_[@a (*n_it).index()] == old iud_[back].
   *
   * Remove entry at @a (*n_it).index() from iud_ vector and replace with entry at end of vector.
   *
   * Complexity: O(@a n.degree())
   */
  node_iterator remove_node(node_iterator n_it) {
    const Node& n = *n_it;

    auto it = n.edge_begin();
    while (it != n.edge_end())
      it = remove_edge(it);

    size_type n_idx = n.index();
    size_type replacement_uid = i2u_.back();
    i2u_[n_idx] = replacement_uid;
    i2u_.pop_back();
    nodes_[replacement_uid].index = n_idx;
    --num_nodes_;

    return NodeIterator(this, i2u_.begin()+n_idx);
  }

  /** Remove edge given incident nodes.
   * @param[in] a   Node 1;
   * @param[in] b   Node 2;
   * @return true if edge was removed from adj_map_, false otherwise.
   *
   * @pre 0 <= @a a.index() < num_nodes().
   * @pre 0 <= @a b.index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_edges() == old num_edges() - 1.
   * @post new @a a.degree() = old @a a.degree() - 1.
   * @post new @a b.degree() = old @a b.degree() - 1.
   *
   * Remove edge from adjacency map. Update number for edges and degrees of adjacent nodes.
   *
   * Complexity:  Average O(1), Worst Case O(@a a.degree() + @a b.degree()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type uid_a = i2u_[a.index()];
    size_type uid_b = i2u_[b.index()];
    std::unordered_map<size_type,InternalEdge>& adj_map_a = adj_map_[uid_a];
    std::unordered_map<size_type,InternalEdge>& adj_map_b = adj_map_[uid_b];
    bool succ_a = adj_map_a.erase(uid_b);
    bool succ_b = adj_map_b.erase(uid_a);
    if (succ_a && succ_b) {
      num_edges_--;
      nodes_[uid_a].degree--;
      nodes_[uid_b].degree--;
    }
    return succ_a && succ_b;
  }

  /** Remove edge given Edge reference.
   * @param[in] e   Edge.
   * @return true if edge was removed from adj_map_, false otherwise.
   *
   * @pre 0 <= @a e.node1().index() < num_nodes().
   * @pre 0 <= @a e.node2().index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_edges() == old num_edges() - 1.
   * @post new @a e.node1().degree() = old @a e.node1().degree() - 1.
   * @post new @a e.node2().degree() = old @a e.node2().degree() - 1.
   *
   * Remove edge from adjacency map. Update number for edges and degrees of adjacent nodes.
   *
   * Complexity:  Average O(1), Worst Case O(@a e.node1().degree() + @a e.node2().degree()).
   */
  size_type remove_edge(const Edge& e) {
    size_type uid_a = i2u_[e.node1().index()];
    size_type uid_b = i2u_[e.node2().index()];
    std::unordered_map<size_type,InternalEdge>& adj_map_a = adj_map_[uid_a];
    std::unordered_map<size_type,InternalEdge>& adj_map_b = adj_map_[uid_b];
    bool succ_a = adj_map_a.erase(uid_b);
    bool succ_b = adj_map_b.erase(uid_a);
    if (succ_a && succ_b) {
      num_edges_--;
      nodes_[uid_a].degree--;
      nodes_[uid_b].degree--;
    }
    return succ_a && succ_b;
  }

  /** Remove edge given edge iterator.
   * @param[in] e_it   EdgeIterator.
   * @return EdgeIterator following @a e_it.
   *
   * @pre 0 <= @a (*e_it).node1().index() < num_nodes().
   * @pre 0 <= @a (*e_it).node2().index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_edges() == old num_edges() - 1.
   * @post new @a (*e_it).node1().degree() = old @a (*e_it).node1().degree() - 1.
   * @post new @a (*e_it).node2().degree() = old @a (*e_it).node2().degree() - 1.
   *
   * Remove edge from adjacency map. Update number for edges and degrees of adjacent nodes.
   *
   * Complexity:  Average O(1), Worst Case O(@a e.node1().degree() + @a e.node2().degree()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    edge_iterator next_e_it = ++e_it;
    size_type uid_a = i2u_[e.node1().index()];
    size_type uid_b = i2u_[e.node2().index()];
    std::unordered_map<size_type,InternalEdge>& adj_map_a = adj_map_[uid_a];
    std::unordered_map<size_type,InternalEdge>& adj_map_b = adj_map_[uid_b];
    bool succ_a = adj_map_a.erase(uid_b);
    bool succ_b = adj_map_b.erase(uid_a);
    if (succ_a && succ_b) {
      num_edges_--;
      nodes_[uid_a].degree--;
      nodes_[uid_b].degree--;
    }
    return next_e_it;
  }

  /** Remove edge given incident iterator.
   * @param[in] e_it   IncidentIterator.
   * @return IncidentIterator following @a e_it.
   *
   * @pre 0 <= @a (*e_it).node1().index() < num_nodes().
   * @pre 0 <= @a (*e_it).node2().index() < num_nodes().
   * @pre num_nodes() > 0.
   * @post new num_edges() == old num_edges() - 1.
   * @post new @a (*e_it).node1().degree() = old @a (*e_it).node1().degree() - 1.
   * @post new @a (*e_it).node2().degree() = old @a (*e_it).node2().degree() - 1.
   *
   * Remove edge from adjacency map. Update number for edges and degrees of adjacent nodes.
   *
   * Complexity:  Average O(1), Worst Case O(@a e.node1().degree() + @a e.node2().degree()).
   */
  incident_iterator remove_edge(incident_iterator e_it) {
    Edge e = *e_it;
    incident_iterator next_e_it = ++e_it;
    size_type uid_a = i2u_[e.node1().index()];
    size_type uid_b = i2u_[e.node2().index()];
    std::unordered_map<size_type,InternalEdge>& adj_map_a = adj_map_[uid_a];
    std::unordered_map<size_type,InternalEdge>& adj_map_b = adj_map_[uid_b];
    bool succ_a = adj_map_a.erase(uid_b);
    bool succ_b = adj_map_b.erase(uid_a);
    if (succ_a && succ_b) {
      num_edges_--;
      nodes_[uid_a].degree--;
      nodes_[uid_b].degree--;
    }
    return next_e_it;
  }



 private:

  // Private Graph data members.
  size_type num_nodes_;                                                                 // Number of nodes in graph.
  size_type num_edges_;                                                                 // Number of edges in graph.
  std::vector<InternalNode> nodes_;                                                     // Vector of internal nodes indexed by node uid's.
  std::vector<size_type> i2u_;                                                          // Vector of index-to-uid.
  std::unordered_map<size_type,std::unordered_map<size_type,InternalEdge>> adj_map_;    // Adjacency map with following structure - {node_index: {adjacent_node_index: InternalEdge}}.
};

#endif // CME212_GRAPH_HPP
