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
template <typename V>
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
    Point point;
    node_value_type value;
    size_type degree;
  };

  /* Edge struct data container. */
  struct InternalEdge {
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

    /** Return this node's position. */
    const Point& position() const { return graph_->nodes_[index_].point; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return index_; }

    /** Return a reference to this node's value. This is a setter function for node value. */
    node_value_type& value() { InternalNode& node = graph_->nodes_[index_]; return node.value; }

    /** Return a const reference to this node's value. This is a getter function for node value.  */
    const node_value_type& value() const { InternalNode& node = graph_->nodes_[index_]; return node.value; }

    /** Return this node's degree. */
    size_type degree() const { nodes_[index_].degree; }

    /** Return IncidentIterator pointing to the first element in the unordered map container of incident edges to this node. */
    incident_iterator edge_begin() const { return IncidentIterator(graph_, index_, graph_->adj_map_[index_].begin()); }

    /** Return IncidentIterator pointing to the past-the-end element in the unordered map container of incident edges to this node.  */
    incident_iterator edge_end() const { return IncidentIterator(graph_, index_, graph_->adj_map_[index_].end()); }

    /** Test whether this node and @a n are equal. Equal nodes have the same graph and the same index. */
    bool operator==(const Node& n) const { return (index_ == n.index_) && (graph_ == n.graph_); }

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
     *    - this.graph_ == n.graph_ && this.index_ < n.index_
     */
    bool operator<(const Node& n) const {
      auto g = reinterpret_cast<std::uintptr_t>(graph_);
      auto n_g = reinterpret_cast<std::uintptr_t>(n.graph_);
      return (g < n_g) || (graph_ == n.graph_ && index_ < n.index_);
    }


   private:
    friend class Graph;         // Allow Graph to access Node's private member data and functions.
    
    // Node private members
    graph_type* graph_;         // Pointer back to the Graph container - 8 bytes in 64 bit.
    size_type index_;           // This Node's index                   - 4 bytes as unsigned integer.

    // Private constructor accessible by friend class Graph.
    Node(const graph_type* graph, size_type index) : graph_(const_cast<graph_type*>(graph)), index_(index) {}
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
    InternalNode internal_node;
    internal_node.point = position;
    internal_node.value = value;
    internal_node.degree = 0;

    nodes_.emplace(num_nodes_, internal_node);
    ++num_nodes_;
    
    return Node(this, num_nodes_-1);
  }

  /** Determine if a Node belongs to this Graph.
   * @param[in] n   Node to search.
   * @return  True if @a n is currently a Node of this Graph
   *
   * True if same graph and valid index.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const { return (n.graph_ == this && n.index_ < num_nodes_); }

  /** Return the node with index @a i.
   * @param[in] i   Node index.
   * @return  Constructed Node.
   *
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, i); }


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
    Node node1() const { return graph_->node(graph_->edges_[index_].node1_index); }

    /** Return the other node of this Edge */
    Node node2() const { return graph_->node(graph_->edges_[index_].node2_index); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Equal if they contain the same graph and same index.
     */
    bool operator==(const Edge& e) const { return (graph_ == e.graph_ && index_ == e.index_); }

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
      return (g < n_g) || (graph_ == e.graph_ && index_ < e.index_);
    }


   private:
    friend class Graph;       // Allow Graph to access Edge's private member data and functions.
    
    // Edge private members
    graph_type* graph_;       // Pointer back to the Graph container - 8 bytes in 64 bit.
    size_type index_;         // This edge's index                   - 4 bytes as unsigned integer.

    // Private constructor accessible by friend class Graph.
    Edge(const graph_type* graph, size_type index) : graph_(const_cast<graph_type*>(graph)), index_(index) {}
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
   * Complexity: O(1).
   */
  Edge edge(size_type i) const { return Edge(this, i); }

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
    auto adj_nodes_to_a = adj_map_.find(a.index());
    if (adj_nodes_to_a == adj_map_.end()) {
      return false;
    }
    else {
      auto edge_to_b = (adj_nodes_to_a->second).find(b.index());
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

    // get node indices
    size_type index_a = a.index();
    size_type index_b = b.index();

    // pointers to check if adjacency map contains nodes a, b
    auto adj_nodes_to_a = adj_map_.find(index_a);
    auto adj_nodes_to_b = adj_map_.find(index_b);

    // conditions on existence of nodes in adjacency maps, with appropriate insertions
    if (adj_nodes_to_a != adj_map_.end()) {
      auto edge_to_b = (adj_nodes_to_a->second).find(index_b);
      if (edge_to_b != (adj_nodes_to_a->second).end()) {
        return edge(edge_to_b->second);
      } else {
        (adj_nodes_to_a->second).insert({index_b, num_edges_});

        if (adj_nodes_to_b != adj_map_.end()) {
          (adj_nodes_to_b->second).insert({index_a, num_edges_});
        } else {
          std::unordered_map<size_type,size_type> adj_map_b = {{index_a, num_edges_}};
          adj_map_.insert({index_b, adj_map_b});
        }
      }
    } else {
      std::unordered_map<size_type,size_type> adj_map_a = {{index_b, num_edges_}};
      adj_map_.insert({index_a, adj_map_a});

      if (adj_nodes_to_b != adj_map_.end()) {
        (adj_nodes_to_b->second).insert({index_a, num_edges_});
      } else {
        std::unordered_map<size_type,size_type> adj_map_b = {{index_a, num_edges_}};
        adj_map_.insert({index_b, adj_map_b});
      }
    }

    // Add InternalEdge object with given data to edges_ vector
    InternalEdge internal_edge;
    internal_edge.node1_index = index_a; internal_edge.node2_index = index_b;
    edges_.emplace(num_edges_, internal_edge);

    ++num_edges_;

    // Increase degree of node a and node b
    nodes_[index_a].degree++;
    nodes_[index_b].degree++;

    return Edge(this, num_edges_-1);
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
    edges_.clear();
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
    Node operator*() const { return Node(graph_, iter_->first); }
    
    /** Pre-increment operator. */
    NodeIterator& operator++() { ++iter_; return *this; }

    /** Equality operator. */
    bool operator==(const NodeIterator& iit) const { return iter_ == iit.iter_; }

   private:
    friend class Graph;         // Allow Graph to access NodeIterator's private member data and functions. 

    // NodeIterator private data members.
    graph_type* graph_;                                                           // Pointer to graph.  
    typename std::unordered_map<size_type,InternalNode>::const_iterator iter_;    // Iterator for node_ unordered map.

    // Private constructor accessible by friend class Graph.
    NodeIterator(const graph_type* graph, typename std::unordered_map<size_type,InternalNode>::const_iterator iter) : graph_(const_cast<graph_type*>(graph)), iter_(iter) {}
  };

  /** Return NodeIterator pointing to the first element in the unordered map container of internal nodes. */
  node_iterator node_begin() const { return NodeIterator(this, nodes_.begin()); }
  
  /** Return NodeIterator pointing to the past-the-end element in the unordered map container of internal nodes. */
  node_iterator node_end() const { return NodeIterator(this, nodes_.end()); }


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
    Edge operator*() const {
      size_type edge_index = iter_->second;
      InternalEdge& edge = graph_->edges_[edge_index];
      edge.node1_index = node_index_;
      edge.node2_index = iter_->first;
      return Edge(graph_, iter_->second);
    }
    
    /** Pre-increment operator. */
    IncidentIterator& operator++() { ++iter_; return *this; }

    /** Equality operator. */
    bool operator==(const IncidentIterator& iit) const { return iter_ == iit.iter_; }

   private:
    friend class Graph;         // Allow Graph to access IncidentIterator's private member data and functions.

    // IncidentIterator private data members
    graph_type* graph_;                                                       // Pointer to graph.   
    size_type node_index_;                                                    // Index of outgoing node.
    typename std::unordered_map<size_type,size_type>::const_iterator iter_;   // Iterator for incident edges unordered map.

    // Private constructor accessible by friend class Graph.
    IncidentIterator(const graph_type* graph, size_type node_index, typename std::unordered_map<size_type,size_type>::const_iterator iter)
     : graph_(const_cast<graph_type*>(graph)), node_index_(node_index), iter_(iter) { }

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
    Edge operator*() const { return Edge(graph_, iter_->first); }

    /** Pre-increment operator. */
    EdgeIterator& operator++() { ++iter_; return *this; }

    /** Equality operator. */
    bool operator==(const EdgeIterator& iit) const { return iter_ == iit.iter_; }

   private:
    friend class Graph;         // Allow Graph to access EdgeIterator's private member data and functions.

    graph_type* graph_;                                                           // Pointer to graph.
    typename std::unordered_map<size_type,InternalEdge>::const_iterator iter_;    // Iterator for edges_ unordered map.

    // Private constructor accessible by friend class Graph
    EdgeIterator(const graph_type* graph, typename std::unordered_map<size_type,InternalEdge>::const_iterator iter)
      : graph_(const_cast<graph_type*>(graph)), iter_(iter) { }
  };

  /** Return EdgeIterator pointing to the first element in the unordered map container of internal edges. */
  edge_iterator edge_begin() const { return EdgeIterator(this, edges_.begin()); }
  
  /** Return EdgeIterator pointing to the past-the-end element in the unordered map container of internal edges. */
  edge_iterator edge_end() const { return EdgeIterator(this, edges_.end()); }


 private:

  // Private Graph data members.
  size_type num_nodes_;                                                             // Number of nodes in graph.
  size_type num_edges_;                                                             // Number of edges in graph.
  std::unordered_map<size_type,InternalNode> nodes_;                                // Unordered_map of internal nodes with following structure - {node_index: InternalNode}.
  std::unordered_map<size_type,InternalEdge> edges_;                                // Unordered_map of internal edges with following structure - {edge_index: InternalEdge}.
  std::unordered_map<size_type,std::unordered_map<size_type,size_type>> adj_map_;   // Adjacency map with following structure - {node_index: {adjacent_node_index: edge_index}}.
};

#endif // CME212_GRAPH_HPP
