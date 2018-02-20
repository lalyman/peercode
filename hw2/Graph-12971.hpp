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

  // Predeclare the internal structs
  struct internal_node;
  struct internal_edge;
  struct internal_adj;  

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

    /** Return this node's position, read-only. */
    const Point& position() const {
      return fetch_node().position;
    }

    /** Return this node's position and possibly modify it. */
    Point& position() {
      return fetch_node().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

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
      auto pos = graph_->adj_.find(index_);
      if (pos == graph_->adj_.end()) {
        return 0;
      }
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

    /** Return this edge's index. */
    size_type index() const {
      return index_;
    }

    /** Return the length of this Edge. */
    double length() const {
      return norm(graph_->node(nid1_).position() - graph_->node(nid2_).position());
    }

    /** Return this edge's value, read-only. */
    const edge_value_type& value() const {
      return graph_->edges_[index_].value;
    }

    /** Return this edge's value and possibly modify it. */
    edge_value_type& value() {
      return graph_->edges_[index_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges belong to the same graph and represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->graph_ == e.graph_) && (this->nid1_ == e.nid1_ && this->nid2_ == e.nid2_) || (this->nid1_ == e.nid2_ && this->nid2_ == e.nid1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // first compare graph indices, then compare edge indices
      return (this->graph_ < e.graph_) || (this->graph_ == e.graph_ && this->index_ < e.index_);
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
        if ((*it).nid2 == nid2) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    size_type nid1 = a.index();
    size_type nid2 = b.index();
    if (has_edge(a, b)){
      auto it = adj_.at(nid1).begin();
      for (; it != adj_.at(nid1).end(); ++it) {
        if ((*it).nid2 == nid2) {
          return Edge(this, (*it).index, nid1, nid2); // valid edge index
        }
      } 
    } else {
      size_type index = num_edges();  
      
      internal_edge new_edge {nid1, nid2, value};
      edges_.push_back(new_edge);      

      internal_adj new_adj1 {nid2, index};
      adj_[nid1].push_back(new_adj1); // update adjacency list: n1 -> n2
      internal_adj new_adj2 {nid1, index};
      adj_[nid2].push_back(new_adj2); // update adjacency list: n2 -> n1
      
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

  /** Remove a node in the graph.
   * @param[in] n The node that we want to remove.
   * @return The index of the removed node if @a n exists in the graph and has been  
   *         removed successfully. Return size() if @a n does not exist in the graph. 
   * @post new size() == old size() - 1 if @a n has been removed successfully.
   *       new size() == old size() if @a n does not exist in the graph.
   *       The last node in the container has been moved to the position of the removed node.
   *       The following outstanding objects might be invalidated by this function
   *       1. Node with the same index as the removed node;
   *       2. Node with the same index as the last node in the container;   
   *       3. Edge with the removed node as one of its endpoint;
   *       4. The graph's NodeIterator and the IncidentIterator of the removed node.  
   * Complexity: O(d^2) < O(num_nodes()), 
   *             where d is the degree of the node and assume the graph is sparse.   
   */ 
  size_type remove_node(const Node& n) {
   
    // Check whether the node exists.
    if(!has_node(n))
      return size();

    // Remove all edges incident to this node.
    auto ei = n.edge_begin();
    while (ei != n.edge_end()) {
      remove_edge(*ei);
    }

    size_type idx = n.index(); // index of the removed node: e.g. 3 
    // Move the last node to the current position.
    swap_helper(nodes_, idx);  

    size_type old_idx = size(); // index of the last node in the container: e.g. 9
   
    // If no edges associated with the moved node, DONE. No need to update its adjacency list.
    auto pos = adj_.find(old_idx); 
    if (pos == adj_.end()) {
      return idx;
    }      
   
    auto it = adj_.at(old_idx).begin();
    for (; it != adj_.at(old_idx).end(); ++it) {
        size_type idx2 = (*it).nid2; // 6
        size_type edge_idx = (*it).index; 

        // Add new adjacency for the moved node.
        internal_adj new_adj1 {idx2, edge_idx};
        adj_[idx].push_back(new_adj1); // add 3->6 
        internal_adj new_adj2 {idx, edge_idx};
        adj_[idx2].push_back(new_adj2); // add 6->3          

        // Remove old adjacency for the moved node.
        auto it2 = adj_.at(idx2).begin();
        for (; it2 != adj_.at(idx2).end(); ++it2) {
          if ((*it2).nid2 == old_idx) {
            adj_.at(idx2).erase(it2); // erase 6->9
            break;
          } 
        }

        // Update the endpoint of associated edge.   
        if (edges_[edge_idx].nid1 == old_idx) // set (9,6) to (3,6) or (6,9) to (6,3)
          edges_[edge_idx].nid1 = idx;
        else
          edges_[edge_idx].nid2 = idx; 
    }
    // Remove old adjacency for the moved node.
    adj_.erase(old_idx); // erase 9->6
    
    return idx;
  }

  /** Remove a node in the graph.
   * @param[in] it The iterator pointing to the node that we want to remove.
   * @return The iterator pointing to the position of the removed node.
   *         It is supposed to point to the node that we moved from the end of the container. 
   * @post 1. new size() == old size() - 1 if @a *it has been removed successfully.
   *       2. new size() == old size() if @a *it does not exist in the graph.
   *       3. The last node in the container has been moved to the position of removed node.
   *       4. The following outstanding objects might be invalidated by this function
   *       (1) Node with the same index as the removed node;
   *       (2) Node with the same index as the last node in the container;   
   *       (3) Edge with the removed node as one of its endpoint;
   *       (4) The graph's NodeIterator and the IncidentIterator of the removed node.  
   * Complexity: O(d^2) < O(num_nodes()), 
   *             where d is the degree of the node and assume the graph is sparse.   
   */ 
  NodeIterator remove_node(NodeIterator it) {
    remove_node(*it);
    return it;
  }

  /** Remove an edge in the graph.
   * @param[in] a, b The endpoints of the edge that we want to remove.
   * @return The index of the removed edge if @a e exists in the graph and has been  
   *         removed successfully. Return num_edges() if @e n does not exist in the graph. 
   * @post 1. new num_edges() == old num_edges() - 1 if @e n has been removed successfully.
   *       2. new size() == old size() if @a n does not exist in the graph.
   *       3. The last edge in the container has been moved to the position of the removed edge.
   *       4. The following outstanding objects might be invalidated by this function
   *       (1) Edge with the same index as the removed edge;
   *       (2) Edge with the same index as the last edge in the container;   
   *       (3) The graph's EdgeIterator and the IncidentIterator of the removed edge's endpoints.  
   * Complexity: O(d), where d is the degree of the edge's endpoints.   
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type nid1 = a.index();
    size_type nid2 = b.index();
    auto it = adj_.at(nid1).begin();
      for (; it != adj_.at(nid1).end(); ++it) {
        if ((*it).nid2 == nid2) {
          Edge e = Edge(this, (*it).index, nid1, nid2);    
          return remove_edge(e);
        }
      }
  }
  
  /** Remove an edge in the graph.
   * @param[in] e The edge that we want to remove.
   * @return The index of the removed edge if @a e exists in the graph and has been  
   *         removed successfully. Return num_edges() if @e n does not exist in the graph. 
   * @post 1. new num_edges() == old num_edges() - 1 if @e n has been removed successfully.
   *       2. new size() == old size() if @a n does not exist in the graph.
   *       3. The last edge in the container has been moved to the position of removed edge.
   *       4. The following outstanding objects might be invalidated by this function
   *       (1) Edge with the same index as the removed edge;
   *       (2) Edge with the same index as the last edge in the container;   
   *       (3) The graph's EdgeIterator and the IncidentIterator of the removed edge's endpoints.  
   * Complexity: O(d), where d is the degree of the edge's endpoints.   
   */ 
  size_type remove_edge(const Edge& e) {
    size_type idx = e.index(); 
    size_type old_nid1 = e.node1().index(); // e.g. 1
    size_type old_nid2 = e.node2().index(); // e.g. 3

    // Check whether the edge exist.
    if (!has_edge(node(old_nid1), node(old_nid2))) {
      return num_edges();
    }

    // Remove this edge from adjacency list.
    auto it1 = adj_.at(old_nid1).begin();
    for (; it1 != adj_.at(old_nid1).end(); ++it1) {
        if ((*it1).nid2 == old_nid2) {
          adj_.at(old_nid1).erase(it1); // erase 1->3
          break;
        }
    }

    // Remove this edge (w/ reversed direction) from adjacency list. 
    auto it2 = adj_.at(old_nid2).begin();
    for (; it2 != adj_.at(old_nid2).end(); ++it2) {
        if ((*it2).nid2 == old_nid1) {
          adj_.at(old_nid2).erase(it2); // erase 3->1
          break;
        }
    }

    // Move the last edge to the current position.
    swap_helper(edges_, idx);
    
    size_type new_nid1 = edges_[idx].nid1; // e.g 2
    size_type new_nid2 = edges_[idx].nid2; // e.g 4

    // Update edge index in adjacency list for the moved edge.
    auto it3 = adj_.at(new_nid1).begin();
    for (; it3 != adj_.at(new_nid1).end(); ++it3) {
        if ((*it3).nid2 == new_nid2) {
          (*it3).index = idx; // associate (2,4) with new edge_idx 
          break;
        }
    }

    // Update edge index in adjacency list for the moved edge (w/ reserved direction). 
    auto it4 = adj_.at(new_nid2).begin();
    for (; it4 != adj_.at(new_nid2).end(); ++it4) {
        if ((*it4).nid2 == new_nid1) {
          (*it4).index = idx; // associate (4,2) with new edge_idx
          break;
        }
    }

    return idx;
  }

  /** Remove an edge in the graph.
   * @param[in] ei The iterator pointing to th edge that we want to remove.
   * @return The iterator pointing to the position of the removed node. 
   *         It is supposed to point to the edge that we moved from the end of the container. 
   * @post 1. new num_edges() == old num_edges() - 1 if @e n has been removed successfully.
   *       2. new size() == old size() if @a n does not exist in the graph.
   *       3. The last edge in the container has been moved to the position of removed edge.
   *       4. The following outstanding objects might be invalidated by this function
   *       (1) Edge with the same index as the removed edge;
   *       (2) Edge with the same index as the last edge in the container;   
   *       (3) The graph's EdgeIterator and the IncidentIterator of the removed edge's endpoints.  
   * Complexity: O(d), where d is the degree of the edge's endpoints.   
   */ 
  EdgeIterator remove_edge(EdgeIterator ei) {
    remove_edge(*ei);
    return ei;
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
    // Pointer back to the Graph
    Graph* graph_;
    // Current position in the Graph's Node list
    size_type count_;
    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type count) 
        : graph_(const_cast<Graph*>(graph)), count_(count) {
    }
  };

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
    
    /** Dereference operator that returns the Edge at current position. */
    Edge operator*() const {
      size_type nid1 = index_;
      size_type nid2 = graph_->adj_.at(index_).at(count_).nid2;
      //return Edge(graph_, 0, nid1, nid2); // invalidate edge index
      size_type edge_idx = graph_->adj_.at(index_).at(count_).index; // valid edge index
      return Edge(graph_, edge_idx, nid1, nid2); 
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
    // Pointer back to the Graph
    Graph* graph_;
    // Current position in the Graph's Edge list
    size_type count_; 
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type count) 
        : graph_(const_cast<Graph*>(graph)), count_(count) {
    }
    
  };
  
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

  // swap the element to be erased with the last element in the container
  template <typename T>
  bool swap_helper(std::vector<T>& v, const size_type index) {
   bool swapped = true;
   if (index == v.size()-1) 
     swapped = false;
   else
    v[index] = v.back();
   v.pop_back();
   return swapped;
  } 

  struct internal_node {
    Point position; // The node's position
    node_value_type value; // The node's value 
  };

  struct internal_edge {
    size_type nid1, nid2; // Indices of the edge's vertices
    edge_value_type value; // The edge's value
  };

  struct internal_adj {
    size_type nid2; // The index of the adjacent node 
    size_type index; // The edge's index 
  };
  
  // Node container
  std::vector<internal_node> nodes_;
  // Edge container
  std::vector<internal_edge> edges_;
  // Adjacency list
  std::map<size_type, std::vector<internal_adj>> adj_;

};

#endif // CME212_GRAPH_HPP
