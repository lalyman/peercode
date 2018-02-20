#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>

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
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of the node values. */
  using node_value_type = V;

  /** Type of the edge values. */
  using edge_value_type = E;

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
    Node() {}

    /** Return a constant reference to this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_].position;
    }

    /** Return a modifiable reference to this node's position. */
    Point& position() {
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx;
    }

    /**
    * @brief Return a reference to this node's value
    *
    * @return a reference to an object of type node_value_type
    */
    node_value_type& value() {
      return graph_->nodes_[uid_].value;
    }

    /**
    * @brief Return a reference to this node's value
    *
    * @return a reference to a constant object of type node_value_type
    * @post this will not be modified
    */
    const node_value_type& value() const {
      return const_cast<node_value_type& >(value());
    }

    /**
    * @brief Return this node's degree (number of incident edges)
    *
    * @return the number of edges incident to the current node in its graph
    * @post this will not be modified
    */
    size_type degree() const {
      return graph_->edges_[uid_].size();
    }

    /**
    * @brief Return an iterator for the position belonging to the first edge incident to this node
    *
    * @return an iterator of type incident_iterator
    * @post this will not be modified
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /**
    * @brief Return an iterator for the position one past the final edge incident to this node
    *
    * @return an iterator of type incident_iterator
    * @post this will not be modified
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph_ == graph_ && n.index() == index();
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
      if (n == *this)
        return false;
      if (n.graph_ == graph_)
        return (index() < n.index());
      return (this < &n);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // pointer back to the graph container
    Graph* graph_;
    // This node's unique identifier
    size_type uid_;

    /** Constructor used by the graph class to construct valid nodes. */
    Node(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid) {};

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
   * @post result_node.value() == value
   * @post result_node.position() == position
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type new_uid = nodes_.size();
    size_type new_idx = size();
    internal_node n = {position, value, new_idx};
    nodes_.push_back(n);
    edges_.push_back(std::vector<internal_edge>{});
    i2u_.push_back(new_uid);
    return Node(this, new_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.uid_ >= nodes_.size())
      return false;
    if (n.index() >= size())
      return false;
    if (i2u_[n.index()] != n.uid_)
      return false;
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i >= size()){
      return Node();
    };
    return Node(this, i2u_[i]);
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Copy constructor */
    //Edge(const Edge& e):graph_(e.graph_), node1_uid_(e.node1_uid_), node2_uid_(e.node2_uid_) {};

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(graph_->nodes_[node1_uid_].idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(graph_->nodes_[node2_uid_].idx);
    }

    /** Return the value of this Edge
    *
    * Complexity: O(d)
    */
    edge_value_type& value() {
      for(size_type i = 0; i < node1().degree(); i ++){
        if (graph_->edges_[node1_uid_][i].node2_uid == node2_uid_){
          return graph_->edges_[node1_uid_][i].value;
        }
      }
      return graph_->edges_[0][0].value;
    }

    /** Return a constant reference to the value of this Edge */
    const edge_value_type& value() const {
      return const_cast<edge_value_type&>(value());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (node1() == e.node1() && node2() == e.node2())
        return true;
      if (node1() == e.node2() && node2() == e.node1())
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (e == *this)
        return false;
      if (e.graph_ == graph_ && e.node1_uid_ != node1_uid_)
        return (node1_uid_ < e.node1_uid_);
      if (e.graph_ == graph_ && e.node2_uid_ != node2_uid_)
        return (node2_uid_ < e.node2_uid_);
      return (this < &e);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_uid_;
    size_type node2_uid_;

    /** Constructor used by the graph class to construct valid edges. */
    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
    : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i >= num_edges()){
      return Edge();
    };
    auto it = edge_begin();
    for (size_type j = 0; j < i; ++ j){
      ++it;
    }
    return *(it);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return adj_index(a,b) != -1;
  }

  /** Find the index of node b in node a's adjacency list
   * @pre @a a and @a b are valid nodes of this graph
   * @return j if there exists an edge between a and b,
   *    and edges_[a.index()][j].node2_index = b.index(),
   *    otherwise return -1
   */
  int adj_index(const Node&a, const Node& b) const {
    if (!has_node(a) || !has_node(b)){
      return -1;
    }
    for (size_type i = 0; i < edges_[a.uid_].size(); ++ i){
      if (edges_[a.uid_][i].node2_uid == b.uid_)
        return i;
    }
    return -1;
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
    if (!has_node(a) || !has_node(b) || a == b){
      return Edge();
    }
    if (!has_edge(a,b)){
      edges_[a.uid_].push_back(internal_edge{b.uid_, value});
      edges_[b.uid_].push_back(internal_edge{a.uid_, value});
      num_edges_++;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_ = std::vector<internal_node>();
    edges_ = std::vector<std::vector<internal_edge>>();
    i2u_ = std::vector<size_type>();
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    /** @brief Return the Node at the iterator's current position
    *
    * @return a Node
    * @post this will not be modified
    */
    Node operator*() const {
      return graph_->node(index_);
    }

    /** @brief Return a reference to a NodeIterator at one past the iterator's current position
    *
    * @return a reference to a NodeIterator
    */
    NodeIterator& operator++() {
      index_++;
      return *this;
    }

    /** @brief Test whether this node iterator and @a n are equal
    *
    * @return a boolean value; true only if both iterators belong to the same graph and have the same node position
    * @post this will not be modified
    */
    bool operator==(const NodeIterator& n) const {
      return n.graph_ == graph_ && n.index_ == index_;
    }

   private:
    friend class Graph;

    NodeIterator(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {};

    Graph* graph_;
    size_type index_;

  };

  /** @brief Return an interator for the position of the node at index 0
  *
  * @return an interator of type node_iterator
  * @post this will not be modified
  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** @brief Return an interator for one past the final node in this graph (index num_nodes())
  *
  * @return an interator of type node_iterator
  * @post this will not be modified
  */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

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

    /** @brief Return the Edge at the iterator's current position
    *
    * @return an Edge
    * @post this will not be modified
    */
    Edge operator*() const {
      return Edge(graph_, node1_uid_, graph_->edges_[node1_uid_][node2_index_].node2_uid);
    }

    /** @brief Return a reference to an IncidentIterator at one past the iterator's current position
    *
    * @return a reference to an IncidentIterator
    */
    IncidentIterator& operator++(){
      node2_index_++;
      return *this;
    }

    /** @brief Test whether this incident iterator and @a i are equal
    *
    * @return a boolean value; true only if both iterators belong to the same graph and have the same position
    * @post this will not be modified
    */
    bool operator==(const IncidentIterator& i) const {
      return i.graph_ == graph_ && i.node2_index_ == node2_index_ && i.node1_uid_ == node1_uid_;
    }

   private:
    friend class Graph;

    IncidentIterator(const Graph* graph, size_type node1_uid, size_type node2_index)
    : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_index_(node2_index) {};

    Graph* graph_;
    size_type node1_uid_;
    size_type node2_index_;

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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

    /** @brief Return the Edge at the iterator's current position
    *
    * @return an Edge
    * @post this will not be modified
    */
    Edge operator*() const {
      return Edge(graph_, graph_->i2u_[node1_index_], graph_->edges_[graph_->i2u_[node1_index_]][node2_index_].node2_uid);
    }

    /** @brief Return an EdgeIterator for the position one past this iterator's position
    *
    * @return a reference to an EdgeIterator
    * @post if result != edge_end(), (*result).node1().index() < (*result).node2().index()
    */
    EdgeIterator& operator++() {
      ++node2_index_;
      findNextEdge();
      return *this;
    }

    /** @brief Test whether this EdgeIterator and @a ei are equal
    *
    * @return a boolean value; true only if both edge iterators belong to the same graph and have the same position
    * @post this will not be modified
    */
    bool operator==(const EdgeIterator& ei) const {
      return graph_ == ei.graph_ && node2_index_ == ei.node2_index_ && node1_index_ == ei.node1_index_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_index_;
    size_type node2_index_;

    EdgeIterator(const Graph* graph, size_type node1_index, size_type node2_index)
    : graph_(const_cast<Graph*>(graph)), node1_index_(node1_index), node2_index_(node2_index) {
      findNextEdge();
    };

    void findNextEdge(){
      while (node1_index_ < graph_->size()) {
        while (node2_index_ < graph_->edges_[graph_->i2u_[node1_index_]].size()){
          // only return the current edge if node1 < node2, otherwise, keep searching.
          // This is to avoid returning duplicate edges
          if (node1_index_ < graph_->nodes_[graph_->edges_[graph_->i2u_[node1_index_]][node2_index_].node2_uid].idx){
            return;
          }
          ++node2_index_;
        }
        // Increment the node iterator once we get to the end of the current node's incident edges
        // As long as the node iterator is still pointing to a valid node, reset the incident iterator and continue
        ++node1_index_;
        if (node1_index_ < graph_->size()){
          node2_index_ = 0;
        }
      }
      return;
    }

  };

  /** @brief Return an iterator for the edge at position 0 in the graph
  *
  * @return an edge_iterator
  * @post this will not be modified
  */
  edge_iterator edge_begin() const{
    if(num_edges() == 0){
      return edge_end();
    }
    return EdgeIterator(this, 0, 0);
  }

  /** @brief Return an iterator for the position one past the last edge in the graph
  *
  * @return an edge_iterator
  * @post this will not be modified
  */
  edge_iterator edge_end() const{
    if(num_edges() == 0){
      return EdgeIterator(this, num_nodes(), 0);
    }
    return EdgeIterator(this, num_nodes(), edges_[i2u_[num_nodes()-1]].size());
  }


  /** @brief Remove the edge connecting Node @node1 and Node @node2 from the graph
  *
  * Does not invalidate outstanding Edge objects, except for those equal to the Edge being deleted.
  * May invalidate outstanding Edge indices (new edge(i) != old edge(i)) and Edge/Incident Iterators
  *
  * @param[in] Node node1, one of the edge's endpoints
  * @param[in] Node node2, the other endpoint of the edge
  * @return 0 if the edge was not in the graph, 1 if the edge was successfully removed
  *
  * @pre has_node(node1) && has_node(node2)
  * @pre has_edge(node1, node2)
  *
  * @post !(has_edge(node1, node2))
  * @post new num_edges() = old num_edges() - 1
  * @post new node1.degree() = old node1.degree() - 1
  * @post new node2.degree() = old node2.degree() - 1
  *
  * Complexity: O(d)
  */
  size_type remove_edge(const Node& node1, const Node& node2){
    // find the index of node2 in node1's adjacency list
    int node1_index = adj_index(node1, node2);

    // if missing, return 0 (There was no edge between node1 and node2)
    if (node1_index == -1)
      return 0;

    // find the index of node1 in node2's adjacency list
    int node2_index = adj_index(node2, node1);

    //remove both adjacency list entries
    edges_[node1.uid_].erase(edges_[node1.uid_].begin() + node1_index);
    edges_[node2.uid_].erase(edges_[node2.uid_].begin() + node2_index);

    // decrease the number of edges
    num_edges_--;

    return 1;
  }


  /** @brief Remove Edge @e from the graph
  *
  * Does not invalidate outstanding Edge objects, except for those equal to the Edge being deleted.
  * May invalidate outstanding Edge indices (new edge(i) != old edge(i)) and Edge/Incident Iterators
  *
  * @param[in] Edge @e, the Edge to remove
  * @return 0 if @e was not in the graph, 1 if @e was successfully removed
  *
  * @pre has_node(e.node1()) && has_node(e.node2())
  * @pre has_edge(e.node1(), e.node2())
  *
  * @post !(has_edge(e.node1(), e.node2()))
  * @post new num_edges() = old num_edges() - 1
  * @post new e.node1().degree() = old e.node1().degree() - 1
  * @post new e.node2().degree() = old e.node2().degree() - 1
  *
  * Complexity: O(d)
  */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /** @brief Remove the Edge at position e_it from the graph
  *
  * Does not invalidate outstanding Edge objects, except for those equal to the Edge being deleted.
  * May invalidate outstanding Edge indices (new edge(i) != old edge(i)) and Edge/Incident Iterators
  *
  * @param[in] Edge @e, the Edge to remove
  * @return 0 if @e was not in the graph, 1 if @e was successfully removed
  *
  * @pre has_node(e.node1()) && has_node(e.node2())
  * @pre has_edge(e.node1(), e.node2())
  *
  * @post !(has_edge(e.node1(), e.node2()))
  * @post new num_edges() = old num_edges() - 1
  * @post new e.node1().degree() = old e.node1().degree() - 1
  * @post new e.node2().degree() = old e.node2().degree() - 1
  *
  * Complexity: O(d)
  */
  edge_iterator remove_edge(edge_iterator e_it){
    auto node1_index = e_it.node1_index_;
    auto node2_index = e_it.node2_index_;
    //remove the Edge
    remove_edge(*e_it);

    return EdgeIterator(this, node1_index, node2_index);
  }

  /** @brief Remove Node n (and all of its edges) from the graph
  *
  * Does not invalidate outstanding Node objects, except for those equal to the Node being deleted
  * May invalidate outstanding Node indices (new node(i) != old node(i)) and Node iterators
  *
  * @param[in] Node n, the node to remove
  * @returns 0 if the node was not in the graph, 1 if the node was successfully removed
  *
  * @pre has_node(n)
  *
  * @post num_nodes() = old num_nodes() - 1
  * @post num_edges() = old num_edges() - n.degree()
  * @post !(has_node(n))
  * @post !(has_edge(n,node(i))), 0 <= i < num_nodes()
  *
  * Complexity: O(N) + O(d^2), which is assumed to be O(N) since N >> d
  */
  size_type remove_node(const Node& n){

    // if n isn't a node of this graph, return 0
    if(!has_node(n))
      return (size_type)0;

    // make a list of Edges to remove (all Edges incident to n)
    std::vector<Edge> edgesToRemove {};
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      edgesToRemove.push_back(*it);
    }

    // Remove all Edges incident to the node
    for(auto it = edgesToRemove.begin(); it != edgesToRemove.end(); ++it){
      remove_edge(*it);
    }

    size_type idx = n.index();

    // Remove the Node's uid from the i2u mapping
    i2u_.erase(i2u_.begin() + idx);

    // Update all Nodes' indices that were shifted after the erase()
    for(; idx < i2u_.size(); idx++){
      nodes_[i2u_[idx]].idx = idx;
    }

    return (size_type)1;
  }

  /** @brief Remove the Node at the current position of n_it
  *
  * Does not invalidate outstanding Node objects, except for those equal to the Node being deleted
  * May invalidate outstanding Node indices (new node(i) != old node(i)) and Node iterators
  *
  * @param[in] node_iterator n_it, an iterator to the Node to be removed
  * @returns a node_iterator pointing to the next valid Node in the graph
  *
  * @pre has_node(n)
  *
  * @post num_nodes() = old num_nodes() - 1
  * @post num_edges() = old num_edges() - (*n_it).degree()
  * @post !(has_node(*n_it))
  * @post !(has_edge(*n_it,node(i))), 0 <= i < num_nodes()
  *
  * Complexity: O(N) + O(d^2), which is assumed to be O(N) since N >> d
  */
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
    return n_it;
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
    Point position;
    node_value_type value;
    size_type idx;
  };

  struct internal_edge{
    size_type node2_uid;
    edge_value_type value;
  };

  // Internal representation of nodes, indexed by uid
  std::vector<internal_node> nodes_;

  // Internal representation of edges.
  // This is an adjacency list representation
  // edges[i][j] = k means that there is an edge between nodes_[i] & nodes_[k]
  std::vector<std::vector<internal_edge>> edges_;

  // Indexed by i, look up to the uid of the node with index i
  std::vector <size_type> i2u_;

  // Quick reference to the number of edges in the graph
  size_type num_edges_ = 0;

};

#endif // CME212_GRAPH_HPP
