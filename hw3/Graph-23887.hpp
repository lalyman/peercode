#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V,  typename E>
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
  /** Specifying node value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Specifying edge value type. */
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

  /** Construct an empty graph.
   * Uses default constructors for each of the member variables. List
   * initialization is used.
   */
  Graph() : nodes_(), edges_() {}

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
     *
     * Here, pointer graph_ is used as a test of validity. Value will be
     * nullptr, unless the Node is initialized by class Graph, and graph_
     * pointer points toward the comprising Graph object.
     */
    Node() : nodeID_(0), graph_(nullptr) {}

    /** Return this node's position.
     *
     * @pre Node must be valid !(graph_ == nullptr)
     *          This is because data is stored in underlying Graph class, not
     *          in the node itself.
     */
    Point& position() const {
      assert(graph_);   // Check for validity
      return graph_->nodes_[nodeID_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     *
     * Will return default value of 0 if called on invalid node.
     */
    size_type index() const {
      return nodeID_;
    }

    /** Return Lvalue reference to node valule.
     *
     * @pre !(graph_ ==  nullptr) Node must be valid
     *
     * Will return default value of V if called on invalid node.
     */
    node_value_type& value() {
      return graph_->nodes_[nodeID_].val;
    }

    /** Return Lvalue reference to const node value.
     *
     * @pre !(graph_ == nullptr) Node must be valid
     *
     * Will return default value of V if called on invalid node.
     */
    const node_value_type& value() const {
      return graph_->nodes_[nodeID_].val;
    }

    /** Return number of nodes connected to a given node.
     *
     * @pre !(graph_ == nullptr) Node must be valid
     *
     */
    size_type degree() const {
      return graph_->nodes_[nodeID_].incident.size();
    }

    /** Return iterator to beginning of incident nodes.
     *
     * @pre !(graph_ == nullptr) Node must be valid
     *
     * See implementation details below.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_->nodes_[nodeID_].incident.begin(),
                              nodeID_,
                              graph_);
    }

    /** Return iterator to end of incident nodes.
     *
     * @pre !(graph_ == nullptr) Node must be valid
     *
     * See implementation details below.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_->nodes_[nodeID_].incident.end(),
                              nodeID_,
                              graph_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (nodeID_ == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning. Here, it is
     * implemented by comparing addresses.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x and
     *  y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return &(graph_->nodes_[nodeID_]) < &(n.graph_->nodes_[n.nodeID_]);
    }

   private:
    /** Allow Graph to access Node's private member data and functions. */
    friend class Graph;

    /** Node ID Number */
    size_type nodeID_;

    /** Pointer to graph */
    Graph* graph_;

    /** Constructor for use by Graph object.
     * @param[in] ID Node ID number
     * @param[in[ graph pointer to graph that node will be part of
     *
     * The cast here is required becase we do not want graph_ to be const, but
     * this ctor is used in certain instances where the second argument is
     * const. I am not actually changing graph_ at all here. I am simply
     * storing it as a non-const member variable.
     */
    Node(size_type ID, const Graph* graph): nodeID_(ID), graph_((Graph*) graph)
    {}
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
   * @param[in] pos The new node's position
   * @param[in] val Value with which the node is initialized
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& pos,
                const node_value_type& val = node_value_type()) {
    internalNode n;
    n.p = pos;
    n.val = val;
    n.incident = std::unordered_map<size_type, size_type> {};
    nodes_.push_back(n);
    Node returnNode(nodes_.size() - 1, this);
    return returnNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * No need to check nodeID < size(), since nodeID is set appropriately in
   * Node constructor at the same time graph_ is set. Thus, the test used here
   * is both necessary and sufficient. All nodes in graph have valid ID number.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert((i < size()) && i >= 0);
    Node returnNode(i, this);
    return returnNode;
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
    /** Construct an invalid Edge.
     *
     * graph_ pointer is used as a test of validity. Value will be nullptr
     * unless Edge is constructed by class Graph, and graph_ points to the
     * comprising Graph object.
     */
    Edge() : edgeID_(0), graph_(nullptr) {}

    /** Return a node of this Edge
     * 
     * @pre Edge must be valid !(graph_ == nullptr)
     *          This is because data is stored in graph object, not in the Edge
     *          object itself. 
     */
    Node node1() const {
      return graph_->node(graph_->edges_[edgeID_].nodes[0]);
    }

    /** Return the other node of this Edge
     *
     * @pre Edge must be valid !(graph_ == nullptr)
     *          This is because data is stored in graph object, not in the Edge
     *          object itself.
     */
    Node node2() const {
      return graph_->node(graph_->edges_[edgeID_].nodes[1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * @pre !(graph_ == nullptr) Edge must be valid
     *
     * Equal edges represent the same undirected edge between two nodes. Since
     * already checking if in same graph, it is sufficient to check that
     * indices are the same.
     */
    bool operator==(const Edge& e) const {
      bool sameGraph = graph_ == e.graph_;
      return (edgeID_ == e.edgeID_) && (sameGraph);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning. Here, we compare
     * using addresses.
     */
    bool operator<(const Edge& e) const {
      return &(graph_->edges_[edgeID_]) < &(e.graph_->edges_[e.edgeID_]);
    }

    /** Reorder the edge nodes.
     *
     * @pre !(graph_ == nullptr) Edge must be valid
     * @post edge.reorder().node1() == edge.node2()
     * @post edge.reorder().node2() == edge.node1()
     *
     * rvalue reference is not used to make the swap because size_type is very
     * lightweight, so creating the extra copy is not a huge deal. This is
     * somewhat ugly, but I unfortunately did not have enough time to change
     * this part of the implementation yet. Potential workaround: store
     * variable representing node1() in each Edge. This allows for lookup of
     * node2() without having to store it (eliminates half the data in egdes_.
     * Still unsure about how to get rid of that last half though... without
     * losing a lot of speed
     */
    void reorder() const {
      size_type temp = graph_->edges_[edgeID_].nodes[0];
      graph_->edges_[edgeID_].nodes[0] = graph_->edges_[edgeID_].nodes[1];
      graph_->edges_[edgeID_].nodes[1] = temp;
      return;
    }

    /** Return Lvalue reference to edge value.
     *
     * @pre !(graph_ ==  nullptr) Edge must be valid
     *
     * Will return default value of E if called on invalid node.
     */
    edge_value_type& value() {
      return graph_->edges_[edgeID_].val;
    }

    /** Return Lvalue reference to edge value.
     *
     * @pre !(graph_ ==  nullptr) Edge must be valid
     *
     * Will return default value of E if called on invalid node.
     */
    const edge_value_type& value() const {
      return graph_->edges_[edgeID_].val;
    }

   private:
    /** Allow Graph to access Edge's private member data and functions. */
    friend class Graph;

    /** Variable for edge ID. */
    size_type edgeID_;

    /** Pointer to graph. */
    Graph* graph_;

    /** Constructor for use by Graph object
     * @param[in] ID ID number of edge to be created
     * @param[in] graph Graph to which edge belongs
     *
     */
    Edge(size_type ID, Graph* graph) : edgeID_(ID), graph_(graph) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1).
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1).
   */
  Edge edge(size_type i) const {
    assert((0 <= i) && (i < num_edges()));
    return Edge(i, (Graph*) this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1) since the lookup in an unordered map is that
   * fast.
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (has_node(a) && has_node(b)) {
      auto it = nodes_[a.index()].incident.find(b.index());
      return it != nodes_[a.index()].incident.end();
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @param[in] a       First Node at end of edge
   * @param[in] n       Second Node at other end of edge
   * @param[in] val     Value associated with the edge, default if unspec
   * @param[out]        Edge    Edge object just created, or already present
   *
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects. Note:
   * The order of node1() and node2() may be swapped in existing Edge objects.
   * This is assumed not to be a problem since edges are symmetric.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. In
   * this case, I achieve O(1) using the has_edge function written above. The
   * details of this scheme are explained further below.
   */
  Edge add_edge(const Node& a,
                const Node& b,
                const edge_value_type& val = edge_value_type()) {
    assert(!(a == b));                  // Check distinctness
    assert(has_node(a) && has_node(b)); // Check validity and belonging
    if (has_edge(a, b)) {
      Edge e = edge(nodes_[a.index()].incident[b.index()]);
      if (e.node1() != a) {e.reorder();}
      return e;
    } else {

      // Add to each node association list
      nodes_[a.index()].incident[b.index()] = num_edges();
      nodes_[b.index()].incident[a.index()] = num_edges();

      // Add to edges_register
      internalEdge entry {};
      entry.val = val;
      entry.nodes = std::vector<size_type> {a.index(), b.index()};
      edges_.push_back(entry);
      return edge(num_edges()-1);
    }
  }

  //
  // REMOVAL METHODS
  //

  /** Remove specified node from Graph.
   *
   * @param[in]         n       Node to be removed from graph
   * @param[out]        size_type       1 if Node successfully removed
   *
   * @pre this->has_node(@a n)  Node is valid and contained in graph.
   * @post new *it == old *it for all elements but n and the last element
   * @post new node(i) == old node(i) for i != num_nodes-1, n.index()
   * @post node(i).index() == i for 0 <= i < num_nodes()
   * @post node(n.index()) == n for 0 <= i < num_nodes()
   * @post new this->size()  = old this->size() - 1
   * @post new this->num_edges()  = old this->num_edges() - @a n.degree()
   * @post new degree() = old degree() - 1 for each incident node
   *
   * The last Node is invalidated, and the removed node points to old last node
   * Iterators are preserved for all elements but the removed node (new node in
   * this position is previous end node) and the end node (effectively
   * removed.) node_end() is invalidated.
   *
   * All edges connected to n are invalidated, meaning all edges with ID
   * given by values of nodes_[n.index()].incident are removed. Edge iterators
   * after the earliest-removed edge are invalidated. edge_end() is invalidated
   * as well.
   *
   * Incident Edge iterators for unaffected nodes are preserved. Incident edge
   * iterators for nodes that were once connected to n are only preserved for
   * iterators up to the edge that was removed. Incident iterators after this
   * point are invalidated, including edge_end().
   *
   * Fuction is O(degree), which we assume to be small.
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      // Remove incident edges first
      auto inc = nodes_[n.index()].incident;
      for (auto it = inc.begin(); it != inc.end(); ++it) {
        remove_edge(n, node(it->first));
      }
      // Swap n with last node and remove last entry, in O(degree)
      // See documentation for this below, near definition
      swap_nodes(n, node(num_nodes()-1));
      nodes_.pop_back();
      return 1;
    }
    return 0;
  }

  /** Remove node specified by iterator from Graph.
   *
   *
   * @param[in]         n_it    Iterator to node thats should be deleted
   * @param[in/out/     n_it    Iterator to new element in position
   *
   * Simply calls the appropriate method for the dereferenced iterator, see
   * documentation above.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;                // Still points to the same position in vector
  }

  /** Swap nodes n1 and n2 in internal memory. 
   *
   * This method exchanges ID numbers of nodes n1 and n2.
   *
   * @param[in] n1      First node to be swapped
   * @param[in] n2      Second node to be swapped
   *
   * @pre !(n1.graph_ == nullptr) n1 is a valid node of this graph
   * @pre !(n2.graph_ == nullptr) n2 is a valid node of this graph
   * @post new n1 == old n2
   * @post new n2 == old n1
   * @post new node(i) == old node(i) for i != n1.index() and i != n2.index()
   * @post new *node_it[i] == old *node_it[i] for i != n1.index(), n2.index()
   * @post new node(j).incident[i] = old node(j).incident[i] for j != n1.index
   *            and j != n2.index(), 0 <= i < node.degree()
   *
   * Users should note that existing iterators might skip/double the swapped
   * entries even though the values of iterators are preserved for nodes that
   * were not swapped.
   * 
   * Edges are preserved for all unaffected edge ID's (edges that do not
   * contain n1 or n2. Edge iterators are preserved for edges that are
   * unaffected, but underlying data in sequence is changed by method.
   * 
   * Incident iterators are invalidated for affected edges, after the
   * lowest-indexed change in the iterator list.
   *
   * Runs in O(degree), which we assume to be small.
   */
  void swap_nodes(const Node& n1, const Node& n2) {
    // Account for edges incident to n1
    auto inc = nodes_[n1.index()].incident;
    for (auto it = inc.begin(); it != inc.end(); ++it) {
      auto& nVect = edges_[it->second].nodes;
      nVect[0] = nVect[0] == (n1.index()) ? n2.index() : nVect[0];
      nVect[1] = nVect[1] == (n1.index()) ? n2.index() : nVect[1];
      nodes_[it->first].incident[n2.index()] = it->second;
      if (n1 != n2) {
        nodes_[it->first].incident.erase(n1.index());
      }
    }

    // Account for edges incident to n2
    inc = nodes_[n2.index()].incident;
    for (auto it = inc.begin(); it != inc.end(); ++it) {
      auto& nVect = edges_[it->second].nodes;
      nVect[0] = nVect[0] == (n2.index()) ? n1.index() : nVect[0];
      nVect[1] = nVect[1] == (n2.index()) ? n1.index() : nVect[1];
      nodes_[it->first].incident[n1.index()] = it->second;
      if (n1 != n2) {
        nodes_[it->first].incident.erase(n2.index());
      }
    }

    // Actually swap the internalNodes now
    std::swap(nodes_[n1.index()], nodes_[n2.index()]);
  }

  /** Remove edge between two specified nodes
   *
   * @param[in]         n1      First node specifying edge to be deleted.
   * @param[in]         n2      Second node specifying edge to be deleted.
   * @param[out]        size_type    1 if Edge successfully removed
   *
   * @post if !has_edge(n1,  n2), then nothing happens and returns 0.
   * Say that edge e does exist between nodes n1 and n2
   * @post new node(i) = old node(i) for 0 <= i < num_nodes()
   * @post new NodeIterator = old NodeIterator, still valid
   * @post new edge(i) = old edge(i) for i != e.edgeID_, i != num_edges()-1
   * @post new *it = old *it for edges that are not e, edge(num_edges(0-1)
   * @post edge(i).edgeID_ = i for 0 <= i < num_edges()
   * @post edge(e.edgeID_) = e for 0 <= i < num_edges()
   * @post new n1.degree() = old n1.degree() - 1;
   * @post new n2.degree() = old n2.degree() - 1;
   * @post new size() = old size()
   * @post new num_edges = old num_edges - 1;
   *
   * Node iterators and nodes are all still valid
   *
   * Old edge(num_edgdes()-1) is now in e.edheID_ position. Other edges are
   * preserved. Similarly, edge itertors are preserved for edges apart from the
   * last edge and the removed edge.
   *
   * Incident edge iterators for unaffected nodes are preserved. Incident edge
   * iterators for nodes that were once a part of the removed edge are only
   * preserved up to the edge that was removed. After this point, EdgeIterators
   * including edge_end() are invalidated.
   *
   * Runs in O(degree), which we assume to be small.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      // Fix the nodes for old back() before swapping with e
      auto& backNodes = edges_.back().nodes;
      size_type edgeNum = nodes_[n1.index()].incident[n2.index()];
      nodes_[backNodes[0]].incident[backNodes[1]] = edgeNum;
      nodes_[backNodes[1]].incident[backNodes[0]] = edgeNum;

      // Swap edge with back() in O(d) and erase last element
      std::swap(edges_[edgeNum], edges_.back());
      edges_.pop_back();

      // Erase edge from association registers
      nodes_[n1.index()].incident.erase(n2.index());
      nodes_[n2.index()].incident.erase(n1.index());
      return 1;
    }
    return 0;
  }

  /** Remove specified edge from graph
   *
   * @param[in]         e       Edge to be removed.
   * @param[out]        size    New num_edges(), after edge is removed.
   *
   * Simlply calls appropriate method for edge between nodes e.node1() and
   * e.node2(). See documentation above.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove specified edge from graph
   *
   * @param[in]         e_it    Iterator to edge to be removed.
   * @param[out]        size    New num_edges(), after edge is removed.
   *
   * Simply calls appropriate method for dereferenced iterator. See
   * documentation above.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;                // Still points to the same positio in vector
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
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
    NodeIterator() : itNum_(0), graph_(nullptr) {}

    /** Return Node object that is being pointed to by the iterator. */
    Node operator*() const {
      return graph_->node(itNum_);
    }

    /** Increment iterator to point to next Node object in graph. */
    NodeIterator& operator++() {
      itNum_++;
      return *this;
    }

    /** Test for equality, iterators must be in the same graph. */
    bool operator==(const NodeIterator& NI) const {
      return (itNum_ == NI.itNum_) && (graph_ == NI.graph_);
    }

    /** Constructor for use by Graph class, to return begin and end.
     * @param[in] itNum The node number to which the iterator should point
     * @param[in] graph Pointer to the graph in which the node is stored
     */
    NodeIterator(size_type itNum, Graph* graph) : itNum_(itNum), graph_(graph)
      {}

   private:
    friend class Graph;

    /** Iterator is represented by index number of node (nodeID). */
    size_type itNum_;

    /* Iterator needs to know which graph it is a part of. */
    Graph* graph_;

  };

  /** Return beginning iterator of Nodes in graph. */
  node_iterator node_begin() const {
    node_iterator returnNI {0, (Graph*) this};
    return returnNI;
  }

  /** Return end iterator of Nodes in graph. */
  node_iterator node_end() const {
    node_iterator returnNI {num_nodes(), (Graph*) this};
    return returnNI;
  }

  /** Return iterator to specified node in graph.
   *
   * @pre n.graph_ == this, Node must belong to this graph.
   *
   * This functionality is useful in applications, shortest_path in particular.
   */
  node_iterator node_it(const Node& n) const {
    node_iterator returnNI {n.index(), (Graph*) this};
    return returnNI;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   *
   * Wrapper of std::unordered_map<size_type, size_type>::iterator*/
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : it_(), nodeNum_(0), graph_(nullptr) {}

    /** Overload dereferencing operator. */
    Edge operator*() const {
      Edge e = graph_->edge(it_->second);
      if (e.node1() !=  graph_->node(nodeNum_)) {e.reorder();}
      return e;
    }

    /** Increment iterator by one. */
    IncidentIterator& operator++() {
      ++it_;
      return *this;
    }

    /** Check equality of iterators. */
    bool operator==(const IncidentIterator& II) const {
      bool sameGraph = graph_ == II.graph_;
      bool sameIT = it_ == II.it_;
      bool sameNode = nodeNum_ == II.nodeNum_;
      return (sameGraph) && (sameIT) && (sameNode);
    }

   private:
    /** Allow for access by class Graph. */
    friend class Graph;

    /** Iterator is actually just a wrapper around this container iterator. */
    std::unordered_map<size_type, size_type>::iterator it_;

    /** Keep track of which node the iterator is looking around. */
    size_type nodeNum_;

    /** Keep track of which graph is being explored. */
    Graph* graph_;

    /** Constructor for use, by Node class.
     * @param[in] it Iterator to map that stores adjacency information
     * @param[in] nodeNum The index number of the node for which iterator is
     *            requested.
     * @param[in] graph Pointer to the graph in which nodes are contained
     *
     * @pre node(it->first).graph_ == graph, Node must be in graph
     *
     */
    IncidentIterator(std::unordered_map<size_type, size_type>::iterator it,
                     size_type nodeNum,
                     Graph* graph) : it_(it), nodeNum_(nodeNum), graph_(graph)
    {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : itNum_(0), graph_(nullptr) {}

    /** Overload dereferencing operator. */
    Edge operator*() const {
      return graph_->edge(itNum_);
    }

    /** Overload increment operator. */
    EdgeIterator& operator++() {
      itNum_++;
      return *this;
    }

    /** Test for equality of iterators. */
    bool operator==(const EdgeIterator& EI) const {
      bool sameGraph = graph_ == EI.graph_;
      bool sameNum = itNum_ == EI.itNum_;
      return sameGraph && sameNum;
    }

   private:
    friend class Graph;
    /** Keep track of position along vector of edges. */
    size_type itNum_;

    /** Pointer to graph that edges are contained in. */
    Graph* graph_;

    /** Constructor for use by Graph.
     * @param[in] itNum Position in edge vector that iterator will point to
     * @param[in] graph Pointer to graph that edges belong to.
     */
    EdgeIterator(size_type itNum, Graph* graph) : itNum_(itNum), graph_(graph)
    {}
  };

  /** Specify beginning edge iterator. */
  edge_iterator edge_begin() const {
    edge_iterator returnEI {0, (Graph*) this};
    return returnEI;
  }

  /** Specify end edge iterator. */
  edge_iterator edge_end() const {
    edge_iterator returnEI {num_edges(), (Graph*) this};
    return returnEI;
  }

 private:
  /** structure used to contain information about nodes. Each node is given a
   * Point to represent position, and a value_type to represent that node's
   * value. The map (incident) is used to keep track of nodes adjacent to each
   * node in the graph. This allows for O(1) implementation of has_edge and
   * add_edge methods.
   */
  struct internalNode {
    /** Position. */
    Point p;

    /** Value. */
    node_value_type val;

    /**The key represents the node ID of an adjacent node. The value represents
     * the edge ID of the edge that connects the two nodes. Say edge with ID
     * number x connects nodes with ID numbers y and z. The internalNode for
     * node y has a map with one entry: KEY is z, VALUE is x. This allows for
     * O(1) lookup of edge number given two nodes.
     */
    std::unordered_map<size_type, size_type> incident;
  };

  /** Data container for nodes. Each node is represented as an internalNode
   * struct. Vector is indexed by Node ID numbers (size_type). This allows for
   * O(1) return of Point based on NodeID. Vector allows us to add Nodes in
   * O(1) amortized.
   */
  std::vector<internalNode> nodes_;

  /** Data structure for internal storage of edges, with values. Each edge
   * contatins a vector that specifies the nodes at either edge. Additionally,
   * each edge contains a value associated with it.
   */
  struct internalEdge {
    /** Value. */
    edge_value_type val;

    /** Nodes connected to edge. Represented using size_type node ID numbers.*/
    std::vector<size_type> nodes;
  };

  /** Data container for edges. Each edge is represented as an internalEdge
   * object, which is defined above. The edge ID is used to index into this
   * vector. The edge ID can be found O(1) given two nodes, using the map
   * inside each internalNode object. Node data is stored twice, but this is
   * necessary for fastest lookup. (Basically, I didn't have enough time to
   * address this and I'll try to take care of it for the last assignment.)
   */
  std::vector<internalEdge> edges_;
};

#endif // CME212_GRAPH_HPP
