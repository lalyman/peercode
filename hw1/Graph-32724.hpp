#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <cassert>
#include <map>
#include <vector>

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
    const Point& position() const {
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
     * @pre !(graph_ == nullptr) Node mustb be valid
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
     * implemented by comparing ID numbers, which are unique to Node objects
     * in a graph.
     *
     * @pre graph_ == n.graph_ is required for comparison to be valid
     *          Otherwise, consider two nodes in different graphs but with the
     *          same ID number. This would not obey trichotomy (below)
     *
     * The node ordering relation must obey trichotomy: For any two nodes x and
     *  y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      assert(graph_ == n.graph_);
      return nodeID_ < n.index();
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
     * */
    Node(size_type ID, Graph* graph) : nodeID_(ID), graph_(graph) {}
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
    internalNode n {pos, val};
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
    return this == (const Graph*)n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert((i < size()) && i >= 0);
    Node returnNode(i, (Graph*) this);
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
      return graph_->node(graph_->edges_[edgeID_][0]);
    }

    /** Return the other node of this Edge
     * 
     * @pre Edge must be valid !(graph_ == nullptr)
     *          This is because data is stored in graph object, not in the Edge
     *          object itself.
     */
    Node node2() const {
      return graph_->node(graph_->edges_[edgeID_][1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * @pre !(graph_ == nullptr) Edge must be valid
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool sameGraph = graph_ == e.graph_;
      auto nodePair = graph_->edges_[edgeID_];
      auto eNodePair = e.graph_->edges_[e.edgeID_];
      bool match = nodePair[0] == eNodePair[0] && nodePair[1] == eNodePair[1];
      bool cMatch = nodePair[1] == eNodePair[0] && nodePair[0] == eNodePair[1];
      return (match || cMatch) && (sameGraph);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning. Here, we compare
     * using edgeID numbers, which are unique within a single graph.
     *
     * @pre graph_ == e.graph_ Edges must belong to same graph
     *          Otherwise, it is possible to violate trichotomy. Consider two
     *          edges in two different graphs, but with the same ID number.
     * @pre !(graph_ == nullptr) Edge must be valid.
     */
    bool operator<(const Edge& e) const {
      assert(graph_ == e.graph_);
      return edgeID_ < e.edgeID_;
    }

    /** Reorder the edge nodes.
     *
     * @pre !(graph_ == nullptr) Edge must be valid
     * @post edge.reorder().node1() == edge.node2()
     * @post edge.reorder().node2() == edge.node1()
     *
     * rvalue reference is not used to make the swap because size_type is very
     * lightweight, so creating the extra copy is not a hude deal.
     */
    void reorder() const {
      size_type temp = graph_->edges_[edgeID_][0];
      graph_->edges_[edgeID_][0] = graph_->edges_[edgeID_][1];
      graph_->edges_[edgeID_][1] = temp;
      return;
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
   * this case, I achieve O(1).
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    auto it = nodes_[a.index()].incident.find(b.index());
    return it != nodes_[a.index()].incident.end();
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
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
  Edge add_edge(const Node& a, const Node& b) {
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
      std::vector<size_type> entry {a.index(), b.index()};
      edges_.push_back(entry);
      return edge(num_edges()-1);
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
   * Implemented as a wrapper of std::map<size_type, size_type>::iterator*/
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
    std::map<size_type, size_type>::iterator it_;

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
    IncidentIterator(std::map<size_type, size_type>::iterator it,
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
    internalNode(const Point& pIn, const V& valIn) :
      p(pIn),
      val(valIn),
      incident() {}

    /** Position. */
    Point p;

    /** Value. */
    node_value_type val;

    /**The key represents the node ID of an adjacent node. The value represents
     * the edge ID of the edge that connects the two nodes. Say edge with ID
     * number x connects nodes with ID numbers y and z. The internalNode for
     * node y has a map with one entry: KEY is z, VALUE is x.
     */
    std::map<size_type, size_type> incident;
  };

  /** Data container for nodes. Each node is represented as an internalNode
   * struct. Vector is indexed by Node ID numbers (size_type). This allows for
   * O(1) return of Point based on NodeID. Vector allows us to add Nodes in
   * O(1) amortized.
   */
  std::vector<internalNode> nodes_;

  /** Data container for edges. Each edge is represented as a vector with two
   * size_type elements. Each size_type corresponds to the ID number of a node.
   * The ID number of the vector is used to index into the (outer) vector. This
   * ID number can be found O(1), since a map is used in the internalNode
   * struct. See above.
   */
  std::vector<std::vector<size_type>> edges_;
};

#endif // CME212_GRAPH_HPP
