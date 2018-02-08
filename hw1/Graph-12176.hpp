#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
 
/** @file Graph.hpp
 * @brief An undirected graph type
 * 
 * These methods were implemented with helpful collaboration from:
 * Joe Higgins
 * Amel Awadelkarim
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 * 
 * Template allows nodes to support user-specified value of type node_value_type
 */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct My_node;
  struct My_edge;

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
  /** Synonym for V, the value of a node (following STL conventions). */
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
  Graph() {
    // HW0: YOUR CODE HERE
    // blank because node_list, edge_list, and adj_list are empty by default
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
    Node() : g_(nullptr), ind_(-1) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->node_list[ind_].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return g_->node_list[ind_].ind;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Getter for the value stored in node_value_type object */
    const node_value_type& value() const {
      return g_->node_list[ind_].val;
    }

    /** Setter for the value stored in node_value_type object */
    node_value_type& value() {
      return const_cast<node_value_type&>(static_cast <Node const &>(*this).value()); 
    }

    /** Return the degree of the node (number of incident edges */
    size_type degree() const {
      return g_->adj_list[ind_].size();
    }

    /** Return incident_iterator to the first incident edge in adj_list[ind_] */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, ind_, 0);
    }

    
    /** Return incident_iterator to the last incident edge in adj_list[ind_] */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, ind_, (*this).degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return this->g_ == n.g_ && this->ind_ == n.ind_;
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
      // HW0: YOUR CODE HERE
      if (this->ind_ == n.ind_) {
        return this->g_ < n.g_;
      } else {
        return this->ind_ < n.ind_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* g_;
    size_type ind_;
    Node(const Graph* g, size_type ind) 
      : g_(const_cast<Graph*>(g)), ind_(ind) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_list.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return (*this).size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    node_list.push_back(My_node(position, node_list.size(), value));
    std::vector<std::tuple<size_type, size_type>> v {};
    adj_list.push_back(v);
    return Node(this, node_list.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1). 
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return node_list[n.index()].pos == n.position();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
    Edge() : g_(nullptr), ind_(-1) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(g_, g_->edge_list[ind_].a); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g_, g_->edge_list[ind_].b);
    }
    
    /** Return the index of this Edge */
    size_type index() const {
      return this->ind_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->node1() == e.node1() && this->node2() == e.node2()) ||
             (this->node1() == e.node2() && this->node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return this->ind_ < e.ind_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* g_;
    size_type ind_;
    Edge(const Graph* g, size_type ind) 
      : g_(const_cast<Graph*>(g)), ind_(ind) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i);
  }
 
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    const std::vector<std::tuple<size_type, size_type>> &temp {adj_list[a.index()]};

    for (size_type i = 0; i < temp.size(); i++) {
      if (std::get<0>(temp[i]) == b.index()) {
        return true;
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
   * Complexity: O(num_nodes() + num_edges())
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    const std::vector<std::tuple<size_type, size_type>> &temp {adj_list[a.index()]};
    
    /* First look for existing edge */
    for (size_type i = 0; i < temp.size(); i++) {
      if (std::get<0>(temp[i]) == b.index()) {
        return edge(std::get<1>(temp[i]));
      }
    }

    /* If edge does not already exist, make a new edge */
    My_edge new_edge = My_edge(num_edges(), a.index(), b.index());
    edge_list.push_back(new_edge);
    adj_list[a.index()].push_back(std::make_tuple(b.index(), num_edges()-1));
    adj_list[b.index()].push_back(std::make_tuple(a.index(), num_edges()-1));
    return edge(num_edges()-1);
  }

  /** Returns the adjacency list */ 
  std::vector<std::vector<size_type>> get_edges() {
    return adj_list;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    edge_list.clear();
    adj_list.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereference operator
     *  @return The current Node pointed to
     */ 
    Node operator*() const {      
      return g_->node(curr);
    }
 
    /** Increment operator
     *  @return NodeIterator to sequential node (by index)
     */
    NodeIterator& operator++() {
      this->curr++;
      return *this;
    }

    /** Two iterators are equal if their graph and indices are equal */
    bool operator==(const NodeIterator& rhs) const {
      return this->g_ == rhs.g_ && this->curr == rhs.curr;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* g_;
    size_type curr;
    NodeIterator(const Graph* g, size_type n) : 
      g_(const_cast<Graph*>(g)), curr(n) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /** Construct node_iterator to point to the first node (by index)
   *  Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Construct node_iterator to point to the last node (by index)
   *  Complexity: O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(this, (*this).num_nodes());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereference operator
     *  @return Edge object being pointed to
     */ 
    Edge operator*() const {      
      return g_->edge(std::get<1>(g_->adj_list[n_ind][curr]));
    }
 
    /** Increment operator
     *  @return IncidentIterator to sequential incident edge
     *          by the order stored within adj_list[n_ind]
     */
    IncidentIterator& operator++() {
      curr++;
      return *this;
    }

    /** Two IncidentIterators are equal if their graphs, node indices,
     *  and incident edge index are equal.
     */
    bool operator==(const IncidentIterator& rhs) const {
      return (this->g_ == rhs.g_ && this->curr == rhs.curr)
       && (this->n_ind == rhs.n_ind);
    }
    
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* g_;
    size_type n_ind; // Index of current node
    size_type curr;  // Index within list of incident edges
    IncidentIterator(const Graph* g, size_type i, size_type n) : 
      g_(const_cast<Graph*>(g)), n_ind(i), curr(n) {
    }
    
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
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    /** Dereference operator
     *  @return Edge being pointed to
     */ 
    Edge operator*() const {      
      return g_->edge(curr);
    }
 
    /** Increment operator
     *  @return EdgeIterator pointing to the next edge (by index)
     */
    EdgeIterator& operator++() {
      this->curr++;
      return *this;
    }

    /** Two EdgeIterators are equal if their graphs and indices are equal */
    bool operator==(const EdgeIterator& rhs) const {
      return this->g_ == rhs.g_ && this->curr == rhs.curr;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* g_;
    size_type curr;
    EdgeIterator(const Graph* g, size_type n) : 
      g_(const_cast<Graph*>(g)), curr(n) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns edge_iterator to the first edge with index 0
   *  Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Returns edge_iterator to one past the last edge
   *  Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, (*this).num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Struct to store node information that proxy will refer to
   *  @param pos   Point object storing the node's (x,y,z) coordinates
   *  @param ind   Index of the node
   *  @param val   Value stored by the node (template type)
   */
  struct My_node {
    Point pos;
    size_type ind;
    node_value_type val;
    // Constructor
    My_node(Point p, size_type i, node_value_type v) : pos(p), ind(i), val(v) {
    }
  };

  /** Struct to store edge information that proxy will refer to
   *  @param ind   Index of the edge
   *  @param a     Index of node1() of the edge
   *  @param b     Index of node2() of the edge
   */
  struct My_edge {
    size_type ind;
    size_type a;
    size_type b;
    // Constructor
    My_edge(size_type i, size_type first, size_type second) :
      ind(i), a(first), b(second) {
    }
  };

  std::vector<My_node> node_list; // Store all My_node objects of the graph
  std::vector<My_edge> edge_list; // Store all My_edge objects of the graph
  /** Adjacency list representing the connections in the graph
   *  adj_list[i] contains a vector of tuples (j, k) such that edge(k)
   *  connects node(i) to node(j).
   */
  std::vector<std::vector<std::tuple<size_type, size_type>>> adj_list;
};

#endif  // CME212_GRAPH_HPP
