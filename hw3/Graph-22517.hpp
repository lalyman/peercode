#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <functional>
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
 private:
  class SecretNode;
  class AdjEdge;
  std::vector<SecretNode> nodes_;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Synonym for V (following STL conventions). */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for E (following STL conventions). */
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
  Graph() {
  	//vectors already initialized to empty vectors
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
      g_ = nullptr;
      ind_ = size_type(-1);
    }

    /** Return this node's position. */
    Point& position() {
      return g_->nodes_[g_->i2u_[ind_]].p_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return g_->nodes_[g_->i2u_[ind_]].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return ind_;
    }

    /** Return this node's value. */
    node_value_type& value() {
      return g_->nodes_[g_->i2u_[ind_]].v_;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return g_->nodes_[g_->i2u_[ind_]].v_;
    }

    /** Return this node's degree, the number of nodes it is adjacent to. */
    size_type degree() const {
      return g_->nodes_[g_->i2u_[ind_]].adj_.size();
    }

    /** Return an iterator starting over edges incident to this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, ind_, size_type(0));
    }

    /** Return an iterator ending over edges incident to this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, ind_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->ind_ == n.ind_) && (this->g_ == n.g_));
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
      if (*this == n) {
      	//same node so not less than
      	return false;
      }
      else if (this->ind_ == n.ind_) {
      	//same index, different graph
      	std::less<Graph*> graph_ptr_less;
      	return graph_ptr_less(this->g_, n.g_);
      }
      else{
      	return (this->ind_ < n.ind_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* g_;
    size_type ind_;

    Node(const Graph* g, size_type ind) {
      g_ = const_cast<Graph*>(g);
      ind_ = ind;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& val = node_value_type()) {
    size_type u = nodes_.size();
    i2u_.push_back(u);
    size_type i = i2u_.size() -1;
    nodes_.push_back(SecretNode(position, val, i));
    return Node(this, i);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= i2u_.size())
      return false;
    return (this->node(n.ind_) == n);
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

  /** Remove a node from the graph, returning the index of the removed node.
   * @param[in] n The node to remove
   * @post new num_nodes() == old num_nodes() - 1
   * @post Removes all edges incident to n
   * @post Nodes node(i) are preserved for n.index() < i < num_nodes()
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      for (auto efirst = n.edge_begin(); efirst != n.edge_end(); ++efirst) {
          auto adjacent_edge = *efirst;
          one_way_remove_edge(adjacent_edge.node2(), n);
      }
      nodes_[i2u_[n.index()]].adj_.clear();
      i2u_[n.index()] = i2u_[num_nodes() -1];
      i2u_.pop_back();
      nodes_[i2u_[n.index()]].i_ = n.index(); 
    }
    return n.index();
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
      g_ = nullptr;
      n1id_ = size_type(-1);
      n2id_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, n1id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, n2id_);
    }

    /** Return length of Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return this edge's value. */
    edge_value_type& value() {
      size_type i = 0;
      for( ; i < node1().degree(); ++i) {
        if (g_->nodes_[g_->i2u_[n1id_]].adj_[i].nid_ == g_->i2u_[n2id_])
          break;
      }
      return g_->nodes_[g_->i2u_[n1id_]].adj_[i].ev_;
    }

    /** Return this edge's value. */
    const node_value_type& value() const {
      size_type i = 0;
      for( ; i < node1().degree(); ++i) {
        if (g_->nodes_[g_->i2u_[n1id_]].adj_[i].nid_ == g_->i2u_[n2id_])
          break;
      }
      return g_->nodes_[g_->i2u_[n1id_]].adj_[i].ev_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1() && node2() == e.node2()) 
      	|| (node2() == e.node1() && node1() == e.node2()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(*this == e) {
      	//same edge
      	return false;
      }
      else {
        if (node1() < e.node1()) {
          return true;
        }
        else if (node2() < e.node2()) {
          return true;
        }
        else{
          return false;
        }
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* g_;
    size_type n1id_;
    size_type n2id_;

    Edge (const Graph* g, size_type n1id, size_type n2id) {
      g_ = const_cast<Graph*>(g);
      n1id_ = n1id;
      n2id_ = n2id;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type sum_deg = 0;
    for(unsigned i = 0; i < num_nodes(); ++i)
      sum_deg += node(i).degree();
    return sum_deg/2;

  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *ssss
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    EdgeIterator eit = edge_begin();
    size_type j;
    for (j = 0; j < i; j++){
      ++eit;
    }
    return *eit;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for(unsigned i = 0; i < a.degree(); ++i) {
      if(nodes_[i2u_[a.index()]].adj_[i].nid_ == i2u_[b.index()])
        return true;
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
    if (has_edge(a, b)){
      return Edge(this, a.ind_, b.ind_);      
    }
    else{
      nodes_[i2u_[a.ind_]].adj_.push_back({i2u_[b.ind_], edge_value_type()});
      nodes_[i2u_[b.ind_]].adj_.push_back({i2u_[a.ind_], edge_value_type()});
      return Edge(this, a.ind_, b.ind_);
    }
  }

   /** Remove an edge from the graph, connecting the specified nodes.
   *   Return size_type for success (1) or failure (0).
   * @param[in] a One end point of edge to remove
   * @param[in] a One end point of edge to remove
   * @post if edge exists, new num_edges() == old num_edges() - 1
   * @post if edge exists, new a.degree() = old a.degree() - 1
   * @post if edge exists, new b.degree() = old b.degree() - 1
   * @post has_edge(a,b) != true
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a,b)) {
      one_way_remove_edge(a, b);
      if (has_edge(b,a)) {
        one_way_remove_edge(b, a);   
        return size_type(1);     
      }
    }
    return size_type(0);
  }

  /** Remove an edge from the graph, connecting the specified nodes.
   *   Return size_type for success (1) or failure (0).
   * @param[in] a One end point of edge to remove
   * @param[in] a One end point of edge to remove
   * @post if edge exists, new num_edges() == old num_edges() - 1
   * @post has_edge(e.node1(),e.node2()) != true
   */
  size_type remove_edge(const Edge& e) {
    remove_edge(e.node1(), e.node2());
    return size_type(0);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
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
      g_ = nullptr;
      curind_ = size_type(-1);
    }

    /** Return object this iterator represents. */
    Node operator*() const {
      return g_->node(curind_);
    }

    /** Update this iterator to represent the next node in the graph. 
      * Return a reference to the resulting iterator. */
    NodeIterator& operator++() {
      if (curind_ < g_->num_nodes())
        curind_++;
      return *this;
    }

    /** Check if two node iterators are equivalent.
      * Return true if they are. */
    bool operator==(const NodeIterator& ni) const {
      return(g_ == ni.g_) && (curind_ == ni.curind_);
    }

   private:
    friend class Graph;
    Graph* g_;
    size_type curind_;

    NodeIterator(const Graph* g, size_type ind) {
      g_ = const_cast<Graph*>(g);
      curind_ = ind;
    }
  };

  /** Return an iterator starting over all nodes in the graph. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an iterator ending over all nodes in the graph. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  /** Remove a node from the graph, returning the index of the removed node.
   * @param[in] n_t The iterator representing position of node to remove
   * @post new num_nodes() == old num_nodes() - 1
   * @post Removes all edges incident to *n_it
   * @post Nodes node(i) are preserved for (*n_it).index() < i < num_nodes()
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type nid = remove_node(*n_it);
    return NodeIterator(this, nid);
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

    /** Return edge that the iterator represents. */
    Edge operator*() const {
      return Edge(g_, nid_, 
        g_->nodes_[g_->nodes_[g_->i2u_[nid_]].adj_[curind_].nid_].i_);
    }

    /** Increment the iterator.
      * Return a reference to the updated iterator. */
    IncidentIterator& operator++() {
      if (curind_ != g_->nodes_[g_->i2u_[nid_]].adj_.size())
        curind_++;
      return *this;
    }

    /** Check if two incident iterators are equivalent.
      * Return true if they are. */
    bool operator==( const IncidentIterator& ii) const {
      return(g_ == ii.g_) && (curind_ == ii.curind_) && (nid_ == ii.nid_);
    }

   private:
    friend class Graph;
    Graph* g_;
    size_type nid_;
    size_type curind_;

    IncidentIterator(const Graph* g, size_type n, size_type ind) {
      g_ = const_cast<Graph*>(g);
      nid_ = n;
      curind_ = ind;
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

    /** Return the edge that this iterator represents. */
    Edge operator*() const {
      if(curind_ < g_->nodes_[g_->i2u_[nid_]].adj_.size())
        return Edge(g_, nid_, 
          g_->nodes_[g_->nodes_[g_->i2u_[nid_]].adj_[curind_].nid_].i_);
      else
        return Edge();
    }

    /** Increment the edge iterator. Will not increment past end.
      * Return a reference to the updated iterator. */
    EdgeIterator& operator++() {
      do {
        if (curind_ == 0 && nid_ == g_->num_nodes())
          return *this; //end it
        curind_++;
        if (curind_ == (g_->nodes_[g_->i2u_[nid_]].adj_.size())) {
          curind_ = 0;
          nid_++;
        }
      }
      while (nid_ < g_->num_nodes() && g_->i2u_[g_->nodes_[g_->i2u_[nid_]].adj_[curind_].nid_] < g_->i2u_[nid_]);
      return *this;
    }

    /** Check if two edge iterators are equivalent. 
      * Return true if they are. */
    bool operator==(const EdgeIterator& ei) const {
      return(g_ == ei.g_) && (curind_ == ei.curind_) && (nid_ == ei.nid_);
    }

   private:
    friend class Graph;
    Graph* g_;
    size_type nid_;
    size_type curind_;

    EdgeIterator(const Graph* g, size_type nid, size_type ind) {
      g_ = const_cast<Graph*>(g);
      nid_ = nid;
      curind_ = ind;

      while (nid_ < g_->num_nodes() && g_->nodes_[g_->i2u_[nid_]].adj_.size() == 0)
      {
        nid_++;
        curind_ = 0;
      }
    }
  };

  /** Return an iterator starting over all edges in the graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }

  /** Return an iterator ending over all edges in the graph. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }

  /** Remove an edge from the graph, corresponding to  the specified nodes.
   *   Return edge_iterator representing first edge in graph.
   * @param[in] a One end point of edge to remove
   * @param[in] a One end point of edge to remove
   * @post if edge exists, new num_edges() == old num_edges() - 1
   * @post if edge exists, new a.degree() = old a.degree() - 1
   * @post if edge exists, new b.degree() = old b.degree() - 1
   * @post has_edge(a,b) != true
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return edge_begin();
  }

 private:
  /** RI: for all i in 0 < i < i2u_.size(), nodes_[i2u_[i]].i_ = i. */
  std::vector<size_type> i2u_;

  /** Helper function to remove edges. Undirected edges are internally 
  * stored as 2-way connections between nodes. This function removes a one-way
  * connection. Two calls to this function will remove the the undirected edge.
  */
  void one_way_remove_edge(const Node& a, const Node& b) {
    size_type i = 0;
    for( ; i < a.degree(); ++i) {
      if (nodes_[i2u_[a.index()]].adj_[i].nid_ == i2u_[b.index()])
        break;
    }
    //swap & pop back

    nodes_[i2u_[a.index()]].adj_[i] = nodes_[i2u_[a.index()]].adj_[a.degree()-1];
    nodes_[i2u_[a.index()]].adj_.pop_back();
  }

  //Class to hold information for incident edges (and the node endpoint)
  class AdjEdge {
    private:
      friend class Graph;
      size_type nid_;
      edge_value_type ev_;
      AdjEdge(const size_type id, const edge_value_type val) {
        nid_ = id;
        ev_ = val;
      }
  };
  //Class to represent Nodes internally
  class SecretNode {
    private:
      friend class Graph;
      Point p_;
      node_value_type v_;
      //invariant: i2u_[nodes_[i].index()] = i
      size_type i_;
      std::vector<AdjEdge> adj_;

      SecretNode(const Point& position, const node_value_type val, size_type i) {
        p_ = position;
        v_ = val;
        i_ = i;
      }
  };
};

#endif // CME212_GRAPH_HPP
