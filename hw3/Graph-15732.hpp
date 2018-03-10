#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#include <iostream>

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

   using node_value_type = V;
   using edge_value_type = E;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node {
    Point position_;
    node_value_type value_;
  };

  struct internal_edge {
    unsigned n_uid1_, n_uid2_;
    edge_value_type value_;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<std::vector<std::pair<unsigned, unsigned> > > adj_;

  // The vectors nodes_ and edges_ hold internal_node and internal_edge structs
  //   containing node_value_type and edge_value_type variables, respectively.
  //   Since these variables could be anything and any size, we want to avoid
  //   removing elements from these vectors to not have to make copies of them,
  //   which potentially could be very expensive. Instead, the vectors nodes_i2u_
  //   and edges_i2u_ hold the uid's pointing to the non-removed elements of the
  //   nodes_ and edges_ vectors. The expenses of deleting from teh i2u_ vectors
  //   are limited because those vectors only contain unsigned variables.
  std::vector<unsigned> nodes_i2u_;
  std::vector<unsigned> edges_i2u_;

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

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

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
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODE ITERATOR
  //

  /** @class Graph::node_iterator
   * @brief Class representing iterator on the graph's nodes.
   *
   * node_iterator objects are used to access information about the Graph's nodes.
   */
  class node_iterator : private totally_ordered<node_iterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    node_iterator() {
    }

    // HW1 #2: YOUR CODE HERE

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    // Dereference operator, returns the node belonging to the node_iterator.
    Node operator*() const {
      return Node(g_, n_idx_);
    }

    // Increment operator, increasing the node id.
    node_iterator& operator++() {
      ++n_idx_;
      return *this;
    }

    // Smaller operator, comparing the node ids of two node_iterators.
    bool operator<(const node_iterator& nit) const {
      if (g_ != nit.g_) {
        if (g_ < nit.g_) return true;
        return false;
      }
      if (n_idx_ < nit.n_idx_) return true;
      return false;
    }

    // Equal operator, comparing graphs and node ids of two iterators.
    bool operator==(const node_iterator& nit) const {
      if (g_ == nit.g_ && n_idx_ == nit.n_idx_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // The private member variables define the graph and the node id.
    const Graph* g_;
    unsigned n_idx_;

    // Construct a valid node_iterator.
    node_iterator(const Graph* g, unsigned n_idx) : g_(g), n_idx_(n_idx) {}
  };

  // Returns node iterator pointing to first node of graph.
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  // Returns node iterator pointing to last node of graph.
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }

  //
  // EDGE ITERATOR
  //

  class edge_iterator : private totally_ordered<edge_iterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    edge_iterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    // Dereference operator, returns the edge belonging to the edge_iterator.
    Edge operator*() const {
      unsigned e_uid_  = g_->edges_i2u_[e_idx_];
      unsigned n_idx1_ = g_->n_u2i(g_->edges_[e_uid_].n_uid1_);
      unsigned n_idx2_ = g_->n_u2i(g_->edges_[e_uid_].n_uid2_);
      return Edge(g_, e_uid_, n_idx1_, n_idx2_);
    }

    // Increment operator, increasing the edge id.
    edge_iterator& operator++() {
      ++e_idx_;
      return *this;
    }

    // Smaller operator, comparing the edge ids of two edge_iterator.
    bool operator<(const edge_iterator& eit) const {
      if (g_ != eit.g_) {
        if (g_ < eit.g_) return true;
        return false;
      }
      if (e_idx_ < eit.e_idx_) return true;
      return false;
    }

    // Equal operator, comparing graphs and edge ids of two iterators.
    bool operator==(const edge_iterator& eit) const {
      if (g_ == eit.g_ && e_idx_ == eit.e_idx_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    // The private member variables define the graph and the edge id.
    const Graph* g_;
    unsigned int e_idx_;

    // Construct a valid edge_iterator.
    edge_iterator(const Graph* g, unsigned e_idx) : g_(g), e_idx_(e_idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  // Returns the edge_iterator pointing to first edge of graph.
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  // Returns the edge_iterator pointing to last edge of graph.
  edge_iterator edge_end() const {
    return edge_iterator(this, edges_i2u_.size());
  }

  //
  // INCIDENT ITERATOR
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class incident_iterator : private totally_ordered<incident_iterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    incident_iterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    // Dereference operator, returns the edge belonging to the incident_iterator.
    Edge operator*() const {
      unsigned e_uid_  = g_->adj_[root_][iid_].second;
      // unsigned n_idx1_ = g_->n_u2i(g_->edges_[e_uid_].n_uid1_);
      // unsigned n_idx2_ = g_->n_u2i(g_->edges_[e_uid_].n_uid2_);
      unsigned n_idx1_ = g_->n_u2i(root_);
      unsigned n_idx2_ = g_->n_u2i(g_->adj_[root_][iid_].first);
      return Edge(g_, e_uid_, n_idx1_, n_idx2_);
    }

    // Increment operator, increasing the incident id.
    incident_iterator& operator++() {
      ++iid_;
      return *this;
    }

    // Decrement operator, decreasing the incident id.
    incident_iterator& operator--() {
      if (iid_>0) --iid_;
      return *this;
    }

    // Minus operator, decreasing the incident id by the int argument.
    incident_iterator& operator-(unsigned i) {
      if (iid_ > i) {
        iid_ -= i;
      } else {
        iid_ = 0;
      }
      return *this;
    }

    // Smaller operator, comparing the incident ids of two incident_iterators.
    bool operator<(const incident_iterator& iit) const {
      if (g_ != iit.g_) {
        if (g_ < iit.g_) return true;
        return false;
      }
      if (root_ != iit.root_) {
        if (root_ < iit.root_) return true;
        return false;
      }
      if (iid_ < iit.iid_) {
        return true;
      }
      return false;
    }

    // Equal operator, comparing graphs and incident ids of two iterators.
    bool operator==(const incident_iterator& iit) const {
      if (g_ == iit.g_ && root_ == iit.root_ && iid_ == iit.iid_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    // The private member variables define the graph, the root,
    //  and the incident id.
    const Graph* g_;
    unsigned int iid_, root_;

    // Construct a valid incident_iterator.
    incident_iterator(const Graph* g, unsigned iid, unsigned root)
        : g_(g), iid_(iid), root_(root) {}
  };


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
      // HW0: YOUR CODE HERE
    }

    // Number of incident edges from the adjacency vector.
    size_type degree() const {
      return g_->adj_[g_->nodes_i2u_[n_idx_]].size();
    }

    // Creates an incident_iterator for the first incident edge.
    incident_iterator edge_begin() const {
      return incident_iterator(g_, 0, g_->nodes_i2u_[n_idx_]);
    }

    // Creates an incident_iterator for the last incident edge.
    incident_iterator edge_end() const {
      return incident_iterator(g_, g_->adj_[g_->nodes_i2u_[n_idx_]].size(), g_->nodes_i2u_[n_idx_]);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->nodes_[g_->nodes_i2u_[n_idx_]].position_;
    }

    /** Return this node's position in a modifiable way. */
    Point& position() {
      return g_->nodes_[g_->nodes_i2u_[n_idx_]].position_;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return g_->nodes_[g_->nodes_i2u_[n_idx_]].value_;
    }

    /** Return this node's value in a modifiable way. */
    node_value_type& value() {
      return g_->nodes_[g_->nodes_i2u_[n_idx_]].value_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return n_idx_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (g_ == n.g_ && n_idx_ == n.n_idx_) {
        return true;
      }
      return false;
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
      if (g_ != n.g_) {
        if (g_ < n.g_) return true;
        return false;
      }
      if (n_idx_ < n.n_idx_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* g_;
    unsigned n_idx_;
    Node(const Graph* g, unsigned int n_idx)
        : g_(const_cast<Graph*>(g)), n_idx_(n_idx) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }


  /*Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    internal_node new_internal_node_;
    new_internal_node_.position_ = position;
    Node          new_node_      = Node(this, size());
    nodes_.push_back(new_internal_node_);
    std::vector<unsigned> adj_node_;
    adj_.push_back(adj_node_);
    return new_node_;
  }*/

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {

    // Create new internal node and new entry for adjacency vector.
    internal_node new_internal_node_;
    new_internal_node_.position_ = position;
    new_internal_node_.value_    = value;
    Node new_node_               = Node(this, nodes_i2u_.size());
    std::vector<std::pair<unsigned, unsigned> > adj_node_;

    // Push back to nodes_, nodes_i2u_, and adj_ vector. Then return the new node instance.
    nodes_i2u_.push_back(nodes_.size());
    nodes_.push_back(new_internal_node_);
    adj_.push_back(adj_node_);
    return new_node_;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.g_ && size() > n.n_idx_) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    Node new_node = Node(this, i);
    return new_node;
  }


  // Return the current node index of the node with the passed node uid.
  size_type n_u2i(const size_type n_uid_) const {
    // Assume that not super many nodes have been deleted. Then it is likely
    //   much more efficient to search backwards from the uid position (or the
    //   back of the vector if that comes earlier).
    size_type n_idx_ = std::min(n_uid_, num_nodes()-1);
    while (nodes_i2u_[n_idx_] != n_uid_) --n_idx_;
    return n_idx_;
  }


  /** Remove a node and all edges incident to it from the graph and return the
   *    new graph size.
   * @param[in] n The node that is to be removed.
   * @pre @n is a node of the graph.
   * @return The graph size after the removal of the node.
   * @post The node @n is invalid.
   * @post All nodes with indices (that's the n_idx_ member variable) larger
   *         larger than the one of the removed node are invalid.
   * @post All edges incident to the removed node @n have been removed as well
   *         from the graph.
   * @post Can invalidate edge indices (that's the e_idx_ member variables).
   *
   * Complexity: Complexity: O(num_edges()*degree + degree^3)
   */
  size_type remove_node(const Node& n) {

    // Removing all edges incident to the node n.
    for (auto iit_=(n.edge_end()); iit_>n.edge_begin(); ) {
      remove_edge(*(--iit_));
    }

    // Removing the node n itself and returning the number of nodes in the graph.
    nodes_i2u_.erase(nodes_i2u_.begin()+n.index());
    return nodes_i2u_.size();
  }


  /** Remove the node associated with the passed node_iterator and all edges
   *    incident to it from the graph and return the node_iterator associated
   *    with the next node.
   * @param[in] n_it The node_iterator to the node n that is to be removed.
   * @pre The associated node n is a node of the graph.
   * @return The node_iterator associated with the next node, or the end
   *           iterator if there is no next node.
   * @post The node n is removed.
   * @post All nodes with indices (that's the n_idx_ member variable) larger
   *         larger than the one of the removed node are invalid.
   * @post All edges incident to the removed node @n have been removed as well
   *         from the graph.
   * @post Can invalidate edge indices (that's the e_idx_ member variables).
   *
   * Complexity: Complexity: O(num_edges()*degree + degree^3)
   */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it != node_end()) {
      remove_node(*n_it);
    }
    return n_it;
  }

  /** Remove the edge between the two passed nodes from the graph.
   * @param[in] a, b The two nodes between which we want to remove the edge.
   * @pre @a and @b are nodes of the graph.
   * @return The number of edges in the graph after removing the edge.
   * @post The edge is removed: It is not in the adjacency vector adj_ and it
   *         is not in the edges_i2u_ vector.
   * @post All edges with indices (that's the e_idx_ member variable) larger
   *         larger than the one of the removed edge are invalid.
   *
   * Complexity: O(num_edges() + degree^2)
   */
  size_type remove_edge(const Node& a, const Node& b) {

    // Only try to remove the edge if it exists.
    if (has_edge(a, b)) {

      // Get the uid's of the two connected nodes.
      auto n_uid1_ = nodes_i2u_[a.index()];
      auto n_uid2_ = nodes_i2u_[b.index()];

      // Remove the connection from the adjacency vector, remembering the edge uid.
      unsigned e_uid_;
      for (auto it_=adj_[n_uid1_].begin(); it_!=adj_[n_uid1_].end(); ++it_) {
        if ((*it_).first == n_uid2_) {
          e_uid_ = (*it_).second;
          adj_[n_uid1_].erase(it_);
          break;
        }
      }
      for (auto it_=adj_[n_uid2_].begin(); it_!=adj_[n_uid2_].end(); ++it_) {
        if ((*it_).first == n_uid1_) {
          adj_[n_uid2_].erase(it_);
          break;
        }
      }

      // Remove the edge from the edges_i2u_ container.
      auto it_ = edges_i2u_.begin();
      while (it_ != edges_i2u_.end()) {
        if (*it_ == e_uid_) {
          break;
        }
        ++it_;
      }
      edges_i2u_.erase(it_);
    }

    // // Return the container size.
    return edges_i2u_.size();
  }

  /** Remove the edge from the graph.
   * @param[in] e The edge to be removed from the graph.
   * @pre @e is an edge of the graph.
   * @return The number of edges in the graph after removing the edge.
   * @post The edge is removed: It is not in the adjacency vector adj_ and it
   *         is not in the edges_i2u_ vector.
   * @post All edges with indices (that's the e_idx_ member variable) larger
   *         larger than the one of the removed edge are invalid.
   *
   * Complexity: O(num_edges() + degree^2)
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  /** Remove the edge associated with the passed edge_iterator from the graph.
   * @param[in] e_it The edge_iterator associated with the edge to be removed
   *                   from the graph.
   * @pre @e_it is a valid edge_iterator of the graph.
   * @return The edge_iterator pointing to the next edge.
   * @post The edge is removed: It is not in the adjacency vector adj_ and it
   *         is not in the edges_i2u_ vector.
   * @post All edges with indices (that's the e_idx_ member variable) larger
   *         larger than the one of the removed edge are invalid.
   *
   * Complexity: O(num_edges() + degree^2)
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end()) {
      remove_edge(*e_it);
    }
    return e_it;
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(g_, n_idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g_, n_idx2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (g_ == e.g_ && e_uid_ == e.e_uid_) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (g_ != e.g_) {
        if (g_ < e.g_) return true;
        return false;
      }
      if (e_uid_ < e.e_uid_) return true;
      return false;
    }

    // Return edge length.
    double length() const {
      return norm(node1().position() - node2().position());
    }


    // Returning the edge value of type edge_value_type.
    edge_value_type& value() {
      return g_->edges_[e_uid_].value_;
    }

    // Returning the edge value of type edge_value_type.
    const edge_value_type& value() const {
      return g_->edges_[e_uid_].value_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* g_;
    unsigned e_uid_, n_idx1_, n_idx2_;

    Edge(const Graph* g, unsigned e_uid, unsigned n_idx1, unsigned n_idx2)
        : g_(const_cast<Graph*>(g)), e_uid_(e_uid), n_idx1_(n_idx1), n_idx2_(n_idx2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    unsigned e_uid_  = edges_i2u_[i];
    unsigned n_idx1_ = n_u2i(edges_[e_uid_].n_uid1_);
    unsigned n_idx2_ = n_u2i(edges_[e_uid_].n_uid2_);
    return Edge(this, e_uid_, n_idx1_, n_idx2_);
    // Edge return_edge = Edge(this, e_uid_, edges_[e_uid_].n_idx1_, edges_[e_uid_].n_idx2_);
    // return return_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // std::vector<std::pair<unsigned, unsigned> > adj_a_ = adj_[nodes_i2u_[a.nid_]];
    // unsigned idx(0), adj_a_size_(adj_a_.size());

    auto adj_a_ = adj_[nodes_i2u_[a.index()]];
    unsigned n_uidb_ = nodes_i2u_[b.index()];
    for (unsigned i=0; i<adj_a_.size(); ++i) {
      if (adj_a_[i].first == n_uidb_) return true;
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
    // HW0: YOUR CODE HERE

    // Check if edge already exists, and if so return it.
    for (auto iit_=a.edge_begin(); iit_!=a.edge_end(); ++iit_) {
      auto c = (*iit_).node2();
      if (b==c) {
        return (*iit_);
      }
    }

    // Create a new edge and return it.
    unsigned e_uid_   = edges_.size();
    std::pair<unsigned, unsigned> pb_pair_;
    pb_pair_.first    = nodes_i2u_[b.n_idx_];
    pb_pair_.second   = e_uid_;
    adj_[nodes_i2u_[a.n_idx_]].push_back(pb_pair_);
    pb_pair_.first    = nodes_i2u_[a.n_idx_];
    adj_[nodes_i2u_[b.n_idx_]].push_back(pb_pair_);

    internal_edge e_;
    e_.n_uid1_ = nodes_i2u_[a.n_idx_];
    e_.n_uid2_ = nodes_i2u_[b.n_idx_];
    edges_i2u_.push_back(e_uid_);
    edges_.push_back(e_);
    return Edge(this, e_uid_, a.n_idx_, b.n_idx_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    adj_.clear();
    nodes_i2u_.clear();
    edges_i2u_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
