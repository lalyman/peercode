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
template <typename V>
class Graph {
 public:
 
   using node_value_type = V;
   
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
  std::vector<internal_node> nodes_;
  std::vector<std::vector<std::pair<unsigned, unsigned> > > adj_;
  std::vector<std::pair<unsigned, unsigned> > edges_;

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
      return Node(g_, nid_);      
    }
    
    // Increment operator, increasing the node id.
    node_iterator& operator++() {
      ++nid_;
      return *this;
    }
    
    // Smaller operator, comparing the node ids of two node_iterators.
    bool operator<(const node_iterator& nit) const {
      if (g_ == nit.g_ && nid_ < nit.nid_) {
        return true;
      }
      return false;
    }
    
    // Equal operator, comparing graphs and node ids of two iterators.
    bool operator==(const node_iterator& nit) const {
      if (g_ == nit.g_ && nid_ == nit.nid_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    // The private member variables define the graph and the node id.
    const Graph* g_;
    unsigned int nid_;
    
    // Construct a valid node_iterator.
    node_iterator(const Graph* g, unsigned nid) : g_(g), nid_(nid) {}
  };
  
  // Returns node iterator pointing to first node of graph.
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }
  
  // Returns node iterator pointing to last node of graph.
  node_iterator node_end() const {
    return node_iterator(this, size());
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
      return Edge(g_, eid_, g_->edges_[eid_].first, g_->edges_[eid_].second);      
    }
    
    // Increment operator, increasing the edge id.
    edge_iterator& operator++() {
      ++eid_;
      return *this;
    }
    
    // Smaller operator, comparing the edge ids of two edge_iterator.
    bool operator<(const edge_iterator& eit) const {
      if (g_ == eit.g_ && eid_ < eit.eid_) {
        return true;
      }
      return false;
    }
    
    // Equal operator, comparing graphs and edge ids of two iterators.
    bool operator==(const edge_iterator& eit) const {
      if (g_ == eit.g_ && eid_ == eit.eid_) {
        return true;
      }
      return false;
    }
    
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    // The private member variables define the graph and the edge id.
    const Graph* g_;
    unsigned int eid_;
    
    // Construct a valid edge_iterator.
    edge_iterator(const Graph* g, unsigned eid) : g_(g), eid_(eid) {}
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
    return edge_iterator(this, edges_.size());
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
      return Edge(g_, g_->adj_[root_][iid_].second, root_, g_->adj_[root_][iid_].first);
    }
    
    // Increment operator, increasing the incident id.
    incident_iterator& operator++() {
      ++iid_;
      return *this;
    }
    
    // Smaller operator, comparing the incident ids of two incident_iterators.
    bool operator<(const incident_iterator& iit) const {
      if (g_ == iit.g_ && root_ == iit.root_ && iid_ < iit.iid_) {
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
      return adj_[nid_].size();
    }
    
    // Creates an incident_iterator for the first incident edge.
    incident_iterator edge_begin() const {
      return incident_iterator(g_, 0, nid_);
    }
  
    // Creates an incident_iterator for the last incident edge.
    incident_iterator edge_end() const {
      return incident_iterator(g_, g_->adj_[nid_].size(), nid_);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->nodes_[nid_].position_;
    }
    
    node_value_type& value() {
      return g_->nodes_[nid_].value_;
    }
    
    const node_value_type& value() const {
      return g_->nodes_[nid_].value_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (g_ == n.g_ && nid_ == n.nid_) {
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
      if (nid_ < n.nid_) {
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
    unsigned int nid_;
    Node(const Graph* g, unsigned int nid)
        : g_(const_cast<Graph*>(g)), nid_(nid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_internal_node_;
    new_internal_node_.position_ = position;
    new_internal_node_.value_    = value;
    Node new_node_ = Node(this, size());
    nodes_.push_back(new_internal_node_);
    std::vector<std::pair<unsigned, unsigned> > adj_node_;
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
    if (this == n.g_ && size() > n.nid_) {
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
      return Node(g_, nid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g_, nid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (g_ == e.g_ && eid_ == e.eid_) {
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
      if (eid_ < e.eid_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* g_;
    unsigned eid_, nid1_, nid2_;

    Edge(const Graph* g, unsigned eid, unsigned nid1, unsigned nid2)
        : g_(const_cast<Graph*>(g)), eid_(eid), nid1_(nid1), nid2_(nid2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    Edge return_edge = Edge(this, i, edges_[i].first, edges_[i].second);
    return return_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::vector<std::pair<unsigned, unsigned> > adj_a_ = adj_[a.nid_];
    unsigned idx(0), adj_a_size_(adj_a_.size());
    while (idx < adj_a_size_) {
      if (adj_a_[idx].first == b.nid_) {
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    std::vector<std::pair<unsigned, unsigned> > adj_a_ = adj_[a.nid_];
    unsigned idx(0), adj_a_size_(adj_a_.size());
    while (idx < adj_a_size_) {
      if (adj_a_[idx].first == b.nid_) {
        Edge return_edge_ = Edge(this, adj_a_[idx].second, a.nid_, b.nid_);
        return return_edge_;
      }
      ++idx;
    }
    std::pair<unsigned, unsigned> pb_pair_;
    pb_pair_.first    = b.nid_;
    pb_pair_.second   = num_edges();
    adj_[a.nid_].push_back(pb_pair_);
    pb_pair_.first    = a.nid_;
    adj_[b.nid_].push_back(pb_pair_);
    pb_pair_.second   = b.nid_;
    edges_.push_back(pb_pair_);
    Edge return_edge_ = Edge(this, num_edges(), a.nid_, b.nid_);
    return return_edge_;
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
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
