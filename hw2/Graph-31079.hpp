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
template <typename V, typename E>
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
  /** Synonum for E, the value of an edge */
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
      return g_->node_list[g_->i2un_[ind_]].pos;
    }

    Point& position() {
      return const_cast<Point&>(static_cast <Node const &>(*this).position());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return g_->node_list[g_->i2un_[ind_]].ind;
      return this->ind_;
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
      return g_->node_list[g_->i2un_[ind_]].val;
    }

    /** Setter for the value stored in node_value_type object */
    node_value_type& value() {
      return const_cast<node_value_type&>(static_cast <Node const &>(*this).value()); 
    }

    /** Return the degree of the node (number of incident edges */
    size_type degree() const {
      return g_->adj_list[g_->i2un_[ind_]].size();
    }

    /** Return incident_iterator to the first incident edge in adj_list[uid] */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, g_->i2un_[ind_], 0);
    }

    
    /** Return incident_iterator to the last incident edge in adj_list[uid] */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, g_->i2un_[ind_], (*this).degree());
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
    return i2un_.size();
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
    node_list.push_back(My_node(position, i2un_.size(), value));

    //DEBUG
    //std::cout << "Node ID map size: " << i2un_.size() << std::endl;
    //std::cout << "New node value is: " << value << std::endl;

    std::vector<std::tuple<size_type, size_type>> v {};
    adj_list.push_back(v);
    i2un_.push_back(node_list.size() - 1);
    return Node(this, i2un_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1). 
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return node_list[i2un_[n.index()]].pos == n.position()
        && node_list[i2un_[n.index()]].ind == n.ind();
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
      return Node(g_, g_->edge_list[g_->i2ue_[ind_]].a); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g_, g_->edge_list[g_->i2ue_[ind_]].b);
    }

    /** Getter for the value stored in edge_value_type object */
    const edge_value_type& value() const {
      return g_->edge_list[g_->i2ue_[ind_]].val;
    }

    /** Setter for the value stored in edge_value_type object */
    edge_value_type& value() {
      return const_cast<edge_value_type&>(static_cast <Edge const &>(*this).value()); 
    }
    
    /** Return the length of this Edge (Euclidean distance between nodes) */
    double length() const {
      return norm(node1().position() - node2().position());
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
      if (this->ind_ == e.ind_) {
        return this->g_ < e.g_;
      } else {  
        return this->ind_ < e.ind_;
      }
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
    return i2ue_.size();
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
    const std::vector<std::tuple<size_type, size_type>> &temp {adj_list[i2un_[a.index()]]};

    for (size_type i = 0; i < temp.size(); i++) {
      if (std::get<0>(temp[i]) == i2un_[b.index()]) {
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
    const std::vector<std::tuple<size_type, size_type>> &temp {adj_list[i2un_[a.index()]]};
    
    /* First look for existing edge */
    for (size_type i = 0; i < temp.size(); i++) {
      if (std::get<0>(temp[i]) == i2un_[b.index()]) {
        return edge(edge_list[std::get<1>(temp[i])].ind);
      }
    }

    /* If edge does not already exist, make a new edge */
    My_edge new_edge = My_edge(num_edges(), a.index(), b.index());
    edge_list.push_back(new_edge);
    i2ue_.push_back(edge_list.size() - 1);
    adj_list[i2un_[a.index()]].push_back(std::make_tuple(i2un_[b.index()], edge_list.size() - 1));
    adj_list[i2un_[b.index()]].push_back(std::make_tuple(i2un_[a.index()], edge_list.size() - 1));
    return edge(num_edges() - 1);
  }

  // DEBUG PURPOSES
  /** Returns the adjacency list */ 
  std::vector<std::vector<std::tuple<size_type, size_type>>> get_adj_list() {
    return adj_list;
  }

  /** Returns the node list */
  std::vector<My_node> get_node_list() {
    return node_list;
  }

  /** Returns the i2u node list */
  std::vector<size_type> get_i2u_nodes() {
    return i2un_;
  }

  /** Returns the edge list */
  std::vector<My_edge> get_edge_list() {
    return edge_list;
  }

  /** Returns the i2u edge list */
  std::vector<size_type> get_i2u_edges() {
    return i2ue_;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    i2un_.clear();
    edge_list.clear();
    i2ue_.clear();
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
      curr++;
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
      return g_->edge(g_->edge_list[std::get<1>(g_->adj_list[n_uid][curr])].ind);
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
       && (this->n_uid == rhs.n_uid);
    }
    
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* g_;
    size_type n_uid; // Unique index of current node
    size_type curr;  // Index within list of incident edges
    IncidentIterator(const Graph* g, size_type i, size_type n) : 
      g_(const_cast<Graph*>(g)), n_uid(i), curr(n) {
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
    size_type curr;       // Index of current Edge
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
  
  // HW2 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // size_type remove_node(const Node&);
  // node_iterator remove_node(node_iterator n_it);
  // size_type remove_edge(const Node&, const Node&);
  // size_type remove_edge(const Edge&);
  // edge_iterator remove_edge(edge_iterator e_it);  


  /** Remove an edge given its two nodes as input 
   *  @param[in] n1, n2  Two node objects
   *  @return  An index to i2ue_ where old edge was removed
   * 
   *  @pre   @a n1 is a valid node where 0 <= n1_uid < adj_list.size()
   *  @post  @new num_edges() == @old num_edges() - 1 (size of i2ue_ decreased by 1)
   *  @post  The size of edge_list is unchanged, as the edge uid isn't erased, just inaccessible
   *  @post  If edge is found, it is removed twice from adj_list,
   *        once under n1 and a second time under n2
   *  Complexity: O(num_nodes() + num_edges()) for adjacency list traversal
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    size_type n1_uid = i2un_[n1.index()];
    size_type n2_uid = i2un_[n2.index()];
    size_type ind = 0;

    for (unsigned k = 0; k < adj_list[n1_uid].size(); k++) {
      if (std::get<0>(adj_list[n1_uid][k]) == n2_uid) {
        ind = edge_list[std::get<1>(adj_list[n1_uid][k])].ind;
        size_type uid = i2ue_[ind];
        
        // We are removing the edge with 'last' index in i2ue_, so just pop
        if (ind == num_edges() - 1) {
          i2ue_.pop_back();
        } else {
          // We are removing from within i2ue_, so replace current ind with
          // final index in i2ue_ and update index of last edge accordingly
          i2ue_[ind] = i2ue_[num_edges() - 1];
          i2ue_.pop_back();
          edge_list[i2ue_[ind]].ind = ind;
        }

        // remove connection from adjacency list
        for (unsigned i = 0; i < adj_list.size(); ++i) {
          for (unsigned j = 0; j < adj_list[i].size(); ++j) {
            if (std::get<1>(adj_list[i][j]) == uid) {
              adj_list[i].erase(adj_list[i].begin() + j);
              --j;
            }
          }
        }
        break;
      }
    }
    return ind;
  }

  /** Remove an edge given that edge as input 
   *  @param[in] e  An Edge object
   *  @return  An index to i2ue_ where old edge was removed
   *
   *  @pre   Edge object has two Nodes that are part of this graph
   *  @post  This Edge @a e  can no longer be accessed by calling edge(@a e.index())
   *  Complexity: O(num_nodes() + num_edges()) as it calls the above remove_edge method
   */
  size_type remove_edge(const Edge& e) {
    //std::cout << "The two nodes have indices " << e.node1().index() << " and " << e.node2().index() << std::endl; // DEBUG
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge given an edge iterator as input 
   *  @param e_it  Iterator to Edge object
   *  @return  An EdgeIterator to the location where *(@a e_it) was removed
   *
   *  @pre  @a e_it is a valid EdgeIterator to this graph
   *  @post The edge *(@a e_it) is no longer accessible by edge(*(@a e_it)).index()
   *  Complexity: O(num_nodes() + num_edges()) as it calls the above remove_edge method
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type ind = remove_edge(*e_it);
    return EdgeIterator(this, ind);
  }

// DEBUG PURPOSES
size_type count_edges() {
  size_type count = 0;
  for (unsigned i = 0; i < this->adj_list.size(); i++) {
    count += this->adj_list[i].size();
  }
  return count;
}


  /** Remove a node given that node as input 
   *  @param[in] n  A Node object
   *  @return  An index to location in i2un where node was removed
   *
   *  @pre   @a n has index within bounds of i2un_ (0 <= n.index() < num_nodes())
   *  @post  This Node @a n can no longer be accessed by calling node(@a n.index())
   *  Complexity: O(num_nodes()) assuming the graph is sufficiently sparse
   */
  size_type remove_node(const Node& n) {
    size_type ind = n.index();
    size_type uid = i2un_[ind];
 
    // remove node information
    if (ind == num_nodes() - 1) { 
      // node has final index in i2un_ array, so just pop
      i2un_.pop_back();
    } else { 
      // we replace current index with final index and update Nodes accordingly
      i2un_[ind] = i2un_[num_nodes() - 1];
      i2un_.pop_back();
      node_list[i2un_[ind]].ind = ind;
    }

    // remove associated edges in adjacency list
    std::vector<std::tuple<size_type, size_type>> &incident_edges = adj_list[uid];
    while (!incident_edges.empty()) {
      std::tuple<size_type, size_type> &curr = incident_edges.back();
      size_type uid2 = std::get<0>(curr);
      size_type uide = std::get<1>(curr);
      
      // remove from neighbor adjacency list
      std::vector<std::tuple<size_type, size_type>> &other_adj_list = adj_list[uid2];
      for(unsigned i = 0; i < other_adj_list.size(); ++i) {
        if (std::get<0>(other_adj_list[i]) == uid) {
          other_adj_list.erase(other_adj_list.begin() + i);
          break;
        }
      }

      // renumber the indices if necessary
      size_type inde = edge_list[uide].ind;
        if (inde == num_edges() - 1) {
          i2ue_.pop_back();
        } else {
          i2ue_[inde] = i2ue_[num_edges() - 1];
          i2ue_.pop_back();
          edge_list[i2ue_[inde]].ind = inde;
        }
      
      incident_edges.pop_back();
    }

    return ind;
  }
  
  /** Remove a node given a node iterator as input
   *  @param n_it  Iterator to Node object
   *  @return  An NodeIterator to the location in i2un where node was removed
   *
   *  @pre  @a n_it is a valid NodeIterator to this graph
   *  @post The node *(@a n_it) is no longer accessible by node(*(@a n_it)).index()
   *  Complexity: O(num_nodes()) as it calls the above remove_node method
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type ind = remove_node(*n_it);
    return NodeIterator(this, ind);
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
   *  @param val   Value stored by the edge (template type)
   */
  struct My_edge {
    size_type ind;
    size_type a;
    size_type b;
    edge_value_type val;
    // Constructor
    My_edge(size_type i, size_type first, size_type second) :
      ind(i), a(first), b(second) {
    }
  };

  std::vector<My_node> node_list; // Store all My_node objects of the graph
                                  // Indexed by node uid
  std::vector<size_type> i2un_;   // Indexed by node ind
  std::vector<My_edge> edge_list; // Store all My_edge objects of the graph
                                  // Indexed by edge uid
  std::vector<size_type> i2ue_;   // Indexed by edge ind
  /** Adjacency list representing the connections in the graph
   *  adj_list[i] contains a vector of tuples (j, k) such that edge(k)
   *  connects node(i) to node(j).
   */
  std::vector<std::vector<std::tuple<size_type, size_type>>> adj_list;
};

#endif  // CME212_GRAPH_HPP
