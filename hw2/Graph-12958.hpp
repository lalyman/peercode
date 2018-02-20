#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
//HW2_CJ - HW1 update
template <typename V, typename E>
class Graph {
 private:

  /*** Predeclare internal structs for node and edge ***/
  struct internal_node;
  //struct internal_edge;

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

  //HW1_CJ
  using node_value_type = V;
  //HW2_CJ
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    //HW2_CJ - HW1 update
    : nodes_(), adj_(), idx2uid_() { 
   /*** Graph constructed with zero size ***/
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
  class Node: private totally_ordered<Node> {
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

    /* Creates an invalid Node object */
    Node() {
      /*** Invalid Node; must be created by calling Graph methods  ***/
    }

    /*** Return this node's position. ***/
    const Point& position() const {
      /*** access this node's position (via its uid) through graph_ ptr ***/
      return graph_->nodes_[node_uid_].position;
    }

    //HW2_CJ
    /*** Return modifiable Node position***/
    Point& position() {
      /*** access this node's position (via its uid) through graph_ ptr ***/
      return graph_->nodes_[node_uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      /*** access this node's index (uid) through graph_ ptr ***/
      return this->node_uid_;
    }

    //HW1_CJ
    /** @brief Return this Node's value
    *   @post return internal_node.node_val data member; 
    */ 
    node_value_type& value() {
      return graph_->nodes_[node_uid_].node_val;
    }
    /** @brief Return this Node's value
    *   @post return internal_node.node_val data member; 
    */ 
    const node_value_type& value() const {
      return graph_->nodes_[node_uid_].node_val;
    }
    /** @brief Return this Node's degree
    *   @pre 0 <= @a i <= num_edges()
    *   @post result_node.degree() == i
    *   Degree represents number of edges adjacent to the node
    */ 
    size_type degree() const {
      size_type node_idx = this->node_uid_;
      return graph_->adj_[node_idx].size();
    }
    /** @brief Return first iterator for Edge element 
    *   @pre 0 <= @a i <= num_edges()
    *   @pre must be implemented on Graph::Node object
    *   @post node.edge_begin() == 0
    */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, node_uid_, 0);
    }

    /** @brief Return last iterator for Edge element 
    *   @pre 0 <= @i <= num_edges()
    *   @pre 0 <= node_uid_ < size()
    *   @pre must be implemented on Graph::Node object
    *   @post node.edge_end() == adj_[node_uid_].size()
    *   Use adjacency list to access a given node in the adj_ list
    */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, node_uid_, graph_->adj_[node_uid_].size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /*** Test if this node's ID and graph ptr are equivalent to n's ***/
      if ((n.index() == this->node_uid_) && (n.graph_ == this->graph_)) {
        return true;
      }
      else {
        return false;
      }
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
      /*** Test if nodes are in same graph ***/
      if (n.graph_ != this->graph_) {
        //implement global ordering
        return this->graph_ < n.graph_;
      }
      /*** Compare node uids for ordering purposes ***/
      if (n.index() > this->node_uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    /*** Private data members and methods for Node ***/
    //Pointer back to the Graph container
    Graph* graph_;
    //This node's unique identification number
    size_type node_uid_;
    /* Private constructor for Graph::Node object */
    Node(const Graph* graph, size_type node_uid) 
      : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid) {
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    /*** size_ == num_nodes() ***/
    //HW2_CJ
    return this->nodes_.size();
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

  /*** Add new internal_node to nodes_ vector of internal_nodes ***/
  Node add_node(const Point& position, const node_value_type& new_node_value = node_value_type()) {
    //HW2_CJ
    //instantiate an internal node to add to nodes_ vector
    internal_node new_node;
    //assign data attributes to new_node
    new_node.position = position;
    new_node.node_val = new_node_value; //HW1_CJ
    new_node.index = nodes_.size(); //HW2_CJ
    //append new_node to nodes_ vector
    nodes_.push_back(new_node);
    idx2uid_.push_back(new_node.index); //HW2_CJ
    /*** update adj_ list ***/
    adj_.emplace_back(std::vector<std::pair<size_type, edge_value_type>>());
    //update number of nodes and next node ID 
    return Node(this, this->nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /*** same graph if ptr to memory for node is the same as graph_ ptr 
         AND if n's index is less than num_nodes() of graph ***/ 
    if ((this == n.graph_) && (num_nodes() > n.index())) {
      return true;
    }
    else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    /*** return a proxy object for node @a i ***/
    assert(i < num_nodes());
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      /*** Invalid Edge; must be created by calling Graph methods ***/
    }

    /** Return a node of this Edge */
    Node node1() const {
      /*** access this edge's node1 member through graph_ ptr ***/
      return graph_->node(node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      /*** access this edge's node2 member through graph_ ptr ***/
      return graph_->node(node2_uid_); 
    }

    //HW2_CJ
    /** @brief Return this Edge's values
    *   @pre match_idx < adj_[node1.index()].size()
    *   @pre edge.graph_ == this.graph_
    *   @post return edge value stored in sub-vector pair data member 
    */
    edge_value_type& value() {
      Node node1 = this->node1();
      Node node2 = this->node2();
      unsigned int match_idx = 0;
      //implement IncidentIterator to iterate though adjacent edges
      for (auto iit = node1.edge_begin(); iit != node1.edge_end(); ++iit) {
        //test if node uids match current Edge
        Edge e = *iit;
        if (node2 == e.node1()|| node2 == e.node2()) {
          return (graph_->adj_[node1.index()][match_idx]).second;
        }
        ++match_idx;
      }
      return (graph_->adj_[node1.index()][match_idx]).second;
    }

    /** @brief Return this Edge's values
    *   @pre match_idx < adj_[node1.index()].size()
    *   @pre edge.graph_ == this.graph_
    *   @post return edge value stored in sub-vector pair data member 
    */
    const edge_value_type& value() const {
      Node node1 = this->node1();
      Node node2 = this->node2();
      unsigned int match_idx = 0;
      //implement IncidentIterator to iterate though adjacent edges
      for (auto iit = node1.edge_begin(); iit != node1.edge_end(); ++iit) {
        //test if node uids match current Edge
        Edge e = *iit;
        if (node2 == e.node1()|| node2 == e.node2()) {
          return (graph_->adj_[node1.index()][match_idx]).second;
        }
        ++match_idx;
      }
      return (graph_->adj_[node1.index()][match_idx]).second;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //this Edge and Edge e must have same nodes
      if (((e.node1() == this->node1()) && (e.node2() == this->node2())) 
          || ((e.node1() == this->node2()) && (e.node2() == this->node1()))) {
        //this Edge and Edge e must be in same graph
        if (e.graph_ == this->graph_) {
          return true;
        } 
        else {
          return false;
        }
      }  
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      /*** Test if Edges are in same Graph ***/
      if (e.graph_ != this->graph_) {
        //implement global ordering
        return this->graph_ < e.graph_;
      }
      /*** Compare node uids within Edge for ordering purposes ***/
      if (e.node1() > this->node1() || (this->node1() == e.node1() && this->node2() < e.node2())) {
        return true;
      }
      else {
        return false;
      }
    }
    
    //HW2_CJ
    /* @brief Determine the Euclidean distance between two Nodes in an Edge
    *  @pre node1.graph_ == node2.graph_
    *  @pre has_edge(node1, node2) == true
    *  @pre node1.degree() > 1
    *  @pre node2.degree() > 1
    */
    double length() const {
      //get both Nodes for this Edge
      auto node1_position = node1().position();
      auto node2_position = node2().position();
      //difference between each Cartesian coordinate for nodes
      double del_x = node1_position.x - node2_position.x;
      double del_y = node1_position.y - node2_position.y;
      double del_z = node1_position.z - node2_position.z;
      //Euclidean distance
      double length = sqrt(pow(del_x,2) + pow(del_y,2) + pow(del_z,2));
      return length;
    }

   //HW2_CJ
   /** Set the value of this Edge to a given input parameter 
   *   @pre idx < adj_[node1_uid_].size()
   *   @pre idx < adj_[node2_uid_].size()
   *   result: Edge's edge_val data member is updated for both Nodes in Edge
   */
   void set_edge_value(edge_value_type edge_val) {
     //set edge value for first node in Edge pair
     for (unsigned int idx = 0; idx < graph_->adj_[node1_uid_].size(); ++idx) {
       //compare node uids
       if ((graph_->adj_[node1_uid_][idx]).first == node2_uid_) {
         (graph_->adj_[node1_uid_][idx]).second = edge_val;
       }
     }
     //set edge value for other node in Edge pair
     for (unsigned int idx = 0; idx < graph_->adj_[node2_uid_].size(); ++idx) {
       //compare node uids
       if ((graph_->adj_[node2_uid_][idx]).first == node1_uid_) {
         (graph_->adj_[node2_uid_][idx]).second = edge_val;
       }
     }
   }

   private:
    /*** Private data members and methods for Edge ***/
    // Pointer back to the Graph container
    Graph* graph_;
    //HW2_CJ - HW1 update
    /* Unique IDs of Nodes connected by this Edge */
    size_type node1_uid_;
    size_type node2_uid_;
    /* private constructor Graph::Edge */
    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
      : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  //HW2_CJ - HW1 update
  size_type num_edges() const {
    ///calculate total number of connections in adjacency list
    unsigned int total_connects = 0;
    for (auto it = adj_.begin(); it != adj_.end(); ++it) {
      total_connects += (*it).size();
    }
    //return number of unique connections (no double counting)
    return total_connects/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  //HW2_CJ - HW1 update
  Edge edge(size_type i) const {
    /*** return a proxy object for edge @ i ***/
    assert(i < num_edges());
    auto ith_edge = std::next(edge_begin(), i);
    return *ith_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    /*** Verify that a and b nodes are witin this graph_ ***/
    if (this == a.graph_ && this == b.graph_) {
      //HW2_CJ - HW1 update
      for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
        //Check if two nodes are connected by Edge via adj_ list
        auto edge = *it;
        if (edge.node1() == b || edge.node2() == b) {
          return true;
        } 
      }
    }
    //return false if a and b are not within same Edge
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
  //HW2_CJ - HW1 update
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& new_edge_value = edge_value_type()) {
    //test if edge already exists via has_edge()
    if (has_edge(a,b)) {
      //return same edge
      return Edge(this, a.index(), b.index());
    }
    //Add edge if not already existent
    /*** Update adj_ list corresponding to both nodes ***/
    //add pair of (connecting nodes, corresponding edge value)
    adj_[a.node_uid_].push_back(std::make_pair(b.node_uid_, new_edge_value));
    adj_[b.node_uid_].push_back(std::make_pair(a.node_uid_, new_edge_value));

    return Edge(this, a.index(), adj_[a.index()].size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adj_.clear();
    idx2uid_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  //HW1_CJ
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    //HW1_CJ
    /** @brief Dereference current NodeIterator
    *   @pre 0 <= nit_index_ <= num_nodes()
    *   @post Node.graph_ == this->graph_ 
    */ 
    Node operator*() const {
      //return Node corresponding to current iterator position
      return graph_->node(nit_index_);  
    }

    /** @brief Increment current NodeIterator
    *   @pre 0 <= ++nit_index_ <= num_nodes()
    *   @post Node.graph_ == this->graph_ 
    *   @post ++nit_index_ != nullptr
    */ 
    NodeIterator& operator++() {
      ++nit_index_;
      return *this;
    }

    /** @brief Evaluate if current NodeIterator is equal to another
    *   @pre this->graph_ == niter.graph_
    *   @pre 0 <= this->nit_index_ < num_nodes()
    *   @pre 0 <= niter.nit_index_ < num_nodes()
    *   @post niter.graph_ == this->graph_ 
    */ 
    bool operator==(const NodeIterator& niter) const {
      if (this->graph_ == niter.graph_) {
        if (this->nit_index_ == niter.nit_index_) {
          return true;
        }
      }
      return false;
    }
  
   private:
    friend class Graph;
    //HW1_CJ
    Graph* graph_;   
    size_type nit_index_;
    //NodeIterator constructor
    NodeIterator(const Graph* graph, size_type nit_index) 
      :  graph_(const_cast<Graph*>(graph)), nit_index_(nit_index) {}
  };

  //HW1_CJ
  /** @brief Return first iterator for Node element
  *   @pre 0 < @a i <= num_nodes() 
  *   @post graph.node_begin() == 0
  */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /** @brief Return last iterator for Node element
  *   @pre 0 < @a i <= num_nodes() 
  *   @post graph.node_end() == num_nodes()
  */
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    //HW2_CJ - HW1 update
    /** @brief Dereference current IncidentIterator to return an Edge
    *   @pre 0 <= node1_index_ < num_nodes() 
    *   @pre 0 <= iit_index_ < adj_[node1_index_].size()
    *   @post  0<= (graph_->adj_[node1_index_][iit_index_]).first < num_edges()
    */
    Edge operator*() const {
      //return Edge corresponding to current iterator position in adj_ sub-vector
      return Edge(graph_, node1_index_, (graph_->adj_[node1_index_][iit_index_]).first);
    }

    /** @brief Increment IncidentIterator to return next reference
    *   @pre 0 <= iit_index_ < num_edges
    *   @post  0 <= *this < num_edges()
    */
    IncidentIterator& operator++() {
      ++iit_index_;
      return *this;
    }

    /** @brief Evaluate if current IncidentIterator is equal to another
    *   @pre 0 <= iit_index_ < num_nodes()
    *   @pre 0 <= iiter.iit_index_ < num_nodes()
    *   @post iiter.graph_ == this->graph_
    */
    bool operator==(const IncidentIterator& iiter) const {
      //check graph
      if (this->graph_ == iiter.graph_) {
        //check indices
        if (this->iit_index_ == iiter.iit_index_ && this->node1_index_ == iiter.node1_index_) {
          return true;
        }
      }
      return false;
    }

   private:
    friend class Graph;
    //HW1_CJ
    Graph* graph_;
    size_type node1_index_;
    size_type iit_index_;
    //IncidentIterator constructor
    IncidentIterator(const Graph* graph, size_type node1_index, size_type iit_index)
      :  graph_(const_cast<Graph*>(graph)), node1_index_(node1_index), iit_index_(iit_index) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    //HW2_CJ - HW1 update
    /** @brief Dereference current EdgeIterator to return an Edge
    *   @pre 0 <= node1.index() < num_nodes()
    *   @pre 0 <= node2.index() < num_nodes()
    *   @post Edge(graph_,node1.index(),node2.index()).graph_ == this->graph_
    */
    Edge operator*() const {
      //return Edge corresponding to nodes for current iterator position 
      Node node1 = Node(graph_, node1_uid_);
      Node node2 = Node(graph_, (graph_->adj_[node1_uid_][adj_node_idx_]).first);
      return Edge(graph_, node1.index(), node2.index());
    }

    //HW2_CJ - HW1 update
    /** @brief Increment current EdgeIterator to return next EdgeIterator
    *   @pre 0 <= adj_node_idx_ < graph_->node(node1_uid_).degree()
    *   @pre 0 <= node1_uid_ < num_nodes()
    *   @post *this != nullptr
    */
    EdgeIterator& operator++() {
      //increment inner vector counter
      ++adj_node_idx_;
      //iterate through all outer node vectors that make adj_
      size_type outer_vec_size = graph_->adj_.size();
      while (node1_uid_ < outer_vec_size) {
        //iterate through all <adj_node_idx,edge_val> pairs within each subvector
        size_type inner_vec_size = graph_->node(node1_uid_).degree();
        while (adj_node_idx_ < inner_vec_size) {
          auto adj_node_uid = (graph_->adj_[node1_uid_][adj_node_idx_]).first;
          //assert global ordering to avoid double counting
          if (node1_uid_ < adj_node_uid) {
            return *this;
          }
          //progress through inner sub-vector until end is reached
          ++adj_node_idx_;
        }
        //thread through each node1_uid_ vector in adj_
        adj_node_idx_ = 0;
        ++node1_uid_; 
      }    
      return *this;
    }

    /** @brief Evaluate if current EdgeIterator is equal to another
    *   @pre this->graph_ == eiter.graph_
    *   @pre 0 <= this->node1_uid_ < num_nodes()
    *   @pre 0 <= adj_node_idx_ < adj_[node1_uid_].size()
    *   @post eiter.graph_ == this->graph_
    */
    bool operator==(const EdgeIterator& eiter) const {
      //edges must be in same graph
      if (this->graph_ == eiter.graph_) {
        //outer node1 IDs must be the same
        if (this->node1_uid_ == eiter.node1_uid_) {
          //adj_node_idx must be same
          if (this->adj_node_idx_ == eiter.adj_node_idx_) {
            return true;
          }
        }
      }
      return false;
    }

   private:
    friend class Graph;
    //HW1_CJ
    Graph* graph_;   
    //HW2_CJ - HW1 update
    size_type node1_uid_;
    size_type adj_node_idx_;
    //EdgeIterator Constructor
    EdgeIterator(const Graph* graph, size_type node1_uid, size_type adj_node_idx)
      :  graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), adj_node_idx_(adj_node_idx) {}
  };

  //HW2_CJ - HW1 update
  /** @brief Return first iterator for Edge element 
  *   @pre 0 <= @a i < num_edges()
  *   @post node.edge_begin() == 0
  */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0,0);
  }
  
  /** @brief Return last iterator for Edge element 
  *   @pre 0 <= @a i < num_edges()
  *   @pre 0 <= adj_.size() <= num_nodes()
  *   @post graph.edge_end() == adj_.size()
  */
  edge_iterator edge_end() const {
    return edge_iterator(this, this->adj_.size(), 0);
  }

  //
  // REMOVE METHODS
  //
  
  //HW2_CJ
  
  /** @brief Method to remove node from graph 
  *   @param[in] node    Node to be erased
  *   @result node is removed from graph; indices are updated accordingly
  *
  *   @pre 0 <= node.index() < num_nodes()
  *   @pre 0 <= node.index() < adj_.size()
  *   @pre node.graph_ == this->graph_
  *   @post Node uids for all nodes with uid greater than node are invalid
  *   @post 0 <= new adj_.size() < old adj_.size()
  *   @post Iterator range [it, old end()] are invalidated
  *   @post New iterator range [next_it, new end()) is established
  *   
  *   Complexity: O(num_nodes())
  */
 
  size_type remove_node(const Node& node) {
    //verify node is in graph
    if (has_node(node) == false) {
      return 0;
    } 

    //remove all edges with connections to node
    for (auto it = idx2uid_.begin(); it != idx2uid_.end(); ++it) {
      Node cur_node = Node(this, (*it));
      if (has_edge(node,cur_node)) {
        remove_edge(node, cur_node); 
      }
    }

    //remove Node from nodes_ vector
    nodes_.erase(nodes_.begin() + node.index());
    //remove node's subvector from adj_ list using STL iterator
    adj_.erase(adj_.begin() + node.index());
    //iterate through adjacent nodes to check for uids larger than deleted node
    for (auto out_it = adj_.begin(); out_it != adj_.end(); ++out_it) {
      //iterate through subvector to check individual node ids
      for (auto in_it = (*out_it).begin(); in_it != (*out_it).end(); ++in_it) {  
        if ( (*in_it).first > node.index()) {
          //decrememnt node id to preserve continuous order
          --(*in_it).first;
        }
      }
    } 

    //node deleted
    return 1;
  }

  /** @brief Method to remove node from graph 
  *   @param[in,out] n_it    Node iterator used to remove corresponding Node
  *   @result n_it  iterator at formerly next element in relation to @a n_it[in]
  *   
  *   @pre (*n_it).graph_ == this->graph_
  *   @pre 0 <= @a n_it <= old end()
  *   @post All iterators in range [old n_it, old end()) invalidated
  *   @post New iterator range [new n_it, new end()) established
  *   @post 0 <= @a result <= new end() 
  *  
  *   Complexity: O(num_nodes())
  */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** @brief Method to remove Edge from graph 
  *   @param[in] node1   First node in Edge pair
  *   @param[in] node2   Second node in Edge pair
  *   @result Edge between node1 and node2 is erased
  *
  *   @pre node1.graph_ == node2.graph_
  *   @pre node1.degree() > 0 && node2.degree() > 0
  *   @pre 0 <= node1.index() < num_nodes() && 0<= node2.index() < num_nodes()
  *   @pre node1.index() != node2.index()
  *   @post 0 <= new node1.degree() < old node1.degree()
  *   @post 0 <= new node2.degree() < old node2.degree()
  *   @post Iterators, Nodes, Edges containing uid greater than node1 and node2 uids are invalid
  * 
  *   Complexity: O(num_nodes() + num_edges())
  *
  *    NOTE: Attempted (multiple ways) to perform this using STL iterator over sub-vector
  *          with std::iter_swap & push_back() erase method (no vector size updating?)
  */
  size_type remove_edge(const Node& node1, const Node& node2) {
    //verify edge is in graph
    if (has_edge(node1,node2) == false) {
      return 0;
    }

    /*** adj_<node1> subvector ***/
    //iterate through node1's subvector to remove connection and update subvector size
    for(unsigned int subvec_idx = 0; subvec_idx < node1.degree(); ++subvec_idx) {
      //access id in <pair> to find connection between node1 and node2
      if (node2.index() == adj_[node1.index()][subvec_idx].first) {
        //implement efficient erase method in constant time
        auto last_elem = adj_[node1.index()][node1.degree()-1];
        std::swap(adj_[node1.index()][subvec_idx], last_elem);
        adj_[node1.index()].pop_back();
      }
    }

    /*** adj_<node2> subvector ***/
    //iterate through node2's subvector to remove connection and update subvector size
    for(unsigned int subvec_idx = 0; subvec_idx < node2.degree(); ++subvec_idx) {
      //access id in <pair> to find connection between node1 and node2
      if (node1.index() == adj_[node2.index()][subvec_idx].first) {
        //implement efficient erase method in constant time
        auto last_elem = adj_[node2.index()][node2.degree()-1];
        std::swap(adj_[node2.index()][subvec_idx], last_elem);
        adj_[node2.index()].pop_back();
      }
    }

    //edge deleted
    return 1 ;
  }

  /** @brief Method to remove edge
  *   @param[in] edge   Edge to be removed
  *   @result edge is removed from graph
  *   
  *   @pre edge.node1() != edge.node2()
  *   @pre node1.graph_ == node2.graph_
  *   @pre edge.graph_ == this->graph_
  *   @post 0 <= new node1.degree() < old node1.degree()
  *   @post 0 <= new node2.degree() < old node2.degree()
  *   @post Iterators, Nodes, Edges containing uid greater than node1 and node2 uids are invalid
  *  
  *   Complexity: O(num_nodes() + num_edges())
  */

  size_type remove_edge(const Edge& edge) {
    Node node1 = edge.node1();
    Node node2 = edge.node2();
    remove_edge(node1,node2);
    return 1;
  }

  /** @brief Method to remove edge
  *   @param[in,out] e_it   Edge iterator used to remove corresponding Edge
  *   @result e_it  iterator at formerly next element in relation to @a e_it[in]
  *   
  *   @pre (*e_it).graph_ == this->graph_
  *   @pre 0 <= @a e_it <= old end()
  *   @post All iterators in range [old e_it, old end()) invalidated
  *   @post New iterator range [new e_it, new end()) established
  *   @post 0 <= @a result <= new end() 
  *  
  *   Complexity: O(num_nodes() + num_edges())
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge edge = *e_it;
    Node node1 = edge.node1();
    Node node2 = edge.node2();
    remove_edge(node1,node2);
    return e_it;
  }


 private:
  /*** Graph class internals for Node and Edge proxys ***/
  //HW1_CJ
  struct internal_node {
    Point position; 
    node_value_type node_val;
    //HW2_CJ
    size_type index;
  };

//HW2_CJ - HW1 update
/*
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    edge_value_type edge_val;
  };
*/
  /*** Graph class vectors of internals ***/
  std::vector<internal_node> nodes_;
  //std::vector<internal_edge> edges_;

  //HW2_CJ - HW1 update 
  /*** Adjacency list for looping efficiency **/
  //outer vector index: Node ID
  //inner vector: pairs of <connecting node ID, corresponding edge value>
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> adj_; 
  //HW2_CJ
  //map between indices and unique node identifiers
  std::vector<size_type> idx2uid_;

  /*** Graph class data members ***/ 
  //size_type next_node_uid_;
  //size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP
