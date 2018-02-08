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
//HW1_CJ
template <typename V>
class Graph {
 private:

  /*** Predeclare internal structs for node and edge ***/
  struct internal_node;
  struct internal_edge;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), edges_(), adj_(), next_node_uid_(0), next_edge_uid_(0) {
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
      /*** Verify Nodes are in same Graph ***/
      if (n.graph_ != this->graph_) {
        return false;
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
    //instantiate an internal node to add to nodes_ vector
    internal_node new_node;
    //assign data attributes to new_node
    new_node.position = position;
    //HW1_CJ
    new_node.node_val = new_node_value;
    //append new_node to nodes_ vector
    nodes_.push_back(new_node);
    //HW1_CJ
    /*** update adj_ list ***/
    adj_.emplace_back(std::vector<std::pair<size_type, size_type>>());
    //update number of nodes and next node ID 
    ++next_node_uid_;
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
    if ((this == n.graph_) && (num_nodes() < n.index())) {
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
    /*** return a proxy object for node @ i ***/
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
      size_type node1_index = graph_->edges_[edge_uid_].node1_uid;
      return graph_->node(node1_index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      /*** access this edge's node2 member through graph_ ptr ***/
      size_type node2_index = graph_->edges_[edge_uid_].node2_uid;
      return graph_->node(node2_index);
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
      /*** Verify Edges are in same Graph ***/
      if (e.graph_ != this->graph_) {
        return false;
      }
      /*** Compare edge uids for ordering purposes ***/
      if (e.edge_uid_ > this->edge_uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    /*** Private data members and methods for Edge ***/
    // Pointer back to the Graph container
    Graph* graph_;
    //This Edge's unique ID number
    size_type edge_uid_;
    /* Unique IDs of Nodes connected by this Edge */
    size_type node1_uid_;
    size_type node2_uid_;
    /* private constructor Graph::Edge */
    Edge(const Graph* graph, size_type edge_uid)
      : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    /*** number of edges present in the graph ***/
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /*** return a proxy object for edge @ i ***/
    assert(i < num_edges());
    return Edge(this, i);
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
      //HW1_CJ
      /*** Iterate through adj_ list to see if an edge contains a and b as nodes ***/
      std::vector<std::pair<size_type, size_type>> adj_elem = adj_[a.node_uid_];
      unsigned int adj_elem_size = adj_elem.size();
      for (size_type i=0; i < adj_elem_size; ++i) {
        //Check if two nodes are connect by Edge via adj_ list
        if (b.index() == adj_elem[i].first) {
          return true;
        } 
      }
    }
    //return false if a and b are not within an Edge in edges_
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
  /*** Add new internal_edge to edges_ vector of internal_edges ***/
  Edge add_edge(const Node& a, const Node& b) {
    //test if edge already exists via has_edge()
    if (has_edge(a,b)) {
      /*** Iterate through adj_ list to see if an edge contains a and b as nodes ***/
      std::vector<std::pair<size_type, size_type>> adj_elem = adj_[a.node_uid_];
      unsigned int adj_elem_size = adj_elem.size();
      for (size_type i = 0; i < adj_elem_size; ++i) {
        //Check if Node b's ID is in adj_elem; if so, return the corresponding Edge ID 
        if (b.index() == adj_elem[i].first) {
          return edge(adj_elem[i].second); 
        }
      }
    }
    /*** Add edge if not already existent ***/
    //instantiate an internal edge to add to edges_ vector
    internal_edge new_edge;
    //assign data attributes to new_edge
    new_edge.node1_uid = a.index();
    new_edge.node2_uid = b.index();
    //append new_edge to edges_ vector
    edges_.push_back(new_edge);
    //HW1_CJ
    /*** Update adj_ list corresponding to both nodes ***/
    //add pair of (connecting nodes, corresponding edge ID)
    adj_[a.node_uid_].push_back(std::make_pair(b.node_uid_, next_edge_uid_));
    adj_[b.node_uid_].push_back(std::make_pair(a.node_uid_, next_edge_uid_));
    //update number of edges and next edge ID
    ++next_edge_uid_;

    return Edge(this, this->edges_.size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    /*** clear nodes_ and edges_ vectors ***/ 
    nodes_.clear();
    edges_.clear();
    adj_.clear();
    //reset ID numbers
    next_node_uid_ = 0;
    next_edge_uid_ = 0;
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

    //HW1_CJ
    /** @brief Dereference current IncidentIterator to return an Edge
    *   @pre 0 <= node1_index_ < num_nodes() 
    *   @pre 0 <= iit_index_ < num_edges()
    *   @post  0<= pair_idx.second < num_edges()
    */
    Edge operator*() const {
      //return Edge corresponding to current iterator position in adj_
      std::pair<size_type,size_type> pair_idx = graph_->adj_[node1_index_][iit_index_];
      //edge ID is the second element in pair
      return graph_->edge(pair_idx.second);
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
    *   @pre 0 <= iit_index_ < num_edges()
    *   @pre 0 <= iiter.iit_index_ < num_edges()
    *   @post iiter.graph_ == this->graph_
    */
    bool operator==(const IncidentIterator& iiter) const {
      if (this->graph_ == iiter.graph_) {
        if (this->iit_index_ == iiter.iit_index_) {
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

    //HW1_CJ
    /** @brief Dereference current EdgeIterator to return an Edge
    *   @pre 0 <= eit_index < num_edges()
    *   @post edge(eit_index_).graph_ == this->graph_
    */
    Edge operator*() const {
      //return Edge corresponding to current iterator position in egdes_ vector
      return graph_->edge(eit_index_);
    }

    /** @brief Increment current EdgeIterator to return next EdgeIterator
    *   @pre 0 <= eit_index < num_edges()
    *   @pre 0 <= ++eit_index < num_edges()
    *   @post ++eit_index != nullptr
    */
    EdgeIterator& operator++() {
      ++eit_index_;
      return *this;
    }

    /** @brief Evaluate if current EdgeIterator is equal to another
    *   @pre this->graph_ == eiter.graph_
    *   @pre 0 <= this->eit_index_ < num_edges()
    *   @pre 0 <= eiter.eit_index_ < num_edges()
    *   @post eiter.graph_ == this->graph_
    */
    bool operator==(const EdgeIterator& eiter) const {
      //edges must be in same graph
      if (this->graph_ == eiter.graph_) {
        //edge IDs must be the same
        if (this->eit_index_ == eiter.eit_index_) {
          return true;
        }
      }
      return false;
    }

   private:
    friend class Graph;
    //HW1_CJ
    Graph* graph_;   
    size_type eit_index_;
    //EdgeIterator constructor
    EdgeIterator(const Graph* graph, size_type eit_index)
      :  graph_(const_cast<Graph*>(graph)), eit_index_(eit_index) {}
  };

  //HW1_CJ
  /** @brief Return first iterator for Edge element 
  *   @pre 0 <= @a i <= num_edges()
  *   @post node.edge_begin() == 0
  */
  edge_iterator edge_begin() const {
      return edge_iterator(this, 0);
  }
  
  /** @brief Return last iterator for Edge element 
  *   @pre 0 <= @i <= num_edges()
  *   @post graph.edge_end() == num_edges()
  */
  edge_iterator edge_end() const {
      return edge_iterator(this, num_edges());
  }

 private:
  /*** Graph class internals for Node and Edge proxys ***/
  //HW1_CJ
  struct internal_node {
    Point position; 
    //HW1_CJ
    node_value_type node_val;
  };

  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
  };

  /*** Graph class vectors of internals ***/
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  //HW1_CJ
  /*** Adjacency list for looping efficiency **/
  //outer vector index: Node ID
  //inner vector: pairs of (connecting node ID, corresponding edge id)
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_; 

  /*** Graph class data members ***/ 
  size_type next_node_uid_;
  size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP
