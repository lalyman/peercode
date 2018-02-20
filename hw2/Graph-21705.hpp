#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <iterator>
#include <vector>
#include <tuple>
#include <cassert>

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

  // predeclare internal structs that will hold data for nodes and edges
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

  /** Type of node value, a user-specified value */
  using node_value_type = V;
 
  /** Type of edge value, a user-specified value */
  using edge_value_type = E; 

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    :  nodes_(), adj_() {
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
    }

    /** Return this node's position, allow modification. */
    Point& position() const{
      return graph_->nodes_[uid_].position; //-> to access private member
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /* Return this node's value */
    node_value_type& value() {
      return graph_->nodes_[uid_].val;
    }

    /* Return this node's value */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].val;
    }

    /* Returns number of edges incident to this Node */
    size_type degree() const {
      return graph_->adj_[uid_].size();
    }

    /* Returns the first edge iterator for an IncidentIterator */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }
    /* Returns the last edge iterator for an IncidentIterator */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check whether uid's are the same
      // Node can access another node's graph_ because Graph is a friend class
      return (this->uid_ == n.index() and this->graph_ == n.graph_);
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
      // check whether nodes are in the same graph	
      // compare uid's, which are ordered by construction
      if (this->graph_ == n.graph_) {
        return (this->uid_ < n.index());
      }
      else {
        return (this->graph_ < n.graph_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_; //pointer to graph object
    size_type uid_; // This element's unique identification number

    /**Private Constructor accessed by Graph to construct valid Node objects*/
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //returns size of nodes_ vector
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {

    //create a new internal_node - instantiating a struct
    internal_node next_node;

    //set up attributes of new node
    next_node.position = position;
    next_node.val = value;

    nodes_.push_back(next_node);

    //push a placeholder vector in adjacency matrix
    adj_.emplace_back(std::vector<std::tuple<size_type, edge_value_type>>());

    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //return true if this is the same as graph_ (both are memory addresses)
    return (n.graph_ == this and n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // check that index is within size of graph (i.e. it already exists)
    assert(i < size());
  
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      //return node at location of first node
      return graph_->node(n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return node at this location of second node
      return graph_->node(n2_uid_);
    }

    /** Return length of edge, i.e. Euclidian distance between two nodes */
    double length() const {
      double delx = node1().position().x - node2().position().x;
      double dely = node1().position().y - node2().position().y;
      double delz = node1().position().z - node2().position().z;
      return sqrt(pow(delx, 2.0) + pow(dely, 2.0) + pow(delz, 2.0));
    }


    /** Set this edge's value to @a val 
     * @param[in] val The object to be set to this edge's value
     *
     * Complexity: max O(d), where d is largest degree of a node
     */
    void set_value(edge_value_type val) {
      //NOTE: not straightforward to use IncidentIterator, since we need to keep
      // track of how many iterations forward we've moved 
      //set val of first instance in adj_
      for (size_type i = 0; i < graph_->adj_[n1_uid_].size(); i++) {
        // check first element of tuple (i.e. uid for node)
        if (std::get<0>(graph_->adj_[n1_uid_][i]) == n2_uid_) { 
          std::get<1>(graph_->adj_[n1_uid_][i]) = val;
        }
      }
      //set val of second instance in adj_
      for (size_type i = 0; i < graph_->adj_[n2_uid_].size(); i++) {
        // check first element of tuple (i.e. uid for node)
        if (std::get<0>(graph_->adj_[n2_uid_][i]) == n1_uid_) { 
          std::get<1>(graph_->adj_[n2_uid_][i]) = val;
        }
      }
    }

    /** Return this edge's value 
     * Complexity: max O(d), where d is largest degree of a node
     */
    edge_value_type& value() {
      //NOTE: this method is unnecessary with set_value implemented
      for (size_type i = 0; i < graph_->adj_[n1_uid_].size(); i++) {
        // check first element of tuple (i.e. uid for node)
        if (std::get<0>(graph_->adj_[n1_uid_][i]) == n2_uid_) { 
          return std::get<1>(graph_->adj_[n1_uid_][i]);
        }
      }
      return std::get<1>(graph_->adj_[0][0]); //this will never be reached, just to silence warning
    }

    /** Return a const of this edge's value
     * Complexity: max O(d), where d is largest degree of a node
     */
    const edge_value_type& value() const {

      for (size_type i = 0; i < graph_->adj_[n1_uid_].size(); i++) {
        // check first element of tuple (i.e. uid for node)
        if (std::get<0>(graph_->adj_[n1_uid_][i]) == n2_uid_) { 
          return std::get<1>(graph_->adj_[n1_uid_][i]);
        }
      }
      return std::get<1>(graph_->adj_[0][0]); //this will never be reached, just to silence warning
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // check whether uid and graph are the same
      // Edge can access another Edge's member attributes b/c Graph is friend
      if (graph_ == e.graph_) {
        return ((e.node1() == node1() and e.node2() == node2())
                or (e.node1() == node2() and e.node2() == node1()));
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // check that the edges are in same graph, and compare uid's
      if (graph_ == e.graph_) {
        if (node1() == e.node1() and node2() == e.node2()) {
          return false;
        }
        return (node1() < e.node1() or (node1() == e.node1() and node2() < e.node2()));
      }
      else {
        return (this->graph_ < e.graph_);
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_; // Pointer to Graph container
    size_type n1_uid_;// Unique identification number of one node
    size_type n2_uid_;// Unique identification number of the other node
    /** Private constructor accessed by Graph to construct valid Edge objects */
    Edge(const Graph* graph, size_type n1_uid, size_type n2_uid)
      : graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_uid_(n2_uid)  {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    //use stl iterator over adj_
    size_type sum = 0;
    for (auto it = adj_.begin(); it != adj_.end(); ++it) {
      sum += (*it).size();
    }
    return sum/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    //use std method to find ith step in edge iterator
    auto ith_iter = std::next(edge_begin(), i);
    return *ith_iter;
  }

  /** Test whether two nodes are connected by an edge.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check that a and b are in the graph
    if (this != a.graph_ or this != b.graph_) {
      return false;
    }
    //check internal vector of adjacency matrix to see if b is in a's vector
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
      // check first element of tuple (i.e. uid for node)
      Edge e = *it;
      if (e.node1() == b or e.node2() == b) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {

    //check whether edge already exists
    if (has_edge(a, b)) {
      return Edge(this, a.index(), b.index());
    }

    //add entries to adjacency matrix
    adj_[a.index()].push_back(std::make_tuple(b.index(), value));
    adj_[b.index()].push_back(std::make_tuple(a.index(), value));

    //return a new edge object
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // use vector clear to remove all data from nodes_ and adj_
    nodes_.clear();
    adj_.clear();
  }



  /** Remove edge addressed by iterator @a e_it
   * @param e_it an iterator pointing to Edge to be erased
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove Edge @a e
   * @param e the Edge to be removed
   */
  size_type remove_edge(const Edge& e) {
    remove_edge(e.node1(), e.node2());
  }


  /** Remove edge connecting Node @a a and Node @a b
   * @param a a Node to which the edge belongs
   * @param a the other Node to which the edge belongs
   */
  size_type remove_edge(const Node& a, const Node& b) {
    //if graph does not contain edge, don't remove!
    if (!(has_edge(a, b))) {
      return 0;
    }

    size_type a_uid = a.index();
    size_type b_uid = b.index();

    //TODO! if removing an edge leaves one (or both) of its nodes
    // unconnected, we need to remove that node (or both nodes)!
    // Node is unconnected if its inner vector is empty
    // Commented out code below attempted to do this..., but kept
    // failing tests 

    //remove entry in a's inner vector in adj_ 
    for	(size_type i = 0; i < adj_[a_uid].size(); i++) {
      if (std::get<0>(adj_[a_uid][i]) == b_uid) {

        //if (false) {
        /*
        //remove node a if this edge is it's only connection
        if (adj_[a_uid].size() == 1) {
          std::cout << "Empty a -------------------" << std::endl;

          //remove entries from b's vector

          for (size_type i = 0; i < adj_[b_uid].size(); i++) {
            if (std::get<0>(adj_[b_uid][i]) == a_uid) {
              adj_[b_uid].erase(adj_[b_uid].begin() + i);
            }
          }        
  
          // 2) Remove uid of Node a from nodes_ vector
          nodes_.erase(nodes_.begin() + a_uid);

          // 3) Remove Node a's inner vector from adj_
          adj_.erase(adj_.begin() + a_uid);
  
          // 4) Correct ordering in adj_
          // Loop through adj_, if a uid in inner vector is greater than a_uid, decrement
          for (auto it = adj_.begin(); it != adj_.end(); it++) {
            //(*it) is inner vector
            for (size_type i = 0; i < (*it).size(); i++) {
              //sanity check
              if (std::get<0>((*it)[i]) == a_uid) {
                std::cout << "Hi!" << std::endl;
                std::cout << "This should have been removed!!" << std::endl;
                std::cout << "a_uid: " << a_uid << std::endl;
              }
              if (std::get<0>((*it)[i]) > a_uid) {
                std::get<0>((*it)[i])--; //decrement
              }
            }
          }

          // need now check whether b has any connections
          if (b.degree() == 0) {
   	    std::cout << "b has no connections" << std::endl;
            remove_node(b);
          }

          return 1; //exit, we've already removed b's entries
        
        */
        //}

        //else just erase element in inner vector
        //else {
          adj_[a_uid].erase(adj_[a_uid].begin() + i);
        //} 

      }
    }


    //remove entry in b's inner vector in adj_
    for	(size_type i = 0; i < adj_[b_uid].size(); i++) {
      if (std::get<0>(adj_[b_uid][i]) == a_uid) {

        //if (false) {
        /*
        //if this is it's only connection
        if (adj_[b_uid].size() == 1) {
          std::cout << "Empty b -------------------" << std::endl;

          //remove entries from a's vector
          for (size_type i = 0; i < adj_[a_uid].size(); i++) {
            if (std::get<0>(adj_[a_uid][i]) == b_uid) {
              adj_[b_uid].erase(adj_[b_uid].begin() + i);
            }
          }        

          // 2) Remove uid of Node b from nodes_ vector
          nodes_.erase(nodes_.begin() + b_uid);

          // 3) Remove Node b's inner vector from adj_
          adj_.erase(adj_.begin() + b_uid);
  
          // 4) Correct ordering in adj_
          // Loop through adj_, if b uid in inner vector is greater than b_uid, decrement
          for (auto it = adj_.begin(); it != adj_.end(); it++) {
            //(*it) is inner vector
            for (size_type i = 0; i < (*it).size(); i++) {
              //sanity check
              if (std::get<0>((*it)[i]) == b_uid) {
                std::cout << "This should have been removed!!" << std::endl;
                std::cout << "b_uid: " << a_uid << std::endl;
              }
              if (std::get<0>((*it)[i]) > b_uid) {
                std::get<0>((*it)[i])--; //decrement
              }
            }
          }

          // need now check whether b has any connections
          if (a.degree() == 0) {
   	    std::cout << "a has no connections" << std::endl;
            remove_node(a);
          }

          return 1; //exit, we've already removed b's entries
          */
        //}


        //erase element in inner vector
        //else {
          adj_[b_uid].erase(adj_[b_uid].begin() + i);
        //}   

      }
    }

    return 1; //for test_edges
  }

 

  /** Remove node addressed by iterator @a n_it
   * @param n_it an iterator pointing to Edge to be erased
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

 
  /** Remove Node @a a and all edges it is connected to 
   * @param a the Node to remove
   */
  void remove_node(const Node& a) {
    //if graph does not contain node, don't remove!
    if (!(has_node(a))) {
      return;
    }

    //save index (for safety)
    size_type a_uid = a.index();

    // 1) Remove all connections to Node a
    std::vector<size_type> a_connections; //temp structure to store connected nodes
    // iterate through a's inner vector and store connections
    for	(auto it = adj_[a_uid].begin(); it != adj_[a_uid].end(); it++) {
      // *it is a tuple
      a_connections.push_back(std::get<0>(*it));
    }
    //call remove_edge on each connection
    for (auto it = a_connections.begin(); it != a_connections.end(); it++) {
      remove_edge(a, Node(this,*it));
    }

    // 2) Remove uid of Node a from nodes_ vector
    nodes_.erase(nodes_.begin() + a_uid);

    // 3) Remove Node a's inner vector from adj_
    adj_.erase(adj_.begin() + a_uid);
  
    // 4) Correct ordering in adj_
    // Loop through adj_, if a uid in inner vector is greater than a_uid, decrement
    for (auto it = adj_.begin(); it != adj_.end(); it++) {
      //*it is inner vector
      for (size_type i = 0; i < (*it).size(); i++) {
        //sanity check
        if (std::get<0>((*it)[i]) == a_uid) {
          std::cout << "This should have been removed!!" << std::endl;
          std::cout << "a_uid: " << a_uid << std::endl;
        }
        if (std::get<0>((*it)[i]) > a_uid) {
          std::get<0>((*it)[i])--; //decrement
        }
      }
    }

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

    /** Dereferences the node iterator to return a Node object*/
    Node operator*() const {
      return graph_->node(index_);
    }
    /** Increments the node iterator*/
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }
    /** Checks whether this NodeIterator equals @a x
     * @param[in] x The NodeIterator to compare to this 
     */ 
    bool operator==(const NodeIterator& x) const {
      return (this->graph_ == x.graph_ and this->index_ == x.index_);
    }

   private:
    friend class Graph;
    Graph* graph_; //pointer to graph object
    size_type index_;// counter for node iterator class
    /** Private constructor for node_iterator class */
    NodeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  /** Returns the first node iterator  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Returns the last node iterator */
  node_iterator node_end() const {
    return NodeIterator(this, this->size());
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

    /** Dereferences the incident iterator to return an Edge object */ 
    Edge operator*() const {
      size_type n2_uid = std::get<0>(graph_->adj_[node1_uid_][index_]);
      return Edge(graph_, node1_uid_, n2_uid);
    }

    /** Increments the incident iterator */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Check whether this IncidentIterator equals @a x
     * @param x The IncidentItertor to compare to this
     */
    bool operator==(const IncidentIterator& x) const {
      return (this->graph_ == x.graph_ and this->index_ == x.index_ and this->node1_uid_ == x.node1_uid_);
    }

   private:
    friend class Graph;
    Graph* graph_; //pointer to graph object
    size_type node1_uid_;// uid of node that spawned the incident iterator
    size_type index_;// counter for incident iterator class
    /* Private constructor for incident_iterator class */
    IncidentIterator(const Graph* graph, size_type node1_uid, size_type index)
      : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), index_(index) {
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

    /** Dereferences the edge iterator to return an Edge object */
    Edge operator*() const {
      Node n1 = Node(graph_, uid1_);
      Node n2 = Node(graph_, std::get<0>(graph_->adj_[uid1_][neighbor_index_])); 
      return Edge(graph_, n1.index(), n2.index());
    }

    /** Increments the edge iterator */
    EdgeIterator& operator++() {
      // Moves down inner vector, if it gets to the end moves the next index in the outer
      // Returns an iterator only if nodeA < nodeB to avoid double counting edges
      ++neighbor_index_;
      while (uid1_ < graph_->adj_.size()) {
    
	if (neighbor_index_ == graph_->adj_[uid1_].size()) {
	  ++uid1_;
	  neighbor_index_ = 0;
          //return end iterator
          if (uid1_ == graph_->num_nodes()) {
            return *this;
          }
        }

        size_type uid2 = std::get<0>(graph_->adj_[uid1_][neighbor_index_]);
        
        if (graph_->node(uid1_) < graph_->node(uid2)) { //use bool operator of Node
          return *this;
        }
        ++neighbor_index_;
      }
      return *this;
    }

    /** Checks whether this EdgeIterator is the same as x 
     * @param x The EdgeIterator to compare to this
     */
    bool operator==(const EdgeIterator& x) const {
      return (this->graph_ == x.graph_ and this->uid1_ == x.uid1_ 
              and this->neighbor_index_ == x.neighbor_index_);
    }

   private:
    friend class Graph;
    Graph* graph_; //pointer to graph object
    size_type uid1_;//counter for outer vector of adj_
    size_type neighbor_index_;//counter for inner vector of adj_
    /* Private constructor for EdgeIterator class */
    EdgeIterator(const Graph* graph, size_type uid1, size_type neighbor_index)
      : graph_(const_cast<Graph*>(graph)), uid1_(uid1), neighbor_index_(neighbor_index) {
    }
  };

  /** Returns the first edge iterator */
  edge_iterator edge_begin() const {
    //return EdgeIterator(this, 0);
    return EdgeIterator(this, 0, 0);
  }
  /** Returns the last edge iterator */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }

 private:

  /** Structure containing data attributes that define a node */
  struct internal_node {
    Point position;
    node_value_type val;
  };

  // STL containers for internal_nodes and adj_, which contains the
  // connections of each node
  std::vector<internal_node> nodes_;
  std::vector<std::vector<std::tuple<size_type, edge_value_type>>> adj_;

};

#endif // CME212_GRAPH_HPP
