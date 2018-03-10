#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

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
  using node_value_type = V;
  using edge_value_type = E;

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
  };

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
  class Node : private totally_ordered<Node>{
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
    Node()
      : graph_(nullptr){
    }

    /** Return this node's position. */
    const Point& position() const {
      size_type uid = graph_->i2u_[idx_];
      return graph_->nodes_[uid].location;
    }

    Point& position() {
      size_type uid = graph_->i2u_[idx_];
      return graph_->nodes_[uid].location;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->idx_;
    }

    /** Return this node's value, of templated type V */
    node_value_type& value(){
      size_type uid = graph_->i2u_[idx_];
      return graph_->nodes_[uid].value;
    }

    /** Return this node's value, of templated type V */
    const node_value_type& value() const{
      size_type uid = graph_->i2u_[idx_];
      return graph_->nodes_[uid].value;
    }

    /** Return the number of nodes connected by an edge to this node */
    size_type degree() const{
      size_type uid = graph_->i2u_[idx_];
      return graph_->nodes_[uid].adj_nodes.size();
    };

    /** Create a beginning iterator for the edges incident to this node*/
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, idx_, 0);
    };

    /** Create an ending iterator for the edges incident to this node*/
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, idx_, this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        if(graph_ == n.graph_)
          if(idx_ == n.idx_)
            return true;
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
      if(this->graph_ < n.graph_)
        return true;
      if(idx_ < n.idx_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    //Pointer back to the graph container
    Graph* graph_;
    size_type idx_;

    Node(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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

    size_type next_uid = nodes_.size();

    internal_node_element node_to_add = {
      position,
      value,
      std::vector<internal_edge_element>(),
      next_uid //a node's uid, safe since we don't delete from nodes_
    };

    nodes_.push_back(node_to_add);
    i2u_.push_back(next_uid);

    return node(nodes_.size()-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(this != n.graph_)
      return false;
    if(n == node(n.idx_))
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    //i refers to an INDEX not a uid
    assert(i < nodes_.size());
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
    Edge()
      : graph_(nullptr){
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(n1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2_idx_);
    }

    edge_value_type& value(){
      size_type n1_uid = graph_->i2u_[n1_idx_];
      size_type n2_uid = graph_->i2u_[n2_idx_];
      auto adj_list = graph_->nodes_[n1_uid].adj_nodes;
      size_type i = 0;
      for(i = 0; i < adj_list.size(); i++){
        if(n2_uid == adj_list[i].node_id){
          break;
        }
      }
      return graph_->nodes_[n1_uid].adj_nodes[i].value;
    }

    const edge_value_type& value() const{
      size_type n1_uid = graph_->i2u_[n1_idx_];
      size_type n2_uid = graph_->i2u_[n2_idx_];
      auto adj_list = graph_->nodes_[n1_uid].adj_nodes;
      size_type i = 0;
      for(i = 0; i < adj_list.size(); i++){
        if(n2_uid == adj_list[i].node_id){
          break;
        }
      }
      return graph_->nodes_[n1_uid].adj_nodes[i].value;
    }

    double length(){
      Point p1 = graph_->node(n1_idx_).position();
      Point p2 = graph_->node(n2_idx_).position();
      return norm_2(p1 - p2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */

    bool operator==(const Edge& e) const {

      //if both edges empty, true
      if(this->graph_ == nullptr && e.graph_ == nullptr)
        return true;

      //if either edge empty (but not both) throw false
      if(this->graph_ == nullptr && e.graph_ != nullptr)
        return false;
      if(this->graph_ != nullptr && e.graph_ == nullptr)
        return false;

      //check global before index
      if(this->graph_ != e.graph_)
        return false;

      //else go ahead and access the edges, checking node equality
      if(this->node1().index() == e.node1().index()
      && this->node2().index() == e.node2().index())
        return true;
      if(this->node1().index() == e.node2().index()
      && this->node2().index() == e.node1().index())
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(this->graph_ < e.graph_)
        return true;
      if(this->node1().index() < e.node1().index())
        return true;
      if(this->node1().index() == e.node1().index()
      && this->node2().index() <  e.node2().index())
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    //Pointer back to the graph container
    Graph* graph_;
    size_type n1_idx_;
    size_type n2_idx_;

    Edge(const Graph* graph, size_type n1_idx_, size_type n2_idx_)
      : graph_(const_cast<Graph*>(graph)),
        n1_idx_(n1_idx_),
        n2_idx_(n2_idx_){
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type edge_counter = 0;
    for(auto ei = this->edge_begin(); !(ei == this->edge_end()); ++ei){
      edge_counter++;
    }
    return edge_counter;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  Edge edge(size_type i) const {
    assert(i < this->num_edges());
    auto ei = std::next(this->edge_begin(), i);
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(edge_finder(a,b) == Edge()){
      return false;
    } else {
      return true;
    }
  }

  /** Function that finds an edge already stored in the graph,
   *  and returns a null edge if it doesn't find it
  */
  Edge edge_finder(const Node& a, const Node& b) const {
    size_type n1_uid = i2u_[a.index()];
    size_type n2_uid = i2u_[b.index()];
    auto adj_nodes = nodes_[n1_uid].adj_nodes;

    for(size_type i = 0; i < adj_nodes.size(); i++){
      size_type adj_node_uid = adj_nodes[i].node_id;
      if(adj_node_uid == n2_uid)
        return Edge(this, nodes_[n1_uid].idx_, nodes_[n2_uid].idx_);
    }

    return Edge();
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
    //find the edge if it already exists
    Edge edge_if_exists = edge_finder(a,b);

    //add edge if it doesn't exist already
    if(edge_if_exists == Edge()){
      //add an edge to nodes_ adjecency list, each way...
      nodes_[a.index()].adj_nodes.push_back({b.index(), edge_value_type()});
      nodes_[b.index()].adj_nodes.push_back({a.index(), edge_value_type()});
      return Edge(this, a.index(), b.index());

    //otherwise, return the edge
    }else{
      return edge_if_exists;
    }

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
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    NodeIterator(const Graph* graph, size_type idx)
       : graph_(const_cast<Graph*>(graph)), current_idx(idx) {
    }

    /** De-references a node iterator to the current node*/
    Node operator*() const {
      //This is OK since if the user de-references g.node_end() there
      //should be an error, so don't need to do any checking for that
      return graph_->node(this->current_idx);
    }

    /** Defines equality for node_iterator */
    bool operator==(const NodeIterator& x) const {
      if(graph_ == x.graph_)
        if(current_idx == x.current_idx)
          return true;
      return false;
    }

    /** Defines incrementing for node_iterator */
    NodeIterator& operator++() {
      current_idx++;
      return *this;
    }

   private:
     friend class Graph;
     Graph* graph_;
     size_type current_idx;
  };

  /** Creates beginning iterator for iterating through a graph's nodes */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Creates ending iterator for iterating through a graph's nodes */
  node_iterator node_end() const {
    return NodeIterator(this, i2u_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** De-references an incident iterator to the current edge  */
    Edge operator*() const{
      size_type n1_uid = graph_->i2u_[n1_idx_];
      auto adj_nodes = graph_->nodes_[n1_uid].adj_nodes;
      size_type n2_uid = adj_nodes[iterator_location_].node_id;
      size_type n2_idx =  graph_->nodes_[n2_uid].idx_;
      return Edge(graph_, n1_idx_, n2_idx);
    }

    /** Defines incrementing for incident_iterator */
    IncidentIterator& operator++(){
      iterator_location_++;
      return *this;
    }

    /** Defines equality for incident_iterator */
    bool operator==(const IncidentIterator& x) const{
      if(graph_ == x.graph_)
        if(n1_idx_ == x.n1_idx_)
          if(iterator_location_ == x.iterator_location_)
            return true;
      return false;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type n1_idx_;
    size_type iterator_location_;

   IncidentIterator(const Graph* graph, size_type n1_idx_, size_type iterator_location_)
     : graph_(const_cast<Graph*>(graph)),
       n1_idx_(n1_idx_),
       iterator_location_(iterator_location_){
   }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator> {
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

    /** De-references an edge_iterator to its current edge */
    Edge operator*() const {
      size_type n1_uid = graph_->i2u_[n1_idx_];

      auto adj_list = graph_->nodes_[n1_uid].adj_nodes;
      size_type n2_uid = adj_list[adj_list_pos_].node_id;
      size_type n2_idx = graph_->nodes_[n2_uid].idx_;

      Edge edge_to_return = graph_->edge_finder(graph_->node(n1_idx_),
                                                graph_->node(n2_idx));
      return edge_to_return;
    }

    /** Defines incrementing for edge_iterator */
    EdgeIterator& operator++() {
      size_type n1_uid = graph_->i2u_[n1_idx_];
      auto adj_list = graph_->nodes_[n1_uid].adj_nodes;

      //go until we get a non-trivial entry of the adjacency list
      //this is critical for the first iteration, and should only enter on first
      //if the first entry of the adjacency list is empty
      while(adj_list.size() == 0){
        if(n1_idx_ == graph_->i2u_.size()-1){
          n1_idx_++;
          return *this;
        }

        n1_idx_++;
        n1_uid = graph_->i2u_[n1_idx_];
        adj_list = graph_->nodes_[n1_uid].adj_nodes;
        adj_list_pos_ = 0;

        if(adj_list.size() != 0)
          return *this;
      }

      //increment our place in the adj list by 1, or go to next node
      if(adj_list_pos_ < adj_list.size() - 1){
        adj_list_pos_++;
      }else if(adj_list_pos_ == adj_list.size() - 1){
        n1_idx_++;
        n1_uid = graph_->i2u_[n1_idx_];
        adj_list = graph_->nodes_[n1_uid].adj_nodes;
        adj_list_pos_ = 0;
      }

      //get information on second node to check if our next step will be valid
      size_type n2_uid = adj_list[adj_list_pos_].node_id;
      size_type n2_idx = graph_->nodes_[n2_uid].idx_;

      //check that this next step is valid, iterate until the end or the next
      //valid step it found
      while(n1_idx_ > n2_idx || adj_list.size() == 0){
        //check if we are at the final edge in the adjacency matrix
        if(n1_idx_ == graph_->i2u_.size()-1 && \
           adj_list_pos_ == adj_list.size()-1){
          n1_idx_++;
          adj_list_pos_++;
          return *this;
        //if there are no adjacent nodes, to n1, go to the next n1
        }else if(n1_idx_ == graph_->i2u_.size()-1 && \
                 adj_list.size() == 0){
          n1_idx_++;
          adj_list_pos_ = 0;
          return *this;
        //if there are no adjacent nodes, to n1, go to the next n1
        }else if(adj_list.size() == 0){
          n1_idx_++;
          adj_list_pos_ = 0;
        //if the index of the adjacent node is < n1, go to next adj_node
        }else if(n1_idx_ > n2_idx && \
                 adj_list_pos_ < adj_list.size() - 1){
          adj_list_pos_++;
        //if we're at the end of the current n1's adjacency nodes
        }else if(n1_idx_ > n2_idx && \
                 n1_idx_ < graph_->i2u_.size()-1 && \
                 adj_list_pos_ == adj_list.size()-1){
          n1_idx_++;
          adj_list_pos_ = 0;
        }

        //update the info for our checks
        n1_uid = graph_->i2u_[n1_idx_];
        adj_list = graph_->nodes_[n1_uid].adj_nodes;
        n2_uid = adj_list[adj_list_pos_].node_id;
        n2_idx = graph_->nodes_[n2_uid].idx_;
      }

      return *this;
    }

    /** Defines equality for incident_iterator */
    bool operator==(const EdgeIterator& x) const{
      if(graph_ == x.graph_)
        if(n1_idx_ == x.n1_idx_ && adj_list_pos_ == x.adj_list_pos_)
          return true;
      return false;
    }
    
   private:
    friend class Graph;
    Graph* graph_;
    size_type n1_idx_;
    size_type adj_list_pos_;

    EdgeIterator(const Graph* graph, size_type n1_idx_, size_type adj_list_pos_)
      : graph_(const_cast<Graph*>(graph)),
        n1_idx_(n1_idx_),
        adj_list_pos_(adj_list_pos_){
    }

  };

  /** Creates beginning iterator to go through all edges of a graph */
  edge_iterator edge_begin() const {

    //find the first non-empty adjacency list, or if there are no edges
    int first_nonempty_node_idx = 0;
    int total_edges = 0;
    for(auto ni = this->node_begin(); !(ni == this->node_end()); ++ni){
      auto node = *ni;
      for(auto ii = node.edge_begin(); ii != node.edge_end(); ++ii)
        total_edges++;

      size_type node_uid = i2u_[node.idx_];
      if(nodes_[node_uid].adj_nodes.size() > 0)
        break;
      first_nonempty_node_idx++;
    }

    //handle the case where all adjacency lists are empty, and return end iter
    if(total_edges == 0){
      return EdgeIterator(this,
                          i2u_.size(),
                          nodes_[i2u_[i2u_.size()-1]].adj_nodes.size());
    }

    //otherwise return first non-empty node index
    return EdgeIterator(this, first_nonempty_node_idx, 0);
  }
  /** Creates ending iterator to go through all edges of a graph */
  edge_iterator edge_end() const{
    return EdgeIterator(this,
                        i2u_.size(),
                        nodes_[i2u_[i2u_.size()-1]].adj_nodes.size());
  }

  /**
    @param[in] n, the node to remove
    @result    success code: 1 on success, 0 if graph doesn't contain node

    Removes a node from the graph. Also removes all edges incident to the node.

    @pre:
      @a n is a node object
      NOTE: Function is able to handle nodes not contained in graph
    @post:
      g.node(i).index() == i for all i with 0 <= i < g.num nodes()
      g.node(n.index()) == n
      new num_nodes() = old num_nodes() - 1
      new num_edges() = old num_edges() - n.degrees()
      For all i where i are adjacent nodes to graph, invalidates edge(n,i) and
        edge(i,n)
      For all i, s.t. n.index() < i: new g.node(i) = old g.node(i+1)

    Complexity:
      O(num_nodes()). Graph assumed sparse: n.degree() << g.num_nodes(),
      so calls depending on n.degree() are not counted toward the complexity
      requirement
  */
  size_type remove_node(const Node& n){
    if(!has_node(n))
      return 0;

    //remove edges incident to the node
    while(nodes_[i2u_[n.idx_]].adj_nodes.size() > 0){
      size_type adj_node_uid = nodes_[i2u_[n.idx_]].adj_nodes[0].node_id;
      remove_edge(node(n.idx_), node(nodes_[adj_node_uid].idx_));
    }

    //erase the node
    i2u_.erase(i2u_.begin() + n.idx_);

    //update the indices of all nodes with higher index
    for(size_type i = n.idx_; i < i2u_.size(); ++i){
      this->nodes_[i2u_[i]].idx_ = i;
    }

    return n.idx_;
  }

  /**
    @param[in] n_it, a node iterator
    @result    a node iterator

    Removes a node from the graph that the @a n_it iterator points to the node
    with index i. Returns an iterator pointing to the new node now occupying
    index i

    @pre
      @a n_it is a node iterator object
      NOTE: Function is able to handle nodes not contained in graph
    @post, complexity requirements all equivalent to:
      size_type remove_node(const Node& n)
  */
  node_iterator remove_node(node_iterator n_it){
    Node n = *n_it;
    size_type idx = remove_node(n); (void) idx;
    return n_it;
  }

  /**
    @param[in] a, the first node of an edge(a,b)
    @param[in] b, the second node of an edge(a,b)
    @result    success code: 1 on success, 0 if graph doesn't contain edge(a,b)

    Removes an edge from the graph. Removes both edge(a,b) and edge(b,a) from
    the graph. Invalidates edges that were instantiated on (a,b).

    @pre:
      @a e is an edge object
      Function is able to handle edges not contained in graph
    @post:
      new num_nodes() = old num_nodes()
      new num_edges() = old num_edges() - 1
      new a.degree() = old a.degree() - 1
      new b.degree() = old b.degree() - 1

    Complexity:
      O(a.degrees() + b.degrees()), which is < O(num_nodes() + num_edges())
      because graph assumed sparse.
  */
  size_type remove_edge(const Node& a, const Node& b){
    //Exit with error code 0 if edge doesnt exist
    if(!has_edge(a,b))
      return 0;

    //delete the edge(a,b) by going through a's adjacency list
    //We use the full path for the adj_nodes vector so that we can manipulate
    //its data in the graph
    for(auto it =  nodes_[i2u_[a.index()]].adj_nodes.begin(); \
             it != nodes_[i2u_[a.index()]].adj_nodes.end(); ++it){

      auto internal_edge_element = *it;
      size_type adj_node_uid = internal_edge_element.node_id;

      if(i2u_[b.index()] == adj_node_uid){
        *it = *(nodes_[i2u_[a.index()]].adj_nodes.end() - 1);
        nodes_[i2u_[a.index()]].adj_nodes.pop_back();
        break;
      }
    }

    //delete the edge(b,a) by going through b's adjacency list
    for(auto it =  nodes_[i2u_[b.index()]].adj_nodes.begin(); \
             it != nodes_[i2u_[b.index()]].adj_nodes.end(); ++it){

      auto internal_edge_element = *it;
      size_type adj_node_uid = internal_edge_element.node_id;

      if(i2u_[a.index()] == adj_node_uid){
        *it = *(nodes_[i2u_[b.index()]].adj_nodes.end() - 1);
        nodes_[i2u_[b.index()]].adj_nodes.pop_back();
        break;
      }
    }

    //Exit with success code
    return 1;
  }

  /**
    @param[in] e, an edge
    @result    success code: 1 on success, 0 if graph doesn't contain edge(a,b)

    Removes an edge from the graph. Removes both edge(e.node1(),node2())
    and edge (e.node2(),e.node1()) from the graph.
    Invalidates edges that were instantiated on @a e.

    @pre:
      @a e is an edge object
      Function is able to handle edges not contained in graph
    @post, complexity requirements all equivalent to:
      size_type remove_edge(const Node& a, const Node& b)
  */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /**
    @param[in] e_it, an edge iterator
    @result    a edge iterator

    Removes an edge from the graph that the @a e_it iterator points to. Returns
    an iterator pointing to the edge that now is in the location that @a e_it
    used to occupy in the order of

    @pre:
      @a e is an edge object
      Function is able to handle edge iterators not contained in graph
    @post, complexity requirements all equivalent to:
      size_type remove_edge(const Node& a, const Node& b)
  */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }

  //Helper function to see what is going on internally in the graph...
  //Useful for making sure that node/edge removal works properly.
  //(it does! all test cases passed in ./test_nodes and ./test_edges)
  void print_internal_structs(){
    std::cout << "VALID INDEXES: [";
    for(size_type i = 0; i < i2u_.size(); i++){
      std::cout << i2u_[i] << ",";
    }
    std::cout << "]" << std::endl << std::endl;

    for(size_type i = 0; i < nodes_.size(); i++){

      std::cout << "uid: " << i << std::endl <<
                   "idx_: " << nodes_[i].idx_ << std::endl <<
                   "location: " << nodes_[i].location << std::endl <<
                   //"value: " << nodes_[i].value << std::endl <<
                   "adj_list" << "[";
      for(size_type j = 0; j < nodes_[i].adj_nodes.size(); j++){
        std::cout << nodes_[i].adj_nodes[j].node_id << ",";
      }
      std::cout << "]" << std::endl << std::endl;
    }

  }

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  struct internal_edge_element{
    size_type node_id;
    edge_value_type value;
  };

  //internal node- like structures
  struct internal_node_element{
    Point location; //The point held by a node
    node_value_type value; //the optional value held by the node
    std::vector<internal_edge_element> adj_nodes; //adjacency list
    size_type idx_; //node's current index
  };

  //an array for node-like structures, holds the data of the graph
  std::vector<internal_node_element> nodes_;

  //an array that maps node_index to node_uid
  //node_index: a value i, 0 <= i < graph_.size() that refers to a valid node
  //node_uid: a value j, 0 <= i < nodes_.size() that refers to an (in)valid node
  std::vector<size_type> i2u_;
};

#endif // CME212_GRAPH_HPP
