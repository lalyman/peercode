#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 *  @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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
class Graph{

  public:
    using node_value_type = V;
    using edge_value_type = E;
    using size_type = unsigned;
    using graph_type = Graph<node_value_type, edge_value_type>;
  
    /** Predeclaration of Node and Edge type. */
    /** Synonym for Node and Edge (following STL conventions). */
    class Node;
    using node_type = Node;
  
    class Edge;
    using edge_type = Edge;
  
    /** Type of node iterators, which iterate over all graph nodes. */
    class NodeIterator;
    using node_iterator = NodeIterator;
  
    /** Type of edge iterators, which iterate over all graph edges. */
    class EdgeIterator;
    using edge_iterator = EdgeIterator;
  
    /** Type of incident iterators, which iterate incident edges to a node. */
    class IncidentIterator;
    using incident_iterator = IncidentIterator;
  
  private:

    struct internal_node {
      Point node_position;
      node_value_type node_val;
      size_type node_idx;
    };

    struct internal_edge {
      edge_value_type edge_val;
      size_type n2_uid;
      size_type edge_uid;
    };

    std::vector<internal_node> node_elements_;
    std::vector<std::vector<internal_edge>> incident_edge_;
    std::vector<size_type> node_valid_;
    std::vector<size_type> edge_valid_;

    size_type edge_count_ = 0;      // track the edge uid
    size_type edge_idx_ = 0;        // track the edge index
    size_type tot_num_edge_ = 0;    // track total number of edges
    bool edge_removed_ = false;     // track if edge_remove() has been called
    
  
  public:
  
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
  
    // CONSTRUCTORS AND DESTRUCTOR
    /** Construct an empty graph. */
    Graph(): node_elements_(), incident_edge_(), node_valid_(), edge_valid_(){ }
  
    /** Default right_uidructor */
    ~Graph() = default;


  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class reXpresenting the graph's nodes.
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
    Node() { }

    /** Return this node's position. */
    Point& position() const { return fetch_node().node_position; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return fetch_node().node_idx; }

    /** Return the value associated with a node. */
    node_value_type& value() { return fetch_node().node_val; }

    /** Return the value associated with a node as a const. */
    const node_value_type& value() const { return fetch_node().node_val; }


    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_ == this->graph_ and n.node_uid_ == this->node_uid_);
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
      // compare graph location 
      if (this->graph_ < n.graph_) return true;

      if (this->node_uid_ < n.node_uid_) return true;
      return false;
    }
    
    /** Returns the number of incident nodes to the node of interest.*/
    size_type degree() const{
      return graph_->incident_edge_[node_uid_].size(); 
    }
    
    /** Constructs the iterator that points to the first element in 
     *  a vector of incident nodes to the current node of interest.*/
    incident_iterator edge_begin() const {
        return incident_iterator(graph_, node_uid_, 0);
    }

    /** Constructs the iterator that points to one past the last element in 
     *  a vector of incident nodes to the current node of interest.*/
    incident_iterator edge_end() const {
        return incident_iterator(graph_, node_uid_, this->degree());
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_uid_;

    Node(const graph_type* graph, size_type node_uid):
        graph_(const_cast<graph_type*>(graph)), node_uid_(node_uid){}

    internal_node& fetch_node() const {
        assert(node_uid_ < graph_->node_elements_.size());
        return graph_->node_elements_[node_uid_];
    }

  };


  /** Return the number of nodes in the graph. */
  size_type size() const {
    return node_valid_.size();
  }
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& node_val_in) {

    auto node_uid_in = node_elements_.size();
    node_valid_.push_back(node_uid_in);

    auto node_idx_in = node_valid_.size() - 1;
    internal_node new_node;
    new_node.node_position = position;
    new_node.node_val = node_val_in;
    new_node.node_idx = node_idx_in;
    
    node_elements_.push_back(new_node);

    // initialize adjacency matrix
    std::vector<internal_edge> v {};
    incident_edge_.push_back(v);

    return Node(this, node_uid_in);

  }

  Node add_node(const Point& position) {
    const node_value_type& node_val = node_value_type();
    return add_node(position, node_val);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this != n.graph_) return false;
    if (n.index() < this -> num_nodes()) return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, node_valid_[i]); 
  }

  /** Invalidate a node by restricting access to the internal_node
   *  @param[in] n The node to be removed
   *  @post new num_nodes() == old num_nodes() - 1
   *  @post new num_edges() == 
   *    old num_edges() - @a incident_edge_[@a n.uid].size()
   *
   *  @return node index of the next valid node in @a node_valid_
   */
  size_type remove_node(const Node& n){
    // get the index of the node to be removed
    auto rm_idx = n.index();
    
    // remove the associated edges
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i ){
       remove_edge(*i); 
    }

    // remove the associated node
    auto new_itr = node_valid_.erase(node_valid_.begin() + rm_idx);

    // update the idx for all internal_nodes
    for (; new_itr < node_valid_.end(); ++new_itr){
        node_elements_[*new_itr].node_idx = new_itr - node_valid_.begin();
    }

    return *new_itr;
  }

  node_iterator remove_node(node_iterator n_it){
    auto new_node = remove_node(*n_it);
    return node_specific(new_node);
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge(){}

    /** Return a node of this Edge */
    Node node1() const { return Node(graph_,uid1_); }

    /** Return the other node of this Edge */
    Node node2() const { return Node(graph_,uid2_); }

    /** Return idx of this Edge */
    size_type index() const {return edge_uid_;}

    /** Return the value associated with an edge. */
    edge_value_type& value() { return fetch_edge().edge_val; }

    /** Return the value associated with an edge. */
    const edge_value_type& value() const { return fetch_edge().edge_val; }

    /** Return edge length*/
    double length() const {
        return norm(this->node1().position() - this->node2().position());
    }

    /** Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if (e.graph_ != graph_) return false;

      if  (e.node1().index() == uid1_  and
           e.node2().index() == uid2_){
          return true;
      } else if
          (e.node1().index() == uid2_  and
           e.node2().index() == uid1_){
          return true;
      } else {
          return false;
      }

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

       // if two edges are the same, return false
       if (*this == e) return false;

       // compare graph location 
       if (this->graph_ < e.graph_) return true;

       // compare node index
       if (uid1_ < e.node1().index()) return true;
       if ((uid1_== e.node1().index())
       and (uid2_ < e.node2().index())) return true;
        return false;
    }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type edge_uid_, uid1_, uid2_;

    Edge(const graph_type* graph, size_type edge_uid, size_type uid1, size_type uid2):
        graph_(const_cast<graph_type*>(graph)),
        edge_uid_(edge_uid), uid1_(uid1), uid2_(uid2) {}

    // loop through all the nodes and the incident edges
    // to find the correct edge to return
    internal_edge& fetch_edge() const {
        if (uid1_ < graph_->tot_num_edge_ and uid2_ < graph_->tot_num_edge_){
            for (size_type i = 0; i < graph_->incident_edge_[uid1_].size(); ++i){
                if(graph_->incident_edge_[uid1_][i].n2_uid == uid2_)
                    return graph_->incident_edge_[uid1_][i];
            }
        }
        assert(false);
    }

  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { 
    return this->edge_valid_.size(); 
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i <= this -> num_edges());
    auto e_uid = this->edge_valid_[i];
    for (size_type j = 0; j < incident_edge_.size(); j++){
        for (size_type k = 0; k < incident_edge_[j].size(); k++){
            if (incident_edge_[j][k].edge_uid == e_uid)
            return Edge(this, incident_edge_[j][k].edge_uid,j,
                        incident_edge_[j][k].n2_uid);
        }
    }
    return Edge(); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
   if (edge_exist(a,b))
     return this->valid(this->edge_count_);
   return false;
  }

  /** Returns true if the edge has been constructed
   *  regardless of the validity of the edge
   */
  bool edge_exist(const Node& a, const Node& b){
    //for (incident_iterator i = a.edge_begin(); i != a.edge_end(); ++i){
    auto n1_uid = node_valid_[a.index()];
    auto n2_uid = node_valid_[b.index()];

    for (size_type i = 0; i < incident_edge_[n1_uid].size(); ++i){
        if (incident_edge_[n1_uid][i].n2_uid == n2_uid){
            this->edge_count_ = incident_edge_[n1_uid][i].edge_uid;
            return (true);
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
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& edge_val_in) {
    // update the connectivity of that node
    size_type n1_idx = a.index();
    size_type n2_idx = b.index();
    size_type n1_uid = node_valid_[n1_idx];
    size_type n2_uid = node_valid_[n2_idx];

    if (edge_exist(a,b)) {
        auto edge_uid = this->edge_count_;

        // check to see if edge removal has happened
        if (this->edge_removed_ == true){
            if ((this->valid(edge_uid)) == false){
                this->edge_valid_.push_back(edge_uid);
            }
        }

        return Edge(this, edge_uid, n1_uid, n2_uid);
    } else {

        // get new edge_uid
        auto edge_uid = this->tot_num_edge_;
        this->edge_valid_.push_back(edge_uid);
        this-> tot_num_edge_++;

        // construct a new internal_edge
        internal_edge new_edge1;
        new_edge1.edge_uid = edge_uid;
        new_edge1.n2_uid = n2_uid;
        new_edge1.edge_val = edge_val_in;
        incident_edge_[n1_uid].push_back(new_edge1);

        internal_edge new_edge2;
        new_edge2.edge_uid = edge_uid;
        new_edge2.n2_uid = n1_uid;
        new_edge2.edge_val = edge_val_in;
        incident_edge_[n2_uid].push_back(new_edge2);

        return Edge(this, edge_uid, n1_uid, n2_uid);
    }
    return Edge();
  }
  Edge add_edge(const Node& a, const Node& b) {
    const edge_value_type& edge_val_in = edge_value_type();
    return add_edge(a,b,edge_val_in);
  }

  /** Invalidate an edge by restricting access to the internal_node
   *  @param[in] n The node to be removed
   *  @post new num_edges() == old num_edgess() - 1
   *
   *  @return edge index of the next valid node in @a edge_valid_
   */
  size_type remove_edge(const Edge& e){
    // get the unique index of the edge to be removed
    auto rm_uid = e.index();
    assert(rm_uid < this->tot_num_edge_);    
    // only remove valid edges
    if(this->valid(rm_uid)){
        auto new_itr = edge_valid_.erase(edge_valid_.begin() + this->edge_idx_);
        this->edge_removed_ = true;
        return *new_itr;
    }
    return 0;
  }

  size_type remove_edge(const Node& a, const Node& b){
    if (has_edge(a,b)){
        auto new_itr = edge_valid_.erase(edge_valid_.begin() + this->edge_idx_);
        this->edge_removed_ = true;
        return *new_itr;
    }
    return 0;
  }

  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return edge_begin();
  }

  /** Takes in an unique id of an edge and checks to see if the edge
   *  is still valid.
   */
  /** Check to see if an edge is valid
   *  @param[in] e_uid Unique id of the edge in question
   *
   *  @return
   *    True  if @a e_uid == edge_valid_[i], where 0<= i < num_edges()
   *    False if @a e_uid != edge_valid_[i], where 0<= i < num_edges()
   */
  bool valid(size_type e_uid){
 
    // unique id needs to be smaller than 
    //the total number of edges ever constructed 
    if (e_uid >= this->tot_num_edge_) return false;
    
    // if there is no valid edges
    if (int(this->edge_valid_.size())< 1){ return false;}

    // if uid exists in the vector of valid edge uids
    // then the edge is valid
    for (size_type i = 0; i < this->edge_valid_.size(); i++){
        if (this->edge_valid_[i] == e_uid) {
            this->edge_idx_ = i;
            return true;
        }
    }

    return false;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_elements_.clear();
    incident_edge_.clear();
    node_valid_.clear();
    edge_valid_.clear();
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

    //size_type p;
    /** Construct an invalid NodeIterator. */
    NodeIterator() { }

    /** De-reference a node_iterator, returns the Node the iterator points to.*/
    Node operator*() const {return graph_->node(p_); }
    
    /** Increase the node_iterator by one*/
    NodeIterator& operator++() { ++p_; return *this; }
    
    /** Compares two node_iterators*/
    bool operator==(const NodeIterator& rhs) const { 
        return (p_ == rhs.p_ and graph_ == rhs.graph_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type p_;
    
    NodeIterator(const graph_type* graph, size_type point_idx):
        graph_(const_cast<graph_type*>(graph)), p_(point_idx){}
  };

  /** Constructs the iterator that points to the first element in 
   *  a vector of nodes.
   */
  node_iterator node_begin() const { return node_iterator(this,0); }
  
  /** Constructs the iterator that points to one past the last element in 
   *  a vector of nodes.
   */
  node_iterator node_end() const {
    return node_iterator(this,node_valid_.size());
  }

  /** Constructs the iterator that points to a specific element in 
   *  a vector of nodes.
   */
  node_iterator node_specific(size_type new_node) const {
    return node_iterator(this,new_node); 
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
    IncidentIterator() { }

    /** De-reference an incident iterator, returns the Edge the iterator points to.*/
    Edge operator*() const {
        return Edge(graph_, 
                    graph_->incident_edge_[n1_uid_][p_].edge_uid, 
                    n1_uid_, 
                    graph_->incident_edge_[n1_uid_][p_].n2_uid);
    }

    /** Increase the incident iterator by one*/
    IncidentIterator& operator++() {
        ++p_;
        return *this;
    }

    /** Compares two incident iterators*/
    bool operator==(const IncidentIterator& iit) const {
        return (graph_ == iit.graph_ and 
                n1_uid_ == iit.n1_uid_ and 
                p_ == iit.p_ );}

   private:
    friend class Graph;
    graph_type* graph_;
    size_type n1_uid_;
    size_type p_;
    
    IncidentIterator(const graph_type* graph, size_type n1_uid, size_type p):
        graph_(const_cast<graph_type*>(graph)), n1_uid_(n1_uid), p_(p){}
  };

 

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
    EdgeIterator() { }

    /** De-reference an edge iterator, returns the Edge the iterator points to.*/
    Edge operator*() const{ return *e_; }

    /** Increase the iterator by one, reference -- peercode-hw1-4943*/
    EdgeIterator& operator++(){
      ++e_;
      while (n_ != graph_-> node_end()) {
        while (e_ != (*n_).edge_end()){
          if ((*e_).node1().index() < (*e_).node2().index()){
            return *this;
          }
          ++e_;
        }
        ++n_;
        if (n_ != graph_->node_end()){
          e_ = (*n_).edge_begin();
        }
      }
      return *this;
    }


    /** Compare two edge iterators */
    bool operator==(const EdgeIterator& eit) const {
        return (e_ == eit.e_ and n_ == eit.n_ and graph_ == eit.graph_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    node_iterator n_;
    incident_iterator e_;

    EdgeIterator(const graph_type* graph, node_iterator n, incident_iterator e):
        graph_(const_cast<graph_type*>(graph)), n_(n), e_(e){}
  };

  /** Constructs the iterator that points to the first element in
   *  a vector of edges. */
  edge_iterator edge_begin() const {
    return edge_iterator(this,node_begin(),node(0).edge_begin());
  }

  /** Constructs the iterator that points to one past the last element in
   *  a vector of edges. */
  edge_iterator edge_end() const {
    return edge_iterator(this,node_end(), node(num_nodes()-1).edge_end());
  }


};

#endif // CME212_GRAPH_HPP




