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
    // The data structure has been updated based on HW2 comments.
    // Reference -- peercode/hw2/Graph-31266.hpp
    // Due to time constraints, grader feedback regarding
    // documentation was not addressed.
    struct internal_node {
      Point node_position;
      node_value_type node_val;
      size_type node_uid;
    };

    struct internal_edge {
      edge_value_type edge_val;
      size_type n2_uid;
    };

    std::vector<internal_node> node_elements_;
    std::vector<std::vector<internal_edge>> incident_edge_;
    std::vector<size_type> node_valid_;

    size_type edge_count_ = 0;      // track the edge uid
    size_type edge_idx_ = 0;        // track the edge index
    
  
  public:
  
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
  
    // CONSTRUCTORS AND DESTRUCTOR
    /** Construct an empty graph. */
    Graph(): node_elements_(), incident_edge_(), node_valid_(){ }
  
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
    size_type index() const { return fetch_node().node_uid; }

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
                const node_value_type& node_val_in = node_value_type()){

    // construct a new internal node
    internal_node new_node;
    new_node.node_position = position;
    new_node.node_val = node_val_in;
    new_node.node_uid = node_valid_.size();

    // add new node to vector of internal nodes
    // and vector of valid modes
    node_valid_.push_back(node_elements_.size());
    node_elements_.push_back(new_node);

    // initialize adjacency matrix
    std::vector<internal_edge> v {};
    incident_edge_.push_back(v);

    return Node(this, node_elements_.size() - 1);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this != n.graph_) return false;
    if (n.index() < this -> num_nodes()) return true;
    return (this == n.graph_ and node_valid_[n.index()] == n.node_uid_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
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

    auto rm_uid = n.node_uid_;
    for(size_type i = 0; i < incident_edge_[rm_uid].size(); ++i){
      size_type uid2 = incident_edge_[rm_uid][i].n2_uid;
      for (size_type j = 0; j < incident_edge_[uid2].size(); j++){
        if(incident_edge_[uid2][j].n2_uid == rm_uid) {
          incident_edge_[uid2][j] = incident_edge_[uid2].back();
          incident_edge_[uid2].pop_back();
        }
      }
    }
    incident_edge_[rm_uid].clear();

    size_type idx = n.index();
    node_elements_[node_valid_[node_valid_.size() - 1]].node_uid = idx;

    node_valid_[idx] = node_valid_.back();
    node_valid_.pop_back();

    return idx;

  }

  node_iterator remove_node(node_iterator n_it){
    auto idx = remove_node(*n_it);
    return NodeIterator(this,idx);
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

    /** Return the value associated with an edge. */
    edge_value_type& value() { return fetch_edge().edge_val; }

    /** Return the value associated with an edge. */
    const edge_value_type& value() const { return value(); }

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
    size_type uid1_, uid2_;

    Edge(const graph_type* graph, size_type uid1, size_type uid2):
        graph_(const_cast<graph_type*>(graph)), uid1_(uid1), uid2_(uid2) {}

    // loop through all the nodes and the incident edges
    // to find the correct edge to return
    internal_edge& fetch_edge() const {
      if(uid1_ < uid2_) {
        for(size_type i = 0; i < graph_->incident_edge_[uid1_].size(); ++i)
          if(graph_->incident_edge_[uid1_][i].n2_uid == uid2_)
            return graph_->incident_edge_[uid1_][i];
      }
      for(size_type i = 0; i < graph_->incident_edge_[uid2_].size(); ++i)
        if(graph_->incident_edge_[uid2_][i].n2_uid == uid1_)
         return graph_->incident_edge_[uid2_][i];
      assert(false);
    }

  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { 
    return std::distance(edge_begin(), edge_end()); 
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i <= this -> num_edges());
    return *(std::next(edge_begin(),i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
    size_type n1_uid = a.node_uid_;
    size_type n2_uid = b.node_uid_;
    if (this == a.graph_ && this == b.graph_){
        for (size_type i = 0; i < incident_edge_[n1_uid].size(); ++i){
            if (incident_edge_[n1_uid][i].n2_uid == n2_uid){
                return true;
            }
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
                const edge_value_type& edge_val_in = edge_value_type()) {
    size_type n1_uid = a.index();
    size_type n2_uid = b.index();

    if (has_edge(a,b)) {
        return Edge(this, n1_uid, n2_uid);
    } else {
        // construct a new internal_edge
        internal_edge new_edge;
        new_edge.n2_uid = n2_uid;
        new_edge.edge_val = edge_val_in;
        incident_edge_[n1_uid].push_back(new_edge);

        new_edge.n2_uid = n1_uid;
        incident_edge_[n2_uid].push_back(new_edge);

        return Edge(this, n1_uid, n2_uid);
    }
    return Edge();
  }

  /** Invalidate an edge by restricting access to the internal_node
   *  @param[in] n The node to be removed
   *  @post new num_edges() == old num_edgess() - 1
   *
   *  @return edge index of the next valid node in @a edge_valid_
   */

  size_type remove_edge(const Node& a, const Node& b){

    // only remove valid edges
    if(has_edge(a,b) == false)
        return a.node_uid_ < b.node_uid_ ? a.degree() : b.degree();
    
    auto uid1 = a.node_uid_;
    auto uid2 = b.node_uid_;
    size_type idx1 = 0;
    size_type idx2 = 0;

    
    bool found_edge = false;
    while(!found_edge && idx1 < incident_edge_[uid1].size()) {
      if(incident_edge_[uid1][idx1].n2_uid == uid2) {
        incident_edge_[uid1][idx1] = incident_edge_[uid1].back();
        incident_edge_[uid1].pop_back();
        found_edge = true;
      }
      ++idx1;
    }

    // now the other other edge
    found_edge = false;
    while(!found_edge && idx2 < incident_edge_[uid2].size()) {
      if(incident_edge_[uid2][idx2].n2_uid == uid1) {
        incident_edge_[uid2][idx2] = incident_edge_[uid2].back();
        incident_edge_[uid2].pop_back();
        found_edge = true;
      }
      ++idx2;
    }
    return uid1 < uid2 ? idx1 : idx2;
  }

  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it){
    auto next_idx1 = remove_edge(*e_it);
    auto next_idx2 = *e_it.node1().node_uid_ < *e_it.node2().node_uid_ ? 
        *e_it.node1().node_uid_ : *e_it.node2().node_uid_;
    return EdgeIterator(this,next_idx2,next_idx1);
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
    Edge operator*() const{ 
        return Edge(graph_, uid1_, graph_->incident_edge_[uid1_][uid2_].n2_uid);
    }

    /** Increase the iterator by one, reference -- peercode-hw1-4943*/
    EdgeIterator& operator++(){
      if(uid1_ < graph_->node_elements_.size())
        uid2_++;
      find_valid_edge();
      return *(this);
    }


    /** Compare two edge iterators */
    bool operator==(const EdgeIterator& eit) const {
      return (graph_ == eit.graph_ && uid1_ == eit.uid1_ && uid2_ == eit.uid2_);
    }

   private:
     Graph* graph_;
     size_type uid1_;
     size_type uid2_;
     EdgeIterator(const Graph* graph, size_type uid1 = 0, size_type uid2 = 0)
       : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
         find_valid_edge();
     }
     friend class Graph;
     void find_valid_edge()
     {
       while(uid1_ < graph_-> incident_edge_.size())
       {
         while(uid2_ < graph_-> incident_edge_[uid1_].size())
         {
           if(uid1_ < graph_-> incident_edge_[uid1_][uid2_].n2_uid)
             return;
           ++uid2_;
         }
         ++uid1_;
         uid2_ = 0;
       }
     }

  };

  /** Constructs the iterator that points to the first element in
   *  a vector of edges. */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  /** Constructs the iterator that points to one past the last element in
   *  a vector of edges. */
  edge_iterator edge_end() const {
    return edge_iterator(this, incident_edge_.size());
  }


};

#endif // CME212_GRAPH_HPP

