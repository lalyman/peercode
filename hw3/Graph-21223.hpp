#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

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

 public:
   typedef V node_value_type;
   typedef E edge_value_type;

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
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node>{
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

    /** return the node value for this node */
    node_value_type& value(){
      // return graph_->node_values_[uid_];
      return graph_->graph_nodes_[uid_].node_value;
    }

    /** return the constant node value for this node */
    const node_value_type& value() const{
      // return graph_->node_values_[uid_];
      return graph_->graph_nodes_[uid_].node_value;
    }

    /** set the node value for this node*/
    // void set_value(node_value_type value_){
    //   graph_->node_values_[uid_] = value_;
    // }

    /** return the degree of this node */
    size_type degree() const{
      return graph_->incidence_[uid_].size();
    }

    /** return an iterator object
    * initially points to the beginning of the list of
    * adjacent edges of the current node */
    incident_iterator edge_begin() const{
      size_type begin_idx = 0;
      return IncidentIterator(graph_, begin_idx, uid_);
    }

    /** return an iterator object
    * points to the end of the list of
    * adjacent edges of the current node */
    incident_iterator edge_end() const{
      size_type end_idx = graph_->incidence_[uid_].size();
      return IncidentIterator(graph_, end_idx, uid_);
    }

    /** Return a reference to this node's position. */
    const Point& position() const {
       return graph_->graph_nodes_[uid_].pos;
    }

    /** Return a reference to this node's position. */
   Point& position(){
       return graph_->graph_nodes_[uid_].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
       return uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        return (uid_ == n.uid_ and graph_ == n.graph_);
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
        bool same_graph = (graph_ == n.graph_ and uid_ < n.uid_);
            return (same_graph || (graph_<n.graph_));
         }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    node_value_type val_;

    // private constructor
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }

    }; // end of node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return graph_nodes_.size();
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
   proxy_node new_node;
   new_node.pos = position;
   new_node.node_idx = graph_nodes_.size();
   new_node.node_value = value;
   graph_nodes_.push_back(new_node);
   // node_values_.push_back(value);

   std::vector<proxy_edge> temp;
   incidence_.push_back(temp);
   return(Node(this, new_node.node_idx));
 }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
   return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** return the value of this edge **/
    edge_value_type& value(){
      return graph_->graph_edges_[uid_].edge_value;
    }

    /** return the value of this edge **/
    const edge_value_type& value() const{
      return graph_->graph_edges_[uid_].edge_value;
    }

    /** return the length of the edge **/
    double length() const{
      return norm(node1().position() - node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      //return (graph_->graph_edges_[uid_].node_1);
      return (graph_->node(node1_));
    }

    /** Return the other node of this Edge */
    Node node2() const {
     //return (graph_->graph_edges_[uid_].node_2);
     return (graph_->node(node2_));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
     return ((node1()==e.node1() && node2()==e.node2())||
          (node2()==e.node1() && node1()==e.node2()));
    }

    size_type get_uid() const{
      return uid_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      bool same_graph = (graph_ == e.graph_ and uid_ < e.uid_);
      return (same_graph || (graph_<e.graph_));
   }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // private class variables
    Graph* graph_;
    size_type uid_;
    size_type node1_;
    size_type node2_;

    // private constructor for the edge class
    Edge(const Graph* graph, size_type uid, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), node1_(node1), node2_(node2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return graph_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
   size_type node1 = graph_edges_[i].node_1.index();
   size_type node2 = graph_edges_[i].node_2.index();
   return Edge(this, i, node1, node2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   if (this != a.graph_ or this != b.graph_){
    return false;
   }

   for (unsigned int i=0; i<incidence_[a.index()].size(); i++){
     if ((incidence_[a.index()].at(i).node_2.index() == b.index())
     ||(incidence_[a.index()].at(i).node_1.index() == b.index())) {
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
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_edge(a,b)){
      for (unsigned int i=0; i<incidence_[a.index()].size(); i++){
        if ((incidence_[a.index()].at(i).node_2.index() == b.index())
      ||(incidence_[a.index()].at(i).node_1.index() == b.index())) {
             return (Edge(this, incidence_[a.index()].at(i).edge_idx, a.index(), b.index()));
        }
      }
    }

     proxy_edge new_edge_a;
     new_edge_a.node_1 = a;
     new_edge_a.node_2 = b;
     new_edge_a.edge_idx = num_edges();
     new_edge_a.edge_value = value;
     graph_edges_.push_back(new_edge_a);

     proxy_edge new_edge_b;
     new_edge_b.node_1 = a;
     new_edge_b.node_2 = b;
     new_edge_a.edge_value = value;
     new_edge_b.edge_idx = num_edges()-1;

     // add the new edge to incidence_
     incidence_[a.index()].push_back(new_edge_a);
     incidence_[b.index()].push_back(new_edge_b);

     return (Edge(this, new_edge_a.edge_idx, a.index(), b.index()));
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    graph_nodes_.clear();
    graph_edges_.clear();
    incidence_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
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

    /** return an NodeIterator object that points to the
    * node with @a uid_*/
    Node operator*() const {
     return graph_->node(uid_);
    }

    /** return an NodeIterator object that points to the
    * node with @a uid_+1 */
    node_iterator& operator++(){
     uid_ = uid_+1;
     return *this;
    }

    /** check if this NodeIterator object is equal to the
    input NodeIterator object_*/
    bool operator==(const node_iterator& x) const{
     return (uid_ == x.uid_ and graph_ == x.graph_);
    }

   private:
    friend class Graph;
    Graph* graph_; // pointer to a graph object
    size_type uid_;

    // constructor for the node iterator
    NodeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  }; // end of node iterator class


  /** return an NodeIterator object that points to the
  * node with index 0 */
  node_iterator node_begin() const{
   size_type begin_idx = 0;
   return NodeIterator(this, begin_idx);
  }

  /** return an NodeIterator object that points to the
  * node with index equal to the number of nodes in the graph */
  node_iterator node_end() const{
   size_type end_idx = graph_nodes_.size();
   return NodeIterator(this, end_idx);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator>{
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

    Edge operator*() const{
      assert(nid_< graph_->num_nodes());
      if (graph_->incidence_[nid_].at(uid_).node_1.index() == nid_){
        return Edge(graph_, graph_->incidence_[nid_].at(uid_).edge_idx, nid_, graph_->incidence_[nid_].at(uid_).node_2.index());
      }
      else{
        return Edge(graph_, graph_->incidence_[nid_].at(uid_).edge_idx, nid_, graph_->incidence_[nid_].at(uid_).node_1.index());
       //return graph_->edge(graph_->incidence_[nid_].at(uid_).edge_idx);
    }
  }

    /** return an IncendentIterator object that points to the
    * node with index uid_+1 */
    incident_iterator& operator++(){
      ++uid_;
      return *this;
    }

    /** check if two incident iterator objects are equal */
    bool operator==(const incident_iterator& iit) const{
      return (uid_ == iit.uid_ and graph_ == iit.graph_ and nid_ == iit.nid_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type uid_; //iterator
    size_type nid_; //node id

    // constructor for the IncidentIterator
    IncidentIterator(const Graph* graph, size_type uid, size_type nid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), nid_(nid) {
    }
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

    /** return an EdgeIterator object that points to the
    * edge with index uid_ */
    Edge operator*() const{
      return graph_->edge(graph_->graph_edges_[uid_].edge_idx);
    }

    /** return an EdgeIterator object that points to the
    * edge with index uid_+1 */
    EdgeIterator& operator++(){
      ++uid_;
      return *this;
    }

    /** check if two EdgeIterator objects are equal */
    bool operator==(const EdgeIterator& x) const{
      return (x.graph_==graph_ && x.uid_==uid_);
    }


   private:
    friend class Graph;
    Graph* graph_;
    size_type uid_;

    /** Constructor for the EdgeIterator object */
    EdgeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** return an EdgeIterator object that points to the
  * edge with index 0 */
  edge_iterator edge_begin() const{
    size_type begin_idx = 0;
    return edge_iterator(this, begin_idx);
  }

  /** return an EdgeIterator object that points to the
  * edge with index equal to total number of edges in graph */
  edge_iterator edge_end() const{
    size_type end_idx = num_edges();
    return edge_iterator(this, end_idx);
  }

  /**
   * @brief Remove the edge (@a a,@a b) from the graph
   *
   * @param a, node in the graph
   * @param b, node in the graph
   * @return  Edge index of the removed edge or the number of edges
   *
   * @pre has_edge(a,b)
   * @post new_graph_edges_.size() = old_graph_edges.size() - 1
   * @post incidence_[a.index()] and incidence_[b.index()] does not contain the
   *       index of the edge removed
   * @post new graph_edges_[index] = old graph_edges_.back()
   *
   * The edge connecting nodes a and b is invalidated.
   *
   * The complexity of remove_edge is O(num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b){
    if (has_edge(a,b)){
      size_type edge_index = 0;

      for (auto iter = graph_edges_.begin(); iter < graph_edges_.end(); ++iter){
        auto temp = *iter;
        if ((a.index() == temp.node_1.index() && b.index() == temp.node_2.index()) ||
          (b.index() == temp.node_1.index() && a.index() == temp.node_2.index())){
            edge_index = iter - graph_edges_.begin();
            break;
          }
        }

      size_type num_edges = graph_edges_.size();
      // replace the "removed edge" with the last edge, and pop the last edge
      graph_edges_[edge_index] = graph_edges_.back();
      // graph_edges_[edge_index].edge_idx = graph_edges_.back().edge_idx;
      graph_edges_[edge_index].edge_idx = edge_index;

      graph_edges_.pop_back();

      auto& adj_a = incidence_[a.index()];
        for (unsigned int i = 0; i != adj_a.size(); ++i){
          auto current = adj_a[i];
            if ((a == current.node_1 && b == current.node_2)||
              (b == current.node_1 && a == current.node_2)) {
                incidence_[a.index()][i] = incidence_[a.index()].back();
                incidence_[a.index()].pop_back();
                // incidence_[a.index()].erase(adj_a.begin()+i);
              // for (unsigned int j = 0; j!=adj_a.size(); j++){
              //   if (adj_a[j].edge_idx>adj_a[i].edge_idx){
              //     adj_a[j].edge_idx -= 1;
              //   }
              //   }
                break;
            }
          }

      auto& adj_b = incidence_[b.index()];
          for (unsigned int i = 0; i != adj_b.size(); ++i){
            auto current = adj_b[i];
            if ((a == current.node_1 && b == current.node_2)||
                (b == current.node_1 && a == current.node_2)){
                  incidence_[b.index()][i] = incidence_[b.index()].back();
                  incidence_[b.index()].pop_back();
                  // incidence_[b.index()].erase(adj_b.begin()+i);
                  // for (unsigned int j = 0; j!=adj_b.size(); j++){
                  //   if (adj_b[j].edge_idx>adj_b[i].edge_idx){
                  //     adj_b[j].edge_idx -= 1;
                  //   }
                  // }
                break;
              }
            }

      for (auto it=incidence_[graph_edges_[edge_index].node_1.index()].begin(); it !=incidence_[graph_edges_[edge_index].node_1.index()].end(); ++it) {
        if ((*it).edge_idx == num_edges-1) {
          (*it).edge_idx = edge_index;
          break;
        }
      }

      for (auto it=incidence_[graph_edges_[edge_index].node_2.index()].begin(); it !=incidence_[graph_edges_[edge_index].node_2.index()].end(); ++it) {
        if ((*it).edge_idx == num_edges-1) {
          (*it).edge_idx = edge_index;
          break;
        }
      }


      // change the edge indexs in the incidence_ vector after removing an edge
      // for (unsigned int i=0; i!=incidence_.size(); ++i){
      //   for (unsigned int j=0; j!=incidence_[i].size(); ++j){
      //     if (incidence_[i][j].edge_idx == num_edges-1) {
      //       incidence_[i][j].edge_idx = edge_index;
      //     }
      //   }
      // }


      return edge_index;
      }
        else {
          return num_edges();
        }
    }

  /**
   * @brief Remove the edge @e from the graph
   *
   * @param e, edge in the graph
   * @return  Edge index of the removed edge or the number of edges
   *
   * @pre has_edge(e)
   * @post new_graph_edges_.size() = old_graph_edges_.size() - 1
   * @post incidence_[a.index()] and incidence_[b.index()] does not contain the
   *       index of the edge removed
   * @post new graph_edges_[index] = old graph_edges_.back()
   *
   * The edge connecting nodes a and b is invalidated.
   *
   * The complexity of remove_edge is O(num_edges())
   */

  size_type remove_edge(const Edge& e){
    size_type idx = remove_edge(e.node1(), e.node2());
    return idx;
  }

  /**
   * @brief Remove the edge @a *e_it from the graph
   *
   * @param e_it, an edge iterator
   * @return  Edge index of the removed edge or the number of edges
   *
   * @pre has_edge(*e_it)
   * @post new_graph_edges_.size() = old_graph_edges_.size() - 1
   * @post incidence_[a.index()] and incidence_[b.index()] does not contain the
   *       index of the edge removed
   * @post new graph_edges_[index] = old graph_edges_.back()
   *
   * The edge connecting nodes a and b is invalidated.
   *
   * The complexity of remove_edge is O(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }


  /**
   * @brief Remove the edge @a nfrom the graph
   *
   * @param n, node in the graph
   * @return  node index of the removed edge or the number of edges
   *
   * @pre has_node(n)
   * @post new_incidence_[n.index()] = old_incidence_[n.index()].back();
   * @post new incidence_.size() = old_incidence_.size()-1
   * @post new graph_nodes_[index] = old graph_nodes_.back()
   *
   * The edge connecting nodes a and b is invalidated.
   *
   * The complexity of remove_edge is O(num_edges())
   */
  size_type remove_node(const Node& n){
    if (!has_node(n)){
      return num_nodes();
    }

  // remove all edges associated with n in the incidence_ vector
  // auto temp = incidence_[n.index()];
  // for (unsigned int iter = 0;
  //     iter != temp.size(); ++iter) {
  while (! incidence_[n.index()].empty()){
    // std::cout << incidence_[n.index()].size() << std::endl;
    auto current = incidence_[n.index()].back();
    auto node1 = current.node_1;
    auto node2 = current.node_2;
    remove_edge(node1, node2);
  }


  size_type old_idx = n.index();
  graph_nodes_[old_idx] = graph_nodes_.back();
  unsigned int old_length = graph_nodes_.size();
  // graph_nodes_[old_idx].node_idx = graph_nodes_.back().node_idx;
  graph_nodes_[old_idx].node_idx = n.index();
  graph_nodes_.pop_back();

  incidence_[old_idx] = incidence_.back();
  incidence_.pop_back();

 //  for (unsigned int iter = 0; iter < incidence_[old_idx].size(); ++iter){
 //   auto current = incidence_[old_idx][iter];
 //   if (current.node_1.index() == old_length){
 //     current.node_1 = node(old_idx);
 //   }
 //   else if (current.node_2.index() == old_length){
 //     current.node_2 = node(old_idx);
 //   }
 // }

 // change the edge indexs in the incidence_ vector after removing an edge
 for (unsigned int i=0; i!=incidence_.size(); ++i){
   for (unsigned int j=0; j!=incidence_[i].size(); ++j){
     if (incidence_[i][j].node_1.index() == old_length-1) {
       incidence_[i][j].node_1 = node(old_idx);
     }
     else if (incidence_[i][j].node_2.index() == old_length-1) {
       incidence_[i][j].node_2 = node(old_idx);
     }
   }
 }

 for (unsigned int i=0; i<graph_edges_.size(); ++i){
   if (graph_edges_[i].node_1.index() == old_length-1){
     graph_edges_[i].node_1 = node(old_idx);
   }
   else if (graph_edges_[i].node_2.index() == old_length-1){
     graph_edges_[i].node_2 = node(old_idx);
   }
 }

  return n.index();
}


/**
 * @brief Remove the node @a *n_it from the graph
 *
 * @param n_it, node iterator in the graph
 * @return  node index of the removed edge or the number of edges
 *
 * @pre has_node(*n_it)
 * @post new_incidence_[*n_it.index()] = old_incidence_[*n_it.index()].back();
 * @post new incidence_.size() = old_incidence_.size()-1
 * @post new graph_nodes_[index] = old graph_nodes_.back()
 *
 * The edge connecting nodes a and b is invalidated.
 *
 * The complexity of remove_edge is O(num_edges())
 */
node_iterator remove_node(node_iterator n_it){
  remove_node((*n_it));
  return n_it;
}

  private:
  struct proxy_node{
   Point pos;
   size_type node_idx;
   node_value_type node_value;
  };

  struct proxy_edge{
   Node node_1;
   Node node_2;
   size_type edge_idx;
   edge_value_type edge_value;
  };


  std::vector<proxy_node> graph_nodes_; // stores nodes of the graph as vector of points
  std::vector<proxy_edge> graph_edges_;
  // std::vector<node_value_type> node_values_;
  // the i-th element of this vector contains the edges (i,j) in the edge set
  std::vector<std::vector<proxy_edge>> incidence_;

};

#endif // CME212_GRAPH_HPP
