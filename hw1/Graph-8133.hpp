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
template <typename V>
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
      size_type node_id = this->uid_;
      return graph_->nodes_[node_id].location;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->uid_;
    }

    /** Return this node's value, of templated type V */
    node_value_type& value(){
      return graph_->nodes_[uid_].value;
    }

    /** Return this node's value, of templated type V */
    const node_value_type& value() const{
      return graph_->nodes_[uid_].value;
    }

    /** Return the number of nodes connected by an edge to this node */
    size_type degree() const{
      return graph_->nodes_[uid_].adjacent_node_uids.size();
    };

    /** Create a beginning iterator for the edges incident to this node*/
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    };

    /** Create an ending iterator for the edges incident to this node*/
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(this->graph_ == n.graph_)
        if(this->uid_ == n.uid_)
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
      if(this->uid_ < n.uid_)
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
    size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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

    internal_node_element node_to_add = {
      position,
      value,
      std::vector<size_type>(),
    };

    nodes_.push_back(node_to_add);
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
    if(node(n.uid_) == n)
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
      return graph_->node(node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_uid);
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
      //This now checks for edge inequality GLOBALLY since nodes check equality
      //globally
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
    size_type node1_uid;
    size_type node2_uid;

    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
      : graph_(const_cast<Graph*>(graph)),
        node1_uid(node1_uid),
        node2_uid(node2_uid){
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
    auto ei = this->edge_begin();
    for(size_type counter = 0; counter < i; ++counter)
      ++ei;
    Edge result = *ei;
    return Edge(this, result.node1().index(), result.node2().index());
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
    size_type node1_uid = a.index();
    size_type node2_uid = b.index();
    auto adjacent_node_uids = nodes_[node1_uid].adjacent_node_uids;

    for(size_type i = 0; i < adjacent_node_uids.size(); i++){
      size_type adj_node_uid = adjacent_node_uids[i];
      if(adj_node_uid == node2_uid)
        return Edge(this, node1_uid, node2_uid);
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
      nodes_[a.index()].adjacent_node_uids.push_back(b.index());
      nodes_[b.index()].adjacent_node_uids.push_back(a.index());
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
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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

    NodeIterator(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), current_uid(uid) {
    }

    /** De-references a node iterator to the current node*/
    Node operator*() const {
      //This is OK since if the user de-references g.node_end() there
      //should be an error, so don't need to do any checking for that
      return graph_->node(this->current_uid);
    }

    /** Defines equality for node_iterator */
    bool operator==(const NodeIterator& x) const {
      if(graph_ == x.graph_)
        if(current_uid == x.current_uid)
          return true;
      return false;
    }

    /** Defines inequality for node_iterator */
    bool operator!=(const NodeIterator& x) const {
      return !(*this == x);
    }

    /** Defines incrementing for node_iterator */
    NodeIterator& operator++() {
      current_uid++;
      return *this;
    }

   private:
     friend class Graph;
     Graph* graph_;
     size_type current_uid;
  };

  /** Creates beginning iterator for iterating through a graph's nodes */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Creates ending iterator for iterating through a graph's nodes */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
      auto adjacent_node_uids = graph_->nodes_[node1_uid].adjacent_node_uids;
      size_type node2_uid = adjacent_node_uids[iterator_location];
      return Edge(graph_, node1_uid, node2_uid);
    }

    /** Defines incrementing for incident_iterator */
    IncidentIterator& operator++(){
      iterator_location++;
      return *this;
    }

    /** Defines equality for incident_iterator */
    bool operator==(const IncidentIterator& x) const{
      if(graph_ == x.graph_)
        if(node1_uid == x.node1_uid)
          if(iterator_location == x.iterator_location)
            return true;
      return false;
    }

    /** Defines inequality for incident_iterator */
    bool operator!=(const IncidentIterator& x) const{
      return !(*this == x);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_uid;
    size_type iterator_location;

   IncidentIterator(const Graph* graph, size_type node1_uid, size_type iterator_location)
     : graph_(const_cast<Graph*>(graph)),
       node1_uid(node1_uid),
       iterator_location(iterator_location){
   }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
      auto adj_list = graph_->nodes_[node1_uid].adjacent_node_uids;
      Edge edge_to_return = graph_->edge_finder(graph_->node(node1_uid),
                                 graph_->node(adj_list[adj_list_pos]));
      return edge_to_return;
    }

    /** Defines incrementing for edge_iterator */
    EdgeIterator& operator++() {
      //get adjacency list for current node
      auto adj_list = graph_->nodes_[node1_uid].adjacent_node_uids;

      //if there are more incident nodes to visit, go to the next incident node
      if(adj_list_pos < adj_list.size() - 1){
        adj_list_pos++;
      //else, go to next node and visit that node's first incident node
      }else if(adj_list_pos == adj_list.size() - 1){
        node1_uid++;
        adj_list = graph_->nodes_[node1_uid].adjacent_node_uids;
        adj_list_pos = 0;
      }

      //since we dont want to visit edge(a,b) and edge(b,a), skip the edge
      //if the second node has a lower index than the current node we're at
      //since we have already visited it
      while(adj_list[adj_list_pos] < node1_uid){
        //if we are at the last valid edge, increment to match the end iterator
        //and return so we get the equality to the end iterator
        if(node1_uid == graph_->nodes_.size()-1 && adj_list_pos == adj_list.size()-1){
          node1_uid++;
          adj_list_pos++;
          break;
        }else if(adj_list_pos < adj_list.size() - 1){
          adj_list_pos++;
        }else if(adj_list_pos == adj_list.size() - 1){
          node1_uid++;
          adj_list = graph_->nodes_[node1_uid].adjacent_node_uids;
          adj_list_pos = 0;
        }
      }

      return *this;
    }

    /** Defines equality for incident_iterator */
    bool operator==(const EdgeIterator& x) const{
      if(graph_ == x.graph_)
        if(node1_uid == x.node1_uid && adj_list_pos == x.adj_list_pos)
          return true;
      return false;
    }

    /** Defines inequality for incident_iterator */
    bool operator!=(const EdgeIterator& x) const{
      return !(*this == x);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_uid;
    size_type adj_list_pos;

    EdgeIterator(const Graph* graph, size_type node1_uid, size_type adj_list_pos)
      : graph_(const_cast<Graph*>(graph)),
        node1_uid(node1_uid),
        adj_list_pos(adj_list_pos){
    }

  };

  /** Creates beginning iterator to go through all edges of a graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
  /** Creates ending iterator to go through all edges of a graph */
  edge_iterator edge_end() const{
    return EdgeIterator(this,
                        nodes_.size(),
                        nodes_[nodes_.size()-1].adjacent_node_uids.size());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //internal node- like structures
  struct internal_node_element{
    Point location; //The point held by a node
    node_value_type value; //the optional value held by the node
    std::vector<size_type> adjacent_node_uids; //adjacency list of node_uids
  };

  //an array for node-like structures, holds the data of the graph
  std::vector<internal_node_element> nodes_;
};

#endif // CME212_GRAPH_HPP
