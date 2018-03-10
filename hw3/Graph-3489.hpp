#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
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

  struct internal_node;
  struct internal_edge;

 //private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  typedef V node_value_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  typedef E edge_value_type;

  /** Predeclaration of iterator types. */
  class node_iterator;
  class edge_iterator;
  class incident_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's position. */
    Point& position() {
      // HW0: YOUR CODE HERE
      return fetch().position;
    }

    /** Return this node's unique id. */
    size_type uid() const {
      return uid_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return fetch().index;
    }

    /** Return this node's degree. */
    size_type degree() const {
      return fetch().degree;
    }

    /** Return a reference to this node's value*/
    node_value_type& value() {
      return fetch().value;
    }

    /** Return this node's value*/
    const node_value_type& value() const {
      return fetch().value;
    }

    /** Return whether this node has been removed*/
    bool removed() const {
      if (index() < graph_->num_nodes())
        return graph_->node_ind2uid_[index()] != uid();
      else
        return true;
    }

    /** Return this node's graph pointer*/
    Graph* graph() const {
      return graph_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (uid() == n.uid()) && (graph() == n.graph());
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
      return (index() < n.index());
    }

    /** Return begin incident edge iterator*/
    incident_iterator edge_begin() const {
      incident_iterator it = incident_iterator(this, 0);
      while (it != edge_end())
        if ((*it).removed())
          ++it;
        else break;
      return it;
    }

    /** Return end incident edge iterator*/
    incident_iterator edge_end() const {
      return incident_iterator(this, fetch().incident_edge_uid.size());
    }

   private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;   // Pointer back to the parent Graph
    size_type uid_;  // This node's unique id
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
    /** Helper method to return the appropriate element.*/
    internal_node& fetch() const {
      return graph_->nodes_[uid_];
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class Edge;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_ind2uid_.size();
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
    
    size_type new_uid {(size_type) nodes_.size()};
    size_type new_index {(size_type) num_nodes()};
    nodes_.push_back({position, value, num_nodes()});
    node_ind2uid_.push_back(new_uid);
    return node(new_index);
  }

  /** Remove a node from the graph, returning the new size of the graph.
   * @param[in] n The node to be removed
   * @post new num_nodes == old num_nodes() - 1
   * @post new node(old index) == old node(old index + 1)
   *
   * Complexity: O(n.num_incident_edges).
   */
  size_type remove_node(const Node& n) {
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
      remove_edge(*it);
    
    node_ind2uid_.erase(node_ind2uid_.begin() + n.index());
    for (auto it = node_iterator(this, n.index()); it != node_end(); ++it)
      (*it).fetch().index--;
    return num_nodes();
  }

  /** Remove a node from the graph, returning an iterator pointing to the next
   * non-removed node or returning the end iterator if there is no next node.
   * @param[in] n_it The iterator for the node to be removed
   * @post new num_nodes == old num_nodes() - 1
   * @post new node(old index) == old node(old index + 1)
   *
   * Complexity: O(n.num_incident_edges).
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if ((node(n.index()) == n) && !n.removed())
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
    return Node(this, node_ind2uid_[i]);
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(graph_->nodes_[node1_uid_].index); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(graph_->nodes_[node2_uid_].index); 
    }

    /** Return the distance between the nodes forming this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    size_type uid() const {
      return uid_;
    }
    
    /** Return the edge index */
    size_type index() const {
      return fetch().index;
    }

    /** Return a reference to this edge's value*/
    edge_value_type& value() {
      return fetch().value;
    }

    /** Return a reference to this edge's value*/
    const edge_value_type& value() const {
      return fetch().value;
    }

    /** Return whether this edge has been removed*/
    bool removed() const {
      if (index() < graph_->num_edges())
        return graph_->edge_ind2uid_[index()] != uid();
      else
        return true;
    }

    /** Return this edge's graph pointer*/
    Graph* graph() const {
      return graph_;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.uid() == uid()) && (graph() == e.graph());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (node1() < e.node1()) || (node2() < e.node2());
    }

   private:
    /** Helper method to return the appropriate element.*/
    internal_edge& fetch() const {
      return graph_->edges_[uid()];
    }

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    /** Find the edge uid */
    size_type find_uid() const {
      for (auto it = node1().fetch().incident_edge_uid.begin();
           it != node1().fetch().incident_edge_uid.end(); ++it) {
        size_type test_index {graph_->edges_.at(*it).index};
        size_type test_n1_uid {graph_->edges_.at(*it).node1_uid};
        size_type test_n2_uid {graph_->edges_.at(*it).node2_uid};
        if ((((node1_uid_ == test_n1_uid) && (node2_uid_ == test_n2_uid)) ||
             ((node1_uid_ == test_n2_uid) && (node2_uid_ == test_n1_uid))) &&
            (*it == graph_->edge_ind2uid_[test_index]))
          return *it;
      }

      return -1;
    }

    /** Private constructors */
    Edge(const Node& n1, const Node& n2)
        : graph_(const_cast<Graph*>(n1.graph())), node1_uid_(n1.uid()), node2_uid_(n2.uid()) {
          uid_ = find_uid();
        }
    Edge(const Graph* graph, size_type n1_uid, size_type n2_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(n1_uid), node2_uid_(n2_uid) {
          uid_ = find_uid();
        }
    Edge(const Graph* graph, size_type n1_uid, size_type n2_uid, size_type e_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(n1_uid), node2_uid_(n2_uid), uid_(e_uid) {}

    Graph* graph_;        // Pointer back to the parent Graph
    size_type uid_;       // Edge's unique id
    size_type node1_uid_; // Node 1's unique id
    size_type node2_uid_; // Node 2's unique id

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_ind2uid_.size();
  }

  /** Return the edge with uid @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges_.at(i).node1_uid, edges_.at(i).node2_uid, i);
  }

  /** Return the edge with uid @a i, with the nodes swapped in storage
   * compared to the edge output from edge(size_type i).
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edgerev(size_type i) const {
    return Edge(this, edges_.at(i).node2_uid, edges_.at(i).node1_uid, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(n1.num_incident_edges)
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it)
      if ((*it).node2() == b)
        return true;
      
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
    if (!has_edge(a,b)) {
      size_type new_uid {(size_type) edges_.size()};
      edges_.push_back({a.uid(), b.uid(), num_edges()});
      edge_ind2uid_.push_back(new_uid);
      a.fetch().incident_edge_uid.push_back(new_uid);
      a.fetch().degree++;
      b.fetch().incident_edge_uid.push_back(new_uid);
      b.fetch().degree++;
      return edge(new_uid);
    } else return Edge(a,b);
    
  }

  /** Remove an edge from the graph, returning the new number of graph edges.
   * @param[in] e The edge to be removed
   * @post new num_edges() == old num_edges() - 1
   * @post new edge(old index) == old edge(old index + 1)
   *
   * Complexity: O(1).
   */
  size_type remove_edge(const Edge& e) {
    if (!e.removed()) {
      e.node1().fetch().degree--;
      e.node2().fetch().degree--;
      edge_ind2uid_.erase(edge_ind2uid_.begin() + e.index());
      for (auto it = edge_iterator(this, e.index()); it != edge_end(); ++it)
        (*it).fetch().index--;
    }
    return num_edges();
  }

  /** Remove an edge from the graph, returning the new number of graph edges.
   * @param[in] n1 One node in the edge to be removed
   * @param[in] n2 The other node in the edge to be removed
   * @post new num_edges() == old num_edges() - 1
   * @post new edge(old index) == old edge(old index + 1)
   *
   * Complexity: O(n1.num_incident_edges).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    for (auto it = n1.edge_begin(); it != n1.edge_end(); ++it)
      if ((*it).node2() == n2) {
        remove_edge(*it);
      }
    return num_edges();
  }

  /** Remove an edge from the graph, returning an iterator pointing to the next
   * non-removed edge or returning the end iterator if there is no next edge.
   * @param[in] n1 One node in the edge to be removed
   * @post new num_edges() == old num_edges() - 1
   * @post new edge(old index) == old edge(old index + 1)
   *
   * Complexity: O(1).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    node_ind2uid_.clear();
    edge_ind2uid_.clear();
  }

  /** @class Graph::node_iterator
   * @brief Class for iterating over the graph's nodes.
  */
  class node_iterator : private equality_comparable<node_iterator> {

   public:
    node_iterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {}

    /** Return node pointed to by the iterator. */
    Node operator*() const {
      return graph_->node(index_);
    }

    /** Return iterator pointing to the next node in the graph. */
    node_iterator& operator++() {
      index_++;
      return *this;
    }

    /** Return iterator pointing to the previous node in the graph. */
    node_iterator& operator--() {
      index_--;
      return *this;
    }

    /** Check whether two node_iterators are equivalent. */
    bool operator==(const node_iterator& iter) const {
      return (iter.graph_ == graph_) && (iter.index_ == index_);
    }

   private:

    Graph* graph_;     // Pointer back to the parent Graph
    size_type index_;  // this iterator's current node index, a number in the range [0, graph_size)

  };

  /** Return begin node iterator*/
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /** Return end node iterator*/
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }

  /** @class Graph::incident_iterator
   * @brief Class for iterating over the node's incident edges.
  */
  class incident_iterator : private equality_comparable<incident_iterator> {
   public:
    incident_iterator(const Node* node, size_type index)
        : node_(const_cast<Node*>(node)), index_(index) {}

    /** Return incident edge pointed to by the iterator. */
    Edge operator*() const {
      if (node_->graph()->edge(edge_uid()).node1() == *node_)
        return node_->graph()->edge(edge_uid());
      else
        return node_->graph()->edgerev(edge_uid());
    }

    /** Return iterator pointing to the next incident edge in the graph. */
    incident_iterator& operator++() {
      index_++;

      while (*this != node_->edge_end())
        if (node_->graph()->edge(edge_uid()).removed())
          index_++;
        else break;

      return *this;
    }

    /** Return uid of the edge pointed to by the iterator. */
    size_type edge_uid() const {
      return node_->fetch().incident_edge_uid[index_];
    }

    /** Check whether two incident_iterators are equivalent. */
    bool operator==(const incident_iterator& iit) const {
      return (*(iit.node_) == *node_) && (iit.index_ == index_);
    }

   private:

    Node* node_;     // Pointer back to the parent Node
    size_type index_;  // this iterator's current edge index, a number in the range [0, graph_size)

  };

  /** @class Graph::edge_iterator
   * @brief Class for iterating over the graph's edges.
  */
  class edge_iterator : private equality_comparable<edge_iterator> {
   public:
    edge_iterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {}

    /** Return edge pointed to by the iterator. */
    Edge operator*() const {
      return graph_->edge(graph_->edge_ind2uid_[index_]);
    }

    /** Return iterator pointing to the next edge in the graph. */
    edge_iterator& operator++() {
      index_++;
      return *this;
    }

    /** Check whether two edge_iterators are equivalent. */
    bool operator==(const edge_iterator& iter) const {
      return (iter.graph_ == graph_) && (iter.index_ == index_);
    }
   private:

    Graph* graph_;     // Pointer back to the parent Graph
    size_type index_;  // this iterator's current node index, a number in the range [0, graph_size)

  };

  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }

 private:
  // Internal type for graph nodes
  struct internal_node {
    Point position;         // The position of a node
    node_value_type value;  // The value of a node
    size_type index;        // The index of a node
    size_type degree;       // The degree of a node
    std::vector<size_type> incident_edge_uid;   // List of incident edge uids
    internal_node() : degree(0) {};
    internal_node(Point p, node_value_type v, size_type i)
    : position(p), value(v), index(i), degree(0) {}
  };
  // Internal type for graph edges
  struct internal_edge {
    size_type node1_uid;    // The unique id of one connected node
    size_type node2_uid;    // The unique id of the other connected node
    edge_value_type value;  // The value of an edge
    size_type index;        // The index of an edge
    internal_edge() {};
    internal_edge(size_type n1_uid, size_type n2_uid, size_type i)
    : node1_uid(n1_uid), node2_uid(n2_uid), index(i) {}
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<size_type> node_ind2uid_;
  std::vector<size_type> edge_ind2uid_;

};

#endif // CME212_GRAPH_HPP
