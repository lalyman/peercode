#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are undirected (there exist 
 * two edges for each unique pair of nodes, insensitive to direction).
 */
template<typename V, typename E> //Make Graph a template
class Graph {
 private:
  // Data structures that store information for the Node and Edge class
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of the template parameters. */
  using node_value_type = V;
  using edge_value_type = E;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of NodeIterator type. */
  class NodeIterator;
  /** Synonym for NodeIterator (following STL conventions). */
  using node_iterator = NodeIterator;

  /** Predeclaration of IncidentIterator type. */
  class IncidentIterator;
  /** Synonym for IncidentIterator (following STL conventions). */
  using incident_iterator = IncidentIterator;

  /** Predeclaration of EdgeIterator type. */
  class EdgeIterator;
  /** Synonym for EdgeIterator (following STL conventions). */
  using edge_iterator = EdgeIterator;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
   * Return type of Graph::Node::index(), Graph::num_nodes(),
   * Graph::num_edges(), and argument type of Graph::node(size_type)
   */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
  }

  /** Default destructor */
  ~Graph() {
  }

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

    /** Return this node's position. */
    Point& position() {
      return fetch().position;
    }
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's value. */
    node_value_type& value() {
      return fetch().value;
    }
    const node_value_type& value() const {
      return fetch().value;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this node's degree. */
    size_type degree() const {
      return fetch().incident_edges.size();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph_ == graph_ && (n.index_ == index_);
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
      if (n.graph_ == graph_)
        return n.index_ < index_;
      else
        return n.graph_ < graph_;
    }

    /** Return the Incident Iterators for the first and last incident nodes
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,index_,0);
    }
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,index_,degree());
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date Position information
    // @pre this node is valid (0 <= index_ < graph_->num_nodes())
    internal_node& fetch() const {
      assert(index_ < graph_->num_nodes());
      return graph_->nodes_[index_];
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* graph_;
    size_type index_;

    /** Private Constructor */
    Node(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
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
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes() 
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node(const Point& position,
                     const node_value_type& value = node_value_type()) {
    // Add the new internal node to graph
    nodes_.push_back({position,value,{}});
    // Return a node that associates with the new internal one
    return Node(this, num_nodes()-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    assert(0 <= i && i < num_nodes());
    return Node(this, i);
  }

  /** Remove node @a n from the graph.
   * @param[in] n The node to remove from the graph
   * @post If @a n is in this graph, new num_nodes() = old num_nodes() - 1
   *                                 new num_edges() = old num_edges - n.degree()
   * @return The number of elements removed from the graph (i.e. 0 if @a n 
   *         doesn't exist, and 1 if @a n does exist)
   *
   * Complexity: O(n.degree()^2)
   */
  size_type remove_node(const node_type& n) {
    // Return no removed nodes if this node doesn't exist in graph
    if (!has_node(n))
      return 0;

    size_type i = n.index_;
    // Delete all incident edges from the list of edges
    while (n.degree() > 0)
      remove_edge(edge(nodes_[i].incident_edges.back()));
    
    // Delete this node from the list of nodes
    nodes_[i] = nodes_.back();
    nodes_.pop_back();

    // Update edge information for the swapped node
    if (i < num_nodes())
      update_edge(i);
    return 1;
  }
  node_iterator remove_node(node_iterator niter) {
    auto n = *niter;
    remove_node(n);
    return niter;
  }


  //
  // NODE ITERATOR
  //

  /** @class Graph::NodeIterator
   * @brief Class representing the Node Iterator abstraction.
   *
   * The Node Iterator abstraction is used to efficiently loop over the Nodes
   * in the Graph.
   */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // Intialize all of the iterator traits
    using value_type        = Node;
    using pointer           = Node*;
    using reference         = Node&;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    /** Construct an empty NodeIterator. */
    NodeIterator() {
    }

    /** Return Node with this Node Iterator @a index. */
    Node operator*() const {
      return Node(graph_,index_);
    }
    
    /** Test if this Node Iterator and @a niter are equal. */
    bool operator==(const NodeIterator& niter) const {
      return (graph_ == niter.graph_) && (index_ == niter.index_);
    }
    
    /** Advance the Node Iterator. */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

   private:
    // Allow Graph to access NodeIterator's private member data and functions.
    friend class Graph;

    // Private data members for Node Iterator
    Graph* graph_;
    size_type index_;
    
    /** Private Constructor */
    NodeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
    }
  };

  /** Return the Node Iterators for the first and last nodes
   *
   * Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(this,0);
  }
  node_iterator node_end() const {
    return NodeIterator(this,num_nodes());
  }
   

  //
  // INCIDENT ITERATOR
  //

  /** @class Graph::IncidentIterator
   * @brief Class representing the Incident Iterator abstraction.
   *
   * The Incident Iterator abstraction is used to efficiently loop over the 
   * Edges incident to a Node.
   */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // Intialize all of the iterator traits
    using value_type = Edge;
    using pointer    = Edge*;
    using reference  = Edge&;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    /** Construct an empty IncidentIterator. */
    IncidentIterator() {
    }

    /** Return Edge with this Incident Iterator @a index. */
    Edge operator*() const {
      size_type edge_id = graph_->nodes_[node_].incident_edges[index_];
      size_type node_1 = node_;
      size_type node_2 = graph_->edges_[edge_id].node_2;
      if (node_1 == node_2)
        node_2 = graph_->edges_[edge_id].node_1;
      return Edge(graph_,node_1,node_2,edge_id);
    }
    
    /** Test if this Incident Iterator and @a iiter are equal. */
    bool operator==(const IncidentIterator& iiter) const {
      return (node_ == iiter.node_) && (index_ == iiter.index_);
    }

    /** Return local index. */
    size_type index() const {
      return index_;
    }
    
    /** Advance the Incident Iterator. */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }

   private:
    // Allow Graph to access IncidentIterator's private member data and functions.
    friend class Graph;

    // Private data members for Node Iterator
    Graph* graph_;
    size_type node_;
    size_type index_;
    
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, size_type node, size_type index)
        : graph_(const_cast<graph_type*>(graph)), node_(node), index_(index) {
    }
  };
   

  //
  // EDGE ITERATOR
  //

  /** @class Graph::EdgeIterator
   * @brief Class representing the Edge Iterator abstraction.
   *
   * The Edges Iterator abstraction is used to efficiently loop over all of the
   * Edges in the graph.
   */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // Intialize all of the iterator traits
    using value_type = Edge;
    using pointer    = Edge*;
    using reference  = Edge&;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    /** Construct an empty EdgeIterator. */
    EdgeIterator() {
    }

    /** Return Edge with this Edge Iterator @a index. */
    Edge operator*() const {
      return Edge(graph_,graph_->edges_[index_].node_1,
                  graph_->edges_[index_].node_2,index_);
    }
    
    /** Test if this Edge Iterator and @a eiter are equal. */
    bool operator==(const EdgeIterator& eiter) const {
      return (graph_ == eiter.graph_) && (index_ == eiter.index_);
    }
    
    /** Advance the Edge Iterator. */
    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }

   private:
    // Allow Graph to access EdgeIterator's private member data and functions.
    friend class Graph;

    // Private data members for Edge Iterator
    Graph* graph_;
    size_type index_;
    
    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
    }
  };

  /** Return the Edge Iterators for the first and last nodes
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }


  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-sensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    node_type node1() const {
      return Node(graph_,node_1_);
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      return Node(graph_,node_2_);
    }

    /** Return this edge's value. */
    edge_value_type& value() {
      return fetch().value;
    }
    const edge_value_type& value() const {
      return fetch().value;
    }

    /** Return this edge's length. */
    double length() const{
      return norm_2(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      return (graph_ == e.graph_ && index_ == e.index_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      if (e.graph_ == graph_)
        return e.index_ < index_;
      else
        return e.graph_ < graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type node_1_;
    size_type node_2_;
    size_type index_;

    // Private constructor
    Edge(const graph_type* graph, size_type node_1, 
         size_type node_2, size_type index)
        : graph_(const_cast<graph_type*>(graph)), node_1_(node_1), 
                 node_2_(node_2), index_(index)  {
    }

    // Communicate with the associated Graph object 
    // to obtain up-to-date information
    // @pre this node is valid (0 <= index_ < graph_->num_edges())
    internal_edge& fetch() const {
      assert(index_ < graph_->num_edges());
      return graph_->edges_[index_];
    }

  };

  /** Return the total number of directed edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(0 <= i && i < num_edges());
    return Edge(this, edges_[i].node_1, edges_[i].node_2, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a n1 and @a n2 are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a n1 and @a n2.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& n1, const node_type& n2) const {
    for (auto iiter = n1.edge_begin(); iiter != n1.edge_end(); ++iiter) {
      auto e = *iiter;
      if ((e.node1() == n1 && e.node2() == n2) || 
          (e.node1() == n2 && e.node2() == n1))
        return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a n1 and @a n2 are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a n1 and e.node2() == @a n2
   * @post has_edge(@a n1, @a n2) == true
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges().
   *       Else,                          new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type add_edge(const node_type& n1, const node_type& n2, 
                     const edge_value_type& value = edge_value_type()) {
    assert(this->has_node(n1) && this->has_node(n2));
    // Check if such edge already exists
    if (has_edge(n1, n2))
      return Edge(this, n1.index(), n2.index(), -1);
    // Add the new internal edge to the Graph
    edges_.push_back({n1.index(), n2.index(), value});
    // Add incident edge information to the internal nodes
    nodes_[n1.index()].incident_edges.push_back(num_edges()-1);
    nodes_[n2.index()].incident_edges.push_back(num_edges()-1);
    // Return an edge that associates with the new internal one
    return Edge(this, n1.index(), n2.index(), num_edges()-1);
  }

  /** Remove edge between nodes @a n1 and @a n2 from the graph.
   * @param[in] n1,n2 The nodes defining an edge in the graph
   * @param[in] e     The edge defined by nodes @a n1 and @a n2
   * @post If @a e is in this graph, new num_edges() = old num_edges() - 1
   * @return The number of elements removed from the graph (i.e. 0 if @a e 
   *         doesn't exist, and 1 if @a e does exist)
   *
   * Complexity: O(n1.degree() + n2.degree())
   */
  size_type remove_edge(const node_type& n1, const node_type& n2) {
    // Return no removed edges if this edge doesn't exist in graph
    if (!has_edge(n1,n2))
      return 0;

    // Delete connectivity info in node attributes and find edge index
    size_type i = 0;
    for (auto iiter = n1.edge_begin(); iiter != n1.edge_end(); ++iiter) {
      auto e = *iiter;
      // Delete this node from the list of nodes
      if (e.node2() == n2) {
        nodes_[n1.index()].incident_edges.erase(
          nodes_[n1.index()].incident_edges.begin() + iiter.index());
        break;
      }
    } 
    for (auto iiter = n2.edge_begin(); iiter != n2.edge_end(); ++iiter) {
      auto e = *iiter;
      // Delete this node from the list of nodes
      if (e.node2() == n1) {
        i = nodes_[n2.index()].incident_edges[iiter.index()];
        nodes_[n2.index()].incident_edges.erase(
          nodes_[n2.index()].incident_edges.begin() + iiter.index());
        break;
      }
    } 

    // Delete this edge from the list of edges
    edges_[i] = edges_.back();
    edges_.pop_back();

    // Update node information for the swapped edge
    if (i < num_edges()) {
      update_node(i,edges_[i].node_1);
      update_node(i,edges_[i].node_2);
    }
    return 1;
  }
  size_type remove_edge(const edge_type& e) {
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1,n2);
  }
  edge_iterator remove_edge(edge_iterator eiter) {
    auto e = *eiter;
    remove_edge(e);
    return eiter;
  }


  //
  // DESTRUCTORS
  //

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    Point position;
    node_value_type value;
    std::vector<size_type> incident_edges;
  };

  struct internal_edge {
    size_type node_1;
    size_type node_2;
    edge_value_type value;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

  /** Method to update node information for an edge with a modified index. */
  void update_node(const size_type edge_id, const size_type node_id) {
    for (auto iiter = node(node_id).edge_begin(); iiter != node(node_id).edge_end(); ++iiter) {
      if (nodes_[node_id].incident_edges[iiter.index()] == num_edges()) {
        nodes_[node_id].incident_edges[iiter.index()] = edge_id;
        break;
      }
    }
  }

  /** Method to update edge information for a node with a modified index. */
  void update_edge(const size_type node_id) {
    for (auto iiter = node(node_id).edge_begin(); iiter != node(node_id).edge_end(); ++iiter) {
      size_type edge_id = nodes_[node_id].incident_edges[iiter.index()];
      if (edges_[edge_id].node_1 == num_nodes())
        edges_[edge_id].node_1 = node_id;
      else if (edges_[edge_id].node_2 == num_nodes())
        edges_[edge_id].node_2 = node_id;
    }
  }

};

#endif // CME212_GRAPH_HPP
