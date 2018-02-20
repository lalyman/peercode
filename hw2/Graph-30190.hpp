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
 * The argument of Graph template V is used as the value type of nodes in the graph. 
 * 
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
  
 public:
  /** Type of this graph. */
  using graph_type = Graph;
  
  /** Synonym for V, which is the V argument of the Graph template. */
  using node_value_type = V;
  
  /** Synonym for E, which is the E argument of the Graph template. */
  using edge_value_type = E;
  
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

 private:
  //
  // Predeclaration of internal types and const iterator types.
  //
  struct internal_node;
  struct internal_edge;
  typedef typename std::vector<internal_node>::const_iterator NodeIterator_InternalNodeIterator;
  typedef typename std::vector<size_type>::const_iterator IncidentIterator_InternalEdgelistIterator;
  typedef typename std::vector<internal_edge>::const_iterator EdgeIterator_InternalEdgeIterator;

 public: 
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /* Construct an empty graph. */
  Graph() : nodes_(size_type(0)), edges_(size_type(0)) {
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
    Node() : graph_(NULL), id_(size_type(0)) {
    }
    
    /** Return this node's position. The returning variable are subject to change. */
    Point& position() {
      return graph_->nodes_[id_].point;
    }
    
    /** Return this node's position. The returning type is a const. */
    const Point& position() const {
      return graph_->nodes_[id_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }
    
    /** Return the size of graph that the node is in. */
    size_type graph_size() const {
      return graph_->size();
    }
    
    /** Return this node's value. The returning variable are subject to change. */
    node_value_type& value() {
    	return graph_->nodes_[id_].value;
    }
    
    /** Return this node's value. The returning type is a const. */
    const node_value_type& value() const {
    	return graph_->nodes_[id_].value;
    }
    
    /** Return the degree of this node. */
    size_type degree() const {
      return graph_->nodes_[id_].edgelist.size();
    }
    
    /** Construct the begin iterator for incident iterator. */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, id_, graph_->nodes_[id_].edgelist.begin());
    }
    
    /** Construct the end iterator for incident iterator. */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, id_, graph_->nodes_[id_].edgelist.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (id_ == n.id_);
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
      if (graph_ != n.graph_) {
        return graph_ < n.graph_;
      }
      return (id_ < n.id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    /** Private variables for Node Class. */
    Graph* graph_;
    size_type id_;
    
    /** Private Node constructor. */
    Node(const Graph* graph, size_type id)
        : graph_(const_cast<Graph*>(graph)), id_(id) {
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
   * @param[in] value The new node's value. If not provided, this is filled by default. 
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_node;
    new_node.point = position;
    new_node.value = value;
    nodes_.push_back(new_node);
    return Node(this, size_type(nodes_.size() - 1));
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_) && (n.id_ < this->nodes_.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }
  
  /** Remove the node denoted by @a n. 
   * @param[in] n The Node to be deleted. 
   * @pre G.has_node(n) == true
   * @post for 0 <= @a i < n.index(), new nodes_[i] == old nodes_[i]
   * @post for n.index() < @a i < old nodes_.size(), new nodes_[i] == old nodes_[i]
   * @post new nodes_[n.index] == old nodes_[nodes_.size() - 1]
   * @post new nodes_.size() == old nodes_.size() - 1
   * @post for 0 <= @a i < new nodes_.size(), new nodes_[i].edgelist is consistent 
   *   with edges_. 
   * @return the old n.index()
   *
   * @note This function invalidates all existing Node, Edge, node_iterator, 
   *   edge_iterator, incident_iterator.
   * 
   * Complexity:  
   */
  size_type remove_node(const Node& n) {
    /** Iteratively remove the edges of the node to be removed. */
    auto e_it = n.edge_begin();
    while (e_it != n.edge_end()) {
      remove_edge(*e_it);
      e_it = n.edge_begin();
    }
    
    /** Copy data of the last node to the location of the node to be removed. */
    size_type idr = n.index();
    size_type idl = nodes_.size() - 1;
    nodes_[idr] = nodes_[idl];
    
    /** Iteratively reset the edges in @a edges_ with nodes index @a idl to @a idr. */
    for (auto it = nodes_[idl].edgelist.begin(); it != nodes_[idl].edgelist.end(); ++it) {
      if (edges_[*it].v1 == idl) {
        edges_[*it].v1 = idr;
      } else if (edges_[*it].v2 == idl) {
        edges_[*it].v2 = idr;
      }
    }
    
    /** Remove the last node. */
    nodes_.pop_back();
    return idr;
  }
  
  /** Remove the node pointed by node_iterator @a n_it. 
   * param[in] n_it the node iterator of the node to be deleted
   * return the node iterator pointing to the location of old (*n_it).index()
   * 
   * @pre @post @note and complexity same as "size_type remove_node(const Node& n)"
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type id = remove_node(*n_it);
    return node_iterator(this, this->nodes_.begin() + id);
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
    Edge() : graph_(NULL), id_(size_type(0)), ifflip_(false) {
    }
    
    /** Return a node of this Edge. */
    Node node1() const {
      if (ifflip_ == false) {
        return Node(this->graph_, this->graph_->edges_[this->id_].v1);
      } else {
        return Node(this->graph_, this->graph_->edges_[this->id_].v2);
      }
    }

    /** Return the other node of this Edge. */
    Node node2() const {
      if (ifflip_ == false) {
        return Node(this->graph_, this->graph_->edges_[this->id_].v2);
      } else{ 
        return Node(this->graph_, this->graph_->edges_[this->id_].v1);
      }
    }
    
    /** Return the index of this Edge. */
    size_type index() const {
      return id_;
    }
    
    /** Return the length of this Edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    }
	
	/** Return this edge's value. The returning variable are subject to change. */
    edge_value_type& value() {
    	return graph_->edges_[id_].value;
    }
    
    /** Return this edge's value. The returning type is a const. */
    const edge_value_type& value() const {
    	return graph_->edges_[id_].value;
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) && (id_ == e.id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ != e.graph_) {
        return graph_ < e.graph_;
      }
      return (id_ < e.id_);
    }

   private:
    /** Allow Graph to access Edge's private member data and functions. */
    friend class Graph;
    
    /** Private variables for Edge Class. */
    Graph* graph_;
    size_type id_;
    bool ifflip_;
    
    /** Private constructor for Edge Class. */
    Edge(const Graph* graph, size_type id, bool ifflip = false) 
        : graph_(const_cast<Graph*>(graph)), id_(id), ifflip_(ifflip) {
    }
  };

  /** Return the total number of edges in the graph.
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
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(maximum degree of the graph). 
   */
  bool has_edge(const Node& a, const Node& b) const {
    Node c = a;
    if (a.degree() > b.degree()) {
      c = b;
    }
    for (auto ei = c.edge_begin(); ei != c.edge_end(); ++ei) {
      if ((std::min((*ei).node1().index(), (*ei).node2().index()) 
           == std::min(a.id_, b.id_)) && 
          (std::max((*ei).node1().index(), (*ei).node2().index()) 
           == std::max(a.id_, b.id_))) {
      	return true;
      }
    }
    return false;
  }
  
  /** Find the edge that connects two nodes.
   * @pre @a a and @a b are valid nodes of this graph
   * @return the incident iterator of edgelist of @a, if the edge corresponds to  
   *   @a a and @a b are presented as an edge. Return an invalid incident_iterator 
   *   if it is not an edge. 
   * 
   * Complexity: No more than O(maximum degree of the graph). 
   */
  incident_iterator find_edge(const Node& a, const Node& b) const {
    Node c = a;
    if (a.degree() > b.degree()) {
      c = b;
    }
    for (auto ei = c.edge_begin(); ei != c.edge_end(); ++ei) {
      if ((std::min((*ei).node1().index(), (*ei).node2().index()) 
           == std::min(a.id_, b.id_)) && 
          (std::max((*ei).node1().index(), (*ei).node2().index()) 
           == std::max(a.id_, b.id_))) {
      	return ei;
      }
    }
    return incident_iterator();
  }
  
  /** Remove the edge represented by Edge @a e. 
  * @param[in] e The Edge to be removed. 
  * @pre @a e is a valid edge of this graph. 
  * @post for 0 <= @a i < e.index(), new edge_[i] == old edges_[i]
  * @post for e.index() < @a i < old edges_.size(), new edges_[i] == old edges_[i]
  * @post new edges_[e.index()] == old edges_[edges_.size() - 1]
  * @post new edges_.size() == old edges_.size() - 1
  * @return the old e.index(). 
  * 
  * @note This function will invalidate all Edge, edge_iterator, incident_iterator
  * 
  * Complexity: O(d_max) operations. 
  */
  size_type remove_edge(const Edge& e) {
    size_type id1 = e.node1().index();    // id1 is one node index of the edge to be removed
    size_type id2 = e.node2().index();    // id2 is the other node index of the edge to be removed
    size_type ide = e.index();    // ide is the index of edge to be removed
    size_type le = edges_.size() - 1;    // le is the index of the last edge
	
	/** Remove e.index() from the @a edgelist of the nodes the edge corresponds to. */
    nodes_[id1].edgelist.erase(std::find(nodes_[id1].edgelist.begin(),
        nodes_[id1].edgelist.end(), ide));
    nodes_[id2].edgelist.erase(std::find(nodes_[id2].edgelist.begin(),
        nodes_[id2].edgelist.end(), ide));
    
    size_type idl1 = edges_[le].v1;    // idl1 is one node index of the last edge
    size_type idl2 = edges_[le].v2;    // idl2 is the other node index of the last edge
    
    /** Change the index of the last edge from the @a edgelist of its nodes. */
    *std::find(nodes_[idl1].edgelist.begin(), nodes_[idl1].edgelist.end(), le) = ide;
    *std::find(nodes_[idl2].edgelist.begin(), nodes_[idl2].edgelist.end(), le) = ide;
    
    /** Copy the last edge to the location of the edge to be removed. */
    edges_[ide] = edges_[le];
    /** Delete the last edge. */
    edges_.pop_back();
    return ide;
  }
  
  /** Remove the edge represented by two Nodes. 
  * @param[in] n1, n2 The nodes of the edge to be removed. 
  * @return If the edge represented by n1, n2 is a valid edge of this graph, 
  *    return the old e.index(). If it's not a valid edge, then return 0. 
  *
  * @post, Complexity is the same as "size_type remove_edge(const Edge& e)"
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    incident_iterator pe = find_edge(n1, n2);
    if (!pe.valid()) {
      return size_type(0);
    }
    return remove_edge(*pe);
  }


  /** Remove the edge represented by an edge_iterator. 
  * @param[in] e_it the edge_iterator pointing to the edge to be removed. 
  * @pre (*e_it) is a valid edge of this graph. 
  * @return the edge iterator with index e.index(). 
  *
  * @post, Complexity is the same as "size_type remove_edge(const Edge& e)"
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type ide = remove_edge(*e_it);
    return edge_iterator(this, this->edges_.begin() + ide);
  }


  /** Add an edge to the graph, or return the current edge 
   * if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(maximum degree of the graph). 
   */
  Edge add_edge(const Node& a, const Node& b) {
    if (this->has_edge(a, b) == true) {
	  return Edge(this, 0);
    }
    internal_edge new_edge;
    new_edge.v1 = std::min(a.id_, b.id_);
    new_edge.v2 = std::max(a.id_, b.id_);
    edges_.push_back(new_edge);
    size_type last_edge_id = edges_.size() - 1;
    this->nodes_[a.id_].edgelist.push_back(last_edge_id);
    this->nodes_[b.id_].edgelist.push_back(last_edge_id);
    return Edge(this, last_edge_id);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects. 
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator :  private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), p_() {
    }

    /** Operations for NodeIterator Class. */
    Node operator*() const {
      return Node(graph_, p_ - graph_->node_begin().p_);
    }
    
    NodeIterator& operator++() {
      ++p_; 
      return (*this);
    }

    bool operator==(const NodeIterator& ni) const {
      return (p_ == ni.p_);
    }
    
   private:
    friend class Graph;
    
    /** Private variables for Node Class. */
    NodeIterator_InternalNodeIterator p_;
    Graph* graph_;
    
    /** Private constructor for NodeIterator Class. */
    NodeIterator(const Graph* graph, NodeIterator_InternalNodeIterator p0) 
        : p_(p0), graph_(const_cast<Graph*>(graph)) {
    }
    
  };

  
  /** Construction for the begin iterator for node iterator. */
  node_iterator node_begin() const {
    return(NodeIterator(this, this->nodes_.begin()));
  }

  /** Construction for the end iterator for node iterator. */
  node_iterator node_end() const {
    return(NodeIterator(this, this->nodes_.end()));
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
    IncidentIterator() : graph_(nullptr), node_id_(size_type(0)), p_(){
    }

 	/** Check if IncidentIterator is valid. */
 	bool valid() {
 	  return (graph_ != nullptr);
 	}
 	
    /** Operations for IncidentIterator Class. */    
    Edge operator*() const {
      return Edge(graph_, *p_, graph_->edges_[*p_].v1 != node_id_);
    }
    
    IncidentIterator& operator++() {
      ++p_; 
      return (*this);
    }
    
    bool operator==(const IncidentIterator& ii) const {
      return (p_ == ii.p_);
    }

   private:
    friend class Graph;
    
    /** Private variables for IncidentIterator Class. */
    Graph* graph_;
    size_type node_id_;
    IncidentIterator_InternalEdgelistIterator p_;
    
    /** Private Constructor for IncidentIterator Class. */
    IncidentIterator(const Graph* graph, size_type node_id, 
                     IncidentIterator_InternalEdgelistIterator  p0)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id), p_(p0) {
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
    EdgeIterator() : graph_(nullptr), p_() {
    }

    /** Operations for EdgeIterator Class. */
    Edge operator*() const {
      return Edge(graph_, p_ - graph_->edge_begin().p_);
    }
    
    EdgeIterator& operator++() {
      ++p_;
      return (*this);
    }

    bool operator==(const EdgeIterator& ei) const {
      return (p_ == ei.p_);
    }

   private:
    friend class Graph;
    
    /** Private variables for EdgeIterator Class. */
    Graph* graph_;
    EdgeIterator_InternalEdgeIterator p_;
    
    /** Private constructor for EdgeIterator Class. */
    EdgeIterator(const Graph* graph, EdgeIterator_InternalEdgeIterator p0) 
        : graph_(const_cast<Graph*>(graph)), p_(p0) {
    }
  };

  /** Construction for the begin iterator for edge iterator. */
  edge_iterator edge_begin() const {
    return edge_iterator(this, this->edges_.begin());
  }

  /** Construction for the end iterator for edge iterator. */
  edge_iterator edge_end() const {
    return edge_iterator(this, this->edges_.end());
  }

 private:  
  /** Internal type for node and edge.  */
  struct internal_node {
    Point point;
    node_value_type value;
    std::vector<size_type> edgelist;
  };
  struct internal_edge {
    size_type v1, v2;
    edge_value_type value;
  };
  
  /** Private variables @a nodes_ and @a edges_ for Graph class. */
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
};

#endif // CME212_GRAPH_HPP
