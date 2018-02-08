#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief A directed graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D directed graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are directed (there exist 
 * two edges for each unique pair of nodes, sinsitive to direction).
 */
template<typename V> //Make Graph a template
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

  /** Type of the template parameter. */
  using node_value_type = V;

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
      return graph_->nodes_edge_map_[index_].size();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_->has_node(n) && (n.index_ == index_);
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
      return graph_->has_node(n) && (n.index_ < index_);
    }

    /** Return the Incident Iterators for the first and last incident nodes
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,index_,graph_->nodes_edge_map_[index_].begin());
    }
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,index_,graph_->nodes_edge_map_[index_].end());
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
    // Add the new node object
    node_objs_.push_back(Node(this, num_nodes()));
    // Add the new internal node to graph
    nodes_.push_back({position,value});
    // Return a node that associates with the new internal one
    return node_objs_.back(); 
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
    return node_objs_[i];
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
      return graph_->node_objs_[index_];
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
      return Edge(graph_,iter_->second);
    }
    
    /** Test if this Incident Iterator and @a iiter are equal. */
    bool operator==(const IncidentIterator& iiter) const {
      return (node_ == iiter.node_) && (iter_ == iiter.iter_);
    }
    
    /** Advance the Incident Iterator. */
    IncidentIterator& operator++() {
      ++iter_;
      return *this;
    }

   private:
    // Allow Graph to access IncidentIterator's private member data and functions.
    friend class Graph;

    // Private data members for Node Iterator
    Graph* graph_;
    size_type node_;
    std::map<size_type, size_type>::iterator iter_;
    
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, size_type node,
                     std::map<size_type, size_type>::iterator iter)
        : graph_(const_cast<graph_type*>(graph)), node_(node), iter_(iter) {
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
      return graph_->edge_objs_[index_];
    }
    
    /** Test if this Edge Iterator and @a eiter are equal. */
    bool operator==(const EdgeIterator& eiter) const {
      return (graph_ == eiter.graph_) && (index_ == eiter.index_);
    }
    
    /** Advance the Edge Iterator. */
    EdgeIterator& operator++() {
      index_ += 2;
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
      return graph_->node_objs_[fetch().node1];
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      return graph_->node_objs_[fetch().node2];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same directed edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      return node1() == e.node1() && node2() == e.node2();
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      return node1() < e.node1();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type index_;

    // Private constructor
    Edge(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index)  {
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
    return edge_objs_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    if (nodes_edge_map_.count(a.index_) == 1)
      if (nodes_edge_map_.at(a.index_).count(b.index_) == 1)
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
  edge_type add_edge(const node_type& a, const node_type& b) {
    assert(this->has_node(a) && this->has_node(b));
    // Check if such edge already exists
    size_type s = std::min(a.index_, b.index_);
    size_type l = std::max(a.index_, b.index_);
    if (has_edge(a, b))
      return edge_objs_[nodes_edge_map_[a.index_][b.index_]];
    // Add the new edge with smaller-first-ordering to nodes_edge_map_, add the
    // new edge object, and add the new internal edge to graph
    nodes_edge_map_[s][l] = num_edges();
    edge_objs_.push_back(Edge(this, num_edges()));
    edges_.push_back({s, l});
    // Add the new edge with larger-first-ordering to nodes_edge_map_, add the
    // new edge object, and add the new internal edge to graph
    nodes_edge_map_[l][s] = num_edges();
    edge_objs_.push_back(Edge(this, num_edges()));
    edges_.push_back({l, s});

    // Return an edge that associates with the new internal one
    return edge_objs_[nodes_edge_map_[a.index_][b.index_]];
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

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    Point position;
    node_value_type value;
  };

  struct internal_edge {
    const size_type node1;
    const size_type node2;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<Node> node_objs_;
  std::vector<Edge> edge_objs_;
  std::map<size_type, std::map<size_type, size_type>> nodes_edge_map_;
};

#endif // CME212_GRAPH_HPP
