#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <iostream>
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
  using graph_type = Graph<V>;

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
  Graph():num_nodes_(0),num_edges_(0){}

  /** Default destructor */
  ~Graph() = default;

  using node_value_type = V;

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
    Node(): uid_(0), graph_(nullptr) , value_() {}
    Node(const graph_type* graph, size_type uid, node_value_type value = node_value_type()):
            graph_(const_cast<graph_type *> (graph)), uid_(uid), value_(value){}



    /** Return this node's position.
     * @pre 0 < this->index() < graph node number
     * @return node's position
     */
    const Point& position() const {
      if(this->index() < graph_->num_nodes())
        return graph_->position_[this->index()];
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * @return node id
     */
    size_type index() const {
      return size_type(uid_);
    }

    /** Return this node's graph pointer.
     * @return graph pointer
     */
    graph_type* graph() const {
      return graph_;
    }


    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value by reference.
     *  @return node value
     */
    node_value_type & value () {return value_;};

    /** Return this node's value by constant reference.
     * @return node value
     */
    const node_value_type & value () const {return value_;};


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->graph_ ==  n.graph() && this->uid_ == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * @param[in] Node
     * @return true if two nodes are in the same graph and uid_ < n.index
     */
    bool operator<(const Node& n) const {
       if(this->graph_ !=  n.graph()){
           std::cout <<"Nodes are in different graphs, cannot compare." << std::endl;
           return false;
       }
       return (this->uid_ < n.index());
    }

   /**
    * @return the number of incident edges of the current node
    */
    size_type degree() const {graph_->incident_edges_[uid_].size();};

    /**
    * @return the begin iterator of the node's incident iterator
    */
    incident_iterator edge_begin() const {return IncidentIterator(&graph_->incident_edges_[uid_][0]);};

    /**
    * @return the end iterator of the node's incident iterator
    */
    incident_iterator edge_end() const {return IncidentIterator(&graph_->incident_edges_[uid_][0] + graph_->incident_edges_[uid_].size());};

   private:
    // Allow Graph to access Node's private member data and functions.

      graph_type * graph_;
      size_type uid_;
      node_value_type value_;


    friend class Graph;

  };

  /** @return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value_type The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node ( const Point & position, const node_value_type & value = node_value_type()) {
    num_nodes_ += 1;
    position_.push_back(position);


    Node new_node = Node(this, num_nodes_ - 1, value);
    nodes_.push_back(new_node);

    incident_edges_.push_back( std::vector<Edge>());

    return new_node;
  }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i >= num_edges()) {
      std::cout << "In func node: @pre 0 <= @a i < num_nodes()" <<std::endl;
      return Node();
      }
    else
      return this->nodes_[i];
  }

  /** Return the node with index @a i.
  * @pre 0 <= @a i < num_nodes()
  * @post result_node.index() == i
  *
  * Complexity: O(1).
  */
  Node & node(size_type i) {
    if(i >= num_edges()) {
      std::cout << "In func node: @pre 0 <= @a i < num_nodes()" <<std::endl;
      return Node();
    }
    else
      return this->nodes_[i];
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge(): graph_(nullptr), uid_(0), n1_(0), n2_(0){
    }
    /** Construct a valid Edge. */
    Edge(const graph_type* graph, size_type n1, size_type n2, size_type uid):
            graph_(const_cast<graph_type *>(graph)),uid_(uid), n1_(n1), n2_(n2){}

    /** Return a node of this Edge */
    Node& node1() const {
      return graph_->nodes_[n1_];
    }

    /** Return the other node of this Edge */
    Node& node2() const {
        return graph_->nodes_[n2_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * @return ture if edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.node1()==this->node1() && e.node2()==this->node2())
              || (e.node1()==this->node2() && e.node2()==this->node1()));
    }



    /** @return this edge's index, a number in the range [0, graph_edge_size). */
    size_type index() const {
      return size_type(uid_);
    }


   /** Test whether this edge is less than @a e in a global order.
    *
    * This ordering function is useful for STL containers such as
    * std::map<>. It need not have any interpretive meaning.
    * @return true if the index of the edge <  e.index
    */
    bool operator<(const Edge& e) const {
        return (this->index() < e.index());
    }

   private:
      graph_type * graph_;
      size_type uid_;
      size_type n1_, n2_;
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

  };

  /** @return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** @return  the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i >= num_edges()) {
      std::cout << "In func edge: @pre 0 <= @a i < num_edges()" <<std::endl;
      return Edge();
    }
    else{
      return edges_[i];
    }

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      for(auto it = a.edge_begin(); it != a.edge_end(); ++it)
          if((*it).node2() == b)   return true;
      for(auto it = b.edge_begin(); it != b.edge_end(); ++it)
          if((*it).node2() == a)   return true;
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


      for(auto it = a.edge_begin(); it != a.edge_end(); ++it)
          if((*it).node2() == b)   return *it;
      for(auto it = b.edge_begin(); it != b.edge_end(); ++it)
          if((*it).node2() == a)   return *it;


      num_edges_ += 1;
      Edge new_edge = Edge(this, a.index(), b.index(), num_edges_ - 1);
      edges_.push_back(new_edge);

      incident_edges_[a.index()].push_back(new_edge);
      incident_edges_[b.index()].push_back(Edge(this, b.index(), a.index(), num_edges_ - 1));

      return new_edge;

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_edges_ = 0;
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Construct a NodeIterator.
     *  @param[in] p: Node pointer
     */
    NodeIterator(const Node* p): p_(const_cast<Node*>(p)){}
    /** Dereference a NodeIterator.
    *   @return the constant reference of the node
    */
    const Node &operator *() const { return *p_ ; }
   /** Dereference a NodeIterator.
    *   @return the reference of the node
    */
    Node &operator *() { return *p_ ; }
   /** Compare two iterators.
    *  @param[in] constant reference of a Nodeiterator
    *  @return true, if these two iterator are the same
    */
    bool operator ==( const NodeIterator & x ) const { return p_ == x.p_ ; }
   /** Self increment.
    *  @return next iterator
    */
    NodeIterator& operator ++(){ p_++;  return *this; }


   private:
    friend class Graph;
    Node* p_ ;

  };

  // Supply definitions AND SPECIFICATIONS for:
 /**
  *  @return the head of the node iterator
  */
  node_iterator node_begin() const {return NodeIterator(&nodes_[0]);}
  /**
   *  @return the tail of the node iterator
   */
  node_iterator node_end() const {return NodeIterator(&nodes_[0]+num_nodes_);}

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

   /** Construct an IncidentIterator.
    *  @param[in] p: Edge pointer
    */
    IncidentIterator(const Edge* p): p_(const_cast<Edge*>(p)){}

   /** Dereference.
    *   @return the constant reference of the edge
    */
    Edge &operator*() const{ return *p_;}
   /** Dereference a IncidentIterator.
    *   @return the reference of the node
    */
    IncidentIterator& operator++(){ p_++; return *this;}


   /** Compare two iterators.
    *  @param[in] constant reference of a Incidentiterator
    *  @return true, if these two iterator are the same
    */
    bool operator==(const IncidentIterator& x) const{ return p_ == x.p_;}

   private:
    friend class Graph;
    Edge* p_ ;

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}
   /** Construct an EdgeIterator.
    *  @param[in] p: Edge pointer
    */
    EdgeIterator(const Edge* p): p_(const_cast<Edge*>(p)){}

   /** Dereference an EdgeIterator.
    *   @return the constant reference of the edge
    */

    Edge &operator*() const{ return *p_;}
   /** Self increment.
    *  @return next iterator
    */
    EdgeIterator& operator++(){p_++; return *this;}

   /** Compare two iterators.
    *  @param[in] constant reference of a Edgeiterator
    *  @return true, if these two iterator are the same
    */
    bool operator==(const EdgeIterator& x) const{return p_ == x.p_;}

    private:
    friend class Graph;
    Edge* p_ ;

  };


 /**
  *  @return the head of the edge iterator
  */
  edge_iterator edge_begin() const{ return EdgeIterator(&edges_[0]);}
 /**
  *  @return the head of the edge iterator
  */
  edge_iterator edge_end() const{ return EdgeIterator(&edges_[0] + num_edges_);}

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<Node> nodes_;     //Node vector
  std::vector<Point> position_; //Node position vector
  std::vector<Edge> edges_;     //Edge vector
  std::vector<std::vector<Edge>> incident_edges_;  //Node incident edge vector
  size_type num_nodes_;
  size_type num_edges_;


};

#endif // CME212_GRAPH_HPP
