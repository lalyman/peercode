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
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

// HW1: YOUR CODE HERE
template <typename V>
class Graph {
    
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
    struct internal_nodes; // Structure is introduced to save 'value'
    struct internal_edges;

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

  // HW1: YOUR CODE HERE
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), next_uid_node_(0), edges_(), next_uid_edge_(0) {
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Set this node's value. */
    node_value_type& value() {
      return const_cast<V&>(graph_->nodes_[uid_].value);
    }
    
    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value;
    }
    
    /** Return number of edges incident to this node. */
    size_type degree() const {
        return graph_->adj_[uid_].size();
    }
    
    /** Return iterator that begins iteration of Incident Iteator. */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_,uid_,0);
    }
    
    /** Return iterator that ends iteration of Incident Iteator. */
    incident_iterator edge_end() const {
        return IncidentIterator(graph_,uid_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if ( !(graph_==n.graph_) ) std::cout << "Check Node== : same graph?" << std::endl;
      if ( uid_==n.index() && graph_==n.graph_ )
          return true;
      else
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
      // HW0: YOUR CODE HERE
      (void) n;           // Quiet compiler warning
      if ( !(graph_==n.graph_) ) std::cout << "Check Node< : same graph?" << std::endl;
      if ( uid_ < n.index() )
          return true;
      else
          return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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

  /** No longer needs this add_node(const Point& position)
   *
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    nodes_[next_uid_node_].position = position;
    nodes_[next_uid_node_].value = node_value_type();//check this part
    next_uid_node_++;
    return Node(this,next_uid_node_-1); 
  }
  */

  // HW1: YOUR CODE HERE
  /** Add a node to the graph, returning the added node. **/
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    (void) position, (void) value; // Quiet compiler warning
    local_node_.position = position;
    local_node_.value= value;
    nodes_.push_back(local_node_);
    adj_.push_back(std::vector<size_type>()); //empty vector to store incident edges later 
    next_uid_node_++;
    return Node(this,next_uid_node_-1); 
  }

   
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    if ( !(this==n.graph_) ) std::cout << "Check has_node : same graph?" << std::endl;
    for ( size_type i = 0; i < size(); ++i )
    {
       if ( node(i) == n )
           return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    if (i < num_nodes())
        return Node(this,i);
    else
    {
        std::cout << "Check node : out of node size" << std::endl;
        return Node();    //Invalid node
    }
  }
    // HW1: YOUR CODE HERE
  /** Return the incident node of node @a uid_node with index @a uid_adj. **/
  Node node(size_type uid_node, size_type uid_adj) const {
    (void) uid_node, (void) uid_adj;
    if ( uid_node < num_nodes() && uid_adj < node(uid_node).degree() )
        return Node(this,adj_[uid_node][uid_adj]);
    else
    {
        std::cout << "Check node : out of node size" << std::endl;
        return Node();    //Invalid node
    }
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
    Edge() {
      // HW0: YOUR CODE HERE
    }
    
    // HW0: YOUR CODE HERE
    // This function is added to get axess to uid_ value for edge
    size_type index() const {
      return uid_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(graph_->edges_[uid_].uid1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(graph_->edges_[uid_].uid2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if ( !(graph_==e.graph_) ) std::cout << "Check Edge== : same graph?" << std::endl;
      if ( this->node1()==e.node1() && this->node2()==e.node2() )
          return true;
      if ( this->node1()==e.node2() && this->node2()==e.node1() )
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      (void) e;           // Quiet compiler warning
      if ( !(graph_==e.graph_) ) std::cout << "Check Edge< : same graph?" << std::endl;
      if ( uid_ < e.index() )
          return true;
      else
          return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    Edge(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    if ( i < num_edges() )
        return Edge(this,i);
    else
        std::cout << "Check edge : out of edge size" << std::endl;
        return Edge();        // Invalid Edge
  }
    // HW1: YOUR CODE HERE
  /** Return the edge contains node @a a and node @a b. **/
  Edge edge(const node_type& a, const node_type& b) const {
    (void) a; (void) b;   // Quiet compiler warning
    for(size_type i = 0; i < num_edges(); ++i) {
        if ( edge(i).node1()==a && edge(i).node2()==b )
           return Edge(this,i);
        if ( edge(i).node1()==b && edge(i).node2()==a )
           return Edge(this,i);
    }
    std::cout << "Check edge : no valid node" << std::endl;
    return Edge();        // Invalid Edge
  }
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    if ( !(this==a.graph_ && this==a.graph_) ) std::cout << "Check has_edge : same graph?" << std::endl;
    if (a==b) std::cout << "Check has_edge : same node!" << std::endl;
    for (size_type i = 0; i < num_edges(); ++i)
    {
        if ( edge(i).node1()==a && edge(i).node2()==b )
           return true;
        if ( edge(i).node1()==b && edge(i).node2()==a )
           return true;
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    (void) a, (void) b;   // Quiet compiler warning
    if ( has_edge(a,b) == false )
    {
        local_edge_.uid1 = a.index();   
        local_edge_.uid2 = b.index();
        //edges_[next_uid_edge_] = local_edge_;
        edges_.push_back(local_edge_);
        //HW1 pushing adj
        adj_[a.index()].push_back(b.index());
        adj_[b.index()].push_back(a.index());
        next_uid_edge_++;
        return Edge(this, next_uid_edge_-1);
    }
    else
        return Edge(this, next_uid_edge_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    next_uid_node_ = 0;
    edges_.clear();
    next_uid_edge_ = 0;
    adj_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return node that iterator point. */
    Node operator*() const {
        return graph_->node(uid_);
    }
    
    /** Return next iterator. */
    NodeIterator& operator++() {
        ++uid_;
        return *this;
    }
    
    /** Check whether two iterators are same. */
    bool operator==(const NodeIterator& node_iter) const {
        if ( !(graph_ == node_iter.graph_) )
            std::cout << "Check node interator== : same graph?" << std::endl;
        if ( graph_ == node_iter.graph_)
            return uid_==node_iter.uid_;
        return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_; //id of current iterating node
    NodeIterator(const Graph* graph, const size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return iterator that begins iteration of Node Iteator. */
  node_iterator node_begin() const {
    return NodeIterator(this,0);
  }

  /** Return iterator that ends iteration of Node Iteator. */
  node_iterator node_end() const {
    return NodeIterator(this,num_nodes());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    /** Return edge that iterator point. */
    Edge operator*() const {
        return graph_->edge(graph_->node(uid_node_),graph_->node(uid_node_,uid_adj_));
    }

    /** Return next iterator. */
    IncidentIterator& operator++() {
        ++uid_adj_;
        return *this;
    }
    /** Check whether two iterators are same. */
    bool operator==(const IncidentIterator& incident_iter) const {
        if ( !(graph_ == incident_iter.graph_) )
            std::cout << "Check incident interator== : same graph?" << std::endl;
        if ( graph_ == incident_iter.graph_)
            return (uid_node_==incident_iter.uid_node_ && uid_adj_ == incident_iter.uid_adj_);
        return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE 
    
    Graph* graph_;
    size_type uid_node_; //given node that we are iterating about
    size_type uid_adj_; //other node of incident edge of the given node
    IncidentIterator(const Graph* graph, const size_type uid_node, const size_type uid_adj)
        : graph_(const_cast<Graph*>(graph)), uid_node_(uid_node), uid_adj_(uid_adj) {
        }
  };
  

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    /** Return edge that iterator point. */
    Edge operator*() const {
        return graph_->edge(uid_);
    }

    /** Return next iterator. */
    EdgeIterator& operator++() {
        ++uid_;
        return *this;
    }
    
    /** Check whether two iterators are same. */
    bool operator==(const EdgeIterator& edge_iter) const {
        if ( !(graph_ == edge_iter.graph_) )
            std::cout << "Check edge interator== : same graph?" << std::endl;
        if ( graph_ == edge_iter.graph_)
            return uid_==edge_iter.uid_;
        return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type uid_; //id of current iterating node
    EdgeIterator(const Graph* graph, const size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return iterator that begins iteration of Edge Iteator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }
  
  /** Return iterator that ends iteration of Edge Iteator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }

  
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    struct internal_node {
        Point position;
        node_value_type value;
    };
    struct internal_edge {
        size_type uid1;
        size_type uid2;
    };
  // HW1: YOUR CODE HERE
    std::vector<internal_node> nodes_; // HW0 modified for 'scalable' nodes
  // adj_ stores adjecent nodes uid for each node -> adj_[0] stores vectors of all incident nodes of node(0)
    std::vector<std::vector<size_type>> adj_; 
    size_type next_uid_node_;
    std::vector<internal_edge> edges_;
    size_type next_uid_edge_;
    internal_node local_node_;
    internal_edge local_edge_;
    
public:
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
