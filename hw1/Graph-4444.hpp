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
  Graph(): nodes_(), nsize_(0) {
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
  class Node {
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
      return (graph_->nodes_[id_]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {return graph_->values_[id_];};
    const node_value_type& value() const {return graph_->values_[id_];};
    size_type degree() const {
        return graph_->adj_list_[this->id_].size();
        };
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_,this->id_,0);
        };
    incident_iterator edge_end() const {
        return IncidentIterator(graph_,this->id_,this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->graph_ == n.graph_ and this->id_ == n.id_) {
        return true;
      }
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
      if (this->id_ < n.id_ and this->graph_ == n.graph_) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_ {}; 
    size_type id_ {};
    Node(const Graph* graph, size_type id)
        : graph_(const_cast<Graph*>(graph)), id_(id) {}
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nsize_;
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
  /*Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(position);
    nsize_++; 
    // value_ = 0;  // default value for value_ 
    return Node(this,nsize_-1);        // Invalid node
  }*/

  Node add_node(const Point& position, \
  const node_value_type& value = node_value_type()) {
    nodes_.push_back(position);
    values_.push_back(value);
    nsize_++;
    // std::cout << "value=" << value << std::endl; 
    return Node(this,nsize_-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.id_<nsize_) return true; 
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
    return Node(this,i);        // Invalid node
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // return Node(graph_,graph_->edge_n1s_[id_]);      // Invalid Node
        return Node(graph_,this->n1_id_); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // return Node(graph_,graph_->edge_n2s_[id_]);      // Invalid Node
        return Node(graph_,this->n2_id_); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ == e.graph_ and (this->n1_id_ == e.n1_id_ or \
        this->n2_id_ == e.n2_id_)) return true; 
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->n1_id_ < e.n1_id_ and this->n2_id_ < e.n2_id_ \
        and this->graph_ == e.graph_) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_ {}; 
    size_type n1_id_ {}, n2_id_ {}; 
    Edge(const Graph* graph, size_type n1_id, size_type n2_id)
        : graph_(const_cast<Graph*>(graph)), n1_id_(n1_id), n2_id_(n2_id) {}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return esize_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this,this->edge_n1s_[i],this->edge_n2s_[i]);        
    // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type n1 = a.id_; 
    size_type n2 = b.id_; 
     
    for (size_type i=0; i<esize_; i++) {
        if ((n1==edge_n1s_[i] and n2==edge_n2s_[i]) \
            (n1==edge_n2s_[i] and n2==edge_n1s_[i])) return true;
        {
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    if (adj_list_.size()<=nsize_) {
        adj_list_.resize(nsize_);
    }
    for (size_type i=0; i<(adj_list_[a.id_].size()); i++) {
        if (b.id_==adj_list_[a.id_][i]) {
            return Edge(this,a.id_,b.id_);
        }
    }
    esize_++;
    edge_n1s_.push_back(a.id_);
    edge_n2s_.push_back(b.id_);
    adj_list_[a.id_].push_back(b.id_);
    adj_list_[b.id_].push_back(a.id_);
    return Edge(this,a.id_,b.id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nsize_ = 0; 
    esize_ = 0; 
    nodes_.clear(); 
    edge_n1s_.clear(); 
    edge_n2s_.clear(); 
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

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(graph_,iter_id_);
    };
    NodeIterator& operator++() {
      iter_id_++; 
      return(*this);
    };
    bool operator==(const NodeIterator& iter) const {
      return (this->graph_==iter.graph_ and this->iter_id_==iter.iter_id_);
    };

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; 
    size_type iter_id_; 

    NodeIterator(const Graph* graph, size_type iter_id) : \
    graph_(const_cast<Graph*>(graph)), iter_id_(iter_id) {};
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return (NodeIterator(this,0)); 
  };
  node_iterator node_end() const {
    return (NodeIterator(this,nsize_)); 
  };

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
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
        return (Edge(graph_,node_id_,graph_->adj_list_[node_id_][iter_id_]));
    };
    IncidentIterator& operator++() {
        this->iter_id_++; 
        return(*this);
    };
    bool operator==(const IncidentIterator& iter) const {
        if (this->graph_==iter.graph_ and this->node_id_==iter.node_id_ and \
            this->iter_id_==iter.iter_id_) return true; 
        return false; 
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_ {}; 
    size_type node_id_ {}; 
    size_type iter_id_ {};
    IncidentIterator (Graph* graph, size_type node_id, size_type iter_id):\
    graph_(const_cast<Graph*>(graph)), node_id_(node_id), iter_id_(iter_id) {}
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
    Edge operator*() const {
        return(Edge(graph_,graph_->edge_n1s_[iter_id_],\
            graph_->edge_n2s_[iter_id_]));
    };
    EdgeIterator& operator++() {
        ++iter_id_; 
        return (*this);
    };
    bool operator==(const EdgeIterator& iter) const {
        if (this->graph_==iter.graph_ and this->iter_id_==iter.iter_id_)
            return true; 
        return false; 
    };

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_ {}; 
    size_type iter_id_ {0}; 
    EdgeIterator(const Graph* graph,const size_type iter_id) \
        :graph_(const_cast<Graph*>(graph)), iter_id_(iter_id) {}; 
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    edge_iterator edge_begin() const {
        return EdgeIterator(this,0);
    };
    edge_iterator edge_end() const {
        return EdgeIterator(this,this->esize_);
    };

 private:
  // HW0: YOUR CODE HERE
    std::vector<Point>nodes_ {};
    size_type nsize_ {0};
    std::vector<node_value_type> values_ {}; 
    std::vector<std::vector<size_type>> adj_list_ {}; 
    std::vector<size_type>edge_n1s_ {};
    std::vector<size_type>edge_n2s_ {};
    size_type esize_ {0}; 
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
