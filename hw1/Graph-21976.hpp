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
template <typename V>
class Graph{

  public:
    using node_value_type = V;
    using size_type = unsigned;
    using graph_type = Graph<node_value_type>;
  
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

    struct internal_node {
      Point node_position;
      node_value_type node_val;
      std::vector<size_type> incident_node_idx;
      std::vector<size_type> incident_edge_idx;
    };

    struct internal_edge {
      size_type left_idx;
      size_type right_idx;
    };

    std::vector<internal_node> node_elements_;
    std::vector<internal_edge> edge_elements_;
    size_type edge_count_; // to keep track of the edge index for any edge of interest
    
  
  public:
  
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
  
    // CONSTRUCTORS AND DESTRUCTOR
    /** Construct an empty graph. */
    Graph(): node_elements_(), edge_elements_(){ }
  
    /** Default right_idxructor */
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
    const Point& position() const { return fetch_node().node_position; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return this->node_idx_; }

    /** Return the value associated with a node. */
    node_value_type& value() { return fetch_node().node_val; }

    /** Return the value associated with a node as a const. */
    const node_value_type& value() const { return fetch_node().node_val; }

    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_ == this->graph_ and n.node_idx_ == this->node_idx_);
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
      return (n.graph_ == this->graph_ and n.node_idx_ < this->node_idx_);
    }
    
    /** Returns the number of incident nodes to the node of interest. 
     */
    size_type degree() const{
      return fetch_node().incident_node_idx.size(); 
    }

    /** Constructs the iterator that points to the first element in 
     *  a vector of incident nodes to the current node of interest.
     */
    incident_iterator edge_begin() const {
        return incident_iterator(graph_, 0, this->index());
    }

    /** Constructs the iterator that points to one past the last element in 
     *  a vector of incident nodes to the current node of interest.
     */
    incident_iterator edge_end() const {
        return incident_iterator(graph_, this->degree(), this->index());
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_idx_;

    Node(const graph_type* graph, size_type node_idx):
        graph_(const_cast<graph_type*>(graph)), node_idx_(node_idx){}

    internal_node& fetch_node() const {
        return graph_->node_elements_[node_idx_];
    }



  };

  /** Return the number of nodes in the graph.
   */
  size_type size() const {
    return node_elements_.size();
  }
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
  Node add_node(const Point& position,
                const node_value_type& node_val) {

    internal_node new_node;
    new_node.node_position = position;
    new_node.node_val = node_val;
    node_elements_.push_back(new_node);

    return Node(this, this->num_nodes()-1);

  }
  Node add_node(const Point& position) {
    const node_value_type& node_val = node_value_type();
    return add_node(position, node_val);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  //--design_0
  //--Need to check that the node index is valid.
  bool has_node(const Node& n) const {
    return ( n.index() < this -> num_nodes() );
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
    Node node1() const {
        return Node(graph_,this->fetch_edge().left_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return Node(graph_,this->fetch_edge().right_idx);
    }



    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if  (e.node1().index() == this->node1().index()  and
           e.node2().index() == this->node2().index()){
          return true;
      } else if
          (e.node1().index() == this->node2().index()  and
           e.node2().index() == this->node1().index()) {
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
        return (this->edge_idx_ < e.edge_idx_);
    }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type edge_idx_;

    Edge(const graph_type* graph, size_type edge_idx):
        graph_(const_cast<graph_type*>(graph)),
        edge_idx_(edge_idx){}

    internal_edge& fetch_edge() const {
      return graph_->edge_elements_[edge_idx_];
    }

  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return edge_elements_.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { return Edge(this, i); }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {

    for(size_type i = 0; i < edge_elements_.size(); i++) {
      if  (a.index() == edge_elements_[i].left_idx and
           b.index() == edge_elements_[i].right_idx) {
           this->edge_count_ = i;
           return true;
      } else if
          (a.index() == edge_elements_[i].right_idx and
           b.index() == edge_elements_[i].left_idx) {
           this->edge_count_ = i;
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

    if (has_edge(a, b)) {
      return Edge(this, this->edge_count_);
    } else {

        // construct a new internal_edge
        internal_edge new_edge;
        new_edge.left_idx = a.index();
        new_edge.right_idx = b.index();
        edge_elements_.push_back(new_edge);

        // update the connectivity of that node
        size_type edge_id = this->num_edges()-1;
        this->node_elements_[a.index()].incident_edge_idx.push_back(edge_id);
        this->node_elements_[b.index()].incident_edge_idx.push_back(edge_id);
        this->node_elements_[a.index()].incident_node_idx.push_back(b.index());
        this->node_elements_[b.index()].incident_node_idx.push_back(a.index());

        return Edge(this, edge_id);
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_elements_.clear();
    edge_elements_.clear();
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
    bool operator==(const NodeIterator& rhs) const { return p_ == rhs.p_;}

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
  node_iterator node_end() const { return node_iterator(this,this->size()); }


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
        Node a = graph_->node(n_idx_);
        size_type temp = graph_->node_elements_[n_idx_].incident_node_idx[p_];
        Node b = graph_->node(temp); 
        return graph_->add_edge(a,b);
    }

    /** Increase the incident iterator by one*/
    IncidentIterator& operator++() {
        ++p_; 
        return *this;
    }

    /** Compares two incident iterators*/
    bool operator==(const IncidentIterator& iit) const {return p_ == iit.p_;}

   private:
    friend class Graph;
    graph_type* graph_;
    size_type p_;
    size_type n_idx_;
    
    IncidentIterator(const graph_type* graph, size_type point_idx, size_type node_idx):
        graph_(const_cast<graph_type*>(graph)), p_(point_idx), n_idx_(node_idx){}
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
    EdgeIterator() { }
 
    /** De-reference an edge iterator, returns the Edge the iterator points to.*/
    Edge operator*() const {return graph_->edge(p_);}

    /** Increase the iterator by one*/
    EdgeIterator& operator++() {++p_; return *this;}

    /** Compare two edge iterators */
    bool operator==(const EdgeIterator& eit) const {return p_ == eit.p_;}
 
   private:
    friend class Graph;
    graph_type* graph_;
    size_type p_;
    
    EdgeIterator(const graph_type* graph, size_type point_idx):
        graph_(const_cast<graph_type*>(graph)), p_(point_idx){}
  };

  /** Constructs the iterator that points to the first element in 
   *  a vector of edges.
   */
  edge_iterator edge_begin() const { return edge_iterator(this,0); }

  /** Constructs the iterator that points to one past the last element in 
   *  a vector of edges.
   */
  edge_iterator edge_end() const {return edge_iterator(this,this->num_edges());}

};

#endif // CME212_GRAPH_HPP




