#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
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
 
 struct internal_nodes;
 struct internal_edges;

 public:
    
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of user-specified value. */
  using node_value_type = V;

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
  Graph() 
      :nodes_() ,edges_() {
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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return (this->g_->nodes_[this->uid_].position);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return (this->uid_);
    }

    /** Return this node's value. */
    node_value_type& value() {
      return (this->g_->nodes_[this->uid_].node_val);
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return (this->g_->nodes_[this->uid_].node_val);
    }

    /** Return this node's number of adjacent edges. */
    size_type degree() const {
      return (this->g_->nodes_[this->uid_].incident_edges.size());
    }

    /** Return this node's first IncidentIterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this->g_, this->uid_);
    }

    /** Return this node's last IncidentIterator. */
    incident_iterator edge_end() const {
      return IncidentIterator(this->g_, this->uid_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.g_ == this->g_ && n.uid_ == this-> uid_);
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
      return (n.g_ == this->g_ && this->uid_ < n.uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* g_;
    size_type uid_;    

    Node(const Graph* g, size_type uid) 
        : g_(const_cast<Graph*>(g)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {    
    return this->nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& nval = node_value_type()) {
    internal_nodes node_addition;
    node_addition.position = position;
    size_type nodeid = this->nodes_.size();
    node_addition.node_val = nval;
    nodes_.push_back(node_addition);
    return Node(this, nodeid);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.g_ == this && n.uid_ < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i >= this->nodes_.size()) {return Node();}
    else {return Node(this,i);}      
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->g_, this->node_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->g_, this->node_b_); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->g_ == e.g_ && e.edgeid_ == this-> edgeid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    bool operator<(const Edge& e) const {
      return (this->edgeid_ < e.edgeid_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  
    Graph* g_;
    size_type edgeid_;
    size_type node_a_;
    size_type node_b_;    

    Edge(const Graph* g, size_type edgeid, size_type node_a, size_type node_b) 
        : g_(const_cast<Graph*>(g)), edgeid_(edgeid), node_a_(node_a), node_b_(node_b) {
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
    if (i >= edges_.size()) {return Edge();}
    else {
      size_type na= this->edges_[i].node_a;
      size_type nb= this->edges_[i].node_b;
      return Edge(this,i,na,nb);}
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
      Edge e = *ei;
      if ((e.node1() == a && e.node2() == b) || (e.node1() == b && e.node2() == a)) {
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
    internal_edges edge_addition;
    edge_addition.node_a = a.index();
    edge_addition.node_b = b.index();
    edge_addition.edgeid = num_edges();
    edges_.push_back(edge_addition);
    nodes_[a.index()].incident_edges.push_back(edge_addition.edgeid);
    nodes_[b.index()].incident_edges.push_back(edge_addition.edgeid);
      return Edge(this, edge_addition.edgeid, edge_addition.node_a,
                  edge_addition.node_b);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this->edges_.clear();
    this->nodes_.clear();
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

    /** Overloads the * operator to return the node of this NodeIterator */
    Node operator*() const {
      return (Node(this->iterator_graph, this->graph_position));
    }

    /** Overloads the ++ operator to increment the NodeIterator */
    NodeIterator& operator++() {
      this->graph_position++;
      return *this;
    }

    /** Overloads the == operator and returns true if the passed in NodeIterator 
    is the same as this NodeIterator */
    bool operator==(const NodeIterator& x) const {
      return (this-> iterator_graph == x.iterator_graph &&
            this->graph_position == x.graph_position);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* iterator_graph;
    size_type graph_position;    

    NodeIterator(const Graph* g, const size_type gpos = 0) 
        : iterator_graph(const_cast<Graph*>(g)), graph_position(gpos) {
        }
};
  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

     /** Returns the first NodeIterator. */
     node_iterator node_begin() const {
        return NodeIterator(this);
     }

     /** Returns the last NodeIterator. */
     node_iterator node_end() const {
        return NodeIterator(this, this->nodes_.size());
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
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Overloads the * operator to return the Edge of this IncidentIterator */
    Edge operator*() const {
      size_type thisedgeid = this->iterator_graph->nodes_[start_node_id].incident_edges[this->graph_position];
      size_type first_node = this->start_node_id;
      size_type second_node = this->iterator_graph->edges_[thisedgeid].node_b;
      if (second_node == first_node) {
        second_node = this->iterator_graph->edges_[thisedgeid].node_a;
      }
      return (Edge(this->iterator_graph, thisedgeid, first_node, second_node));
    }
    
    /** Overloads the ++ operator to increment the IncidentIterator */
    IncidentIterator& operator++() {
      this->graph_position++;
      return *this;
    }

    /** Overloads the == operator and returns true if the passed in IncidentIterator 
    is the same as this IncidentIterator */
    bool operator==(const IncidentIterator& x) const {
      return (this-> iterator_graph == x.iterator_graph &&
            this->graph_position == x.graph_position && this->start_node_id == x.start_node_id);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* iterator_graph;
    size_type start_node_id;
    size_type graph_position;

    IncidentIterator(const Graph* g, const size_type snodeid, const size_type gpos = 0) 
        : iterator_graph(const_cast<Graph*>(g)), start_node_id(snodeid), graph_position(gpos) {

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
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Overloads the * operator to return the Edge of this EdgeIterator */
    Edge operator*() const {
      size_type thisedgeid = this->graph_position;
      size_type first_node = this->iterator_graph->edges_[graph_position].node_a;
      size_type second_node = this->iterator_graph->edges_[graph_position].node_b;
      if (second_node == first_node) {
          second_node = this->iterator_graph->edges_[thisedgeid].node_a;
      }
      return (Edge(this->iterator_graph, thisedgeid, first_node, second_node));
    }


    /** Overloads the ++ operator to increment the EdgeIterator */
    EdgeIterator& operator++() {
      this->graph_position++;
      return *this;
    }

    /** Overloads the == operator and returns true if the passed in EdgeIterator 
    is the same as this EdgeIterator */
    bool operator==(const EdgeIterator& x) const {
      return (this-> iterator_graph == x.iterator_graph &&
              this->graph_position == x.graph_position);
    }

   private:
    friend class Graph;
    Graph* iterator_graph;
    size_type graph_position;

    EdgeIterator(const Graph* g, const size_type gpos = 0) 
        : iterator_graph(const_cast<Graph*>(g)), graph_position(gpos) {
        }

    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the first EdgeIterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this);
  }

  /** Return the last EdgeIterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->edges_.size());
  }

 private:
  //Internal type of nodes and edges
  struct internal_nodes {
    Point position;
    node_value_type node_val;
    std::vector<size_type> incident_edges;

  };
  
  struct internal_edges {
    size_type node_a;
    size_type node_b;
    size_type edgeid;
  };

  std::vector<internal_nodes> nodes_;
  std::vector<internal_edges> edges_;

};

#endif // CME212_GRAPH_HPP