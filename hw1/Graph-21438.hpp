#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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
  Graph() {
    num_edges_ = 0;
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
     /* Forward the request to graph_, which looks up the
      * appropriate internal node's position
      */
      return graph_->fetch_node(id_).p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return a non-const reference to the value of a node
    * @pre 0 <= id_ < graph.size()
    * @return A non-const reference to the value of the node
    * 
    * Complexity: O(1)
    */
    node_value_type& value(){
      return graph_->fetch_node(id_).value;

    }

    /** Return a const reference to the value of a node
    * @pre 0 <= id_ < graph.size()
    * @return A const reference to the value of the node
    *
    * Complexity: O(1)
    **/

    const node_value_type& value() const{
      return graph_->fetch_node(id_).value;
    }

    /** Return the degree of a node
    * @pre 0<= id_ < graph.size()
    *
    * Complexity: O(1)
    **/
    
    size_type degree() const{
      return graph_->internal_edges[id_].size();

    }

    /** Return an beginning incident_iterator 
    * @pre Node is either valid or one past the end
    * @post If Node is valid, *result = internal_edges[(*this).get_id()][0]
    *
    * @return A beginning incident_iterator
    *
    **/

    incident_iterator edge_begin() const{
      return IncidentIterator(*this, 0);
    }

    /** Return an ending incident_iterator
    * @pre Node is valid
    * @post result.id_ == this->degree()
    * @return An ending incident_iterator
    */

   incident_iterator edge_end() const {
     return IncidentIterator(*this, this->degree());
   }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.get_graph() == graph_) && (n.get_id() == id_);
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
      return id_ < n.get_id();
    }

    size_type get_id() const {
      return id_;
   }

   const graph_type* get_graph() const{
     return graph_;
   }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class NodeIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A node is characterized by its graph and its index
    graph_type* graph_; 
    size_type id_;

    Node(const graph_type* g, size_type id) : graph_(const_cast<graph_type*>(g)),
                                              id_(id) {}
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return internal_nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val (optional) The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */

  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    InternalNode n = {position, val};
    internal_nodes.push_back(n);
    std::vector<InternalEdge> new_edges {};
    internal_edges.push_back(new_edges);
    return Node(this, this->size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.get_graph() == this) && (n.get_id() < this->size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i && i < this->size());
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Forwards the request to graph_
      // which looks up the appropriate internal
      // edge and returns its first node n1_
      return n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->internal_edges[n1_.get_id()][id_].n2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // First check that the two edges are from the same graph
      if(e.get_graph() != graph_) {
        return false;
      }
      else {
        Node n1_ = this->node1();
        Node n2_ = this->node2();
        
        //Edges are undirected, so check both orientations
        bool ans = (n1_ == e.node1() && n2_ == e.node2());
        ans = ans || (n1_ == e.node2() && n2_ == e.node1());
        return ans; 
      }
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return n1_ < e.node1();
    }

    const graph_type* get_graph() const{
      return graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Edges are characterized by their graph, their first node
    // and their index in the adjacency list at that node

    const graph_type* graph_;
    const node_type n1_;
    size_type id_;    

    Edge(const Node& n1, size_type id) : 
        n1_(n1),
        id_(id) { graph_ = n1.get_graph();}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
      return num_edges_;
    }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    bool has_edge(const Node& a, const Node& b) const {
      // Iterate over all the internal edges, checking each to see
      // if any has endpoints @a a and @a b
      for(auto e : internal_edges[a.get_id()]){
        if(e.n2_ == b){
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
     // Similar to has_edge, iterate over internal_edges
     // returning a proxy associated with the InternalEdge with
     // endpoints a and b if such an InternalEdge exists
     for(size_type i = 0; i < a.degree(); ++i){
        if(internal_edges[a.get_id()][i].n2_ == b){
          return Edge(a, i);
        }
      }
     // If we don't find an edge with these endpoints,
     // we need to add it, increment num_edges_
     // and return a proxy associated with it
  
     internal_edges[a.get_id()].push_back(InternalEdge(a,b));
     internal_edges[b.get_id()].push_back(InternalEdge(b,a));
     ++num_edges_;
     return Edge(a, internal_edges[a.get_id()].size() - 1);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    void clear() {
      internal_edges.clear();
      internal_nodes.clear();
      num_edges_ = 0;
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the node pointed to by this iterator
    * @pre NodeIterator is not one past the end
    *
    * @return A Node pointed to by this iterator
    */
    
    Node operator*() const{
      return Node(graph_, id_);
    }

    /** Increment the NodeIterator 
    * @pre NodeIterator is not one past the end
    *
    * @return The incremented NodeIterator
    */
    node_iterator& operator++(){
      ++id_;
      return *this;
    }

    /** Check equality with another NodeIterator
    *
    * @return true if any only if dereferencing gives equal nodes
    */

    bool operator==(const node_iterator& ni) const{
      return *(*this) == *ni;
    }

   const graph_type* get_graph() const {
     return graph_;
   }

   size_type get_id() const {
     return id_;
   }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type id_; 
    NodeIterator(const graph_type* g, size_type id) : graph_(const_cast<graph_type*>(g)), id_(id) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return a beginning node_iterator
   * 
   * @post *return == this->internal_nodes[0]
   * @return a beginning node_iterator
   */
 
  node_iterator node_begin() const {
     return NodeIterator(this, 0);
  }
  
  /** Return an ending node_iterator
   *
   * @return an ending node_iterator
   */

  node_iterator node_end() const{
    return NodeIterator(this, this->size());
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return an Edge proxy for where this iterator points
     * @pre IncidentIterator is not one past the end
     *
     * @return An Edge proxy for where this iterator points 
     */

    Edge operator*() const {
      return Edge(node_, id_);
    }

    /** Increment the IncidentIterator
     * @pre IncidentIterator is not one past the end
     * @post This iterator points to the next edge in the adjacency list
     * @return The incremented iterator
     */

    incident_iterator& operator++(){
      ++id_;
      return *this;
    }
    
    /** Check equality of IncidentIterators
     *
     * @return true if and only if the two IncidentIterators have
     *          the same node_ and same id_
     */

    bool operator==(const incident_iterator& iit) const {
      return (id_ == iit.get_id()) && (node_ == iit.get_node());        
    }

    size_type get_id() const {
      return id_;
    }
   
    node_type get_node() const {
      return node_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type id_;
    node_type node_;

    IncidentIterator(const node_type node, size_type id = 0) : id_(id), node_(node) {}
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

     /** Return the edge pointed to
      * @pre inc_iter_ is not one off the end
      *
      * @return The Edge pointed to by the iterator
      */
     
     Edge operator*() const {
       return *inc_iter_;   
     }

     /** Increment the EdgeIterator
      * @pre EdgeIterator is not one past the end
      * @brief This increment skips Edges e where e.node1() <= e.node2()   
      *
      * @return The incremented EdgeIterator
      */

     EdgeIterator& operator++() {

       // Regardless of the edge we are currently at, we need to 
       // increment the internal iterators at least once
       internal_increment();

       // If niter_ is at the last node, we are at the end and should stop incrementing
       // Otherwise, we need to increment if inc_iter_ is at the end or if the current edge has second node
       // at least as big as the first

       while((niter_ != (*niter_).get_graph()->node_end())  && (inc_iter_ == (*niter_).edge_end() || ((*inc_iter_).node1() > (*inc_iter_).node2()))) {
         internal_increment();
       }
       return *this;
     }

     /** Check equality of EdgeIterators
      *
      *
      * @return If either iterator is one past the end, true iff their niter_'s are equal
      *         otherwise, true iff both internal iterators are equal
      */

     bool operator==(const EdgeIterator& eiter) const {
       // If either niter_ is at the end, the iterator is one past the end and
       // we don't want to compare inc_iter_
       bool off_the_end = niter_ == (*niter_).get_graph()->node_end();
       off_the_end = off_the_end || (eiter.get_niter() == (*(eiter.get_niter())).get_graph()->node_end());

       if(off_the_end) {
         return (niter_ == eiter.get_niter());
       }
       // Otherwise, both internal iterators must match
       return (inc_iter_ == eiter.get_inc_iter()) && (niter_ == eiter.get_niter());
     }

     
     //Increment the internal iterators to advance one edge in the adjacency list
     void internal_increment() {
       if(inc_iter_ != (*niter_).edge_end()){
         ++inc_iter_;
       }
       else{
         ++niter_;
         inc_iter_ = (*niter_).edge_begin();
       }
     }

     const NodeIterator& get_niter() const {
       return niter_;
     }    

     const IncidentIterator& get_inc_iter() const{
       return inc_iter_;
     }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    NodeIterator niter_;
    IncidentIterator inc_iter_;
    EdgeIterator(NodeIterator niter, IncidentIterator inc_iter) :
      niter_(niter),
      inc_iter_(inc_iter){}
      

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

   /** Return a beginning edge_iterator
    * @pre graph contains at least one edge
    * @post *return_value == internal_edges[0][0]
    * @return A beginning edge_iterator
    */

   edge_iterator edge_begin() const {
     node_iterator first = node_begin();
     return EdgeIterator(first, (*first).edge_begin()); 
   }

   /** Return an ending edge_iterator
    * @return An ending edge_iterator
    */

  edge_iterator edge_end() const {
    node_iterator last = node_end();
    return EdgeIterator(last, (*last).edge_end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
   size_type num_edges_;

   struct InternalNode{
     Point p;
     node_value_type value;
     InternalNode(Point pos, V val) : p(pos), value(val) {}
   };

   struct InternalEdge{
     node_type n1_;
     node_type n2_;

     InternalEdge(Node n1, Node n2) : n1_(n1), n2_(n2) {}

   };
   // We store the internal objects in vectors for simplicity.
   // For the edges, an unordered map keyed by the endpoints
   // would be an alternative that would allow has_edge and
   // add_edge to run in O(1).

   std::vector<InternalNode> internal_nodes;
   std::vector<std::vector<InternalEdge> > internal_edges; //Adjacency list

   const InternalNode& fetch_node(size_type id) const{
     assert (0 <= id && id < this->size());
     return internal_nodes.at(id);
   }

   InternalNode& fetch_node(size_type id) {
     assert (0 <= id && id < this->size());
     return internal_nodes.at(id);
   }

   const InternalEdge& fetch_edge(const Node& n, size_type id) const{
     assert (0 <= id && id < num_edges_);
     return internal_edges[n.get_id()].at(id);
   }

};


#endif // CME212_GRAPH_HPP
