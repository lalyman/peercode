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
template<typename V>
class Graph {
 private:
  struct nodeinfo{
    Point p_;
    V v_;
  };
  std::vector<nodeinfo> Nodes;
  std::vector<std::vector<unsigned>> Edges;
  std::vector<std::vector<unsigned>> Adj;  
  unsigned size0;
  unsigned num_edges0;

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
  using node_value_type = V;
  

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    size0 = 0;
    num_edges0 = 0;
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
  class Node: private totally_ordered<Node> {
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
      uid_n = 0;
      //graph_n = nullptr;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      const Point& result_point = graph_n->Nodes[uid_n].p_;
      return result_point;
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_n;
      //return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return reference of value of the node  
     */
    node_value_type& value(){
      node_value_type& result_v = graph_n->Nodes[uid_n].v_;
      return result_v;
    }

    /** Return a const reference of value of the node
     */ 
    const node_value_type& value() const{
      const  node_value_type& result_v = graph_n->Nodes[uid_n].v_;      
      return result_v;
}

    /** Return number of edges incident to the node */
    size_type degree() const{
      return graph_n->node_degree(uid_n);
    }

    /** Return begin of incident iterator */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_n, uid_n, 0);
    }

    /** Return end of incident iterator of this node */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_n, uid_n, degree());
    } 

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_n == n.graph_n) && (uid_n == n.index()));
      //(void) n;          // Quiet compiler warning
      //return false;
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
      assert(graph_n == n.graph_n);
      return uid_n < n.index();
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type uid_n;
    graph_type* graph_n;

    Node(const graph_type* graph0, size_type uid0): graph_n(const_cast<graph_type*>(graph0)), uid_n(uid0){
    } 
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return number of edges incident to the node with unique id @a id0 */
  size_type node_degree(size_type id0) const{
    assert(id0 < size0);
    return Adj[id0].size();
  }
 
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size0;
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
   *
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    size0 += 1;
    Point result_point= const_cast<Point&>(position);
    Nodes.push_back(result_point);
    return Node(this, size0 - 1);
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
  }*/


  /** Add a node to the graph, returning the added node.
   *  @param[in] position: The new node's position
   *  @param[in] value0: The new node's value. Optional.
   *  @post new num_nodes() == old num_nodes()
   *  @post result_node.index() == old num_nodes()
   *  Complexity: O(1)
   */
  Node add_node(const Point& position, const node_value_type& value0 = node_value_type()){
    size0 += 1;
    nodeinfo result_nodeinfo;
    result_nodeinfo.p_ = const_cast<Point&>(position);
    result_nodeinfo.v_ = const_cast<node_value_type&>(value0);
    Nodes.push_back(result_nodeinfo);
    Adj.push_back({});
    return Node(this, size0 - 1);
}

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.index() < size0);
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((0 <= i) && (i < size0));
    return Node(this, i);
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      node1_id = 0;
      node2_id = 0;
      edge_id = 0;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_e->node(node1_id);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_e->node(node2_id);      // Invalid Node
    }
   
    /** Swap node1 and node2 */
    void swap() { 
      int temp = node1_id;
      node1_id = node2_id;
      node2_id = temp;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1()) && (node2() == e.node2())) || ((node1() == e.node2()) && (node2() == e.node1())) ;
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((node1().index() + node2().index()) < (e.node1().index() + e.node2().index()));
      //(void) e;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_e;
    size_type node1_id;
    size_type node2_id;
    size_type edge_id;
    Edge(const graph_type* graph_e0, size_type node1_id0, size_type node2_id0, size_type edge_id0): graph_e(const_cast<graph_type*>(graph_e0)), node1_id(node1_id0), node2_id(node2_id0), edge_id(edge_id0){
      assert(node1_id0 < graph_e->size0 && node2_id0 < graph_e->size0 && node1_id0 != node2_id0);
    }
    // HW0: YOUR CODE HERE
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
    return num_edges0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges0);
    //if(i < num_edges0) return Edge();
    return Edge(this, Edges[i][0], Edges[i][1], i);
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_n == b.graph_n);
    size_type a0 = a.index();
    size_type b0 = b.index();
    bool test = false;
    size_type i = 0;
    while(i < size0){
      test = ((a0 == Edges[i][0] && b0 == Edges[i][1]) || (a0 == Edges[i][1] && b0 == Edges[i][0]));
      if(test){
        return true;
      }
      i++;
    }
    return false;
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;
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
    assert(a.graph_n == b.graph_n);
    assert(has_node(a) && has_node(b));
    size_type a0 = a.index();
    size_type b0 = b.index();
    bool test = false;
    size_type i = 0;
    while(i < num_edges0){
      test = ((a0 == Edges[i][0] && b0 == Edges[i][1]) || (a0 == Edges[i][1] && b0 == Edges[i][0]));
      if(test){
        return Edge(this, a.index(), b.index(), i);
      }
      i++;
    }
    num_edges0 += 1;
    std::vector<size_type> temp = {a0, b0};
    Edges.push_back(temp);
    Adj[a0].push_back(num_edges0 - 1);
    Adj[b0].push_back(num_edges0 - 1);
    return Edge(this, a0, b0, num_edges0 - 1);
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Nodes.clear();
    Edges.clear();
    size0 = 0;
    num_edges0 = 0; 
    // HW0: YOUR CODE HERE
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>  {
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

    /** Dereference NodeIterator
     *  @post: return the Node pointed by NodeIterator*/
    Node operator*() const{
      //assert(itr_ni < graph_ni->size());
      return graph_ni->node(itr_ni);
    }

    /** Move the iterator forward by 1
     *  @post: return the node_iterator that points to the next node in graph*/ 
    node_iterator& operator++(){
      ++itr_ni;
      return *(this);
    }

    /** Test whether this NodeIterator equals to @a n_itr0.
     *  Equal NodeIterators point to the same Node. */
    bool operator==(const node_iterator& n_itr0) const{
      //return graph_ni->node(itr_ni) == *n_itr0;
      return  graph_ni == n_itr0.graph_ni && itr_ni == n_itr0.itr_ni; 
}



   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_ni;
    size_type itr_ni;
    /** Construct a NodeIterator. @a itr_ni is initialized as 0 by default. */
    NodeIterator(const graph_type* graph0): graph_ni(const_cast<graph_type*>(graph0)), itr_ni(0){
    }
    NodeIterator(const graph_type* graph0, size_type itr0): graph_ni(const_cast<graph_type*>(graph0)), itr_ni(itr0){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return begin of node_iterator in this graph */
  node_iterator node_begin() const{
    return NodeIterator(this);
  }

  /** Return end of node_iterator in this graph */
  node_iterator node_end() const{
    return NodeIterator(this, size0);
  }
  

  //
  // Incident Iterator
  //

  /** A public interface to use @a Adj 
   *  @a n_id0 is index of the first node
   *  @a i is ith neighbour */
  size_type Get_Adj(size_type n_id0, size_type i) const {
    assert(n_id0 < size0);
    assert(i < Adj[n_id0].size());
    return Adj[n_id0][i];
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    
    /** Dereference IncidentIterator */
    Edge operator*() const{
      assert(node_uid < graph_ii->size());
      Node ntemp = graph_ii->node(node_uid);
      assert(itr_ii < ntemp.degree());
      Edge result = graph_ii->edge(graph_ii->Get_Adj(node_uid, itr_ii));
      if(result.node1() != ntemp){
        result.swap();
      }
      return result;
}


    /** Move Incident Iterator forward by one */
    IncidentIterator& operator++() {
      ++itr_ii;
      return *(this);
    }

    /** Return true if this incident iterator equals to @a itr_ii0 */
    bool operator==(const IncidentIterator& itr_ii0) const {
      return (graph_ii == itr_ii0.graph_ii) && (node_uid == itr_ii0.node_uid) && (itr_ii == itr_ii0.itr_ii);
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_ii;
    size_type node_uid;
    size_type itr_ii;
    IncidentIterator(const graph_type* graph_ii0, size_type node_uid0, size_type itr_ii0): graph_ii(const_cast<graph_type*>(graph_ii0)), node_uid(node_uid0), itr_ii(itr_ii0){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>  {
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

    /** Dereference Edge Iterator. Return the edge it points to. */
    Edge operator*() const{
      return graph_ei->edge(itr_ei);    
    }

    /** Move the itrator forward by one */
    EdgeIterator& operator++(){
      ++itr_ei;     
      return *(this);
    }

    /** Return true if this == @a e_itr0 */
    bool operator==(const EdgeIterator& e_itr0) const{
      //return graph_ei->edge(itr_ei) == *e_itr0;
      return graph_ei == e_itr0.graph_ei && itr_ei == e_itr0.itr_ei;
    }  

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_ei;
    size_type itr_ei;
    EdgeIterator(const graph_type* graph_ei0, size_type itr_ei0): graph_ei(const_cast<graph_type*>(graph_ei0)), itr_ei(itr_ei0){
    }
    EdgeIterator(const graph_type* graph_ei0): graph_ei(const_cast<graph_type*>(graph_ei0)), itr_ei(0){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /** Return edge iterator that points to begin */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Return edge iterator that points to end */
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges0);
  }

};

#endif // CME212_GRAPH_HPP
