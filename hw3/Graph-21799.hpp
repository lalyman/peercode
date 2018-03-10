#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

//turn off asserts
//#define NDEBUG

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

template <typename V,typename E>
class Graph  {
 private:


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //


  /** Type of this graph. 
   * Making value_type generalizeable. */
  using node_value_type = V;
  using edge_value_type = E;
  using graph_type = Graph<V,E>;

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

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return m_graph->m_nodes[m_uid].idx;
    }

    /** Return this node's position (can be changed). */
    Point& position() {
      assert(valid());
      return m_graph->m_nodes[m_uid].position;
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(valid());
      return m_graph->m_nodes[m_uid].position;
    }

    /* Return a mutable reference to this node's value */
    node_value_type& value() {
        assert(valid());
        return m_graph->m_nodes[m_uid].value;
    }; 

    /* Return a immutable reference to this node's value */
    const node_value_type& value() const {
        assert(valid());
        return m_graph->m_nodes[m_uid].value;
    }; 

    /* Check how many edges connect this node
     * to other nodes. */
    size_type degree() const {
        assert(valid());
        return m_graph->m_nodes[m_uid].adj_nodes.size();
    }

    /* Return an iterator point to the first edge
     * this node is connected to*/
    incident_iterator edge_begin() const {
        assert(valid());
        return IncidentIterator(m_graph, index(), 0);
    }

    /* Return an iterator pointing one past the
     * last edge this node is connected to. */
    incident_iterator edge_end() const {
        assert(valid());
        return IncidentIterator(m_graph, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ( m_graph == n.m_graph && 
               index() == n.index() );
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     */
    bool operator<(const Node& n) const {
      if (m_graph == n.m_graph)
          return(index() < n.index() ) ;
      else
          return (m_graph < n.m_graph) ;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    size_type m_uid;
    graph_type* m_graph;

    /* For debugging.  Function that checks the representation invariants of node.
     * Can be checked for all functions within Node*/
    bool  valid()  const {
        bool val  (m_uid  >= 0 && 
            m_uid < m_graph->m_nodes.size() && 
            m_graph ->m_nodes[m_uid].idx  < m_graph->m_i2u.size()  &&
            m_graph->m_i2u[m_graph->m_nodes[m_uid].idx] == m_uid  
            );
        if (!val)
            std::cout << "Node UID: " << m_uid <<"; Idx: " << m_graph->m_nodes[m_uid].idx <<std::endl;
        return val;
    }

    // Private constructor for Graph to construct valid Node objects
    Node(graph_type* graph, size_type uid ) {
        if (uid < 0 || uid > graph->m_nodes.size())
            throw std::runtime_error("Trying to construct an invalid node");
        m_graph = graph;
        m_uid = uid;
    }

    // If Node is passed a const graph_type*; need it to be volatile
    // And every node is defined by both a graph and an index, so we
    // need both defined in the constructor
    Node(const graph_type* graph, size_type uid ) {
        m_graph = const_cast<graph_type*> (graph);
        m_uid = uid;
    }

  };//End node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
      return m_i2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new node's value; defaults to 
   *                     type's default value if none given
   * @post new num_node_points() == old num_nodes() + 1
   * @post new num_node_values() == old num_node_values() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    NodeStore node;
    node.position = position;
    node.value = value;
    node.idx = m_i2u.size();
    m_i2u.push_back(m_nodes.size());
    m_nodes.push_back(node);
    return Node(this,m_nodes.size()-1);
  }
  
  /** Remove a node @a n from the graph, returning 0 if 
   * @a n was not in graph or 1 if node was removed.
   * Also removes all edges of @a n.
   * @param[in] n        The node to be removed
   *
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old_num_edges() - n.degree()
   * @post new i.index() for i in range (n,end] == old i.index() - 1 
   *
   * Does not invalidate node or edge instances, but does invalidate all
   * iterators (node, edge and incident).
   *
   * Complexity: O( num_nodes+n.degree() ) amortized operations.
   *
   * Nodes structures are not explicitly deleted, instead an index
   * to uid mapper (m_i2u) and node indices are modified.
   */
  size_type remove_node(const Node& n){
      //Check if node is in graph
      if (has_node(n)==false) {
          //silently return if node to be deleted not in Graph
          return 0;}
      //Remove this node from all adjacent nodes adj_nodes lists
      while(n.edge_begin()!=n.edge_end()) {
           remove_edge(*(n.edge_begin()));}
      //Remove node from i2u mapping list
      m_i2u.erase(m_i2u.begin()+m_nodes[n.m_uid].idx);
      //Finally, decrement all higher nodes' indices
      for (auto i=n.m_uid; i<m_nodes.size(); i++) {
          m_nodes[i].idx--;     }
      return 1;
  }
  
  /** Remove the node of nodeiterator @a nit from the graph, returning 0 if 
   * @a n was not in graph or 1 if node was removed.
   * Also removes all edges of @a nit.
   * @param[in] nit        The node iterator for node to be removed
   *
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old_num_edges() - n.degree()
   * @post new i.index() for i in range (n,end] == old i.index() - 1 
   *
   * Complexity: O( num_nodes+(*nit).degree() ) amortized operations.
   */
  node_iterator remove_node(node_iterator nit ) {
      remove_node(*nit);
      return nit;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const  {
      return ( n.m_graph == this && 
              n.index() < size() );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const  {
      assert( i >= 0 && i < size() );
      return Node(this,m_i2u[i]);  
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

          //Is there a better way to write these const and non-const value() functions
          //such that the same code is not written twice?
          /* Returns this edge's value */
          edge_value_type& value() {
              //always access values in sorted order
              std::vector<size_type> pair = m_graph->sort(m_uid1,m_uid2);
              //searching over adj_nodes list of lower index node ( O(n.degree()) )
              auto it_start{ m_graph->m_nodes[pair[0]].adj_nodes.begin() };
              auto it_end{ m_graph->m_nodes[pair[0]].adj_nodes.end() };
              ptrdiff_t idx = find( it_start, it_end, pair[1]) - it_start;
              return m_graph->m_nodes[pair[0]].adj_edge_values[idx];
          }

          const edge_value_type& value() const {
              //always access values in sorted order
              std::vector<size_type> pair = m_graph->sort(m_uid1,m_uid2);
              auto it_start{ m_graph->m_nodes[pair[0]].adj_nodes.begin() };
              auto it_end{ m_graph->m_nodes[pair[0]].adj_nodes.end() };
              ptrdiff_t idx = find( it_start, it_end, pair[1]) - it_start;
              return m_graph->m_nodes[pair[0]].adj_edge_values[idx];
          }
          /* Return the length of this edge, ie euclidean distance between 
           * connecting nodes. */
          double length() const {
              Point P1 {node1().position()};
              Point P2 {node2().position()};
              return norm(P1-P2);
          }

          /** Return a node of this Edge */
          Node node1() const {
              return Node(m_graph,m_uid1);
          }

          /** Return the other node of this Edge */
          Node node2() const {
              return Node(m_graph,m_uid2);
          }

          /** Test whether this edge and @a e are equal.
           *
           * Equal edges represent the same undirected edge between two nodes.
           */
          bool operator==(const Edge& e) const {
            auto ids_this = m_graph->sort(m_uid1,m_uid2);
            auto ids = m_graph->sort(e.m_uid1,e.m_uid2);
            return (m_graph  == e.m_graph  &&
                    ids_this == ids        );
            }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
          if (m_graph == e.m_graph)
          {
              auto ids_this = m_graph->sort(m_uid1,m_uid2);
              auto ids = m_graph->sort(e.m_uid1,e.m_uid2);
              if (ids_this[0] == ids[0])
                  return (ids_this[1] < ids[1]);
              else
                  return (ids_this[0] < ids[0]);
          }
          else
              return ( m_graph < e.m_graph );
        }

       private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
        graph_type* m_graph;
        size_type m_uid1;
        size_type m_uid2;

        //Private Constructor
        Edge(graph_type* graph, size_type uid1, size_type uid2) {
            m_uid1 = uid1;
            m_uid2 = uid2;
            m_graph = graph;
        }

      
  }; // End Edge Class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return m_edge_ctr;
  }

  /** Return the edge with node indices @a i and @a j.
   * @pre 0 <= @a i < size()
   * @pre 0 <= @a j < size()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Complexity: O(1)
   */
  Edge edge(size_type i,size_type j) const {
    assert ( 0<= i && i < size() &&
             0<= j && j < size() &&
             i != j );
    return Edge(this,m_i2u[i],m_i2u[j]);
  }
  
  /** Return the edge with index idx.
   * Least preferred way to access edge, since in this representation
   * edges do not have indices and are instead defined by the indices 
   * of their nodes.
   * @pre 0 <= @a idx < num_edges()
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  Edge edge(size_type idx) const {
    assert( 0 <= idx && idx < num_edges() );
    auto eit = edge_begin();
    for (size_type i=0;i<=idx;i++)
    {
        ++eit;
    }
    return *eit;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert( a.m_graph == this &&
            b.m_graph == this  ); 
    std::vector<size_type> adj_nodes = m_nodes[a.m_uid].adj_nodes;
    if ( find( adj_nodes.begin(), adj_nodes.end(), b.m_uid) == adj_nodes.end() )
        return false;
    else
        return true;
  }

  bool has_edge(const Edge& e) const {
    node_type a{Node(this,e.m_uid1)};
    node_type b{Node(this,e.m_uid2)};
    return has_edge(a,b);
  }

/*  bool has_edge(size_type i, size_type j) const {
    node_type a{Node(this,i)};
    node_type b{Node(this,j)};
    return has_edge(a,b);
  }*/

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
  Edge add_edge(const Node& a, const Node& b,const edge_value_type& value = edge_value_type()) {
    assert( a.m_graph == this &&
            b.m_graph == this &&
            a.m_uid != b.m_uid );
    if (this->has_edge(a,b)==false)
    {
        //Add secondary node to primary nodes adjacent nodes list
        this->m_nodes[a.m_uid].adj_nodes.push_back(b.m_uid);
        this->m_nodes[a.m_uid].adj_edge_values.push_back(value);

        //Symmetric; Add primary node to secondary nodes adjacent nodes list
        this->m_nodes[b.m_uid].adj_nodes.push_back(a.m_uid);
        this->m_nodes[b.m_uid].adj_edge_values.push_back(value);
        
        //Keep track of added edges
        m_edge_ctr++;
    }
    return Edge(this,a.m_uid,b.m_uid);
  }
  
  /** Remove an edge @a e from the graph, returning 0 if @a e did not exist or 1
   * if @a e was removed successfully.
   * 
   * @post has_edge(@a e) == false
   * @post If old has_edge(@a e), new num_edges() == old num_edges().
   *       Else,                  new num_edges() == old num_edges() - 1.
   *
   * Does not invalidate edge or node instances, but does invalidate
   * edge iterators and incident iterators.
   *
   * Complexity: No more than O( degrees() of edge nodes 1 and 2 )
   */

  size_type remove_edge(const Edge& e) {
      
      if (has_edge(e)==false)
          return 0;

      //Remove edge from adjacency lists of both of its nodes
      std::vector<size_type>& adj_nodes = m_nodes[e.m_uid1].adj_nodes;
      size_type adj_idx = find(adj_nodes.begin(), adj_nodes.end(), e.m_uid2) - adj_nodes.begin();
      adj_nodes.erase(adj_nodes.begin()+adj_idx);
      m_nodes[e.m_uid1].adj_edge_values.erase(m_nodes[e.m_uid1].adj_edge_values.begin()+adj_idx);

      std::vector<size_type>& adj_nodes2 = m_nodes[e.m_uid2].adj_nodes;
      size_type adj_idx2 = find(adj_nodes2.begin(), adj_nodes2.end(), e.m_uid1) - adj_nodes2.begin();
      adj_nodes2.erase(adj_nodes2.begin()+adj_idx2);
      m_nodes[e.m_uid2].adj_edge_values.erase(m_nodes[e.m_uid2].adj_edge_values.begin()+adj_idx2);
      
      //Decrememnt edge counter
      m_edge_ctr--;

      return 1;
  }

  /* See specification for remove_edge(const Edge& e) above*/
  edge_iterator remove_edge(edge_iterator eit){
      remove_edge(*eit);
      return ++eit;
  }

  /* See specification for remove_edge(const Edge& e) above*/
  size_type remove_edge(const Node& a, const Node& b) {
      remove_edge(Edge(this,a.m_uid,b.m_uid));
      return 0;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    m_nodes.clear();
    m_i2u.clear();
    m_edge_ctr = 0;
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

    
    /* Returns dereferenced iterator*/
    Node operator*() const {
        return Node(m_graph, m_graph->m_i2u[m_index]);
    }

    /* Increments node iterator*/
    node_iterator& operator++() {
        m_index++;
        return *this;
    }   
    
    /* Decrement node iterator*/
    node_iterator& operator--() {
        m_index--;
        return *this;
    }

    /* Checks equality of this node iterator with another, i */
    bool operator==(const node_iterator& i) const {
        return (i.m_index == m_index &&
                i.m_graph == m_graph );
    }

   private:
    friend class Graph;
    size_type m_index;   
    const graph_type* m_graph;
    NodeIterator(const graph_type* graph, size_type index) {
        m_index = index;
        m_graph = graph;
    }
  };

  /* Returns first node iterator in a global order*/
  node_iterator node_begin() const {
      return NodeIterator(this,0);
  }

  /* Returns last node iterator in a global order*/
  node_iterator node_end() const {
      return NodeIterator(this, size());
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

    
    /* Returns edge this incident iterator is pointing to*/
    edge_type operator*() const {
        size_type uid1,uid2;
        uid1 = m_graph->m_i2u[m_nidx];
        uid2 = m_graph->m_nodes[m_graph->m_i2u[m_nidx]].adj_nodes[m_iidx] ;
        return Edge(m_graph,uid1 ,uid2 );
    }

    /* Increments incident iterator*/
    incident_iterator& operator++() {
        m_iidx ++ ;
        return  *this;
    }

    /* Checks equality between two incident iterators */
    bool operator==(const incident_iterator& iit) const {
        return (iit.m_nidx == m_nidx &&
                iit.m_iidx == m_iidx);
    }

   private:
    friend class Graph;
    graph_type* m_graph;
    size_type m_nidx;
    size_type m_iidx;
    IncidentIterator(graph_type* graph, size_type nidx, size_type iidx) { 
        m_graph = graph;
        m_nidx = nidx;
        m_iidx = iidx;

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

    
    /* return edge this edge iterator is pointing to*/
    edge_type operator*() const {
        return *m_iit;
    }

    /* Increment edge iterator using node and incident iterator functionality
     * Continues to increment until 2 conditions are satisfied:
     * a) incident iterator does not point to node.edge_end()
     * b) the edge given by (*eit) has node2.index() < node1.index()
     *    so as to skip duplicate edges
     * */
    EdgeIterator& operator++() {
      ++m_iit;
      while (true) {
          if (m_iit==(*m_nit).edge_end())
          {  
              ++m_nit;
              m_iit = (*m_nit).edge_begin();
              if ( m_nit == m_graph->node_end() )
              { return *this; }
          }
          else if ( (*m_nit).index() < (*m_iit).node2().index() )
          {  ++m_iit;}
          else
          { 
            return *this;
          }
      }
    }

    /* Check whether two edge iterators are equivalent*/
    bool operator==(const EdgeIterator& eit) const {
        return (m_nit == eit.m_nit &&
                m_iit == eit.m_iit );
    }

   private:
    friend class Graph;
    const graph_type* m_graph;
    NodeIterator m_nit;
    IncidentIterator  m_iit;
    EdgeIterator(const graph_type* graph, NodeIterator nit, IncidentIterator iit) {
        m_graph = graph;
        m_nit = nit;
        m_iit = iit;
    }

  };

  /* Start of edges defined as start of nodes and start of that
   * node's adjacency list*/
  EdgeIterator edge_begin() const {
    NodeIterator nit (this->node_begin());
    IncidentIterator iit ((*nit).edge_begin());
    return EdgeIterator(this,nit,iit);
  }

  /* End of edges defined as one past the last node in the graph (node_end)
   * and an incident iterator for that "node" set to its beginning */
  EdgeIterator edge_end() const {
    NodeIterator nit (this,size());
    IncidentIterator iit ((*nit).edge_begin());
    return EdgeIterator(this,nit,iit);
  }

 private:
//  DOCUMENTATION
//
// The data is represented as a vector of node data structures.
// Edges are not defined explicitly but are given by a member 
// list of nodes adjacent to a given node.  Edge values are stored
// similarly, with the index of the adjacent node equaling the 
// index of that edge's value.
//
// In this way, all edges are stored exactly twice.
//
// An edge counter member counts the number of edges for fast
// access.
// 
// A m_i2u member maps node indices to unique identifiers, so that 
// node structures are not deleted.
//
//
  struct NodeStore {
      node_value_type              value;
      Point                        position;
      size_type                    idx;
      std::vector<size_type>       adj_nodes;
      std::vector<edge_value_type> adj_edge_values;
  };

  std::vector<NodeStore> m_nodes;
  std::vector<size_type> m_i2u;
  size_type m_edge_ctr=0;
  
  //Sorts values in ascending order
  std::vector<size_type> sort(size_type i, size_type j) const {
      std::vector<size_type> pair;
      if (i < j)
      {
         pair.push_back(i);
         pair.push_back(j);
      }
      else 
      {
        pair.push_back(j);
        pair.push_back(i);
      }
      return pair;
  }

  //Sorts values in descending order
  std::vector<size_type> bsort(size_type i, size_type j) const {
      std::vector<size_type> pair;
      if (i > j)
      {
         pair.push_back(i);
         pair.push_back(j);
      }
      else 
      {
        pair.push_back(j);
        pair.push_back(i);
      }
      return pair;
  }


};

#endif // CME212_GRAPH_HPP
