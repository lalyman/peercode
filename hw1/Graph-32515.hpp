#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <list>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V>
class Graph  {
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


  /** Type of this graph. 
   * Making value_type generalizeable. */
  using node_value_type = V;
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
  Graph() {
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
      // HW0: YOUR CODE HERE
    }


    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return m_graph->m_node_points[m_index];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert( m_index < m_graph->size() ); //Is it correct to use an assert here? Is an if statement better?
                                           //How should an error be thrown if this is not met?
      // HW0: YOUR CODE HERE
      return m_index;
    }

    // HW1: YOUR CODE HERE
    /* Return a reference to this node's value 
     * Value returned here can be modified */
    node_value_type& value() {
        return m_graph->m_node_values[m_index];
    }; 

    /* Return a constant reference to this node's value, 
     * s.t. it cannot be modified */
    const node_value_type& value() const {
        return m_graph->m_node_values[m_index]; //Do I need to change the return type to reflect function?
    }; 

    // size_type degree() const;
    /* Check how many edges connect this node
     * to other nodes. */
    size_type degree() const {
        return m_graph->m_node_edges[m_index].size();
    }

    // incident_iterator edge_begin() const;
    /* Return an iterator point to the first edge
     * this node is connected to*/
    incident_iterator edge_begin() const {
        return IncidentIterator(m_graph, m_index, 0);
    }

    // incident_iterator edge_end() const;
    /* Return an iterator pointing one past the
     * last edge this node is connected to. */
    incident_iterator edge_end() const {
        return IncidentIterator(m_graph, m_index, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Checking whether the two graphs are equal
      return ( m_graph == n.m_graph && //Do I need to define operators == and < for Graphs now too?
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
      // HW0: YOUR CODE HERE
      if (m_graph == n.m_graph)
          return(index() < n.index() ) ;
      else
          return (m_graph < n.m_graph) ;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    size_type m_index;
    graph_type* m_graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Private constructor for Graph to construct valid Node objects
    Node(graph_type* graph, size_type index ) {
        m_graph = graph;
        m_index = index;
    }

    // If Node is passed a const graph_type*; need it to be volatile
    Node(const graph_type* graph, size_type index ) {
        m_graph = const_cast<graph_type*> (graph);
        m_index = index;
    }

  };//End node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return m_node_points.size();
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
    // HW0: YOUR CODE HERE
    m_node_points.push_back(position); 
    m_node_values.push_back(value);
    return Node(this,size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const  {
    // HW0: YOUR CODE HERE
    return (n.m_graph == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const  {
    // HW0: YOUR CODE HERE
    assert( i >= 0 && i < size() );
    return Node(this,i);  
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
      return Node(m_graph,m_index1);
      //return m_graph->get_edge_node(m_index,0); ///FIX: no function; just index and graph and Node()
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(m_graph,m_index2);
      //return m_graph->get_edge_node(m_index,1); ///FIX
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      auto ids_this = m_graph->sort(m_index1,m_index2);
      auto ids = m_graph->sort(e.m_index1,e.m_index2);
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
          auto ids_this = m_graph->sort(m_index1,m_index2);
          auto ids = m_graph->sort(e.m_index1,e.m_index2);
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
    // HW0: YOUR CODE HERE
    graph_type* m_graph;
    size_type m_index1;
    size_type m_index2;

    //Private Constructor
    Edge(graph_type* graph, size_type index1, size_type index2) {
        m_index1 = index1;
        m_index2 = index2;
        m_graph = graph;
    }

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

  }; // End Edge Class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return m_edges.size();
  }

  /** Return the edge with node indices @a i and @a j.
   * @pre 0 <= @a i < size()
   * @pre 0 <= @a j < size()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i,size_type j) const {
    // HW0: YOUR CODE HERE
    assert ( 0<= i && i < size() &&
             0<= j && j < size() &&
             i != j );
    std::vector< size_type > pair;
    pair = sort(i,j);
    return Edge(this,pair[0],pair[1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert( a.m_graph == this &&
            b.m_graph == this &&
            a.m_index != b.m_index ); // Edge cannot connect same nodes
    
    std::vector< size_type > pair;
    pair = sort(a.index(),b.index());
    auto i = find(m_edges.begin(), m_edges.end(), pair);
    if (i == m_edges.end())
        return false;
    else
        return true;
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
    assert( a.m_graph == this &&
            b.m_graph == this );
    std::vector<size_type> pair,bpair;// index1, index2;
    pair = sort(a.index(),b.index());
    bpair = bsort(a.index(),b.index());
    //this if/else ensures that edges are always written
    //in the convention of the smaller node index first 
    if (this->has_edge(a,b) == false)
    {   
        this->m_edges.push_back(pair);
        this->m_node_edges[a.index()].push_back(b.index());
        this->m_node_edges[b.index()].push_back(a.index());
    }
    return Edge(this,pair[0],pair[1]);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    m_edges.clear();
    m_node_points.clear();
    m_node_values.clear();
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
    Node operator*() const {
        return Node(m_graph, m_index);
    }

    // NodeIterator& operator++()
    node_iterator& operator++() {
        m_index++;
        return *this;
    }

    // bool operator==(const NodeIterator&) const
    bool operator==(const node_iterator& i) const {
        return (i.m_index == m_index &&
                i.m_graph == m_graph );
    }

   private:
    friend class Graph;
    size_type m_index;   // Knows node index for easy defining
    const graph_type* m_graph;
    NodeIterator(const graph_type* graph, size_type index) {
        m_index = index;
        m_graph = graph;

    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // node_iterator node_begin() const
  node_iterator node_begin() const {
      return NodeIterator(this,0);
  }

  // node_iterator node_end() const
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    edge_type operator*() const {
        return Edge(m_graph, m_node_id, m_graph->m_node_edges[m_node_id][m_it_id]);
    }

    incident_iterator& operator++() {
        m_it_id ++ ;
        return  *this;
    }

    bool operator==(const incident_iterator& iit) const {
        return (iit.m_it_id == m_it_id);
    }

   private:
    friend class Graph;
    graph_type* m_graph;
    size_type m_node_id;
    size_type m_it_id;
    IncidentIterator(graph_type* graph, size_type node_id, size_type it_id) { 
        //graph and index of spawning node
        m_graph = graph;
        m_node_id = node_id;
        m_it_id = it_id;

    }
    // HW1 #3: YOUR CODE HERE
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
    
    edge_type operator*() const {
        return *m_id2;
    }

    EdgeIterator& operator++() {
      ++m_id2;
      //While loop continues if the edge has been done
      //or if the incidentiterator has come to an end.
      while ( done.find(*m_id2) != done.end() ||
              m_id2 == (*m_id1).edge_end()  )
      {  
        //Go to next primary node if iit is at end  
        if (m_id2 == (*m_id1).edge_end())
        {   
            //if the new primary node is the end, need
            //to exit here (end of edge iterators defined as
            //the end of incident iterators of the last primary node
            if ( m_id1 == NodeIterator(m_graph,m_graph->size()-1) )
                {return *this;}
            else 
            {   
                ++m_id1;
                m_id2 = (*m_id1).edge_begin();
            }
        }
        //incrememnt incident iterator if edge already
        //has been done.
        if (done.find(*m_id2) != done.end() )
        {   
            ++m_id2;
        }
      } 
      done.insert(*m_id2);  
      return *this;  
    }

    bool operator==(const EdgeIterator& eit) const {
        return (m_id1 == eit.m_id1 &&
                m_id2 == eit.m_id2 );
    }

   private:
    friend class Graph;
    const graph_type* m_graph;
    NodeIterator m_id1;
    IncidentIterator  m_id2;
    std::set< Edge > done;
    EdgeIterator(const graph_type* graph, NodeIterator id1, IncidentIterator id2) {
        m_graph = graph;
        m_id1 = id1;
        m_id2 = id2;
    }

    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  //
  EdgeIterator edge_begin() const {
    //return EdgeIterator();
    NodeIterator id1 (this->node_begin());
    Node node1 (*id1);
    IncidentIterator id2 (node1.edge_begin());
    //graph_type* g (this);
    return EdgeIterator(this,id1,id2);
  }

  EdgeIterator edge_end() const {
    //return EdgeIterator();
    //NodeIterator id1 (this->node_end());//make 1 less than the node end
    //--id1;
    NodeIterator id1 (this,size()-1);//make 1 less than the node end
    Node node1 (*id1);
    IncidentIterator id2 (node1.edge_end());
    return EdgeIterator(this,id1,id2);
  }

 private:
//  std::vector<Point> m_nodes;
//  DOCUMENTATION
//
//
//
  std::vector<Point> m_node_points;
  std::vector<node_value_type> m_node_values;
  std::vector< std::vector<size_type> > m_edges;
  std::map< const size_type, std::vector<size_type> > m_node_edges;
  //size_type m_index;

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


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
