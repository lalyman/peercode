#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
  using node_value_type = V;
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
      return graph->actualnodes[nodeindex];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nodeindex;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    /* return value of current node
    * @post some value of node_value_type 
    */
    node_value_type& value() {
      return graph->values[nodeindex];
    }
    // const node_value_type& value() const;
    /*
    * return value of current node 
    * @pre nodeindex < number of nodes  
    * @post some value of node_value_type 
    */
    const node_value_type& value() const {
      if (nodeindex < graph->size())
        return graph->values[nodeindex];
      assert(false);
    }
    // size_type degree() const;
    /*
    * return degree of node
    * @post result is a value of size_type, equal to the degree of node 
    */
    size_type degree() const{
      return graph->adj_lists[nodeindex].size();
    }
    //incident_iterator edge_begin() 
    /*
    * return iterator at first edge
    * @post result is an iterator, pointing to the first incident edge 
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph, nodeindex, size_type(0)); 
    } 
    //incident_iterator edge_end()
    /*
    * return iterator at last edge
    * @post result is an iterator, pointing to the last incident edge -> next 
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph, nodeindex, degree());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph == n.graph && nodeindex == n.index());
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
      return (graph == n.graph && nodeindex < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph;
    size_type nodeindex;

    Node(const Graph* g, size_type nindex) : graph(const_cast<Graph*>(g)), nodeindex(nindex) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return actualnodes.size();
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
  Node add_node(const Point& someposition, const node_value_type& somevalue = node_value_type(0)) {
    // HW0: YOUR CODE HERE
    actualnodes.push_back(someposition);  
    adj_lists[num_nodes() - 1];    //initialization 
    values.push_back(somevalue); 

    return Node(this, num_nodes()-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= 0 && i < num_nodes())
      return Node(this, i); 
    assert(false);        
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

    /*Edge(Graph* g, size_type a, size_type b, size_type l) {
      graph = g;
      nodeindex_1 = a;
      nodeindex_2 = b;
      edgeindex = l;
    }*/
    Edge(const Graph* g, size_type someindex, size_type nindex1, size_type nindex2) : graph(const_cast<Graph*>(g)), edgeindex(someindex), nodeindex_1(nindex1), nodeindex_2(nindex2) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph->node(nodeindex_1);    // get node at index      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(nodeindex_2);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2() || 
        node1() == e.node2() && node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (graph == e.graph && edgeindex < e.edgeindex);    // check that graph are the same 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph;
    size_type nodeindex_1;
    size_type nodeindex_2;   
    size_type edgeindex; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= 0 && i < num_edges())    // add pre-conditions 
      return Edge(this, i, edges[i][0], edges[i][1]); 
    assert(false);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b)) {
        size_type i1 = 0;
        size_type i2 = 0;
        for (size_type i = 0; i < num_nodes(); ++i) {
            if (a == node(i)) {
              i1 = i;
            } else if (b == node(i)) {
              i2 = i;
            }
        }
        for (size_type i = 0; i < num_edges(); ++i) {
          if (edges[i][0] == i1 && edges[i][1] == i2) {
            return true;
          }
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
    if (has_node(a) && has_node(b)) {
      std::vector<size_type> adj_list_a = adj_lists.at(a.index());    // search thru adj list
      for (size_type i = 0; i < adj_list_a.size(); ++i) {
        Edge e = edge(adj_list_a[i]);
        if (((a == e.node1()) && (b == e.node2())) 
           || ((a == e.node2()) && (b == e.node1()))) 
           return Edge(this, i, a.index(), b.index());
      }

      std::vector<size_type> some_edge;    // new edge 
      some_edge.push_back(a.index());
      some_edge.push_back(b.index());
      edges.push_back(some_edge);

      adj_lists[a.index()].push_back(num_edges()-1);
      adj_lists[b.index()].push_back(num_edges()-1);
      return edge(num_edges()-1);
    }

    assert(false);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    actualnodes.clear();
    values.clear();
    edges.clear();
    
    adj_lists.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    NodeIterator(const Graph* g, size_type someindex)
        : graph(const_cast<Graph*>(g)), nodeindex(someindex) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    /*
    * dereference
    * @post result is an element of type Node, that the iterator points to 
    */
    Node operator*() const {
      return graph->node(nodeindex);   
    }

    // NodeIterator& operator++()
    /*
    * iterator points to the next node 
    * @post result is an element of type NodeIterator pointing to the next node 
    */
    NodeIterator& operator++() {
      nodeindex++;
      return *this;
    }
    // bool operator==(const NodeIterator&) const
    /*
    * equality check 
    * @post result is bool, with true meaning that the iterators equal  
    */
    bool operator==(const NodeIterator& itr1) const {
      return (graph == itr1.graph && nodeindex == itr1.nodeindex);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph;
    size_type nodeindex;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  /*
  * return the iterator pointing to the first node
  * @post result is an element of type node_iterator, pointing to the first node 
  */
  node_iterator node_begin() const {
    return NodeIterator(this, size_type(0));
  }
  // node_iterator node_end() const
  /* 
  * return the iterator pointing to the end
  * @post result is an node_iterator, pointing to the end (not the last node)
  */
  node_iterator node_end() const {
    return NodeIterator(this, size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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


    IncidentIterator(const Graph* g, size_type index1, size_type index2) : graph(const_cast<Graph*>(g)), nodeindex(index1), adjindex(index2) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    /*
    * return an Edge that the incident iterator poitns to 
    * @pre nodeindex < number of nodes in graph 
    * @pre adjindex < number of incident edges of the node
    * @post result is an Edge, corresponding to the adjindex-th element of the adjacency edges of node corresponding to nodeindex 
    */
      Edge operator*() const {
        size_type edge_index = graph->adj_lists[nodeindex][adjindex];
        size_type node_index_1 = graph->edges[edge_index][0];
        size_type node_index_2 = graph->edges[edge_index][1];
        
        //if (nodeindex == node_index_1) {
          return Edge(graph, edge_index, nodeindex, node_index_2);    //directed - -!
        //} else {
        //  return Edge(graph, edge_index, node_index_2, nodeindex); 
        //}
      }

    // IncidentIterator& operator++()
      /*
      * move to the next
      * @pre nodeindex < number of nodes in graph 
      * @pre adjindex < number of incident edges of the node
      * @post result is an IncidentIterator
      */
      IncidentIterator& operator++() {
        adjindex++;
        return *this;
      }

    // bool operator==(const IncidentIterator&) const
      /*
      * equality check 
      * @para[in] i an object of IncidentIterator type
      * @post result is an boolean variable, equals true when the iterator i equals the current iterator 
      */
      bool operator==(const IncidentIterator& i) const {
         return ((graph == i.graph)
         && (nodeindex == i.nodeindex) 
         && (adjindex == i.adjindex));  
      }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
      Graph* graph;
      size_type nodeindex;
      size_type adjindex;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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
     // Edge operator*() const
    /*
    * return an Edge that the iterator points to 
    * @post result is an Edge that the iterator points to 
    */
    Edge operator*() const {
      return graph->edge(edgeindex); 
    }

    // EdgeIterator& operator++()
    /*
    * return an EdgeIterator
    * @post result is the next EdgeIterator
    */
    EdgeIterator& operator++() {
      edgeindex++;
      return *this;
    }

    // bool operator==(const EdgeIterator&) const
    /*
    * equality check 
    * @para[in] i  of type EdgeIterator
    * @post result is an boolean, and is true when the two iterators equal 
    */
    bool operator==(const EdgeIterator& i) const {
      return (graph == i.graph && edgeindex == i.edgeindex); 
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph;
    size_type edgeindex; 


    EdgeIterator(const Graph* g, size_type index)
        : graph(const_cast<Graph*>(g)), edgeindex(index) {
    }
  };  

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  /*
  * return an edge iterator pointing to the first edge
  * @post result is an edge_iterator, pointing to the first edge 
  */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, size_type(0));  
  }

  // edge_iterator edge_end() const
  /*
  * return an edge iterator pointing to the end
  * @post result is an edge iterator, pointing to the end (not the last edge, after that)
  */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges()); 
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<Point> actualnodes;
  std::vector<node_value_type> values;
  std::vector<std::vector<size_type>> edges;    // two indices on nodes 
  std::map<size_type, std::vector<size_type>> adj_lists;   // node index ~ edge index 
};

#endif // CME212_GRAPH_HPP
