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

template < typename V >
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the value stored into the nodes. */
  using node_value_type = V ;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    this->point_map = std::unordered_map< size_type, Point* >();
    this->node_map = std::unordered_map< size_type, Node* >();
    this->value_map = std::unordered_map< size_type, node_value_type* >();

    this->edge_map = std::unordered_map< size_type, std::unordered_map< size_type, Edge* > >();
    this->edge_indexation = std::unordered_map< size_type, Edge* >();
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
   * Nodes implement operators == and <, and get totally ordered through the totally ordered template
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
      this->graph = nullptr;
      this->node_index = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      return *(this->graph->point_map.at(this->node_index));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->node_index;
    }

    /** Return a modifiable reference to a Node's value. */
    node_value_type & value (){
      return *(this->graph->value_map.at(this->node_index));
    }

    /** Return a const reference to a Node's value. This cannot modify a node's value. */
    const node_value_type & value () const{
      return *(this->graph->value_map.at(this->node_index));
    }

    /** Return the degree of a node. O(1) operations. */
    size_type degree() const{
      return this->graph->edge_map[this->node_index].size();
    }

    /** Return an object able to iterate over the neighbours of n, pointing on the first neighbor. */   
    incident_iterator edge_begin() const{
      incident_iterator iterator;
      iterator.it= this->graph->edge_map[this->node_index].begin();
      return iterator;
    }

    /** Return an iterator over the end of the neighbors of n. */
    incident_iterator edge_end() const{
      incident_iterator iterator;
      iterator.it= this->graph->edge_map[this->node_index].end();
      return iterator;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->graph == n.graph && this->node_index == n.node_index;
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
      bool inferior;
      // if both Nodes belong to the same graph, check the order of their index
      // else, cast pointer graph to an integer, and compare the adresses as integers
      if(n.graph == this->graph){
        inferior = this->node_index < n.node_index;
      }
      else{
        auto i = reinterpret_cast<std::uintptr_t>(this->graph);
        auto j = reinterpret_cast<std::uintptr_t>(n.graph);

        inferior = i < j;
      }
      return inferior;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type node_index;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->node_map.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value The new node's value
   * 
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == node_value
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type & node_value = node_value_type ()) {
    Node* new_node = new Node();
    new_node->graph = this;
    new_node->node_index = this->size();

    // Since this->size() returns the size of node_map, we have to add the new node to this map last. 
    this->value_map[this->size()] = new node_value_type(node_value);
    this->point_map[this->size()] = new Point(position.x, position.y, position.z);
    this->node_map[this->size()] = new_node;

    return *new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    bool result = (this->node_map.count(n.index()) > 0);
    return result;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    Node* n_out;
    try{
      n_out = this->node_map.at(i);
    }
    catch(const std::out_of_range e){
      throw "Index not in graph";
    }
    return *n_out;
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   * Nodes implement operators == and <, and get totally ordered through the totally ordered template
   */
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      node_1 = nullptr;
      node_2 = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return *(this->node_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return *(this->node_2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (*(this->node_1) == *(e.node_2) && *(this->node_2) == *(e.node_2)) || (*(this->node_1) == *(e.node_2) && *(this->node_2) == *(e.node_1));
    }

    bool operator<(const Edge& e) const {
      bool inferior;
      if (*(this->node_1) < *(this->node_2)){
        if (*(e.node_1) < *(e.node_2)){
          inferior = (*(this->node_1) < *(e.node_1)) || (*(this->node_1) == *(e.node_1) && *(this->node_2) < *(e.node_2));
        }
        else{
          inferior = (*(this->node_1) < *(e.node_2)) || (*(this->node_1) == *(e.node_2) && *(this->node_2) < *(e.node_1));
        }
      }
      else{
        if (*(e.node_1) < *(e.node_2)){
          inferior = (*(this->node_2) < *(e.node_1)) || (*(this->node_2) == *(e.node_1) && *(this->node_1) < *(e.node_2));
        }
        else{
          inferior = (*(this->node_2) < *(e.node_2)) || (*(this->node_2) == *(e.node_2) && *(this->node_1) < *(e.node_1));
        }
      }
      return inferior;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const Node* node_1;
    const Node* node_2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edge_indexation.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    Edge* edg_out;
    try{
      edg_out = this->edge_indexation.at(i);
    }
    catch(const std::out_of_range e){
      throw "Edge not in graph";
    }
    return *edg_out;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool result;
    try{
      result = (this->edge_map.at(a.index()).count(b.index()) > 0);
    }
    catch(const std::out_of_range e){
      result = false;
    }
    return result;
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
    if (this->has_edge(a, b)){
      return *(this->edge_map.at(a.index()).at(b.index()));
    }
    else{
    Edge* new_edge = new Edge();
    new_edge->node_1 = &a;
    new_edge->node_2 = &b;

    this->edge_map[a.index()][b.index()] = new_edge;
    this->edge_map[b.index()][a.index()] = new_edge;

    this->edge_indexation[this->num_edges()] = new_edge;
    return *new_edge;
    }

    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // map.clear will delete all the pointers contained in the map

    this->edge_indexation.clear();
    this->node_map.clear();
    this->point_map.clear();
    this->value_map.clear();

    this->point_map = std::unordered_map< size_type, Point* >();
    this->node_map = std::unordered_map< size_type, Node* >();
    this->value_map = std::unordered_map< size_type, Node* >();
    this->edge_map = std::unordered_map<size_type, std::unordered_map<size_type, Edge* > >();
    this->edge_indexation = std::unordered_map< size_type, Edge* >();
  }

  //
  // Node Iterator
  //

  // Define the type of an iterator through a map
  typedef typename std::unordered_map<size_type, Node*>::const_iterator node_map_iterator;
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   * This will basically be equivalent to an iterator within graph.node_map.
   * Iterating through a unordered_map is basically O(n), with each operation being constant in time
   * */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

      /** The iterator is an iterator over the graph.node_map. Taking `second' yields the node. */
      Node operator *() const{
        return *(this->it->second);
      }
      /** Move the iterator forward. This is equivalent as calling ++ on the underlying iterator
       * over the graph.node_map. */
      node_iterator & operator ++(){
        (this->it)++;
        return *this;
      }
      /** This implicitely defines operator !=. */
      bool operator ==( const node_iterator & rhs) const{
        return this->it == rhs.it;
      }


    private:
    friend class Graph;
    // This is an iterator over a map similar to graph.node_map.
    node_map_iterator it;
  };

  /** Return an iterator over the nodes of the graph, pointing over the `first' element of graph.node_map. */
  node_iterator node_begin() const{
    node_iterator iterator;
    iterator.it= this->node_map.begin();
    return iterator;
  }

  /** Return an iterator pointing to the end of graph.node_map. */ 
  node_iterator node_end() const{
    node_iterator iterator;
    iterator.it = this->node_map.end();
    return iterator;
  }

  //
  // Incident Iterator
  //
  // Define the type of an iterator through a map
  typedef typename std::unordered_map<size_type, Edge*>::const_iterator edge_map_iterator;
  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. 
   * This will basically be equivalent to an iterator within graph.edge_map[n].*/
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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
      /** Returns the edge to which the iterator refers.
       * Makes sure that node_1 is the origin node.
       * this->it->second is the Edge to which it points
       * this->it->first is the index of the destination node */
      Edge operator *() const{
        if (this->it->second->node_1->index() == this->it->first){
          const Node* temp = this->it->second->node_1;
          this->it->second->node_1 = this->it->second->node_2;
          this->it->second->node_2 = temp;
        } 
        return *(this->it->second);
      }

      /** Increments the underlying pointer in the unordered_map. */
      IncidentIterator & operator ++(){
        (this->it)++;
        return *this;
      }

      /** This implicitely defines operator !=. */
      bool operator ==( const IncidentIterator & rhs) const{
        return this->it == rhs.it;
      }

   private:
    friend class Graph;
    edge_map_iterator it;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * This will basically be equivalent to an iterator within graph.edge_indexation*/
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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
    /** Returns the edge to which the iterator refers.
     * Does not care which node is the origin */
    Edge operator *() const{
        return *(this->it->second);
      }

      /** Increments the underlying pointer in the unordered_map. */
      EdgeIterator & operator ++(){
        (this->it)++;
        return *this;
      }
      
      /** This implicitely defines operator !=. */
      bool operator ==( const EdgeIterator & rhs) const{
        return this->it == rhs.it;
      }

   private:
    friend class Graph;
    edge_map_iterator it;

  };
     /** Return an iterator over the nodes of the graph, pointing over the `first' element of graph.edge_indexation. */
    edge_iterator edge_begin() const{
    edge_iterator iterator;
    iterator.it= this->edge_indexation.begin();
    return iterator;
  }

    /** Return an iterator pointing to the end of graph.edge_indexation. */ 
    edge_iterator edge_end() const{
    edge_iterator iterator;
    iterator.it= this->edge_indexation.end();
    return iterator;
  }

 private:

  std::unordered_map<size_type, Point* > point_map;
  std::unordered_map<size_type, Node* > node_map;
  std::unordered_map<size_type, node_value_type* > value_map;

  std::unordered_map<size_type, std::unordered_map<size_type, Edge*>> edge_map;
  std::unordered_map<size_type, Edge*> edge_indexation;
};

#endif // CME212_GRAPH_HPP
