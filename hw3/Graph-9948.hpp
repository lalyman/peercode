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

template <typename V, typename E>
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
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

  /** Syonym for type of vaiable a node stores */
  using node_value_type = V;

  /** Syonym for type of vaiable a edge stores */
  using edge_value_type = E;

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
    // HW0: YOUR CODE HERE
    nNodes = 0;
    nEdges = 0;

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
        ind = -1;
        parentGraph = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      return parentGraph->nodePosition(ind);
    }

    /** Return a modifiable position. */
    Point& position() {
      return parentGraph->node_list[ind].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return ind;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**
     * @brief Returns a modifiable value for the node.
     *
     * @return returns modifiable value for the node
     * @pre this is a valid node
     * @post node value set to modified value
     *
     */
    node_value_type& value(){
      return parentGraph->node_list[ind].value;
    }

    /**
     * @brief Returns a unchangeable value for the node.
     *
     * @return returns const variable value for the node
     * @pre this is a valid node
     * @post The value of the node is the same after function call
     *
     */
    const node_value_type& value() const{
      return parentGraph->node_list[ind].value;
    }

    /**
     * @brief Access to the number of edges connecting this node to others
     *
     * @return returns the degree of the node as a const variable
     * @pre this is a valid node
     * @post The value of the degree is the same after function call
     *
     */
    size_type degree() const{
      return parentGraph->adjacency_list[ind].size();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (parentGraph == n.parentGraph && n.index()==this->index());
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
      return ind<n.index();
    }

    /**
     * @brief Provide an iterator for the edges connecting to this node in a graph
     *
     * @return Iterator that points to the first edge connecting this node to another
     * @pre parentGraph != nullptr
     * @pre ind < nEdges
     * @post returned iterator has an index of 0
     *
     * Ordering of edges enforces strict weak ordering
     */
    incident_iterator edge_begin() const{
      assert( parentGraph != nullptr);
      return incident_iterator(parentGraph, ind, 0);
    }

    /**
     * @brief Provide the end iterator for the edges in a graph
     *
     * @return Iterator that points to the last edge in the graph
     * @pre parentGraph != nullptr
     * @pre ind < nEdges
     * @post returned iterator has an index of node.degree()
     *
     * Ordering of edges enforces strict weak ordering
     */
    incident_iterator edge_end() const{
      return incident_iterator(parentGraph,ind, this->degree());
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* parentGraph;
    size_type ind;

    /**
     * @brief Private constructor to be used by the Graph class to construct nodes
     * @param g Pointer to graph that is creating the node.
     * @param i index of node being created
     *
     * @pre @a g!= nullptr
     * @post new node created for graph @a g with indexx @a ind
     */

    Node(const Graph* g, size_type i) : parentGraph(const_cast<Graph*> (g)), ind(i) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return adjacency_list.size();
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

  Node add_node ( const Point & position, const node_value_type & val= node_value_type ()){
    node_type newNode (this, nNodes);
    nNodes++;
    Point p;
    p.x=position.x;
    p.y= position.y;
    p.z= position.z;
    internal_node temp {p,val};
    node_list.push_back(temp);
    std::vector<adjacencyInfo> emptyList;
    adjacency_list.push_back(emptyList);
    return newNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    return (n.index() < size() && n.parentGraph == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i<nNodes);
    //return node_list[i].node;
    return node_type(this,i);
  }

  /**
   * @brief Remove node @a n from the graph and delete all connecting edges
   * @param n Node that needs to be removed
   * @return number of nodes removed
   *
   * @pre node @a is in this graph
   * @post new graph.num_nodes() = old graph.num_nodes-1
   * @post new graph.num_edges() = old.graph.num_edges - @a n.degree()
   * @post invalidates all node iterators in the range [0, old graph.num_nodes]
   *        and edge iterators in the range [0, old graph.num_edges]
   * @post Invalidates and edges that are connected to @a n
   *
   * This has time complexity  of < O(num_nodes()) assuming sparse graph
   */
  size_type remove_node(const Node& n){
      assert(has_node(n));
      size_type ind = n.index();
      node_list[ind] = node_list.back();
      node_list.pop_back();

      for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
        remove_edge(*ei);
      }

      adjacency_list[ind] = adjacency_list.back();
      size_type prevInd = adjacency_list.size()-1;

      for(size_type i = 0; i < adjacency_list[ind].size(); ++i){
          size_type con = adjacency_list[ind][i].node2_index;
          for(size_type j = 0; j < adjacency_list[con].size(); ++j){
              if (adjacency_list[con][j].node2_index == prevInd){
                  adjacency_list[con][j].node2_index = ind;
              }
          }
      }
      adjacency_list.pop_back();
      --nNodes;
      return 1;

  }

  /**
   * @brief Remove node pointed at @a n_it from the graph and delete all connecting edges
   * @param n_it Iterator pointing to node that needs to be removed
   * @return number a new valid iterator that points to the next node
   *
   * @pre n_it is a valid iterator and != graph.node_end()
   * @post new graph.num_nodes() = old graph.num_nodes-1
   * @post new graph.num_edges() = old.graph.num_edges - @a n.degree()
   * @post invalidates all node iterators in the range [0, old graph.num_nodes]
   *        and edge iterators in the range [0, old graph.num_edges]
   * @post Invalidates and edges that are connected to @a n
   *
   * This has time complexity  of < O(num_nodes()) assuming sparse graph
   */
  node_iterator remove_node(node_iterator n_it){
      assert(n_it != node_end());
      remove_node(*n_it);
      return n_it;
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
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return (e.node1()==n1 && e.node2() == n2) || (e.node2()==n1 && e.node1() == n2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {


      size_type aMin = std::min(this->node1().index(),this->node2().index());
      size_type aMax = std::max(this->node1().index(),this->node2().index());
      size_type bMin = std::min(e.node1().index(),e.node2().index());
      size_type bMax = std::max(e.node1().index(),e.node2().index());

      if (n1.parentGraph != e.node1().parentGraph){
          return n1.parentGraph < e.node1().parentGraph;
      }
      if (aMin == bMin){
          return aMax<bMax;
      }
      else{
          return aMin<bMin;
      }
    }

    /** return the length of the edge */
    double length(){
      return norm(n1.position() - n2.position());
    }

    /** return a modifiable value for this edge */
    edge_value_type& value(){
      size_type index;
      for (size_type i = 0; i<n1.degree(); ++i){
        if (n1.parentGraph->adjacency_list[n1.index()][i].node2_index == n2.index()){
          index = i;
        }
      }
      return n1.parentGraph->adjacency_list[n1.index()][index].edge_value;
    }

    /** return a constant value for this edge */
    const edge_value_type& value() const{
      size_type index;
      for (size_type i = 0; i<n1.degree(); ++i){
        if (n1.parentGraph->adjacency_list[n1.index()][i].node2_index == n2.index()){
          index = i;
        }
      }
      return n1.parentGraph->adjacency_list[n1.index()][index].edge_value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    node_type n1;
    node_type n2;
    Edge(node_type node1, node_type node2) : n1(node1), n2(node2){

    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return nEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {

    assert(i < num_edges());
    auto it = edge_begin();
    for(size_type ind = 0; ind < i; ind++ ){
        ++it;
    }
    return *it;        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    for (size_type i =0; i < a.degree(); i++){
      if (adjacency_list[a.index()][i].node2_index == b.index()){
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
      assert(a.index()!=b.index());
    edge_type newEdge (a,b);
    if (this->has_edge(a,b)){
        for (size_type i =0; i < a.degree(); i++){
          if (adjacency_list[a.index()][i].node2_index == b.index()){
              return newEdge;
          }
        }
    }
    edge_value_type edgeValue = edge_value_type();
    adjacencyInfo n1 {b.index(), edgeValue};
    adjacencyInfo n2 {a.index(), edgeValue};
    adjacency_list[a.index()].push_back(n1);
    adjacency_list[b.index()].push_back(n2);
    nEdges++;
    return newEdge;

  }

  /**
   * @brief Remove edge between nodes @a a and @a b
   * @param a First node of edge that needs to be removed
   * @param b Second node of edge that needs to be removed
   * @return number of edges removed
   *
   * @pre Nodes @a a and b are valid nodes in this graph
   * @post If graph.has_edge(a,b), new graph.num_edges() = old graph.num_edges-1
   * @post If graph.has_edge(a,b), invalidates all edge iterators in the range [0, old graph.num_edges]
   * @post Else no change to graph
   *
   * It has time complexity O(a.degree() + b.degree())
   */
  size_type remove_edge(const Node& a, const Node& b){
      assert(has_node(a));
      assert(has_node(b));
      if (has_edge(a,b)){
          //std::cout<<"Is it getting here?"<<std::endl;
          size_type aInd = a.index();
          size_type bInd = b.index();
          for(size_type i = 0; i < adjacency_list[aInd].size(); ++i){
              if (adjacency_list[aInd][i].node2_index == bInd){
                  adjacency_list[aInd][i] = adjacency_list[aInd].back();
                  adjacency_list[aInd].pop_back();
              }
          }

          for(size_type i = 0; i < adjacency_list[bInd].size(); ++i){
              if (adjacency_list[bInd][i].node2_index == aInd){
                  adjacency_list[bInd][i] = adjacency_list[bInd].back();
                  adjacency_list[bInd].pop_back();
              }
          }
          --nEdges;
          return 1;
      }
      return 0;
  }

  /**
   * @brief Remove edge @a e
   * @param e Edge that needs to be removed
   * @return number of edges removed
   *
   * @pre @a e is a valid edge
   * @post If graph contains @ e, new graph.num_edges() = old graph.num_edges-1
   * @post If graph contains @ e, invalidates all edge iterators in the range [0, old graph.num_edges]
   * @post Else no change to graph
   *
   * It has time complexity O(e.node1().degree() + e.node2().degree())
   */
  size_type remove_edge(const Edge& e){
      return remove_edge(e.node1(), e.node2());
  }

  /**
   * @brief Remove edge pointed at by @a e_it
   * @param e_it iterator that points to edge to be removed
   * @return edge iterator that points to the next edge
   *
   * @pre @a e_it != graph.edge_end()
   * @post If graph contains @ *e_it, new graph.num_edges() = old graph.num_edges-1
   * @post If graph contains @ *e_it, invalidates all edge iterators in the range [0, old graph.num_edges]
   * @post Else no change to graph
   *
   * It has time complexity O((*e_it).node1().degree() + (*e_it).node2().degree())
   */
  edge_iterator remove_edge(edge_iterator e_it){
      remove_edge(*e_it);
      return e_it;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nEdges=0;
    nNodes=0;
    adjacency_list.clear();
    node_list.clear();
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

    /**
     * @brief operator * used to dereference the node in the list being pointed too
     * @return returns the node being pointed too
     */
    Node operator*() const{
      return parentGraph->node(ind);
    }

    /**
     * @brief operator ++ used to increment the iterator by 1
     * @return same iterator with a index incremented by 1, unless the iterator is at
     * the end of the list
     */
    NodeIterator& operator++(){
        if (ind<parentGraph->size()){
            ++ind;
        }
        return *this;
    }

    /**
     * @brief operator == used to check the equality between two iterators
     * @param test: the iterator to be compared too
     * @return result of comparison
     *
     * Two iteratrso are equal if they point to the same graph and to the same index
     */
    bool operator==(const NodeIterator& test) const{
        return (test.parentGraph == parentGraph) && (test.ind == ind);
    }

   private:
    friend class Graph;

    graph_type* parentGraph;
    size_type ind;
    /**
     * @brief NodeIterator private constructor for the node_itreator
     * @param g pointer to graph that initialized this iterator
     * @param i index of the node this iterator is currenlty pointing too
     */
    NodeIterator(const Graph* g, size_type i) : parentGraph(const_cast<Graph*> (g)), ind(i) {}

  };

  // Supply definitions AND SPECIFICATIONS for:
  /**
   * @brief node_begin returns an iterator that points to beginning of the list of nodes
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /**
   * @brief node_end returns an iterator that points to end of the list of nodes
   */
  node_iterator node_end() const{
    return NodeIterator(this,this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<incident_iterator>{
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

    /**
     * @brief operator * used to dereference the connecting edge that this iterator
     * is pointing too
     * @return The connecting edge being pointed to by this iterator
     */
    Edge operator*() const{
      size_type node2_ind = parentGraph->adjacency_list[node_ind][ind].node2_index;
      //return parentGraph->edge_list[edge_ind];
//      return Edge(parentGraph->node_list[node_ind].node,
//                  parentGraph->node_list[node2_ind].node);
      return Edge(node_type(parentGraph, node_ind), node_type(parentGraph, node2_ind));
    }

    /**
     * @brief operator ++ used to increment the iterator by 1
     * @return same iterator with an ind = ind+1;
     */
    IncidentIterator& operator++(){
      if (ind < parentGraph->adjacency_list[node_ind].size()){
          ind = ind+1;
      }
      return *this;
    }

    /**
     * @brief operator == checks the equality of two iterators
     * @param test The iterator to be checked against
     * @return te result of equality comparison
     *
     * These edges enforce a strict weak ordering
     */
    bool operator==(const IncidentIterator& test) const{
      return (parentGraph == test.parentGraph)
            && (node_ind == test.node_ind)
            && (ind == test.ind);
    }

   private:
    friend class Graph;
    graph_type* parentGraph;
    size_type node_ind;
    size_type ind;

    /**
     * @brief IncidentIterator private constructor to iterate over all the edge connected
     * to this node
     * @param g pointer to graph that constructed this iterator
     * @param ni index of the node whose edges we are iterating over
     * @param i refers to the @a i'th edge connecting to this node
     */
    IncidentIterator(const Graph* g, size_type ni, size_type i)
        : parentGraph(const_cast<Graph*> (g)), node_ind (ni), ind(i) {}


    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<edge_iterator>{
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

    /**
     * @brief operator * Used to derefernce the edge being pointed too
     * @return The edge that th eiterator is currently pointing to
     */
    Edge operator*() const{
      return edge_type(parentGraph->node(ind1),parentGraph->node(parentGraph->adjacency_list[ind1][ind2].node2_index));
    }

    /**
     * @brief operator ++ Used to increment the iterator by 1
     * @return same iterator with ind = ind+1
     */
    EdgeIterator& operator++(){
      ++ind2;

      // find next edge by going through adjacency info
      while (ind1<parentGraph->adjacency_list.size()){

          while (ind2 < parentGraph->adjacency_list[ind1].size()){
              if (ind1 < parentGraph->adjacency_list[ind1][ind2].node2_index){
                  return *this;
              }
              ++ind2;
          }
          ++ind1;
          ind2 = 0;
      }

      return *this;
    }

    /**
     * @brief operator == Used to check the equality between two iterators
     * @param test The iterator that is beinng compared too
     * @return result of equality comparison
     *
     * Ordering of edges enforces a strict weak ordering
     */
    bool operator==(const EdgeIterator& test) const{
        return (test.parentGraph == parentGraph) && (test.ind1 == ind1) && (test.ind2 ==ind2);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* parentGraph;
    size_type ind1;
    size_type ind2;

    /**
     * @brief EdgeIterator private constructor to iterate over all the edge

     * @param g pointer to graph that constructed this iterator
     * @param i refers to the @a i'th edge of the graph
     */
    EdgeIterator(const Graph* g, size_type n1i, size_type n2i) : parentGraph(const_cast<Graph*> (g)), ind1(n1i), ind2(n2i) {}

  };

  /** Return an iterator that points to the first edge in the graph
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  /** Return an iteratort hat points to the last edge in the graph
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this,adjacency_list.size(), 0);
  }

 private:

  size_type nNodes; // Numbder of nodes in the graph
  size_type nEdges; // Number of edges in the graph

  /**
   * @brief The internal_node struct contains a node's position and value
   */
  struct internal_node {
    Point position;
    node_value_type value;
  };

  /**
   * @brief The adjacencyInfo struct contains the connectivity information of the nodes
   * Each struct contains the edge index of the connecting edge, and the node index of
   * the connected node
   */
  struct adjacencyInfo{
    size_type node2_index;
    edge_value_type edge_value;
  };

  std::vector<internal_node> node_list; // contains all the nodes
  std::vector<std::vector<adjacencyInfo>> adjacency_list; // contains connectivity info


  /** Get the position of node[i]
   */
  const Point& nodePosition(size_type i) const {
    return node_list[i].position;
  }

};

#endif // CME212_GRAPH_HPP
