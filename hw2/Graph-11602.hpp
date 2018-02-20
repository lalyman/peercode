#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <cmath>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
	typedef V node_value_type;
 	typedef E edge_value_type;
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // Adding template usage
  //using node_value_type = V;

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
  Graph() 
    : nodes(), edges() {
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
    }

    /** Return this node's position. */
    Point& position() {
      return graphPtr->nodes[idx].xy;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->idx;
    }

    // Set the value type
    const node_value_type& value() const{
    	return graphPtr->nodes[idx].value;
    }
    node_value_type& value(){
      return graphPtr->nodes[idx].value;
    }

    // Find the number of edges incident to node
    size_type degree() const{
      return (graphPtr->connections.at(idx).size());
    }

    // Create the iterators for the incident edges
    incident_iterator edge_begin() const{
      return IncidentIterator(graphPtr,graphPtr->connections[this->index()].begin());
    }
    incident_iterator edge_end() const{
      return IncidentIterator(graphPtr,graphPtr->connections[this->index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (idx == n.idx && graphPtr == n.graphPtr){
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
      // Make sure they're in the same graph, then use < idx
      if (graphPtr != n.graphPtr){
        return false;
      }
      return (idx < n.idx);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Add special privat members for proxy to graph
    Graph* graphPtr;
    size_type idx;
    Node(const Graph* graph, size_type index)
        : graphPtr(const_cast<Graph*>(graph)), idx(index){
        }


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
  Node add_node(const Point& position) {
    
    // Create the nodes struct
    node_single newNode;
    newNode.xy = position;
    newNode.idx = nodes.size()+1;
    nodes.push_back(newNode);

    // Add node to the connections (with empty connection)
    connections.insert(std::pair<size_type,std::map<size_type,size_type>>(newNode.idx,{}));

    return Node(this,newNode.idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // See if the node is a part of this graph
    if (this != n.graphPtr){
      return false;
    }

    // If it is, see if the index is within range
    return (nodes.size() < n.idx);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,nodes[i].idx);       
  }

  /** Remove a node from the graph.
   * @pre @a a is a distinct valid node of this graph
   * @return an index to the next node
   * @post nodes now have one less index
   * @post edges incident to node are now gone
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_node(const Node& n){
  	nodes.erase(nodes.begin()+n.idx);

  	//Remove the edges too
  	for (a=n.edge_begin(); a != n.edge_end(); a++)
  		edges.erase(connections[((*a).idx)]);
  	return n.idx;
  };
  node_iterator remove_node(node_iterator n_it){
  	nodes.erase(n_it);
  	return n_it;
  };


  /** Remove an edge from the graph.
   * @pre @a a is a distinct valid edge of this graph
   * @return an index to the next edge
   * @post edges now have one less index
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node& n1, const Node& n2){
  	auto edgeidx = connections[n1.idx][n2.idx];
  	edges.erase(edges.begin()+edgeidx);
  	return edgeidx;
  };
  size_type remove_edge(const Edge& e){
  	edges.erase(edges.begin()+e.idx);
  	return e.idx;
  };
  edge_iterator remove_edge(edge_iterator e_it){
  	edges.erase(e_it);
  	return e_it;
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

    edge_value_type& value(){
  		return graphPtr->edges[this->idx].value;
  	}
  	const edge_value_type& value() const{
  		return graphPtr->edges[this->idx].value;
  	}

    /** Return a node of this Edge */
    Node node1() const {
      return graphPtr->edges[idx].a;     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graphPtr->edges[idx].b;    
    }

    // Return the length of the edge
    double length() const {
    	auto n1 = this->node1().position();
    	auto n2 = this->node2().position();
    	return std::sqrt(std::pow(n1.x-n2.x,2) + std::pow(n1.y-n2.y,2) + std::pow(n1.z-n2.z,2));
    }

    /** Test whether this edge and @a


     e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (idx < e.idx);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (idx < e.idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Add needed variables
    Graph* graphPtr;
    size_type idx;
    Edge(const Graph* graph, size_type index)
        : graphPtr(const_cast<Graph*>(graph)), idx(index){
        }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this,i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    //Search through the edges for a connection
    if (connections.at(a.idx).count(b.idx) == 1){
      return true;
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
    
    // Search through the edges and return if one exists
    if (has_edge(a,b) == true){
      return Edge(this,connections.at(a.idx).at(b.idx));
    }

    // Otherwise make a new edge
    edge_single newEdge;
    newEdge.a = a;
    newEdge.b = b;
    newEdge.idx = edges.size();
    edges.push_back(newEdge);
    
    // Add the edge to the connections
    connections[a.idx].insert(std::pair<size_type,size_type>(b.idx,newEdge.idx));
    connections[b.idx].insert(std::pair<size_type,size_type>(a.idx,newEdge.idx));

    return Edge(this,newEdge.idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    connections.clear();
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

    // Return a Node when dereferencing
    Node operator*() const{
      return Node(this->graphPtr,this->idx);
    };

    // Move to the next node when iterating
    NodeIterator& operator++(){
      this->idx++;
      return *this;
    };

    // They are equal if the index is equal (assumed they are from same graph)
    bool operator==(const NodeIterator& n) const{
      return (this->idx == n.idx);
    };

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graphPtr;
    size_type idx;
    NodeIterator(const Graph* graph, size_type index)
        : graphPtr(const_cast<Graph*>(graph)), idx(index){
        }
  };

  // Set the start and end point based on size of the vector
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  };
  node_iterator node_end() const{
    return NodeIterator(this,nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    // Return an edge when dereferencing
    Edge operator*() const{
      return Edge(graphPtr,idxMap->second);
    };

    // Index the map via built-in iterator - return result
    IncidentIterator& operator++(){
      ++idxMap;
      return *this;
    };

    // Compare index (assumed theyre from the same graph)
    bool operator==(const IncidentIterator& n) const{
      return (this->idxMap->second == n.idxMap->second);
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graphPtr;
    std::map<size_type, size_type>::iterator idxMap;
    IncidentIterator(const Graph* graph, std::map<size_type, size_type>::iterator index)
        : graphPtr(const_cast<Graph*>(graph)), idxMap(index){
        }
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

    // Return Edge when dereferencing
    Edge operator*() const{
      return Edge(graphPtr,idx);
    };

    // Index to next element of vector
    EdgeIterator& operator++(){
      this->idx++;
      return *this;
    };

    // Compare index number (assume from same graph)
    bool operator==(const EdgeIterator& n) const{
      return (this->idx == n.idx);
    };


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graphPtr;
    size_type idx;
    EdgeIterator(const Graph* graph, size_type index)
        : graphPtr(const_cast<Graph*>(graph)), idx(index){
        }
  };

  // Set start and end point of iter based on size of edge vector
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  };
  edge_iterator edge_end() const{
    return EdgeIterator(this,edges.size());
  };

 private:

  // Define the node and edge structs
  struct node_single
  {
    Point xy;
    size_type idx;
    node_value_type value;
  };
  struct edge_single
  {
    Node a;
    Node b;
    size_type idx;
    edge_value_type value;
  };

  // Create the vectors to hold the nodes and edges
  std::vector<node_single> nodes;
  std::vector<edge_single> edges;

  // Create a map to hold the connections and idx
  std::map<size_type, std::map <size_type, size_type>> connections;

};

#endif // CME212_GRAPH_HPP
