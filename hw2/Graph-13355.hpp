
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
*/

#include <algorithm>
#include <vector>
#include <cassert>
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
 private:

  // Predeclaring the interal structures. 
  struct internal_node; 
  struct interal_edge;

 public:
  /**Synonym for node type value*/
  using node_value_type = V;

  /**Synonym for node type value*/
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
  Graph() 
      : nodes(), edges(), index2uid(), edges2uid(){  
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
    const Point& position() const {
     return (graph->nodes[this-> node_id_]).position; 
    } 

    /** Return this node's position. */
    Point& position() {
      return (graph->nodes[this-> node_id_]).position; 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
     return(this->graph->nodes[this->node_id_].node_index);
    }

    /** Return this node's value. */
    node_value_type& value() {
      return (this ->graph -> nodes[this -> node_id_].node_value);
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return (this ->graph -> nodes[this -> node_id_].node_value);
    }

    /** Return this node's degree, the number of edges connected to it. */
    size_type degree() const {
      return(this-> graph -> nodes[this -> node_id_].adjacent_edges.size());
    }

    /** Return where to start iterating the adjacent edges to a node. */
    incident_iterator edge_begin() const {
        return (IncidentIterator(this -> graph, this-> node_id_)); 
    }

    /** Return where to end iterating the adjacent edges to a node. */
    incident_iterator edge_end() const {
      return (IncidentIterator(this -> graph, this -> node_id_, this -> degree()));  
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    
    bool operator==(const Node& n) const {
      return (n.graph == this -> graph && n.node_id_ == this -> node_id_);
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
      if (n.graph == this -> graph) {
      return (node_id_ < n.node_id_);
      }
      else {
        return (this -> graph < n.graph); 
      }
    }
    

  private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Allow Graph to access Node's private member data and functions.
    friend class Graph; 
    Graph* graph;
    size_type node_id_;

    // Private Node constructor.
    Node(const Graph* g, size_type id)
      : graph(const_cast<Graph*>(g)), node_id_(id) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index2uid.size();
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

 Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_node;
    new_node.position = position; 
    new_node.node_value = value;
    new_node.node_index = this-> index2uid.size();
    nodes.push_back(new_node); // Adds node to internal node struct.
    index2uid.push_back(nodes.size()-1); // Adds index to i2u. 
    return Node(this,nodes.size()-1); 
  }


  /** Remove a node to the graph. 
   * This removes the nodes index from the i2u vector but does not
   * change the interal struct of nodes. This also calls the
   * remove_edge function in order to remove all adjacent edges.  
   * @param[in] n The node to be removed. 
   * @post g.node(i).index() = i for i < g.num_nodes()
   * @post g.node(n.index()) == n
   * @post index2uid.size() = num_nodes 
   * 
   * Complexity: O(degree + degree). (degree is the number
   * of edges connected to a node)
   */
  size_type remove_node(const Node& n) {
    //Checks to see if the graph has the node. 
    if(!has_node(n)) {
      return (0);
    }
    std::vector<Edge> edges_remove;
    
    //Removes the adjacent edges to each node. 
    for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge def_ei = *ei;
      edges_remove.push_back(def_ei);
    }
    for(size_type i = 0; i < edges_remove.size(); i++) {
      remove_edge(edges_remove[i]);
    }
    //Removes the node from the internal structure and the i2u list. 
    index2uid[n.index()] = index2uid.back();
    index2uid.pop_back();
    nodes[index2uid[n.index()]].node_index = n.index();
    return 1; 
  }

    /** Remove a node to the graph. 
   * This removes the nodes index from the i2u vector but does not
   * change the interal struct of nodes. This also calls the
   * remove_edge function in order to remove all adjacent edges.  
   * @param[in] n_it The iterator to the node to be removed. 
   * @post g.node(i).index() = i for i < g.num_nodes()
   * @post g.node(n.index()) == n
   * @post index2uid.size() = num_nodes 
   * 
   * Complexity: O(degree + degree). 
   */
   node_iterator remove_node(node_iterator n_it) {
    Node def_n_it = *n_it; 
    this->remove_node(def_n_it);
    return(n_it);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return(n.graph == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(this -> index2uid.size() <= i) {
      return Node();
    } else {
      return Node(this,index2uid[i]);
    } 
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
      return(Node(this->graph, this-> first_node_id));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return(Node(this->graph, this->second_node_id));
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.edge_id_ == this->edge_id_ && this->graph == e.graph);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    
    bool operator<(const Edge& e) const {
      if(this->graph == e.graph) {
        return(this->edge_id_ < e.edge_id_);
      } else {
        if(this->edge_id_ != e.edge_id_) {
          return(this->edge_id_ < e.edge_id_);
        } else {
          return(true);
        }
      }
    }

    /**Returns the length between two nodes using the norm*/
    double length() const{
        return(norm(this->graph->nodes[first_node_id].position -this->graph->nodes[second_node_id].position));
    }

    /** Return this edges value. */
    edge_value_type& value() {
       return (this ->graph -> edges[this -> edge_id_].edge_value);
    }

    /** Return this edges value. */
    const edge_value_type& value() const {
       return (this ->graph -> edges[this -> edge_id_].edge_value);
    }

   private:
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type edge_id_;
    size_type first_node_id;
    size_type second_node_id;
    // Private Edge constructor.
    Edge(const Graph* g, size_type id, size_type first_node, size_type second_node)
      : graph(const_cast<Graph*>(g)), edge_id_(id), first_node_id(first_node), 
      second_node_id(second_node) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges2uid.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges2uid[i], this -> edges[edges2uid[i]].first_node, 
      this -> edges[edges2uid[i]].second_node);
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
        if((e.node1() == a && e.node2() == b) ||(e.node1() == b && e.node2() == a)) {
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

 Edge add_edge(const Node& a, const Node& b,const edge_value_type value = edge_value_type()) {
   // Checking if the current edge already exists.
   if (has_edge(a,b)) {
      for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
      Edge e = *ei; 
        if((e.node1() == a && e.node2() == b) ||(e.node1() == b && e.node2() == a)) {
          return Edge(e); 
        }
      }
    } 
   // If distinct and does not exist, then create a new node and add it.
      internal_edge new_edge;
      new_edge.first_node = a.index(); 
      new_edge.second_node = b.index();
      new_edge.edge_id = this -> edges2uid.size();
      new_edge.edge_value = value;
      this -> edges.push_back(new_edge);
      this -> edges2uid.push_back(edges.size()-1); // Adds edge index to e2u. 
      this -> nodes[index2uid[a.index()]].adjacent_edges.push_back(edges2uid.size()-1); 
      this -> nodes[index2uid[b.index()]].adjacent_edges.push_back(edges2uid.size()-1); 
      return Edge(this,new_edge.edge_id, new_edge.first_node,
      new_edge.second_node);
    }


  /** Remove an edge between two nodes.  
   * This removes the edges index from the e2u vector but does not
   * change the interal struct of edges.
   * @param[in] a, b The nodes of that the edge is connected to. 
   * @post edges2uid.size() = num_edges 
   * 
   * Complexity: No more than O(num_edges + degree + degree). 
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Checks to see if the edge exists. 
    if(!has_edge(a, b)) return(0);

    //Finds the edge id of of the edge connected to a and b. 
    size_type id;
    for(auto it = edge_begin(); it != edge_end(); ++it) {
      Edge e = *it;
      if((e.node1() == a && e.node2() == b) || 
         (e.node1() == b && e.node2() == a)) {
          id = e.edge_id_;
          break;
      }
    }

    // Removes the edge id from the adjacency lists of node a and b. 
    for (size_type i =0; i < this-> nodes[index2uid[a.index()]].adjacent_edges.size(); i++) {
        if (this-> nodes[index2uid[a.index()]].adjacent_edges[i] == id) {
          nodes[index2uid[a.index()]].adjacent_edges.erase(this-> nodes[index2uid[a.index()]].adjacent_edges.begin() + i); 
        }
    }
    for (size_type i =0; i < this-> nodes[index2uid[b.index()]].adjacent_edges.size(); i++) {
        if (this-> nodes[index2uid[b.index()]].adjacent_edges[i] == id) {
          nodes[index2uid[b.index()]].adjacent_edges.erase(this-> nodes[index2uid[b.index()]].adjacent_edges.begin() + i); 
        }
    }

    // Removes the edge index from the e2u vector. 
    edges2uid[id] = edges2uid.back();
    edges2uid.pop_back();
    this -> edges[edges2uid[id]].edge_id = id;
    return(1);
  }


  /** Remove an edge. 
   * This removes the edges index from the e2u vector but does not
   * change the interal struct of edges.
   * @param[in] e The edge to be removed.
   * @post edges2uid.size() = num_edges 
   * 
   * Complexity: No more than O(num_edges + degree + degree). 
   */
  size_type remove_edge(const Edge& e) {
    return(remove_edge(e.node1(), e.node2()));

  }

  /** Remove an edge. 
   * This removes the edges index from the e2u vector but does not
   * change the interal struct of edges.
   * @param[in] e_it The iterator of the edge to be removed. 
   * @post edges2uid.size() = num_edges 
   * 
   * Complexity: No more than O(num_edges + degree + degree). 
   */
  edge_iterator remove_edge (edge_iterator e_it) {
    Edge def_e_it = *e_it;
    return (EdgeIterator(this, remove_edge(def_e_it)));
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this -> index2uid.clear();
    this -> edges2uid.clear();
    this -> nodes.clear();
    this -> edges.clear();
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

    /** Return a dereferenced version of the node iterator. */
    Node operator*() const {
      return (Node(this -> iterator_graph, this -> graph_position));
    }

    /** Return the node iterator after incrementing it. */
    NodeIterator& operator++() {
      this -> graph_position++; 
      return *this; 
      
    }

    /** Tests whether this NodeIterator and @n are the same. */
    bool operator==(const NodeIterator& n) const {
      return (this -> graph_position == n.graph_position && this-> iterator_graph 
        == n.iterator_graph); 
    }


   private:
    // Private data members and methods for NodeIterator that will not be visible to users, 
    // but may be useful within Graph.
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* iterator_graph;
    size_type graph_position; 

    // Private NodeIterator constructor.
    NodeIterator(const Graph* g, const size_type position = 0)
      : iterator_graph(const_cast<Graph*>(g)), graph_position(position){
    }
  };

  /**Returns where to start the node iterator.*/
  node_iterator node_begin() const {
    return node_iterator(this); 
  }

  /**Returns where to end the node iterator.*/
  node_iterator node_end() const {
    return node_iterator(this, this -> nodes.size()); 
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


    /** Return a dereferenced version of the incident iterator. */
    Edge operator*() const {
      size_type edge_id_start = this -> iterator_graph -> 
      nodes[this ->iterator_graph -> index2uid[node_start_id]].adjacent_edges[this -> graph_position]; 
      size_type adjecent_node_id = this -> iterator_graph -> edges[this ->iterator_graph -> 
        edges2uid[edge_id_start]].second_node;
      if (adjecent_node_id == node_start_id) {
        adjecent_node_id = this -> iterator_graph -> edges[this -> iterator_graph 
          -> edges2uid[edge_id_start]].first_node; 
      }

      return (Edge(this -> iterator_graph, edge_id_start, this -> node_start_id, adjecent_node_id)); 
    }

    /** Return the incident iterator after incrementing it. */
    IncidentIterator& operator++() {
      this -> graph_position++; 
      return (*this);
    }

    /** Tests whether this IncidentIterator and @n are the same. */
    bool operator==(const IncidentIterator& n) const {
      return (this -> graph_position == n.graph_position && this-> iterator_graph 
        == n.iterator_graph && this -> node_start_id == n.node_start_id); 
    }


   private:
    // Private data members and methods for IncidentIterator that will not be visible to users, 
    // but may be useful within Graph.
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* iterator_graph;
    size_type node_start_id; 
    size_type graph_position; 

    // Private NodeIterator constructor.
    IncidentIterator(const Graph* g, const size_type node_given, const size_type position = 0)
      : iterator_graph(const_cast<Graph*>(g)), node_start_id(node_given), graph_position(position){
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

    /** Return a dereferenced version of the edge iterator. */
  Edge operator*() const{
      size_type node_one = this -> iterator_graph -> edges[this -> iterator_graph 
        -> edges2uid[graph_position]].first_node; 
      size_type node_two = this -> iterator_graph -> edges[this -> iterator_graph 
        -> edges2uid[graph_position]].second_node;
      if (node_one == node_two){
        node_two = this -> iterator_graph -> edges[this -> iterator_graph -> edges2uid[graph_position]].first_node; 
      }
      return (Edge(this -> iterator_graph, this -> graph_position, node_one, node_two));
    }

    /** Return the edge iterator after incrementing it. */
    EdgeIterator& operator++() {
      this -> graph_position++; 
      return *this; 
    }

    /** Tests whether this EdgeIterator and @e are the same. */
    bool operator==(const EdgeIterator& e) const {
      return (this -> graph_position == e.graph_position && this-> iterator_graph 
      == e.iterator_graph); 
    }

   private:
    // Private data members and methods for EdgeIterator that will not be visible to users, 
    // but may be useful within Graph.
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* iterator_graph;
    size_type graph_position; 

    //Private EdgeIterator constructor
    EdgeIterator(const Graph* g, const size_type position = 0) 
      : iterator_graph(const_cast<Graph*>(g)), graph_position(position){
    }
  };

  /**Returns where to start the edge iterator.*/
  edge_iterator edge_begin() const {
   return  edge_iterator(this); 
  }

  /**Returns where to end the edge iterator.*/
  edge_iterator edge_end() const {
    return edge_iterator(this, this -> edges2uid.size()); 
  }

 private:
  // Internal type node and edge.
  struct internal_node {
    Point position; // Position of the node.
    V node_value; // The identification for the node type.
    size_type node_index; //The identification for the node. 
    std::vector<size_type> adjacent_edges; // All of the incident edges. 
  };

  struct internal_edge {
  size_type first_node; // First node of the edge.
  size_type second_node;  // Second node of the edge.
  size_type edge_id; // The identification for the edge. Each is unique.
  E edge_value; // The identiication for the edge value. 
  };

  // Initializing the vectors of the interal structs. 
  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;

  //Initializign the vectors that allow removal of edges and nodes.
  std::vector<size_type> index2uid; 
  std::vector<size_type> edges2uid; 
};

#endif // CME212_GRAPH_HPP
