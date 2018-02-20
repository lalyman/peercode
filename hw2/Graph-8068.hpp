#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

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
   
   struct internal_nodes;
   struct internal_edges;
   

 public:

  /** Type of this graph. */
  using graph_type = Graph;
  //using node_value_type = V;

  typedef V node_value_type;
  typedef E edge_value_type;

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

    /** HW2 #1: Return this node's position. */
    Point& position() {
      return graph_->nodes_vec[this->id_].pos;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_vec[this->id_].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the node value */
    node_value_type& value(){
        return graph_->nodes_vec[this->id_].node_value;
    }

    /** Return the node value */
    const node_value_type& value() const{
        return graph_->nodes_vec[this->id_].node_value;
    }

    /** Return the degree of a node */
    size_type degree() const{
        return graph_->nodes_vec[this->id_].incident_edge_index.size();
    }

    /** Return the beginning incident iterator */
    incident_iterator edge_begin() const{
        return IncidentIterator(graph_, id_, 0);
    }

    /** Return the ending incident iterator */
    incident_iterator edge_end() const{
        return IncidentIterator(graph_, id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((n.graph_ == graph_) && (n.id_ == id_));
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
      bool same_graph = ((graph_ == n.graph_) && (id_ < n.id_));
      return (same_graph || (graph_ < n.graph_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type id_;

    // Private constructor
    Node(const Graph* g, size_type id)
      : graph_(const_cast<Graph*>(g)), id_(id){
    }
   };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_vec.size();
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
     internal_nodes new_node;
       new_node.pos = position;
       new_node.node_value = value;
       new_node.node_id = size();
       nodes_vec.push_back(new_node);
       return Node(this, new_node.node_id);
    
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((n.graph_ == this) && (n.id_ < size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < size());
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,this->node_1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,this->node_2_id);  
    }

    /** Return the value of this edge */
    edge_value_type& value() {
        return graph_->edges_vec[this->edge_id].edge_value;
    }

    /** Return the value of this edge */
    const edge_value_type& value() const {
        return graph_->edges_vec[this->edge_id].edge_value;
    }

    double length() const {
      return norm(node1().position() - node2().position());
    } 

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((node1() == e.node1())&&(node2() == e.node2())) || 
        ((node1() == e.node2())&&(node2() == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      bool same_graph = ((graph_ == e.graph_) && (edge_id < e.edge_id));
      return ((graph_ < e.graph_) || same_graph);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type edge_id;
    size_type node_1_id;
    size_type node_2_id;

    // Private constructor
    Edge(const Graph* g, size_type id, size_type node_1_, size_type node_2_)
      : graph_(const_cast<Graph*>(g)), edge_id(id), node_1_id(node_1_) 
        ,node_2_id(node_2_){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return Edge(this, i, edges_vec[i].node_one, edges_vec[i].node_two);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool contains_edge;
    for(size_type i : nodes_vec[a.index()].incident_edge_index){
      contains_edge = ((edge(i).node1() == a && edge(i).node2() == b) ||
                      (edge(i).node1() == b && edge(i).node2() == a));
      if (contains_edge){
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
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& value = edge_value_type()) {
      if (has_edge(a,b)){
        for (size_type i = 0; i < edges_vec.size(); ++i){
          if ((edges_vec[i].node_one == a.index() && 
               edges_vec[i].node_two == b.index()) || 
              (edges_vec[i].node_one == b.index() && 
              edges_vec[i].node_two == a.index())){
            return Edge(this, i, a.index(), b.index()); 
          }
        }   
      }
      internal_edges new_edge; 
      new_edge.node_one = a.index(); 
      new_edge.node_two = b.index();
      new_edge.edge_value = value;
      edges_vec.push_back(new_edge);
      this->nodes_vec[a.index()].incident_edge_index.push_back(num_edges()-1);
      this->nodes_vec[b.index()].incident_edge_index.push_back(num_edges()-1); 
      
      return Edge(this, edges_vec.size()-1, new_edge.node_one, new_edge.node_two); 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_vec.clear();
    edges_vec.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator> {
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

    /** Return the node with a specified index */
    Node operator*() const{
        return Node(graph_, this->iter_id);
    }

    /** Return the next node iterator */
    NodeIterator& operator++(){
        iter_id += 1;
        return *this;
    }

    /** Compare if two node iterators are equal.
   * @param[in] n A node iterator
   * @return a boolean checking whether the two node iterators
   *  are in the same graph and their indices are the same
   *
   * Complexity: O(1)
   */
    bool operator==(const NodeIterator& n) const{
        return ((graph_ == n.graph_) && (iter_id == n.iter_id));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type iter_id;

    // Private constructor
    NodeIterator(const Graph* g, size_type id)
      : graph_(const_cast<Graph*>(g)), iter_id(id){
    }
   };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Return the beginning node iterator */
    node_iterator node_begin() const{
        return NodeIterator(this, 0);
    }

    /** Return the ending node iterator */
    node_iterator node_end() const{
        return NodeIterator(this, size());
    }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge that is connected to the node */
    Edge operator*() const{
      size_type e_id = graph_->nodes_vec[this->node_id].incident_edge_index[edge_id];
      if (graph_->edges_vec[e_id].node_one == node_id){
            return Edge(graph_, e_id, node_id, graph_->edges_vec[e_id].node_two);
        }
        else{
            return Edge(graph_, e_id, node_id, graph_->edges_vec[e_id].node_one);
        }
    }

    /** Return the next incident iterator */
    IncidentIterator& operator++(){
        edge_id += 1;
        return *this;
    }

    /** Compare if two incident iterators are equal.
   * @param[in] n An incident iterator
   * @return a boolean checking whether the two incident iterators
   *  are in the same graph, their node indices are the same, and
   *  their indices are the same
   *
   * Complexity: O(1)
   */
    bool operator==(const IncidentIterator& n) const{
        return ((graph_ == n.graph_) && (node_id == n.node_id) 
             && (edge_id == n.edge_id));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_id;
    size_type edge_id;

    // Private constructor
    IncidentIterator(const Graph* g, size_type node_id_, size_type edge_)
      : graph_(const_cast<Graph*>(g)), node_id(node_id_), edge_id(edge_){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge that is associated with the edge index */
    Edge operator*() const{
        internal_edges int_edge = graph_->edges_vec[this->edge_];
        return Edge(graph_, edge_, int_edge.node_one, int_edge.node_two);
    }

    /** Return the next edge iterator */
    EdgeIterator& operator++(){
        edge_+= 1;
        return *this;
    }

    /** Compare if two edge iterators are equal.
   * @param[in] n An edge iterator
   * @return a boolean checking whether the two edge iterators
   *  are in the same graph and their indices are the same
   *
   * Complexity: O(1)
   */
    bool operator==(const EdgeIterator& n) const{
        return ((graph_==n.graph_) && (edge_==n.edge_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_;

    // Private constructor
    EdgeIterator(const Graph* g, size_type edge_id)
      : graph_(const_cast<Graph*>(g)), edge_(edge_id){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the beginning edge iterator */
  edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
  }

  /** Return the ending edge iterator */
  edge_iterator edge_end() const{
      return EdgeIterator(this, edges_vec.size());
  }
  
  /**
   * @brief Remove an edge for two given nodes
   *
   * @param a A node
   * @param b A node
   * @return  Edge index of the edge removed or the number of edges
   *
   * @pre has_edge(a,b)
   * @post new_edges_vec.size() = old_edges_vec.size() - 1
   * @post nodes_vec[a.index()].incident_edge_index does not contain the
   *       index of the edge removed
   * @post nodes_vec[b.index()].incident_edge_index does not contain the
   *       index of the edge removed
   * @post new_edges_vec[edge] = old_edges_vec.back()
   *
   * The edge connecting nodes a and b is invalidated. 
   *
   * The complexity of remove_edge is O(num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b){
    //Delete the edge in edges_vec and incident_edge_index for a and b
    //Reorder

    /** Find the edge index, delete the edge in edges_vec,
        delete the edge in incident edges for a and b */
    if (has_edge(a,b)) {
      size_type index = 0;
      for (auto i = edges_vec.begin(); i != edges_vec.end(); ++i){
        if ((a.index()==(*i).node_one && b.index()==(*i).node_two) || 
            (a.index()==(*i).node_two && b.index()==(*i).node_one)){
          index = i - edges_vec.begin();
          break;      
        }
      }
      
      edges_vec[index] = edges_vec.back();
      edges_vec.pop_back();

      for (auto i = nodes_vec[a.index()].incident_edge_index.begin(); 
                i != nodes_vec[a.index()].incident_edge_index.end(); ++i){
        if ((*i)==index){
          (*i) = nodes_vec[a.index()].incident_edge_index.back();
          nodes_vec[a.index()].incident_edge_index.pop_back();
          break;
        }
      }

      for (auto i = nodes_vec[b.index()].incident_edge_index.begin(); 
                i != nodes_vec[b.index()].incident_edge_index.end(); ++i){
        if ((*i)==index){
          (*i) = nodes_vec[b.index()].incident_edge_index.back();
          nodes_vec[b.index()].incident_edge_index.pop_back();
          break;
        }
      }

      /** Update the edges in incident edges for the new
          edge in that position */
      if (index < num_edges()) {
        size_type new_one = edges_vec[index].node_one;
        size_type new_two = edges_vec[index].node_two;

        for (auto i = nodes_vec[new_one].incident_edge_index.begin();
                  i != nodes_vec[new_one].incident_edge_index.end(); ++i){
          if ((*i) == num_edges()){
            (*i) = index;
            break;
          }
        }

        for (auto i = nodes_vec[new_two].incident_edge_index.begin();
                  i != nodes_vec[new_two].incident_edge_index.end(); ++i){
          if ((*i) == num_edges()){
            (*i) = index;
            break;
          }
        }
      }
      return index;
    }
    else
      return num_edges();
  }

  /**
   * @brief Remove a given edge
   *
   * @param e An edge
   * @return  Edge index of the edge removed or the number of edges
   *
   * @pre has_edge(e.node1(), e.node2())
   * @post new_edges_vec.size() = old_edges_vec.size() - 1
   * @post nodes_vec[e.node1().index()].incident_edge_index does not 
   *       contain the index of the edge removed
   * @post nodes_vec[e.node2().index()].incident_edge_index does not 
   *       contain the index of the edge removed
   * @post new_edges_vec[edge] = old_edges_vec.back()
   *
   * The edge e is invalidated. 
   *
   * The complexity of remove_edge is O(num_edges())
   */
  size_type remove_edge(const Edge& e){
    size_type result = remove_edge(e.node1(), e.node2());
    return result;
  }

  /**
   * @brief Remove an edge given its iterator
   *
   * @param e_it An edge iterator
   * @return e_it
   *
   * @pre has_edge((*e_it).node1(), (*e_it).node2())
   * @post new_edges_vec.size() = old_edges_vec.size() - 1
   * @post nodes_vec[(*e_it).node1().index()].incident_edge_index does 
   *       not contain the index of the edge removed
   * @post nodes_vec[(*e_it).node2().index()].incident_edge_index does 
   *       not contain the index of the edge removed
   * @post new_edges_vec[edge] = old_edges_vec.back()
   *
   * The edge *e_it is invalidated. 
   *
   * The complexity of remove_edge is O(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }

  /**
   * @brief Remove a given node
   *
   * @param n A node
   * @return the node index or the number of nodes
   *
   * @pre has_node(n)
   * @post new_nodes_vec.size() = old_nodes_vec.size() - 1
   * @post old_nodes_vec[n.index()].incident_edge_index.size() = 0
   * @post new_nodes_vec[n.index()] = old_nodes_vec.back()
   *
   * The node n is invalidated and all incident edges associated with it
   * are invalidated.
   *
   * The complexity of remove_node is O(num_nodes())
   */
  size_type remove_node(const Node& n){
    if (!has_node(n)){
      return num_nodes();
    }

    /** Delete incident edges associated with n */ 
    while (nodes_vec[n.index()].incident_edge_index.size() != 0)
      remove_edge(edge(nodes_vec[n.index()].incident_edge_index.back()));

    /** Delete the node in nodes_vec */
    size_type old_n_idx = n.index(); 
    nodes_vec[old_n_idx] = nodes_vec.back();
    nodes_vec.pop_back();
    
    /** Update indices associated with the new node at n.index() */
    if (old_n_idx < num_nodes()){
      nodes_vec[old_n_idx].node_id = old_n_idx;
      for (auto i = nodes_vec[old_n_idx].incident_edge_index.begin(); 
              i != nodes_vec[old_n_idx].incident_edge_index.end(); ++i){
          if (edges_vec[*i].node_one == num_nodes()) {
              edges_vec[*i].node_one = old_n_idx;
          }
          if (edges_vec[*i].node_two == num_nodes()) {
              edges_vec[*i].node_two = old_n_idx;
          }
      }
    }
    return n.index();
  }

  /**
   * @brief Remove a node given its iterator
   *
   * @param n_it A node iterator
   * @return n_it
   *
   * @pre has_node(n)
   * @post new_nodes_vec.size() = old_nodes_vec.size() - 1
   * @post old_nodes_vec[(*n_it).index()].incident_edge_index.size() = 0
   * @post new_nodes_vec[(*n_it).index()] = old_nodes_vec.back()
   *
   * The node n is invalidated, all incident edges associated with it
   * are invalidated, and the node iterator n_it is invalidated.
   *
   * The complexity of remove_node is O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it){
    remove_node((*n_it));
    return n_it;
  }
  
 
 private:
  
  /** Internal type for nodes */
  /** pos - a Point object specifying the position
      node_id - index number for that node
      node_value - value of that node
      incident_edge_index - a vector containing edge indices
                            incident to that node
  */
  struct internal_nodes {
    Point pos;
    size_type node_id;
    node_value_type node_value;
    std::vector<size_type> incident_edge_index;
  };

  /** Internal type for edges */
  /** node_one - the index associated with the first node in that edge
      node_two - the index associated with the second node in that edge
      edge_value - value of that edge
  */
  struct internal_edges {
    size_type node_one;
    size_type node_two;
    edge_value_type edge_value;
  };

  /** vector containing internal nodes */
  std::vector<internal_nodes> nodes_vec;
  /** vector containing internal edges */
  std::vector<internal_edges> edges_vec; 
};

#endif // CME212_GRAPH_HPP