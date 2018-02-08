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
      if (id_ < graph_->nodes_vec.size())
        return graph_->nodes_vec[id_];
      assert(false);
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
        if (id_ < graph_->nodes_vec.size())
            return graph_->value_vec[id_];
        assert(false);
    }

    /** Return the node value */
    const node_value_type& value() const{
        if (id_ < graph_->nodes_vec.size())
            return graph_->value_vec[id_];
        assert(false);
    }

    /** Return the degree of a node */
    size_type degree() const{
        return graph_->e_map[id_].size();
    }

    /** Return the beginning incident iterator */
    incident_iterator edge_begin() const{
        return IncidentIterator(graph_, id_, size_type(0));
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
    node_value_type value_;

    // Private constructor
    Node(const Graph* g, size_type id, 
    	node_value_type value = node_value_type())
      : graph_(const_cast<Graph*>(g)), id_(id), value_(value){
    }
   };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
    // HW0: YOUR CODE HERE
    // Add one element (position) to the vector
    // of nodes. Also, create a new key in the
    // map of edges. 
    nodes_vec.push_back(position);
    value_vec.push_back(value);
    e_map[nodes_vec.size()-1];
    return Node(this, nodes_vec.size()-1, value);     
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
      // HW0: YOUR CODE HERE
      if (edge_id < graph_->edge_vec.size())
      	  return graph_->node(node_1_id);
      assert(false); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (edge_id < graph_->edge_vec.size())
      	  return graph_->node(node_2_id);
      assert(false);   
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
    return edge_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return Edge(this, i, edge_vec[i][0], edge_vec[i][1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Keys of the edges' map are nodes' ids and values associated
    // with each key are the edges' ids. Hence, we loop over every
    // edge id related to node a. If b is in the value vector,
    // then there exists an edge between a and b.
    if (has_node(a) && has_node(b)){
        std::vector<size_type> current = e_map.at(a.index());
        for (size_type i = 0; i < current.size(); ++i){
          Edge curr_e = edge(current[i]);
          if (((curr_e.node1() == a) && (curr_e.node2() == b)) || 
            ((curr_e.node1() == b) && (curr_e.node2() == a)))
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
    // First, check a and b are valid nodes and they are not
    // equal. Then, check if there exists an edge between a and b.
    // If it does, find out that edge by looping over the value vector
    // associated with the key a. Otherwise, create a new vector
    // containing the indices of a and b and add to the vector
    // of edges. Also update the map of edges.
    
    if (has_node(a) && has_node(b) && !(a == b)){
          std::vector<size_type> current = e_map.at(a.index());
          if (has_edge(a,b)){           
            for (size_type i = 0; i < current.size(); ++i){
                Edge curr_e = edge(current[i]);
                if ((curr_e.node1() == a && curr_e.node2() == b) || 
                  (curr_e.node1() == b && curr_e.node2() == a))
                    return (Edge(this, i, a.index(), b.index()));
            }
          }
         
        else{
        std::vector<size_type> new_edge;
        new_edge.push_back(a.index());
        new_edge.push_back(b.index());
        edge_vec.push_back(new_edge);
        e_map.at(a.index()).push_back(edge_vec.size()-1);
        e_map.at(b.index()).push_back(edge_vec.size()-1);
        return Edge(this, edge_vec.size()-1, a.index(), b.index());
        }     
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
    nodes_vec.clear();
    edge_vec.clear();
    e_map.clear();
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
        if(iter_id < graph_->size())
            return graph_->node(iter_id);
        assert(false);
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
        return NodeIterator(this, size_type(0));
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
    	size_type e_id = graph_->e_map[node_id][edge_id];
        if (graph_->edge_vec[e_id][0] == node_id){
            return Edge(graph_, edge_id, node_id, graph_->edge_vec[e_id][1]);
        }
        else{
            return Edge(graph_, edge_id, node_id, graph_->edge_vec[e_id][0]);
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
        if (edge_ < graph_->edge_vec.size()){
            return graph_->edge(edge_);
        }
        assert(false);
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
      return EdgeIterator(this, size_type(0));
  }

  /** Return the ending edge iterator */
  edge_iterator edge_end() const{
      return EdgeIterator(this, edge_vec.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.
  std::vector<Point> nodes_vec;
  std::vector<node_value_type> value_vec;
  std::vector<std::vector<size_type>> edge_vec;
  std::map<size_type, std::vector<size_type>> e_map;
};

#endif // CME212_GRAPH_HPP