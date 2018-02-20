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
template <typename V, typename E>
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

  /** Type of node value */
  //using node_value_type = V;

  typedef V node_value_type;
  typedef E edge_value_type;

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
  Graph(): nodes(), edges(), node_edges() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() = default;

  /** Remove a graph's nodes
   * @param[in] n      Node that is to be removed
   * @return The size_type 1 indicating node has been removed
   *         The size_type 0 indicating node has not been removed
   *
   * @post new num_nodes() = old num_nodes() - 1
   * @post edges and adjacency edges that contain @a n.node_id_
   *       are removed.
   * @post adjacency edges are updated with new adjacent node position
   *
   * Complexity: No more than O(num_nodes()), should be less  
   */
  size_type remove_node(const Node& n) {

    // Check if node is available
    if (!has_node(n)) {
      return 0;
    }
    // Record down the id of the node that is to be removed and replaced
    size_type mid_node_ID = n.node_id_;
    size_type old_ori_node_ID = nodes.size()-1;

    // Remove all the adjacent edges of node n using remove_edge function
    while (!(node_edges[mid_node_ID].empty())) {
      auto n1 = node(mid_node_ID);
      auto n2 = node(node_edges[mid_node_ID].back().first);   
      remove_edge(n1,n2);   
    }

    // Remove the node n from the node_edges that conains all the adjacent nodes
    node_edges[mid_node_ID] = node_edges[old_ori_node_ID];
    node_edges.pop_back();

    // 1. Update all the edges with the new node ID
    // 2. Update adjacent node ID with the new node ID  
    for (auto it = node_edges[mid_node_ID].begin(); it != node_edges[mid_node_ID].end(); ++it) {
      // Purpose 1
      if (edges[(*it).second].first == old_ori_node_ID) {
        edges[(*it).second].first = mid_node_ID;
      } else if (edges[(*it).second].second == old_ori_node_ID) {
        edges[(*it).second].second = mid_node_ID; 
      }           
      for (auto it2 = node_edges[(*it).first].begin(); it2 != node_edges[(*it).first].end(); ++it2) {
        if ((*it2).first == old_ori_node_ID) {
          // Purpose 2
          (*it2).first = mid_node_ID;

          // Purpose 1
          if (edges[(*it2).second].first == old_ori_node_ID) {
            edges[(*it2).second].first = mid_node_ID;
          } else if (edges[(*it2).second].second == old_ori_node_ID) {
            edges[(*it2).second].second = mid_node_ID;            
          }
        }
      }
    }

    // Delete n from nodes and replace n with the last node
    nodes[mid_node_ID] = nodes[old_ori_node_ID];
    nodes.pop_back();
    return 1;
  }

  /** Remove a graph's nodes
   * @param[in] n_it    Node iterator that points to the node to be removed
   * @return The node iterator that points to the next node
   *
   * @pre 0 <= @a *it.node_id_ < @a old num_nodes()
   * @post new num_nodes() = old num_nodes() - 1
   * @post edges and adjacency edges that contain @a *it.node_id_
   *       are removed.
   * @post adjacency edges are updated with new adjacent node position
   *  
   * Complexity: No more than O(num_nodes()), should be less 
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** Remove a graph's edge
   * @param[in] n1    First node of the edge
   * @param[in] n2    Second node of the edge
   * @return The size_type 1 indicating node has been removed
   *         The size_type 0 indicating node has not been removed
   *
   * @pre 0 <= @a n1.node_id_ < @a old num_nodes()
   * @pre 0 <= @a n2.node_id_ < @a old num_nodes()
   * @post new num_edges() = old num_edges() - 1
   * @post adjacency edges are updated with new edge ID
   *  
   * Complexity: No more than O(num_nodes()+num_edges()), should be less 
   */

  size_type remove_edge(const Node& n1, const Node& n2) {

    // Check if edge exists
    if (!has_edge(n1, n2)) {
      return 0;
    }

    // Initialize a variabe to keep track of the edge_ID to be removed.
    size_type mid_edge_ID = 0;

    // Remove the n2 from n1's adjacent nodes
    for (auto it = node_edges[n1.node_id_].begin(); it !=node_edges[n1.node_id_].end(); ++it) {
      if ((*it).first == n2.node_id_) {
        // Record down the ID of the edge to be removed
        mid_edge_ID = (*it).second;
        node_edges[n1.index()].erase(it);
        break;
      }
    }

    // Remove the n1 from n2's adjacent nodes
    for (auto it = node_edges[n2.node_id_].begin(); it !=node_edges[n2.node_id_].end(); ++it) {
      if ((*it).first == n1.node_id_) {
        node_edges[n2.index()].erase(it);
        break;
      }
    }

    // Remove the edge and edge_value;
    size_type old_ori_edge_ID = edges.size()-1;
    edges[mid_edge_ID] = edges[old_ori_edge_ID];
    edges.pop_back();

    edge_values[mid_edge_ID] = edge_values.back();
    edge_values.pop_back();

    size_type new_n1_ID = edges[mid_edge_ID].first;
    size_type new_n2_ID = edges[mid_edge_ID].second;

    // Update the node_edge (adjecent node) with the new edge ID
    for (auto it = node_edges[new_n1_ID].begin(); it !=node_edges[new_n1_ID].end(); ++it) {
      if ((*it).first == new_n2_ID) {
        (*it).second = mid_edge_ID;
        break;
      }
    }

    for (auto it = node_edges[new_n2_ID].begin(); it !=node_edges[new_n2_ID].end(); ++it) {
      if ((*it).first == new_n1_ID) {
        (*it).second = mid_edge_ID;
        break;
      }
    }
    return 1; 
  }

  /** Remove a graph's edge
   * @param[in] e    Edge to be removed
   * @return The size_type 1 indicating node has been removed
   *         The size_type 0 indicating node has not been removed
   *
   * @pre 0 <= @a e.node1_id_ < @a old num_nodes()
   * @pre 0 <= @a e.node2_id_ < @a old num_nodes()
   * @post new num_edges() = old num_edges() - 1
   * @post adjacency edges are updated with new edge ID
   *  
   * Complexity: No more than O(num_nodes()+num_edges()), should be less 
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = node(e.node1_id_);
    Node n2 = node(e.node2_id_);
    if (remove_edge(n1, n2)) {
      return 1;
    } 
    return 0;
  }

  /** Remove a graph's edge
   * @param[in] e_it    Edge iterator that points to the edge to be removed
   * @return The edge iterator that points to the next edge to be removed
   *
   * @pre 0 <= @a *e_it.node1_id_ < @a old num_nodes()
   * @pre 0 <= @a *e_it.node2_id_ < @a old num_nodes()
   * @post new num_edges() = old num_edges() - 1
   * @post adjacency edges are updated with new edge ID
   *  
   * Complexity: No more than O(num_nodes()+num_edges()), should be less 
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

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
      return graph_->nodes[node_id_].first;
    }

    /** Return this node's position with modification capability. */
    Point& position() {
      return graph_->nodes[node_id_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_id_;
    }

    // HW1: YOUR CODE HERE

    /**
     * @return The value of the non-constant node_value_type node
     */
    node_value_type& value() {
    	return graph_->nodes[node_id_].second;
    }

    /**
     * @return The value of the constant node_value_type node
     */
    const node_value_type& value() const {
    	return graph_->nodes[node_id_].second;
    }

    /**
     * @return The number of adjacent nodes to the current node
     */
    size_type degree() const {
    	return graph_->node_edges[node_id_].size();
    }

    /**
     * @return An incident iterator that points to the first edge of all the
     * 				 adjacent edges to the node
     */
    incident_iterator edge_begin() const {
    	return incident_iterator(graph_, node_id_, 0);
    }

    /**
     * @return An incident iterator that points to the last edge of all the
     * 				 adjacent edges to the node
     */
    incident_iterator edge_end() const {
    	return incident_iterator(graph_, node_id_, graph_->node_edges[node_id_].size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ && node_id_ == n.node_id_) {
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
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_) {
        if (node_id_ < n.node_id_) {
          return true;
        }
      } else if (graph_ < n.graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type node_id_;
    Node(const Graph* graph, size_type node_id)
      : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   *
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes.push_back(std::make_pair(position, value));

    // Initialize empty vector of adjacent edges for the node
    node_edges.push_back(std::vector<std::pair<size_type, size_type>> ());
    return Node(this, size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.node_id_ < this->nodes.size()) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_id_);
    }

    /** Return the initial length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_) {
        if ((node1_id_ == e.node1_id_ and node2_id_ == e.node2_id_) or 
          (node1_id_ == e.node2_id_ and node2_id_ == e.node1_id_)) {
          return true;
        }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ and node_id_ < e.node_id_) {
      	return true;
      } else if (graph_ < e.graph_) {
        return true;
      }
      return false;
    }

    edge_value_type& value() {
    	return graph_->edge_values[node_id_];
    }

    const edge_value_type& value() const {
    	return graph_->edge_values[node_id_];
    }

   private:
   	// Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node_id_;
    size_type node1_id_;
    size_type node2_id_;

    Edge(const Graph* graph , size_type node_id , size_type node1_id , size_type node2_id)
    : graph_(const_cast<Graph*>(graph)) , node_id_(node_id) , node1_id_(node1_id) , node2_id_(node2_id) {
    }
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
    return Edge(this , i , edges[i].first , edges[i].second);
  }

  /** Test whether two nodes are connected by an edge by iterating
   *  through only the adjecent nodes.
   *
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *         Scenario 1: @a a> @a b, then loop through all the first
   *                     element of adjacent nodes of @a a, if @a b is
   *                     found, then return true
   *         Scenario 2: @a b > @a a, then loop through all the first
   *                     element of adjacent nodes of @a b, if @a a is
   *                     found, then return true
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (unsigned i = 0; i < node_edges[a.node_id_].size(); i++) {
      if (this->node_edges[a.node_id_][i].first == b.node_id_) {
        return true;
      }
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true for two scenarios. Please check
   *       has_edge(@a a, @a b) documentation
   *
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
  	// HW0: YOUR CODE HERE

    // has_edge(@a a, @a b) but return the edge
    if (a.node_id_ > b.node_id_) {
      for (unsigned i = 0; i < node_edges[a.node_id_].size(); i++) {
        if (this->node_edges[a.node_id_][i].first == b.node_id_) {
          return Edge(this, node_edges[a.node_id_][i].second, a.node_id_, b.node_id_);
        }
      }
      // Push back new edge to the back of the vector
      // if same edge is not found
      edges.push_back(std::make_pair(a.node_id_ , b.node_id_));
      edge_values.push_back(value);
    }
    else {
      for (unsigned i = 0; i < node_edges[b.node_id_].size(); i++) {
        if (this->node_edges[b.node_id_][i].first == a.node_id_) {
          return Edge(this, node_edges[b.node_id_][i].second, b.node_id_, a.node_id_);
        }
      }
      // Push back new edge to the back of the vector
      // if same edge is not found
      edges.push_back(std::make_pair(b.node_id_ , a.node_id_));
      edge_values.push_back(value);
    }

    // Push back the adjacent node and the identity ID for that node
    node_edges[a.node_id_].push_back(std::make_pair(b.node_id_, edges.size()-1));
    node_edges[b.node_id_].push_back(std::make_pair(a.node_id_, edges.size()-1));

    // Return the last edge
    return Edge(this, edges.size()-1 , a.node_id_ , b.node_id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    edges.clear();
    node_edges.clear();
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

    /** Constructor for NodeIterator
     * @param[in] graph The type Graph object
     * @param[in] node_size The size_type value
     *
     * @pre 0 <= @a node_size < num_nodes()
     * @post A NodeIterator that points to the node @a node_size of the @a graph
     */
    NodeIterator(const Graph* graph, size_type node_size) : graph_(graph), node_size_(node_size) {
    }

    /** Dereferencing the target NodeIterator
     * @post The node that the NodeIterator points to
     */
    Node operator*() const {
    	return Node(graph_, node_size_);
    }

    /** Increment the target NodeIterator
     * @pre  node_size < num.node()
     *
     * @post NodeIterator points to the next Node if (node_size < num.node())
     *			 NodeIterator points to the nullptr if !(node_size < num.node())
     */
    NodeIterator& operator++() {
    	node_size_++;
    	return *this;
    }

    /** Equality of the current NodeIterator and the target NodeIterator
     * @param[in] node_iterator The target NodeIterator to be compared
     *
     * @return true if the graph that both NodeIterators point to is the same and
     *         the values of the node_size_ of both NodeIterators are the same
     */
    bool operator==(const NodeIterator& node_iterator) const {
    	if (graph_ == node_iterator.graph_ && node_size_ == node_iterator.node_size_) {
    		return true;
    	}
    	return false;
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    const Graph *graph_;
    size_type node_size_;

  };

  // HW1 #2: YOUR CODE HERE

  /** node_iterator of the first node
   * @return An node_iterator that points to the first node of all the nodes
   */
  node_iterator node_begin() const {
  	return NodeIterator(this, 0);
  }

  /** node_iterator of the last node
   * @return An node_iterator that points to the last node of all the nodes
   */
  node_iterator node_end() const {
  	return NodeIterator(this, nodes.size());
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

    /** Constructor for IncidentIterator
     * @param[in] graph The type Graph object
     * @param[in] curr_node The size_type value of the current node
     * @param[in] curr_edge The size_type value of the current edge
     *
     * @pre 0 <= @a curr_node < num_nodes()
     * @pre 0 <= @a curr_edge < degree()
     * @post An IncidentIterator that points to the @a curr_edge of the adjacent
     *       @a curr_node of the @a graph
     */
    IncidentIterator(const Graph* graph, size_type curr_node, size_type curr_edge):
    graph_(graph), curr_node_(curr_node), curr_edge_(curr_edge) {
    }

    /** Dereferencing the target IncidentIterator
     * @post The Edge that the IncidentIterator points to
     */
    Edge operator*() const {
    	return Edge(graph_ , graph_->node_edges[curr_node_][curr_edge_].second ,curr_node_, graph_->node_edges[curr_node_][curr_edge_].first);
    }

    /** Increment the target IncidentIterator
    	*
     * @pre  curr_edge < degree()
     * @post IncidentIterator points to the next Node if (curr_edge < degree())
     *			 NodeIterator points to the nullptr if !(curr_edge < degree())
     */
    IncidentIterator& operator++() {
    	curr_edge_++;
    	return *this;
    }

    /** Equality of the current IncidentIterator and the target IncidentIterator
     * @param[in] incident_iterator The target IncidentIterator to be compared
     *
     * @return true if the graph that all of the conditions are met
     *         1. Both IncidentIterator point to the same graph
     *         2. Both IncidentIterator point to the same node
     *         3. Both IncidentIterator point to the same edge of the same node
     */
    bool operator==(const IncidentIterator& incident_iterator) const {
    	if (incident_iterator.graph_ == graph_ && incident_iterator.curr_node_ == curr_node_
    		&& incident_iterator.curr_edge_ == curr_edge_) {
    		return true;
    	}
    	return false;
    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph *graph_;
    size_type curr_node_;
    size_type curr_edge_;

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

    /** Constructor for EdgeIterator
     * @param[in] graph The type Graph object
     * @param[in] curr_edge The size_type value of the current edge in global order
     *
     * @pre 0 <= @a curr_edge < num.edges()
     * @post An EdgeIterator that points to the @a curr_edge the @a graph
     */
    EdgeIterator(const Graph* graph, size_type curr_edge):
    graph_(graph), curr_edge_(curr_edge) {
    }

    /** Dereferencing the target EdgeIterator
     * @post The Edge that the EdgeIterator points to
     */
    Edge operator*() const {
    	return Edge(graph_, curr_edge_, graph_->edges[curr_edge_].first, graph_->edges[curr_edge_].second);
    }

    /** Increment the target EdgeIterator
     *
     * @pre  curr_edge < num.edges()
     * @post EdgeIterator points to the next Node if (curr_edge < num.edges())
     *			 EdgeIterator points to the nullptr if !(curr_edge < num.edges())
     */
    EdgeIterator& operator++() {
    	curr_edge_++;
    	return *this;
    }

    /** Equality of the current EdgeIterator and the target EdgeIterator
     * @param[in] edge_iterator The target EdgeIterator to be compared
     *
     * @return true if the graph that all of the conditions are met
     *         1. Both EdgeIterator point to the same graph
     *         2. Both EdgeIterator point to the same edge of the same node
     */
    bool operator==(const EdgeIterator& edge_iterator) const {
    	if (curr_edge_ == edge_iterator.curr_edge_ && graph_ == edge_iterator.graph_) {
    		return true;
    	}
    	return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    size_type curr_edge_;
  };

  // HW1 #5: YOUR CODE HERE

  /** edge_iterator of the first edge
   * @return An edge_iterator that points to the first edge of all the edges
   */
  edge_iterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }

  /** edge_iterator of the last edge
   * @return An edge_iterator that points to the last edge of all the edges
   */
  edge_iterator edge_end() const {
  	return EdgeIterator(this, edges.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

    /* A vector of points.
     * x, y and z coordinates are stored at each point.
     * The index number on the vector is same as the node number.
     */
    std::vector<std::pair<Point, node_value_type>> nodes;

    /* A vector of edges.
     * In each edge, there is a pair of node stating which nodes are connected.
     */
    std::vector<std::pair<size_type, size_type>> edges;

    /* A vector of vector of adjacent nodes and edge ID
		 * The outer vector contains the ID of the node
		 * The inner vector contains the identity of the adjacent node
		 * The idetity of the adjacent node is represented by a pair of node ID and
		 * edge ID
     */
    std::vector<std::vector<std::pair<size_type, size_type>>> node_edges;

    std::vector<E> edge_values;

};

#endif // CME212_GRAPH_HPP
