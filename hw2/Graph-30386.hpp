#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <utility>
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
// template <typename V>
class Graph {
 private:

  /** Predeclaring internal structs */
  struct internal_edge;
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  std::vector<internal_node> nodes_;
  std::vector<std::vector<std::pair<unsigned, internal_edge>>> adj_;
  std::vector<unsigned> i2u_; // map between indices (i) and unique identifiers (uids).

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Specifies the type of the nodal value - provided by user. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Specifies the type of the edge value - also provided by user. */
  using edge_value_type = E;

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

    /** Return this node's position. */
    const Point& position() const {
      return node_graph->nodes_[node_graph->i2u_[node_index]].p;
    }

    /** Return a mutable reference to this node's position. */
    Point& position() {
      return node_graph->nodes_[node_graph->i2u_[node_index]].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index; // IS THIS RIGHT?!
    }

    /** Method which returns the mutable value of the node. */
    node_value_type& value(){
        return node_graph->nodes_[node_graph->i2u_[node_index]].value;
    };

    /** Const method which returns an IMMUTABLE reference to the node value. */
    const node_value_type& value() const {
        return node_graph->nodes_[node_graph->i2u_[node_index]].value;
    };

    /** Returns the number of incident edges to this node. */
    size_type degree() const{
        return node_graph->adj_[node_index].size();
    };

    /** Returns pointer to this nodes "first" incident edge. */
    incident_iterator edge_begin() const {
        return IncidentIterator(node_graph, node_index, 0);
    };

    /** Returns pointer one-past this nodes "last" incident edge. */
    incident_iterator edge_end() const {
        return IncidentIterator(node_graph, node_index, this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return node_graph == n.node_graph && node_index == n.node_index;
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
      if (node_index < n.node_index){
        return true;
      } else {
        return (node_index == n.node_index && node_graph < n.node_graph);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* node_graph;
    size_type node_index;

    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : node_graph(const_cast<Graph*>(graph)), node_index(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value =
      node_value_type()) {

    internal_node new_node;
    new_node.p = position;
    new_node.value = node_value;
    i2u_.push_back(nodes_.size());
    nodes_.push_back(new_node);
    std::vector<std::pair<size_type, internal_edge>> new_vec;
    adj_.push_back(new_vec);

    return Node(this, i2u_.size() - 1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      for (auto nit = node_begin(); nit != node_end(); ++nit){
          if ((*nit).position() == n.position()){
              return (this == n.node_graph && n.index() == (*nit).index());
          }
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
    return Node(this, i);        // Invalid node
  }

  /**  If @n is a node in graph: Removes @n and all its incident edges
   * from graph
   * Else: Returns 0.
   *
   * @pre @n is a valid node.
   * @post new_num_nodes() = old_num_nodes() - 1
   * @post new_num_edges() = old_num_edges() - n.degree()
   * @return 0
   *
   * Complexity: O(num_edges())
   */
  size_type remove_node(const Node& n){
      if (has_node(n)){
          // Removing node from i2u_ vector
          i2u_[n.index()] = i2u_.back();
          i2u_.pop_back();

          // Removing nodes incident edges from adj_
          adj_[n.index()] = adj_.back();
          adj_.pop_back();

          // Removing all other instances of node in adj_
          for (size_type i = 0; i < adj_.size(); i++){
              for (size_type j = 0; j < adj_[i].size(); j++){
                  if (adj_[i][j].first == n.index()){
                      adj_[i].erase(adj_[i].begin() + j);
                  }
              }
          }
          // Decreasing edge count
          edge_count -= n.degree();

          return 0;
      } else {
          return 0;
      }
  };

  /** If edge exists: Removes edge defined by @n1 and @n2 from graph.
   * Else: Returns 0.
   *
   * @pre @n1 and @n2 are valid nodes
   * @post new_num_edges() = old_num_edges() - 1
   *
   * Complexity: O(num_edges)
   */
  size_type remove_edge(const Node& n1, const Node& n2){
      if (has_edge(n1,n2)){
          // Remove edge from n1 adjacency list.
          for (size_type i = 0; i < adj_[n1.index()].size(); i++){
              if (adj_[n1.index()][i].first == n2.index()){
                  adj_[n1.index()][i] = adj_[n1.index()].back();
                  adj_[n1.index()].pop_back();
              }
          }
          // Remove edge from n2 adjacency list.
          for (size_type i = 0; i < adj_[n2.index()].size(); i++){
              if (adj_[n2.index()][i].first == n1.index()){
                  adj_[n2.index()][i] = adj_[n2.index()].back();
                  adj_[n2.index()].pop_back();
              }
          }
          // Decrease total number of edges by 1.
          edge_count -= 1;
          return 0;
      } else {
          return 0;
      }
  };

  /** Calls main remove_edge function (above). Same specification. */
  size_type remove_edge(const Edge& e){
      Node n1 = e.node1();
      Node n2 = e.node2();
      return remove_edge(n1,n2);
  };

  /** If edge is in graph: Removes edge at position marked by iterator
   * @e_it
   * Else: Returns null iterator.
   *
   * @pre @e_it is a valid edge iterator
   * @post new_num_edges() = old_num_edges() - 1
   * @return edge_iterator pointing to new_edge at old_edge position
   *
   * Complexity: O(num_nodes())
   */
  edge_iterator remove_edge(edge_iterator e_it){
      Edge e = *e_it;
      if (has_edge(e.node1(),e.node2())){
          size_type index1 = e.node1().index();
          size_type index2 = e.node2().index();
          for (unsigned i = 0; i < adj_[index1].size(); i++){
              if (adj_[index1][i].first == index2){
                  adj_[index1].erase(adj_[index1].begin() + i);
                  return EdgeIterator(this, index1, i);
              }
          }
          edge_count -= 1;
      } else {
          return EdgeIterator(this, 0, 0);
      }
  };


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
      return Node(edge_graph, node1_index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(edge_graph, node2_index);
    }

    /** Returns the value of this Edge by reference */
    edge_value_type& value() {
        // Value is assigned to edge (n1, n2) where n1.index() < n2.index().
        if (node1_index < node2_index){
            for (size_type i = 0; i < edge_graph->adj_[node1_index].size(); i++){
                if (edge_graph->adj_[node1_index][i].first == node2_index){
                    return edge_graph->adj_[node1_index][i].second.value;
                }
            }
        } else {
            // If n1.index() > n2.index(), then find the value associated with edge (n2, n1).
            for (size_type i = 0; i < edge_graph->adj_[node2_index].size(); i++){
                if (edge_graph->adj_[node2_index][i].first == node1_index){
                    return edge_graph->adj_[node2_index][i].second.value;
                }
            }
        }
        return edge_graph->adj_[node1_index][0].second.value;
    }

    /** Returns an IMMUTABLE value of this Edge by reference
     *
     * Same format as above function.
     */
    const edge_value_type& value() const{
        if (node1_index < node2_index){
            for (size_type i = 0; i < edge_graph->adj_[node1_index].size(); i++){
                if (edge_graph->adj_[node1_index][i].first == node2_index){
                    return edge_graph->adj_[node1_index][i].second.value;
                }
            }
        } else {
            for (size_type i = 0; i < edge_graph->adj_[node2_index].size(); i++){
                if (edge_graph->adj_[node2_index][i].first == node1_index){
                    return edge_graph->adj_[node2_index][i].second.value;
                }
            }
        }
        return edge_graph->adj_[node1_index][0].second.value;
    }

    /** Returns Euclidean length of the edge */
    double length() const {
        Point& pos1 = edge_graph->nodes_[edge_graph->i2u_[node1_index]].p;
        Point& pos2 = edge_graph->nodes_[edge_graph->i2u_[node2_index]].p;
        return norm(pos1 - pos2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      Node node1 = Node(edge_graph, node1_index);
      Node node2 = Node(edge_graph, node2_index);
      return (node1 == e.node1() && node2 == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      Node node1 = Node(edge_graph, node1_index);
      Node node2 = Node(edge_graph, node2_index);

      if (node1 < e.node1()){
        return true;
      } else if (node1 == e.node1() && node2 < e.node2()) {
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* edge_graph;
    size_type node1_index;
    size_type node2_index;

    /** Private Constructor */
    Edge(const Graph* graph, const size_type node1, const size_type node2)
        : edge_graph(const_cast<Graph*>(graph)), node1_index(node1),
        node2_index(node2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    EdgeIterator eit_ = this->edge_begin();
    std::advance(eit_,i);
    return (*eit_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      for (size_type i = 0; i < adj_[a.index()].size(); i++){
          if (adj_[a.index()][i].first == b.index()){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value = edge_value_type()) {
    if (this->has_edge(a,b)){
        return Edge(this, a.index(), b.index());
    } else {
        internal_edge newedge;
        newedge.node1_index = a.index();
        newedge.node2_index = b.index();
        newedge.value = edge_value;

        // Adding pair corresponding to new edge to @a's and @b's adjacency list
        std::pair<size_type, internal_edge&> n1_pair (b.index(), newedge);
        std::pair<size_type, internal_edge&> n2_pair (a.index(), newedge);
        adj_[a.index()].push_back(n1_pair);
        adj_[b.index()].push_back(n2_pair);

        // Increasing edge count by 1
        edge_count += 1;

        return Edge(this, a.index(), b.index());        // Invalid Edge
        }
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_count = 0;
    nodes_.clear();
    i2u_.clear();
    adj_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   */
  class NodeIterator {
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

    /** Dereferencing operation. Dereferenced NodeIterator returns a Node. */
    value_type operator*() const {
        return Node(niterator_graph, id);
    };

    /** Increment operation. Incrementing a NodeIterator should return a
     * NodeIterator&.
     */
    node_iterator& operator++() {
        ++id;
        return *this;
    }

    /** Checks if this node_iterator is equal to input. */
    bool operator==(const node_iterator& nit) const {
        return (niterator_graph == nit.niterator_graph && id == nit.id);
    };

    /** Checks if this node_iterator is NOT equal to input. */
    bool operator!=(const node_iterator& nit) const {
        return !(niterator_graph == nit.niterator_graph && id == nit.id);
    };

   private:
    friend class Graph;

    Graph* niterator_graph;
    size_type id;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type index)
        : niterator_graph(const_cast<Graph*>(graph)), id(index) {}
  };

  /** Points to the first element in the collection of nodes. */
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  };

  /** Points to one past the last element in the collection of nodes. */
  node_iterator node_end() const{
      return NodeIterator(this, i2u_.size());
  };

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   */
  class IncidentIterator {
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

    /** Dereferencing operation. Dereferenced IncidentIterator returns an Edge.
     * The Edge returned has node_1 = current node, and node_2 =
     * neighbor.
     */
    value_type operator*() const {
        size_type node2_index = iiterator_graph->adj_[node_index][id].first;
        return Edge(iiterator_graph, node_index, node2_index);
    };

    /** Incrementing the incident edge. Returns IncidentIterator. */
    IncidentIterator& operator++() {
        ++id;
        return *this;
    };

    /** Checks if current iterator is equal to input iterator. Returns a
     * boolean value.
     */
    bool operator==(const IncidentIterator& iit) const {
        return (iiterator_graph == iit.iiterator_graph && node_index == iit.node_index && id == iit.id);
    };

    /** Checks if this incident_iterator is NOT equal to input. */
    bool operator!=(const IncidentIterator& iit) const {
        return !(iiterator_graph == iit.iiterator_graph && node_index == iit.node_index && id == iit.id);
    };

   private:
    friend class Graph;
    Graph* iiterator_graph;
    size_type node_index;
    size_type id;

    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type node, size_type index)
        : iiterator_graph(const_cast<Graph*>(graph)), node_index(node), id(index) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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

    /** Dereferencing operator returns an Edge. */
    value_type operator*() const {
        size_type node2_index = eiterator_graph->adj_[outer_id][inner_id].first;
        return Edge(eiterator_graph, outer_id, node2_index);
    };

    /** Increment operator returns an EdgeIterator. */
    EdgeIterator& operator++() {
        // Keep doing this loop indefinitely :|
        // (Probably not the best)
        while (true) {
            // Handles the case where inner_increment reaches the end of the
            // inner_vector
            while ((inner_id + 1) >= (eiterator_graph->adj_[outer_id].size())){
                inner_id = 0;
                ++outer_id;
                // Handles case where we reach the end of the adjacency list
                if (outer_id == eiterator_graph->adj_.size()) {
                    return *this;
                }
            }
            ++inner_id;
            if (outer_id < eiterator_graph->adj_[outer_id][inner_id].first){
                // Only return iterators to edges where e.node1().index() <
                // e.node2().index()
                return *this;
            }
        }
        return *this; // Should never reach this return statement
    };



    /** Checks if input EdgeIterator is equal to current iterator. Graphs and ids must be equal.
     */
    bool operator==(const EdgeIterator& eit) const {
        return (eiterator_graph == eit.eiterator_graph && outer_id == eit.outer_id && inner_id == eit.inner_id);
    };

    /** Checks if this edge_iterator is NOT equal to input. */
    bool operator!=(const EdgeIterator& eit) const {
        return !(eiterator_graph == eit.eiterator_graph && outer_id == eit.outer_id && inner_id == eit.inner_id);
    };

   private:
    friend class Graph;
    Graph* eiterator_graph;
    size_type outer_id;
    size_type inner_id;
    // size_type id;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type out, size_type in)
        : eiterator_graph(const_cast<Graph*>(graph)), outer_id(out), inner_id(in) {}
  };

  // Points to the beginning of the adjacency list
  edge_iterator edge_begin() const {
      size_type i = 0;
      while (adj_[i].size() == 0){
          if (i+1 < adj_.size()){
              ++i;
          } else {
              return EdgeIterator(this, adj_.size(), 0);
          }
      }
      return EdgeIterator(this, i, 0);
  };

  // Points to one-past the end of the adjacency list.
  edge_iterator edge_end() const {
      return EdgeIterator(this, adj_.size(), 0);
  };

 private:
    struct internal_edge {
        size_type node1_index;
        size_type node2_index;
        edge_value_type value;
    };

    struct internal_node {
        Point p;
        node_value_type value;
    };

    size_type edge_count;
};

#endif // CME212_GRAPH_HPP
