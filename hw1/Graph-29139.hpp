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

  /** Predeclaring internal structs */
  struct internal_edge;
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  std::vector<std::vector<unsigned>> map_of_nodes;
  std::vector<internal_edge> edges;
  std::vector<internal_node> nodes;

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
      return node_graph->nodes[node_index].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index;
    }

    /** Method which returns the mutable value of the node. */
    node_value_type& value(){
        return node_graph->nodes[node_index].value;
    };

    /** Const method which returns an IMMUTABLE reference to the node value. */
    const node_value_type& value() const {
        return node_graph->nodes[node_index].value;
    };

    /** Returns the number of incident edges to this node. */
    size_type degree() const{
        return map_of_nodes[node_index].size();
    };

    /** Returns pointer to this nodes "first" incident edge. */
    incident_iterator edge_begin() const {
        return IncidentIterator(node_graph, node_index, 0);
    };

    /** Returns pointer one-past this nodes "last" incident edge. */
    incident_iterator edge_end() const {
        return IncidentIterator(node_graph, node_index,
            node_graph->map_of_nodes[node_index].size());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      (void) n;          // Quiet compiler warning
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
      (void) n;           // Quiet compiler warning
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
  Node add_node(const Point& position, const node_value_type& node_value =
      node_value_type()) {

    internal_node new_node;
    new_node.p = position;
    new_node.index = num_nodes();
    new_node.value = node_value;
    nodes.push_back(new_node);

    // Initialize a spot in map_of_nodes for this node
    std::vector<unsigned> node_neighbors;
    map_of_nodes.push_back(node_neighbors);

    return Node(this, new_node.index);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.node_graph && n.node_index < num_nodes();
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

    // /** Private Constructor */
    // Edge(const Graph* graph, size_type index)
    //     : edge_graph(const_cast<Graph*>(graph)), edge_index(index) {
    // }

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
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    (void) i;           // Quiet compiler warning
    for (int j = 0; j != edges.size(); j++){
        if (edges[j].edge_index == i){
            return Edge(this, edges[j].node1_index, edges[j].node2_index);
        }
    }
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      if (std::find(map_of_nodes[a.node_index].begin(),
      map_of_nodes[a.node_index].end(), b.node_index) !=
      map_of_nodes[a.node_index].end()){
          return true;
      } else {
          return false;
      }
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
    if (this->has_edge(a,b)){
        return Edge(this, a.node_index, b.node_index);
    } else {
        internal_edge newedge;
        newedge.node1_index = a.node_index;
        newedge.node2_index = b.node_index;
        newedge.index = num_edges();
        edges.push_back(newedge);

        map_of_nodes[a.node_index].push_back(b.node_index);
        map_of_nodes[b.node_index].push_back(a.node_index);
        (void) a, (void) b;   // Quiet compiler warning
        return Edge(this, a.node_index, b.node_index);        // Invalid Edge
        }
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    map_of_nodes.clear();
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
  node_iterator node_end() const {
      return NodeIterator(this, nodes.size());
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
     * Specifically, the Edge returned has node_1 = current node, and node_2 =
     * neighbor.
     */
    Edge operator*() const {
        return Edge(iiterator_graph, node_index,
            iiterator_graph->map_of_nodes[node_index][id]);
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

    /** Dereferencing operator returns an Edge from an EdgeIterator. */
    Edge operator*() const {
        return Edge(eiterator_graph, eiterator_graph->edges[id].node1_index, eiterator_graph->edges[id].node2_index);
    };

    /** Increment operator returns an EdgeIterator from an EdgeIterator. */
    EdgeIterator& operator++() {
        ++id;
        return *this;
    };

    /** Checks if input EdgeIterator is equal to current iterator. Graphs and
     * ids must be equal.
     */
    bool operator==(const EdgeIterator& eit) const {
        return (eiterator_graph == eit.eiterator_graph && id == eit.id);
    };

    /** Checks if this edge_iterator is NOT equal to input. */
    bool operator!=(const EdgeIterator& eit) const {
        return !(eiterator_graph == eit.eiterator_graph && id == eit.id);
    };

   private:
    friend class Graph;
    Graph* eiterator_graph;
    size_type id;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type index)
        : eiterator_graph(const_cast<Graph*>(graph)), id(index) {}
  };

  /** Points to the beginning of vector of edges. */
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  };

  /** Points to one-past the end of edges vector. */
  edge_iterator edge_end() const {
      return EdgeIterator(this, edges.size());
  };

 private:
    struct internal_edge {
        size_type node1_index;
        size_type node2_index;
        size_type index;
    };

    struct internal_node {
        Point p;
        size_type index;
        node_value_type value;
    };
};

#endif // CME212_GRAPH_HPP
