#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <list>
#include <vector>
#include <map>
#include <iterator>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#include <fstream>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph : private totally_ordered<Graph<V,E>> {
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
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  using node_value_type = V;
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
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
  Graph() : nodes(0, InternalNode()), edges(0, InternalEdge()) {
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
  class Node : private totally_ordered<Graph::Node> {
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
    Node() : idx(-1), graph(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph->nodes[idx].pos;
    }

    /** Return this node's position as a reference. */
    Point& position() {
      return graph->nodes[idx].pos;
    }
    /** Return this node's index, a number in the range [0, graph_size). */

    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /*
     * @return a rvalue of node_value_type. For example:
     *
     * @code
     * Node n
     * node_value_type v = n.value();
     * @endcode
     */
    node_value_type& value() {
      return graph->nodes[idx].val;
    }

    /*
     * @return a lvalue of node_value_type. For example:
     *
     * @code
     * Node n;
     * n.value() = node_value_type();
     * @ endcode
     */
    const node_value_type& value() const {
      return graph->nodes[idx].val;
    }

    /*
     * @return number of Edge object adjcent to this node.
     *
     * If the node is invalid, attempting to access it will result in segmentation fault.
     */
    size_type degree() const {
      size_type degree = 0;
      for (auto i = 0; i < graph->adjcent[idx].size(); ++i) {
        auto edge_idx = graph->adjcent[idx][i];
        if (graph->edges[edge_idx].valid) ++degree;
      }
      return degree;
    }

   /*
    * @return An Iterator to the first valid Edge object, that is adjcent to the node.
    *
    * If there are no edges, edge_begin() == edge_end().
    */
    incident_iterator edge_begin() const {
      //TODO: need to return an incident_iterator type
      auto nbr = graph->adjcent[idx];
      if (graph->edges[nbr[0]].valid) return IncidentIterator(idx, 0, graph);
      for(auto i = 0; i < nbr.size(); ++i) {
        auto edge_idx = nbr[i];
        if (graph->edges[edge_idx].valid)
          return IncidentIterator(idx, i, graph);
      }
      return IncidentIterator(idx, nbr.size(), graph);
    }

  /*
   * @return An Iterator to the element following the last valid 
   * Edge object adjcent to this node.
   *
   * This element acts as a placeholder; attempting to access it results in undefined behaviour.
   */
    incident_iterator edge_end() const {
      //TODO: need to return an incident_iterator type
      auto nbr = graph->adjcent[idx];
      return IncidentIterator(idx, nbr.size(), graph);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return n.graph == graph and n.index() == index();
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
      return graph == n.graph ? idx < n.idx : (void *) graph < (void *) n.graph;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type idx;
    Graph* graph;
    Node(size_type _idx, const Graph* _g) : idx(_idx), graph(const_cast<Graph *>(_g)) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    size_type N = 0;
    for(auto n : nodes) {
      if (n.valid) ++N;
    }
    return N;
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
  /*
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    Node new_node = Node();
    new_node.idx = this->size();
    new_node.pos = position;
    new_node.value() = node_value_type();
    new_node.incident = std::vector<Edge>(0, Edge());
    this->nodes.push_back(new_node);
    return new_node;
  }
  */

  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW1: allow node_value_type
    // update this->nodes and this->adjcent accordingly
    size_type idx = this->nodes.size();
    InternalNode new_node(position, val);
    this->nodes.push_back(new_node);
    this->adjcent[idx] = std::vector<size_type>(0);
    return Node(idx, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.index() < nodes.size() and n.graph == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(i, this);
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
  class Edge : private totally_ordered<Graph::Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() : i(-1), j(-1), idx(-1), graph(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(i, graph);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(j, graph);
    }
    edge_value_type& value() {
      return graph->edges[idx].val;
    }

    const edge_value_type& value() const {
      return graph->edges[idx].val;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph == e.graph and idx == e.idx;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return graph == e.graph ? this->idx < e.idx : (void *) graph < (void *) e.graph;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph;
    size_type i;
    size_type j;
    size_type idx;

    Edge(size_type _i, size_type _j, size_type _idx, const Graph* _g) :
      i(_i), j(_j), idx(_idx), graph(const_cast<Graph *>(_g)) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type N = 0;
    for (auto e : edges) {
      if (e.valid) ++N;
    }
    return N;
  }

  /** Return the edge with index @a idx.
   * @pre 0 <= @a idx < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type idx) const {
    // HW0: YOUR CODE HERE
    auto i = edges[idx].i;
    auto j = edges[idx].j;
    return Edge(i, j, idx, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto i = a.index();
    auto j = b.index();
    return find_edge(i, j) < edges.size();
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
    // updated for HW2
    auto i = a.index();
    auto j = b.index();
    size_type idx = find_edge(i, j);
    if (idx < this->edges.size()) {
      return Edge(i, j, idx, this);
    }
    InternalEdge e(i, j);
    // insert at the end
    // update this->edges and this->adj accordingly
    idx = this->edges.size();
    this->edges.push_back(e);
    this->adjcent[i].push_back(idx);
    this->adjcent[j].push_back(idx);
    return Edge(i, j, idx, this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->adjcent.clear();
    this->edges.clear();
    this->nodes.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<Graph::NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : idx(-1), graph(nullptr) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /*
     * @return A valid Node object in the graph.
     */
    Node operator*() const {
      if (idx == graph->nodes.size() || !graph->nodes[idx].valid) return Node();
      return Node(idx, graph);
    }

    /* Increment the iterator
     * @return *this.
     */
    NodeIterator& operator++() {
      ++idx;
      if (idx >= graph->nodes.size()) {
        idx = graph->nodes.size();
      }
      return *this;
    }

    /* Equality check with another iterator
     * @return true or false.
     */
    bool operator==(const NodeIterator& p) const {
      return graph == p.graph and idx == p.idx;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph;
    size_type idx;
    NodeIterator(size_type _idx, const Graph* _g) : 
      idx(_idx), graph(const_cast<Graph *>(_g)) {};
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /*
   * @return An Iterator to the first valid Node object added to the graph.
   *
   * If there are no nodes, node_begin() == node_end().
   */
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }
  /*
   * @return An Iterator to the element following the last valid 
   * Node object added to the graph
   *
   * This element acts as a placeholder; attempting to access it results in undefined behaviour.
   */
  node_iterator node_end() const {
    return NodeIterator(this->size(), this);
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
    IncidentIterator() : idx(-1), graph(nullptr) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /*
     * @return A valid Edge object, adjcent to a Node.
     */
    Edge operator*() const {
      auto nbr = graph->adjcent[node_idx];
      if (idx == nbr.size()) return Edge();
      auto e_idx = nbr[idx];
      auto e = graph->edges[e_idx];
      auto j = e.i == node_idx ? e.j : e.i;
      return Edge(node_idx, j, e_idx, graph);
    }

    /* Increment the iterator
     * @return *this.
     */
    IncidentIterator& operator++() {
      auto nbr = graph->adjcent[node_idx];
      ++idx;
      if(idx >= nbr.size())
        idx = nbr.size();
      return *this;
    }

    /* Equality check with another iterator
     * @return true or false.
     */
    bool operator==(const IncidentIterator& p) const {
      return p.graph == graph and p.node_idx == node_idx
        and p.idx == idx;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type idx;
    size_type node_idx;
    Graph* graph;

    IncidentIterator(size_type _node_idx, size_type _idx,
        const Graph* _g): node_idx(_node_idx), idx(_idx),
    graph(const_cast<Graph*>(_g)) { }
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
    EdgeIterator() : idx(-1), graph(nullptr) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /*
     * @return A valid Edge object in the graph
     */
    Edge operator*() const {
      if (idx == graph->edges.size()) return Edge();
      auto i = graph->edges[idx].i;
      auto j = graph->edges[idx].j;
      return Edge(i, j, idx, graph);
    }
    /* Increment the iterator
     * @return *this.
     */
    EdgeIterator& operator++() {
      ++idx;
      if (idx >= graph->edges.size()) {
        idx = graph->edges.size();
      }
      return *this;
    }

    /* Equality check with another iterator
     * @return true or false.
     */
    bool operator==(const EdgeIterator& p) const {
      return p.graph == graph and p.idx == idx;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph;
    size_type idx;

    EdgeIterator(size_type _idx, const Graph* _g) : idx(_idx),
    graph(const_cast<Graph *>(_g)) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for
  /*
   * @return An Iterator to the first valid Edge object added to the graph.
   *
   * If there are no edges, edge_begin() == edge_end().
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this);
  }
  /*
   * @return An Iterator to the element following the last valid 
   * Edge object added to the graph
   *
   * This element acts as a placeholder; attempting to access it results in undefined behaviour.
   */
  edge_iterator edge_end() const {
    size_type idx = edges.size();
    return EdgeIterator(idx, this);
  }
  //
  // HW2: public interface of remove node and edge
  //
  size_type remove_node(const Node& n) {
    auto idx = n.index();
    return delete_node(idx);
  }

  node_iterator remove_node(node_iterator n_it) {
    auto n = *n_it;
    auto next_idx = remove_node(n);
    return NodeIterator(next_idx, this);
  }

  size_type remove_edge(const Node& a, const Node& b) {
    auto i = a.index();
    auto j = b.index();
    auto e_idx = find_edge(i, j);
    if (e_idx == edges.size()) return false;
    return delete_edge(e_idx);
  }

  size_type remove_edge(const Edge& e) {
    auto e_idx = e.index();
    return delete_edge(e_idx);
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    auto e_idx = remove_edge(e);
    return EdgeIterator(e_idx, this);
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct InternalEdge {
    bool valid = true;
    size_type i;
    size_type j;
    edge_value_type val;
    InternalEdge() : i(-1), j(-1), val(edge_value_type()) {};
    InternalEdge(size_type _i, size_type _j) {
      i = std::min(_i, _j);
      j = std::max(_i, _j);
      val = edge_value_type();
    };
    bool operator<(const InternalEdge& e) const {
      return i == e.i ? j < e.j : i < e.i;
    }
  };

  struct InternalNode {
    bool valid = true;
    Point pos;
    node_value_type val;
    InternalNode() : pos(0), val(node_value_type()) {};
    InternalNode(Point _pos, const node_value_type& _val) : pos(_pos), val(_val) {}
  };

  size_type find_edge(size_type i, size_type j) {
    auto n = adjcent[i].size() < adjcent[j].size() ? i : j;
    auto ii = std::min(i, j);
    auto jj = std::max(i, j);
    for (auto i = 0; i < adjcent[n].size(); ++i) {
      auto e_idx = adjcent[n][i];
      auto e = edges[e_idx];
      if (e.i == ii && e.j == jj) return e_idx;
    }
    return this->edges.size();
  }

  void update_adjcency_matrix(size_type i, size_type j, size_type e) {
    adjcent[i].push_back(e);
    adjcent[j].push_back(e);
  };

  std::vector<InternalNode> nodes;
  std::vector<InternalEdge> edges;
  std::map<size_type, std::vector<size_type> > adjcent;
  //
  // find_edge functionality
  //
  /**
   * find if an edge with node indices a and b exists.
   * @return : the edge index if found; else the number of edges.
   */
  size_type find_edge(const size_type a, const size_type b) const {
    if(adjcent.at(a).empty() or  adjcent.at(b).empty()) return edges.size();
    auto u = adjcent.at(a).size() < adjcent.at(b).size() ? a : b;
    auto aa = std::min(a, b);
    auto bb = std::max(a, b);
    auto& nbr = adjcent.at(u);
    for (auto e_idx : nbr) {
      auto e = edges[e_idx];
      if (e.i == aa and e.j == bb) return e_idx;
    }
    return edges.size();
  }

  //
  // HW2: remove node and edge functionality
  //

  size_type delete_node(size_type node_idx) {
    // 1. remove all associated edges till its degree is zero
    auto nbr = adjcent[node_idx];
    for(auto e : nbr) {
      delete_edge(e);
    }
    // 2. swap this zero degree node with last node, update
    // this->edges and this->adj
    auto last_idx = nodes.size() - 1;
    if(node_idx < last_idx) swap_last_node(node_idx);
    // 3. update this->nodes
    nodes.pop_back();
    return node_idx;
  }

  // swap uid of node to be deleted with last nodes in this->nodes
  // update this->edges and this->adj accordingly
  void swap_last_node(size_type node_idx) {
    auto last_idx = nodes.size() - 1;
    assert(last_idx > node_idx);
    auto last_node_nbr = adjcent[last_idx];
    assert(adjcent[node_idx].size() == 0);
    // update this->edges
    for (auto e_idx : last_node_nbr) {
      auto& e = edges[e_idx];
      assert(e.j == last_idx);
      if (e.i == node_idx) continue; // case where the node to be 
      //swapped has an edge with last node
      auto ii = std::min(e.i, node_idx);
      auto jj = std::max(e.i, node_idx);
      assert(ii < jj);
      e.i = ii;
      e.j = jj;
    }
    // update this->adj
    adjcent[node_idx] = last_node_nbr;
    adjcent[last_idx].clear();
  }

  /**
   * remove the edge with index == edge_idx
   * It performs the operation by first swapping with the last
   * edge in the list and then removing the last edge
   */
  size_type delete_edge(size_type edge_idx) {
    auto last_idx = edges.size() - 1;
    if (last_idx != edge_idx) swap_last_edge(edge_idx);
    delete_last_edge();
    return edge_idx;
  }

  /**
   * swap the last edge with edge with index == edge_idx
   * It updates this->edges and this->adjcent accordingly.
   */
  void swap_last_edge(size_type edge_idx) {
    auto i = edges[edge_idx].i;
    auto j = edges[edge_idx].j;
    auto ii = edges.back().i;
    auto jj = edges.back().j;
    auto last_idx = edges.size() - 1;
    // update this->edges
    std::swap(edges[edge_idx], edges.back());
    // update this->adjcent
    substitute(adjcent[i], edge_idx, last_idx + 1);
    substitute(adjcent[j], edge_idx, last_idx + 1);
    substitute(adjcent[ii], last_idx, edge_idx);
    substitute(adjcent[jj], last_idx, edge_idx);
    substitute(adjcent[i], last_idx + 1, last_idx);
    substitute(adjcent[j], last_idx + 1, last_idx);
  }

  /**
   * remove the edge with index == this.num_edges() - 1
   */
  void delete_last_edge() {
    auto last_idx = edges.size() - 1;
    auto i = edges[last_idx].i;
    auto j = edges[last_idx].j;
    // update this->adjcent
    remove(adjcent[i], last_idx);
    remove(adjcent[j], last_idx);
    // update this->edges
    edges.pop_back();
  }

  /**
   * given an unsorte array @a V, substitute the entry V[i] == old with sub
   */
  void substitute(std::vector<size_type>& v, size_type old, size_type sub) {
    auto it = std::find(v.begin(), v.end(), old);
    if (it != v.end()) *it = sub;
  }

  /**
   * given an unsorted array @a V, remove the entry V[i] == val
   */
  void remove(std::vector<size_type>& v, size_type val) {
    auto it = std::find(v.begin(), v.end(), val);
    std::swap(*it, v.back());
    v.pop_back();
  }
};

#endif // CME212_GRAPH_HPP
