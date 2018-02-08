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
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  // PUBLIC TYPE DEFINITIONS

    /** Type of this graph. */
    using graph_type = Graph<V>;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    /** Type of node attribute. */
    using node_value_type = V;

  // Predeclarations 
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

  /** Type of vector of size_type's, useful in storing the adjacency list of each node */
  using szvec = std::vector<size_type>;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : numEdges(0) {
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
    Node() : G(NULL), key(size_type(-1)) {
      // HW0(Done): YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0(Done): YOUR CODE HERE
      assert(G != NULL);
      assert(G->has_node(*this));
      return G->points[G->nodeIds.at(key)];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0(Done): YOUR CODE HERE
      assert(G->has_node(*this));
      return G->nodeIds.at(key);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the value of attribute of this node. */
    node_value_type& value() {
      assert(G->has_node(*this));
      return const_cast<node_value_type&>(G->attrs[G->nodeIds.at(key)]);
    }
    const node_value_type& value() const {
      assert(G->has_node(*this));
      return G->attrs[G->nodeIds.at(key)];
    }


    /** Return the degree of node, which is the number of incident edges. */
    size_type degree() const {
      assert(G->has_node(*this));
      return G->ajLists[G->nodeIds.at(key)].size();
    }


    /** Creates an incident_iterator pointing to the first nbr of this node (nbr is smallest index(). */
    incident_iterator edge_begin() const {
      return incident_iterator(*this, 0);
    }

    /** Creates an incident_iterator pointing to the place after the last nbr of this node. */
    incident_iterator edge_end() const {
      return incident_iterator(*this, this->degree());      
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0(Done): YOUR CODE HERE
      return this->G == n.G && this->key == n.key;
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
      // HW0(Done): YOUR CODE HERE
      if (this->key < n.key) return true;
      if (this->key > n.key) return false;
      return this->G < n.G;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    /** A graph pointer to the graph this node belongs to. */
    const graph_type* G;
    /** The key of this node. Note that this is not the same as the index. 
     *  See the explanations at the definition of private variables of Graph
     */
    size_type key;

    /** Private Node constructor that can only be initialized by a graph
     * to create a valid node;
     *  It assigns the necessary value of a valid node according to the @ */
    Node(const graph_type* _G, size_type _key) : G(_G), key(_key) {}
  };

  // Node-related graph functions
    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
      // HW0(Done): YOUR CODE HERE
      return points.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
      return size();
    }

    /* To generate a new available key for nodes.
     * This key is generated at random, except to make sure it is different
     *  from existing node keys.
     */
    size_type getAvailableKey() {
      long key = rand();
      while (nodeIds.count(key)) key = rand();
      return key;
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] val The new node's modifiable attribute value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
      // HW0(Done): YOUR CODE HERE
      long newKey = getAvailableKey();
      this->nodeIds[newKey] = this->size();
      this->keys.push_back(newKey);
      this->points.push_back(position);
      this->ajLists.push_back(szvec(0));
      this->deg_small.push_back(0);
      this->attrs.push_back(val);
      return Node(this, newKey);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
      // HW0(Done): YOUR CODE HERE
      return n.G == this && this->nodeIds.count(n.key);
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
      // HW0(Done): YOUR CODE HERE
      assert(i < this->num_nodes());
      return Node(this, keys[i]);
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
    Edge() : G(NULL), key1(size_type(-1)), key2(size_type(-1)) {
      // HW0(Done): YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0(Done): YOUR CODE HERE
      return Node(G, key1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0(Done): YOUR CODE HERE
      return Node(G, key2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return this->G == e.G && this->key1 == e.key1 && this->key2 == e.key2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->key1 < e.key1) return true;
      if (this->key1 > e.key1) return false;

      if (this->key2 < e.key2) return true;
      if (this->key2 > e.key2) return false;

      return this->G < e.G;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects


    /** A graph pointer to the graph this node belongs to. */
    const graph_type* G;
    /** The keys of two nodes related to this edge.
     *  Note that the node key is not the same as the node index. 
     *  See the explanations at the definition of private variables of Graph.
     */
    size_type key1, key2;

    /** Private Edge constructor that can only be initialized by a graph
     * to create a valid edge;
     *  It assigns the necessary value of a valid edge according to the @ */
    Edge(const graph_type* _G, size_type _key1, size_type _key2) 
      : G(_G), key1(_key1), key2(_key2) {};
  };

  // Edge-related graph functions
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
      // HW0(Done): YOUR CODE HERE
      return numEdges;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: [goal] No more than O(num_nodes() + num_edges()), hopefully less
     *             [real] O(1)
     */
    Edge edge(size_type i) const {
      // HW0(Done): YOUR CODE HERE
      assert(i < this->num_edges());
      size_t nodeId1 = 0;
      while (i >= deg_small[nodeId1]) {
        i -= deg_small[nodeId1];
        nodeId1 ++;
      }
      return Edge(this, keys[nodeId1], keys[ajLists[nodeId1][i]]);
    }


    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: [goal] No more than O(num_nodes() + num_edges()), hopefully less
     *             [real] O(log(deg_max))
     */
    bool has_edge(const Node& a, const Node& b) const {
      // HW0(Done): YOUR CODE HERE
      size_type aidx = a.index();
      size_type bidx = b.index();
      auto it = lower_bound(ajLists[aidx].begin(), ajLists[aidx].end(), bidx); 
      // Bug!! Will be a bug if I use szvec::iterator instead of auto. No idea why ...
      if (it == ajLists[aidx].end()) return false;
      return *it == bidx;
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
     * Complexity: [Goal] No more than O(num_nodes() + num_edges()), hopefully less
     *             [real] O(log(deg_max))
     */
    Edge add_edge(const Node& a, const Node& b) {
      // HW0(Done): YOUR CODE HERE
      assert(!(a == b));
      size_type aidx = a.index(), bidx = b.index();
      if (!this->has_edge(a, b)) {
        szvec::iterator it = lower_bound(ajLists[aidx].begin(), ajLists[aidx].end(), bidx);
        ajLists[aidx].insert(it, bidx);
        it = lower_bound(ajLists[bidx].begin(), ajLists[bidx].end(), aidx);
        ajLists[bidx].insert(it, aidx);
        numEdges ++;
        deg_small[std::max(aidx, bidx)] ++;
      }
      return Edge(this, keys[aidx], keys[bidx]);
    }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    points.clear();
    ajLists.clear();
    deg_small.clear();
    numEdges = 0;
    keys.clear();
    nodeIds.clear();
    attrs.clear();
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
    NodeIterator() : G(NULL), nodeId(size_type(-1)) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** De-reference a node_iterator */
    Node operator*() const {
      assert(G != NULL && nodeId < G->size());
      return Node(this->G, this->G->keys[nodeId]);
    }

    /** 
     * Increase the node_iterator by 1 (to point to the next node). 
     * @pre: the node iterator does not point to the end of graph nodes.
     */
    node_iterator& operator++() {
      assert(nodeId < G->size());
      nodeId ++;
      return *this;
    }

    /** 
     * Definition of equality between two node_iterators:
     *  they points to the same node in the same graph.
     */
    bool operator==(const NodeIterator& NI) const {
      return this->G == NI.G && this->nodeId == NI.nodeId;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /** A graph pointer to the graph this node_iterator points to. */
    const graph_type* G;
    /** The index of node this node_iterator points to (or G.size() to mark the end of node_iterator). */
    size_type nodeId;

    /** Private node_iterator constructor that can only be initialized by a graph
     * to create a valid node_iterator;
     *  It assigns the necessary value of a valid node_iterator according to the @ */
    NodeIterator(const graph_type* _G, size_type _nodeId) : G(_G), nodeId(_nodeId) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Creates a node_iterator pointing to the first node (nodes with index 0). */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }
  /** Creates a node_iterator pointing to the place after the last node (nodes with highest index). */
  node_iterator node_end() const {
    return node_iterator(this, this->size());    
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
    IncidentIterator() : i(size_type(-1)) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** De-reference incident_iterator */
    Edge operator*() const {
      assert(i < n.degree());
      return Edge(n.G, n.key, n.G->keys[n.G->ajLists[n.index()][i]]);
    }

    /** 
     * Increase the incident_iterator by 1 (to point to the next neighbor). 
     * @pre: the incident_iterator does not point to the end of list of neighbors.
     */
    incident_iterator& operator++() {
      assert(i < n.degree());
      i ++;
      return *this;
    }

    /** 
     * Definition of equality between two incicent_iterators:
     *  they points to the same nbr of same node in the same graph.
     */
    bool operator==(const IncidentIterator& II) const {
      return this->n == II.n && this->i == II.i; 
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    friend class Node;

    /** The node where this incicent_iterators comes from. */
    Node n;
    /** counter of neighbor this incicent_iterators points to */
    size_type i;    

    /** Private incicent_iterators constructor that can only be initialized by a node
     * to create a valid incicent_iterators;
     *  It assigns the necessary value of a valid incicent_iterators according to the @ */
    IncidentIterator(const Node& _n, size_type _i) :
     n(_n), i(_i) {}
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
    EdgeIterator() : G(NULL), nodeId(size_type(-1)), i(size_type(-1)) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** 
     * De-reference am edge_iterator.
     * @post: the returning edge E has E.node1() > E.node2();
     */
    Edge operator*() const {
      assert(G != NULL && nodeId < this->G->size() && i < this->G->deg_small[nodeId]);
      return Edge(this->G, this->G->keys[nodeId], this->G->keys[this->G->ajLists[nodeId][i]]);
    }

    /** 
     * Iteratively update the edge_iterator until either of :
     *  (1) @var it_ points to an edge E satisfying E.node1() > E.node2();
     *  (2) @var it_ points to the end of "edge_lists"
     */
    edge_iterator& settle() {
      while (nodeId < G->size() && i == G->deg_small[nodeId]) {
        i = 0;
        nodeId ++;
      }    
      return *this;  
    }

    /** 
     * Increase the incident_iterator by 1 (to point to the next neighbor). 
     * @pre: the incident_iterator does not point to the end of "edge_lists".
     */
    edge_iterator& operator++() {
      assert(nodeId < G->size());
      i ++;
      return this->settle();
    }

    /** 
     * Definition of equality between two edge_iterators:
     *  they points to the same edge in the same graph.
     */
    bool operator==(const EdgeIterator& EI) const {
      return this->G == EI.G && this->nodeId == EI.nodeId && this->i == EI.i; 
    }

   private:
    friend class Graph;
    // HW1 #5 YOUR CODE HERE

    /** A graph pointer to the graph this edge_iterator points to. */
    const graph_type* G;
    /** The index of node1() this edge_iterator comes from. */
    size_type nodeId;
    /** counter of neighbor of node1() this edge_iterator points to */
    size_type i;    

    /** Private edge_iterator constructor that can only be initialized by a graph
     * to create a valid edge_iterator;
     *  It assigns the necessary value of a valid edge_iterator according to the @ */
    EdgeIterator(const graph_type* _G, size_type _nodeId, size_type _i) :
     G(_G), nodeId(_nodeId), i(_i) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Creates a edge_iterator pointing to the first edge. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0).settle();
  }

  /** Creates a edge_iterator pointing to the place after the last edge. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->size(), 0);
  }

 private:
  /*
   * In this graph implementation, each node has two id's: the index and the key.
   *  index: same as in the instruction, which takes value in {0, 1, ..., n-1} where n is 
   *         the number of nodes in the graph. Note: values are consecutive!
   *  key: values are not consecutive, and can take any value of long.
   * We use consecutive index not only due to the requirement in assignment requirement, but also
   *  from the motivation of storing node attributes (point, adjacency list, deg_small, values, ...)
   * We also use non-consecutive key to make it convenient to change the index of every outstanding
   *  node object.
   * Each node (and edge) stores only the key, and its index may change due to possible node
   *  removal.
   * The graph maintains the mapping between the index and key using the following two variables.
   */
  szvec keys;       // keys[index] = the key of node with index
  std::unordered_map<size_type, size_type> nodeIds; // nodeIds.at(key) = the index of node with key

  /** 
   * We store the edge information using adjacency lists. 
   * ajLists[nodeId] = list of neighbor indices of node nodeId (in increasing order) 
   */
  std::vector<szvec> ajLists;
  size_type numEdges;

  /** Attribute of each node */
  /** Position each each node */
  std::vector<Point> points;
  /** node attribute */
  std::vector<node_value_type> attrs;
  /** deg_small[nodeId] = # of neighbors of node nodeId whose index is smaller than nodeId. */
  szvec deg_small;


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
