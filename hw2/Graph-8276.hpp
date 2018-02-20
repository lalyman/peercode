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
template <typename V = float, typename E = float>
class Graph {

 public:

  // PUBLIC TYPE DEFINITIONS

    /** Type of this graph. */
    using graph_type = Graph<V, E>;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    /** Type of node attribute. */
    using node_value_type = V;
    using edge_value_type = E;

    static const size_type INVALID_SIZETYPE = size_type(-1);
    static const size_type INVALID_KEY = size_type(-1);

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



 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

  /* An internal structure to store the attributes of each node */
  struct NodeInfo {
    // The key of this node.
    size_type key;

    /** Position each each node */
    Point point;

    /** node attribute */
    node_value_type attr;

    /** Number of neighbors whose index is smaller than this nodeId. */
    size_type deg_s;

    NodeInfo(size_type _key, const Point& _point, const node_value_type& _attr, size_type _deg_s = 0) 
    : key(_key), point(_point), attr(_attr), deg_s(_deg_s) {}
  };

  /** 
   * An internal structure to store the attributes of each edge.
   * This structure is used in the ajList of type vector< vector< EdgeInfo> > 
   *   ajList[nodeId][i] stores the edge information of the edge between nodeId and its ith neighbor
   * The neighbor node id is stored in the dstId, and edge attribute is stored in attr;
   */
  class EdgeInfo : private totally_ordered<EdgeInfo> {
   public:
    EdgeInfo(const size_type& _dstId = INVALID_SIZETYPE, const edge_value_type& _attr = edge_value_type()) 
     : dstId(_dstId), attr(_attr) {}

    /**
     * We define two EdgeInfo to be equal/less_than based on their dstId, which is useful 
     *   when we examine if a node is a neighbor of current node.
     */
    bool operator==(const EdgeInfo& E2) const {
      return this->dstId == E2.dstId;
    }
    bool operator<(const EdgeInfo& E2) const {
      return this->dstId < E2.dstId;
    }

   private:
    friend class Graph;
    /** node2() of this edge */
    size_type dstId;

    /** edge attribute */
    edge_value_type attr;
  };

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : numEdges(0) {
    // nodeIds[INVALID_KEY] = INVALID_SIZETYPE; // No idea why the linker fails if I use INVALID_KEY here...
    nodeIds[size_type(-1)] = INVALID_SIZETYPE;
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
    Node() : G(NULL), key(INVALID_KEY) {
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(G != NULL);
      assert(G->has_node(*this));
      return G->NodesInfo[G->nodeIds.at(key)].point;
    }

    Point& position() {
      assert(G != NULL);
      assert(G->has_node(*this));
      return const_cast<Point&>(G->NodesInfo[G->nodeIds.at(key)].point);      
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(G->has_node(*this));
      return G->nodeIds.at(key);
    }

    /** Return the value of attribute of this node. */
    node_value_type& value() {
      assert(G->has_node(*this));
      return const_cast<node_value_type&>(G->NodesInfo[G->nodeIds.at(key)].attr);
      // return G->NodesInfo[G->nodeIds.at(key)].attr;
    }
    
    const node_value_type& value() const {
      assert(G->has_node(*this));
      return G->NodesInfo[G->nodeIds.at(key)].attr;
    }


    /** Return the degree of node, which is the number of incident edges. */
    size_type degree() const {
      assert(G->has_node(*this));
      return G->ajLists[G->nodeIds.at(key)].size();
    }


    /** Creates an incident_iterator pointing to the first nbr of this node (nbr is smallest index(). */
    incident_iterator edge_begin() const {
      return incident_iterator(G, index(), 0);
    }

    /** Creates an incident_iterator pointing to the place after the last nbr of this node. */
    incident_iterator edge_end() const {
      return incident_iterator(G, index(), this->degree());      
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
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
      if (this->key < n.key) return true;
      if (this->key > n.key) return false;
      return this->G < n.G;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

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
      return NodesInfo.size();
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
      size_type key = rand();
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
      long newKey = getAvailableKey();
      this->nodeIds[newKey] = this->size();
      this->NodesInfo.push_back(NodeInfo(newKey, position, val));
      this->ajLists.push_back(std::vector<EdgeInfo>(0));
      return Node(this, newKey);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
      return n.G == this && n.key != INVALID_KEY && this->nodeIds.count(n.key);
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
      assert(i < this->num_nodes());
      return Node(this, NodesInfo[i].key);
    }

    /** 
     * Remove a node to the graph as well as all its incident edges. If the node
     *     to remove is not the last node (node with index = num_nodes() - 1), this
     *     function will change the last node's id to the removed node index to maintain
     *     consecutive storage
     * @param n   The node to remove
     * @return    number of nodes in the graph after node removal
     *
     * @pre  n is a node of the graph
     * @post new num_nodes() == old num_nodes() - 1
     * @post new num_edges() == old num_edges() - @a n 's degree, 
     * @post invalidate 
     *        1) any outstanding incident_iterator starting at @a n or n's neighbor;
     *        2) any outstanding edge_iterator;
     *        3) any outstanding node denoting @a n;
     *        4) any outstanding edge s.t. @a n is one of its endpoints.
     *        5) any outstanding node_iterator pointing after @a n.
     *
     * Complexity: O(deg_max ^2) =(under sparsity assumption)= o(num_nodes()).
     */
    size_type remove_node(const Node& n) {
      assert(has_node(n));
      size_type nid = n.index();
      // Remove nid from its neighbors ajList
      for (auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it) {
        size_type nbrId = (*e_it).node2().index();
        removeFromOrderedVec(ajLists[nbrId], EdgeInfo(nid));
        if (nid < nbrId) {
          NodesInfo[nbrId].deg_s --;    
        }
      }

      // update the total number of edges in the graph
      numEdges -= n.degree();

      // if nid is not the last node, put the last node at the position of nid to 
      //    maintain consecutive storage.
      size_type end_nid = num_nodes() - 1;
      if (nid < end_nid) {
        // Move the node information to the right place
        nodeIds[NodesInfo[end_nid].key] = nid;
        NodesInfo[nid] = NodesInfo[end_nid];
        ajLists[nid] = ajLists[end_nid];

        // update the ajList of end_nid's neighbors
        for (EdgeInfo e : ajLists[end_nid]) {
          size_type nbrId = e.dstId;
          e.dstId = nid;
          ajLists[nbrId].pop_back(); // Note that end_nid must be the last element in ajLists[nbrId]
          addToOrderedVec(ajLists[nbrId], e);
          if (nid < nbrId) {
            NodesInfo[nbrId].deg_s ++;
            NodesInfo[nid].deg_s --;
          }
        }
      }

      // erase and pop
      nodeIds.erase(n.key);
      ajLists.pop_back();
      NodesInfo.pop_back();
      return num_nodes();
    }

    /** 
     * Remove a node to the graph as well as all its incident edges. If the node
     *     to remove is not the last node (node with index = num_nodes() - 1), this
     *     function will change the last node's id to the removed node index to maintain
     *     consecutive storage
     * @param n_it  A node iterator s.t. *n_it is the node to remove
     * @return      A node iterator pointing to the same memory as the @a n_it, i.e., it
     *              either points to the old last node which has been shift to @a n_it, or
     *              points to node_end().
     *
     * @pre  n is a node of the graph
     * @post new num_nodes() == old num_nodes() - 1
     * @post new num_edges() == old num_edges() - @a n 's degree, 
     * @post invalidate 
     *        1) any outstanding incident_iterator starting at @a n or n's neighbor;
     *        2) any outstanding edge_iterator;
     *        3) any outstanding node denoting @a n;
     *        4) any outstanding edge s.t. @a n is one of its endpoints.
     *        5) any outstanding node_iterator pointing after @a n.
     *
     * Complexity: O(deg_max ^2) =(under sparsity assumption)= o(num_nodes()).
     */
    node_iterator remove_node(node_iterator n_it) {
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
    Edge() : G(NULL), key1(INVALID_KEY), key2(INVALID_KEY) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(G, key1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
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

    // To obtain the length (L2-norm) of an edge
    double length() const {
      return norm(node1().position() - node2().position());
    }

    // To obtain the attribute of the edge
    // Complexity: O(log deg_max).
    edge_value_type& value() {
      Node n1 = node1(), n2 = node2();
      assert(G->has_edge(n1, n2));
      size_type nodeId1 = n1.index(), nodeId2 = n2.index();
      auto it = lower_bound(G->ajLists[nodeId1].begin(), G->ajLists[nodeId1].end(), EdgeInfo(nodeId2)); 
      return const_cast<edge_value_type&>(it->attr);
      // return G->NodesInfo[G->nodeIds.at(key)].attr;
    }
    
    const edge_value_type& value() const {
      Node n1 = node1(), n2 = node2();
      assert(G->has_edge(n1, n2));
      size_type nodeId1 = n1.index(), nodeId2 = n2.index();
      auto it = lower_bound(G->ajLists[nodeId1].begin(), G->ajLists[nodeId1].end(), EdgeInfo(nodeId2)); 
      return it->attr;
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

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
      return numEdges;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: [goal] No more than O(num_nodes() + num_edges()), hopefully less
     *             [real] O(num_nodes())
     * 
     * Explanation: 
     *   When we iterate over edges using edge id @a i, to avoid visiting each edge twice 
     *   (i.e., both (nodeId1, nodeId2) and (nodeId2, nodeId1)), we only output edges such  
     *   that nodeId1 > nodeId2. To implement this, note that NodesInfo[nodeId].deg_s is the number 
     *   of neighbors with id smaller than nodeId, and ajLists[nodeId] is a sorted (increasing) 
     *   vector of EdgeInfo containing neighbor Id as dstId. Therefore, neighbors with smaller id are stored on the 
     *   0, 1, ..., NodesInfo[nodeId].deg_s-1 positions in ajLists[nodeId], which is the half of
     *   each adjacency list we look at.
     */
    Edge edge(size_type i) const {
      assert(i < this->num_edges());
      size_t nodeId1 = 0;
      while (i >= NodesInfo[nodeId1].deg_s) {
        i -= NodesInfo[nodeId1].deg_s;
        nodeId1 ++;
      }
      return Edge(this, NodesInfo[nodeId1].key, NodesInfo[ajLists[nodeId1][i].dstId].key);
    }


    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: [goal] No more than O(num_nodes() + num_edges()), hopefully less
     *             [real] O(log(deg_max))
     */
    bool has_edge(const Node& a, const Node& b) const {
      assert(has_node(a) && has_node(b));
      size_type aidx = a.index();
      size_type bidx = b.index();
      auto it = lower_bound(ajLists[aidx].begin(), ajLists[aidx].end(), EdgeInfo(bidx)); 
      if (it == ajLists[aidx].end()) return false;
      return (*it).dstId == bidx;
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
      assert(!(a == b));
      size_type aidx = a.index(), bidx = b.index();
      if (!this->has_edge(a, b)) {
        addToOrderedVec(ajLists[aidx], EdgeInfo(bidx));
        addToOrderedVec(ajLists[bidx], EdgeInfo(aidx));
        numEdges ++;
        NodesInfo[std::max(aidx, bidx)].deg_s ++;
      }
      return Edge(this, NodesInfo[aidx].key, NodesInfo[bidx].key);
    }

    /** 
     * Remove an edge in the graph. Do nothing if the edge does not exist.
     * @param n1, n2   The two endpoints of the edge to remove.
     * @return    number of edges in the graph after edge removal
     *
     * @post new num_edges() == old num_edges() - 1 if the edge to remove exists 
     * @post invalidate any outstanding incident_iterator starting at n1 or n2, and any edge_iterator,
     *       and any outstanding edge denoting this to-remove edge 
     *
     * Complexity: O(deg_max).
     */
    size_type remove_edge(const Node& n1, const Node& n2) {
      if (!has_edge(n1, n2)) {
        return numEdges;
      }
      size_type idx1 = n1.index();
      size_type idx2 = n2.index();
      removeFromOrderedVec(ajLists[idx1], EdgeInfo(idx2));
      removeFromOrderedVec(ajLists[idx2], EdgeInfo(idx1));
      NodesInfo[std::max(idx1, idx2)].deg_s --;
      numEdges --;
      return numEdges;
    }

    /** 
     * Remove an edge in the graph. Do nothing if the edge does not exist.
     * @param e   The edge to remove.
     * @return    number of edges in the graph after edge removal
     *
     * @post new num_edges() == old num_edges() - 1 if the edge to remove exists 
     * @post invalidate any outstanding incident_iterator starting at n1 or n2, and any edge_iterator,
     *       and any outstanding edge denoting this to-remove edge
     *
     * Complexity: O(deg_max).
     */
    size_type remove_edge(const Edge& e) {
      return remove_edge(e.node1(), e.node2());
    }

    /** 
     * Remove an edge in the graph. Do nothing if the edge does not exist.
     * @param e_it   The edge_iterator pointing to the edge to remove.
     * @return    An edge_iterator pointing to the next edge (i.e., the edge pointed
     *            by ++e_it before the edge removal) or edge_end();
     *
     * @post new num_edges() == old num_edges() - 1 if the edge to remove exists 
     * @post invalidate any outstanding incident_iterator starting at n1 or n2, and any edge_iterator,
     *       and any outstanding edge denoting this to-remove edge
     *
     * Complexity: O(deg_max).
     */
    edge_iterator remove_edge(edge_iterator e_it) {
      Edge e = *e_it;
      remove_edge(e);
      return e_it.settle();
    }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeIds.clear();
    ajLists.clear();
    NodesInfo.clear();
    numEdges = 0;
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
    NodeIterator() : G(NULL), nodeId(INVALID_SIZETYPE) {
    }


    /** De-reference a node_iterator */
    Node operator*() const {
      assert(G != NULL && nodeId < G->size());
      return Node(G, G->NodesInfo[nodeId].key);
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

    /** A graph pointer to the graph this node_iterator points to. */
    const graph_type* G;
    /** The index of node this node_iterator points to (or G.size() to mark the end of node_iterator). */
    size_type nodeId;

    /** Private node_iterator constructor that can only be initialized by a graph
     * to create a valid node_iterator;
     *  It assigns the necessary value of a valid node_iterator according to the @ */
    NodeIterator(const graph_type* _G, size_type _nodeId) : G(_G), nodeId(_nodeId) {}
  };

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
    IncidentIterator() : G(NULL), nodeId(INVALID_SIZETYPE), i(INVALID_SIZETYPE) {}

    // Supply definitions AND SPECIFICATIONS for:

    /** De-reference incident_iterator */
    Edge operator*() const {
      assert(i < G->node(nodeId).degree());
      size_type nbrId = G->ajLists[nodeId][i].dstId;
      return Edge(G, G->NodesInfo[nodeId].key, G->NodesInfo[nbrId].key);
    }

    /** 
     * Increase the incident_iterator by 1 (to point to the next neighbor). 
     * @pre: the incident_iterator does not point to the end of list of neighbors.
     */
    incident_iterator& operator++() {
      assert(i < G->node(nodeId).degree());
      i ++;
      return *this;
    }

    /** 
     * Definition of equality between two incicent_iterators:
     *  they points to the same nbr of same node in the same graph.
     */
    bool operator==(const IncidentIterator& II) const {
      return this->G == II.G && this->nodeId == II.nodeId && this->i == II.i; 
    }

   private:
    friend class Graph;

    /** 
     * The Graph and nodeId where this incicent_iterators comes from. 
     */
    const graph_type* G;
    size_type nodeId;
    /** counter of neighbor this incicent_iterators points to */
    size_type i;    

    /** Private incicent_iterators constructor that can only be initialized by a node
     * to create a valid incicent_iterators;
     *  It assigns the necessary value of a valid incicent_iterators according to the @ */
    IncidentIterator(const graph_type* _G, size_type _nodeId, size_type _i) :
     G(_G), nodeId(_nodeId), i(_i) {}
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
    EdgeIterator() : G(NULL), nodeId(INVALID_SIZETYPE), i(INVALID_SIZETYPE) {
    }


    /** 
     * De-reference am edge_iterator.
     * @post: the returning edge E has E.node1() > E.node2();
     */
    Edge operator*() const {
      assert(G != NULL && nodeId < this->G->size() && i < G->NodesInfo[nodeId].deg_s);
      return Edge(G, G->NodesInfo[nodeId].key, G->NodesInfo[G->ajLists[nodeId][i].dstId].key);
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

    /** 
     * Iteratively update the edge_iterator until either of :
     *  (1) @var it_ points to an edge E satisfying E.node1() > E.node2();
     *  (2) @var it_ points to the end of "edge_lists"
     */
    edge_iterator& settle() {
      while (nodeId < G->size() && i == G->NodesInfo[nodeId].deg_s) {
        i = 0;
        nodeId ++;
      }    
      return *this;  
    }

  };

  
  /** Creates a edge_iterator pointing to the first edge. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0).settle();
  }

  /** Creates a edge_iterator pointing to the place after the last edge. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->size(), 0);
  }







 private:
  /**
   * In this graph implementation, each node has two id's: the index and the key.
   *  index: same as in the instruction, which takes value in {0, 1, ..., n-1} where n is 
   *         the number of nodes in the graph. Note: index values are consecutive!
   *  key: values are not consecutive, and can take any value of size_type (except INVALID_KEY which
   *         is saved for invalid_node).
   * We use consecutive index not only due to the requirement in assignment requirement, but also
   *  from the motivation of storing node attributes (NodeInfo, adjacency list)
   * We also use non-consecutive key to make it convenient to change the index of every outstanding
   *  node object.
   * Each node (and edge) stores only the key, and its index may change due to possible node
   *  removal.
   * The graph uses the following unordered_map to store the mapping key->index, while store the
   *  mapping index->key in NodeInfo[].key.
   */
  std::unordered_map<size_type, size_type> nodeIds; // nodeIds.at(key) = the index of node with key

  /** 
   * We store the edge information using adjacency lists. 
   * ajLists[nodeId] = list of neighbor indices of node nodeId (in increasing order) 
   */
  std::vector<std::vector<EdgeInfo> > ajLists;
  size_type numEdges;


  /** Attribute of each node */
  std::vector<NodeInfo> NodesInfo;


  // Helper functions useful in our graph implementation.

    /** 
     * Add an element to an ordered list (i.e., to put it in the right position)
     * This function is useful in adding edges. Note that our graph stores edge information
     *  in an adjacent list, where each node corresponds to an increasing list of neighbor ids.
     * @para vec: an vector of elements in increasing order
     * @para val: a new element to be add
     *
     * @pre: elements in @a vec are put in increasing put
     * @post: @a vec been added the new element @a val and still in increasing
     */
    template<typename T>
    void addToOrderedVec(std::vector<T>& vec, T val) {
      auto it = lower_bound(vec.begin(), vec.end(), val);
      vec.insert(it, val);
    }

    /** 
     * Remove an element from an ordered list (i.e., find and remove)
     * This function is useful in removing edges. Note that our graph stores edge information
     *  in an adjacent list, where each node corresponds to an increasing list of neighbor ids.
     * @para vec: an vector of elements in increasing order
     * @para val: an element contained in the list to be removed
     *
     * @pre: elements in @a vec are put in increasing order, and @a val is equal to one and only
     *        one element in @a vec.
     * @post: @a val is removed from @a vec, while keeping @a vec increasing and consecutive storage.
     */
    template<typename T>
    void removeFromOrderedVec(std::vector<T>& vec, T val) {
      auto it = lower_bound(vec.begin(), vec.end(), val);
      assert(*it == val);
      vec.erase(it);
    }

};

#endif // CME212_GRAPH_HPP
