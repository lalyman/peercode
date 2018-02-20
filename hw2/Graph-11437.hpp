
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;


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

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
      graph_ = NULL;
      nodeIndex_ = 0;
    }

    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> nodesList[nodeIndex_].pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_ -> nodesList[nodeIndex_].index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      assert(graph_);
      return graph_ -> nodesList[nodeIndex_].val_;
    }

    const node_value_type& value() const {
      assert(graph_);
      return graph_ -> nodesList[nodeIndex_].val_;
    }
    
    size_type degree() const {
      assert(graph_);
      return graph_ -> nodesList[nodeIndex_].adjList_.size();
    }

    incident_iterator edge_begin() const {
      return incident_iterator(graph_, nodeIndex_, graph_ -> nodesList[nodeIndex_].adjList_.begin());
    }

    incident_iterator edge_end() const {
      return incident_iterator(graph_, nodeIndex_, graph_ -> nodesList[nodeIndex_].adjList_.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ && n.graph_) {
        return graph_ == n.graph_ && nodeIndex_ == n.nodeIndex_;
      } else {
        return false;
      }
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
      return tie(graph_, nodeIndex_) < tie(n.graph_, n.nodeIndex_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    size_type nodeIndex_;
    graph_type* graph_;
    Node(const graph_type* graph, size_type nodeIndex) {
      graph_ = const_cast<graph_type*>(graph);
      nodeIndex_ = nodeIndex;
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nodeVal = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type nodeIndex = nodesList.size();
    size_type index = nodes.size();
    nodes.push_back(nodeIndex);
    nodesList.push_back(nodeInfo(position, nodeVal, index));
    return node(index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_ && nodesList[n.nodeIndex_].index_ != -1;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, nodes[i]);     
  }


	/** Remove a designated node if exists and all connected edges
	* @param[in] n  a valid node
	* @return the number of nodes removed
	* 
	* @pre @a n is a node to be removed
	* @post if the node exists in the graph, remove the node and all connected edges
	* @return the number of nodes removed, 0 or 1
	* 
	* Complexity: O(num_nodes()).
	*/
	size_type remove_node(const Node& n) {
		if (!has_node(n)) {
			return 0;		
		}
		while (n.degree()) {
			remove_edge(*(n.edge_begin()));
		}
		if (n.index() != nodes.size() - 1) {
			nodes[n.index()] = nodes[nodes.size() - 1];
			nodesList[nodes[n.index()]].index_ = n.index();
		}
		nodes.pop_back();
		nodesList[n.nodeIndex_].index_ = -1;
		return 1;
	}
	
	/** Remove a designated node if exists and all connected edges
	* @param[in] n  a valid node
	* @return the iterator of the node
	* 
	* @pre @a n_it is a node iterator pointing to the node to be removed
	* @post if the node that @a n_it points to exists in the graph, remove the node and all connected edges
	* @return the iterator
	* 
	* Complexity: O(num_nodes()).
	*/
	node_iterator remove_node(node_iterator n_it) {
		remove_node(*(n_it));
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
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node1_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ && e.graph_) {
        return graph_ == e.graph_ && node1_ == e.node1_ && node2_ == e.node2_;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return tie(graph_, node1_, node2_) < tie(e.graph_, e.node1_, e.node2_);
    }
		// the length of the edge
    double length() const {
      return norm(node1().position() - node2().position());
    }

    edge_value_type& value() {
      return graph_ -> nodesList[node1_].adjList_[node2_];
    }

    const edge_value_type value() const {
      return graph_ -> nodesList[node1_].adjList_[node2_];
    }

		// the other direction of edge
		edge_type dir() const {
			return edge_type(graph_, node2_, node1_);
		}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      
      size_type node1_;
      size_type node2_;
      graph_type* graph_;

      Edge(const graph_type* graph, size_type n1, size_type n2) {
        graph_ = const_cast<graph_type*>(graph);
        node1_ = n1;
        node2_ = n2;
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_nodes_index.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return *next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (b.nodeIndex_ < a.nodeIndex_) {
      return has_edge(b, a);
    }
    if (num_edges() == 0 || !has_node(a) || !has_node(b)) {
      return false;
    }
    return edge_nodes_index.count(pair<size_type, size_type>(a.nodeIndex_, b.nodeIndex_));
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE
    assert(has_node(a) && has_node(b));
    if (a.nodeIndex_ > b.nodeIndex_) {
      return add_edge(b, a, val);
    }
    if (has_edge(a, b)) {
      return edge_type(this, a.nodeIndex_, b.nodeIndex_);
    }
    nodesList[a.nodeIndex_].adjList_[b.nodeIndex_] = val;
    nodesList[b.nodeIndex_].adjList_[a.nodeIndex_] = val;
    edge_nodes_index.insert(pair<pair<size_type, size_type>, size_type>(pair<size_type, size_type>(a.nodeIndex_, b.nodeIndex_), edge_nodes_index.size()));
    edge_index_nodes.push_back(pair<size_type, size_type>(a.nodeIndex_, b.nodeIndex_));
    return edge_type(this, a.nodeIndex_, b.nodeIndex_);
  }

  Edge edge(size_type i) {
    return edge_type(this, edge_index_nodes[i].first, edge_index_nodes[i].second);
  }


  

	/** Remove a designated edge if exists
	* @param[in] n1, n2  two distinct valid edges
	* @return the number of edges removed
	* 
	* @pre @a n1 and @a n2 are two nodes
	* @post if there is an edge between @a n1 and @a n2, remove the edge
	* @return the number of edgs removed, 0 or 1
	* 
	* Complexity: O(num_nodes() + num_edges())).
	*/
	size_type remove_edge(const Node& n1, const Node& n2) {
		
		if (n1.nodeIndex_ > n2.nodeIndex_) {
			return remove_edge(n2, n1);
		} 
		if (!has_edge(n1, n2)) {
			return 0;
		}
		auto index1 = n1.nodeIndex_;
		auto index2 = n2.nodeIndex_;
	
		nodesList[index1].adjList_.erase(index2);
		nodesList[index2].adjList_.erase(index1);
		auto index = edge_nodes_index[pair<size_type, size_type>(index1, index2)];
		
		if (index != edge_index_nodes.size() - 1) {
			edge_index_nodes[index] = edge_index_nodes[edge_index_nodes.size() - 1];
			edge_nodes_index[edge_index_nodes[index]] = index;
		}	
		edge_index_nodes.pop_back();
		edge_nodes_index.erase(pair<size_type, size_type>(index1, index2));
		return 1;
	}
	
	/** Remove a designated edge if exists
	* @param[in] e a valid edge
	* @return the number of edges removed
	* 
	* @pre @a e is a valid edge
	* @post if the edge exists, remove the edge
	* @return the number of edgs removed, 0 or 1
	* 
	* Complexity: O(num_nodes() + num_edges()).
	*/
	size_type remove_edge(const Edge& e) {
		return remove_edge(e.node1(), e.node2());
	}
	

	/** Remove a designated edge if exists
	* @param[in] e_it an edge iterator that points to the edge to be removed
	* @return the iterator of the edge removed
	* 
	* @pre @a e_it is an edge iterator
	* @post if the edge that @a e_it points to exists, remove the edge
	* @return the iterator
	* 
	* Complexity: O(num_nodes() + num_edges()).
	*/	
	edge_iterator remove_edge(edge_iterator e_it) {
		remove_edge(*(e_it));
		return e_it;
	}
	
  /** Remove all nodes and edges from this graph.
  * @post num_nodes() == 0 && num_edges() == 0
  *
  * Invalidates all outstanding Node and Edge objects.
  */
	void clear() {
    // HW0: YOUR CODE HERE
    for (auto itr = nodesList.begin(); itr != nodesList.end(); ++itr) {
      itr -> adjList_.clear();
    }
    nodesList.clear();
    nodes.clear();
    edge_nodes_index.clear();
    edge_index_nodes.clear();
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
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return value_type(graph_, graph_ -> nodes[index_]);
    }

    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

    bool operator==(const NodeIterator& nodeItr) const {
      if (graph_ && nodeItr.graph_) {
        return graph_ == nodeItr.graph_ && index_ == nodeItr.index_;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

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
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return edge_type(graph_, nodeIndex_, itr_ -> first);
    }

    IncidentIterator& operator++() { 
      ++itr_++;
      return *this;
    }

    bool operator==(const IncidentIterator& itr) const {
      if (graph_ && itr.graph_) {
        return graph_ == itr.graph_ && nodeIndex_ == itr.nodeIndex_ && itr_ == itr.itr_;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type nodeIndex_;
    typename unordered_map<size_type, edge_value_type> :: iterator itr_;
    IncidentIterator(const Graph* graph, size_type nodeIndex, typename unordered_map<size_type, edge_value_type> :: iterator itr) {
      graph_ = const_cast<Graph*>(graph);
      nodeIndex_ = nodeIndex;
      itr_ = itr;
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return value_type(graph_, graph_ -> edge_index_nodes[index_].first, graph_ -> edge_index_nodes[index_].second);
    }

    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }

    bool operator==(const EdgeIterator& itr) const {
      if (graph_ && itr.graph_) {
        return graph_ == itr.graph_ && index_ == itr.index_;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;
    EdgeIterator(const graph_type* graph, size_type index) {
      graph_ = const_cast<graph_type*>(graph);
      index_ = index;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_index_nodes.size());
  }
	


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct nodeInfo {
    Point pos_;
    node_value_type val_;
    int index_;
    unordered_map<size_type, edge_value_type> adjList_;
    nodeInfo(const Point& pos, const node_value_type& val, int index) {
      pos_ = pos;
      val_ = val;
      index_ = index;
    }
  };

  vector<nodeInfo> nodesList;
  vector<size_type> nodes;  // for hw2
  vector<pair<size_type, size_type>> edge_index_nodes;
  map<pair<size_type, size_type>, size_type> edge_nodes_index;

};

#endif // CME212_GRAPH_HPP
