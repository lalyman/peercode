#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
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
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return get_node().node_pt; //Get the node and A node's position
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
        //internal_node& node = graph_->nodes_[uid_];
        return get_node().node_value; //???
    }
    const node_value_type& value() const {
        //internal_node& node = graph_->nodes_[uid_];
        return get_node().node_value; //???
    }
    //Look in the adjacency list and find how many incident nodes/edges there are
    size_type degree() const {
        return graph_->adj_list[uid_].size();
    }
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_,uid_,0); ///????
    }
    incident_iterator edge_end() const {
        return IncidentIterator(graph_,uid_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        if ((n.graph_ == this->graph_) && (n.uid_ == this -> uid_)) {
            return true;
        }
        else {
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
        if ((n.graph_ == this->graph_) && (this->uid_ < n.uid_)) {
            return true;
        }
        else {
            return false;
        }
    }


   private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Allow Graph to access Node's private member data and functions.
    friend class Graph<V>;
    // HW0: YOUR CODE HERE

    // Pointer back to the Graph
    Graph* graph_;

    //Node's unique identification number
    size_type uid_;

    /** Private Constructor **/
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {

        }

        //Helper function gets the internal node referenced in nodes_ vector
    internal_node& get_node() const {
        //Check that node is in the graph
        if (!((uid_ >= 0) && uid_ < graph_ -> size())) {
            std::cout << "The uid_ is: " << uid_ <<'\n';
            std::cout << "The size is: " << graph_ -> size() <<'\n';
        }
        assert(uid_ >= 0 && uid_ < graph_ -> size());
        return graph_ -> nodes_[uid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
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
   //???
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node newNode;
    newNode.node_pt = position;
    newNode.node_index = nodes_.size();
    newNode.node_value = value;
    nodes_.push_back(newNode);

    std::vector<size_type> no_neigh; //A new node has no new neighbors
    adj_list.push_back(no_neigh);   //so add an empty vector in at that position
    return Node(this,nodes_.size() - 1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.index < this->num_nodes);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this,i);        // Invalid node
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
      return Node(graph_,n1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,n2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (((n1_ == e.n1_ && n2_ == e.n2_) ||
      (n1_ == e.n2_ && n2_ == e.n1_) )&& graph_ == e.graph_) {
          return true;
      }
      else {
          return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if ((uid_ < e.uid_) && graph_ == e.graph_) {
            return true;
        }
        else {
            return false;
        }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph<V>;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph
    Graph* graph_;

    // Edge's unique identification number
    size_type uid_;
    size_type n1_;
    size_type n2_;
    /** Private Constructor **/
    Edge(const Graph* graph, size_type n_1, size_type n_2)
        : graph_(const_cast<Graph*>(graph)), n1_(n_1), n2_(n_2) {

        }

    internal_edge& get_edge() const {
        std::cout << "uid: " <<uid_ <<'\n';
        assert((uid_ >= 0) && uid_ < graph_->num_edges());
        return graph_->edges_.at(uid_);
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if(i < num_edges()) {

        return Edge(this,edges_[i].node_1,edges_[i].node_2);
    }
    else {
        return Edge(); // Invalid Edge
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE If node A is connected to B, return true
    bool flag;
    for(size_type i : adj_list[a.index()]) {
        flag = ((edge(i).node1() == a) && edge(i).node2() == b);
    };
    if (flag == True) {
        std::cout << "Has edge" << '\n';
        return true;
    }
    else {
        return false;
    }
}


  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.n1() == @a a and e.node2() == @a b
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
    // if(has_edge(a,b) == true) {
    //     return Edge(this,a.index(),b.index());
    // }
        internal_edge newEdge;
        newEdge.node_1 = a.index();
        newEdge.node_2 = b.index();
        newEdge.uid_ = edges_.size();
        edges_.push_back(newEdge);

        //unordered_map implementation
        //adj_list[newEdge.node_1][newEdge.node_2] = this->edges_.size() -1;
        //adj_list[newEdge.node_2][newEdge.node_1] = this->edges_.size() -1;

        adj_list[newEdge.node_1].push_back(newEdge.uid_);
        adj_list[newEdge.node_2].push_back(newEdge.uid_);
        return Edge(this,a.index(),b.index());        // Invalid Edge

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      nodes_.clear();
      edges_.clear();
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

    //Return the node corresponding to the given node index
    Node operator*() const {
        return graph_->node(uid_);
    }

    //If valid node, increment the node index and return it
    NodeIterator& operator++() {
        if (uid_ < graph_->size()) {
            ++uid_;
            return *this;
        }
        else {
            return *this;
        }
    }

    //Ensure that the graphs are the same and the node indices are the same
    bool operator==(const NodeIterator& x) const {
        return (x.graph_ == graph_ && x.uid_ == uid_);
    }

   private:
    friend class Graph<V>;
    //Private Constructor of NodeIterator Class

    // HW1 #2: YOUR CODE HERE

    //Pointer to graph and index
    graph_type* graph_;
    size_type uid_;
    NodeIterator(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  //Node iterator beings at the beginning of the internal node vector
  node_iterator node_begin() const {
      return NodeIterator(this,0);
  }
  //Node iterator ends at the end of the internal node vector
  node_iterator node_end() const {
      return NodeIterator(this,size());
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

    //Dereference an incident operator by ensuring that the proper node
    // is selected, then return the Edge with its two node indices
    Edge operator*() const {
        if (node_idx_ != graph_->edges_[iter_].node_2)
        {
            return Edge(graph_,node_idx_,graph_->edges_[iter_].node_2);
        }
        else {
            return Edge(graph_,graph_->edges_[iter_].node_2, node_idx_);
        }

    }
    IncidentIterator& operator++() {
        if(iter_ < graph_->adj_list[node_idx_].size()) {
            ++iter_;
            return *this;
        }
        else {
            return *this;
        }
    }
    bool operator==(const IncidentIterator& x) const {
        return (x.graph_ == graph_ && x.node_idx_ == node_idx_ &&
                x.iter_ == iter_);
    }

   private:
    friend class Graph<V>;
    // HW1 #3: YOUR CODE HERE
    //Private Constructor of Incident Iterator
    graph_type* graph_; // Pointer to the graph so we can use edge()
    size_type node_idx_; //id of the node that we're looking for incidence with
    size_type iter_; //edges incident to a node are stored in a vector at the index of node in
                    // so we can iterator over index of the edges vector
    IncidentIterator(const graph_type* graph, size_type node_idx, size_type iter)
        : graph_(const_cast<graph_type*>(graph)), node_idx_(node_idx),
            iter_(iter) {
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
        return graph_->edge(uid_);
    }
    EdgeIterator& operator++() {
        if (uid_ < graph_->size()) {
            ++uid_;
            return *this;
        }
        else {
            return *this;
        }
    }
    bool operator==(const EdgeIterator& x) const {
        return (x.graph_ == graph_ && x.uid_ == uid_);
    }

   private:
    friend class Graph<V>;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type uid_;
    EdgeIterator(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
      return EdgeIterator(this,0);
  }
  edge_iterator edge_end() const {
      return EdgeIterator(this,size());
  }
//Compare(const Point& p) : p {} {
    //return (n1()-p > node2()-p)

 private:
     struct internal_node {
         Point node_pt;
         int node_index;
         node_value_type node_value;
     };

     //Indices of two nodes that the edges touches
     struct internal_edge {
         size_type uid_;
         size_type node_1;
         size_type node_2;
     };
  // HW0: YOUR CODE HERE
  std::vector<internal_node> nodes_; //Pair is (Point, Node's Value)
  std::vector<internal_edge> edges_; //Each index is an edge, at each index is the vector of nodes that edge is incident to
  //std::unordered_map<size_type,std::unordered_map<size_type,size_type>> adj_list;
  std::vector<std::vector<size_type>> adj_list; //adjacency list
  std::map<size_type, std::set<size_type>> nodeEdges;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
