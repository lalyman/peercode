#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <unordered_map>

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
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  typedef V node_value_type;
  typedef E edge_value_type;

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  //using node_value_type = V;
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

  //** Type of data structure that records data about edges, a nested unordered_set
  //** in an unordered_map */
  using adjacency_map_type = std::unordered_map<size_type,
      std::unordered_map<size_type,size_type>>;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
        : nodes_(), edges_(), adj_map(), num_edges_(0)
  {

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
      return (graph_->nodes_[uid_]).position; //Get the node and A node's position
    }

    /** Node Position Made Modifiable **/
    Point& position() {
        return (graph_->nodes_[uid_]).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
        //internal_node& node = graph_->nodes_[uid_];
        return (graph_->nodes_[uid_]).value; //???
    }
    const node_value_type& value() const {
        //internal_node& node = graph_->nodes_[uid_];
        return (graph_->nodes_[uid_]).value; //???
    }
    //Look in the adjacency list and find how many incident nodes/edges there are
    size_type degree() const {
        //if the node is in the graph, return the number of edges it has
        if(graph_->adj_map.find(uid_) != graph_->adj_map.end()) {
            return (size_type) (graph_->adj_map[uid_]).size();
        }
        else {
            return 0;
        }
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,uid_,(graph_->adj_map[uid_]).begin());
    }

    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, (graph_->adj_map[uid_]).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
         return (graph_ == n.graph_ && uid_ == n.uid_);
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
        return (graph_ == n.graph_ && uid_ < n.uid_) || graph_ < n.graph_;
    }


   private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Allow Graph to access Node's private member data and functions.
    friend class Graph<V,E>;
    // HW0: YOUR CODE HERE

    // Pointer back to the Graph
    graph_type* graph_;

    //Node's unique identification number
    size_type uid_;

    /** Private Constructor **/
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {

        }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return (size_type) i2u_.size();
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

  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE

    //Add a new internal node with appropriate position, value, and index
    nodes_.push_back({position,value,size(), false});
    size_type uid = nodes_.size() - 1;
    i2u_.push_back(uid);
    return Node{this,uid};
}


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //????
    return (nodes_[n.index()].removed != true && (n.graph_ == this));
    //return nodes_[i2u_[n.index()]] == n.index()
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node{this,i2u_[i]};        // Invalid node
  }

  /** Removes a node from a graph and the edges incident to that node
   * @param[in]     n      Const Node&, the node we want to remove
   * @return 1 when node is sucessfully removed
   *
   * @tparam n: Node @n must be in the graph
   */
  size_type remove_node(const Node& n) {
      //Remove if has node
      //if (has_node(n)) {
          //Remove incident edges
          for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
              this->remove_edge(*iit);
          }
          i2u_[n.index()] = i2u_.back();
          i2u_.pop_back();
          nodes_[n.index()].removed = true;
          nodes_[i2u_[n.index()]].index = n.index();
          return 1;
      //}
      //else {
          //return 0;
      //}
}

/** Removes a node from a graph and the edges incident to that node
 *    when given a node iterator by defrencing the node iterator and
 *    calling remove_node() on the dereferenced node
 * @param[in]     n_it      The node iterator, the node we want to remove
 * @returns n_it
 *
 * @tparam n_it: Iterator must point to a valid in the graph
 */
node_iterator remove_node(node_iterator n_it) {
    this->remove_node(*n_it);
    return n_it;
}

/** Removes an edge from the graph
 * @param[in]     a      Const Node&, the edge's root node
 * @param[in]     b      Const Node&, the edge's other node
 * @return 1 if edge is in graph and successfully removed , else return 0
 */
size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a,b) ) {
        internal_edge e = (edges_.back());
        size_type e_idx = adj_map.at(a.index()).at(b.index());
        size_type aidx = a.index();
        size_type bidx = b.index();
        edges_[e_idx] = edges_.back();
        edges_.pop_back();

        adj_map[e.node1][e.node2] = e_idx;
        adj_map[e.node2][e.node1] = e_idx;
        adj_map.at(a.index()).erase(b.index());
        --num_edges_;
        return 1;
    }

    else {
        return 0;
    }

}

/** Removes an inputted edge from the graph by calling remove_edge() on the
 *  edge's two incident nodes
 * @param[in]     e      Const Edge&, the edge to remove
 * @param[in]     b      Const Node&, the edge's other node
 * @return 1 if edge is in graph and successfully removed , else return 0
 */
size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),e.node2());
}

/** Removes an  edge from the graph when passing an edge iterator by calling
 *  remove_edge() on the dereferenced edge iterator
 *  edge's two incident nodes
 * @param[in]     e_it      edge_iterator
 * @returns e_it
 */
edge_iterator remove_edge(edge_iterator e_it) {
    this->remove_edge(*e_it);
    return e_it;

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

    /** Length of edge is Euclidean Distance between its nodes' positions **/

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node2_);      // Invalid Node
    }
    double length() const {
        return norm(node1().position() - node2().position());
    }

  edge_value_type& value() {
      //internal_node& node = graph_->nodes_[uid_];
      return graph_->edges_[graph_->adj_map[node1_][node2_]].value;
  }
  const edge_value_type& value () const {
      return graph_->edges_[graph_->adj_map[node1_][node2_]].value;
  }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if (e.graph_ != graph_) {
            return false;
        }
        return ( (node1_ == e.node1_ && node2_ == e.node2_) ||
            (node1_ == e.node2_ && node2_ == e.node1_) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if (e.graph_ != graph_) {
            return graph_ > e.graph_;
        }
        return ( (node1_ < e.node1_) ||
                (node1_ == e.node1_ && node2_ < e.node2_) );
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph<V,E>;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph
    graph_type* graph_;

    // Edge's unique identification number
    size_type node1_;
    size_type node2_;
    /** Private Constructor **/
    Edge(const graph_type* graph, size_type n_1, size_type n_2)
        : graph_(const_cast<graph_type*>(graph)), node1_(n_1), node2_(n_2) {

        }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;

    // size_type count_edges = 0;
    // for (unsigned k = 0; k < num_nodes(); ++k) {
    //   for (unsigned j = k+1; j < num_nodes(); ++j) {
    //     if (has_edge(node(k), node(j)))
    //       ++count_edges;
    //   }
    // }
    // return count_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge{this, edges_[i].node1, edges_[i].node2};
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE If node A is connected to B, return true
        if (this != a.graph_ || this != b.graph_) {
            return false;
        }
        if (nodes_[a.index()].removed == true &&
            nodes_[b.index()].removed == true) {
                return false;
            }
        if (nodes_[a.index()].removed == true ||
            nodes_[b.index()].removed == true)  {
                return ( ((adj_map.find(a.index()) != adj_map.end())
                      && (adj_map.at(a.index()).find(b.index()) )
                         != adj_map.at(a.index()).end()) ||

                         ((adj_map.find(b.index()) != adj_map.end())
                               && (adj_map.at(b.index()).find(a.index())
                                  != adj_map.at(b.index()).end())));
        }
        else {
            return ( (adj_map.find(a.index()) != adj_map.end())
                  && (adj_map.at(a.index()).find(b.index())
                     != adj_map.at(a.index()).end()))
                     &&

                     ((adj_map.find(b.index()) != adj_map.end())
                           && (adj_map.at(b.index()).find(a.index())
                              != adj_map.at(b.index()).end()));

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

  Edge add_edge(const Node& a, const Node& b, edge_value_type value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if(!has_edge(a, b)) {
     adj_map[a.index()][b.index()] = num_edges_;
     adj_map[b.index()][a.index()] = num_edges_;
     edges_.push_back({a.index(), b.index(), num_edges_,value});
     num_edges_ = num_edges_ + 1;
    }

    return Edge{this, a.index(), b.index()};

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      nodes_.clear();
      i2u_.clear();
      edges_.clear();
      num_edges_ = 0;
      adj_map.clear();
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
        return graph_->node(idx_);
    }

    //If valid node, increment the node index and return it
    NodeIterator& operator++() {
        idx_++;
        return *this;
    }

    //Ensure that the graphs are the same and the node indices are the same
    bool operator==(const NodeIterator& x) const {
        return (x.idx_ == idx_ && (x.graph_ == graph_));
    }

   private:
    friend class Graph<V,E>;
    //Private Constructor of NodeIterator Class

    // HW1 #2: YOUR CODE HERE

    //Pointer to graph and index
    graph_type* graph_;
    size_type idx_;
    NodeIterator(const graph_type* graph, size_type idx)
        : graph_(const_cast<graph_type*>(graph)), idx_(idx) {
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
      return NodeIterator(this,this->num_nodes());
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
        return Edge{graph_, node_, (iter_->first)};
    }

    IncidentIterator& operator++() {
        iter_++;
        return *this;
    }
    bool operator==(const IncidentIterator& x) const {
        return node_ == x.node_ && iter_ == x.iter_ && graph_ == x.graph_;
    }

   private:
    friend class Graph<V,E>;
    // HW1 #3: YOUR CODE HERE
    //Private Constructor of Incident Iterator
    graph_type* graph_; // Pointer to the graph so we can use edge()
    size_type node_; //id of the node that we're looking for incidence with
    std::unordered_map<size_type,size_type>::iterator iter_;

    IncidentIterator(const graph_type* graph, size_type node,
                    std::unordered_map<size_type,size_type>::iterator iter)
        : graph_(const_cast<graph_type*>(graph)), node_(node),
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
        ++uid_;
        return *this;
    }
    bool operator==(const EdgeIterator& x) const {
        return uid_ == x.uid_ && graph_ == x.graph_;
    }

   private:
    friend class Graph<V,E>;
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
      return EdgeIterator(this,num_edges());
  }
//Compare(const Point& p) : p {} {
    //return (n1()-p > node2()-p)

 private:
     struct internal_node {
         Point position;
         node_value_type value;
         size_type index;
         bool removed;
     };

     //Indices of two nodes that the edges touches
     struct internal_edge {
         size_type node1;
         size_type node2;
         size_type uid_;
         edge_value_type value;
     };
  // HW0: YOUR CODE HERE

  // Vector of internal_node structs //
  std::vector<internal_node> nodes_;

  //Vector of internal_edge structs //
  std::vector<internal_edge> edges_;

  //Index to unique id
  std::vector<size_type> i2u_;
  std::vector<size_type> e2u_;
  //Adjacency map that has node has key and set of incident edges as the value //
  adjacency_map_type adj_map; //adjacency list

  size_type num_edges_;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
