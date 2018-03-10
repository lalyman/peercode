#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include<iterator>
#include<list>
#include <functional>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "VectorList.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {

 public:

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */

   using size_type = unsigned;

 private:
 struct AdjInfo;
 struct InternalEdge;
 struct InternalNode;

 std::vector<InternalNode> internal_nodes;
 VectorList<size_type> node_mask;
 std::list<InternalEdge> internal_edges;
 std::list<std::list<AdjInfo>> adjacency;

 public:
 //
 // PUBLIC TYPE DEFINITIONS
 //

 /** Type of this graph. */

  using adj_info_iter = typename std::list<AdjInfo>::iterator;
  using e_iter = typename std::list<InternalEdge>::iterator;
  using adj_iter = typename std::list<std::list<AdjInfo>>::iterator;

  using graph_type = Graph;
  typedef V node_value_type;
  typedef E edge_value_type;
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
     /* Forward the request to graph_, which looks up the
      * appropriate internal node's position
      */
      return graph_->fetch_node(id_).p;
    }

    Point& position() {
      return graph_->fetch_node(id_).p;
    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    /** Return a non-const reference to the value of a node
    * @pre 0 <= id_ < graph.size()
    * @return A non-const reference to the value of the node
    * 
    * Complexity: O(1)
    */
    node_value_type& value(){
      return graph_->fetch_node(id_).value;

    }

    /** Return a const reference to the value of a node
    * @pre 0 <= id_ < graph.size()
    * @return A const reference to the value of the node
    *
    * Complexity: O(1)
    **/

    const node_value_type& value() const{
      return graph_->fetch_node(id_).value;
    }

    /** Return the degree of a node
    * @pre 0<= id_ < graph.size()
    *
    * Complexity: O(1)
    **/
    
    size_type degree() const{
      return (*(graph_->internal_nodes[id_].adj_ref)).size();
    }

    /** Return an beginning incident_iterator 
    * @pre Node is either valid or one past the end
    * @post If Node is valid, *result = internal_edges[(*this).get_id()][0]
    *
    * @return A beginning incident_iterator
    *
    **/

    incident_iterator edge_begin() const{
      auto i_node = graph_->fetch_node(id_);
      
      return IncidentIterator((*i_node.adj_ref).begin());
    }

    /** Return an ending incident_iterator
    * @pre Node is valid
    * @post result.id_ == this->degree()
    * @return An ending incident_iterator
    */

   incident_iterator edge_end() const {
      auto i_node = graph_->fetch_node(id_);
     return IncidentIterator((*i_node.adj_ref).end());
   }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.get_graph() == graph_) && (n.get_id() == id_);
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
      auto less = std::less<const graph_type*>();
      return less(graph_, n.get_graph()) || (graph_ == n.get_graph() && id_ < n.get_id());
    }

    size_type get_id() const {
      return id_;
   }

   const graph_type* get_graph() const{
     return graph_;
   }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class NodeIterator;

    // A node is characterized by its graph and its index
    graph_type* graph_; 
    size_type id_;

    Node(const graph_type* g, size_type id) : graph_(const_cast<graph_type*>(g)),
                                              id_(id) {}
  };

  size_type unmask(size_type uid) {
    return node_mask[uid];
  }

  /** Remove a node and all of its adjacent edges from the graph
  * @param[in] n the node to be removed
  *
  * @post new size() = old size() - 1 
  * @post new num_edges() = old num_edges() - degree(n)
  * @post All outstanding Node objects with id greater than @a n.get_id()
  *       are invalidated
  * @post All outstading Edge objects with id greater than the smallest id
  *       of an edge adjacent to @a n are invalidated
  * @post All edge_iterators not pointing at a removed edge remain valid
  * @post All incident_iterators not pointing at a removed edge remain valid
  * @return @a n.get_id(), which is now the id of the node previously follwowing @a n
  * Complexity: O(N) where N is the number of nodes in the graph
  **/
  size_type remove_node(const Node& n){
    size_type id = n.get_id();
    InternalNode i_node = internal_nodes[id];
    
    for(AdjInfo e : *i_node.adj_ref){
      auto n2_ = (*e.edge_).node2();
      if (n2_ == n){
        n2_ = (*e.edge_).node1();
      }

      InternalNode i_n2_ = internal_nodes[n2_.get_id()];
      //Remove the second iterator in the adjacency list that points at this edge
      (*(i_n2_.adj_ref)).erase(e.alt_edge_);

      internal_edges.erase(e.edge_);
    }
    //Remove this node's entry in the adjacency list
    adjacency.erase(i_node.adj_ref);

    //Remove this node from internal_nodes
    internal_nodes.erase(internal_nodes.begin() + id);


    //Update node_mask by decrementing all nodes still in the graph that have
    // a uid_ greater than the uid_ of the remvoed node
    for(auto it = node_mask.iter(i_node.uid_); it != node_mask.end(); ++it) {
      --((*it).get_val());
    }

    node_mask.erase(node_mask.iter(i_node.uid_));
    return id;
  }

   /** Remove a node and its adjacenct edges from the graph
    * param[in] A node_iterator pointing at the node to be removed
    * Identical to remove_node(size_type) in every respect expect its return type and input type
    * @return A node iterator pointing at the node following the removed node @a *n_it
    **/
   node_iterator remove_node(node_iterator n_it){
     size_type id = remove_node(*n_it);
     return node_iterator(this, id);
   }

  /** Remove a node and all of its adjacent edges from the graph
  * @param[in] n the node to be removed
  *
  * @post new num_edges() = old num_edges() - 1
  *
  * @post All outstading Edge objects with id greater than the id of the
  *       removed edge are invalidated
  * @post All edge_iterators not pointing at the removed edge remain valid
  * @post All incident_iterators not pointing at the removed edge remain valid
  * @return 1 if the edge was in the graph and then removed. 0 otherwise. 
  * Complexity: O(d) where d = @a a.degree(). Worst case: O(size())
  **/

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return internal_nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val (optional) The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */

  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    InternalNode n = {position, val, static_cast<size_type>(node_mask.underlying_size())};
    internal_nodes.push_back(n);
    node_mask.push_back(size()-1);

    std::list<AdjInfo> new_edges {};
    adjacency.push_back(new_edges);

    internal_nodes[internal_nodes.size() - 1].adj_ref = std::prev(adjacency.end(), 1);
    
    return Node(this, this->size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this) && (n.get_id() < this->size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i && i < this->size());
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return (*iter_).node1();
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return (*iter_).node2();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // First check that the two edges are from the same graph
      if(e.graph_ != graph_) {
        return false;
      }
      else {
        Node n1_ = node1();
        Node n2_ = node2();
        
        //Edges are undirected, so check both orientations
        bool ans = (n1_ == e.node1() && n2_ == e.node2());
        ans = ans || (n1_ == e.node2() && n2_ == e.node1());
        return ans; 
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1() == e.node1()) {
        return node2() < e.node2();
      }
      return node1() < e.node1();
    }

    const graph_type* get_graph() const{
      return graph_;
    }

   const edge_value_type& value() const {
     return (*iter_).value;
   }

   edge_value_type& value() {
     return (*iter_).value;
   }

   double length() const {
     return norm(this->node1().position() - this->node2().position());
   }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;

    // Edges are characterized by their graph, their first node
    // and their index in the adjacency list at that node

    const graph_type* graph_;
    e_iter iter_;    

    Edge(e_iter iter) :
        iter_(iter)  {graph_ = (*iter_).node1().get_graph();}
   };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
      return internal_edges.size();
    }

    Edge edge(size_type id) {
      auto iter_ = this->edge_begin();
      for(size_type i = 0; i < id; ++i) {
        ++iter_;
      }
      return *iter_;
    }
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    bool has_edge(const Node& a, const Node& b) const {
      // Iterate over all the internal edges, checking each to see
      // if any has endpoints @a a and @a b
      for(auto e : internal_edges){
        if(e.has_endpoints(a,b)){
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
     // Similar to has_edge, iterate over internal_edges
     // returning a proxy associated with the InternalEdge with
     // endpoints a and b if such an InternalEdge exists

     auto n1_ = fetch_node(a.get_id());
     auto n2_ = fetch_node(b.get_id());

     for(AdjInfo e : *n1_.adj_ref) {
        InternalEdge i_edge = *e.edge_;
        if(i_edge.has_endpoints(a,b)) { 
          return Edge(e.edge_);
        }
      }
     
     internal_edges.push_back(InternalEdge(this, n1_.uid_,n2_.uid_));

     e_iter e_tail = std::prev(internal_edges.end(), 1);
 
     AdjInfo adj1 {e_tail};
     AdjInfo adj2 {e_tail};

     (*n1_.adj_ref).push_back(adj1);
     (*n2_.adj_ref).push_back(adj2);
     
     (*n1_.adj_ref).back().alt_edge_ = std::prev((*n2_.adj_ref).end(), 1);
     (*n2_.adj_ref).back().alt_edge_ = std::prev((*n1_.adj_ref).end(), 1);

     return Edge(e_tail);
    }

  size_type remove_edge(const Node& a, const Node& b){
    InternalNode i_n1 = internal_nodes[a.get_id()];
    InternalNode i_n2 = internal_nodes[b.get_id()];

    //Look for the edge in the adjacency list
    //then erase it from the adjacency list (in both places)
    //and from internal_edges
    for(auto it = (*i_n1.adj_ref).begin(); it != (*i_n1.adj_ref).end(); ++it) {
      auto e = *it;
      if ((*e.edge_).has_endpoints(a,b)) {
        internal_edges.erase(e.edge_);
        (*(i_n2.adj_ref)).erase(e.alt_edge_);
        (*(i_n1.adj_ref)).erase(it);
        return 1;
      }
    }
    return 0;
  }

  /** Remove an edge from the graph
   * @param[in] The Edge to be removed
   * Identical to remove_edge(const Node& a, const Node& b) in every respect
   * expect for its input type
   **/
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph and return an edge_iterator to the next edge
   * @param[in] An edge_iterator pointing to the edge to be removed
   * Identical to remove_edge(const Node& a, const Node& b) except for its input and return types
   * @return An edge_iterator pointing at the next edge after @a *e_it pre-removal.
   **/
  edge_iterator remove_edge(edge_iterator e_it) {
    //Very similar to remove_edge(const Node&, const Node&)
    //but we need to keep track of the iterator to return a new edge_iterator
    Node a = (*e_it).node1();
    Node b = (*e_it).node2();

    InternalNode i_n1 = internal_nodes[a.get_id()];
    InternalNode i_n2 = internal_nodes[b.get_id()];

    for(auto it = (*i_n1.adj_ref).begin(); it != (*i_n1.adj_ref).end(); ++it) {
      auto e = *it;
      if ((*e.edge_).has_endpoints(a,b)) {
        (*(i_n2.adj_ref)).erase(e.alt_edge_);
        auto next = internal_edges.erase(e.edge_);
        (*(i_n1.adj_ref)).erase(it);
        return edge_iterator(next);
      }
    }

  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    void clear() {
      internal_edges.clear();
      internal_nodes.clear();
      adjacency.clear();
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

    /** Return the node pointed to by this iterator
    * @pre NodeIterator is not one past the end
    *
    * @return A Node pointed to by this iterator
    */
    
    Node operator*() const{
      return Node(graph_, id_);
    }

    /** Increment the NodeIterator 
    * @pre NodeIterator is not one past the end
    *
    * @return The incremented NodeIterator
    */
    node_iterator& operator++(){
      ++id_;
      return *this;
    }

    /** Check equality with another NodeIterator
    *
    * @return true if any only if dereferencing gives equal nodes
    */

    bool operator==(const node_iterator& ni) const{
      return (graph_ == ni.get_graph()) && (id_ == ni.get_id());
    }

   const graph_type* get_graph() const {
     return graph_;
   }

   size_type get_id() const {
     return id_;
   }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type id_; 
    NodeIterator(const graph_type* g, size_type id) : graph_(const_cast<graph_type*>(g)), id_(id) {}

  };

  /** Return a beginning node_iterator
   * 
   * @post *return == this->internal_nodes[0]
   * @return a beginning node_iterator
   */
 
  node_iterator node_begin() const {
     return NodeIterator(this, 0);
  }
  
  /** Return an ending node_iterator
   *
   * @return an ending node_iterator
   */

  node_iterator node_end() const{
    return NodeIterator(this, this->size());
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

    /** Return an Edge proxy for where this iterator points
     * @pre IncidentIterator is not one past the end
     *
     * @return An Edge proxy for where this iterator points 
     */

    Edge operator*() const {
      return Edge((*iter_).edge_);
    }

    /** Increment the IncidentIterator
     * @pre IncidentIterator is not one past the end
     * @post This iterator points to the next edge in the adjacency list
     * @return The incremented iterator
     */

    incident_iterator& operator++(){
      ++iter_;
      return *this;
    }
    
    /** Check equality of IncidentIterators
     *
     * @return true if and only if the two IncidentIterators have
     *          the same node_ and same id_
     */

    bool operator==(const incident_iterator& iit) const {
      return (iter_ == iit.get_iter());        
    }

   const adj_info_iter& get_iter() const {
     return iter_;
   }

   private:
    friend class Graph;
    adj_info_iter iter_;

    IncidentIterator(adj_info_iter iter) : iter_(iter) {}
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


     /** Return the edge pointed to
      * @pre inc_iter_ is not one off the end
      *
      * @return The Edge pointed to by the iterator
      */
     
     Edge operator*() const {
       return Edge(iter_);
     }

     /** Increment the EdgeIterator
      * @pre EdgeIterator is not one past the end
      * @brief This increment skips Edges e where e.node1() <= e.node2()   
      *
      * @return The incremented EdgeIterator
      */

     EdgeIterator& operator++() {
         ++iter_;
       return *this;
     }

     /** Check equality of EdgeIterators
      *
      *
      * @return If either iterator is one past the end, true iff their niter_'s are equal
      *         otherwise, true iff both internal iterators are equal
      */

     bool operator==(const EdgeIterator& eiter) const {
       return iter_ == eiter.iter_;
     }

     const e_iter& get_iter() const {
       return iter_;
     }
     
   private:
    friend class Graph;
    e_iter iter_;
    EdgeIterator(e_iter iter) : iter_(iter) {}

  };

   /** Return a beginning edge_iterator
    * @pre graph contains at least one edge
    * @post *return_value == internal_edges[0][0]
    * @return A beginning edge_iterator
    */

   edge_iterator edge_begin() {
     return EdgeIterator(internal_edges.begin()); 
   }

   /** Return an ending edge_iterator
    * @return An ending edge_iterator
    */

  edge_iterator edge_end() {
    return EdgeIterator(internal_edges.end());
  }

 
 private:

   struct AdjInfo {
         
     e_iter edge_;
     adj_info_iter alt_edge_;
     AdjInfo(e_iter edge) : edge_(edge) {}
   };


   struct InternalNode{
     Point p;
     node_value_type value;
     adj_iter adj_ref;
     size_type uid_;

     InternalNode(Point pos, node_value_type val, size_type uid) : p(pos), value(val), uid_(uid) {}

   };
   struct InternalEdge{
     size_type uid1_;
     size_type uid2_;
     edge_value_type value;
     graph_type* graph_;

     InternalEdge(graph_type* graph, size_type uid1, size_type uid2) : uid1_(uid1), uid2_(uid2), graph_(graph) {}

     InternalEdge(graph_type* graph, size_type uid1, size_type uid2, edge_value_type val) : InternalEdge(graph, uid1, uid2), value(val) {}

     Node node1() const {
       return Node(graph_, graph_->unmask(uid1_));
     }

     Node node2() const {
       return Node(graph_, graph_->unmask(uid2_));
     }

     bool has_endpoints(const Node& a, const Node& b){
       return ((node1() == a) && (node2() == b)) || ((node1() == b) && (node2() == a));
     }
   };

   const InternalNode& fetch_node(size_type id) const{
     assert (0 <= id && id < this->size());
     return internal_nodes[id];
   }

   InternalNode& fetch_node(size_type id) {
     assert (0 <= id && id < this->size());
     return internal_nodes[id];
   }

   const InternalEdge& fetch_edge(const Edge& e) const{
     for(auto adj : *(e.node1()).adj_ref) {
       if (*adj.edge_.has_endpoints(e.node1(),e.node2())) {
         return e;
       }
     }
   }
   
   InternalEdge& fetch_edge(const Edge& e) {
     for(auto adj : *(e.node1()).adj_ref) {
       if (*adj.edge_.has_endpoints(e.node1(),e.node2())) {
         return e;
       }
     }
   }
};
#endif // CME212_GRAPH_HPP
