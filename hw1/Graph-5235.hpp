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
template <typename V>
class Graph {
 // private:
 //  // HW0: YOUR CODE HERE
 //  struct node_element{
 //    Point point;
 //    size_type uid;
 //  }

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
  class Node:private totally_ordered<Node> {
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
      return graph_->node_vec[uid_].point;
      // return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->node_vec[uid_].uid;
      // return size_type(-1);
    }

    // HW1: YOUR CODE HERE
	
    // Supply definitions AND SPECIFICATIONS for:
  /** return the value of this node value  */
     node_value_type& value(){
      return graph_->node_vec[uid_].val;
    };
    /** return the value of this node value  */
      const node_value_type& value() const{
      return graph_->node_vec[uid_].val;
    };
  /** return the number of nodes adjacent to this node */
     size_type degree() const{
      return graph_->adj_vec[uid_].size();
     };
  /** return an incidentIterator pointing to the fist edge adjacent to this node */
     incident_iterator edge_begin() const{
      return incident_iterator(graph_,uid_,0);
     };
  /** return an incidentIterator pointing to the end edge adjacent to this node */
     incident_iterator edge_end() const{
      return incident_iterator(graph_,uid_,degree());
     };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
     return graph_ == n.graph_ and uid_==n.uid_;

      //  (void) n;          // Quiet compiler warning
      //return false;
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
      return (graph_<n.graph_)or(graph_==n.graph_ and uid_<n.uid_);
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;
    // private constructor
    Node(const Graph* graph,size_type uid)
        :graph_(const_cast<Graph*>(graph)),uid_(uid) {
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //return size_;  //define at the private end 
    return node_vec.size();
    // return 0;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes() ???????????
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
     node_element new_node;
     new_node.point = position;
     new_node.uid = num_nodes();
     new_node.val = value;
     node_vec.push_back(new_node);
     size_ = node_vec.size();
    // adj_element new_adj;
    // new_adj.adj_uid =-1;
    // new_adj.adj_eid =-1;
     std::vector<adj_element> temp;
     adj_vec.push_back(temp);
     return Node(this,node_vec.size()-1);
    // (void) position;      // Quiet compiler warning
    // return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.uid_<=size();
    // (void) n;            // Quiet compiler warning
    // return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this,i);
    // (void) i;             // Quiet compiler warning
    // return Node();        // Invalid node
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,graph_->edge_vec[eid_].n1);
      // return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,graph_->edge_vec[eid_].n2);
      // return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((graph_==e.graph_)and(eid_==e.eid_));
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
     if(graph_==e.graph_){
       if(eid_<e.eid_) return true;
     }
     return false;
     // (void) e;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_;
    size_type eid_;
    size_type uid1_;
    size_type uid2_;
    //Edge(const Graph* graph,size_type eid) :graph_(const_cast<Graph*>(graph)),eid_(eid) {}
    Edge(const Graph* graph,size_type eid, size_type uid1, size_type uid2):graph_(const_cast<Graph*>(graph)),eid_(eid),uid1_(uid1),uid2_(uid2){}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    //return edge_size_;
    return edge_vec.size();
    // return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const{
    // HW0: YOUR CODE HERE
    return Edge(this,i,edge_vec[i].n1,edge_vec[i].n2);
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    // HW0: YOUR CODE HERE
     for (unsigned int i =0;i<adj_vec[a.uid_].size();i++){
       if(b.uid_==adj_vec[a.uid_][i].adj_uid){
       return true; 
       }
      return false;
     }
   // for(auto& e:edge_vec){
    // if((e.n1==a.index()&&e.n2==b.index())||(e.n2==a.index()&&e.n1==b.index()))
     // return true;
     // }
     //return false;
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;
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
 //   if(has_edge(a,b)){
  //    for(auto& e:edge_vec){
   //     if((e.n1 == a.index() && e.n2 == b.index()) || (e.n2 == a.index() && e.n1 == b.index())) 
    //      return Edge(this,adj_vec[a.uid_][e[]].eid,a.uid_,b.uid_);
     //   }
   // }
    for(unsigned int i = 0; i<adj_vec[a.uid_].size();i++){
      if(has_edge(a,b)) return Edge(this,adj_vec[a.uid_][i].adj_eid,a.uid_,b.uid_);
      }
    edge_element new_edge;
    new_edge.n1 = a.index();
    new_edge.n2 = b.index();
    new_edge.eid = num_edges();
    edge_vec.push_back(new_edge);
    adj_element new_adj1;
    adj_element new_adj2;
    new_adj1.adj_uid = b.index();
    new_adj1.adj_eid = edge_vec.size()-1;
    new_adj2.adj_uid = a.index();
    new_adj2.adj_eid = edge_vec.size()-1;
    adj_vec[a.uid_].push_back(new_adj1);
    adj_vec[b.uid_].push_back(new_adj2);
    return Edge(this,edge_vec.size()-1,new_edge.n1,new_edge.n2);
    // (void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edge_vec.clear();
    node_vec.clear();
    adj_vec.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator> {
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
  /** return the Node NodeIterator pointing to */
    Node operator*() const{
      return Node(graph_,index_p);
    }
   /** return the NodeIterator pointing to next node or nullptr */
    NodeIterator& operator++(){
     index_p++;
     return *this;
    }
    /** return true if this node and the node pointed by nodeIterator are the same */
    bool operator==(const NodeIterator& nodeIterator)const{
     return (nodeIterator.graph_==graph_)&&(nodeIterator.index_p==index_p); 
    }
    

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph *graph_;
    size_type index_p;
     NodeIterator(const Graph* graph,size_type i):graph_(graph),index_p(i){
    }


    };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  /**return NodeIterator pointing to the fisrt Node in graph*/ 
    node_iterator node_begin() const{
    return NodeIterator(this,0);
     }
  // node_iterator node_end() const
/** return NodeIterator pointing to the last Node in graph**/
    node_iterator node_end() const{
    return NodeIterator(this,node_vec.size());
    }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
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
/** return the edge incidentIterator pointing to   */
    Edge operator*() const{
       return Edge(graph_,graph_->adj_vec[node_id_][inc_id_].adj_eid,node_id_,graph_->adj_vec[node_id_][inc_id_].adj_uid);
    }
/** return IncidentIterator pointing to the next adjacent node or null_ptr*/
    IncidentIterator& operator++(){
     inc_id_++;
     return *this;
    }
/** return true if this node is the same as the adjacent node pointed at   */
     bool operator==(const IncidentIterator& iit) const{
           return((iit.graph_==graph_)&&(iit.node_id_==node_id_)&&(iit.inc_id_==inc_id_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_;
    size_type inc_id_;
    IncidentIterator(Graph *graph,size_type node_id,size_type inc_id):graph_(graph),node_id_(node_id),inc_id_(inc_id){}


  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
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
    /**return Edge this pointing to   */
     Edge operator*() const{
     return Edge(graph_,edge_id_,graph_->edge_vec[edge_id_].n1,graph_->edge_vec[edge_id_].n2);
    }
    /** return EdgeIterator pointing to next edge or nullptr */
     EdgeIterator& operator++(){
     edge_id_++;
     return *this;    
     }
    /** return true if this edge is the same as edge pointing to*/
     bool operator==(const EdgeIterator& eit) const{
     return (graph_==eit.graph_)&&(eit.edge_id_==edge_id_);
     }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    size_type edge_id_;
    EdgeIterator(const Graph *graph,size_type edge_id):graph_(graph),edge_id_(edge_id){} 
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
/**return EdgeIterator pointing to the first edge */
   EdgeIterator edge_begin() const{
    return  EdgeIterator(this,0);
   }
 /** return Edge Iterator pointing to the last edge*/
   EdgeIterator edge_end() const{
   return EdgeIterator(this,edge_vec.size());    
}

 private:
  // HW0: YOUR CODE HERE
   struct node_element{
    Point point;
    size_type uid;
    node_value_type val;
   };
   struct edge_element{
     size_type n1;
     size_type n2;
     size_type eid;
   };
   struct adj_element{
     size_type adj_uid;
     size_type adj_eid;
   };

   size_type size_;
   size_type edge_size_;
   std::vector<node_element> node_vec;
   std::vector<edge_element> edge_vec;
   std::vector<std::vector<adj_element>> adj_vec;  
 //std::map<std::map<>>
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
