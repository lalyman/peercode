#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <sstream>
#include <string>
#include <stdlib.h>

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
 //struct to save node info
  struct NodeInfo{
      Point position_;
     
      V value_ = 0;// default value to zero
  };
  using node_info=NodeInfo;
  std::vector<node_info> internal_nodes; // a vector to save node information
  unsigned size_,edge_size;
  std::vector<std::vector<std::pair<unsigned,unsigned>>> connectivity; // a 3-dim vector to save edges and connectivities
  //std::map<unsigned, std::unordered_set<unsigned>> connectivity ;
 // connectivity.insert(std::pair<unsigned, unordered_set<unsigned> (0,{}));
  
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  using node_value_type = V;
  /** Type of this graph. */
  using graph_type = Graph<V>;

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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph():size_(0),edge_size(0){}
    // HW0: YOUR CODE HERE 

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
    Node():node_id(0), graph_(){}
      // HW0: YOUR CODE HERE
    /** Return reference to this node value*/
    node_value_type& value(){
       return graph_->internal_nodes[node_id].value_;
    }
    /** Return this node's value.*/
    const  node_value_type& value() const{
       return graph_->internal_nodes[node_id].value_;
    }
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      
      return graph_->internal_nodes[node_id].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
	return node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
      /** Return number of edges connected to this node*/
      size_type degree() const{
         return graph_->connectivity[node_id].size() - 1;
      }


    /** Return the iterator to the first edge connected to this node.
   * @post this node.node_id = (result.ptr)->node_1
   * @post conncting node.node_id = (result.ptr)->node_2
   * @post connectivity[this node.node_id][connecing node.node_id].second = result.ptr
   * Complexity: O(1) amortized operations.
   */

    // incident_iterator edge_begin() const;
      incident_iterator edge_begin() const{
          
            incident_iterator e_begin(graph_->connectivity[node_id][1].second,graph_,this);
           
         return e_begin; 
     }

    /** Return the iterator to the last edge connected to this node.
   * @post this node.node_id = (result.ptr)->node_1
   * @post conncting node.node_id = (result.ptr)->node_2
   * @post connectivity[this node.node_id][connecing node.node_id].second = result.ptr
   * Complexity: O(num(connecting edges to this node))
   */
    // incident_iterator edge_end() const;
      incident_iterator edge_end() const{
          size_type csize = graph_->connectivity[node_id].size();
          
          size_type last_id = graph_->connectivity[node_id][csize-1].second;
          
         incident_iterator e_end(last_id,graph_,this); 
         return e_end; 
  
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
       */
    bool operator==(const Node& n) const {
     
       return (n.index()==this->index() && n.graph_==this->graph_);
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
      return (n.index() < this->index() && n.graph_ == this->graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type node_id;
    graph_type* graph_;
    Node(size_type node_id_ ,const graph_type* graph):node_id(node_id_),graph_(const_cast<graph_type*>(graph)){}
    Node add_node(const Point&);
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
  	
    // HW0: YOUR CODE HERE
    return size_;
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
       //incrementing size
       size_++;
       size_type new_node_id = size_ - 1;
       //pushing back node_info to the internal_nodes
       internal_nodes.push_back(node_info());
       
       internal_nodes[new_node_id].position_ = position;
       //pre-allocating memory for connectivity vector
       std::vector<std::pair<size_type,size_type>> temp = {std::make_pair(new_node_id,0)};
       connectivity.push_back(temp);

       node_type new_node(new_node_id,this);
    return new_node;        // Invalid node
  }
 Node add_node(const Point& position,const node_value_type value)
 {
	size_++;
	size_type new_node_id = size_ -1;
	internal_nodes.push_back(node_info());
        //std::cout << "pushing back node info with value" <<std::endl;
        internal_nodes[new_node_id].position_ = position;
        internal_nodes[new_node_id].value_ = value;
        std::vector<std::pair<size_type,size_type>> temp = {std::make_pair(new_node_id,0)};
        connectivity.push_back(temp);
	node_type new_node(new_node_id,this);
	return new_node;
 }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
      return (n.node_id < size() && n.graph_ == this->graph_);
    }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
        assert(i < size());
        node_type new_node(i,this);
      return new_node;
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(): edge_id(0),graph_e(){}
      // HW0: YOUR CODE HERE

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      node_type node1_(node1_id,graph_e);
      return node1_;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
     	node_type node2_(node2_id,graph_e);
      return node2_;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	return (this->node2() == e.node2() && this->node1() == e.node1() && e.graph_e == this->graph_e);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
       return (this->edge_id < e.edge_id && e.graph_e == this->graph_e);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type node1_id, node2_id, edge_id;
    graph_type* graph_e;
    Edge(size_type edge_id_ ,const graph_type* graph_e_)
        :edge_id(edge_id_),graph_e(const_cast<graph_type*>(graph_e_)){}
    Edge edge(const Node&,const Node&);
    // HW0: YOUR CODE HERE
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
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    edge_type new_edge(i,this);
    //searchin for edge with id =i
    for (size_type j=0; j < size(); ++j){
        //searching through each node
        for(size_type ii=1 ; ii < connectivity[j].size(); ++ii){
           if(connectivity[j][ii].second == i){
               //checking which node has smaller id 
               if(j < connectivity[j][ii].first){
                   new_edge.node1_id = j; 
                   new_edge.node2_id = connectivity[j][ii].first;  
                }
                else{
                   new_edge.node2_id = j; 
                   new_edge.node1_id = connectivity[j][ii].first;
                }
            }
       }
    }
    
    return new_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	assert(a.node_id < size() && b.node_id < size() && a.graph_ == b.graph_);
        size_type n1,n2;
        if (a.node_id < b.node_id) {n1=a.node_id; n2=b.node_id;} 
	else	{n1=b.node_id; n2=a.node_id;}
        for (size_type j=1 ; j < connectivity[n1].size();++j){ 
	    if (connectivity[n1][j].first ==n2)  return true;
         }
   return false;		
       // Quiet compiler warning
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
	assert(a.node_id < size() && b.node_id < size() && a.graph_ == b.graph_);
        size_type n1,n2;
        if (a < b) {n1 = a.node_id; n2= b.node_id;}
        else{ n2= a.node_id; n1= b.node_id;}
        size_type new_edge_id = edge_size - 1;
        if(!has_edge(a,b)){
		edge_size++;
                new_edge_id++;
                
                connectivity[n1].push_back(std::make_pair(n2,new_edge_id));
                connectivity[n2].push_back(std::make_pair(n1,new_edge_id));
	}
        
	edge_type new_edge(new_edge_id,this);
        new_edge.node1_id = n1;
        new_edge.node2_id = n2;
    return new_edge;        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_type edge_(0,this);
    node_type node_(0,this);
    size_ = node_.node_id;
    edge_size = edge_.edge_id;
    internal_nodes.erase(internal_nodes.begin(),internal_nodes.end());
    connectivity.erase(connectivity.begin(),connectivity.end());
    
    // HW0: YOUR CODE HERE
 }
//
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator():ptr_id(0),graph_ptr() {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the node which the iterator is pointing to.
   * @post  (result) = this node
   * 
   * Complexity: O(1) amortized operations.
   */
    // Node operator*() consT
     Node operator*() const{
	node_type this_node(ptr_id,graph_ptr);
        return this_node;
     }
    // NodeIterator& operator++()
     /** Return the refernce iterator to the next node.
   * 
   * Complexity: O(1) amortized operations.
   */
     node_iterator& operator++(){

        //std::cout << "ptr->node_id"<< ptr->node_id<<std::endl;
        
         ++(ptr_id);
        //std::cout << "ptr->node_id:" << ptr_id<<std::endl;
     return *(this);
     }
    // bool operator==(const NodeIterator&) const
    /** check whether two pointers points to the same node .
   * 
   * Complexity: O(1) amortized operations.
   */
    bool operator ==(const NodeIterator& node_iter1) const {
        return (this-> ptr_id == node_iter1.ptr_id && this->graph_ptr == node_iter1.graph_ptr );
    }
   private:
    friend class Graph;
    size_type ptr_id;
    graph_type* graph_ptr;
    NodeIterator(size_type ptr_id_ ,const graph_type* graph_ptr_)
        :ptr_id(ptr_id_),graph_ptr(const_cast<graph_type*>(graph_ptr_)){}
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
    /** Return the iterator to the firs node.
   * @post *(result.ptr) = first node
   * Complexity: O(1) amortized operations.
   */
    node_iterator node_begin() const{ 
        node_iterator nbegin(0,this);
        
     return nbegin;
   }
  // node_iterator node_end() const
   /** Return the iterator to the last node.
   * @post *(result.ptr) = last node
   * Complexity: O(1) amortized operations.
   */
   node_iterator node_end() const{
     
     node_iterator nend(size()-1,this);
     return nend;
   }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator():ptr_id_e(0),graph_ptr_e(),node_ptr() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
     /** Return the edge pointed by the incidentiterator.
   * 
   * Complexity: O(1) amortized operations.
   */
     Edge operator*() const{
         //initializing edge with this node.
         edge_type this_edge(ptr_id_e,graph_ptr_e);
         this_edge.node1_id = node_ptr->index();
         //initializing node1
         node_type n1(this_edge.node1_id,graph_ptr_e);
         size_type csize = graph_ptr_e->connectivity[n1.index()].size();
         //searching through edges connected to the node to find node2
         for(size_type i =1; i< csize;++i){
             if(graph_ptr_e->connectivity[n1.index()][i].second == ptr_id_e){
                   this_edge.node2_id = graph_ptr_e->connectivity[n1.index()][i].first;
              }
         }
         return this_edge;
     }
    // IncidentIterator& operator++()
     /** Return the iterator the next edge connected to the node.
   * 
   * Complexity: O(num of edges connected to node1) amortized operations.
   */
     IncidentIterator& operator++(){
       size_type id1  = node_ptr->index();
           //going through correspondinh edges to check the id
           for(size_type j = 1; j< graph_ptr_e->connectivity[id1].size() - 1;j++){
             if (graph_ptr_e->connectivity[id1][j].second == ptr_id_e )
                {//if edge id is found, we shift to the next
                 ptr_id_e = graph_ptr_e->connectivity[id1][j+1].second;
                 return *(this);
                 }
          }
     return *(this);
     }
    // bool operator==(const IncidentIterator&) const

     /** Check whether  the pointer points to the same edge
   * Complexity: O(1) amortized operations.
   */
     bool operator==(const IncidentIterator& edge_iter) const{
         return (this->ptr_id_e ==  edge_iter.ptr_id_e && this->graph_ptr_e == edge_iter.graph_ptr_e);

    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
   
    size_type ptr_id_e;
    graph_type* graph_ptr_e;
    node_type* node_ptr;
     IncidentIterator(size_type ptr_id_e_ ,const graph_type* graph_ptr_e_,const node_type* node_ptr_)
        :ptr_id_e(ptr_id_e_),graph_ptr_e(const_cast<graph_type*>(graph_ptr_e_)) ,node_ptr(const_cast<node_type*>(node_ptr_)){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator():ptr_id_e(0),graph_ptr_e() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
      /** Return the edge pointed by the edgeiterator.
   * 
   * Complexity: O(1) amortized operations.
   */
      Edge operator*() const{
          edge_type new_edge(ptr_id_e,graph_ptr_e);
         return new_edge;
      }
    // EdgeIterator& operator++()
     /** Return the iterator the next edge.
   * 
   * Complexity: O(num of edges connected to node1) amortized operations.
   */

     EdgeIterator& operator++(){
         ++(ptr_id_e);
         return *(this);
     }
    // bool operator==(const EdgeIterator&) const
     /** Check whether  the pointer points to the same edge
   * Complexity: O(1) amortized operations.
   */
      bool operator==(const EdgeIterator& edge_iter) const{
         return  (this->ptr_id_e == edge_iter.ptr_id_e && this->graph_ptr_e == edge_iter.graph_ptr_e);
     }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
     size_type ptr_id_e;
     graph_type* graph_ptr_e;
     EdgeIterator(size_type ptr_id_e_ ,const graph_type* graph_ptr_e_)
        :ptr_id_e(ptr_id_e_),graph_ptr_e(const_cast<graph_type*>(graph_ptr_e_)){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
   /** Return the iterator to thefirst edge.
   * @post *(result.ptr_id_e) = first edge
   * Complexity: O(1) amortized operations.
   */
    edge_iterator edge_begin() const{
        edge_iterator first_edge(0,this);
     return first_edge;
    }
  // edge_iterator edge_end() const
    /** Return the iterator to the last edge.
   * @post *(result.ptr_id_e) = last edge
   * Complexity: O(1) amortized operations.
   */

    edge_iterator edge_end() const{
    edge_iterator last_edge(edge_size,this);
    return last_edge;
   }

 private:
  
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
