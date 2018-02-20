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
template <typename V,typename E>
class Graph {
 private:
 //struct to save node info
  struct NodeInfo{
      Point position_;
      V value_;// default value to zero
  };
  struct EdgeInfo{
     unsigned edge_id;
      E edge_values_;
  };
  using edge_info = EdgeInfo;
  using node_info = NodeInfo;

  std::vector<node_info> internal_nodes; // a vector to save node information
  std::vector<std::vector<std::pair<unsigned,edge_info>>> connectivity; // a 3-dim vector to save edges and connectivities
  std::vector<unsigned> node_id_;
  std::vector<unsigned> eid_;
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  //template <typename E>  
  using node_value_type = V;
  using edge_value_type = E;
  /** Type of this graph. */
  using graph_type = Graph<V,E>;

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
    Node():idx_(0), graph_(){}
      // HW0: YOUR CODE HERE
    /** Return reference to this node value*/
    node_value_type& value(){
       return graph_->internal_nodes[graph_->node_id_[idx_]].value_;
    }
    /** Return this node's value.*/
    const  node_value_type& value() const{
       return graph_->internal_nodes[graph_->node_id_[idx_]].value_;
    }
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //std::cout << "position :"<< graph_->internal_nodes[graph_->node_id_[idx_]].position_ <<std::endl;
      return graph_->internal_nodes[graph_->node_id_[idx_]].position_;
    }
    Point& position(){
    //std::cout << "position :"<< graph_->internal_nodes[graph_->node_id_[idx_]].position_;
       
    return graph_->internal_nodes[graph_->node_id_[idx_]].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
	return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
      /** Return number of edges connected to this node*/
      size_type degree() const{
         return graph_->connectivity[graph_->node_id_[idx_]].size();
      }


    /** Return the iterator to the first edge connected to this node.
   * @post this node.node_id = (result.ptr)->node_1
   * @post conncting node.node_id = (result.ptr)->node_2
   * @post connectivity[this node.node_id][connecing node.node_id].second = result.ptr
   * Complexity: O(1) amortized operations.
   */

    // incident_iterator edge_begin() const;
      incident_iterator edge_begin() const{
            //std::cout << "edge_begin()"<< graph_->connectivity[graph_->node_id_[idx_]][0].second.edge_id<<std::endl;
            //incident_iterator e_begin(0,graph_,this);
	//incident_iterator e_begin(this->graph_->connectivity[graph_->node_id_[idx_]][0].second.edge_id,graph_,this);
        incident_iterator e_begin(0,graph_,this); 
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
          //size_type csize = graph_->connectivity[node_id].size();
          
          //size_type last_id = graph_->c:onnectivity[node_id][csize-1].second;
         //int last_id = graph_->connectivity[node_id][0].second.edge_id;
         //std::cout << "edge_size"<< graph_->edge_size<<std::endl;
         incident_iterator e_end(this->degree(),graph_,this);
         //incident_iterator e_end(graph_->num_edges(),graph_,this); 
         return e_end; 
  
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
       */
    bool operator==(const Node& n) const {
     
       return (this->index()==n.index() && n.graph_==this->graph_);
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
      //return (this->index() < n.index() && n.graph_ == this->graph_);
       if(this->graph_ == n.graph_) return (this->index() < n.index());
       else return (this->graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type idx_;
    graph_type* graph_;
    Node(size_type idx ,const graph_type* graph):idx_(idx),graph_(const_cast<graph_type*>(graph)){}
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
       size_type new_idx_ = size_ - 1;
       //pushing back node_info to the internal_nodes
       internal_nodes.push_back(node_info());
       //std::cout << "internal_nodes.size:"<< internal_nodes.size()<<std::endl;
       node_id_.push_back(internal_nodes.size()-1);
       //std::cout << "node_id:"<< internal_nodes.size()-1<<std::endl;
       V default_value{};
       //edge_info new_edge{};
       
       //E default_edge_values;
       internal_nodes[node_id_[new_idx_]].value_ = default_value;
       internal_nodes[node_id_[new_idx_]].position_ = position;
       //pre-allocating memory for connectivity vector
       std::vector<std::pair<size_type,edge_info>> temp ;
       connectivity.push_back(temp);

       //std::cout << "node"<< new_node_id<<std::endl;
       node_type new_node(new_idx_,this);
    return new_node;        // Invalid node
  }
 Node add_node(const Point& position,const node_value_type value)
 {
	size_++;
	size_type new_idx_ = size_ -1;
	internal_nodes.push_back(node_info());
        //std::cout << "node_id:"<< new_idx_<<std::endl;
        node_id_.push_back(internal_nodes.size()-1);
        //std::cout << "pushing back node info with value" <<std::endl;
        internal_nodes[node_id_[new_idx_]].position_ = position;
        internal_nodes[node_id_[new_idx_]].value_ = value;
        //edge_info new_edge{};
       
        //internal_nodes[new_node_id].edge_info.push_back({});
        std::vector<std::pair<size_type,edge_info>> temp;
        connectivity.push_back(temp);
	node_type new_node(new_idx_,this);

	return new_node;
 }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
      return (n.idx_ < size() && n.graph_ == this);
    }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
        //assert(i < size())
        //std::cout << "here"<<std::endl;
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
    Edge(): edge_id_(0),graph_e_(){}
      // HW0: YOUR CODE HERE

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      size_type idx;
      //auto p = std::find(graph_e_->node_id_.begin(),graph_e_->node_id_.end(),node1_id_);
      for(size_type i=0; i < graph_e_->size();++i){
	if(graph_e_->node_id_[i] == node1_id_)  idx = i;
	}
      node_type node1_(idx,graph_e_);
      return node1_;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type idx;
      //auto p = std::find(graph_e_->node_id_.begin(),graph_e_->node_id_.end(),node1_id_);
      for(size_type i=0; i < graph_e_->size();++i){ 
        if(graph_e_->node_id_[i]==node2_id_)  idx = i;
        }
        //auto p = std::find(graph_e_->node_id_.begin(),graph_e_->node_id_.end(),node2_id_);
     	node_type node2_(idx,graph_e_);
      return node2_;      // Invalid Node
    }
    double length() const{
      Point distance{Point(0,0,0)};
      
      distance += node1().position();
      distance -= node2().position();
      double distance_norm= norm_2(distance);
      
    return distance_norm;
    }
    edge_value_type& value(){
      size_type id1 = graph_e_->node_id_[node1().index()];
      for(size_type j=0; j <  node1().degree() ;++j){ 
         if (graph_e_->connectivity[id1][j].first == graph_e_->node_id_[node2().index()]){
               return graph_e_->connectivity[id1][j].second.edge_values_;         
        }
      }
      edge_info default_values{};
     return default_values.edge_values_;
    }
    const edge_value_type& value() const{
      size_type id1 = graph_e_->node_id_[node1().index()];
      for(size_type j=0; j< node1().degree();++j){ 
         if (graph_e_->connectivity[id1][j].first == graph_e_->node_id_[node2().index()]){
               return graph_e_->connectivity[id1][j].second.edge_values_;         
        }
      }
      edge_info default_edge_values{};
     return default_edge_values.edge_values_;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	bool comp_nodes = (this->node2() == e.node2() && this->node1() == e.node1());	
       return (comp_nodes && e.graph_e_ == this->graph_e_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
       if(e.graph_e_ == this->graph_e_) return (this->edge_id_ < e.edge_id_);
       else return (this->graph_e_ < e.graph_e_);
   }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type node1_id_, node2_id_, edge_id_;
    graph_type* graph_e_;
    Edge(size_type edge_id ,const graph_type* graph_e)
        :edge_id_(edge_id),graph_e_(const_cast<graph_type*>(graph_e)){}
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
     size_type n1,n2;
  
    for (size_type idx=0; idx < size(); ++idx){
        //searching through each node
        size_type id1 = node_id_[idx];
        for(size_type ii=0 ; ii < connectivity[id1].size(); ++ii){
            if(connectivity[id1][ii].second.edge_id == eid_[i]){
		n1 = id1;
		n2 = connectivity[id1][ii].first;
		break;
	   }
	   }
       }
		
               //checking which node has smaller id 
               		if(n1 < n2){
                 		new_edge.node1_id_ = n1; 
                  		new_edge.node2_id_ = n2;  
                	}
                	else{
                   		new_edge.node2_id_ = n1; 
                   		new_edge.node1_id_ = n2;
                	} 
			//id1 : " << new_edge.node1_id_<< " id2 : " << new_edge.node2_id_<<std::endl;
    
    
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
	//assert(a.node_id < size() && b.node_id < size() && a.graph_ == b.graph_);
        size_type n1,n2;
        if (a < b) {n1=a.index(); n2=b.index();} 
	else	{n1=b.index(); n2=a.index();}
        size_type id1 = node_id_[a.index()];
        size_type id2 = node_id_[b.index()];
        for (size_type j=0 ; j < a.degree();++j){ 
	    if (connectivity[id1][j].first == id2) {
			return true;
	   }
         }
        for (size_type j=1 ; j < b.degree();++j){ 
            if (connectivity[id2][j].first == id1)  return true;
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
	//assert(a.node_id < size() && b.node_id < size() && a.graph_ == b.graph_);
        size_type n1,n2,new_edge_id;
        if (a < b) {n1 = a.index(); n2= b.index();}
        else{ n2= a.index(); n1= b.index();}
        size_type id1 = node_id_[n1];
        size_type id2 = node_id_[n2];
	//check if the edge exists or not
        if(!has_edge(a,b)){
                	connectivity[id1].push_back(std::make_pair(id2,edge_info()));
			connectivity[id2].push_back(std::make_pair(id1,edge_info()));
		        eid_.push_back(edge_size);
			++edge_size;
            		connectivity[id1].back().second.edge_id = eid_.back();
        		connectivity[id2].back().second.edge_id = eid_.back();
			new_edge_id  = eid_.back();
	}
	//if exists find the corresponding id
	else{
		for(size_type j=0;j<connectivity[id1].size() ;++j)
		{
			if(connectivity[id1][j].first== id2) {
				for(size_type jj=0; jj< eid_.size(); ++jj){
					if(eid_[jj] == connectivity[id1][j].second.edge_id){
						new_edge_id = jj;
					}
				}
			}
		}
	}
        edge_type new_edge(new_edge_id,this);
	new_edge.node1_id_ = n1;
        new_edge.node2_id_ = n2;
	return new_edge;        
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_type edge_(0,this);
    node_type node_(0,this);
    size_ = node_.idx_;
    edge_size = edge_.edge_id_;
    internal_nodes.erase(internal_nodes.begin(),internal_nodes.end());
    connectivity.erase(connectivity.begin(),connectivity.end());
    node_id_.erase(node_id_.begin(),node_id_.end());
    // HW0: YOUR CODE HERE
 }
/**To remove node @a n with @a n.index() 
* @param[in] n  Node to remove from @a graph
* @param[out] size_ number of remaining nodes
*
* @post new @a graph has one less node, all nodes @a index() > @a n.index() are invalid aftre.
* nodes are shifted down by one. 
*       all edges connected to @a a are removed. edge_vector (@a eid_) has @a a.degree() less elements.
*       edges with index larger than min(edge indecis connecting to @a n are all invalid.  
* complexity O(num nodes connecting to this node+ number of edges connecting to adjacent nodes)
*/
  size_type remove_node(const Node& n){
        if(has_node(n)){
	
        size_type n1 = node_id_[n.index()];
	auto last = n.edge_end();
	auto it = n.edge_begin();
                while(it != last){
		auto e = *(it);
		//looping through connecting nodes to delete node @a n
		for(size_type j=0; j < e.node2().degree();++j){
                        size_type id2 = node_id_[e.node2().index()];
			if(connectivity[id2][j].first == n1){
				--edge_size;
				//deleting corresponding edge
                                for(size_type eid = 0; eid < eid_.size();++eid){
					if(eid_[eid]==connectivity[id2][j].second.edge_id){
						eid_.erase(eid_.begin()+eid);
  						break;
					}
				}
				connectivity[id2].erase(connectivity[id2].begin()+j);
                		break;
			}
         		
		}
               ++(it);
;
	}
	connectivity[n1].erase(connectivity[n1].begin(),connectivity[n1].end());
       	 node_id_.erase(node_id_.begin()+n.index());      
       size_ = size_ - 1;
       
       }
   return size_;
  }
/**To remove edge with @a a and @a b 
* @param[in] a  first node
* @param[in] b second node
* @param[out] edge_size number of remaining edges
*
* @post new @a graph has one less edge connecting @a a and @a b, 
*       @a a and @2 b has one less edge.
*       all edges with @a edge_id > @a edge_id connecting @a and @b are invalid.
*       nodes are all valid.
* element in edge_vector(@a eid_) are shifted down by one. 
* complexity O(num nodes connecting to @a a+ @a b)
*/
size_type remove_edge(const Node& a, const Node& b){
        //checking whether the edge exists or not
        if(has_edge(a,b)){
        size_type id1 = node_id_[a.index()];
        size_type id2 = node_id_[b.index()];
        //loop through each node to delete edge connecting @a a and @2 b.
	for(size_type j=0; j < a.degree() ; ++j){
		if(connectivity[id1][j].first == id2){
			//deleting edge from eid_
			for(size_type eid = 0; eid < eid_.size(); ++eid){
                                        if(eid_[eid]==connectivity[id1][j].second.edge_id){
                                                eid_.erase(eid_.begin()+eid);
                                        }
                                }
 			connectivity[id1].erase(connectivity[id1].begin()+j);
                 break;
		}
	}
	//second loop for the second node
	for(size_type j=0; j < b.degree() ; ++j){
                if(connectivity[id2][j].first == id1){
                        connectivity[id2].erase(connectivity[id2].begin()+j);		
                 break;
                }
        }

	edge_size = edge_size - 1;
        }
 return edge_size;
}
/**To remove edge with @a a and @a b 
* @param[in] e  edge to removed
* @param[out] edge_size number of remaining edges
*
* @post new @a graph has one less edge connecting @a a and @a b, 
*       all edges with @a edge_id > @a edge_id connecting @a and @b are invalid.
*       nodes are all valid.
* element in edge_vector(@a eid_) are shifted down by one. 
* complexity O(num nodes connecting to two nodes of @a e)
*/

size_type remove_edge(const Edge& e){
	Node a = e.node1();
	Node b = e.node2();
	size_type  new_edge_size = remove_edge(a,b);
	return new_edge_size;
}
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
    NodeIterator():id_ptr_(0),graph_ptr_() {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the node which the iterator is pointing to.
   * @post  (result) = this node
   * 
   * Complexity: O(1) amortized operations.
   */
    // Node operator*() consT
     Node operator*() const{
	node_type this_node(id_ptr_,graph_ptr_);
        //std::cout << "this_node_id"<< id_ptr_<<std::endl;
        return this_node;
     }
    // NodeIterator& operator++()
     /** Return the refernce iterator to the next node.
   * 
   * Complexity: O(1) amortized operations.
   */
     node_iterator& operator++(){
         ++(id_ptr_);
     return *(this);
     }
    // bool operator==(const NodeIterator&) const
    /** check whether two pointers points to the same node .
   * 
   * Complexity: O(1) amortized operations.
   */
    bool operator ==(const NodeIterator& node_iter1) const {
        return (this-> id_ptr_ == node_iter1.id_ptr_ && this->graph_ptr_ == node_iter1.graph_ptr_ );
    }
   private:
    friend class Graph;
    size_type id_ptr_;
    graph_type* graph_ptr_;
    NodeIterator(size_type id_ptr ,const graph_type* graph_ptr)
        :id_ptr_(id_ptr),graph_ptr_(const_cast<graph_type*>(graph_ptr)){}
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
     node_iterator nend(size(),this);
     return nend;
   }
 /**To remove node @a n with @a n.index() 
* @param[in] n_it pointer to the node
* @param[out] n_it pointer to the next node
*
* @post new @a graph has one less node, all nodes @a index() > @a n_it->index() are invalid aftre.
* nodes are shifted down by one. 
*       all edges connected to @a *(n_it) are removed. edge_vector (@a eid_) has @a n_it->degree() less elements.
*       edges with index larger than min(edge indecis connecting to @a *(n_it) are all invalid.  
* complexity O(num nodes connecting to this node+ number of edges connecting to adjacent nodes)
*/
  node_iterator remove_node(node_iterator n_it){
        auto n = *(n_it);
        ++(n_it);
        remove_node(n);
  return n_it;
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
    IncidentIterator():edge_id_ptr_(0),graph_ptr_e_(),node_ptr_() {
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
	
          size_type id1 = graph_ptr_e_->node_id_[node_ptr_->index()];
          size_type id2 = graph_ptr_e_->connectivity[id1][edge_id_ptr_].first;
          size_type eid = graph_ptr_e_->connectivity[id1][edge_id_ptr_].second.edge_id; 
         
        edge_type this_edge(eid,graph_ptr_e_);
         this_edge.node1_id_ = graph_ptr_e_->node_id_[node_ptr_->index()];
	this_edge.node2_id_ = id2;
         
         return this_edge;
     }
    // IncidentIterator& operator++()
     /** Return the iterator the next edge connected to the node.
   * 
   * Complexity: O(num of edges connected to node1) amortized operations.
   */
     IncidentIterator& operator++(){
          ++(edge_id_ptr_);
     return *(this);
     }
    // bool operator==(const IncidentIterator&) const

     /** Check whether  the pointer points to the same edge
   * Complexity: O(1) amortized operations.
   */
     bool operator==(const IncidentIterator& edge_iter) const{
         return (this->edge_id_ptr_ ==  edge_iter.edge_id_ptr_ &&
                 this->graph_ptr_e_ == edge_iter.graph_ptr_e_ && 
                  this->node_ptr_ == edge_iter.node_ptr_);

    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
   
    size_type edge_id_ptr_;
    graph_type* graph_ptr_e_;
    node_type* node_ptr_;
     IncidentIterator(size_type ptr_ ,const graph_type* graph_ptr_e,const node_type* node_ptr)
        :edge_id_ptr_(ptr_),graph_ptr_e_(const_cast<graph_type*>(graph_ptr_e)) ,node_ptr_(const_cast<node_type*>(node_ptr)){}
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
    EdgeIterator():edge_itr_ptr_(0),graph_ptr_e_itr_() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
      /** Return the edge pointed by the edgeiterator.
   * 
   * Complexity: O(1) amortized operations.
   */
      Edge operator*() const{
          edge_type new_edge = graph_ptr_e_itr_->edge(edge_itr_ptr_);
         return new_edge;
      }
    // EdgeIterator& operator++()
     /** Return the iterator the next edge.
   * 
   * Complexity: O(num of edges connected to node1) amortized operations.
   */

     EdgeIterator& operator++(){
        
        ++(edge_itr_ptr_);
         return *(this);
     }
    // bool operator==(const EdgeIterator&) const
     /** Check whether  the pointer points to the same edge
   * Complexity: O(1) amortized operations.
   */
      bool operator==(const EdgeIterator& edge_iter) const{
         return  (this->edge_itr_ptr_ == edge_iter.edge_itr_ptr_
			&& this->graph_ptr_e_itr_ == edge_iter.graph_ptr_e_itr_);
     }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
     size_type edge_itr_ptr_;
     graph_type* graph_ptr_e_itr_;
     EdgeIterator(size_type ptr_ ,const graph_type* graph_ptr_e)
        :edge_itr_ptr_(ptr_),graph_ptr_e_itr_(const_cast<graph_type*>(graph_ptr_e)){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
   /** Return the iterator to thefirst edge.
   * @post *(result.ptr_id_e) = first edge
   * Complexity: O(1) amortized operations.
   */
    edge_iterator edge_begin() const{
        //edge_iterator first_edge(0,this);
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
/**To remove edge *(e_it) 
* @param[in] e_it pointer to the edge
* @param[out] e_it iterator to the next edge
*
* @post new @a graph has one less edge connecting , 
*       all edges with @a edge_id > @a e_it->edge_id_.
*       nodes are all valid.
* element in edge_vector(@a eid_) are shifted down by one. 
* complexity O(num nodes connecting to two nodes of @a *(e_it))
*/

edge_iterator remove_edge(edge_iterator e_it){
	auto e = *(e_it);
        ++e_it;
	remove_iterator(e);
        return e_it;




}
 private:
  size_type size_,edge_size;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
