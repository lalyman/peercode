#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>


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

  /** Synonym for Value */
  using node_value_type = V;
  using edge_value_type = E;

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
  class Node : private totally_ordered<Node>{
   public:
    
    Node() {
    }

    const Point& position() const {
    	if(nid_ < graph_->nodes_.size())
        return graph_->nodes_[nid_].position;

		  assert(false);
    }

    Point& position() {
      if(nid_ < graph_->nodes_.size())
        return graph_->nodes_[nid_].position;

      assert(false);
    }

    size_type index() const {
    	if(nid_ < graph_->nodes_.size() )
	    	return nid_; 

		  assert(false);
    }

    /** 
    *  Get the user-specified value of each node
	  */
    node_value_type& value() {
    	if( nid_ < graph_->nodes_.size() )
        return graph_->nodes_[nid_].nodevalue; 

		  assert(false);
    }

    const node_value_type& value() const{
    	if( nid_ < graph_->nodes_.size() )
        return graph_->nodes_[nid_].nodevalue; 

      assert(false);
    }

    size_type degree() const{
    	return graph_->incidence_map.at(nid_).size();
    }

    incident_iterator edge_begin() const{
    	return IncidentIterator(graph_, 0, nid_);
    }

    incident_iterator edge_end() const{
    	return IncidentIterator(graph_, degree(), nid_);
    }

    bool operator==(const Node& n) const {
    	return (graph_ == n.graph_ and nid_ == n.index());
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
    	return ((nid_ < n.nid_ and graph_ == n.graph_) || graph_<n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type nid_;

    // private constructor
    Node(const Graph* graph, size_type nid)
	        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
          }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& pos, 
                const node_value_type& val = node_value_type()) {
  	proxy_node n {pos, val};
    nodes_.push_back(n);
  	//add a key in connect mapping and link it to an empty vector
  	std::vector<size_type> empty;
  	incidence_map[nodes_.size()-1] = empty;
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // if the node is not in this graph
    return (n.graph_== this and n.nid_<size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < size())
      return Node(this, i);

    assert(false);
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (eid_ < graph_->edges_.size())
        return graph_->node(aidx_);

      assert(false);
    }

    /** Return the other node of this Edge */
    Node node2() const {
    	if (eid_ < graph_->edges_.size())
        return graph_->node(bidx_);

      assert(false);
    }

    edge_value_type& value() {
      if (eid_ < graph_->edges_.size())
        return graph_->edges_[eid_].edgevalue;

      assert(false);
    }

    const edge_value_type& value() const {
      if (eid_ < graph_->edges_.size())
        return graph_->edges_[eid_].edgevalue;

      assert(false);
    }

    // Returns the euclidean length of the edge
    double length() const {
      return norm_2(node1().position()-node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.node1()==node1() && e.node2()==node2()) or 
              (e.node2()==node1() && e.node1()==node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	return ((eid_ < e.eid_ and graph_ == e.graph_ ) || graph_<e.graph_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions. 
    friend class Graph;
    Graph* graph_;
    size_type eid_;
    size_type aidx_;
    size_type bidx_;

    // private constructor
    Edge(const Graph* graph, size_type eid, size_type aidx, size_type bidx)
	        : graph_(const_cast<Graph*>(graph)), eid_(eid), 
	        aidx_(aidx), bidx_(bidx) {
	    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i < edges_.size()) {
      return Edge(this, i, edges_[i].n1_id, edges_[i].n2_id);
    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (this->has_node(a) and this->has_node(b)) {
      for (auto e : incidence_map.at(a.index())) {
        if (edges_[e].n1_id==b.index() or edges_[e].n2_id==b.index())
        return true;
      }
      // auto innermap = incidence_map.at(a.index()); // all adjacent edge-node map
      // for (auto it = innermap.begin(); it!=innermap.end(); ++it) {
      //   return it->second == b.index();
      // }
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
  	if (a == b) 
  		return Edge();
  	
  	if(has_edge(a, b)){
		  for (auto e : incidence_map.at(a.index())) {
        if (edges_[e].n1_id==a.index() and edges_[e].n2_id==b.index()) 
    			return edge(e);
    		if (edges_[e].n1_id==b.index() and edges_[e].n2_id==a.index())
    			return Edge(this, e, edges_[e].n2_id, edges_[e].n1_id);
		  }
      // auto innermap = incidence_map.at(a.index());
      // for (auto it = innermap.begin(); it!=innermap.end(); ++it) {
      //   if (it->second == b.index())
      //     return Edge(this, it->first, a.index(), b.index());
      // }
  	}

  	proxy_edge newedge {a.index(), b.index(), val};
  	edges_.push_back(newedge);

    incidence_map[a.index()].push_back(edges_.size()-1);
    incidence_map[b.index()].push_back(edges_.size()-1);
    // std::map<size_type,size_type> a_temp, b_temp;
    // a_temp[edges_.size()-1] = b.index();
    // b_temp[edges_.size()-1] = a.index();
  	// incidence_map.insert(
   //    std::pair<size_type,std::map<size_type,size_type>>(a.index(),a_temp)
   //    );
  	// incidence_map.insert(
   //    std::pair<size_type,std::map<size_type,size_type>>(b.index(),b_temp)
   //    );

  	return edge(edges_.size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
  	nodes_.clear();
  	edges_.clear();
  	incidence_map.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    /** 
  	* Dereference the node iterator
  	* @return The node with current node_id
    */
    Node operator*() const{
    	if (iid_ < graph_->size())
	    	return graph_->node(iid_);

	    assert(false);
    }

    /** 
	  * Increment the current node_id and returns the node iterator
    */
    NodeIterator& operator++() {
    	iid_++;
    	return *this;
    }

    /** 
  	* Check if two edge iterators are equal by comparing their edge_id 
  	* and if they are in the same graph.
    */
    bool operator==(const NodeIterator& i) const {
    	return (iid_ == i.iid_ and graph_ == i.graph_);
    }

    private:
    friend class Graph;
    Graph* graph_;
    size_type iid_;
    // private constructor
    NodeIterator(const Graph* graph, size_type iid)
	        : graph_(const_cast<Graph*>(graph)), iid_(iid) {
          }
  };

  /** 
  * Returns a node iterator with node id 0, the begin of a node iterator
  */
  node_iterator node_begin() const{
  		return NodeIterator(this, 0);
  }
  /** 
  * Returns the end of node iterator
  */
  node_iterator node_end() const{
  		return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. 
   */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** 
  	* Dereference the incident iterator
  	* @return The edge with one end connected to node with node_id nidx
    */
    Edge operator*() const {
    	if (incid_ < graph_->node(nidx_).degree()) {
        size_type temp_eid = graph_->incidence_map.at(nidx_)[incid_];
      	return graph_->edge(temp_eid);
    	}

    	assert(false);
    }

    /**
	* Increment the current edge_id 
	* @return The incident iterator representing the incident edge
    */
    IncidentIterator& operator++() {
    	incid_++;
    	return *this;
    }

    /** 
	* Check if two incident iterators are equal by comparing their host node id, 
	* positions in the mapping vector and if they are in the same graph
    */
    bool operator==(const IncidentIterator& i) const {
    	return (incid_ == i.incid_ and graph_ == i.graph_ 
    		and nidx_ == i.nidx_); 
    }

   private:
    friend class Graph; 
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type incid_; // index/position in the vector of edge_id 
    size_type nidx_;

    // private constructor
    IncidentIterator(const Graph* graph, size_type incid, size_type nidx)
		: graph_(const_cast<Graph*>(graph)), incid_(incid), nidx_(nidx) {
	}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()  {
    }

	/** 
	* Dereference the edge iterator
	* @return The edge with edge_id eiid_
  */
    Edge operator*() const {
    	if (eiid_<graph_->num_edges()){
	    	return graph_->edge(eiid_);
      }

	    assert(false);
    }

	/** 
	* Increment the current edge_id 
	* @return The edge iterator 
    */
    EdgeIterator& operator++() {
    	eiid_++;
    	return *this;
    }

    /** 
  	* Check if two edge iterators are equal by comparing their node_id, 
  	* and if they are in the same graph
    */
    bool operator==(const EdgeIterator& e) const {
    	return (graph_ == e.graph_ and eiid_ == e.eiid_);
    }

  private:
    friend class Graph;
    Graph* graph_;
    size_type eiid_; 

    // private constructor
    EdgeIterator(const Graph* graph, size_type eiid)
     : graph_(const_cast<Graph*>(graph)), eiid_(eiid) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** 
  * Returns an edge iterator with edge id 0, the begin of the edge iterator
  */
  EdgeIterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }
  /** 
  * Returns the end of edge iterator
  */
  EdgeIterator edge_end() const {
  	return EdgeIterator(this,num_edges());
  }

  //--------------Removal of edges--------------


  size_type remove_edge (const Node & a, const Node & b) {
    if(!has_edge(a,b)) return num_edges();

    size_type rm_idx, e;
    // loop through the vector of incidence to find the key for the edge
    for (size_type i=0; i<a.degree();i++){
      e = incidence_map.at(a.index())[i]; //edge index
      if (edges_[e].n1_id==b.index() or edges_[e].n2_id==b.index()) {
        rm_idx = e; //find the edge index
        incidence_map.at(a.index())[i]=incidence_map.at(a.index()).back();
        incidence_map.at(a.index()).pop_back();
        break;
      }
    } 

    for (size_type i=0; i<b.degree();i++){
      e = incidence_map.at(b.index())[i];
      if (edges_[e].n1_id==a.index() or edges_[e].n2_id==a.index()) {
        incidence_map.at(b.index())[i]=incidence_map.at(b.index()).back();
        incidence_map.at(b.index()).pop_back();
        break;
      }
    } 

    // update edge id in the proxy vector
    edges_[rm_idx]=edges_.back();
    edges_.pop_back();

    if (rm_idx!=num_edges()){
      //update the new edge_id rm_idx's end nodes
      size_type n1 = edges_[rm_idx].n1_id;
      size_type n2 = edges_[rm_idx].n2_id;

      // update node a index
      for (auto i=incidence_map[n1].begin();i!=incidence_map[n1].end();++i){
        if ((*i)==num_edges()) {
          (*i) = rm_idx;
          break;
        }
      }

      // update node b index
      for (auto i=incidence_map[n2].begin();i!=incidence_map[n2].end();++i){
        if ((*i)==num_edges()) {
          (*i) = rm_idx;
        }
      } 
    }
    return rm_idx;
  }

  size_type remove_edge (const Edge& e){
    return remove_edge(e.node1(),e.node2());
  }

  edge_iterator remove_edge (edge_iterator e_it){
    size_type rm_idx = remove_edge(*e_it);
    if (e_it==this->edge_end())
      return edge_end();
    else return EdgeIterator(this,rm_idx);
  }


  size_type remove_node (const Node & n){
    assert(has_node(n));
    size_type n_idx = n.index();
    size_type e;
    // std::cout<<"num edge = "<<num_edges()<<"\n";
    // remove all incident edge
    while(n.degree()>0){
      //current edge id
      e = incidence_map.at(n_idx)[0]; 
      
      remove_edge(edge(e));
    }
    
    // update node id in the proxy vector
    nodes_[n_idx]=nodes_.back();
    nodes_.pop_back();

    // update edges with node with original id = nodes_.size()-1;
    // it should now have edge id n_idx in the inc_map.
    size_type old_id = incidence_map.size()-1;
    incidence_map[n_idx]=incidence_map[old_id];
    for (auto e : incidence_map.at(n_idx)) {
      if (edges_[e].n1_id==old_id) edges_[e].n1_id=n_idx;
      if (edges_[e].n2_id==old_id) edges_[e].n2_id=n_idx;
    }
    incidence_map.erase(old_id);
    
    return n_idx;
  }

  node_iterator remove_node (node_iterator n_it) {
    size_type rm_idx = remove_node(*n_it);
    if (n_it==this->node_end())
      return node_end();
    else return NodeIterator(this,rm_idx);
  }


 private:

  struct proxy_node{
    Point position; //std::vector<Point> node_;
    node_value_type nodevalue; //std::vector<node_value_type> nodevalue_;
  };

 	struct proxy_edge{
    size_type n1_id;
    size_type n2_id; //std::vector<std::vector<size_type>> edge_;
    
    edge_value_type edgevalue; //std::vector<edge_value_type> edgevalue_;
  };

 	std::vector<proxy_node> nodes_;
  std::vector<proxy_edge> edges_;
  
  // Mapping from one node (node_id) to the set of 
  // edges (edge_id) that has one end on this node.
  std::map<size_type, std::vector<size_type>> incidence_map;
};

#endif // CME212_GRAPH_HPP

