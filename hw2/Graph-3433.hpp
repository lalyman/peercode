#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <functional>

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
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  //
  // PUBLIC TYPE DEFINITIONS
  //
 public:
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
  
  //Type of the value assosciated with nodes and edges.
  using node_value_type=V;
  using edge_value_type=E;
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

  //This iterator class is used to iterate through the Nodes of a graph.
  class NodeIterator{
    public:
	 // These type definitions let us use STL's iterator_traits.
	using value_type        = Node;                     // Element type
	using pointer           = Node*;                    // Pointers to elements
	using reference         = Node&;                    // Reference to elements
	using difference_type   = std::ptrdiff_t;           // Signed difference
	using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
    //Construct an invalid NodeIterator, that does not belong to any graph.
    NodeIterator(){
        graph=nullptr;
        ptr=0;
    }
    
    //Derefernce the NodeIterator, will return a Node object.
    Node operator *() const{
        return graph->node(this->ptr);
    };
    
    //Increment the NodeIterator and return it.
    NodeIterator& operator ++(){
        ptr++;
        return *this;
    };

    //Check equality by checking if two iterators belong to the same graph and point to the same
    //position.
    bool operator ==( const NodeIterator& iter ) const{
        if(iter.graph!=this->graph){
             return false;
        }
        return (iter.ptr==this->ptr);
    };

    //The negation of the equality operator.
    bool operator !=(const node_iterator& iter ) const{
        return !(iter==*this);
    };

    //The private attributes of the iterator are a pointer to the graph it
    //belongs to. The position of the pointer, and a constructor that can be
    //called by the Graph class.
    private:
        friend class Graph;
        size_type ptr;
        Graph* graph;
        NodeIterator(const Graph* id_graph,size_type index){
            this->graph=const_cast<Graph*>(id_graph);
            this->ptr=index;
        };

  };

  //Create a NodeIterator that points to the beginning of the NodeInfo
  //vector.
  NodeIterator node_begin() const{
    return node_iterator(this,0);
  }

  //Create a NodeIterator that points to one index past the end of the NodeInfo
  //vector.
  NodeIterator node_end() const{
    return node_iterator(this,this->num_nodes());
  }
  
  //This class iterates through all the adjacent nodes to a root node.
  class IncidentIterator{
    public:
	using value_type        = Edge;                     // Element type
   	using pointer           = Edge*;                    // Pointers to elements
   	using reference         = Edge&;                    // Reference to elements
   	using difference_type   = std::ptrdiff_t;           // Signed difference
   	using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
    //Construct an invalid IncidentIterator, that does not belong to any graph.
    IncidentIterator(){
        ptr_=0;
        nid_=0;
        graph=nullptr;
    }

    //Derefernce the IncidentIterator, will return an Edge object.
	Edge operator *() const{
            return Edge(graph,nid_,graph->Points[nid_].adj_[ptr_]);
    	}

    //Increment the IncidentIterator and return it.
    IncidentIterator& operator ++(){
        ptr_++;
        return(*this);
    }
    
    //check if two iterators are equal by checking if they belong to the same
    //graph and have the same indices for both the adjacency list and the root
    //node.
    bool operator ==( const incident_iterator & iit ) const{
        if(iit.graph!=this->graph){
            return false;
        }
        return (iit.ptr_==this->ptr_&&iit.nid_==this->nid_);
    };
    
    //The negation of the equality operator.
    bool operator !=(const incident_iterator& iter ) const{
        return !(iter==*this);
    };

    //The private attributes of the iterator are a pointer to the graph it
    //belongs to. The index of the root node, the position in the adjaceny list, 
    //and a constructor that can be called by the Graph and Node classes.
    private:
        size_type ptr_;
        size_type nid_;
        friend class Graph;
        friend class Node;
        Graph* graph;
        IncidentIterator(const Graph* id_graph,size_type index,size_type nid){
            this->graph=const_cast<Graph*>(id_graph);
            this->ptr_=index;
            this->nid_=nid;
        };
  };
  
  //This class is used to iterate through all the edges belonging to a graph.
  class EdgeIterator{
    public:
    	using value_type        = Edge;                     // Element type
    	using pointer           = Edge*;                    // Pointers to elements
    	using reference         = Edge&;                    // Reference to elements
    	using difference_type   = std::ptrdiff_t;           // Signed difference
    	using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        //Construct an invalid edge operator that does not belong to a graph,
        EdgeIterator(){
            nid1_=0; 
            nid2_=0; 
            graph_=nullptr;
        }

        //Derefernce the iterator, this returns an Edge object.
    	Edge operator *() const{
            return(Edge(graph_,nid1_,graph_->Points[nid1_].adj_[nid2_]));
        }
    
        //Increment the operator. If we are not at the end of the adjaceny list
        //of the first node we increment the second position counter.
        //Otherwise we increment the first position counter and set the second
        //one to zero.
        //To avoid double counting edges we keep incrementing until
        //node1.index()>node2.index().
        EdgeIterator& operator ++(){
            nid2_++;
            ValidateEdge();
            return(*this);
        }

        //Two iterators are equal if they belong to the same graph and their
        //id's are identical.
        bool operator ==( const edge_iterator & eit ) const {
            if(eit.graph_!=this->graph_){
                return false;
            }
            return (eit.nid1_==this->nid1_ && eit.nid2_==this->nid2_);
        }
        
        //The negation of the equality operator.
        bool operator !=( const edge_iterator & eit ) const {
            return!(eit==*this);
        }
  
      //The private attributes of the iterator are a pointer to the graph it
      //belongs to. The index of the node1, the position in the adjaceny list, 
      //and a constructor that can be called by the Graph class.
      private:
        friend class Graph;
        size_type nid1_;
        size_type nid2_;
        Graph* graph_;
        EdgeIterator(const Graph* graph,size_type nid1,size_type nid2){
            this->nid1_=nid1;
            this->nid2_=nid2;
            this->graph_=const_cast<Graph*>(graph);
            ValidateEdge();
		}

        /* Check if the iterator points at a valid edge.
         * If we are at the end of the adjacency list of the first node we move
         * over to the next node.
         * If we are at the end of Nodes vector we exit as this is edge_end().
         * If index of node 2 is larger than node one we contnue. This prevents
         * the double counting of edges.
         * If none of the conditions above hold we exit the function.
         * This function is called in the constructor to see if the first edge
         * is valid and in the incement operator to see if the new edge is
         * valid.
         */
        void ValidateEdge()
        {
            while(true){
              if(nid1_==graph_->num_nodes()&&nid2_>=0){
                  break;
              }
              if(nid2_==(graph_->node(nid1_).degree())){
                nid1_++;
                nid2_=0;
                continue;
              }
              if(graph_->Points[nid1_].adj_[nid2_]<nid1_){
                nid2_++;
                continue;
              }
              break;
            }
        }
  };
  
  //This function returns an EdgeIterator to the first edge. This is defined as
  //the edge that connects node(0) and the first node in its adjaceny list.
  EdgeIterator edge_begin() const {
    return(edge_iterator(this,0,0));
  }

  //This function returns an EdgeIterator to the end of the edge list. This is
  //define as one past the node indices for the first node and zero for the
  //second node.
  EdgeIterator edge_end() const{
    return(edge_iterator(this,this->num_nodes(),0));
  }


  class Node: private totally_ordered<Node> {
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
        this->id=size_type(-1);
        this->graph=nullptr;
    }

    /** Return a const reference this node's position. */
    const Point& position() const {
      return graph->Points[id].point_;
    }

    /** Return a reference this node's position. */
    Point& position() {
      return graph->Points[id].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        return (n.graph==this->graph&&n.index()==this->id);
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
        if(this->id<n.index()){
            return true;
        }
        if(this->id==n.index()){
            std::less<Graph*> comp;
            return comp(this->graph,n.graph);
        }
      return false;
    }

    //Return the value assosciated with a node.
    //The first function allows the value to be modified, the second declares it
    //to be a const.
    node_value_type& value (){
        return graph->Points[id].node_value_;
    }
    const node_value_type & value () const{
        return graph->Points[id].value_;
    }

    //Return the number of incident edges of a node.
    size_type degree() const{
        return graph->Points[id].adj_.size();
    }

    //Return an Incident iterator to the start of the adjacency list. 
    IncidentIterator edge_begin() const{
        return IncidentIterator(graph,0,id); 
    }

    //Return an Incident iterator to one past the end of the adjacency list.
    IncidentIterator edge_end() const{
        return IncidentIterator(graph,this->degree(),id); 
    
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    size_type id;
    Graph* graph;
    Node(size_type i,const Graph* id_graph){
        this->graph=const_cast<Graph*>(id_graph);
        this->id=i;
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
    return Points.size();
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

  Node add_node ( const Point& position , const node_value_type & value = node_value_type()){
    NodeInfo new_node= NodeInfo(position,value);
    Points.push_back(new_node);
    return this->node(this->Points.size()-1);        
  }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    bool in_graph;
    if(n.index()>this->size()){
        in_graph=false;
    }
    in_graph=(this==n.graph);
    return in_graph;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  
  Node node(size_type i) const {
    if (i>=this->Points.size()){
        std::cout<< "Index: "<<i<<" out of bounds."<<std::endl; 
        return Node();
    }
    else{
        return Node(i,this);        
    }
  }
 /* Dereference the node iterator and call remove_node with that node. 
  */
  NodeIterator remove_node(NodeIterator& nit){
    Node n = *nit;
    remove_node(n);
    return nit;
  }


  /* Remove a node by going through the following steps.
   * Remove all the incident edges on the node by repeatedly calling
   * edge_remove on the adjaceny list until it is empty.
   * 
   * Get the index of the node that is removed and the last Node.
   * 
   * Copy the last node to the position of the node to be removed.
   * 
   * Go through the adjacency list of the last node and for each of the adjacent
   * nodes update their adjacency list such that they point to the new index
   * instead of the last node.
   * Copy the edge_values with the correct index into both the new node's
   * map and the second nodes map. 
   * 
   * Remove the node at the end.
   *
   *If the node is not present return 1 of the node is removed succesfully return 0.
   */
  size_type remove_node(const Node& n){
    if(!has_node(n)){
        return 0;
    }
    for(unsigned i = 0; i<n.degree();){
        this->remove_edge(*n.edge_begin());
    }
    size_type new_index = n.index();
    size_type old_index = this->num_nodes()-1;
    this->Points[new_index]=this->Points[old_index];
    Node old_node = this->node(new_index);
    for(unsigned i =0; i<old_node.degree();i++){
        Node n2 = node(Points[old_index].adj_[i]);
        for(unsigned j=0;j<n2.degree();j++){
            if(Points[n2.index()].adj_[j]==old_index){
                Points[n2.index()].adj_[j]=new_index;
            }
        }
            Points[n2.index()].edge_values_[new_index]=Points[n2.index()].edge_values_[old_index];
            Points[new_index].edge_values_[n2.index()]=Points[n2.index()].edge_values_[new_index];
    } 

    /*
    for(unsigned i = 0; i<old_node.degree();i++){
        this->remove_edge(*old_node.edge_begin());
    }*/
    Points.pop_back();
    return 1;
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
        nid1_=size_type(-1);
        nid2_=size_type(-1);
        graph=nullptr;

    }


    /** Return a node of this Edge */
    Node node1() const {;
        return graph->node(nid1_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return graph->node(nid2_);
    }
    
    /* Return the value assosciated with an edge.
     * These values are stored in a map stored in NodeInfo.
     * To ensure that the storage is consistent we always find the edge value by
     * finding it in the map assosciated with the smaller node index of the
     * edge's nodes.
     */
    const edge_value_type& value() const{
        if(nid1_>nid2_){
            return graph->Points[nid2_].edge_values_[nid1_];
        }
        return graph->Points[nid1_].edge_values_[nid2_];
    }
    edge_value_type& value(){
        if(nid1_>nid2_){
            return graph->Points[nid2_].edge_values_[nid1_];
        }
        return graph->Points[nid1_].edge_values_[nid2_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if(this->node1()==e.node1()&&this->node2()==e.node2()){
            return true;
        }       
        if(this->node1()==e.node2()&&this->node2()==e.node1()){
            return true;
        }       
        return false;
    }
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        //Order edges by first comparing node1 of both edges and only if they
        //are equal use node2 for comparison.
        if(this->node1()<e.node1()){
            return true;
        }
        if(this->node1()==e.node1()){
            return (this->node2()<node2());
        }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type nid1_;
    size_type nid2_;
    Graph* graph;
    Edge(const Graph* id_graph,size_type nid1,size_type nid2){
        this->graph=const_cast<Graph*>(id_graph);
        this->nid1_=nid1;
        this->nid2_=nid2;
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
	  if(i>=num_edges()){
		return Edge();
	  }
      edge_iterator eit=this->edge_begin();
      for(size_type k=0;k<i;k++){
        ++eit;
      }
      Edge e= *eit;
      return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   // Check if an edge is present by checking if one node of the edge is present
   // in the other node's adjaceny list.
   if(std::find(Points[a.index()].adj_.begin(),Points[a.index()].adj_.end(),b.index())!=Points[a.index()].adj_.end()){
        return true;
      }
    else{
        return false;
      }
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *       edge.value() has the value given by the default constructor of
   *       edge_value_type.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    //Check if an edge already is present in the graph.
    //If it isn't update the adjaceny lists of the nodes that the edge connects.
    if(has_edge(a,b)){
        return Edge(this,a.index(),b.index());
    }

    edge_num++;
    this->Points[a.index()].adj_.push_back(b.index());
    this->Points[b.index()].adj_.push_back(a.index());
    this->Points[a.index()].edge_values_[b.index()]= edge_value_type {};
    this->Points[b.index()].edge_values_[a.index()]= edge_value_type {};
    return Edge(this,a.index(),b.index());
  }
  
/* Call remove_edge() with the two nodes that define an edge.
 * 
 */
  size_type remove_edge(const Edge& e){
      return remove_edge(e.node1(),e.node2());
  }

/* Dereference the EdgeIterator call remove_edge() with its two nodes.
 * To return a valid iterator we set the iterator to the next valid edge if the
 * current node1() of the edge has no incident edges.
 * return the iterator.
 */
  EdgeIterator remove_edge(EdgeIterator& eit){
    Edge e =*eit;
    remove_edge(e.node1(),e.node2());
    if(e.node1().degree()==0){
        eit.nid1_++;
        eit.nid2_=0;
        eit.ValidateEdge();
    }
    return eit;
  }
/* Remove an edge connecting two nodes. The edge does not have to exist but the
 * nodes must be valid.
 * Find the index of the second node in the adjaceny list of the first.
 * If the index is in the list remove it by copying the back to the position and
 * poping the back.
 *
 * Do the same for the index of the first node in the adjacency list of the
 * second.
 * Decrement the number of edges counter and return 0.
 *
 * If the edge does not exist neither of the if statements will be true and
 * therefor we go straight to the return statement at the end and return 1.
 *
 */
  size_type remove_edge(const Node& n1, const Node& n2){
    size_type nid1 = n1.index();
    size_type nid2 = n2.index();
    auto rm_it1=std::find(this->Points[nid1].adj_.begin(),this->Points[nid1].adj_.end(),nid2);
    if(rm_it1!=this->Points[nid1].adj_.end()){
        *rm_it1=(*(--this->Points[nid1].adj_.end()));
        this->Points[nid1].adj_.pop_back();
        this->Points[nid1].edge_values_.erase(nid2);
    } 
    auto rm_it2=std::find(this->Points[nid2].adj_.begin(),this->Points[nid2].adj_.end(),nid1);
    if(rm_it2!=this->Points[nid2].adj_.end()){
        *rm_it2=(*(--this->Points[nid2].adj_.end()));
        this->Points[nid2].adj_.pop_back();
        this->Points[nid2].edge_values_.erase(nid1);
        edge_num--;
        return 1;
    } 
    return 0;
  }
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    edge_num=0;
  }

 private:
  // This class contains the actual information of a node. It stores the position,
  // value and adjaceny list for each node. To store the edge values each
  // NodeInfo has a map structure. The keys are the second nodes that the edge
  // connects the first to. This means that the edge values are stored twice.
  // the node connects to while the value is the node 
  // The Node class is effectively a 
  //pointer to an object of NodeInfo.
    class NodeInfo{
        private:
            friend class Graph;
            Point point_;
            std::vector<size_type> adj_; 
            node_value_type node_value_;
            std::map<size_type,edge_value_type> edge_values_;
            NodeInfo(Point point, node_value_type node_value){
                this->point_=point;
                this->node_value_=node_value;
            }
    };
    //The Graph class stores a vector of NodeInfo objects and a counter to keep
    //track of the number of edges in the graph.
    std::vector<NodeInfo> Points;
    size_type edge_num=0;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
