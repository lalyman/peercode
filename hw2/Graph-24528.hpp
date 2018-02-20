#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


template <typename V, typename E>
class Graph {

  // PRIVATE TYPE DEFINITIONS
  private:
  /** declare the internal struct of Node. **/
    struct Internal_Nodes;
    struct Internal_Edges;

  // PUBLIC TYPE DEFINITIONS
  public:

    /** Node value type **/
    using node_value_type = V;

    /** Edge value type **/
    using edge_value_type = E;

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

    /** Type of indexes and sizes.
    Return type of Graph::Node::index(), Graph::num_nodes(),
    Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;


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


    // CONSTRUCTORS AND DESTRUCTOR

    /** Construct an empty graph. */
    /* The constructor initialises two empty vectors which store the information about nodes and edges of the graph.
       It also initializes the count of number of nodes and edges in the graph to be zero. */
    Graph() :
      nodes_(),i2u_(),e2u_(), num_nodes_(0),num_edges_(0),counter_euid(0),incident_edges_() { } 
      //i2u_ keeps track of accessible nodes
      //e2u_ keeps track of accessible edges
      //counter to keep track of 
    /** Default destructor */
      ~Graph() = default;

    // NODES

    /** @class Graph::Node
    * @brief Class representing the graph's nodes.
    * Node objects are used to access information about the Graph's nodes.*/


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
        Node() { } //Invalid Node constructor when called by the user directly

        /** Hw 2 implementation of function position() which returns a
        * a non-const reference to node position**/ 
        Point& position(){
          assert(valid());
          return const_cast<Graph*> (graph_) -> nodes_[uid_].position;
        }

        /** Return this node's position. */
        const Point& position() const {
          assert(valid());
          return graph_ -> nodes_[uid_].position; 
          //Access the position in the nodes vector of the graph containing the node at uid_
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
          assert(valid());
          return graph_ -> nodes_[uid_].idx; 
        } 

        /** Test whether this node and @a n are equal.
        * Equal nodes have the same graph and the same index.
        */
        bool operator==(const Node& n) const {
          assert(valid());
          assert(n.index()<n.graph_ -> num_nodes());
          return (n.uid_ == uid_ && n.graph_ == graph_); //Checking if the index and graph pointers of two nodes is the same
        }

        /** Test whether this node is less than @a n in a global order.
        * This ordering function is useful for STL containers such as
        * std::map<>. It need not have any geometric meaning.
        * The node ordering relation must obey trichotomy: For any two nodes x
        * and y, exactly one of x == y, x < y, and y < x is true.
        */

        bool operator<(const Node& n) const {
          assert(valid());
          assert(n.index() < n.graph_ -> num_nodes());
          return (uid_ < n.uid_); //  graph_ == n.graph_;
          // I guess it's not important to check if graphs are same
        }

        /** This function returns the value associated with the Node. **/

                // HW1: YOUR CODE HERE
                // Supply definitions AND SPECIFICATIONS for:
                // node_value_type& value();
                // const node_value_type& value() const;
                // size_type degree() const;
                // incident_iterator edge_begin() const;
                // incident_iterator edge_end() const;

        /** This method returns the value stored in node by reference.
        * The value can be altered by another methods.**/
        node_value_type& value() {
          assert(valid());
          return const_cast<Graph*>(graph_) -> nodes_[uid_].value ;
        }

        /** This method also returns the value stored in node but as a const.**/
        const node_value_type& value() const {
          assert(valid());
          return graph_ -> nodes_[uid_].value ;
        }

       /** This function calculates the degree of a node. **/
        size_type degree() const {
          assert(valid());
          return graph_ -> incident_edges_[uid_].size();
        } 

        /** this function constructs incident iterator for the node
        * object by calling the constructor.**/
        incident_iterator edge_begin() const {
         // return incident_iterator(graph_ -> nodes_[uid_].idx,0,graph_);
          return incident_iterator(uid_,0,graph_);
        }
        /** this function also constructs the incident iterator for this node
        * object  and index of edge passed is equal to the size of incident edges
        * vector for this node. **/ 
        incident_iterator edge_end() const {
         // return incident_iterator(graph_ -> nodes_[uid_].idx,graph_->incident_edges_[uid_].size(), graph_);
          return incident_iterator(uid_,graph_->incident_edges_[uid_].size(),graph_);
        } 

        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node o
      private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        size_type uid_; // Works like a uid for a node
        const Graph* graph_; // Pointer to a Graph to which the node belongs
        
        bool valid() const {
          return uid_ >=0 && uid_ < graph_ -> nodes_.size() 
                 && graph_ -> nodes_[uid_].idx < graph_ -> i2u_.size()
                 && graph_ -> i2u_[graph_->nodes_[uid_].idx] == uid_;
        }

        // Private Node constructor which takes in index and pointer to a graph
        Node(size_type uid, const Graph* graph) : uid_(uid),graph_(graph) {}
    };


    /** Return the number of nodes in the graph.
    *
    * Complexity: O(1).
    */
    size_type size() const {
      return num_nodes_; //Private data member of  graph class
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
      Internal_Nodes newNode;
      newNode.position = position;
      newNode.value = value;
      newNode.idx = num_nodes_;
      nodes_.push_back(newNode); //nodes_ vector is a vector of struct Internal Nodes 
      num_nodes_++; //The number of nodes is incremented by one
      size_type uid = nodes_.size()-1;
      i2u_.push_back(uid);//Pushing the uid of nodes in i2u_
      return Node(uid,this);// The new Node is returned by calling the Node constructor
    }

    /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
    */
    bool has_node(const Node& n) const {
      return (n.graph_ == this && 0 <= n.index() < num_nodes()) ;
     // return (n.graph_ == this && n.index() < i2u_.size() ); 
    }

    /** Return the node with index @a i.
        * @pre 0 <= @a i < num_nodes()
        * @post result_node.index() == i
        *
        * Complexity: O(1).
    */
    Node node(size_type i) const {
      assert( i < num_nodes_); 
      return Node(i2u_[i], this); //Return the node by calling the private constructor  of Node
    }

/** This method invalidates node passed to it and all incident edges.
  * @pre 0<= @a.index() < num_nodes()
  * @post result = Number of nodes removed
  * @post all incident edges to Node @a a are invalidated.
  * @post the indexes of all the valid nodes in graph change
  * @post the indexes of all the valid edges in graph change
  * @post new num_nodes() <= old num_nodes()
  * @post new num_edges() <= old num_edges()
  * Complexity: O(num_nodes())
  */

    size_type remove_node(const Node& a){
      assert(a.index() < num_nodes());
      size_type EdgesRemoved = 0;
      size_type incident_edges = incident_edges_[a.uid_].size();
      for(size_type i =0; i< incident_edges; i++){
        size_type incident_uid = incident_edges_[a.uid_][i].uidNode;
        if(std::find(i2u_.begin(),i2u_.end(),incident_uid)!=i2u_.end())
          EdgesRemoved+=remove_edge(a,Node(incident_uid,this));
      }
      auto it = i2u_.begin() + a.index();
      i2u_.erase(it);
      nodes_[a.uid_].idx = -1;//Invalidate node index
      num_nodes_= num_nodes_ -1;
      for(size_type i =0; i< i2u_.size(); i++){
        nodes_[i2u_[i]].idx = i;
      } 
      return 1; //the number of nodes removed
    }


    // EDGES

    /** @class Graph::Edge
    * @brief Class representing the graph's edges.
    * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
    * are considered equal if they connect the same nodes, in either order.
    */
    class Edge : private totally_ordered<Edge> {
      public:
        
        /** Construct an invalid Edge. */
        Edge() {    } // An invalid edge constructor
    
        /** Return a node of this Edge */
        Node node1() const {
            assert( 0 <= index1_ < graph_ -> num_nodes());
            return Node(graph_ -> i2u_[index1_],graph_);//get the uid for the node
          //Node constructor called by taking the uid and graph pointer of node 
        }

        /** Return the other node of this Edge */
        Node node2() const {
            assert( 0 <= index2_ < graph_ -> num_nodes());
            return Node(graph_ -> i2u_[index2_],graph_);//get the uid for the node 
        }

        /** Test whether this edge and @a e are equal.
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
        /** Check if the two nodes of two edges are same by 
        calling the node1() and node2() methods of edge objects. **/
          if((e.node1() == node1() && e.node2() == node2()) || 
             (e.node1() ==  node2() && e.node2() == node1()))
            return true;
          else
            return false;
        }

        /** Test whether this edge is less than @a e in a global order.
        *
        * This ordering function is useful for STL containers such as
        * std::map<>. It need not have any interpretive meaning.
        */
        bool operator<(const Edge& e) const {
        /** Compare the edge_index of two edges which is a 
        private attribute of edge object. **/
          return (graph_ == e.graph_ && uid_edge_ < e.uid_edge_);// && uid_edge_ < e.uid_edge_);
        }

        // Method to access edge value
        edge_value_type& value(){
          assert( 0 <= index1_ < graph_ -> num_nodes());
          assert( 0 <= index2_ < graph_ -> num_nodes()); 
          //Explore the incident edges of any of the two nodes of the edge
          for(size_type i = 0; i < graph_ -> node(index1_).degree();i++){
            if(graph_ -> incident_edges_[graph_ -> i2u_[index1_]][i].uidNode 
                == graph_->i2u_[index2_]){
              return const_cast<Graph*>(graph_) -> incident_edges_[graph_->i2u_[index1_]][i].
                      value;
            }
          } 
        //  return edge_value_type();
        }

        const edge_value_type& value() const {
          return value();
        }

        double length() const{
          return norm(graph_ -> node(index1_).position() - graph_ -> node(index2_).position());
        }

    private:
    // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      size_type index1_;// index associated with  node 1 of this edge object
      const Graph* graph_; // graph pointer associated with both nodes of this edge object
      size_type index2_; //index associated with node 2 of this edge object
      size_type uid_edge_; //index associated with this edge object

      /** Private edge constructor. Takes in index and graph pointers 
      * of two nodes and an edge index . **/

      Edge(size_type index1,const Graph* graph,size_type index2,size_type uid_edge): 
          index1_(index1),graph_(graph),index2_(index2), uid_edge_(uid_edge) { }
      // Use this space to declare private data members and methods for Edge
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Edge objects
    };

    /** Return the total number of edges in the graph.
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
    */
    size_type num_edges() const {
      return num_edges_;// Private data member of Graph class
    }

    /** Return the edge with index @a i.
        * @pre 0 <= @a i < num_edges()
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
    */
    Edge edge(size_type k) const {
      assert( k < num_edges());
      size_type euid = e2u_[k];//get me the edge uid from the index 
      for( size_type i = 0; i < incident_edges_.size(); i++){
        for(size_type j = 0; j < incident_edges_[i].size(); j++){
          if(incident_edges_[i][j].uidEdge ==  euid){//compare the uid of edge
            //find the index of other node
            size_type indexNode2=std::find(i2u_.begin(),i2u_.end(),incident_edges_[i][j].uidNode) -
                                i2u_.begin();
            size_type indexNode1= std::find(i2u_.begin(),i2u_.end(),i)-i2u_.begin();
            if(indexNode1 < indexNode2)
              return Edge(indexNode1,this,indexNode2,euid); 
            else
              return Edge(indexNode2,this,indexNode1,euid);
              //Construct edge objects  using edge uids
          }
        }
      }
    return Edge();
    }

    /** Test whether two nodes are connected by an edge.
        * @pre @a a and @a b are valid nodes of this graph
        * @return True if for some @a i, edge(@a i) connects @a a and @a b.
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
    */
    bool has_edge(const Node& a, const Node& b) const {
    /** call the has_node() method of Graph class to check 
    if the nodes passed are valid nodes of the graph. **/                
      assert ( has_node (a)  &&  has_node (b) );
      size_type uid_edge;
      bool edgeExists=false;
      if(incident_edges_.size()==0) return false;
      for ( size_type i = 0; i < incident_edges_[i2u_[a.index()]].size() ; i++){ 
      //iterate over incident_edges of node a of graph
        if( incident_edges_[i2u_[a.index()]][i].uidNode == i2u_[b.index()] ){
          uid_edge = incident_edges_[i2u_[a.index()]][i].uidEdge;
          edgeExists = true;
        }
      }
      if(edgeExists)
        return(std::find(e2u_.begin(),e2u_.end(),uid_edge)!=e2u_.end());//edge exists and valid
      else
       return false;//edge does not exist
    }

    /** Add an edge to the graph, or return the current edge if it already exists.
        * @pre @a a and @a b are distinct valid nodes of this graph
        * @return an Edge object e with  e.node1() == @a a and e.node2() == @a b
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
      assert ( has_node(a) && has_node(b) && !(a==b)); // valid nodes of graph and distinct
      /** New implementation of add_edges. **/
      size_type n1, n2;
      if ( a.index()< b.index()) { n1 = a.index(); n2 = b.index();}
      else { n1 = b.index(); n2 = a.index();}
      incident_edges_.resize(num_nodes());
      /** Iterating over the incident edges of node a and 
          checking if the other end has node b.**/
      size_type uid_edge;
      for(size_type i= 0; i< incident_edges_[i2u_[a.index()]].size();i++){
        if(incident_edges_[i2u_[a.index()]][i].uidNode == i2u_[b.index()]){
          uid_edge = incident_edges_[i2u_[a.index()]][i].uidEdge;
          if(std::find(e2u_.begin(),e2u_.end(),uid_edge)!=e2u_.end()){
              //Edge exists and is valid
              return Edge(a.index(),this,b.index(),uid_edge);
          }
          else{//Edge exists and it is invalidated
            e2u_.push_back(uid_edge);
            num_edges_ = num_edges_ +1;
            return Edge(a.index(),this,b.index(),uid_edge);
          }
        }
      }
     //Edge does not exist
      num_edges_++; // if no edge exists then create a new edge and increment the counter
      counter_euid ++;
      Internal_Edges newEdge;
      newEdge.uidNode=i2u_[b.index()];
      newEdge.uidEdge=counter_euid-1;//num_edges_-1;
      newEdge.value=edge_value_type(); 
      incident_edges_[i2u_[a.index()]].push_back(newEdge);
      newEdge.uidNode = i2u_[a.index()];
      incident_edges_[i2u_[b.index()]].push_back(newEdge);
      e2u_.push_back(counter_euid -1); 
      return Edge(a.index(),this,b.index(), counter_euid -1); //return the last added edge
    }

/** This method is used to invalidate edge between two nodes.
  * @pre 0<=@a a.index() < num_nodes()
  * @pre 0<=@a b.index() < num_nodes()
  * @pre @a a.valid() and @a b.valid()
  * @pre edge between Node @a a and @a b is edge e
  * @post the result gives the number of edges invalidated by the method
  * @post The edge between @a Node a and @a b is invalidated
  * @post @a e2u_ vector .find(e) == e2u_.end()
  * @post the indexes of all the valid edges in graph change
  * @post new num_edges() <= old num_edges()
  * Complexity: O(num_nodes() + num_edges()) assuming max degree() of Node @a a and @a b is small
  */
    size_type remove_edge(const Node& a, const Node& b){
      assert(has_node(a) && has_node(b));
      if(!has_edge(a,b)){return 0;}//if no edge between the two nodes
      size_type euid;
      for(size_type i =0; i < incident_edges_[i2u_[a.index()]].size(); i++){
        if(incident_edges_[i2u_[a.index()]][i].uidNode == i2u_[b.index()]){
          euid= incident_edges_[i2u_[a.index()]][i].uidEdge;
          auto it = std::find(e2u_.begin(),e2u_.end(),euid);
          if(it!=e2u_.end()){
            e2u_.erase(it);
            num_edges_ = num_edges_ -1 ;//decrease the number of edges by 1
          }
        }
      }
      return 1; //number of edges removed by the function
    }

    /** This method also invalidates edge passed to it.
    * @pre Edge @a a is valid i.e. find(e2u_.begin(),e2u_.end(),a.uid_edge_)!=e2u_.end()
    * @post e2u_.find(@a a) == e2u_.end()
    * @post returns the number of edges invalidated by the method
    * @post the indexes of all the valid edges in graph change
    * @post new num_edges() <= old num_edges()
    * Complexity: O(num_nodes() + num_edges())
    */
    size_type remove_edge(const Edge& a){
      assert(std::find(e2u_.begin(),e2u_.end(),a.uid_edge_)!=e2u_.end());
      return remove_edge(a.node1(),a.node2());
    }

    /** This method invalidates edge pointed by edge iterator @a e_it.
  * @pre @a *e_it is a valid edge i.e. find(e2u_.begin(),e2u_.end(),*e_it.uid_edge_)!=e2u_.end()
  * @post result gives the number of edges invalidated by method.
  * @post all edge iterators in the range [it, end()] are invalidated
  * @post the indexes of all the valid edges in graph change
  * @post new num_nodes() <= old num_nodes()
    @Complexity: O(num_nodes() + num_edges())
    */
    size_type remove_edge(edge_iterator e_it){
      return remove_edge(*e_it);
    }
    /** Remove all nodes and edges from this graph.
        * @post num_nodes()  == 0 && num_edges() == 0
        *
        * Invalidates all outstanding Node and Edge objects.
    */
    void clear() {
     // edges_.clear(); //clears all the vectors and assigns 0 to all the counters
      nodes_.clear();
      num_nodes_=0;
      num_edges_=0;
      incident_edges_.clear();
      i2u_.clear();
      e2u_.clear();
      counter_euid = 0;
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
      NodeIterator() {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** This function is used for derefencing node iterator pointer. **/
      Node operator*() const {
      
        if (indexIt_ < graphIt_ -> num_nodes())
          return Node(graphIt_ -> i2u_[indexIt_],graphIt_); //calling the node constructor
        return Node();// return an invalid node if index is greater than num_nodes
      }
      //indexIt_ is used to iterate over i2u_ vector
    /** This function is used to implement ++ operator for node iterator. **/
      node_iterator& operator++() {
        if(indexIt_ == graphIt_ -> num_nodes())
          return *this; //do not increment the index
        // Increment the index by 1 and call the nodeIterator constructor
        indexIt_ = indexIt_ + 1;
        return *this;
      } 

    /** This function checks if two iterators are equal. **/
      bool operator==(const node_iterator& secondIt) const {
        return ((graphIt_ == secondIt.graphIt_) && (indexIt_ == secondIt.indexIt_));
      }

      private:
        //Private Node Iterator constructor
        NodeIterator(const Graph* graph, size_type index) : graphIt_(graph), indexIt_(index) {}
        const Graph* graphIt_;//private attributes
        size_type indexIt_;
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
    };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
    node_iterator node_begin() const {
      return NodeIterator(this, 0); 
      //call the node iterator with this graph object and index 0 of node
    }

    node_iterator node_end() const {
      return NodeIterator(this, num_nodes());
    }//call the node iterator with this graph and index equal to num of nodes

    /** This method invalidates the Node object @a *n_it.
    *@pre  Valid iterator @a n_it i.e. points to a valid Node object
    *@post Invalidates all iterators in the range[it, end())
    *@post the result gives the number of nodes invalidated by the method
    *@post the edges incident to @a *n_it are invalidated
    * @post the indexes of all the valid nodes in graph change
    * @post the indexes of all the valid edges in graph change
    * @post new num_nodes() <= old num_nodes()
    * @post new num_edges() <= old num_edges()
    * Complexity: O(num_nodes())
    */

    node_iterator remove_node (node_iterator n_it){
      assert(*n_it.index() < num_nodes() && n_it != i2u_.end());
      remove_node(*n_it);
      return n_it;//dont know if i can return n_it
    }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator : private totally_ordered<IncidentIterator>{
      public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

      /** Construct an invalid IncidentIterator. */
      IncidentIterator() {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
  /** This method is used for derefencing incident iteartor.**/
      Edge operator*() const {
        // the idea is go to the uid_nodeIIt and index_edgeIIt and get the uid .
        //then get the index and call the edge method
        if(index_edgeIIt_ < graphIIt_ -> incident_edges_[uid_nodeIIt_].size()){
          unsigned edge_uid =  graphIIt_ -> incident_edges_[uid_nodeIIt_][index_edgeIIt_].uidEdge;
          auto it = std::find(graphIIt_ -> e2u_.begin(),graphIIt_ ->e2u_.end(),edge_uid);
          if( it != graphIIt_->e2u_.end())
            return graphIIt_ -> edge(size_type(it-graphIIt_->e2u_.begin()));
        }
        return Edge();
      }

  /** This method is used for incrementing the incident iterator. **/
      incident_iterator& operator++() {
        if(index_edgeIIt_ == graphIIt_ -> incident_edges_[uid_nodeIIt_].size())
          return *this; //not incrementing
        else
          index_edgeIIt_ = index_edgeIIt_ + 1;
        return *this;
      }
  /** This method compares two incident iterators.**/
      bool operator==(const incident_iterator& iit) const {
        return (graphIIt_ == iit.graphIIt_ && uid_nodeIIt_ == iit.uid_nodeIIt_
                && index_edgeIIt_ == iit.index_edgeIIt_);
      }
      private:
        friend class Graph;
        //Private incident iterator takes in node index,edge index and graph
        IncidentIterator(size_type uid_nodeIIt, size_type index_edgeIIt, const Graph* graph):
        uid_nodeIIt_(uid_nodeIIt), index_edgeIIt_(index_edgeIIt), graphIIt_(graph) {}
        size_type uid_nodeIIt_;
        size_type index_edgeIIt_;
        const Graph* graphIIt_;
    // HW1 #3: YOUR CODE HERE
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
  /** Method to deference edge iterator.**/
      Edge operator*()const {
        if(index_ < graph_ -> num_edges())
          return graph_ -> edge(index_);
        else
          return Edge(); //return an invalid edge object
      }
  /** Method to increment edge iterator.**/
      EdgeIterator& operator++() {
        if(index_ == graph_ -> num_edges()) //edges_.size())
          return *this; //not incrementing end iterator
        index_ = index_ + 1;
        return *this;
      }
  /** this method checks for equality of two edge iterators. **/
      bool operator==(const EdgeIterator& second) const{
        return (graph_ == second.graph_ && index_ == second.index_);
      }

      private:
        friend class Graph;
    // HW1 #5: YOUR CODE HERE
    //edge iterator has two private attributes edge index and graph pointer
        const Graph* graph_;
        size_type index_;
        EdgeIterator(const Graph* graph, size_type index): graph_(graph), index_(index) {} 
    };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    // edge iterator to the first element of edges_ vector
    edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
    }
    //edge iterator to the last element of edges_ vector
    edge_iterator edge_end() const{
      return EdgeIterator(this, num_edges_);
    }
    // PRIVATE TYPE DEFINITIONS
  private:
    // defining the internal structure of Node
    struct Internal_Nodes{
      size_type idx; //index associated with the node
      Point position;
      node_value_type value;
    };
    // defining the internal structure of Edge
    struct Internal_Edges{
      size_type uidNode; //this should be the uid of node index of incident node
      size_type uidEdge; //index of incident edge
      edge_value_type value;
    };
    // New definition of nodes_;
    std::vector<Internal_Nodes> nodes_;
    std::vector<size_type> i2u_;
    std::vector<size_type> e2u_;
    size_type num_nodes_; // counter for number of nodes in graph
    size_type num_edges_; // counter for number of edges in graph
    size_type counter_euid;
    std::vector<std::vector<Internal_Edges>> incident_edges_;
    // this vector stores the indexes of incident edges for all the nodes
};

#endif





