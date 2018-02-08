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

  // PRIVATE TYPE DEFINITIONS
  private:
  /** declare the internal struct of Node. **/
    struct Internal_Nodes;


  // PUBLIC TYPE DEFINITIONS
  public:

    /** Node value type **/
    using node_value_type = V;

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
      nodes_(),edges_(), num_nodes_(0),num_edges_(0), incident_edges_(){ }               

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

        /** Return this node's position. */
        const Point& position() const {
          return graph_ -> nodes_[index_].position; 
          //Access the position in the nodes vector of the graph containing the node at index_
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
          return index_; //This returns the index associated with Node object
        }

        /** Test whether this node and @a n are equal.
        * Equal nodes have the same graph and the same index.
        */
        bool operator==(const Node& n) const {
          return (n.index_ == index_  && n.graph_ == graph_); //Checking if the index and graph pointers of two nodes is the same
        }

        /** Test whether this node is less than @a n in a global order.
        * This ordering function is useful for STL containers such as
        * std::map<>. It need not have any geometric meaning.
        * The node ordering relation must obey trichotomy: For any two nodes x
        * and y, exactly one of x == y, x < y, and y < x is true.
        */

        bool operator<(const Node& n) const {
          return (index_ < n.index()); //  graph_ == n.graph_;
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
          return const_cast<Graph*>(graph_) -> nodes_[index_].value ;
        }

        /** This method also returns the value stored in node but as a const.**/
        const node_value_type& value() const {
          return graph_ -> nodes_[index_].value ;
        }

       /** This function calculates the degree of a node. **/
        size_type degree() const {
          return graph_ -> incident_edges_[index_].size();
        } 

        /** this function constructs incident iterator for the node
        * object by calling the constructor.**/
        incident_iterator edge_begin() const {
          return incident_iterator(index_,0,graph_);
        }
        /** this function also constructs the incident iterator for this node
        * object  and index of edge passed is equal to the size of incident edges
        * vector for this node. **/ 
        incident_iterator edge_end() const {
          return incident_iterator(index_,graph_->incident_edges_[index_].size(), graph_);
        } 

        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node o
      private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        size_type index_; // Works like a uid for a node
        const Graph* graph_; // Pointer to a Graph to which the node belongs

        // Private Node constructor which takes in index and pointer to a graph
        Node(size_type index, const Graph* graph) : index_(index),graph_(graph) {}
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
      nodes_.push_back(newNode); //nodes_ vector is a vector of struct Internal Nodes 
      num_nodes_++; //The number of nodes is incremented by one
      return Node(num_nodes_-1,this);// The new Node is returned by calling the Node constructor
    }

    /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
    */
    bool has_node(const Node& n) const {
      return (n.graph_ == this && n.index() < num_nodes()) ;
    }

    /** Return the node with index @a i.
        * @pre 0 <= @a i < num_nodes()
        * @post result_node.index() == i
        *
        * Complexity: O(1).
    */
    Node node(size_type i) const {
      assert(i < size());
      return Node(i, this); //Return the node by calling the private constructor  of Node
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
          return Node(index1_,graph1_);  
          //Node constructor called by taking the index and graph pointer of node 
        }

        /** Return the other node of this Edge */
        Node node2() const {
          return Node(index2_,graph2_);      // call the Node constructor
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
          return (index_edge_ < e.index_edge_);
        }

    private:
    // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      size_type index1_;// index associated with  node 1 of this edge object
      const Graph* graph1_; // graph pointer associated with node 1 of this edge object
      size_type index2_; //index associated with node 2 of this edge object
      const Graph* graph2_; // graph pointer associated with node 2 of this edge object
      size_type index_edge_; //index associated with this edge object

      /** Private edge constructor. Takes in index and graph pointers 
      * of two nodes and an edge index . **/

      Edge(size_type  index1, const Graph* graph1, size_type index2, const Graph* graph2, 
            size_type index_edge): index1_(index1),graph1_(graph1), index2_(index2), 
            graph2_(graph2), index_edge_(index_edge) { }
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
    Edge edge(size_type i) const {
      assert( i < num_edges());
      return edges_[i]; //edges_ is a vector containing the "Edges" objects of the graph 
      // I return the Edge object at index i stored in edges_vector of Graph object
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
      for ( size_type i = 0; i < num_edges() ; i++){ 
      //iterate over edges_ vector of graph
        if( (edges_[i].node1() == a && edges_[i].node2() == b) || 
            (edges_[i].node1() == b && edges_[i].node2() == a)) // call the node methods
          return true;
      }
      return false;// the nodes are present in graph but they don't have an edge
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
      /** Old implementation  of add_edges. **/
     /** for(size_type i = 0; i< num_edges(); i++){ iterate over edges to find if edge exists_
        if (edges_[i].node1() == a && edges_[i].node2() == b )//
          return Edge(a.index(),a.graph_,b.index(),b.graph_,i);
        else if (edges_[i].node1() == b && edges_[i].node2() == a ) {
          return Edge(b.index(),b.graph_,a.index(),a.graph_,i); ;//return the edge object
        }
      }**/
      /** New implementation of add_edges. **/
      incident_edges_.resize(num_nodes());
      /** Iterating over the incident edges of node a and 
          checking if the other end has node b.**/
      for(size_type i= 0; i< incident_edges_[a.index()].size();i++){
        size_type edgeIndex=incident_edges_[a.index()][i];//edge index
        if(edge(edgeIndex).node1()==b)
          return Edge(b.index(),b.graph_,a.index(),a.graph_,edgeIndex);
        else if(edge(edgeIndex).node2() ==b)
          return Edge(a.index(),a.graph_,b.index(),b.graph_,edgeIndex);
      }

      num_edges_++; // if no edge exists then create a new edge and increment the counter
      // edges_ vector holds the edge objects which are created using private constructor
      Edge newEdge = Edge(a.index_,a.graph_,b.index_,b.graph_,num_edges_-1);
      edges_.push_back(newEdge);
      //update the incident_edges_ vector
      incident_edges_[a.index()].push_back(num_edges_ - 1); 
      // pushing the index of incident edge
      incident_edges_[b.index()].push_back(num_edges_ -1); 
      // pushing the index of incident edge
      return newEdge; //return the last added edge
    }

    /** Remove all nodes and edges from this graph.
        * @post num_nodes()  == 0 && num_edges() == 0
        *
        * Invalidates all outstanding Node and Edge objects.
    */
    void clear() {
      edges_.clear(); //clears all the vectors and assigns 0 to all the counters
      nodes_.clear();
      num_nodes_=0;
      num_edges_=0;
      incident_edges_.clear();
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
          return Node(indexIt_,graphIt_); //calling the node constructor
        return Node();// return an invalid node if index is greater than num_nodes
      }

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
        if(index_edgeIIt_ < graphIIt_ -> incident_edges_[index_nodeIIt_].size()){
          unsigned edge_index =  graphIIt_ -> incident_edges_[index_nodeIIt_][index_edgeIIt_];
          return graphIIt_ -> edge(edge_index); //call the method to return Edge from index
        }
        else 
          return Edge(); 
          //return an invalid edge if edge index is more than num of incident edges
      }

  /** This method is used for incrementing the incident iterator. **/
      incident_iterator& operator++() {
        if(index_edgeIIt_ == graphIIt_ -> incident_edges_[index_nodeIIt_].size())
          return *this; //not incrementing
        else
          index_edgeIIt_ = index_edgeIIt_ + 1;

        return *this;
      }
  /** This method compares two incident iterators.**/
      bool operator==(const incident_iterator& iit) const {
        return (graphIIt_ == iit.graphIIt_ && index_nodeIIt_ == iit.index_nodeIIt_
                && index_edgeIIt_ == iit.index_edgeIIt_);
      }
      private:
        friend class Graph;
        //Private incident iterator takes in node index,edge index and graph
        IncidentIterator(size_type index_nodeIIt, size_type index_edgeIIt, const Graph* graph):
        index_nodeIIt_(index_nodeIIt), index_edgeIIt_(index_edgeIIt), graphIIt_(graph) {}
        size_type index_nodeIIt_;
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
      Edge operator*() const{
        if(index_ < graph_ -> num_edges())
          return graph_ -> edges_[index_];
        else
          return Edge(); //return an invalid edge object
      }
  /** Method to increment edge iterator.**/
      EdgeIterator& operator++() {
        if(index_ == graph_ -> edges_.size())
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
      Point position;
      node_value_type value;
    };
    // New definition of nodes_;
    std::vector<Internal_Nodes> nodes_;
    std::vector<Edge> edges_; //edges_ vector to hold all information about edges i.e. Edges
    size_type num_nodes_; // counter for number of nodes in graph
    size_type num_edges_; // counter for number of edges in graph
  public:
    std::vector<std::vector<size_type>> incident_edges_;
    // this vector stores the indexes of incident edges for all the nodes
};

#endif





