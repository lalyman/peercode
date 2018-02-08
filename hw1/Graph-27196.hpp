#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/// @file Graph.hpp
/// @brief An undirected graph type
///
/*Note: HW0 content used for HW1 was taken from peer code repository*/

#include <algorithm>
#include <cassert>
#include <set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/// @class Graph
/// @brief A template for 3D undirected graphs.
/// 
/// Users can add and retrieve nodes and edges. Edges are unique (there
/// is at most one edge between any pair of distinct nodes).
/// 

template <typename V>
class Graph {
  
  private:

    // Pre-declare structs.
    struct internal_node;
    struct internal_edge;

  public:

    /// Type of this graph.
    using graph_type = Graph;

    /// Predeclaration of Node type.
    class Node;
    /// Synonym for Node (following STL conventions).
    using node_type = Node;

  using node_value_type = V;

    /// Predeclaration of Edge type.
    class Edge;
    /// Synonym for Edge (following STL conventions).
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

    /// Type of indexes and sizes.
    /// Return type of Graph::Node::index(), Graph::num_nodes(),
    /// Graph::num_edges(), and argument type of Graph::node(size_type)
    using size_type = unsigned;


    // CONSTRUCTORS AND DESTRUCTOR

    /// Construct an empty graph.
    Graph()
        : num_nodes_(0), num_edges_(0),
          nodes_(std::vector<internal_node *>()),
          edges_(std::vector<internal_edge *>()) {
    }

    /// Default destructor
    ~Graph() {

      // Since edges_ and nodes_ store pointers to structs on the heap,
      // we need to free that memory.

      for (size_type i = 0; i < num_nodes_; i++) {
        delete nodes_[i];
      }
      for (size_type i = 0; i < num_edges_; i++) {
        delete edges_[i];
      }
    }

    // NODES

    /// @class Graph::Node
    /// @brief Class representing the graph's nodes.
    /// 
    /// Node objects are used to access info about the Graph's nodes.
    /// 
  class Node : private totally_ordered<Node> {
      
      public:
        
        /// Construct an invalid node.
        /// 
        /// Valid nodes are obtained from the Graph class, but it is
        /// occasionally useful to declare an @i invalid node, and
        /// assign a valid node to it later. For example:
        //
        /// @code
        /// Graph::node_type x;
        /// if (...should pick the first node...)
        ///   x = graph.node(0);
        /// else
        ///   x = some other node using a complicated calculation
        /// do_something(x);
        /// @endcode
        /// 
        Node() {
        }

        /// Return this node's position.
        const Point &position() const {
          return graph_->nodes_[index_]->position;
        }

        /// Return this node's index, a number in the range [0, graph_size).
        size_type index() const {
          return index_;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;
    
    /**
     * @brief   Returns the value assigned to a node
     * @return  Value assigned to parent node
     * @pre     a value has been assigned to parent node
     */
    node_value_type& value() const {
      return graph_->nodes_[index_]_->nodeVal;
    }

    /**
     * @brief   Returns the value assigned to a node
     * @return  Value assigned to parent node
     * @pre     a value has been assigned to the parent node
     */
    const node_value_type& value() const {
      return graph_->nodes_[index_]_->nodeVal;
    }

    /**
     * @brief   Returns the number of edges connected to a node
     * @return  Number of edges connected to parent node
     * @pre     Edges have been added to parent graph
     */
    size_type degree() const { // number of edges incident to a node
      size_type nConnect = 0;
      for (size_type ii = 0; ii = graph_->size(); ii++) {
	if graph_->has_edge(index_,ii) {
	    nConnect = nConnect + 1;
	  }
      }
      return nConnect; //problem 3
    }

    /**
     * @brief   Returns iterator to the first edge connected to a node
     * @return  Iterator to the first edge connected to a node
     * @pre     Edges have been added to parent graph
     */    
    incident_iterator edge_begin() const {
      size_type nConnect = 0;
      while !(graph_->has_edge(index_,nConnect){
	  nConnect = nConnect + 1;
	}
	return incident_iterator(graph_->find_edge(index_,nConnect)); //problem 3
    }

     /**
     * @brief   Returns iterator to one past the last edge connected to a node
     * @return  nullptr
     * @pre     Edges have been added to parent graph
     */
    incident_iterator edge_end() const{
	return incident_iterator(nullptr); //problem 3
    }
    
      
        /// Test whether this node and @a n are equal.
        /// 
        /// Equal nodes have the same graph and the same index.
        /// 
        bool operator==(const Node &n) const {
          // If the corresponding entries in nodes_ for n's index and
          // this's index point to the same internal_node object, then
          // both their graphs and indices must be the same.
          return graph_->nodes_[index_] == graph_->nodes_[n.index()];
        }

        /// Test whether this node is less than @a n in a global order.
        /// 
        /// This ordering function is useful for STL containers such as
        /// std::map<>. It need not have any geometric meaning.
        /// 
        /// The node ordering relation must obey trichotomy: For any
        /// two nodes x and y, exactly one of x == y, x < y, and y < x
        /// is true.
        ///
        bool operator<(const Node &n) const {
          return index_ < n.index();
        }

      private:
        
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Pointer back to the Graph to which the node belongs.
        Graph *graph_;

        // Node's unique index in the set of nodes.
        size_type index_;

        /// Private constructor.
        Node(const Graph *graph, size_type ind)
            : graph_(const_cast<Graph *>(graph)), index_(ind) {
        }
    };

    /// Return the number of nodes in the graph.
    /// 
    /// Complexity: O(1).
    /// 
    size_type size() const {
      return num_nodes_;
    }

    /// Synonym for size().
    size_type num_nodes() const {
      return size();
    }

    /// Add a node to the graph, returning the added node.
    /// @param[in] position The new node's position
    /// @post new num_nodes() == old num_nodes() + 1
    /// @post result_node.index() == old num_nodes()
    /// 
    /// Complexity: O(1) amortized operations.
    /// 
  Node add_node(const Point &position, const node_value_type& nodeVal = node_value_type()) {
      // Just add a new internal_node to the end of nodes_ and
      // increment num_nodes_ after assigning the new node the old
      // value of num_nodes_ for its index.
    nodes_.push_back(new internal_node(position, nodeVal));
    return Node(this, num_nodes_++);
    }

    /// Determine if a Node belongs to this Graph
    /// @return True if @a n is currently a Node of this Graph
    /// 
    /// Complexity: O(1).
    /// 
    bool has_node(const Node& n) const {
      // Checking equality ensures that their Graphs are the same.
      return n == node(n.index());
    }

    /// Return the node with index @a i.
    /// @pre 0 <= @a i < num_nodes()
    /// @post result_node.index() == i
    /// 
    /// Complexity: O(1).
    /// 
    Node node(size_type i) const {
      // Simply create a new instance of the proxy.
      return Node(this, i);
    }

    // EDGES

    /// @class Graph::Edge
    /// @brief Class representing the graph's edges.
    /// 
    /// Edges are order-insensitive pairs of nodes. Two Edges with the
    /// same nodes are considered equal if they connect the same nodes,
    /// in either order.
    /// 
  class Edge : private totally_ordered<Edge> {
      
      public:
        
        /// Construct an invalid Edge.
        Edge() {
        }

        // There's no need to store the node indices in the proxy
        // object since we can easily get them from the internal_edge
        // object defined below.

        /// Return a node of this Edge.
        Node node1() const {
          return graph_->node(graph_->edges_[index_]->one);
        }

        /// Return the other node of this Edge.
        Node node2() const {
          return graph_->node(graph_->edges_[index_]->two);
        }

        /// Test whether this edge and @a e are equal.
        /// 
        /// Equal edges are the same undirected edge between two nodes.
        /// 
        bool operator==(const Edge &e) const {
          // Checking equality ensures that their Graphs are the same.
          return ((node1() == e.node1() and node2() == e.node2())
            or (node1() == e.node2() and node2() == e.node1()));
        }

        /// Test whether this edge is less than @a e in a global order.
        /// 
        /// This ordering function is useful for STL container such as
        /// std::map<>. It need not have any interpretive meaning.
        ///
        bool operator<(const Edge &e) const {
          return index_ < graph_->find_edge(e.node1(), e.node2());
        }

      private:
      
        // Allow Graph to access Edge's private members.
        friend class Graph;

        // Pointer back to the Graph to which the node belongs.
        Graph *graph_;

        // Edge's unique index in the set of edges.
        size_type index_;

        /// Private constructor.
        Edge(const Graph *graph, size_type ind)
            : graph_(const_cast<Graph *>(graph)), index_(ind) {
        }
    };

    /// Return the total number of edges in the graph.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    size_type num_edges() const {
      return num_edges_;
    }

    /// Return the edge with index @a i.
    /// @pre 0 <= @a i < num_edges()
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    Edge edge(size_type i) const {
      return Edge(this, i);
    }

    /// Test whether two nodes are connected by an edge.
    /// @pre @a a and @a b are valid nodes of this graph
    /// @return True if for some @a i, edge(@a i) connects @a a and @a b.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    bool has_edge(const Node &a, const Node &b) const {
      return find_edge(a, b) < num_edges_;
    }

    /// Add an edge to the graph, or return the current edge if it
    /// already exists.
    /// @pre @a a and @a b are distinct valid nodes of this graph
    /// @return an Edge object e with e.node1() == @a a and 
    ///                               e.node2() == @a b
    /// @post has_edge(@a a, @a b) == true
    /// @post If old has_edge(@a a, @a b),
    ///         new num_edges() == old num_edges().
    ///       Else
    ///         new num_edges() == old num_edges() + 1.
    /// 
    /// Can invalidate edge indexes -- in other words, old edge(@a i)
    /// might not equal new edge(@a i). Must not invalidate outstanding
    /// Edge objects.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    Edge add_edge(const Node &a, const Node &b) {

      size_type edge_index = find_edge(a, b);
      if (edge_index < num_edges_) {
        return Edge(this, edge_index);
      }

      // If edge does not exist yet, wire it up in internal structures.
      edges_.push_back(new internal_edge(a.index(), b.index()));
      nodes_[a.index()]->incident_edges->insert(num_edges_);
      nodes_[b.index()]->incident_edges->insert(num_edges_);
      return Edge(this, num_edges_++);
    }

    /// Remove all nodes and edges from this graph.
    /// @post num_nodes() == 0 && num_edges() == 0
    /// 
    /// Invalidates all outstanding Node and Edge objects.
    /// 
    void clear() {
      nodes_.clear();
      edges_.clear();
      num_nodes_ = 0; num_edges_ = 0;
    }

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  template <typename V>
  class NodeIterator {
    NodeIterator(Node<V>* p) : p{p} {}
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * @ brief   Overloads * to dereference the iterator
     * @ return  Node that the iterator points to
     */
    Node<v> operator*() const {
      return graph_->Node(p->index());
    }

     /**
     * @ brief   Overloads ++ to increment the iterator
     * @ return  Node with next highest index
     */
    NodeIterator<V>& operator++() {
      size_type nodeInd = p->index() + 1;
      p = graph_->Node(nodeInd);
      return *this;
    }

     /**
     * @ brief   Overloads == to check equality of node iterators
     * @ return  Bool indicating equality (true) or inequality (false)
     */
    bool operator==(const NodeIterator& x) const {
      return p == x.p;
    }

   private:
    friend class Graph<V>;
    friend class Node<V>;
    // HW1 #2: YOUR CODE HERE
    Graph *graph_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
  * @ brief   Returns iterator that points to first node in graph
  * @ return  Node with index 0
  */
  NodeIterator node_begin() const {
    Node<V>* head = this->node(0);
    return NodeIterator<V>(head);
  }

  /**
  * @ brief   Returns iterator that points one past last node in graph
  * @ return  nullptr
  */
  node_iterator node_end() const {
    return NodeIterator<V>(nullptr);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
    IncidentIterator() : p{nullptr} {}
    Edge* p;
    p1* = p.node1(); // common node
    p2* = p.node2(); // changing node
    //store vector of all nodes that "connect" to p1?  How?
    
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * @ brief   Overloads * to dereference the iterator
     * @ return  Edge that the iterator points to
     */
    Edge operator*() const {
      return graph_->edge(graph_, graph_->find_edge(p1,p2));
    }

     /**
     * @ brief   Overloads ++ to increment the iterator
     * @ return  Connected edge with next highest index
     */
    IncidentIterator& operator++() {
      p2 = graph_->node(p2->index()+1);
      while !(graph_->has_edge(p1, p2) || p2->index() >= graph_->size()) {
	  p2 = graph_->node(p2->index()+1);
	}
      p = graph_->find_edge(p1,p2);
      return *this; //go to the next "p2"
    }
    
     /**
     * @ brief   Overloads == to check equality of incident iterators
     * @ return  Bool indicating equality (true) or inequality (false)
     */
    bool operator==(const IncidentIterator& x) const {
      return p == x.p;
    }

   private:
    friend class Graph;
    friend class Edge;
    friend class Node;

    IncidentIterator(Edge* p) : p{p} {}
    Graph *graph_;
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
    EdgeIterator() : itEdge{nullptr} {}
    
    Edge* itEdge;
    p1 = itEdge->node1();
    p2 = itEdge->node2();
    
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
    
    /**
     * @ brief   Overloads * to dereference the iterator
     * @ return  Edge that the iterator points to
     */
    Edge operator*() const {
      return Edge(p1,p2);
    }
    
     /**
     * @ brief   Overloads ++ to increment the iterator
     * @ return  Edge with next highest index
     */
    EdgeIterator& operator++() {
      itEdge = Edge(itEdge->index()+1);
      return *this;
    }
    
     /**
     * @ brief   Overloads == to check equality of edge iterators
     * @ return  Bool indicating equality (true) or inequality (false)
     */
    bool operator==(const EdgeIterator& x) const {
      return itEdge == x.itEdge;
    }

   private:
    friend class Graph;
    friend class Edge;
    friend class Node;

    Graph* graph_;
    EdgeIterator(Edge* itEdge) : itEdge{itEdge} {}
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    
    /**
     * @brief   Returns iterator to the first edge in the graph
     * @return  Iterator to the edge with index 0
     * @pre     Edges have been added to parent graph
     */   
    edge_iterator edge_begin() const {
      Edge* firstEdge = Edge(graph_,0)
      return edge_iterator(firstEdge);
    }
    
    /**
     * @brief   Returns iterator to the first edge connected to a node
     * @return  nullptr
     * @pre     Edges have been added to parent graph
     */   
    edge_iterator edge_end() const {
      return edge_iterator(nullptr);
    }

  private:

    // Helper function to get an edge's index in edges_.
    // 
    // Searches through the edges incident to Node a for one whose
    // other endpoint is b. If such a node is found, return its index
    // in edges_, otherwise return num_edges_.
    //
    size_type find_edge(const Node &a, const Node &b) const {

      // Get the indices of the edges incident to a.
      std::set<size_type> *a_edges = nodes_[a.index()]->incident_edges;

      // Now search these indices to find the edge whose other endpoint
      // is b using the std::find_if function with a custom predicate.
      // References for find_if with custom predicates:
      //   http://www.cplusplus.com/reference/algorithm/find_if/
      //   https://stackoverflow.com/q/6679096/902812
      std::set<size_type>::iterator it = std::find_if(
        a_edges->begin(), a_edges->end(), 

        // Use a lambda function for the custom predicate so that we
        // can access the node b and the member edges_. References for
        // using lambda functions:
        //   https://stackoverflow.com/q/7627098/902812
        [this, b](size_type edge_index) {
          return (
            // Either endpoint can be b, so check both. Since we can
            // assume a and b are from this graph, we just need to
            // check indices.
            edges_[edge_index]->one == b.index()
            or edges_[edge_index]->two == b.index()); } );
      
      // Check whether such a node was found.
      if (it == a_edges->end()) {
        return num_edges_;
      }
      // Can dereference an iterator of std::set<size_type> to get a 
      // size_type (as though it were a size_type *). Reference for
      // using iterators:
      //   https://www.cprogramming.com/tutorial/stl/iterators.html
      return *it;
    }

    // This graph's internals are based on the description in the last
    // bullet of this section from Wikipedia on adjacency lists:
    // https://en.wikipedia.org/wiki/Adjacency_list#Implementation_details

    // Internal struct to bundle node data (Point) with a set of
    // indices of edges to which it's connected.
  template <typename V>
    struct internal_node {

      // Constructor makes pointer to set that's stored on the heap.
      internal_node(const Point &pos, V &nVal)
	: position(pos), nodeVal(nVal), incident_edges(new std::set<size_type>()) {
      }

      // Free the heap memory.
      ~internal_node() {
        delete incident_edges;
      }

      // Store node data.
      const Point position;
      V  nodeVal;

      // Store a pointer to the set instead of the entire set in order
      // to keep the amount of memory on the stack small.
      std::set<size_type> *incident_edges; 
    };

    // Internal struct to store an edge.
    struct internal_edge {

      internal_edge(size_type ind1, size_type ind2)
          : one(ind1), two(ind2) {
      }
      
      // Keep pointers to node endpoints (as indices).
      size_type one;
      size_type two;
    };
    
    // Explicitly track num_nodes_ and num_edges_ with member variables
    // to avoid casting results from nodes_.size() or edges_.size().
    size_type num_nodes_, num_edges_;

    // Positions associated with each node (access by index).
    std::vector<internal_node *> nodes_;
    std::vector<internal_edge *> edges_;
};

#endif // CME212_GRAPH_HPP
