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


  /** @class Graph
   * @brief A template for 3D undirected graphs.
   *
   * Users can add and retrieve nodes and edges. Edges are unique (there is at
   * most one edge between any pair of distinct nodes).
   */
  class Graph {
   private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)

    /* Predeclaration of the internal representation of a node. */
    struct internal_node;

    /* Predeclaration of the internal representation of an edge. */
    struct internal_edge;

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

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() : nodes_(), edges_() {
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
    class Node {
     public:
      /** Construct an invalid node.
       * Invalid nodes are represented by a null pointer to the graph
       * and an index of -1. Since index is an unsigned int, -1 is actually
       * the maximum value of the size_type attribute. This should be fine,
       * unless we have a graph with the maximum amount of nodes. In that 
       * case, the graph could overflow anyway. 
       * @endcode
       */
      Node() :  graph_(nullptr), index_(-1) {
        // HW0: YOUR CODE HERE
      }

      /** Return this node's position. */
      const Point& position() const {
        // HW0: YOUR CODE HERE
        // Delegate the query to the Graph class via fetchNode()
        return fetchNode().location;
      }

      /** Return this node's index, a number in the range [0, graph_size). */
      size_type index() const {
        // HW0: YOUR CODE HERE
        return this->index_;
      }

      /** Test whether this node and @a n are equal.
       *
       * Equal nodes have the same graph and the same index.
       */
      bool operator==(const Node& n) const {
        // HW0: YOUR CODE HERE
        // Compare graph pointers and indexes.
        return this->graph_ == n.graph_ && this->index_ == n.index_;
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
        return this->index_ < n.index_;
      }

     private:
       // Allow Graph to access Node's private member data and functions.
       friend class Graph;
       // HW0: YOUR CODE HERE
       // Use this space to declare private data members and methods for Node
       // that will not be visible to users, but may be useful within Graph.
       // i.e. Graph needs a way to construct valid Node objects

       /* Pointer to graph to which this node belongs. */
       Graph* graph_;

       /* Unique index identifying this node. */    
       size_type index_;

       /** Private constructor to create Node objects.
        * 
        * Node objects should be created by the user with the add_node method 
        * of the Graph class.
        */
       Node(const Graph* newgraph, size_type id) 
         : graph_(const_cast<Graph*>(newgraph)), index_(id) {
       }

       /** Find the node data represented by this Node.
        *       
        * This method accesses the graph and searches for the 
        * internal representation whose index matches the 
        * identifier for this node.
        */
       internal_node& fetchNode() const {
         // Make sure that the node is in the graph
         assert(index_ < graph_->size());
         
         // Return the node whose index matches the caller index
         return graph_->nodes_[index_];
       }
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
      // HW0: YOUR CODE HERE
      return this->nodes_.size();
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

      // Create a new internal_node with the given position.
      internal_node new_internalnode = {position};

      // Add the internal_node to the node Vector.
      this->nodes_.push_back(new_internalnode);

      // Return a Node object with this graph's pointer and the
      // new Node index.
      return Node(this,this->nodes_.size()-1); 
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
      // HW0: YOUR CODE HERE
      // A Node belongs to this Graph if it points to this Graph
      // and its index is valid, i.e. it's less than the maximum
      // size of the node array.
      return n.graph_ == this && n.index_ < this->nodes_.size();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
      // HW0: YOUR CODE HERE
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
    class Edge {
     public:
      /** Construct an invalid Edge.
       * Invalid nodes are represented by a null pointer to the graph
       * and an index of -1. Since index is an unsigned int, -1 is actually
       * the maximum value of the size_type attribute. This should be fine,
       * unless we have a graph with the maximum amount of nodes. In that 
       * case, the graph could overflow anyway. 
       * @endcode
       */
      Edge() : graph_(nullptr), index_(-1) {
        // HW0: YOUR CODE HERE
      }

      /** Return a node of this Edge */
      Node node1() const {
        // HW0: YOUR CODE HERE
        // Return a Node object pointing to this Graph. To obtain its index,
        // delegate the query to the Graph, which will return the index of one
        // of the nodes this edge links.
        return Node(this->graph_,fetchEdge().node1_id);
      }

      /** Return the other node of this Edge */
      Node node2() const {
        // HW0: YOUR CODE HERE
        // Return a Node object pointing to this Graph. To obtain its index,
        // delegate the query to the Graph, which will return the index of one
        // of the nodes this edge links.
        return Node(this->graph_,fetchEdge().node2_id);
      }

      /** Test whether this edge and @a e are equal.
       *
       * Equal edges represent the same undirected edge between two nodes.
       */
      bool operator==(const Edge& e) const {
        // Equal edges should belong to the same Graph. Also, Edges are 
        // undirected, so we need to check both directions.
        return  this->graph_ == e.graph_ && 
                ((fetchEdge().node1_id == e.fetchEdge().node1_id &&
                  fetchEdge().node2_id == e.fetchEdge().node2_id) ||
                 (fetchEdge().node2_id == e.fetchEdge().node1_id &&
                  fetchEdge().node1_id == e.fetchEdge().node2_id));
      }

      /** Test whether this edge is less than @a e in a global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any interpretive meaning.
       */
      bool operator<(const Edge& e) const {
        return this->index_ < e.index_;
      }

     private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;
      // HW0: YOUR CODE HERE
      // Use this space to declare private data members and methods for Edge
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Edge objects

      /* Pointer to graph to which this node belongs. */
      Graph* graph_;

      /* Index of the edge. */
      size_type index_;

      /** Private constructor to create Edge objects.
       * 
       * Edge objects should be created by the user with the add_edge method 
       * of the Graph class.
       */
      Edge(const Graph* newgraph, size_type id) 
        : graph_(const_cast<Graph*>(newgraph)), index_(id) {
      }

      /** Find the edge data represented by this Edge.
       *       
       * This method accesses the graph and searches for the 
       * internal representation whose index matches the 
       * identifier for this edge.
       */
      internal_edge& fetchEdge() const {
        // Make sure that the edge is in the graph.
        assert(index_ < graph_->num_edges());
        
        // Return the edge whose index matches the caller index.
        return graph_->edges_[index_];
      }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
      // HW0: YOUR CODE HERE
      return this->edges_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
      // HW0: YOUR CODE HERE
      return Edge(this, i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
      // HW0: YOUR CODE HERE

      // For now, we have not found the edge.
      bool result = false;

      // Loop over all edges looking for the one linking the two nodes.
      // Once we find the edge, we can exit the loop. Complexity should 
      // be O(num_edges()).
      for(size_type i = 0; i < num_edges() && result == false; ++i) {
        result =  (edges_[i].node1_id == a.index() &&
                   edges_[i].node2_id == b.index()) ||
                  (edges_[i].node2_id == a.index() &&
                   edges_[i].node1_id == b.index());
      }

      return result;
    }

    /** Add an edge to the graph, or return the current edge if it already 
    *   exists.
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges()+1.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b) {
      // HW0: YOUR CODE HERE

      // Index of the edge, whether it's an existing edge or a new one.
      size_type idx;
      // For now, we have not found the edge.
      bool result = false;

      // Loop over all edges looking for the one linking the two nodes.
      // Once we find the edge, we can exit the loop. We can't use has_edge() 
      // because we also need the edge's index. Complexity should 
      // be O(num_edges()).
      for(idx = 0; idx < num_edges() && result == false; idx++) {
        result =  (edges_[idx].node1_id == a.index() &&
                   edges_[idx].node2_id == b.index()) ||
                  (edges_[idx].node2_id == a.index() &&
                   edges_[idx].node1_id == b.index());
      }

      // If we did not find the edge, create an internal_edge to represent it
      // and add it to the vector.
      if(idx == num_edges() && result == false) {
        internal_edge new_internaledge = {a.index(),b.index()};
      this->edges_.push_back(new_internaledge);
    }

    return Edge(this,idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Empty the cointainer for nodes and edges.
    this->nodes_.clear();
    this->edges_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /* Structure to represent a node internally (Proxy pattern). */
  struct internal_node {
    Point location;
  };

  /* Structure to represent an edge internally (Proxy pattern). */
  struct internal_edge {
    size_type node1_id;
    size_type node2_id;
  };

  /* Vector containing the internal representations of Node objects. */
  std::vector<internal_node> nodes_;
  /* Vector containing the internal representations of Edge objects. */
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
