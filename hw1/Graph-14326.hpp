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
template <typename V>
class Graph : private totally_ordered<Graph<V>> {
  private:

    struct internal_nodes;
    struct internal_edges;

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

    /** Allow the Nodes to support a user-specified value. */
    using node_value_type = V;

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
    Graph()  
      : nodes_(), edges_(){
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
          return graph_->nodes_[this->uid_].p;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
          return graph_->nodes_[this->uid_].uid;
        }
  
        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Return this node's value. */
        node_value_type& value(){
          return graph_->nodes_[this->uid_].node_value;
        }

        /** Return this node's value. */
        const node_value_type& value() const{
          return graph_->nodes_[this->uid_].node_value;
        } 

        /* Return the degree of this node */ 
        size_type degree() const{
          return graph_->nodes_[this->uid_].edges_in_node.size();
        }

        incident_iterator edge_begin() const{
          return IncidentIterator(graph_, this->uid_, 0);
        }


        incident_iterator edge_end() const{
          return IncidentIterator(graph_, this->uid_, degree());
        }
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
          return this->graph_ == n.graph_ && this->uid_ == n.uid_;

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
          return (this->graph_ == n.graph_ && this->uid_ < n.uid_) || 
            (this->graph_ < n.graph_ && this->uid_ == n.uid_);
        }

      private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Pointer back to the Graph
        Graph* graph_;
        // This nodes's unique identification number
        size_type uid_;
        /** Private Constructor */
        Node(const Graph* graph, size_type uid)
          : graph_(const_cast<Graph*>(graph)), uid_(uid) {
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
     * @param[in] user-specified value for the node
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */

    Node add_node(const Point& position, const node_value_type& node_value = node_value_type()){
      internal_nodes new_node; 
      new_node.p = position;
      new_node.node_value = node_value;
      new_node.uid = size();
      nodes_.push_back(new_node);
      return Node(this, new_node.uid); 
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
      return n.graph_ == this; 
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
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
    class Edge : private totally_ordered<Edge> {
      public:
        /** Construct an invalid Edge. */
        Edge() {
        }

        /** Return a node of this Edge */
        Node node1() const {
          return Node(graph_, this->n1_uid_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
          return Node(graph_, this->n2_uid_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
          return (this->node1() == e.node1() && this->node2() == e.node2() && this->graph_ == e.graph_) ||
            (this->node1() == e.node2() && this->node2() == e.node1() && this->graph_ == e.graph_);  
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
          return this->node1() < e.node1() || (this->node1() == e.node1() && this->node2() < e.node2()) ||
            (this->node1() == e.node1() && this->node2() == e.node2() && this->graph_ < e.graph_);
        }

      private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        // Pointer back to the Graph
        Graph* graph_;

        // Node1 uid; 
        size_type n1_uid_;

        // Node2 uid:
        size_type n2_uid_; 

        /** Private Constructor */
        Edge(const Graph* graph, size_type n1_uid, size_type n2_uid)
          : graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_uid_(n2_uid){
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
      edge_iterator ei =  edge_begin(); 
      for(size_type  k = 0; k < i; k++) {
        ++ei;
      }
      return *ei;

    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const { 
      bool contains_edge;
      for(size_type i : nodes_[a.uid_].edges_in_node) {
        contains_edge = (edge(i).node1() == a && edge(i).node2() == b);
        if (contains_edge){
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
      if (has_edge(a,b)){
        return Edge(this, a.uid_, b.uid_); 
      }
      /* std::cout << "Here" << std::endl; */ 
      internal_edges new_edge; 
      new_edge.n1_uid = a.uid_; 
      new_edge.n2_uid = b.uid_;
      this->nodes_[a.uid_].edges_in_node.push_back(num_edges());
      this->nodes_[b.uid_].edges_in_node.push_back(num_edges()); 
      edges_.push_back(new_edge);
      return Edge(this, new_edge.n1_uid, new_edge.n2_uid); 
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
      nodes_.clear();
      edges_.clear(); 
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
        NodeIterator() {}

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /* Return the node for the iterator */
        Node operator*() const{
          return Node(graph_, this->uid_);
        }

        /* Increment the node iterator */
        NodeIterator& operator++(){
          this->uid_ ++;
          return *this;
        }

        /** Test whether this node iterator  and @a node_iter are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const NodeIterator& node_iter) const{
          return this->graph_ == node_iter.graph_ && this->uid_ == node_iter.uid_;
        }

      private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE

        // Pointer back to the Graph
        Graph* graph_;
        // This nodes's unique identification number
        size_type uid_;
        /** Private Constructor for NodeIterator. */ 
        NodeIterator(const Graph* graph, size_type uid) 
          : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
    }; 

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* Returns the beggining of the iterator for nodes. */
    node_iterator node_begin() const{ 
      return NodeIterator(this, 0);  
    }

    /* Returns the end of the iterator for nodes. */
    node_iterator node_end() const{
      return NodeIterator(this, size()); 
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

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
  
        /* Return the incident edge for the iterator */
        Edge operator*() const{
          return graph_->edge(graph_->nodes_[this->uid_].edges_in_node[edge_position_]);
        }; 

        /* Increment the incident edge iterator */
        IncidentIterator& operator++(){
          this->edge_position_ ++;
          return *this;
        }

        /** Test whether this incident edge iterator and @a incident_iter are equal.
         *
         * Equal incident edges have the same graph and the same index.
         */
        bool operator==(const IncidentIterator& incident_iter) const{
          return this->graph_ == incident_iter.graph_ && this->uid_ == incident_iter.uid_ && this->edge_position_ == incident_iter.edge_position_;
        }

      private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE

        // Pointer back to the Graph
        Graph* graph_;
        // This nodes's unique identification number
        size_type uid_;

        // The position of the edge we are iterating over. 
        size_type edge_position_; 


        /** Private Constructor for IncidentIterator. */ 
        IncidentIterator(const Graph* graph, size_type uid, size_type edge_position) 
          : graph_(const_cast<Graph*>(graph)), uid_(uid), edge_position_(edge_position) {}
    };


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

        /* Returns the edge for the iterator */
        Edge operator*() const{
          internal_edges int_edge = graph_->edges_[this->uid_];
          return Edge(graph_, int_edge.n1_uid, int_edge.n2_uid);
        }

        /* Increments the node iterator */
        EdgeIterator& operator++(){
          this->uid_ ++;
          return *this;
        }

        /** Test whether this edge iterator  and @a edge_iter are equal.
         *
         * Equal edges have the same graph and the same index.
         */
        bool operator==(const EdgeIterator& edge_iter) const{
          return this->graph_ == edge_iter.graph_ && this->uid_ == edge_iter.uid_;
        }

      private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE

        // Pointer back to the Graph
        Graph* graph_;
        // This edge's unique identification number
        size_type uid_;
        /** Private Constructor for EdgeIterator. */ 
        EdgeIterator(const Graph* graph, size_type uid) 
          : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* Returns the beggining of the iterator for edges. */
    edge_iterator edge_begin() const{ 
      return EdgeIterator(this, 0);  
    }

    /* Returns the end of the iterator for edges. */
    edge_iterator edge_end() const{
      return EdgeIterator(this, edges_.size()); 
    }

  private:

    /* Internal type for nodes */
    struct internal_nodes {
      Point p;   // The point held by a node. 
      size_type uid; // The unique identifcation for a node.
      node_value_type node_value; // Value of the node.
      std::vector<size_type> edges_in_node; // stores the edges in the node;  
    };
  
    /* Internal type for elements */
    struct internal_edges{
      size_type n1_uid; // Node 1 for this edge. 
      size_type n2_uid; // Node 2 for this edge. 
    };

    std::vector<internal_nodes> nodes_; // vector containing the nodes.
    std::vector<internal_edges> edges_; // vector containing the edges. 

};

#endif // CME212_GRAPH_HPP
