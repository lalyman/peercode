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

    struct internal_node;
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
  Graph() {
    // HW0: YOUR CODE HERE
 	std::vector<internal_node> empty_v_n;
	nodes = empty_v_n;
        std::vector<internal_edge> empty_v_e;
	edges = empty_v_e;
         num_nodes_ = 0;
         num_edges_ = 0;
         next_index = 0;
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
graph_pointer = nullptr;
            n_id = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
return fetch().position_;
      
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
 return fetch().n_id;
      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
            if (this->n_id == n.n_id){
                return true;	
            } else {
                return false;
            }
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
     (void) n;           // Quiet compiler warning
            if (this->index() < n.index()){
                return true;	
            } else {
                return false;
            }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
Graph* graph_pointer;
        size_type n_id;
        
        Node(Graph graph, size_type index){
            this->graph_pointer = &graph;
            this->n_id = index;
        }
        
        internal_node& fetch() const {
            for (size_type i = 0; i < graph_pointer->size(); ++i){
                if (graph_pointer->nodes[i].n_id == n_id)
                    return graph_pointer->nodes[i];
                assert(false);
            }
        }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_nodes_;
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
    // Create a new elements array
            internal_node new_node;
		new_node.position_ = position;
		new_node.n_id = next_index;
		nodes.push_back(new_node);
                 ++num_nodes_;
		//this->num_nodes_ = this->num_nodes_ +1;
                ++next_index;
                (void) position;      // Quiet compiler warning
               return Node(*this, next_index-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
            for (size_type i = 0; i < num_nodes_; ++i){
                if (nodes[i].n_id == n.index()){
                    return true;
                }
            }
            (void) n;            // Quiet compiler warning
            return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
            for (size_type i = 0; i < num_nodes_; ++i){
                if (nodes[i].n_id == i){
                    return Node(*this,i);
                }
            }
            (void) i;            // Quiet compiler warning
            return Node();        // Invalid node
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
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
                node_1 = Node();
                node_2 = Node();
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
                return this->node_1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
                return this->node_2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->fetch().e_id < e.fetch().e_id) {
                    return true;
                }
                (void) e;           // Quiet compiler warning
                return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
       if (this->fetch().e_id < e.fetch().e_id) {
                    return true;
                }
                (void) e;           // Quiet compiler warning
                return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
 Node node_1;
            Node node_2;

           Edge(Node n1, Node n2){
                this->node_1 = n1;
                this->node_2 = n2;
            }

            internal_edge& fetch() const {
        Graph* graph_p = this->node_1.graph_pointer;
                for (size_type i = 0; i < graph_p->num_edges(); ++i){
                    if ((graph_p->edges[i].node_1.index() == node_1.index() && graph_p->edges[i].node_2.index() == node_2.index()) || (graph_p->edges[i].node_2.index() == node_1.index() && graph_p->edges[i].node_1.index() == node_2.index())) {
                        return graph_p->edges[i];
                    }
                }
                assert(false);
            }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
            return this->num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
                for (size_type i = 0; i < num_edges_; ++i){
                if (edges[i].e_id == i){
                    return Edge(edges[i].node_1, edges[i].node_2);
                }
            }

            (void) i;             // Quiet compiler warning
            return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
            for (size_type i = 0; i < num_edges_; ++i){
                if ((edges[i].node_1 == a && edges[i].node_2 == b) || (edges[i].node_2 == a && edges[i].node_1 == b) ){
                    return true;
                }
            }
            (void) a; (void) b;   // Quiet compiler warning
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
    // HW0: YOUR CODE HERE
           if (has_edge(a,b)){
                return Edge(a,b);
           } else {
		internal_edge new_edge;
		new_edge.node_1 = a;
		new_edge.node_2 = b;
		new_edge.e_id = next_index;
		edges.push_back(new_edge);
                ++num_edges_;
		// this->num_edges_ = this->num_edges_ +1;
                ++next_index;
                (void) a, (void) b;   // Quiet compiler warning
                return edge(next_index-1);
            }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
            nodes.clear();
	    edges.clear(); 
            num_nodes_ = 0;
            num_edges_ = 0;
            next_index = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
struct internal_node {
            Point position_;   
            size_type n_id;      
        };
        struct internal_edge {
            Node node_1;  
            Node node_2;
            size_type e_id;
        };
        
	size_type num_nodes_;
        size_type num_edges_;
        size_type next_index;
        std::vector<Graph::internal_node> nodes;
        std::vector<Graph::internal_edge> edges;
        
        

};

#endif // CME212_GRAPH_HPP
