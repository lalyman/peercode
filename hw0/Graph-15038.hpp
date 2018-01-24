#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>

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
  std::map<unsigned, Point> internal_nodes;
  unsigned size_,edge_size;
  std::map<unsigned, std::vector<unsigned>> internal_edges;
  std::map<unsigned, std::unordered_set<unsigned>> connectivity ;
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
  Graph():
		size_(0),edge_size(0){
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
    Node() : node_id(0), graph_() {
      // HW0: YOUR CODE HERE
	
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->internal_nodes[node_id];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
	return node_id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
	if(n.node_id==this->node_id && n.graph_==this->graph_ ){
		return true;
	}
      return false;
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
	if(n.index() < this->index()){
		return true;
	}
      return false;
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
    Node(size_type node_id_ ,const graph_type* graph)
	:node_id(node_id_),graph_(const_cast<graph_type*>(graph)){
	}
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
       //size_type new_size = size() + 1;
       size_ ++;
       size_type new_node_id = size_ - 1;
       internal_nodes.insert(std::pair<unsigned,Point> (new_node_id,position));
       connectivity.insert(std::pair<size_type , std::unordered_set<unsigned>> (new_node_id,{}));
       node_type new_node(new_node_id,this);
    return new_node;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
	for(size_type i=0; i<size(); i++){
		if(n.node_id == i){
			return true;
		}
	}
	
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
    //(void) i;             // Quiet compiler warning
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge(): edge_id(0),graph_e(){
      // HW0: YOUR CODE HERE
    }

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
	if(this->node2() == e.node2() && this->node1() == e.node1()){
		return true;
	}
                 // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(this->edge_id > e.edge_id){
		return true;
      }	           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type node1_id, node2_id, edge_id;
    graph_type* graph_e;
    Edge(size_type edge_id_ ,const graph_type* graph_e_)
        :edge_id(edge_id_),graph_e(const_cast<graph_type*>(graph_e_)){
    }
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
    assert(i < num_edges());
    edge_type new_edge(i,this);
    std::vector<size_type> nodes=internal_edges.at(i);
    new_edge.node1_id = nodes[0];
    new_edge.node2_id = nodes[1];
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
	assert(a.node_id < size_ && b.node_id < size_);
	//std::unordered_set<size_type> set1 = {a.node_id,b.node_id};
	//for (size_type i=0; i< size();++i){
		//std::unordered_set<size_type> set2 = internal_edges.at(i);
               
	       if(connectivity.at(b.node_id).count(a.node_id)){
			return true;
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
	assert(a.node_id < size_ && b.node_id < size_);
        size_type new_edge_id = edge_size - 1;
	edge_type new_edge(new_edge_id,this);
		if(!has_edge(a,b)){
			edge_size++;
			new_edge.edge_id++;
                        if (a < b) {
				new_edge.node1_id = a.node_id;
				new_edge.node2_id = b.node_id;
			}
			else{
				 new_edge.node2_id = a.node_id;
 				 new_edge.node1_id = b.node_id;
			}
			std::vector<size_type> nodes={new_edge.node1_id,new_edge.node2_id};
			internal_edges.insert(std::pair<size_type, std::vector<size_type>> (new_edge.edge_id ,  nodes));
                        connectivity[new_edge.node1_id].emplace(new_edge.node2_id);
			connectivity[new_edge.node2_id].emplace(new_edge.node1_id);
		}
	

       // Quiet compiler warning
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
    internal_edges.erase(internal_edges.begin(),internal_edges.end());
    connectivity.erase(connectivity.begin(),connectivity.end());
    // HW0: YOUR CODE HERE
  }


};

#endif // CME212_GRAPH_HPP
