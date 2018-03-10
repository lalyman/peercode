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
template<typename V, typename E>
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

  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_value_type = E;

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
  Graph() {
    // HW0: YOUR CODE HERE
 	
         num_nodes_ = 0;
         num_edges_ = 0;
         
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
      // HW0: YOUR CODE HERE
          graph_pointer = nullptr;
            n_id = 0;
    }

    /** Return this node's position. */
     Point& position() const {
      // HW0: YOUR CODE HERE
return graph_pointer->nodes[n_id].position_;
      
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
 return n_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
	return this->graph_pointer->nodes[n_id].n_value_;

	}

     const node_value_type& value() const {
	return (const node_value_type&)this->graph_pointer->nodes[n_id].n_value_;
	}

    size_type degree() const {
	return this->graph_pointer->adjacency_matrix[this->n_id].size();
	}
    
	incident_iterator edge_begin() const {
	return IncidentIterator(this->graph_pointer, this->n_id, 0);
	}

    incident_iterator edge_end() const {
	return IncidentIterator(this->graph_pointer, this->n_id, this->degree());
	}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
            if (this->graph_pointer == n.graph_pointer && this->n_id == n.n_id){
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
	//assert(this->graph_pointer == n.graph_pointer);
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
        
        Node(Graph* graph, size_type index){
            this->graph_pointer = graph;
            this->n_id = index;
        }
        
        /*internal_node& fetch() const {
            //for (size_type i = 0; i < graph_pointer->size(); ++i){
                //if (graph_pointer->nodes[i].n_id_== this->n_id){
                    return graph_pointer->nodes[n_id];
		//}
            //}
	//assert(false);
        }*/
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
  Node add_node(const Point &position, const node_value_type &value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Create a new elements array
            internal_node new_node;
		new_node.position_ = position;
		new_node.n_id_ = num_nodes_;
		new_node.n_value_ = value; 
		nodes.push_back(new_node);
		adjacency_matrix.push_back(std::vector<std::pair<size_type, edge_value_type>>());
                 ++num_nodes_;     
               return node(num_nodes_-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
	if (n.index()+1 <= num_nodes_ && n.graph_pointer == this) {
	return true;
	} else {
	return false;
	}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
	return Node(const_cast<Graph*>(this), i);     
  }

  
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
	
	double length() const {
	return norm(node_1.position() - node_2.position());
	}

	edge_value_type & value () {
	for (auto it = node_1.edge_begin(); it != node_1.edge_end(); ++it) 
		if ((*it).node2() == node_2) 
			return node_1.graph_pointer->adjacency_matrix[node_1.index()][it.edge_position].second;
	return node_1.graph_pointer->adjacency_matrix[node_1.index()][0].second;	
	}

	const edge_value_type & value () const {
	for (auto it = node_1.edge_begin(); it != node_1.edge_end(); ++it) 
		if ((*it).node2() == node_2) 
			return node_1.graph_pointer->adjacency_matrix[node_1.index()][it.edge_position].second;
	return node_1.graph_pointer->adjacency_matrix[node_1.index()][0].second;	
	}


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->node_1 == e.node_1 && this->node_2 == e.node_2) | (this->node_1 == e.node_2 && this->node_2 == e.node_1)){
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
      if (std::min(this->node_1.n_id, this->node_2.n_id) < std::min(e.node_1.n_id, e.node_2.n_id) )
	return true;
	else if (std::min(this->node_1.n_id, this->node_2.n_id) > std::min(e.node_1.n_id, e.node_2.n_id))
	return false;
	else if (std::max(e.node_1.n_id, e.node_2.n_id) <= std::max(e.node_1.n_id, e.node_2.n_id))
	return true;
	else return false;
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
	size_type counter = 0;
	for (EdgeIterator ni = this->edge_begin(); ni != this->edge_end(); ++ ni) {
	if (ni.node_position < ni.graph_pointer->adjacency_matrix[ni.node_position][ni.edge_position].first) { //to count everything once //
	if (counter == i)
	return Edge(Node(const_cast<Graph*>(ni.graph_pointer), ni.node_position),Node(const_cast<Graph*> (ni.graph_pointer),ni.graph_pointer->adjacency_matrix[ni.node_position][ni.edge_position].first));
	else counter++;
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
    // HW0: YOUR CODE HERE
	assert(a.graph_pointer == b.graph_pointer);
	for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
	if ((*it).node2().index() == b.n_id)
	return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type &value = edge_value_type()) {
    // HW0: YOUR CODE HERE
           if (has_edge(a,b)){
                return Edge(a,b);
           } else {
		assert(a.graph_pointer == b.graph_pointer);
		std::pair<size_type,edge_value_type> bb(b.n_id, value);
		std::pair<size_type,edge_value_type> aa(a.n_id, value);
		a.graph_pointer->adjacency_matrix[a.n_id].push_back(bb);
		b.graph_pointer->adjacency_matrix[b.n_id].push_back(aa);
                ++num_edges_;
                return Edge(a, b);
            }
  }

//HW2 remove functions

/** Remove a node from the graph and all its edges and return the old index of the removed node
   * @param[in]     n node to be removed
   * @pre @a n is a node of this graph
   * @return old index of deleted node
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() = old num_edges() - old n.degree()
   *
   * Complexity: No more than O(num_nodes()), hopefully less, since we assume graph is sparse
   */
size_type remove_node (const Node& n) {
assert(this->has_node(n));
while(!adjacency_matrix[n.index()].empty()) 
	remove_edge(*n.edge_begin());
nodes.erase(nodes.begin() + n.index());
adjacency_matrix.erase(adjacency_matrix.begin() + n.index());
num_nodes_--;
for (size_type i = 0; i < num_nodes_; ++i) {
	if (nodes[i].n_id_ > n.index())
		nodes[i].n_id_--;
	for(size_type j = 0; j < adjacency_matrix[i].size(); ++j) 
		if(adjacency_matrix[i][j].first > n.index())
			adjacency_matrix[i][j].first--;
	}
return n.index();
}

/** Remove a node from the graph and all its edges and return the old index of the removed node
   * @param[in]     n node to be removed
   * @pre @a n is a node of this graph
   * @return node iterator pointing to the same position but different node
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() = old num_edges() - old n.degree()
   *
   * Complexity: No more than O(num_nodes()), hopefully less, since we assume graph is sparse
   */
node_iterator remove_node (node_iterator n_it) {
remove_node(*n_it);
return n_it;
}

/** Remove an edge from the graph and return true if remove worked (i.e. the edge we wanted to remove
   * actually exists, else it returns false
   * @param[in]     n1 node1() of the edge to be removed
   * @param[in]     n2 node2() of the edge to be removed
   * @pre @a n1 and @a n2 are nodes of the graph
   * @return true if has_edge(n1, n2) == true else return false
   * @post new num_edges() = old num_edges() - -1
   *
   * Complexity: O(1) since we are searching incident edges of two nodes and graph is sparse
   */
size_type remove_edge (const Node& n1, const Node& n2) {
for (size_type j = 0; j < adjacency_matrix[n1.index()].size(); ++j) {
	if (adjacency_matrix[n1.index()][j].first == n2.index()) {
		adjacency_matrix[n1.index()].erase(adjacency_matrix[n1.index()].begin() + j);
		break;
		} }
for (size_type j = 0; j < adjacency_matrix[n2.index()].size(); ++j) {
	if (adjacency_matrix[n2.index()][j].first == n1.index()) {
		adjacency_matrix[n2.index()].erase(adjacency_matrix[n2.index()].begin() + j);
		--num_edges_;
		return 1;
		} }
return 0;
}

size_type remove_edge (const Edge& e) {
return remove_edge(e.node1(), e.node2());
}

edge_iterator remove_edge (edge_iterator e_it) {
return remove_edge(*e_it);
}


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
            nodes.clear();
	    adjacency_matrix.clear();
            num_nodes_ = 0;
            num_edges_ = 0;
           
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
	assert(n_position <= graph_pointer->num_nodes_);
	return Node(const_cast<Graph*>(graph_pointer), n_position);  
	}

    NodeIterator& operator++() {
	n_position++;
	return *this;
	}

    bool operator==(const NodeIterator& nodeiter) const {
	if (nodeiter.graph_pointer == this->graph_pointer && nodeiter.n_position == this->n_position)
	return true;
	else return false;
	}

   private:
    friend class Graph;
	// HW1 #2: YOUR CODE HERE
	const Graph *graph_pointer;
	size_type n_position;

// private constructor
	NodeIterator(const Graph *graph_p, size_type n_p) {
	this->graph_pointer = graph_p;
	this->n_position = n_p;
	}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
	return NodeIterator(this,0);
	}

  node_iterator node_end() const {
	return NodeIterator(this,this->num_nodes_);
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
    IncidentIterator() {//
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
	return Edge(Node(const_cast<Graph*>(this->graph_pointer), this->node_id),Node(const_cast<Graph*>(this->graph_pointer), this->graph_pointer->adjacency_matrix[node_id][edge_position].first));

	}

    IncidentIterator& operator++() {
	edge_position++;
	return *this;
	}

    bool operator==(const IncidentIterator& II) const {
	if (this->graph_pointer == II.graph_pointer && this->node_id == II.node_id && this->edge_position == II.edge_position)
	return true;
	else return false;
	}

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
	const Graph *graph_pointer;
	size_type node_id;
	size_type edge_position;

	// private constructor
	IncidentIterator(const Graph *graph_p, size_type nodeid, size_type pos) {
	this->graph_pointer = graph_p;
	this->node_id = nodeid;
	this->edge_position = pos;
	}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {//
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Edge operator*() const {
	return Edge(Node(const_cast<Graph*>(graph_pointer), node_position),Node(const_cast<Graph*>(graph_pointer), graph_pointer->adjacency_matrix[node_position][edge_position].first));
	}
	
     EdgeIterator& operator++() {
	if (edge_position < graph_pointer->adjacency_matrix[node_position].size() - 1)
	edge_position++;
	else {
	node_position++;
	edge_position = 0;
	}
	return *this;
	}
    
	bool operator==(const EdgeIterator& EI) const {
	if (this->graph_pointer == EI.graph_pointer && this->node_position == EI.node_position && this->edge_position == EI.edge_position)
	return true;
	else return false;
	}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	const Graph *graph_pointer;
	size_type node_position;
	size_type edge_position;
	
	EdgeIterator(const Graph* graphp, size_type np, size_type ep) {
	this->graph_pointer = graphp;
	this->node_position = np;
	this->edge_position = ep;
	}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
	return EdgeIterator(this, 0, 0);
	}
	
  edge_iterator edge_end() const {
	return EdgeIterator(this,num_nodes_, 0);
	}

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
struct internal_node {
            Point position_;   
            size_type n_id_; 
	    node_value_type n_value_;    
        };
        
	size_type num_nodes_;
        size_type num_edges_;
        std::vector<Graph::internal_node> nodes;
        std::vector<std::vector<std::pair<size_type, edge_value_type>>> adjacency_matrix;
        
        

};

#endif // CME212_GRAPH_HPP
