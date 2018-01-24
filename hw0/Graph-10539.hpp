#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

 //NOTE: Many of the implementations in this proxy design pattern
 // mimic those found in proxy_example.cpp. Therefore, similarities between
 // proxy_example.cpp may be found.

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
  //Creating the structures for the internals of the graph classes
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

  /** Construct an empty graph.
  In this case, we also initialize the vector of internal_node and the
  vector of internal_edge here with this constructor (private members of
  the Graph class)
  **/
  Graph()
      : graph_nodes(), graph_edges() {
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

    /*In this case, the invalid node simply initializes nothing (meaning
    that graph_ points to nothing and uid_ is not initialized) and thus the
    constructor simply constructs the Node.
    */
    Node() {
    }

    /** Return this node's position. */
    // We use a private helper function, similar to that seen
    // in proxy_example.cpp to fetch the node position.
    const Point& position() const {
      return fetch_node().node_pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //Checking to see if <this> has the same graph and index as @a n
      if (n.graph_ == this->graph_ && n.uid_ == this->uid_)
        return true;

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
      //Checking to see if <this> has the same graph and less than @a n
      //in global order
      if (n.graph_ == this->graph_ && this->uid_ < n.uid_)
        return true;

      return false;
    }

   private:
    //Create the pointer to the Graph class to act as a proxy to the
    //internals of Graph uid_ here is the index or unique identifier for each
    //Node object. Note that in internal_node, there is no index. Node is the
    //one storing index here. The functions above use uid_ as index.
    graph_type* graph_;
    size_type uid_;

    //Constructor for a valid node with two arguments: the graph and a uid
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    //The fetch_node helper function used above in the Node public functions
    //This helper function retrieves the internal node referenced in the
    //graph_nodes vector corresponding to the current Node object by matching
    //uid_ (indices)
    internal_node& fetch_node() const {
      //checking to see if it's in bounds
      assert(uid_ >= 0 && uid_ < graph_->size());

      return graph_->graph_nodes.at(uid_);
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  //Number of nodes is simply the size of the vector of internal nodes
  size_type size() const {
    return graph_nodes.size();
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
    //Using the proxy's position argument, we initialize the internal_node
    //variable with the correct position to add to the internal graph_nodes
    //vector to correctly add this new node.

    internal_node newNode;
    newNode.node_pt = position;
    graph_nodes.push_back(newNode);

    //Once we've added it to the vector, we return the new node wih the correct
    //index
    //NOTE: our nodes are zero indexed
    return Node(this, this->num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //Since n.index() >= 0 is a tautology since n.index() is unsigned, we check
    //to see if the index of the node is within the bounds of the size of
    //the graph
    if(n.index() < this->num_nodes())
      return true;

    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    //Since i >= 0 is a tautology since i is unsigned
    //More or less the same logic as has_node()
    if(i < this->num_nodes())
      return Node(this, i);

    //Invalid node
    return Node();
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      //return a Node object with the correct index. Uses the helper function
      //as defined in the private section
      return Node(graph_, this->fetch_edge().source);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //Same as above
      return Node(graph_, this->fetch_edge().dest);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //Check the current edges' source and destination Node indices match
      //those of e's source and destination Node indices. Check both ways
      if ((this->node1().index() == e.node1().index()  &&
           this->node2().index() == e.node2().index()) ||
          (this->node1().index() == e.node2().index()  &&
           this->node2().index() == e.node1().index()))
          return true;

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //Check global ordering through the edge's uid
      if (this->uid_ < e.uid_)
        return true;

      return false;
    }

   private:

    //This is more or less the same as that of the Node's private members. We
    //have a pointer to the referenced graph, and a uid for each Edge object
    graph_type* graph_;
    size_type uid_;
    // Allow Graph to access Edge's private member data and functions.

    //We have a constructor for two arguments: graph and uid
    Edge(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    //Again, a helper function to fetch the corresponding edge in the internal
    //storage with that of <this> object.
    internal_edge& fetch_edge() const {
      assert(uid_ >= 0 && uid_ < graph_->num_edges());

      return graph_->graph_edges.at(uid_);
    }

    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  //Same as the num_nodes except counting the edges in the graph_edges vector
  size_type num_edges() const {
    return graph_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  //Here we return the Edge that matches the index that is passed in.
  Edge edge(size_type i) const {
    //Validate to see if index is in range and return the Edge object.
    if(i < this->num_edges())
      return Edge(this, i);

    //Otherwise, return an invalid Edge
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool flag = false;

    //Iterate through all of the edges. If we find an edge where the source
    //and destination indices match those of the indices of the nodes passed
    //in, then we have found a match and we toggle the flag to true.
    for(size_type i = 0; i < graph_edges.size(); i++) {
      if ((a.index() == graph_edges[i].source &&
           b.index() == graph_edges[i].dest)  ||
          (b.index() == graph_edges[i].source &&
           a.index() == graph_edges[i].dest))
        flag = true;
    }

    //Returns depending on if the toggle was triggered.
    return flag;
  }

  /** Add an edge to the graph, or return the current edge if it already exists
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

    //This block of code is more or less the same code as the has_edge function
    //If we find an existing edge, simply return it (we can't reuse has_edge)
    //here since we need to know the uid of the object.
    size_type counter;
    for(counter = 0; counter < graph_edges.size(); counter++) {
      if ((a.index() == graph_edges[counter].source &&
           b.index() == graph_edges[counter].dest)  ||
          (b.index() == graph_edges[counter].source &&
           a.index() == graph_edges[counter].dest))
          return Edge(this, counter);
    }

    //If the edge was not found, then we need to add it. We add it by
    //initializing with a new variable, setting the source and dest values
    //and appending it to our graph_edges vector. This way, we update this in
    //memory.
    internal_edge new_edge;
    new_edge.source = a.index();
    new_edge.dest = b.index();
    graph_edges.push_back(new_edge);

    //Ensuring that the edge has been added
    if (has_edge(a, b)) {
      return Edge(this, counter);
    } else {
      return Edge();        // Invalid Edge
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  //Clearing here simply means flushing the vectors of any stored objects.
  //We leave the destruction of these objects to the destructor
  void clear() {
    graph_nodes.clear();
    graph_edges.clear();
  }

 private:
  //internal_node is the internal struct that we will use to store the nodes
  //It contains simply a Point object to store the positions we will be using.
  struct internal_node {
    Point node_pt;
  };

  //internal_edge is the internal struct that we will use to store the edges
  //It contains the Node index of the first node and the index of the
  //second node names source and dest respectively.
  struct internal_edge {
    size_type source;
    size_type dest;
  };

  //Internal STL container that we will use to store the actual nodes
  //and edges.
  std::vector<internal_node> graph_nodes;
  std::vector<internal_edge> graph_edges;
};

#endif // CME212_GRAPH_HPP
