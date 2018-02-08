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
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
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

  using hash_map = std::unordered_map<size_type, size_type>;
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /**
  * @brief Construct an empty graph.
  *
  * @param none
  * @return Graph Object
  *
  * @pre none
  * @post Graph object is created
  *
  * In this case, we also initialize the vector of internal_node and the
  * vector of internal_edge here with this constructor (private members of
  * the Graph class)
  **/
  Graph()
      : graph_nodes(), graph_edges() {
  }

  /**
  * @brief Default destructor
  *
  * @param none
  *
  * @pre Graph object exists
  * @post Graph object is destroyed
  **/
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

   /**
   * @brief Construct an invalid Node.
   *
   * @param none
   * @return Node object
   *
   * @pre none
   * @post Node object is created
   *
   * In this case, the invalid node simply initializes nothing (meaning
   * that graph_ points to nothing and uid_ is not initialized) and thus the
   * constructor simply constructs the Node.
   **/
    Node() {
    }

    /**
    * @brief Return this node's position.
    *
    * @param none
    * @return A Point object of the node with position.
    *
    * @pre Node object exists with valid position
    * @post A Point object is returned
    *
    * We use a private helper function, similar to that seen
    * in proxy_example.cpp to fetch the node position.
    **/
    const Point& position() const {
      return fetch_node().node_pt;
    }

    /**
    * @brief Return this node's index
    *
    * @param none
    * @return A unsigned integer with the index of the node.
    *
    * @pre Node object exists with valid index
    * @post A unsigned integer in the range [0, graph_size) is returned
    **/
    size_type index() const {
      return this->uid_;
    }

    /**
    * @brief Return the reference to the node's value
    *
    * @param none
    * @return A node_value_type value in the Node object
    *
    * @pre Node object exists with valid node_value_type value
    * @post A node_value_type value associated with the node is returned
    *
    * This function acts as a sort of a setter since the value is not const.
    **/
    node_value_type& value() {
      return fetch_node().val;
    }

    /**
    * @brief Return a constant reference to the node's constant value
    *
    * @param none
    * @return A const node_value_type value in the Node object
    *
    * @pre Node object exists with valid node_value_type value
    * @post A const node_value_type value associated with the node is returned
    *
    * This function acts as a sort of a getter since the value is const.
    **/
    const node_value_type& value() const {
      return fetch_node().val;
    }

    /**
    * @brief Computes the degree of the node
    *
    * @param none
    * @return An unsigned integer degree of the node
    *
    * @pre Node object exists
    * @post The unsigned degree of the Node is returned
    *
    * Retrieves the degree of the node using the IncidentIterator class
    **/
    size_type degree() const {
      Node n = *this;
      size_type counter = 0;

      //Using the IncidentIterator class to iterate over
      for(auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        counter++;
      }
      return counter;
    }

    /**
    * @brief Returns an incident iterator pointing to the first edge
    *
    * @param none
    * @return An incident_iterator object pointing to the first element
    *
    * @pre this object exists
    * @post An incident_iterator object is returned pointing to the first edge
    **/
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,
                              graph_->edge_search.find(uid_)->second.begin(),
                              uid_);
    }

    /**
    * @brief Returns an incident iterator pointing to the last edge
    *
    * @param none
    * @return An incident_iterator object pointing to the last element
    *
    * @pre this object exists
    * @post An incident_iterator object is returned pointing to the last edge
    **/
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,
                              graph_->edge_search.find(uid_)->second.end(),
                              uid_);
    }

    /**
    * @brief Test whether this node and @a n are equal.
    *
    * @param[in] n        Node to test equality on
    * @return             true if the two nodes are equal, false otherwise
    *
    * @pre Node object exists
    * @pre n is a valid Node
    * @post If result is true, then this->graph_ == @a n.graph_ and
    *       this->uid_ = @a n.uid_. If the result is false, this does not have
    *       to hold.
    *
    * Equal nodes have the same graph and the same index.
    **/
    bool operator==(const Node& n) const {
      //Checking to see if <this> has the same graph and index as @a n
      if (n.graph_ == this->graph_ && n.uid_ == this->uid_)
        return true;

      return false;
    }

    /**
    * @brief Test whether this node is less than @a n in a global order.
    *
    * @param[in] n        Node to test equality on
    * @return             true if this node is less than @a n, false otherwise
    *
    * @pre Node object exists
    * @pre n is a valid Node
    * @post If result is true, then this->graph_ == n.graph_ and
    *       this->uid_ < n.uid_. If the result is false, this does not have to
    *       hold.
    *
    * This ordering function is useful for STL containers such as
    * std::map<>. It need not have any geometric meaning.
    *
    * The node ordering relation must obey trichotomy: For any two nodes x
    * and y, exactly one of x == y, x < y, and y < x is true.
    **/
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

    /**
    * @brief Constructor for a valid node with two arguments: the graph
    *        and a uid
    *
    * @param[in] graph        Graph object that contains the node
    * @param[in] uid          unsigned integer of node's index/unique id
    * @return                 An Node containing the initialized values
    *
    * @pre graph_ is not a nullptr
    * @pre 0 <= uid < size of the graph
    * @post A Node n such that n.graph_ != nullptr and
    *       0 <= uid_ < size of the graph is returned.
    **/
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    /**
    * @brief Helper function to fetch nodes
    *
    * @param none
    * @return An internal_node object corresponding to the Node
    *
    *
    * @pre graph_ is not a nullptr
    * @pre 0 <= uid_ < size of the graph
    * @pre The vector graph_nodes is not empty
    * @post The result is a valid internal_node struct.
    *
    * The fetch_node helper function used above in the Node public functions
    * This helper function retrieves the internal node referenced in the
    * graph_nodes vector corresponding to the current Node object by matching
    * uid_ (indices)
    **/
    internal_node& fetch_node() const {
      //checking to see if it's in bounds
      assert(uid_ >= 0 && uid_ < graph_->size());

      return graph_->graph_nodes.at(uid_);
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /**
  * @brief Return the number of nodes in the graph.
  *
  * @param none
  * @return the unsigned integer size of the graph (number of nodes)
  *
  *
  * @pre Graph object has been constructed
  * @pre The vector graph_nodes is not empty
  * @post result == the number of nodes in graph_nodes
  *
  * Number of nodes is simply the size of the vector of internal nodes
  * Complexity: O(1).
  **/
  size_type size() const {
    return graph_nodes.size();
  }

  /**
  * @brief Return the number of nodes in the graph.
  *
  * @param none
  * @return the unsigned integer size of the graph (number of nodes)
  *
  *
  * @pre Graph object has been constructed
  * @pre The vector graph_nodes is not empty
  * @post result == the number of nodes in graph_nodes
  *
  * Synonym for size().
  **/
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value  The value stored inside the node
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type ()) {
    //Using the proxy's position argument, we initialize the internal_node
    //variable with the correct position to add to the internal graph_nodes
    //vector to correctly add this new node.

    internal_node newNode;
    newNode.node_pt = position;
    newNode.val = value;
    graph_nodes.push_back(newNode);

    //Once we've added it to the vector, we return the new node wih the correct
    //index
    //NOTE: our nodes are zero indexed
    return Node(this, this->num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @param n   Node to check to see if it belongs in the graph
   * @return True if @a n is currently a Node of this Graph
   *
   * @pre n is a valid Node
   * @post n.index() < this->num_nodes()
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
   * @param[in] Index of node to get
   * @return Node object corresponding to the index passed in
   *
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
  class Edge : private totally_ordered<Edge> {
   public:

    /**
    * @brief Construct an invalid Edge.
    *
    * @param none
    * @return Edge object
    *
    * @pre none
    * @post Edge object is created
    **/
    Edge() {
    }

    /**
    * @brief Return a node of this Edge
    *
    * @param none
    * @return One of the nodes of the edge
    *
    * @pre Edge object exists
    * @post Node that is returned corresponds to one of the nodes in the edge.
    **/
    Node node1() const {
      //return a Node object with the correct index. Uses the helper function
      //as defined in the private section. Direction corresponds to different
      //directions of the edge (either (a,b) or (b,a)).
      if(direction_) {
        return Node(graph_, this->fetch_edge().source);
      }
      else {
        return Node(graph_, this->fetch_edge().dest);
      }
    }

    /**
    * @brief Return the other node of this Edge
    *
    * @param none
    * @return One of the nodes of the edge
    *
    * @pre Edge object exists
    * @post Node that is returned corresponds to the other node in the edge.
    **/
    Node node2() const {
      //Same as above
      if(direction_) {
        return Node(graph_, this->fetch_edge().dest);
      }
      else {
        return Node(graph_, this->fetch_edge().source);
      }
    }

    /**
    * @brief Test whether this edge and @a e are equal.
    *
    * @param e      The edge to compare this object to
    * @return       True if the edges are equal. False otherwise.
    *
    * @pre This object has been initialized
    * @pre @a e is a valid edge
    * @post True if this->node1().index() == e.node1().index()  &&
    *       this->node2().index() == e.node2().index()) and vice versa.
    *       False otherwise
    *
    * Equal edges represent the same undirected edge between two nodes.
    **/
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

    /**
    * @brief Test whether this edge is less than @a e in a global order.
    *
    * @param e      The edge to compare this object to
    * @return       True if this edge is less than @a e. False otherwise.
    *
    * @pre This object has been initialized
    * @pre @a e is a valid edge
    * @post True if this->uid_ < e.uid_. False otherwise
    *
    * This ordering function is useful for STL containers such as
    * std::map<>. It need not have any interpretive meaning.
    **/
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
    bool direction_;
    // Allow Graph to access Edge's private member data and functions.

    /**
    * @brief Constructor for a valid edge with two arguments: the graph
    *        and a uid
    *
    * @param[in] graph        Graph object that contains the edge
    * @param[in] uid          unsigned integer of edge's index/unique id
    * @return                 An Edge containing the initialized values
    *
    * @pre graph_ is not a nullptr
    * @pre 0 <= uid < size of the graph
    * @post A Edge e such that e.graph_ != nullptr and
    *       0 <= uid_ < number of edges is returned.
    **/
    Edge(const graph_type* graph, size_type uid, bool direction)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid),
                 direction_(direction){
    }

    /**
    * @brief Again, a helper function to fetch the corresponding edge
    *        in the internal storage with that of <this> object.
    *
    * @param none
    * @return An internal_edge object corresponding to the Edge
    *
    *
    * @pre graph_ is not a nullptr
    * @pre 0 <= uid_< number of edges
    * @pre The vector graph_edges is not empty
    * @post The result is a valid internal_edge struct.
    **/
    internal_edge& fetch_edge() const {
      assert(uid_ >= 0 && uid_ < graph_->num_edges());

      return graph_->graph_edges.at(uid_);
    }

    friend class Graph;
  };

  /**
   * @brief Return the total number of edges in the graph.
   *
   * @param none
   * @return the unsigned integer number of edges
   *
   *
   * @pre Graph object has been constructed
   * @pre The vector graph_edges is not empty
   * @post result == the number of edges in graph_edges
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return graph_edges.size();
  }

  /**
   * @brief Return the edge with index @a i.
   * @param[in] i Index of the edge to get
   * @return Edge object corresponding to the index passed in.
   *
   * @pre 0 <= @a i < num_edges()
   * @post result is an Edge object corresponding to the index passed in
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Here we return the Edge that matches the index that is passed in.
   */
  Edge edge(size_type i) const {
    //Validate to see if index is in range and return the Edge object.
    if(i < this->num_edges())
      return Edge(this, i, true);

    //Otherwise, return an invalid Edge
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @param[in] a   A Node in the edge
   * @oaram[in] b   The other node in the edge
   * @return True if the edge is in the graph. False otherwise
   *
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool flag = false;

    //Search in the outer map in the edge_search map, the first
    //key corresponds to the first node in the graph. The value of the
    //outer-most map is another map. This map contains all the nodes that @a a
    //connects to with an edge. We can then iterate through this inner map to
    //find the indices of the edge if there is one that connects @a a to @a b.
    //The following code is the nested iteration through the adjacency map.
    auto search = edge_search.find(a.index());
    if(search != edge_search.end()) {
      //
      auto secondSearch = (search->second).find(b.index());
      if(secondSearch != (search->second).end())
        flag = true;
    }

    //Returns depending on if the toggle was triggered.
    return flag;
  }

  /**
   * @brief Add an edge to the graph, or return the current edge
   *        if it already exists
   *
   * @param[in] a   A Node in the edge
   * @oaram[in] b   The other node in the edge
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *
   * @pre @a a and @a b are distinct valid nodes of this graph
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
    //If it has the edge in the graph, return in. Check both directions
    if(has_edge(a, b)) {
      return Edge(this, edge_search[a.index()][b.index()], true);
    }
    if(has_edge(b, a)) {
      return Edge(this, edge_search[a.index()][b.index()], false);
    }
    //If the edge was not found, then we need to add it. We add it by
    //initializing with a new variable, setting the source and dest values
    //and appending it to our graph_edges vector. This way, we update this in
    //memory. In addition, make sure we add it to the edge_search map for ease
    //of search in the future
    internal_edge new_edge;
    new_edge.source = a.index();
    new_edge.dest = b.index();
    graph_edges.push_back(new_edge);

    size_type new_index = graph_edges.size() - 1;
    edge_search[a.index()][b.index()] = new_index;
    edge_search[b.index()][a.index()] = new_index;

    return Edge(this, new_index, true);
  }

  /**
   * @brief Remove all nodes and edges from this graph.
   *
   * @param none
   *
   * @pre Graph object exists
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects. As well as the
   * edge_search map.
   * Clearing here simply means flushing the vectors of any stored objects.
   * We leave the destruction of these objects to the destructor.
   */
  void clear() {
    graph_nodes.clear();
    graph_edges.clear();
    edge_search.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /**
    * @brief Construct an invalid NodeIterator.
    *
    * @param none
    * @return NodeIterator object
    *
    * @pre none
    * @post NodeIterator object is created
    **/
    NodeIterator() {
    }

    /**
    * @brief Operator to deference the iterator and find what node the
    *        iterator is pointing to.
    *
    * @param none
    * @return Node object corresponding to the object the iterator is
    *         pointing to.
    *
    * @pre this object is initialized
    * @post The result is the Node that the iterator is pointing to.
    *
    * We use the iterInd_ which mimics the index of the node to search for the
    * correct node for our purposes.
    **/
    Node operator*() const {
      return Node(graph_, iterInd_);
    }

    /**
    * @brief Operator to increment the iterator
    *
    * @param none
    * @return NodeIterator object that now has an incremented index.
    *
    * @pre this object is initialized
    * @post new iterInd_ = old iterInd_ + 1
    *
    * We use the iterInd_ which mimics the index of the node to search for the
    * correct node for our purposes.
    **/
    NodeIterator& operator++() {
      iterInd_++;
      return *this;
    }

    /**
    * @brief Operator to test equality of the the NodeIterator
    *
    * @param[in] ni the NodeIterator object to compare this to
    * @return True if the NodeIterator objects are the same. False otherwise.
    *
    * @pre this object is initialized
    * @pre ni is a valid NodeIterator object.
    * @post if the result is true, then this->graph_ == ni.graph_ and
    *       this->iterInd_ == ni.iterInd_. If false, then these conditions
    *       do not have to hold.
    **/
    bool operator==(const NodeIterator& ni) const {
      if(this->graph_ == ni.graph_ && this->iterInd_ == ni.iterInd_)
        return true;
      else
        return false;
    }

   private:
    graph_type* graph_;
    size_type iterInd_;

    /**
    * @brief Constructor for a valid NodeIterator with two arguments: the graph
    *        and a uid
    *
    * @param[in] graph        Graph object that contains the node
    * @param[in] iterInd      unsigned integer of index the iterator is
    *                         pointing to
    * @return                 An NodeIterator containing the initialized values
    *
    * @pre graph_ is not a nullptr
    * @pre 0 <= iterInd < size of the graph
    * @post A NodeIterator n such that n.graph_ != nullptr and
    *       0 <= iterInd_ < size of the graph is returned.
    **/
    NodeIterator(const graph_type* graph, size_type iterInd)
        : graph_(const_cast<graph_type*>(graph)), iterInd_(iterInd){
    }

    friend class Graph;
  };

  /**
  * @brief Returns a node_iterator pointing to the first Node in the iterator.
  *
  * @param none
  * @return A node_iterator object pointing to the first element
  *
  * @pre this object exists
  * @post A node_iterator object is returned pointing to the first node
  **/
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
  * @brief Returns a node_iterator pointing to the last Node in the iterator.
  *
  * @param none
  * @return A node_iterator object pointing to the last element
  *
  * @pre this object exists
  * @post A node_iterator object is returned pointing to the last node
  **/
  node_iterator node_end() const {
    return NodeIterator(this, this->size());
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

    /**
    * @brief Construct an invalid IncidentIterator.
    *
    * @param none
    * @return IncidentIterator object
    *
    * @pre none
    * @post IncidentIterator object is created
    **/
    IncidentIterator() {
    }

    /**
    * @brief Operator to deference the iterator and find what node the
    *        iterator is pointing to. Note this only returns incident edges
    *
    * @param none
    * @return Edge object corresponding to the object the iterator is
    *         pointing to.
    *
    * @pre this object is initialized
    * @post The result is the Edge that the iterator is pointing to.
    **/
    Edge operator*() const {
      //For this object, our iterator is a map. Therefore, mapIter is an
      //iterator over all the nodes/keys in the outer map. With each resulting
      //map, we look to see how many of these inner maps have the node with
      //index n_ as the key (meaning the node with index n_ is a destination
      //node). We find the corresponding map and return the edge found.
      if(graph_->graph_edges[mapIter_->second].source == n_) {
        return Edge(graph_, mapIter_->second, true);
      }
      //For the reverse direction
      else {
        return Edge(graph_, mapIter_->second, false);
      }
    }

    /**
    * @brief Operator to increment the iterator
    *
    * @param none
    * @return incident_iterator object that now is pointing to the next object
    *         in the map iterator
    *
    * @pre this object is initialized
    * @post new mapIter_ now points to the next map in the iterator
    **/
    incident_iterator& operator++() {
      mapIter_++;
      return *this;
    }

    /**
    * @brief Operator to test equality of the the IncidentIterator
    *
    * @param[in] iit the IncidentIterator object to compare this to
    * @return True if the IncidentIterator objects are the same.
    *         False otherwise.
    *
    * @pre this object is initialized
    * @pre iit is a valid IncidentIterator object.
    * @post if the result is true, then this->graph_ == ni.graph_ and
    *       this->n_ == ni.n_ and this->mapIter_ == iit.mapIter_.
    *       If false, then these conditions do not have to hold.
    **/
    bool operator==(const incident_iterator& iit) const {
      if(this->graph_ == iit.graph_ && this->mapIter_ == iit.mapIter_ &&
         this->n_ == iit.n_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
     //Note that for this object, we are iterating over maps instead of simple
     //size_type indices. This is so that we can search for the index through
     //iterating over the map. mapIter_ is an iterator over the maps and n_ is
     //the index of the node we are tryiing to find all edges incident to.
     graph_type* graph_;
     hash_map::iterator mapIter_;
     size_type n_;

     /**
     * @brief Constructor for a valid NodeIterator with three arguments: the
     *        graph and the map iterator and the node n that the iterator
     *        is trying to find incident edges to.
     *
     * @param[in] graph   Graph object that contains the node
     * @param[in] mapIter map iterator that will iterator over map objects
     * @param[in] n       index of node we are trying to find all edges
     *                    incident to.
     * @return            An IncidentIterator containing the initialized values
     *
     * @pre graph_ is not a nullptr
     * @pre mapIter is a valid iterator and has been initialized.
     * @pre 0 <= n < size of the graph
     * @post A IncidentIterator iit such that iit.graph_ != nullptr and
     *       0 <= n_ < size of the graph is returned and mapIter is a valid
     *       iterator.
     **/
     IncidentIterator(const graph_type* graph, hash_map::iterator mapIter,
                      size_type n)
         : graph_(const_cast<graph_type*>(graph)), mapIter_(mapIter), n_(n){
     }
    friend class Graph;
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

    /**
    * @brief Construct an invalid EdgeIterator.
    *
    * @param none
    * @return EdgeIterator object
    *
    * @pre none
    * @post EdgeIterator object is created
    **/
    EdgeIterator() {
    }

    /**
    * @brief Operator to deference the iterator and find what edge the
    *        iterator is pointing to.
    *
    * @param none
    * @return Edge object corresponding to the object the iterator is
    *         pointing to.
    *
    * @pre this object is initialized
    * @post The result is the Edge that the iterator is pointing to.
    *
    * We use the iterInd_ which mimics the index of the edge to search for the
    * correct edge for our purposes.
    **/
    Edge operator*() const {
      return Edge(graph_, iterInd_, true);
    }

    /**
    * @brief Operator to increment the iterator
    *
    * @param none
    * @return EdgeIterator object that now has an incremented index.
    *
    * @pre this object is initialized
    * @post new iterInd_ = old iterInd_ + 1
    *
    * We use the iterInd_ which mimics the index of the node to search for the
    * correct node for our purposes.
    **/
    EdgeIterator& operator++() {
      ++iterInd_;
      return *this;
    }

    /**
    * @brief Operator to test equality of the the EdgeIterator
    *
    * @param[in] ni the EdgeIterator object to compare this to
    * @return True if the EdgeIterator objects are the same. False otherwise.
    *
    * @pre this object is initialized
    * @pre @a ni is a valid EdgeIterator object.
    * @post if the result is true, then this->graph_ == ei.graph_ and
    *       this->iterInd_ == ei.iterInd_. If false, then these conditions
    *       do not have to hold.
    **/
    bool operator==(const EdgeIterator& ei) const {
      if(this->graph_ == ei.graph_ && this->iterInd_ == ei.iterInd_) {
        return true;
      }
      else {
        return false;
      }
    }


   private:
     graph_type* graph_;
     size_type iterInd_;

     /**
     * @brief Constructor for a valid EdgeIterator with two arguments: the
     *        graph and a uid
     *
     * @param[in] graph        Graph object that contains the node
     * @param[in] iterInd      unsigned integer of index the iterator is
     *                         pointing to
     * @return                 An EdgeIterator containing the initialized values
     *
     * @pre graph_ is not a nullptr
     * @pre 0 <= iterInd < number of edges
     * @post A EdgeIterator n such that n.graph_ != nullptr and
     *       0 <= iterInd_ < number of edges is returned.
     **/
     EdgeIterator(const graph_type* graph, size_type iterInd)
         : graph_(const_cast<graph_type*>(graph)), iterInd_(iterInd){
     }

    friend class Graph;
  };

  /**
  * @brief Returns a edge_iterator pointing to the first Edge in the iterator.
  *
  * @param none
  * @return A edge_iterator object pointing to the first element
  *
  * @pre this object exists
  * @post A edge_iterator object is returned pointing to the first edge
  **/
  edge_iterator ee_edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
  * @brief Returns a edge_iterator pointing to the last Edge in the iterator.
  *
  * @param none
  * @return A edge_iterator object pointing to the last element
  *
  * @pre this object exists
  * @post A edge_iterator object is returned pointing to the last edge
  **/
  edge_iterator ee_edge_end() const {
    return EdgeIterator(this, num_edges());
  }


 private:
  //internal_node is the internal struct that we will use to store the nodes
  //It contains simply a Point object to store the positions we will be using.
  struct internal_node {
    Point node_pt;
    node_value_type val;
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

  //The nested maps to make edge search faster. Acts as a sort of adjacency
  //matrix that is unordered.
  std::unordered_map<size_type, hash_map> edge_search;
};

#endif
