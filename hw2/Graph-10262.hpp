#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

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

template < typename V, typename E >
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the value stored into the nodes. */
  using node_value_type = V ;
  using edge_value_type = E ;

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

  /** Class containing all the relevant informations about a node.
   * It will be stored in a node_information_map. */
  struct NodeInformations;
  using node_information = NodeInformations;

  /** Class containing all the relevant informations abou a edge.
   * It will be stored in a edge_information_map. */
  struct EdgeInformations;
  using edge_information = EdgeInformations;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    this->node_information_map = std::unordered_map< size_type, node_information* >();

    this->edge_map = std::unordered_map< size_type, std::unordered_map< size_type, Edge* > >();
    this->edge_information_map = std::unordered_map< size_type, edge_information* >();
  }

  //* The destructor has to delete all the heap objects that were created. This is done in clear. */
  ~Graph(){
    this->clear();
  }

  /** This struct will hold the informations for each node.
   * It contains a pointer to a Node object, a pointer to a point object and a pointer to a value object.
   * This is very efficient in memory, since those objects are created once and for all in the heap.
   * All copies will be very fast since NodeInformations only contains three pointers.
   * Moreover, the pointer to the node will be shared by the edges, leading to only one object in the heap. */
  struct NodeInformations{
    Node* node;
    Point* point;
    node_value_type* value;
  
    NodeInformations(){
      this->node = nullptr;
      this->point = nullptr;
      this->value = nullptr;
    }

    NodeInformations(Node* node, Point point, node_value_type value){
      this->node = node;
      this->point = new Point(point.x, point.y, point.z);
      this->value = new node_value_type(value);
    }

    ~NodeInformations(){
      delete this->node;
      delete this->point;
      delete this->value;
    }
  };

  /** This struc will hold the informations about the edges.
   * It will create a pointer to an Edge object and a pointer to a value.
   * The Edge* will be shared by the adjacency map, leading to only one Edge in memory. */ 
  struct EdgeInformations{
    Edge* edge;
    edge_value_type* value;

    EdgeInformations(){
      this->edge = nullptr;
      this->value = nullptr;
    }

    EdgeInformations(Edge* edge){
      this->edge = edge;
      this->value = new edge_value_type();
    }
    ~EdgeInformations(){
      delete this->edge;
      delete this->value;
    }
  };

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   * Nodes implement operators == and <, and get totally ordered through the totally ordered template
   */
  class Node : private totally_ordered<Node>{
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
      this->graph = nullptr;
      this->node_index = -1;
    }

    /** Return this node's position. */
    Point& position() const {
      return *(this->graph->node_information_map.at(this->node_index)->point);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->node_index;
    }

    /** Return a modifiable reference to a Node's value. */
    node_value_type & value (){
      return *(this->graph->node_information_map.at(this->node_index)->value);
    }

    /** Return a const reference to a Node's value. This cannot modify a node's value. */
    const node_value_type & value () const{
      return *(this->graph->node_information_map.at(this->node_index)->value);
    }

    /** Return the degree of a node. O(1) operations. */
    size_type degree() const{
      return this->graph->edge_map[this->node_index].size();
    }

    /** Return an object able to iterate over the neighbours of n, pointing on the first neighbor. */   
    incident_iterator edge_begin() const{
      incident_iterator iterator;
      iterator.it= this->graph->edge_map[this->node_index].begin();
      return iterator;
    }

    /** Return an iterator over the end of the neighbors of n. */
    incident_iterator edge_end() const{
      incident_iterator iterator;
      iterator.it= this->graph->edge_map[this->node_index].end();
      return iterator;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->graph == n.graph && this->node_index == n.node_index;
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
      bool inferior;
      // if both Nodes belong to the same graph, check the order of their index
      // else, cast pointer graph to an integer, and compare the adresses as integers
      if(n.graph == this->graph){
        inferior = this->node_index < n.node_index;
      }
      else{
        auto i = reinterpret_cast<std::uintptr_t>(this->graph);
        auto j = reinterpret_cast<std::uintptr_t>(n.graph);

        inferior = i < j;
      }
      return inferior;
    }

   private:
    /** This is a private consructor of a valid node.
     * Only the graph can access it, so the user cannot create a valid node from outside. */
    Node(Graph* graph, size_type node_index):
    graph(graph), node_index(node_index){}

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type node_index;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->node_information_map.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value The new node's value
   * 
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == node_value
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type & node_value = node_value_type ()) {
    Node* new_node = new Node(this, this->size());

    this->node_information_map[this->size()] = new node_information(new_node, position, node_value);

    return *new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    bool result = (this->node_information_map.count(n.index()) > 0);
    return result;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    Node* n_out;
    try{
      n_out = this->node_information_map.at(i)->node;
    }
    catch(const std::out_of_range e){
      std::cout << "Index not in graph" << std::endl;
      return Node();
    }
    return *n_out;
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   * Nodes implement operators == and <, and get totally ordered through the totally ordered template
   */
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      node_1 = nullptr;
      node_2 = nullptr;
      index = -1;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return *(this->node_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return *(this->node_2);
    }

    /** The edge does not hold any reference to the graph, but its node do. */
    edge_value_type & value () {
      return *(node_1->graph->edge_information_map.at(this->index)->value);
    }

    /** Operator at is overloaded to be const in an unordered_map. */ 
    const edge_value_type & value () const{
      return *(node_1->graph->edge_information_map.at(this->index)->value);
    }

    /** Compute the length of an edge. */
    double length() const{
      return norm(this->node_1->position() - this->node_2->position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (*(this->node_1) == *(e.node_1) && *(this->node_2) == *(e.node_2)) || (*(this->node_1) == *(e.node_2) && *(this->node_2) == *(e.node_1));
    }

    bool operator<(const Edge& e) const {
      bool inferior;
      if (*(this->node_1) < *(this->node_2)){
        if (*(e.node_1) < *(e.node_2)){
          inferior = (*(this->node_1) < *(e.node_1)) || (*(this->node_1) == *(e.node_1) && *(this->node_2) < *(e.node_2));
        }
        else{
          inferior = (*(this->node_1) < *(e.node_2)) || (*(this->node_1) == *(e.node_2) && *(this->node_2) < *(e.node_1));
        }
      }
      else{
        if (*(e.node_1) < *(e.node_2)){
          inferior = (*(this->node_2) < *(e.node_1)) || (*(this->node_2) == *(e.node_1) && *(this->node_1) < *(e.node_2));
        }
        else{
          inferior = (*(this->node_2) < *(e.node_2)) || (*(this->node_2) == *(e.node_2) && *(this->node_1) < *(e.node_1));
        }
      }
      return inferior;
    }

   private:

   /** Private constructor for a valid edge. Only through the graph can a valid edge be built. */
    Edge(Node* node_1, Node* node_2, size_type index):
    node_1(node_1), node_2(node_2), index(index){}
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Node* node_1;
    Node* node_2;
    size_type index;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edge_information_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    Edge* edg_out;
    try{
      edg_out = this->edge_information_map.at(i)->edge;
    }
    catch(const std::out_of_range e){
      throw "Edge not in graph";
    }
    return *edg_out;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    bool result;
    try{
      result = (this->edge_map.at(a.index()).count(b.index()) > 0);
    }
    catch(const std::out_of_range e){
      result = false;
    }
    return result;
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
    if (this->has_edge(a, b)){
      return *(this->edge_map.at(a.index()).at(b.index()));
    }
    else{
    // This Edge* will be stored in the edge_information struct, as well as in the adjacency.
    Edge* new_edge = new Edge(this->node_information_map[a.index()]->node,
                              this->node_information_map[b.index()]->node,
                              this->num_edges());

    this->edge_map[a.index()][b.index()] = new_edge;
    this->edge_map[b.index()][a.index()] = new_edge;

    this->edge_information_map[this->num_edges()] = new edge_information(new_edge);

    return *new_edge;
    }

    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * 
   * We have to iterate through the map to delete all the pointers before clearing the map.
   */
  void clear() {
    for (auto it = this->node_information_map.begin(); it != this->node_information_map.end(); ++it){
      delete (it->second);
    }
    for (auto it = this->edge_information_map.begin(); it != this->edge_information_map.end(); ++it){
      delete (it->second);
    }
    this->edge_information_map.clear();
    this->edge_map.clear();
    this->node_information_map.clear();

  }

  //
  // Node Iterator
  //

  // Define the type of an iterator through a map
  // typedef typename std::unordered_map<size_type, node_information*>::const_iterator node_map_iterator;
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   * This will basically be equivalent to an iterator within graph.node_map.
   * This takes the index as argument, which allows an 'ordered' iteration through an unordered map.
   *  
   * Iterating through a unordered_map is basically O(n), with each operation being constant in time
   * */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

      /** The iterator is an index, we can look in the node_information_map at index it. */
      Node operator *() const{
        return *(this->graph->node_information_map.at(it)->node);
      }
      /** Move the iterator forward. We go to the next index */
      node_iterator & operator ++(){
        this->it += 1;
        return *this;
      }
      /** This implicitely defines operator !=. */
      bool operator ==( const node_iterator & rhs) const{
        return this->it == rhs.it;
      }


    private:
    /** Private constructor, only available in the graph. */
    NodeIterator(const Graph* graph, size_type it):
    it(it), graph(graph){}
    friend class Graph;
    // This is an iterator over a map similar to graph.node_map.
    size_type it;
    const Graph* graph;
  };

  /** Return an iterator over the nodes of the graph, pointing over the `first' element of graph.node_map. */
  node_iterator node_begin() const{
    NodeIterator iterator(this, 0);
    return iterator;
  }

  /** Return an iterator pointing to the end of graph.node_map. */ 
  node_iterator node_end() const{
    NodeIterator iterator(this, this->num_nodes());
    return iterator;
  }

  //
  // Incident Iterator
  //
  // Define the type of an iterator through a map
  typedef typename std::unordered_map<size_type, Edge*>::const_iterator edge_map_iterator;
  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. 
   * This will basically be equivalent to an iterator within graph.edge_map[n].*/
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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
      /** Returns the edge to which the iterator refers.
       * Makes sure that node_1 is the origin node.
       * this->it->second is the Edge to which it points
       * this->it->first is the index of the destination node */
      Edge operator *() const{
        if (this->it->second->node_1->index() == this->it->first){
          Node* temp = this->it->second->node_1;
          this->it->second->node_1 = this->it->second->node_2;
          this->it->second->node_2 = temp;
        } 
        return *(this->it->second);
      }

      /** Increments the underlying pointer in the unordered_map. */
      IncidentIterator & operator ++(){
        (this->it)++;
        return *this;
      }

      /** This implicitely defines operator !=. */
      bool operator ==( const IncidentIterator & rhs) const{
        return this->it == rhs.it;
      }

   private:
    friend class Graph;
    edge_map_iterator it;
  };

  //
  // Edge Iterator
  //
  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * This will basically be equivalent to an iterator within graph.edge_indexation*/
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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
    /** Returns the edge to which the iterator refers.
     * Does not care which node is the origin */
    Edge operator *() const{
        return *(this->graph->edge_information_map.at(it)->edge);
      }

      /** Increments the index count */
      EdgeIterator & operator ++(){
        this->it += 1;
        return *this;
      }
      
      /** This implicitely defines operator !=. */
      bool operator ==( const EdgeIterator & rhs) const{
        return this->it == rhs.it;
      }

   private:
   /** Private constructor, only available through graph class. */
    EdgeIterator(const Graph* graph, size_type it):
    it(it), graph(graph){}
    friend class Graph;
    size_type it;
    const Graph* graph;

  };
     /** Return an iterator over the nodes of the graph, pointing over the edge with index 0 */
    edge_iterator edge_begin() const{
    edge_iterator iterator(this, 0);
    // iterator.it= this->edge_information_map.begin();
    return iterator;
  }

    /** Return an iterator pointing to the last edge */ 
    edge_iterator edge_end() const{
    edge_iterator iterator(this, this->num_edges());
    // iterator.it= this->edge_information_map.end();
    return iterator;
  }

  /** Delete the node n and all the relative edges to n.
   * Ensures that the invariants are respected, i.e that there is no "hole" in the indexation.
   * After this, the node that previously held the last index holds the index of the node to delete.
   * The index of the node to delete now refers to the last index, and this is true for the adjacency maps as well.
   * This ensures that there is no "hole" in the indexation.
   * There is no invariant forcing us to keep track of the "initial index" of each node, so this solution works.
   * 
   * This is done as follows :
   *    swap the node to be deleted by the node with the last index
   *          swap the node informations
   *          swap the adjacency maps
   *          make sure that the edges and the nodes are correctly indexed
   *    delete all edges relative to the node to be deleted
   *    delete the node and its informations
   * 
   * This is done very efficiently, since:
   *    The node informations are stored through a pointer, requiring only swaping two pointers in a map
   *    All edges pointing to a node hold a pointer to this node, which means that changing the index of a node
   *    automatically changes the index in all the edges refering to it.
   *    No memory reallocation is necessary, since we are not storing within a vector but within a map
   * 
   * Thus, we only perform :
   *    A swap of pointer ( O(1) )
   *    An iteration on the two adjacency maps ( O( degree ), much smaller than O( nodes ))
   *    An iteration on the adjacency to delete the vectors (same as above), sine deleting an edge is O(1)
   * 
   * Therefore, it is not less efficient than leaving holes in the indexing and use a mapping vector. */ 
  size_type remove_node ( const Node & n ){

    size_type old_index = n.index();
    size_type old_size = this->num_nodes() - 1;

    if ( old_index != old_size){ 

    // swap adjacencies maps
    this->edge_map[old_index].swap(this->edge_map[old_size]);

    // Swap nodes, values and points
    node_information* temp_node_info = this->node_information_map[old_size];
    this->node_information_map[old_size] = this->node_information_map[old_index];
    this->node_information_map[old_index] = temp_node_info;

    // Swap the indexes of the nodes. This is automatically propagated within the edges.
    this->node_information_map[old_index]->node->node_index = old_index;
    this->node_information_map[old_size]->node->node_index = old_size;

    // Now, old_size refers to the index we want to delete and old_index is the former last vector.
    
    // iterate over the adjacency maps to modify the keys of the map
    for (auto it = this->node_information_map[old_index]->node->edge_begin(); it != this->node_information_map[old_index]->node->edge_end(); ++it){
      
      // For each edge adjacent to the node to delete, go in the other leg of the map to modify the indexation in the keys
      //      If the destination of the edge is also connected to the former last node of the graph, swap the edges
      //      Else, just change the keys. 
      Edge* temp_edge = this->edge_map[(*it).node_2->node_index][old_size];

        if (this->edge_map[(*it).node_2->node_index].find(old_index) != this->edge_map[(*it).node_2->node_index].end()){
        this->edge_map[(*it).node_2->node_index][old_size] = this->edge_map[(*it).node_2->node_index][old_index];
        this->edge_map[(*it).node_2->node_index][old_index] = temp_edge;
        }

        else{
          this->edge_map[(*it).node_2->node_index][old_index] = temp_edge;
          this->edge_map[(*it).node_2->node_index].erase(old_size);
        }
    }

    for (auto it = this->node_information_map[old_size]->node->edge_begin(); it != this->node_information_map[old_size]->node->edge_end(); ++it){
      // Do the same thing with the map of the former last index 
      Edge* temp_edge = this->edge_map[(*it).node_2->node_index][old_index];
        if (!(this->edge_map[(*it).node_2->node_index].find(old_size) != this->edge_map[(*it).node_2->node_index].end())){
          this->edge_map[(*it).node_2->node_index][old_size] = temp_edge;
          this->edge_map[(*it).node_2->node_index].erase(old_index);
        }
    }
    }

    // Delete everything

    // Iterate over the adjacent edges to the node to delete, and delete the edges.
    for (auto it = this->node_information_map[old_size]->node->edge_begin(); it != this->node_information_map[old_size]->node->edge_end();){
      auto it2 = std::next(it);
      this->remove_edge((*it));
      it = it2;
    }

    // Delete the adjacency map
    this->edge_map.erase(old_size);

    // Delete the node and its informations
    delete (this->node_information_map[old_size]);
    this->node_information_map.erase(old_size);

    return 1;
  }

  /** Deletes the node refered to by n_it.
   * It calls function *n_it, which swaps objects and then deletes the object.
   * Since this is done by swapping, the iterator is still valid.
   * If we remove the last one, n_it is now equal to graph.node_end(). */
  node_iterator remove_node ( node_iterator n_it ){
    this->remove_node((*n_it));
    return n_it;
  }

  /** Deletes an edge, refered by its two nodes.
   * Calls the remove_edge function over a single edge, defined below.
   * If the edge does not exist, returns 0. */
  size_type remove_edge ( const Node & a, const Node &b){
    if (this->has_edge(a, b)){
    Edge edge_to_erase = this->add_edge(a, b);
    this->remove_edge(edge_to_erase);
    return 1;
    }
    else{
      return 0;
    }
  }

  /** Deletes the edge e. In order to keep the edge map 'dense', we swap with
   * the last index, then remoe the last from the index list.
   * 
   * We also make sure to clean the adjacency map.
   * Since this is being done by swapping, all edge iterators remain valid. */
  size_type remove_edge ( const Edge & e){

    // Delete the two adjacencies of e
    this->edge_map[e.node_1->node_index].erase(e.node_2->node_index);
    this->edge_map[e.node_2->node_index].erase(e.node_1->node_index);


    size_type old_size = this->num_edges() - 1;
    size_type old_index = e.index;

    // Swap the edge to delete with the last one in the edge list.
    edge_information* temp_edge_information = this->edge_information_map[old_size];
    this->edge_information_map[old_size] = this->edge_information_map[old_index];
    this->edge_information_map[old_index] = temp_edge_information;
    this->edge_information_map[old_index]->edge->index = old_index;
    this->edge_information_map[old_size]->edge->index = old_size;
    // The Edge* stored in the adjacency maps is the same one as in the edge_information,
    // so this change in index is also propagated to the adjacency map

    delete (this->edge_information_map[old_size]);
    this->edge_information_map.erase(old_size);

    return 1;

  }

  /** Delete edge refered to by iterator e_it. Since this is done by swapping, e_it remains valid
   * at the end, pointing to the ormer last element.
   * 
   * If we remove the last element, the iterator is now equal to graph.edge_end(). */
  edge_iterator remove_edge ( edge_iterator e_it ){
    this->remove_edge((*e_it));
    return e_it;
  }

 private:
  std::unordered_map< size_type, node_information* > node_information_map;

  std::unordered_map<size_type, std::unordered_map<size_type, Edge*>> edge_map;
  std::unordered_map< size_type, edge_information* > edge_information_map;

};

#endif // CME212_GRAPH_HPP
