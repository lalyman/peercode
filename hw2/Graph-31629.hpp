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

template <typename V = int,typename E = int>
class Graph {

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

   /** Node values type */
   using node_value_type = V;

   /** Edge values type */
   using edge_value_type = E;


 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  std::unordered_map<size_type,std::pair<Point,node_value_type>>* nodes_;
  std::unordered_map<size_type,std::pair<std::pair<size_type, size_type>,edge_value_type>>* edges_list_;
  std::unordered_map<size_type,std::unordered_map<size_type, size_type>*>* rev_edges_list_;
public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
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



  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    /**
    * General remarks on the datastructure:
    * We use unordered maps to store the values, as it is more scalable than a vector.
    * For the sake of simplicity, we store the edges as an unordered map of
    * unordered maps, corresponding to an adjacency matrix, which allows an easy
    * iteration. However, we also keep a list of all edges as pairs of indexes.
    * This allows to quickly find the edge having a given index, and to
    * quickly iterate over the edges. Since we do not have space constraints in the subject,
    * I have chosen to keep this implementation until there are additional constraints.
    */

    // Nodes
    nodes_ = new std::unordered_map<size_type,std::pair<Point,node_value_type>>();
    // Edges
    edges_list_ = new std::unordered_map<size_type,std::pair<std::pair<size_type, size_type>,edge_value_type>>();
    rev_edges_list_ = new std::unordered_map<size_type,std::unordered_map<size_type, size_type>*>;
  }

  /** Destructor */
  ~Graph() {
    //Clear the graph
   clear();
   //deallocate the member pointers
   delete nodes_;
   delete edges_list_;
   delete rev_edges_list_;
 };

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
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
      graph_ = nullptr;
      uid_ = size_type (-1);
    }

    /** Return this node's position. */
    const Point& position() const {
      const Point& point = (*(graph_->nodes_))[uid_].first;
      return point;
    }

    /** Return this node's position, non const method */
    Point& position() {
      Point& point = (*(graph_->nodes_))[uid_].first;
      return point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return size_type(uid_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /* * Return the value associated with a node
     * @return value associated with the node
     * @pre this is a valid node
     *
     * Non-constant version, can be used to change the value
     */

    node_value_type & value(){
      node_value_type& node_value = (*(graph_->nodes_))[uid_].second;
      return node_value;
    }

    /* * Return the value associated with a node
     * @return value associated with the node
     * @pre this is a valid node
     * @post the graph state cannot be changed
     * Constant version, cannot be used to change the value
     */
    const node_value_type & value() const{
      const node_value_type& node_value =  (*(graph_->nodes_))[uid_].second;
      return node_value;
    }

    /* * Return the degree of a node
     * @return number of neighbors of the node in the graph
     * @pre this is a valid node
     * @pre the rev_edges_list_ map has been properly initialised, i.e each node
     *      index has at least an empty map
     *
     */

    size_type degree() const{
      assert(graph_->rev_edges_list_->count(uid_)>0);
      return (*(graph_->rev_edges_list_))[uid_]->size();
    }

    /* * Return an iterator of neighbors
     * @return iterator of type incident_iterator, iterates over the neighbors of
     *         the node
     * @pre this is a valid node
     * @pre the rev_edges_list_ map has been properly initialised, i.e each node
     *      index has at least an empty map
     * @post the iterator will iterate exactly over the neighbors
     *
     * To iterate over the neighbors, we simply initialize an iterator using the
     * built- in iterators of the standard library over the unordered_maps
     */
    incident_iterator edge_begin() const{
     assert (graph_->rev_edges_list_->count(uid_)>0);
     std::unordered_map<size_type,size_type>::iterator it = (*(graph_->rev_edges_list_))[uid_]->begin();
     incident_iterator ni  = incident_iterator(graph_,it,uid_);
     return ni;
    }

    /* * Return an iterator of neighbors
     * @return iterator correponding to the end of the neighbors
     * @pre this is a valid node
     * @pre the rev_edges_list_ map has been properly initialised, i.e each node
     *      index has at least an empty map
     * @post the iterator is equal to the end of iteration over neighbors
     *
     * To iterate over the neighbors, we simply initialize an iterator using the
     * built- in iterators of the standard library over the unordered_maps
     */
    incident_iterator edge_end() const{
     assert (graph_->rev_edges_list_->count(uid_)>0);
     std::unordered_map<size_type,size_type>::iterator it = (*(graph_->rev_edges_list_))[uid_]->end();
     incident_iterator ni  = incident_iterator(graph_,it,uid_);
     return ni;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool result =(graph_== n.graph_);
      result = result && (uid_ == n.uid_);
      return result;
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
      return index()< n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    const Graph*  graph_;

    size_type  uid_;

    /* This constructor is used to buil valid nodes*/
    Node(const Graph* graph, size_type& uid) {
      graph_ = graph;
      uid_ = uid;
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_->size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] position The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */

  Node add_node(const Point& position,const node_value_type& node_value = node_value_type()) {
    size_type old_size = size();
    nodes_->insert({{old_size,std::pair<Point,node_value_type>(position,node_value)}});
    Node new_node = Node(this,old_size);
    std::unordered_map<size_type,size_type> * new_map = new std::unordered_map<size_type,size_type> ();
    rev_edges_list_->insert({{old_size,new_map}});

    return new_node;
  }

  /**
   * Removes a node from the graph
   * @param  n node to remove
   * @return   0 if a no node has been removed, 1 if a node has been removed
   *
   * @post If has_node(@a n), then after the function is applied, the node is
   *       removed, i.e
   *             - new num_nodes() = old num_nodes()
   *             - new num_edges() = old num_edges() - @a n.degree()
   *             - The nodes in the new graph have exactly the same position as
   *             the nodes in the old graph except the one that is removed
   * @post The graph is still valid
   * @post if old num_nodes()<1 and @a n is a valid node then the graph is cleared
   * @post if old num_nodes()>1 and @a n is a valid node with @a n.index()<num_nodes()-1
   *       let j= @a n.index() and i =num_nodes()-1
   *       Then old node(i) has index j in the new graph
   *       All other nodes have the same index.
   *       - Edges, Nodes objects that do not involve old node(j) or old node(i)
   *       are still valid
   *       - IncidentIterators for nodes that are not old node(j) or old node(i)
   *       and do not have them as neighbors are still valid
   *       - NodeIterators referring to a node of index < min(i,j) are still valid
   *       - EdgeIterators are invalid
   * @post if old num_nodes()>0 and @a n is a valid node with @a n.index()==num_nodes()-1
   *       All remaining nodes have the same index.
   *       - Edges, Nodes objects that do not involve old node(i) are still valid
   *       - IncidentIterators for nodes that are not old node(i) and do not
   *        have it as neighbor are still valid
   *       - NodeIterators referring to a node of index < i are still valid
   *       - EdgeIterators are invalid
   *
   * Given that complexity of all maps operations is O(1), that complexity of
   * remove_edge is O(1), the complexity of this function is
   *  O(@a n.degree() + n1.degree()), which is roughly O(1) if the graph is sparse,
   *  O(num_nodes()) if the graph is not sparse (Worst case)
   *
   */
  size_type remove_node(const Node& n){
    if (has_node(n)){
      if(num_nodes()==1){
        clear();
        return 1;
      }else{
        size_type node_id {n.uid_};
        if (node_id==num_nodes()-1){
         //Case 1: The node to remove has the last index
         //We generate an adjacency list
          std::vector<size_type> neighbors_index;
          for(incident_iterator nii = n.edge_begin();nii!=n.edge_end();++nii){
            neighbors_index.push_back((*nii).node2().index());
          }
          //We remove all incident edges
          for(auto neighbor_index = neighbors_index.begin();neighbor_index != neighbors_index.end();++neighbor_index){
            remove_edge(node(node_id), node(*neighbor_index));
          }
          //We delete all entries corresponding to this node
          delete (*rev_edges_list_)[node_id];
          rev_edges_list_->erase(node_id);
          nodes_->erase(node_id);
        }else{
         //Case 2 : We will switch the node with the node with last index to erase it and keep index invariant

        //set ids
        size_type node_to_remove_id = num_nodes()-1;
        //buffer for neighbors - prevents from invalidating iterators
        std::vector<size_type> neighbors_index;
        //computation of the adj list of the node to remove
        for(incident_iterator nii = n.edge_begin();nii!=n.edge_end();++nii){
          neighbors_index.push_back((*nii).node2().index());
        }

        //removal of all edges
        for(auto neighbor_index = neighbors_index.begin();neighbor_index != neighbors_index.end();++neighbor_index){
          remove_edge(node(node_id), node(*neighbor_index));
        }
        //clear the array
        neighbors_index.clear();
        Node node_to_remove =node(node_to_remove_id);
        // we compute the adj list of the node to swap
        for(incident_iterator nii = node_to_remove.edge_begin();nii!=node_to_remove.edge_end();++nii){
          neighbors_index.push_back((*nii).node2().index());
        }
        //Then we reindex all the edges of the node to move with the new index
        for(auto neighbor_index = neighbors_index.begin();neighbor_index != neighbors_index.end();++neighbor_index){

         size_type edge_to_move_id {(*(*rev_edges_list_)[*neighbor_index])[node_to_remove_id]};
          //update adj list
          (*rev_edges_list_)[*neighbor_index]->insert({{node_id,edge_to_move_id}});
          (*rev_edges_list_)[node_id]->insert({{*neighbor_index,edge_to_move_id}});
          size_type min_idx,max_idx;
          min_idx = (*neighbor_index<node_id)?(*neighbor_index):node_id;
          max_idx =  (*neighbor_index>=node_id)?(*neighbor_index):node_id;
          //update edges list
          (*edges_list_)[edge_to_move_id].first=std::pair<size_type,size_type>(min_idx,max_idx);
          //erase old entries
          (*rev_edges_list_)[*neighbor_index]->erase(node_to_remove_id);
          (*rev_edges_list_)[node_to_remove_id]->erase(*neighbor_index);
        }
        // Now we can swap the actual values of the nodes
        std::swap((*nodes_)[node_id], (*nodes_)[node_to_remove_id]);
        //We now delete the remaining entries for the old index
        delete (*rev_edges_list_)[node_to_remove_id];
        rev_edges_list_->erase(node_to_remove_id);
        nodes_->erase(node_to_remove_id);
      }
        return 1;
    }
  }else{
    return 0;
  }
  }


  /**Iterator version of previous method, returns iterator to next valid node
  *
  * @pre @a n_it!= node_end()
  * @post The return value iterates exactly over the nodes that @a n_it has not
  *       iterated on
  */
  node_iterator remove_node( node_iterator n_it ){
    size_type idx {(*n_it).index()};
    remove_node(*n_it);
    return NodeIterator(this,idx);
  }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_==this)&&(n.uid_ < this->size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    const Node node = Node(this,i);
    return node;        // Invalid node
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
      uid_ = size_type(-1);
      id1_ = size_type(-1);
      id2_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(id1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(id2_)  ;
    }

    /** Return the length of this Edge (Constant complexity) */
    double length() const{
      return norm(node1().position()-node2().position());
    }

    /** Return the value of this Edge (Constant complexity) */
    const edge_value_type & value() const{
      const edge_value_type& edge_value = (*(graph_->edges_list_))[uid_].second;
      return edge_value;
    }

    /** Return the value of this Edge, no const version (Constant complexity) */
    edge_value_type & value(){
      edge_value_type& edge_value = (*(graph_->edges_list_))[uid_].second;
      return edge_value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_==e.graph_ && ((id1_==e.id1_&&id2_==e.id2_)||(id2_==e.id1_&&id1_==e.id2_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (e.graph_<graph_)||((e.graph_==graph_)&&(uid_<e.uid_));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    const Graph * graph_;
    size_type uid_;
    size_type id1_;
    size_type id2_;

    Edge(const Graph* graph, const size_type& uid, const size_type& id1, const size_type& id2) {
      graph_ = graph;
      uid_ = uid;
      id1_ = id1;
      id2_ = id2;
      // HW0: YOUR CODE HERE
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Since we have used an extra map for the edges, this is O(1) complexity
   */
  size_type num_edges() const {
    return edges_list_->size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   * @post result_edge.index() == @a i
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Since we have used an extra map for the edges, this is O(1) complexity
   */
  Edge edge(size_type i) const{
    std::pair<size_type, size_type> edge_pair = (*edges_list_)[i].first;
    const Edge result_edge = Edge(this,i,edge_pair.first,edge_pair.second);
    return result_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * We simply have to perform a lookup in the adjacency table, using the nodes
   * indexes
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type min_idx,max_idx;
    min_idx = (a.uid_<b.uid_)?a.uid_:b.uid_;
    max_idx = (a.uid_>=b.uid_)?a.uid_:b.uid_;
    return (rev_edges_list_->count(min_idx)>0)&&((*rev_edges_list_)[min_idx]->count(max_idx) >0);


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
   * Complexity here is O(1) due to the storage of 2 datastructures
   */
  Edge add_edge(const Node& a, const Node& b,const edge_value_type& edge_value = edge_value_type()) {
    size_type min_idx,max_idx;
    min_idx = (a.uid_<b.uid_)?a.uid_:b.uid_;
    max_idx = (a.uid_>=b.uid_)?a.uid_:b.uid_;
    if (has_edge(a,b)){
     // If we actually have an edge, we just return it
     // To do so, we can use the edge method, that has a constant cost due to
     // the storage of an additional map of pairs. We can retrieve the index of
     // the existing edge in constant time
      return edge((*(*rev_edges_list_)[min_idx])[max_idx]);
    }
    else{
     //Compute the new index
      const size_type old_sz = num_edges();
      const std::pair<size_type, size_type> pair1  = std::pair<size_type, size_type>(min_idx,max_idx);

      //Insert in pairs map
      edges_list_->insert({{old_sz, std::pair<std::pair<size_type, size_type>,edge_value_type>(pair1,edge_value)}});
      //We check that the adjacency lists are properly built
      assert (rev_edges_list_->count(min_idx));
      assert (rev_edges_list_->count(max_idx));

      assert (!((*rev_edges_list_)[min_idx]->count(max_idx)));
      assert (!((*rev_edges_list_)[max_idx]->count(min_idx)));

      //We add the new edge
      (*rev_edges_list_)[min_idx]->insert({{max_idx,old_sz}});
      (*rev_edges_list_)[max_idx]->insert({{min_idx,old_sz}});
      //this method has O(1) cost, hence we can use it to return the result
     return edge(old_sz);
  }
}

/**
 * Removes an edge from the graph
 * @param  edge_index_todelete index of the edge to delete
 * @return index of the deleted edge
 *
 * @post The graph is still valid
 * @post old num_edges()-- = new num_edges()
 * @post The edges of the new graph are exactly the same as the edges of the old graph,
 *       except for the one that is removed
 * @post The edges index after removal are the elements of [0, new num_edges())
 * @post if old num_edges()>1 and @a edge_index_todelete<num_nodes()-1
 *       let j= @a edge_index_todelete and i =num_edges()-1
 *       Then old edge(i) has index j in the new graph
 *       All other edges have the same index.
 *       - Edge objects of indexes different from j,i are still valid
 *       - IncidentIterators for nodes that are not  old edge(j).node1(), old edge(i).node1(),
 *       old edge(j).node2() or old edge(i).node2() are still valid
 *       - EdgeIterators referring to an edge of index < min(i,j) are still valid
 *       - NodeIterators are still valid
 * @post if old num_edges()>0 and @a edge_index_todelete ==num_nodes()-1
 *       let  i = @a edge_index_todelete
 *       All remaining nodes have the same index.
 *       - Edge objects of indexes different from i are still valid
 *       - IncidentIterators for nodes that are not edge(i).node1() or
 *       edge(i).node2() are still valid
 *       - EdgeIterators referring to an edge of index < i are still valid
 *       - NodeIterators are still valid
 *
 * Complexity is O(1), since all the lookups in the map are O(1) on average,
 * erasing operations for unordered map are O(1) on average, and swapping values
 * using std::swap is O(1)
 *
 */
  size_type remove_edge(const size_type edge_index_todelete){
    if (num_edges()==1){
      (*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.first]->erase((*edges_list_)[edge_index_todelete].first.second);
      (*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.second]->erase((*edges_list_)[edge_index_todelete].first.first);
      edges_list_->clear();
    }else{
      if (num_edges()-1==edge_index_todelete){
        (*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.first]->erase((*edges_list_)[edge_index_todelete].first.second);
        (*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.second]->erase((*edges_list_)[edge_index_todelete].first.first);
        edges_list_->erase(edge_index_todelete);
      }else{
    size_type edge_index_toswap  = num_edges()-1;
    //swap in edges list
    std::swap((*edges_list_)[edge_index_toswap], (*edges_list_)[edge_index_todelete]);
    //delete old entry
    (*rev_edges_list_)[(*edges_list_)[edge_index_toswap].first.first]->erase((*edges_list_)[edge_index_toswap].first.second);
    (*rev_edges_list_)[(*edges_list_)[edge_index_toswap].first.second]->erase((*edges_list_)[edge_index_toswap].first.first);

    //update index for swapped entry
    (*(*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.first])[(*edges_list_)[edge_index_todelete].first.second] = edge_index_todelete;
    (*(*rev_edges_list_)[(*edges_list_)[edge_index_todelete].first.second])[(*edges_list_)[edge_index_todelete].first.first] = edge_index_todelete;
    //delete entries in edges list and value
    edges_list_->erase(edge_index_toswap);

  }
}
return edge_index_todelete;
}

/**
 * Remove edge (pair of nodes version)
 * @param  n1 a node of the graph
 * @param  n2 a node of the graph
 * @return    0 if !(has_edge(@a n1,@a n2)), 1 otherwise
 *
 * @post if !(has_edge(@a n1,@a n2)), then the graph is not changed
 *       if has_edge(@a n1,@a n2), then the edge is removed (see docuemtation for
 *       index based implementation for more details on specification)
 *
 * Complexity : O(1)
 */
  size_type remove_edge( const Node& n1, const Node& n2){
    if(has_edge(n1,n2)){
      remove_edge((*(*rev_edges_list_)[n1.uid_])[n2.uid_]);
      return 1;
    }else{
      return 0;
    };
  }

  /**
   * Remove edge (edge object version)
   * @param  edge_todelete an edge of the graph
   * @return    1
   *
   * @pre edge_todelete is a valid edge of the graph
   * @post The edge is removed
   *
   * (see documentation for index based implementation for more details on specification)
   * Complexity : O(1)
   */
  size_type remove_edge(const Edge& edge_todelete){
    remove_edge(edge_todelete.uid_);
    return 1;
  }

  /**
   * Remove edge (edge iterator version)
   * @param  e_it an edge iterator
   * @return  an edge iterator to the next valid edge
   *
   * @pre e_it!= old edge_end()
   * @post The edge *e_it is removed
   * @post the return value is a valid iterator, that iterates over all the
   *       edges of the graph that e_it has not itered on
   *
   * (see documentation for index based implementation for more details on specification)
   * Complexity : O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it){
      size_type index_new_edge{ (*e_it).uid_};
      remove_edge(index_new_edge);
      return EdgeIterator(this, index_new_edge);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   *
   * Complexity : O(num_nodes()+ num_edges()) since we have to destroy the maps of
   * edges and nodes, i.e. O(num_nodes()) if the graph is sparse
   */
  void clear() {
    //destroy items
    nodes_->clear();
    edges_list_->clear();
    for(auto mi = rev_edges_list_->begin();mi!= rev_edges_list_->end();++mi){
     mi->second->clear();
     delete(mi->second);
    }
    rev_edges_list_->clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   *
   * To implement the NodeIterator class, we simply have to store a pointer to
   * the graph and an index. Actually, since the nodes are identified by an
   * index and that this index is in the range [0,num_nodes), this is very
   * easy to iterate over the nodes with all operations in constant time
   */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_ = nullptr;
      index_node = size_type(-1);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /* * Access to the item
    * @return the current iteration item
    * We simply use the current index
    */
    Node operator*() const {
      return graph_->node(index_node);
    }

    /* * Increment the iterator
    * @return the updated iterator
    * Given that we use indexes to iterate, we only have to increment the index
    */
    node_iterator& operator++(){
      index_node++;
      return *this;
    }

    /* * Equality
    * @param ni a node iterator to compare
    * @return true if the iterators are identical, false otherwise
    * With our implementation, 2 iterators are equal iff they point towards the
    * same graph and are at the same index.
    */
    bool operator==(const node_iterator &ni) const{
      return (graph_==ni.graph_)&&(index_node==ni.index_node);
    }

   private:
    friend class Graph;
    const Graph *  graph_;
    size_type index_node;

    /* *Valid constructor for Node iterators  */
    NodeIterator(const Graph * graph,size_type index_nod){
      graph_ = graph;
      index_node = index_nod;
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /* * Return an iterator of nodes
   * @return iterator of type node_iterator, iterates over the nodes of the graph
   * @post the iterator will iterate exactly over the nodes in the graph
   *
   * To iterate over all the nodes, we simply initialize an iterator at the index
   * 0
   */
  node_iterator node_begin() const{
    node_iterator ni  = NodeIterator(this,size_type(0));
    return ni;
  }

  /* * Return an iterator to the end of the nodes
   * @return iterator of type node_iterator, end of iteration over the nodes
   *
   * To iterate over all the nodes, we simply initialize an iterator at the index
   * that is the first one not used in the graph, i.e num_nodes()
   */
  node_iterator node_end() const{
    node_iterator ni  = NodeIterator(this,this->num_nodes());
    return ni;
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   *
   * To implement the IncidentIterator class, we use the fact that we have stored
   * our edges as an adjacency matrix encoded as a map of maps. As a result, the
   * iteration over the neighbors of a given node is simply the iteration over
   * the unordered map at this particular index, i.e over a row of the adjacency
   * list
   * We can therefore use built-in iterators, and store the built in iterator
   * ove the neighbors, the id of the node over whose neighbors we iterate, and
   * a pointer to the graph
   */
  class IncidentIterator :private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graph_ = nullptr;
      it = nullptr;
      source_node_id = size_type(-1);
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const


    /* * Access to the item
    * @return the current iteration item
    * Given that the iterator over the unordered map give the edge id and the id
    * of the neighbor, and that our structure stores the id of the source node,
    * we simply have to build the appropriate Edge object and return it
    *
    * Complexity O(1)
    */
    Edge operator*()const{
      Edge edge = Edge(graph_,it->second,source_node_id,it->first);
      return edge;
    }

    /* * Increment the iterator
    * @return the updated iterator
    * We simply increment the underlying iterator
    *
    * Complexity O(1)
    */
    IncidentIterator& operator++(){
      it++;
      return *this;
    }

    /* * Equality
    * @param nii a incident iterator to compare
    * @return true if the iterators are identical, false otherwise
    * To compare the iterators, we simply compare their members : the raph, the source
    * nodes and the undelying unordered_map iterator
    *
    * Complexity O(1)
    */
    bool operator==(const IncidentIterator& nii) const{
      return graph_==nii.graph_ &&source_node_id==nii.source_node_id && it==(nii.it);
    }

   private:
    friend class Graph;
    const Graph * graph_;
    std::unordered_map<size_type,size_type>::iterator it;
    size_type source_node_id;

    /* *Valid constructor for Incident iterators  */
    IncidentIterator(const Graph * graph,std::unordered_map<size_type,size_type>::iterator ite,size_type source_node_ide){
      graph_ = graph;
      it = ite;
      source_node_id = source_node_ide;
    }
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * Given that we have kept a map of edges, and that we have an invariant on the indexes
   * we can proceed the same way we did for NodeIterators, simply by using an
   * index and a pointer to the graph, that allows to keep all operations in O(1)
   */
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
     graph_ = nullptr;
     index_edge = size_type(-1);
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /* * Access to the item
    * @return the current iteration item
    * We simply use the current index
    *
    * Complexity O(1) due to the extra storage, since the edge method is very efficient
    */
    Edge operator*() const {
      return graph_->edge(index_edge);
    }

    /* * Increment the iterator
    * @return the updated iterator
    * Given that we use indexes to iterate, we only have to increment the index
    */
    edge_iterator& operator++(){
      index_edge++;
      return *this;
    }

    /* * Equality
    * @param ei a edge iterator to compare
    * @return true if the iterators are identical, false otherwise
    * With our implementation, 2 iterators are equal iff they point towards the
    * same graph and are at the same index.
    */
    bool operator==(const edge_iterator &ei) const{
      return (graph_==ei.graph_)&&(index_edge==ei.index_edge);
    }
   private:
    friend class Graph;
    const Graph* graph_;
    size_type index_edge;

    EdgeIterator(const Graph * graph,size_type index_edg){
      graph_ = graph;
      index_edge = index_edg;
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /* * Return an iterator of nodes
   * @return iterator of type edge_iterator, iterates over the edges of the graph
   * @post the iterator will iterate exactly over the edges in the graph
   *
   * To iterate over all the edges, we simply initialize an iterator at the index
   * 0
   */
  edge_iterator edge_begin() const{
   edge_iterator ei = EdgeIterator(this,size_type(0));
   return ei;
  }

  /* * Return an iterator to the end of the edges
   * @return iterator of type edge_iterator, end of iteration over the edges
   *
   * To iterate over all the edges, we simply initialize an iterator at the index
   * that is the first one not used in the graph, i.e num_edges()
   */
  edge_iterator edge_end() const{
   edge_iterator ei = EdgeIterator(this,this->num_edges());
   return ei;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.



};

#endif // CME212_GRAPH_HPP
