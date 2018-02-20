  #ifndef CME212_GRAPH_HPP
  #define CME212_GRAPH_HPP

  /** @file Graph.hpp
   * @brief An undirected graph type
   */

  #include <algorithm>
  #include <vector>
  #include <unordered_map>
  #include <cassert>

  #include "CME212/Util.hpp"
  #include "CME212/Point.hpp"


  /** @class Graph
   * @brief A template for 3D undirected graphs.
   *
   * Users can add and retrieve nodes and edges. Edges are unique (there is at
   * most one edge between any pair of distinct nodes).
   */
  template <typename V, typename E>
  class Graph {
   private:

    /** Predeclaration of the internal representation of a node. */
    struct internal_node;

    /** Predeclaration of the internal representation of an edge. */
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
    
    /** Synonym for the template argument */
    using node_value_type = V;

    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;

    /** Synonym for the template argument */
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
    Graph() : nodes_(), nUid_(), edges_(), eUid_(), adj_() {
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
     * Nodes store a position represented by a Point, and a unique value 
     * defined by the Graph template.
     */
    class Node : private totally_ordered<Node> {
     public:
      /** Construct an invalid node.
       * @post Node.graph_ is a null pointer and the node index 
       *       is max(unsigned int).
       *
       * Invalid nodes are represented by a null pointer to the graph
       * and a unique index of -1. Since index is an unsigned int, -1 is 
       * actually the maximum value of the size_type attribute. This should be 
       *  fine, unless we have a graph with the maximum amount of nodes. In that 
       * case, the graph could overflow anyway. 
       * @endcode
       */
      Node() :  graph_(nullptr), uindex_(-1) {}

      /** Return this node's position.
       * @pre The calling Node has a valid unique index i 
       *      with 0 <= i < nodes_.size()
       *
       * @return Point representing the node position
       */
      const Point& position() const {
        // Delegate the query to the Graph class via fetchNode()
        return fetchNode().location;
      }

      /** Return this node's position.
       * @pre The calling Node has a valid unique index i 
       *      with 0 <= i < nodes_.size()
       *
       * @return Modifiable reference to point representing the node position
       */
      Point& position() {
        // Delegate the query to the Graph class via fetchNode()
        return fetchNode().location;
      }

      /** Return this node's index, a number in the range [0, graph_size).
       * @return Index i of the Node, with 0 <= i < graph_size 
       *         or max(unsigned int) if the Node is invalid.
       */
      size_type index() const {
        return fetchNode().idx;
      }

      /** Get a reference to this Node's value.
       * @pre  The calling Node has a valid unique index i 
       *       with 0 <= i < nodes_.size().
       *
       * @return The calling Node's value reference.
       */
      node_value_type& value() {
        return fetchNode().value;
      }

      /** Get a constant reference to this Node's value.
       * @pre  The calling Node has a valid unique index i 
       *       with 0 <= i < nodes_.size().
       *
       * @return The calling Node's value constant reference.
       */
      const node_value_type& value() const {
        return fetchNode().value;
      }

      /** Get this Node's degree.
       * A Node degree is defined as the number of Edges that are incident
       * on that Node.
       * @pre  The Node is valid within a Graph. It has an entry in the
       *       adjacency matrix.
       *
       * @return The calling Node's degree.
       */
      size_type degree() const {
        // Use the Node unique id to find it in the adjacency matrix
        return graph_->adj_[graph_->nUid_[this->index()]].size();
      }

      /** Get an Incident Iterator pointing to the first incident Edge 
       *  for this Node.
       * An Incident Iterator is an iterator for all the Edges that are
       * incident on the calling Node.
       * @pre  The Node is valid within a Graph. It has an entry in the
       *       adjacency matrix.
       * @post The incident iterator will point to the first Edge in the 
       *       Graph internal Edge vector that has this node in one of its 
       *       vertices.
       *
       * @return An incident iterator to traverse the Edges incident on this
       *         node.
       */
      incident_iterator edge_begin() const {

        //  Make sure the node is inside the adjacency matrix
        assert(graph_->adj_.find(uindex_) != graph_->adj_.end());

        std::unordered_map<size_type,size_type>::iterator it = 
            graph_->adj_.find(uindex_)->second.begin();

        return incident_iterator{graph_,it,uindex_};
      }

      /** Get an Incident Iterator pointing outside of the range of incident 
       *  Edge for this Node.
       * An Incident Iterator is an iterator for all the Edges that are
       * incident on the calling Node.
       * @pre  The Node is valid within a Graph. It has an entry in the
       *       adjacency matrix.
       * @post The incident iterator will point one element past the last
       *       Edge in the Graph internal Edge vector that has this node in 
       *       one of its vertices.
       *
       * @return An incident iterator to mark the limit of Edges incident on 
       *         this node.
       */
      incident_iterator edge_end() const {
        //  Make sure the node is inside the adjacency matrix
        assert(graph_->adj_.find(uindex_) != graph_->adj_.end());

        std::unordered_map<size_type,size_type>::iterator it = 
             graph_->adj_.find(uindex_)->second.end();

        return incident_iterator{graph_,it,uindex_};
      }
  
      /** Test whether this node and @a n are equal.
       *
       * Equal nodes have the same graph and the same index.
       * @param[in] n Node to compare against.
       *
       * @return True if this Node is equal to @a n, false otherwise.
       */
      bool operator==(const Node& n) const {
        // Compare graph pointers and indexes.
        return this->graph_ == n.graph_ && this->uindex_ == n.uindex_;
      }

      /** Test whether this node is less than @a n in a global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any geometric meaning.
       *
       * The node ordering relation must obey trichotomy: For any two nodes x
       * and y, exactly one of x == y, x < y, and y < x is true.
       * @param[in] n Node to compare against.
       *
       * @return True if this Node is less than @a n, false otherwise. 
       */
      bool operator<(const Node& n) const {
        bool result;

        // Test the global order prioritizing graph pointer
        if(this->graph_ < n.graph_)
          result = true;
        else if(this->graph_ == n.graph_)
          result = this->uindex_ < n.uindex_;
        else
          result = false;

        return result;
      }

     private:
       // Allow Graph to access Node's private member data and functions.
       friend class Graph;

       /** Pointer to graph to which this node belongs. */
       Graph* graph_;

       /** Unique index identifying this node. */    
       size_type uindex_;

       /** Private constructor to create Node objects.
        * 
        * Node objects should be created by the user with the add_node method 
        * of the Graph class.
        * @param[in] newgraph Pointer to the Graph this Node belongs to.
        * @param[in] id Unique index identifying Node in the internal 
        *               Node vector.
        * @pre @a newgraph is a Graph that contains the internal_node 
        *         corresponding to this Node in the @a id position of its 
        *         internal_node vector.
        * @post This Node is a valid Node of Graph @a newgraph.
        */
       Node(const Graph* newgraph, size_type id) 
         : graph_(const_cast<Graph*>(newgraph)), uindex_(id) {
       }

       /** Find the node data represented by this Node.
        *       
        * This method accesses the graph and searches for the 
        * internal representation whose unique index matches the 
        * identifier for this node.
        * @pre This Node's unique index is within the range of the Graph.
        *
        * @return A reference to the internal representation of this Node. 
        */
       internal_node& fetchNode() const {
         // Return the node whose index matches the caller unique index
         return graph_->nodes_[uindex_];
       }
    };

    /** Return the number of nodes in the graph.
     *
     * @return Return the number of Nodes in this Graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
      // The number of nodes is determined by the number of
      // entries in the unique key mapping.
      return this->nUid_.size();
    }

    /** Synonym for size().
     *
     * @return Return the number of Nodes in this Graph.
     *
     * Complexity: O(1).
     */
    size_type num_nodes() const {
      return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] value    The new node's value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     * @post The adjacency matrix of the graph has a new entry for the 
     *       added node.
     * @post The unique key mapping of the graph has a new entry for the 
     *       added node.
     *
     * @return Return the newly added node.
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, 
                  const node_value_type& value = node_value_type() ) {
      // The new node index and the new unique index
      size_type newIdx = this->num_nodes();
      size_type newUIdx = this->nodes_.size();

      // Create a new internal_node with the given position.
      internal_node new_internalnode = {newIdx, position, value};

      // Add the new entry into the unique key mapping
      this->nUid_.push_back(newUIdx);
      // Add the internal_node to the node Vector.
      this->nodes_.push_back(new_internalnode);

      // Add the node to the adjacency matrix
      this->adj_[newUIdx];

      // Return a Node object with this graph's pointer and the
      // new Node index.
      return Node(this, newUIdx); 
    }

    /** Determine if a Node belongs to this Graph.
     * @param[in] n Node to check whether it belongs to the Graph or not.
     * 
     * @return True if @a n is currently a Node of this Graph.
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
      // A Node belongs to this Graph if it points to this Graph
      // and its index is valid, i.e. it's less than the maximum
      // size of the node array.
      return n.graph_ == this && n.index() < this->num_nodes();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * @return Node whose index is @a i.
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
      // Use the unique key mapping to construct the Node
      return Node(this, nUid_[i]);
    }

    /** Remove a node from the graph.
     *  This method deletes all the incident Edges of a Node, and then
     *  deletes that Node from the Graph.
     * @post new Graph.num_nodes() == old Graph.num_nodes() - 1
     * @post new Graph.num_edges() == old Graph.num_edges() - @a n.degree()
     * @post for all i such that old n.index() < i < old Graph.num_nodes(),
     *       new node(i).index() == old node(i).index - 1
     * @post The unique index mapping for Nodes is updated.
     * @post All outstanding Nodes equal to @a n are invalid.
     * @post All outstanding Edges incident on Node @a n are invalid.
     * @post All NodeIterators pointing at node(i).index() >= old n.index() are
     *       invalid.
     *
     * @return Number of nodes deleted from the Graph.
     *
     * Complexity: assuming n.degree() << num_edges(), O(num_nodes()).
     */
    size_type remove_node(const Node& n) {
      size_type result = 0;

      // Check that the Node is in the Graph
      if(this->has_node(n)) {
        size_type node_uidx = n.uindex_;

        // To avoid conflicts on iterators when deleting
        // incident Edges, store the pairs of Nodes of each
        // incident Edge in a vector
        std::vector<std::pair<Node,Node>> tempAdj;

        // Add the incident Edge data to temporary vector  
        for(IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it)
          tempAdj.push_back(std::make_pair((*it).node1(),(*it).node2()));
      
        // Delete each incident Edge
        for(auto it = tempAdj.begin(); it != tempAdj.end(); ++it)
          remove_edge(it->first, it->second);

        // Find the node entry in unique index mapping and delete it
        auto it = nUid_.begin() + nodes_[node_uidx].idx;
        it = nUid_.erase(it);
  
        // Update the Node indexes
        for(; it != nUid_.end(); ++it)
          nodes_[*it].idx -= 1;
    
        result = 1;
      }

      return result;
    }

    /** Remove a node from the graph, using iterators.
     *  This method deletes all the incident Edges of a Node, and then
     *  deletes that Node from the Graph.
     * @post Same postconditions as remove_node(const Node&)
     *
     * @return Iterator pointing to the next Node in the Graph.
     *
     * Complexity: O(N).
     */
    node_iterator remove_node(node_iterator n_it) {
      // Remove the Node
      remove_node(*n_it);
      // Since after removing a Node, we update the indexes of all
      // the Nodes that follow it, we do not need to increase the
      // iterator.
      return n_it; 
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
      /** Construct an invalid Edge.
       * @post Edge.graph_ is a null pointer and the edge unique index 
       *       is max(unsigned int). Edge::node1() and Edge::node2()
       *       are also max(unsigned int)
       *
       * Invalid nodes are represented by a null pointer to the graph
       * and a unique index of -1. Since index is an unsigned int, -1 is 
       * actually the maximum value of the size_type attribute. This should be 
       * fine, unless we have a graph with the maximum amount of nodes. In that 
       * case, the graph could overflow anyway. 
       * @endcode
         */
      Edge() : graph_(nullptr), uindex_(-1), node1_(-1), node2_(-1) {}

      /** Return a node of this Edge 
       *
       * @return Return index of node1. 
       */
      Node node1() const {
        // Return a Node object pointing to this Graph.
        return Node(graph_, node1_);
      }

      /** Return the other node of this Edge.
       *
       * @return Return index of node2. 
       */
      Node node2() const {
        // Return a Node object pointing to this Graph.
        return Node(graph_, node2_);
      }

      /** Return the length of the Edge.
       *
       * The length of the Edge is defined as the 2-norm of the
       * difference between the positions of this Edge's Nodes.
       *
       * @return Return the length of the Edge. 
       */
      double length() const{
        return norm(this->node1().position() - this->node2().position());
      }

      /** Get a reference to this Edge's value.
       * @pre  The calling Edge has a valid unique index i
       *       with 0 <= i < edges_.size().
       *
       * @return The calling Edge's value reference.
       */
      edge_value_type& value() {
        return fetchEdge().value;
      }

      /** Get a constant reference to this Edge's value.
       * @pre  The calling Edge has a valid unique index i
       *       with 0 <= i < edges_.size().
       *
       * @return The calling Edge's value constant reference.
       */
      const edge_value_type& value() const {
        return fetchEdge().value;
      }

      /** Test whether this edge and @a e are equal.
       *
       * Equal edges represent the same undirected edge between two nodes.
       * @pre This Edge and @a e belong to valid Graphs.
       *
       * @return True if the Edge is equal to @a e, false otherwise.
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
       * @pre This Edge and @a e belong to valid Graphs.
       *
       * @return True if the Edge is less than @a e, false otherwise.
       */
      bool operator<(const Edge& e) const {
        bool result;

        // Test the global order prioritizing graph pointer
        if(this->graph_ < e.graph_)
          result = true;
        else if(this->graph_ == e.graph_)
          result = this->uindex_ < e.uindex_;
        else
          result = false;

        return result;
      }

     private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      /** Pointer to graph to which this node belongs. */
      Graph* graph_;

      /** Unique index of the edge. */
      size_type uindex_;

      /** Direction of the edge. */
      size_type node1_;
      size_type node2_; 

      /** Private constructor to create Edge objects.
       * 
       * Edge objects should be created by the user with the add_edge method 
       * of the Graph class.
       * @param[in] newgraph Pointer to the Graph this Edge belongs to.
       * @param[in] id Unique index identifying Edge in the internal 
       *               Node vector.
       * @param[in] nd1 Index of node as returned by Edge::node1()
       * @param[in] nd2 Index of node as returned by Edge::node2()
       * @pre @a newgraph is a Graph that contains the internal_edge 
       *         corresponding to this Edge in the @a id position of its 
       *         internal_edge vector.
       * @post This Edge is a valid Edge of Graph @a newgraph.
       */
      Edge(const Graph* newgraph, size_type id, size_type nd1, size_type nd2) 
        : graph_(const_cast<Graph*>(newgraph)), uindex_(id), 
          node1_(nd1), node2_(nd2) {}

      /** Find the edge data represented by this Edge.
       *       
       * This method accesses the graph and searches for the 
       * internal representation whose index matches the 
       * identifier for this edge.
       * @pre This Edge's index is within the range of the Graph's Edges.
       *
       * @return A reference to the internal representation of this Edge. 
       */
      internal_edge& fetchEdge() const {
        // Return the edge whose index matches the caller index.
        return graph_->edges_[uindex_];
      }
    };

    /** Return the total number of edges in the graph.
     *
     * @return Number of Edges in the Graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
      // The number of Edges is determined by the number of
      // entries in the unique key mapping.
      return this->eUid_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * @return Edge with index @a i.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
      return Edge(this, eUid_[i], edges_[eUid_[i]].node1_id, 
                                  edges_[eUid_[i]].node2_id);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     *
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
      // For now, we have not found the edge.
      bool result = false;

      // Look for the first node index in our adjacency structure
      auto adjA = adj_.find(nUid_[a.index()]);

      // If the first node is in the adjacency structure
      if(adjA != adj_.end()) {
        // Look for the second node
        auto adjB = (adjA->second).find(nUid_[b.index()]);

        if(adjB != (adjA->second).end())
          result = true; // The edge connecting a and b is in the Graph
      }

      return result;
    }

    /** Add an edge to the graph, or return the current edge if it already 
     *  exists.
     * @param[in] a First Node of the Edge to add.
     * @param[in] b Second Node of the Edge to add.
     * @pre @a a and @a b are distinct valid nodes of this graph.
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges()+1.
     * @post If not old has_edge(@a a, @a b), add new edge(@a a, @a b) and
     *                                        new edge(@a b, @a a) to adjacency
     *                                        matrix.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
     */
    Edge add_edge(const Node& a, const Node& b, 
                  const edge_value_type& value = edge_value_type()) {
      // For now, we have not found the edge.
      bool result = false;

      // Look in the adjacency struct for the Edge linking the two nodes.
      // Once we find the edge, obtain its index. Complexity should be O(1). 

     // If we do not find the edge, its index is the last available.
      size_type uidx = edges_.size(); 

      // Look for the first node index in our adjacency structure
      auto adjA = adj_.find(nUid_[a.index()]);

      // If the first node is in the adjacency structure
      if(adjA != adj_.end()) {
        // Look for the second node
        auto adjB = (adjA->second).find(nUid_[b.index()]);

        if(adjB != (adjA->second).end()) {
          // The Edge is in the Graph. Retrieve its unique index
          result = true;
          uidx   = adjB->second;
        }
      }

      // If we did not find the edge, create an internal_edge to represent it
      // and add it to the vector and adjacency structure.
      if(uidx == edges_.size() && result == false) {
        size_type idx  = this->num_edges();

        // Add the index to the unique index mapping
        this->eUid_.push_back(idx);

        internal_edge new_internaledge = {idx, nUid_[a.index()],
                                               nUid_[b.index()], value};
        this->edges_.push_back(new_internaledge);

        // Add the two directions of Edge to the adjacency matrix
        adj_[nUid_[a.index()]][nUid_[b.index()]] = uidx;
        adj_[nUid_[b.index()]][nUid_[a.index()]] = uidx;
    }
  
    return Edge(this,uidx,nUid_[a.index()],nUid_[b.index()]);
  }

  /** Remove an edge from the graph.
   *  This method deletes an Edge from the Graph, as specified by
   *  the Nodes it connects.
   * @post new Graph.num_edges() == old Graph.num_edges() - 1
   * @post The unique index mapping for Edges is updated.
   * @post All outstanding Edges equal to the removed Edge are invalid.
   * @post All EdgeIterators pointing at edges from the same Graph
   *       greater than the removed Edge are invalid.
   * @post All IncidentIterators created from any of the nodes linked
   *       by the deleted Edge are invalid.
   *
   * @return Number of Edges deleted from the Graph.
   *
   * Complexity: O(num_edges()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type result = 0;

    // If the Edge belongs to the Graph, delete it
    if(adj_.find(nUid_[a.index()]) != adj_.end() &&
       adj_[nUid_[a.index()]].find(nUid_[b.index()]) 
          != adj_[nUid_[a.index()]].end()) {

      // Get the Edge unique id
      size_type edge_uidx = adj_[nUid_[a.index()]][nUid_[b.index()]];
      // Remove the entries in the adjacency matrix
      adj_[nUid_[a.index()]].erase(nUid_[b.index()]);
      adj_[nUid_[b.index()]].erase(nUid_[a.index()]);

      // Remove the entry in the unique index mapping
      auto it = eUid_.begin() + edges_[edge_uidx].idx;
      it = eUid_.erase(it);

      // Update the rest of the indexes
      for(; it != eUid_.end(); ++it)
        edges_[*it].idx -= 1;

      result = 1;
    }

    return result;
  }

  /** Remove an edge from the graph.
   *  This method deletes an Edge from the Graph.
   * @post This method has the same postconditions as the
   *       remove_edge(const Edge&).
   *
   * @return Number of Edges deleted from the Graph.
   *
   * Complexity: O(num_edges()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  } 

  /** Remove an edge from the graph, as specified by the iterator.
   *  This method deletes an Edge from the Graph.
   * @post This method has the same postconditions as the
   *       remove_edge(const Edge&).
   *
   * @return An iterator pointing to the next edge in the Graph.
   *
   * Complexity: O(num_edges()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge((*e_it).node1(), (*e_it).node2());
    return e_it; 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Empty the cointainer for nodes and edges.
    this->nodes_.clear();
    this->nUid_.clear();
    this->edges_.clear();
    this->eUid_.clear();
    this->adj_.clear();
  }

  //
  // NODE ITERATOR
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

    /** Access the Node pointed by this iterator.
     * @pre This is a valid Node iterator for this Graph.
     * @pre This iterator points to a Node whose index i is 0<= i < graph.size()
     *
     * @return a Node object pointed by this iterator.
     */
    Node operator*() const {
      return Node(graph_, graph_->nUid_[ptrIdx_]);
    }

    /** Post-increment the Node iterator.
     * @pre This is a valid Node iterator for this Graph.
     * @post The new iterator points to the next Node in the Graph.
     *
     * @return the updated iterator.
     */
    NodeIterator& operator++() {
      ptrIdx_++;
      return *this;
    }

    /** Iterator equality comparator.
     * Two Node Iterators are equal if they belong to the same Graph
     * and they point to the same Node.
     * @pre This is a valid Node iterator for this Graph.
     *
     * @return True if the Node Iterator is equal to @a rhs, false otherwise.
     */
    bool operator==(const NodeIterator& rhs) const {
      //return Node(graph_, graph_->nUid_[ptrIdx_]) == *rhs;
      return graph_ == rhs.graph_ && ptrIdx_ == rhs.ptrIdx_;
    }

   private:
    friend class Graph;

    /** Iterator pointer to Graph */
    const Graph* graph_;

    /** Iterator index to element */
    size_type ptrIdx_;

    /** Private constructor to create Node Iterator objects.
     * 
     * Node Iterator objects should be created by the user with 
     * node_begin method of the Graph class.
     * @param[in] graph Pointer to the Graph this Node Iterator belongs to.
     * @param[in] ptrIdx  Index identifying Node in the internal 
     *                    Node vector this iterator points to.
     * @pre @a graph is a Graph that contains the Node
     *         in the @a ptrIdx position of its 
     *         internal_node vector.
     * @post This Node Iterator points to a valid Node of Graph @a graph.
     */
    NodeIterator(const Graph* graph, size_type ptrIdx) 
                    : graph_(graph), ptrIdx_(ptrIdx) {}
  };

  /** Create a Node Iterator pointing to the first Node in the Graph.
   *
   * @pre This Node Iterator belongs to a valid Graph. 
   *
   * @return Return a Node Iterator pointing to the first Node of Graph. The
   *         first Node is the one with the smallest valid index.
   */
  node_iterator node_begin() const {
    return NodeIterator{this, 0};
  }

  /** Create a Node Iterator pointing to the last Node in the Graph.
   *
   * @pre This Node Iterator belongs to a valid Graph. 
   *
   * @return Return a Node Iterator pointing to the last Node of Graph. The
   *         last Node is the one with the largest valid index.
   */
  node_iterator node_end() const {
    return NodeIterator{this, this->size()};
  }

  //
  // INCIDENT ITERATOR
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

    /** Access the incident Edge pointed by this iterator.
     * An incident Edge is an Edge that connects the Node used
     * to create this iterator with any other Nodes.
     * @pre This is a valid Incident Iterator for this Graph.
     * @pre This iterator points to an Edge whose index i i
     *      is 0<= i < graphnum_eges()
     *
     * @return an Edge object pointed by this iterator, which is 
     *         incident to the Node associated with the iterator..
     */
    Edge operator*() const {
      return Edge(graph_, itIdx_->second, orgNodeId_,itIdx_->first);
    }

    /** Post-increment the Incident Iterator.
     * @pre This is a valid Incident iterator for this Graph.
     * @post The new iterator points to the next incident Edge in the Graph.
     *
     * @return the updated iterator.
     */
    IncidentIterator& operator++() {
      itIdx_++;
      return *this;
    }

    /** Iterator equality comparator.
     * Two Incident Iterators are equal if they belong to the same Graph,
     * they were originated with the same Node, and they point to the same Edge.
     * @pre This is a valid Incident iterator for this Graph.
     *
     * @return True if the Incident Iterator is equal to @a rhs, false 
     *         otherwise.
     */
    bool operator==(const IncidentIterator& rhs) const {
      return this->graph_ == rhs.graph_ && this->orgNodeId_ == rhs.orgNodeId_ 
             && this->itIdx_ == rhs.itIdx_;
    }

   private:
    friend class Graph;
    /** Iterator pointer to Graph */
    const Graph* graph_;

    /** Iterator pointing to edges adjacent to this node */
    std::unordered_map<size_type,size_type>::iterator itIdx_;

    /** Id of origin node */
    size_type orgNodeId_;

    /** Private constructor to create Incident Iterator objects.
     * 
     * Incident Iterator objects should be created by the user with 
     * Node.begin_edge() method of the Node class.
     * @param[in] graph Pointer to the Graph this Incident Iterator belongs to.
     * @param[in] itIdx Iterator over the elements adjacent to Node with
     *                  index @a orgNodeId from the adjacency matrix.
     * @param[in] orgNodeId Index of the Node whose Incident edges we want to
     *                      iterate on.
     * @pre @a graph is a Graph that contains the Node
     *         in the @a orgNodeId position of its internal_node vector.
     * @pre @a itIdx is an iterator over the adjacent elements of Node @a 
     *      @a orgNodeId.
     * @post This Incident Iterator points to a valid incident Edge  of 
     *       Graph @a graph.
     */
    IncidentIterator(const Graph* graph, 
                    std::unordered_map<size_type,size_type>::iterator itIdx,
                    size_type orgNodeId) 
                    : graph_(graph), itIdx_(itIdx), orgNodeId_(orgNodeId) {}

  };

  //
  // EDGE ITERATOR
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

    /** Access the Edge pointed by this iterator.
     * @pre This is a valid Edge iterator for this Graph.
     * @pre This iterator points to an Edge whose index i 
     *      is 0<= i < graph.num_edges()
     *
     * @return an Edge object pointed by this iterator.
     */
    Edge operator*() const {
      return Edge(graph_, graph_->eUid_[ptrIdx_], 
                  graph_->edges_[graph_->eUid_[ptrIdx_]].node1_id, 
                  graph_->edges_[graph_->eUid_[ptrIdx_]].node2_id);
    }

    /** Post-increment the Edge iterator.
     * @pre This is a valid Edge iterator for this Graph.
     * @post The new iterator points to the next Edge in the Graph.
     *
     * @return the updated iterator.
     */
    EdgeIterator& operator++() {
      ptrIdx_++;
      return *this;
    }

    /** Iterator equality comparator.
     * Two Edge Iterators are equal if they belong to the same Graph
     * and they point to the same Edge.
     * @pre This is a valid Edge iterator for this Graph.
     *
     * @return True if the Edge Iterator is equal to @a rhs, false otherwise.
     */
    bool operator==(const EdgeIterator& rhs) const {
      return this->graph_ == rhs.graph_ && this->ptrIdx_ == rhs.ptrIdx_; 
    }

   private:
    friend class Graph;

    /** Iterator pointer to Graph */
    const Graph *graph_;

    /** Iterator index to element */
    size_type ptrIdx_;

    /** Private constructor to create Edge Iterator objects.
     * 
     * Edge Iterator objects should be created by the user with 
     * edge_begin method of the Graph class.
     * @param[in] graph Pointer to the Graph this Edge Iterator belongs to.
     * @param[in] ptrIdx  Unique index identifying Edge in the internal 
     *               Edge vector this iterator points to.
     * @pre @a graph is a Graph that contains the Edge
     *         in the @a ptrIdx position of its 
     *         internal_edge vector.
     * @post This Edge Iterator points to a valid Edge of Graph @a graph.
     */
    EdgeIterator(const Graph* graph, size_type ptrIdx) 
                    : graph_(graph), ptrIdx_(ptrIdx) {}
  };

  /** Create an Edge Iterator pointing to the first Edge in the Graph.
   *
   * @pre This Edge Iterator belongs to a valid Graph. 
   *
   * @return Return an Edge Iterator pointing to the first Edge of Graph. The
   *         first Edge is the one with the smallest valid index.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator{this, 0};
  }

  /** Create an Edge Iterator pointing to the last Edge in the Graph.
   *
   * @pre This Edge Iterator belongs to a valid Graph. 
   *
   * @return Return an Edge Iterator pointing to the last Edge of Graph. The
   *         last Edge is the one with the largest valid index.
   */
  edge_iterator edge_end() const {
    return EdgeIterator{this,this->num_edges()};
  }

 private:

  /** Structure to represent a node internally (Proxy pattern). */
  struct internal_node {
    size_type idx;
    Point location;
    node_value_type value;
  };

  /** Structure to represent an edge internally (Proxy pattern). */
  struct internal_edge {
    size_type idx;
    size_type node1_id;
    size_type node2_id;
    edge_value_type value;
  };

  /** Vector containing the internal representations of Node objects. */
  std::vector<internal_node> nodes_;
  /** Vector containing the unique index mapping of Node objects. */
  std::vector<size_type> nUid_;
  /** Vector containing the internal representations of Edge objects. */
  std::vector<internal_edge> edges_;
  /** Vector containing the unique index mapping of Edge objects. 
   * Although technically not necessary, since we have the adjacency matrix, 
   * this struct helps us traverse the Graph faster using EdgeIterators, and 
   * since it contains integer types, it is not too resource consuming.
   */
  std::vector<size_type> eUid_;

  /** Map representing the adjacency matrix of the Graph. */
  std::unordered_map<size_type,std::unordered_map<size_type,size_type>> adj_;
};

#endif // CME212_GRAPH_HPP
