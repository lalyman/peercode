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

template <typename V, typename E>
class Graph {
private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)

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
    /** Synonym for value in Node */
    using node_value_type = V;


    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;
    /** Synonym for value in Edge */
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
        // no initializer for vector is needed
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
            // HW0: YOUR CODE HERE
        }

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            assert(valid());
            return graph_->node_elements_[uid_].position_;
        }

        Point& position() {
            assert(valid());
            return graph_->node_elements_[uid_].position_;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            assert(valid());
            return graph_->node_elements_[uid_].idx_;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        /** Return this node's value */
        node_value_type& value(){
            return graph_->node_elements_[uid_].value_;
        }

        /** Return this node's value */
        const node_value_type& value() const {
            return graph_->node_elements_[uid_].value_;
        }

        /** Return the number of adjacent nodes */
        size_type degree() const{
            return graph_->adj_[uid_].size();
        }

        /** Return an iterator pointing to the first incident edge of this node. */
        incident_iterator edge_begin() const{
            return IncidentIterator(graph_, uid_, 0);
        }

        /** Returns an iterator referring to the past-the-end incident edge of this node. */
        incident_iterator edge_end() const{
            return IncidentIterator(graph_, uid_, degree());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            assert(valid());
            return (n.uid_ == uid_ and n.graph_ == graph_);
            // can access private member of objects of the same class
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
            assert(valid());
            bool cond1 = ((graph_ == n.graph_) && (index() < n.index()));
            bool cond2 = graph_ < n.graph_;
            return cond1 || cond2;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        friend class Edge;
        graph_type* graph_;   // pointer to parent graph
        size_type uid_;     // unique identifier of node
        Node(const graph_type* graph, size_type uid)    // valid constructor
                : graph_(const_cast<Graph*>(graph)), uid_(uid){}

        // test if a node is valid using representation invariant
        bool valid() const {
            return uid_ < graph_->node_elements_.size()   // unsigned > 0 == true
                   && graph_-> node_elements_[uid_].idx_ < graph_->idx_2_uid_node_.size()
                   && graph_->idx_2_uid_node_[graph_->node_elements_[uid_].idx_] == uid_;
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
        return idx_2_uid_node_.size();
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
    Node add_node(const Point& position, const node_value_type& value = node_value_type()){
        internal_node inode;
        inode.position_ = position;
        inode.value_ = value;
        inode.idx_ = idx_2_uid_node_.size();
        node_elements_.push_back(inode);    // add internal node to the list

        idx_2_uid_node_.push_back(node_elements_.size() - 1);   // add node to id map
        adj_.resize(node_elements_.size());
        return Node(this, node_elements_.size() - 1);        // return valid node
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        return n.valid();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(0 <= i and i < idx_2_uid_node_.size());
        return Node(this, idx_2_uid_node_[i]);
    }

    /** Remove the given Node @a n.
     * @param n Node to be removed from graph.
     * @return 0 if remove_node fails, 1 if succeed.
     *
     * @pre @a n is a valid Node in this graph.
     * @post  new num_nodes == old num_nodes() - 1;
     * @post  Edges associated with this node are deleted.
     *        new num_edges() == old num_edges() - n.degree();
     * @post  NodeIterator pointing to the beyond nodes are invalidated. All
     *        NodeIterators pointing to Nodes before the removed one are still valid.
     *
     * Complexity:  O(num_nodes() + num_edges())
     * to TA: for my implementation, remove_node is more complex then remove_edge,
     * remove_edge is much more efficient then required.
     */
    size_type remove_node(const Node& n){
        if (has_node(n) && n.valid()){
          size_type index = n.index();
          std::vector<size_type> e_indices;   // vector to store indices of edges to be removed
          std::vector<std::vector<size_type>>  ie_indices; // incident edges indeces, has the same size as n_uid
          if(n.degree() > 0){
              for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
                  Edge e = *it;   // incident edge
                  e_indices.push_back(e.index());

                  std::vector<size_type> v;
                  v.push_back(e.node2().uid_);

                  for(auto ait = adj_[e.uid_n2_].begin(); ait != adj_[e.uid_n2_].end(); ++ait){
                      if((*ait) == e.uid_){
                          v.push_back(ait - adj_[e.uid_n2_].begin());
                      }
                  }
                  assert(v.size() == 2);
                  ie_indices.push_back(v);
              }
          }
          adj_[n.uid_].clear();
          // erase from adjacency mapping
          for(auto it = ie_indices.begin(); it != ie_indices.end(); ++it){
              std::vector<size_type> v = *it;
              adj_[v[0]].erase(adj_[v[0]].begin() + v[1]);
          }
          // manually change index
          int count = 0;
          for(auto it = e_indices.begin(); it != e_indices.end(); ++it){
              size_type idx = *it;
              // manually change all following edge index
              for(auto eit = edge_elements_.begin(); eit!=edge_elements_.end(); ++eit){
                  if((*eit).idx_ > idx){
                    (*eit).idx_ -= 1;
                  }
              }
              for(auto it = e_indices.begin(); it != e_indices.end(); ++it){
                  if (*it > idx){
                      *it -= 1;
                  }
              }

              //remove from index mapping
              idx_2_uid_edge_.erase(idx_2_uid_edge_.begin() + idx - count);
              count++;
          }

          // manually change index from internal_edge
          for(auto eit = edge_elements_.begin(); eit!=edge_elements_.end(); ++eit){
            if((*eit).idx_n1_ >index){(*eit).idx_n1_ -= 1;}
            if((*eit).idx_n2_ >index){(*eit).idx_n2_ -= 1;}
          }

          // manually change index of all internal nodes following erased node
          for (auto it = idx_2_uid_node_.begin(); it != idx_2_uid_node_.end(); ++it){
              if (node_elements_[*it].idx_ > index){node_elements_[*it].idx_ -= 1;}
          }
          // erase node from mapping
          idx_2_uid_node_.erase(idx_2_uid_node_.begin() + index);
          return 1;
        }

        return 0;
    }
    /** Remove the given Node @a n.
     * @param n_it NodeIterator pointing to the Node to be removed from graph.
     * @return An NodeIterator pointing to the new location of the Node that followed
     *         the removed Node. It is the end if the operation removed the last Node
     *         in the graph.
     *
     * @pre @a n_it pointing to a valid Node in this graph.
     * @post  new num_nodes == old num_nodes() - 1;
     * @post  Edges associated with this node are deleted.
     *        new num_edges() == old num_edges() - n.degree();
     * @post  NodeIterator pointing to the beyond nodes are invalidated. All
     *        NodeIterators pointing to Nodes before the removed one are still valid.
     *
     * Complexity:  O(num_nodes() + num_edges())
     * to TA: for my implementation, remove_node is more complex then remove_edge,
     * remove_edge is much more efficient then required.
     */
    node_iterator remove_node(node_iterator n_it){
      Node n = *n_it;
      assert(n.valid());
      size_type index = n.index();
      std::vector<size_type> e_indices;   // vector to store indices of edges to be removed
      std::vector<std::vector<size_type>>  ie_indices; // incident edges indeces, has the same size as n_uid
      if(n.degree() > 0){
          for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
              Edge e = *it;   // incident edge
              e_indices.push_back(e.index());

              std::vector<size_type> v;
              v.push_back(e.node2().uid_);

              for(auto ait = adj_[e.uid_n2_].begin(); ait != adj_[e.uid_n2_].end(); ++ait){
                  if((*ait) == e.uid_){
                      v.push_back(ait - adj_[e.uid_n2_].begin());
                  }
              }
              assert(v.size() == 2);
              ie_indices.push_back(v);
          }
      }
      adj_[n.uid_].clear();
      // erase from adjacency mapping
      for(auto it = ie_indices.begin(); it != ie_indices.end(); ++it){
          std::vector<size_type> v = *it;
          adj_[v[0]].erase(adj_[v[0]].begin() + v[1]);
      }
      // manually change index
      int count = 0;
      for(auto it = e_indices.begin(); it != e_indices.end(); ++it){
          size_type idx = *it;
          // manually change all following edge index
          for(auto eit = edge_elements_.begin(); eit!=edge_elements_.end(); ++eit){
              if((*eit).idx_ > idx){
                (*eit).idx_ -= 1;
              }
          }
          for(auto it = e_indices.begin(); it != e_indices.end(); ++it){
              if (*it > idx){
                  *it -= 1;
              }
          }

          //remove from index mapping
          idx_2_uid_edge_.erase(idx_2_uid_edge_.begin() + idx - count);
          count++;
      }
      for(auto eit = edge_elements_.begin(); eit!=edge_elements_.end(); ++eit){
        if((*eit).idx_n1_ >index){(*eit).idx_n1_ -= 1;}
        if((*eit).idx_n2_ >index){(*eit).idx_n2_ -= 1;}
      }

      // manually change index of all internal nodes following erased node
      for (auto it = idx_2_uid_node_.begin(); it != idx_2_uid_node_.end(); ++it){
          if (node_elements_[*it].idx_ > index){node_elements_[*it].idx_ -= 1;}
      }
      // erase node from mapping
      idx_2_uid_node_.erase(idx_2_uid_node_.begin() + index);
      return node_begin() + index;
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
            // HW0: YOUR CODE HERE
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return Node(graph_, uid_n1_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return Node(graph_, uid_n2_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            bool cond1 = (e.node1() == this->node1()) && (e.node2() == this->node2());
            bool cond2 = (e.node2() == this->node1()) && (e.node1() == this->node2());
            return (cond1 || cond2);
        }

        /** Return the index of this edge */
        size_type index() const {
            return graph_->edge_elements_[uid_].idx_;
        }

        edge_value_type& value() {
            return graph_->edge_elements_[uid_].value_;
        }



        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
          bool cond1 = ((graph_ == e.graph_) && (index() < e.index()));
          bool cond2 = graph_ <e. graph_;
          return cond1 || cond2;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        graph_type* graph_;   // pointer to the parent graph
        size_type uid_;       // unique identifier of this edge
        size_type uid_n1_;    // unique identifier for node 1
        size_type uid_n2_;    // unique identifier for node 2
        Edge(const graph_type* graph, size_type uid, size_type uid_node1, size_type uid_node2)
                : graph_(const_cast<Graph*>(graph)), uid_(uid), uid_n1_(uid_node1), uid_n2_(uid_node2){}
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
        return idx_2_uid_edge_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(i >= 0 and i < num_edges());
        return Edge(this, idx_2_uid_edge_[i],
                    idx_2_uid_node_[edge_elements_[idx_2_uid_edge_[i]].idx_n1_],
                    idx_2_uid_node_[edge_elements_[idx_2_uid_edge_[i]].idx_n2_]);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        if(has_node(a) && has_node(b)){
            if(a == b){ return false;}
            int i = edge_index(a, b);
            return (i > -1);
        }else{
            return false;
        }
    }

    /** Return the index of an edge in this graph, return -1 if not exist
     * @pree @a a and @b are valid nodes of this graph
     * @return Index of edge if it exists
     *
     * Complexity: No more tha O(num_nodes() + num_edges())
     */
     int edge_index(const Node&a, const Node& b) const {
        assert(has_node(a) && has_node(b));
        if(a != b){
          for(auto it = adj_[a.uid_].begin(); it != adj_[a.uid_].end(); it++){
              if(edge_elements_[*it].idx_n1_ == b.index() || edge_elements_[*it].idx_n2_ == b.index()){
                  return edge_elements_[*it].idx_;
              }
          }
        }
        return -1;
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
        assert(has_node(a) && has_node(b) && (!(a == b)));
        if (has_edge(a, b)) {
            int eid = edge_index(a, b);
            return Edge(this, idx_2_uid_edge_[eid],
                        idx_2_uid_node_[edge_elements_[idx_2_uid_edge_[eid]].idx_n1_],
                        idx_2_uid_node_[edge_elements_[idx_2_uid_edge_[eid]].idx_n2_]);
        }else{
            internal_edge edge;
            edge.idx_n1_ = a.index();
            edge.idx_n2_ = b.index();
            edge.idx_ = idx_2_uid_edge_.size();
            edge_elements_.push_back(edge);     // push back unique edge

            idx_2_uid_edge_.push_back(edge_elements_.size() - 1);   // add map of edge
            adj_[a.uid_].push_back(edge_elements_.size() - 1);
            adj_[b.uid_].push_back(edge_elements_.size() - 1);   // add to adjency vector

            return Edge(this, edge_elements_.size() - 1,
                        idx_2_uid_node_[edge_elements_[edge_elements_.size() - 1].idx_n1_],
                        idx_2_uid_node_[edge_elements_[edge_elements_.size() - 1].idx_n2_]);

        }
    }

    /** Remove the Edge constructed by Node a and Node b.
     * @param a Node of the Edge to be removed from graph.
     * @param b Node of the Edge to be removed from graph.
     * @return 0 if remove_edge fails, 1 if succeed.
     *
     * @pre @a a and @a b are two valid Node in this graph.
     * @pre @a a and @a b have a valid Edge in this graph.
     * @post  new num_edges() == old num_edges() - 1;
     * @post  EdgeIterator pointing to the beyond Edges are invalidated. All
     *        EdgeIterator pointing to Edges before the removed one are still valid.
     *
     * Complexity:  O(num_edges())
     */
    size_type remove_edge(const Node& a, const Node& b){
        int e_idx = edge_index(a, b);
        if(e_idx != -1){
          // remove edge from adj_ mapping
          size_type offset = 0;
          bool flag = false;
          for(auto it = adj_[a.uid_].begin(); it!=adj_[a.uid_].end(); ++it){
              if(*it == idx_2_uid_edge_[e_idx]){
                  offset = it - adj_[a.uid_].begin();
                  flag = true;
                  break;
              }
          }
          assert(flag);
          adj_[a.uid_].erase(adj_[a.uid_].begin() + offset);
          flag = false;
          for(auto it = adj_[b.uid_].begin(); it!=adj_[b.uid_].end(); ++it){
              if(*it == idx_2_uid_edge_[e_idx]){
                  offset = it - adj_[b.uid_].begin();
                  flag = true;
                  break;
              }
          }
          assert(flag);
          adj_[b.uid_].erase(adj_[b.uid_].begin() + offset);

          // manually change index of following edges
          for(auto it = edge_elements_.begin(); it!=edge_elements_.end(); ++it){
              if((*it).idx_ > (unsigned)e_idx){(*it).idx_ -= 1;}
          }
          idx_2_uid_edge_.erase(idx_2_uid_edge_.begin() + e_idx);  // remove from edge mapping
          return 1;
        }
        // return e_idx;
        return 0;
    }

    /** Remove the given Edge from graph.
     * @param e  Edge to be removed from graph.
     * @return 0 if remove_edge fails, 1 if succeed.
     *
     * @pre @a e is a valid Edge in this graph.
     * @post  new num_edges() == old num_edges() - 1;
     * @post  EdgeIterator pointing to the beyond Edges are invalidated. All
     *        EdgeIterator pointing to Edges before the removed one are still valid.
     *
     * Complexity:  O(num_edges())
     */
    size_type remove_edge(const Edge& e){
        return remove_edge(e.node1(), e.node2());
    }

    /** Remove the given Edge from graph.
     * @param e_it  EdgeIterator pointing to the Edge to be removed from graph.
     * @return  An EdgeIterator pointing to the first Edge in the graph.
     *
     * @pre @a e_it pointing to a valid Edge in this graph.
     * @post  new num_edges() == old num_edges() - 1;
     * @post  EdgeIterator pointing to the beyond Edges are invalidated. All
     *        EdgeIterator pointing to Edges before the removed one are still valid.
     *
     * Complexity:  O(num_edges())
     */
    edge_iterator remove_edge(edge_iterator e_it){
        remove_edge((*e_it).node1(), (*e_it).node2());
        return edge_begin();
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        node_elements_.clear();
        idx_2_uid_node_.clear();
        adj_.clear();
        edge_elements_.clear();
        idx_2_uid_edge_.clear();
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

        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        /** Deference of Node Iterator */
        Node operator*() const {
            return graph_->node(idx_);
        }

        /** Increment of Node Iterator */
        NodeIterator& operator++(){
            ++idx_;
            return *this;
        }

        /** Test if two Node Iterators are the same
         *  Same Node Iterators have point the same node in the same graph
         * */
        bool operator==(const NodeIterator& ni) const{
            return (ni.graph_ == graph_) && (ni.idx_ == idx_);
        }

    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        graph_type* graph_;
        size_type idx_;  // the index of node this iterator is pointing to
        NodeIterator(const graph_type* graph, size_type idx)   // constructor through graph
                : graph_(const_cast<Graph*>(graph)), idx_(idx){}
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return an iterator pointing to the first node in the graph. */
    node_iterator node_begin() const{
        return NodeIterator(this, 0);
    }

    /** Returns an iterator referring to the past-the-end node in the graph. */
    node_iterator node_end() const {
        return NodeIterator(this, idx_2_uid_node_.size());
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
        /** Deference this iterator, return the edge this iterator is referring to. */
        Edge operator*() const{
            size_type u1 = graph_->idx_2_uid_node_[graph_->edge_elements_[graph_->adj_[uid_parent_n_][idx_]].idx_n1_];
            size_type u2 = graph_->idx_2_uid_node_[graph_->edge_elements_[graph_->adj_[uid_parent_n_][idx_]].idx_n2_];
            if(u1 == uid_parent_n_){
                return Edge(graph_, graph_->adj_[uid_parent_n_][idx_], u1, u2);
            }else{
                return Edge(graph_, graph_->adj_[uid_parent_n_][idx_], u2, u1);
            }
        }

        /** Increment the iterator by one, return an iterator pointing to the next incident edge. */
        IncidentIterator& operator++(){
            ++idx_;
            return *this;
        }

        /** Test if two incident iterators are the same
         *  Same incident iterators points to the same edge in the same graph.
         * */
        bool operator==(const IncidentIterator& ii) const{
            return (ii.graph_ == graph_) and (ii.idx_ == idx_);
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
        graph_type* graph_;
        size_type uid_parent_n_;   // unique identifier of the parent node
        size_type idx_;   // index of the edge among all incident edges
        IncidentIterator(const graph_type* graph, size_type uid, size_type idx):
                graph_(const_cast<Graph*>(graph)), uid_parent_n_(uid), idx_(idx){}
    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private totally_ordered<EdgeIterator>{

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

        /** Deference an edge iterator. Return the edge its pointing to. */
        Edge operator*() const{
            return Edge(graph_, graph_->idx_2_uid_edge_[idx_],
                             graph_->idx_2_uid_node_[graph_->edge_elements_[graph_->idx_2_uid_edge_[idx_]].idx_n1_],
                             graph_->idx_2_uid_node_[graph_->edge_elements_[graph_->idx_2_uid_edge_[idx_]].idx_n2_]);
        }

        /** Increment the iterator by 1. Return an iterator refering to the next edge in graph. */
        EdgeIterator& operator++(){
            ++idx_;
            return *this;
        }

        /** Test if two edge iterators are the same.
         *  Same edge iterators belong to the same graph and point to the same edge
         * */
        bool operator==(const EdgeIterator& ei) const{
            return (ei.graph_ == graph_) && (ei.idx_ == idx_);
        }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        graph_type* graph_;
        size_type idx_; // index of edge in graph
        EdgeIterator(const graph_type* graph, size_type idx):
                graph_(const_cast<Graph*>(graph)), idx_(idx){}
    };

    // HW1 #5: YOUR CODE HERE
    /** Return an iterator pointing to the first edge in graph. */
    edge_iterator edge_begin() const{
        return EdgeIterator(this, 0);
    }

    /** Returns an iterator referring to the past-the-end edge in graph. */
    edge_iterator edge_end() const{
         return EdgeIterator(this, idx_2_uid_edge_.size());
    }

private:
    // internal type for node and edge
    struct internal_node{
        Point position_;  // 3d position of the node
        node_value_type value_;  // value of the node
        size_type idx_;   // node index in the graph
    };
    struct internal_edge{
        size_type idx_n1_;  // uid of first node
        size_type idx_n2_;  // uid of second node
        size_type idx_;   // edge index in the graph
        edge_value_type value_;
    };

    std::vector<internal_node> node_elements_;   // vector for all nodes (include erased nodes)
    std::vector<size_type> idx_2_uid_node_;   // map from node index to node unique identifier

    std::vector<std::vector<size_type> > adj_;  // vector for all edges indexed by unique nodes
    std::vector<internal_edge> edge_elements_;   // vector for all edges, indexed by unique identifier
    std::vector<size_type> idx_2_uid_edge_;   // map from edge index to edge unique identifier

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
