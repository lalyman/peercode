#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <iostream>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E = double>
class Graph {
private:

    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph<V,E>;

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

    //node value type
    using node_value_type = V;
    //edge value type
    using edge_value_type = E;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph(){}

    /** Default destructor */
    ~Graph() = default;




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
        Node(): graph_(nullptr),uid_(0) {}
        Node(const graph_type* graph, size_type uid):
                graph_(const_cast<graph_type *> (graph)), uid_(uid){}



        /** Return this node's position.
         * @pre 0 < this->index() < graph node number
         * @return node's position
         */
        Point& position() {
            if(this->index() >= graph_->num_nodes())
                std::cout <<"ERROR: index is greater than num_nodes()" << std::endl;
            return graph_->nodes_[uid_].position_;
        }

        /** Return this node's position.
         * @pre 0 < this->index() < graph node number
         * @return node's position
         */
        const Point& position() const {
            if(this->index() < graph_->num_nodes())
                return graph_->nodes_[uid_].position_;
            assert(false);
        }

        /** Return this node's index, a number in the range [0, graph_size).
         * @return node id
         */
        size_type index() const {
            return graph_->nodes_[uid_].idx_;
        }

        /** Return this node's index, a number in the range [0, graph_size).
         *  @return node id
         */
        size_type &index() {
            return graph_->nodes_[uid_].idx_;
        }


        /** Return this node's unique id, a number in the range [0, total number of nodes).
         * @return node uid_
         */
        size_type uid() const {
            return uid_;
        }

        /** Return this node's unique id, a number in the range [0, total number of nodes).
         * @return node uid_
         */
        size_type & uid(){
            return uid_;
        }

        /** Return this node's graph pointer.
         * @return graph pointer
         */
        graph_type* graph() const {
            return graph_;
        }



        /** Return this node's value by reference.
         *  @return node value
         */
        node_value_type & value () {return graph_->nodes_[uid_].value_;}

        /** Return this node's value by constant reference.
         * @return node value
         */
        const node_value_type & value () const {return graph_->nodes_[uid_].value_;};


        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return (this->graph_ ==  n.graph() && this->uid_ == n.uid_);
        }

        /** Test whether this node is less than @a n in a global order.
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         * @param[in] Node
         * @return true if two nodes are in the same graph and uid_ < n.index
         */
        bool operator<(const Node& n) const {
            if(this->graph_ ==  n.graph())
                return (this->index() < n.index());
            else
                return this->graph_ < n.graph();//nodes in different graph, compare the graph pointer
        }

        /**
         * @return the number of incident edges of the current node
         */
        size_type degree() const {
            return graph_->incident_edges_[this->index()].size();
        };

        /**
        * @return the begin iterator of the node's incident iterator
        */
        incident_iterator edge_begin() const {return IncidentIterator(graph_, this->index(), 0);};

        /**
        * @return the end iterator of the node's incident iterator
        */
        incident_iterator edge_end() const {return IncidentIterator(graph_, this->index(), degree());};

    private:
        // Allow Graph to access Node's private member data and functions.

        graph_type * graph_;//graph pointer
        size_type uid_;     //unique node id


        friend class Graph;

    };

    /** @return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return i2n_.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] node_value_type The new node's value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node ( const Point & position, const node_value_type & value = node_value_type()) {
        size_type  num_nodes_ = i2n_.size();
        size_type  num_uids_ = nodes_.size();

        nodes_.push_back(Internal_Node(position, value, num_nodes_));

        i2n_.push_back(num_uids_);

        incident_edges_.push_back( std::vector<Internal_Edge>());

        return Node(this, num_uids_);
    }



    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        return (n.index() < num_nodes());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        if(i >= num_nodes()) {
            std::cout << "In func node: @pre 0 <= @a i < num_nodes()" <<std::endl;
            return Node();
        }
        else
            return Node(this, i2n_[i]);
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
    class Edge: private totally_ordered<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge(): graph_(nullptr), uid_(0) ,  n1_(0), n2_(0){//, value_(edge_value_type()){
        }
        /** Construct a valid Edge. */
        Edge(const graph_type* graph, size_type n1, size_type n2, size_type uid):
                graph_(const_cast<graph_type *>(graph)),uid_(uid), n1_(n1), n2_(n2){
        }

        /** Return a node of this Edge */
        Node node1() const {
            return Node(graph_, n1_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, n2_);
        }

        /** Return the first node number of this Edge */
        size_type& n1()  {
            return n1_;
        }

        /** Return the second node number of this Edge */
        size_type& n2()  {
            return n2_;
        }

        /** Return this edge's graph pointer.
         * @return graph pointer
         */
        graph_type* graph() const {
            return graph_;
        }

        /** Test whether this edge and @a e are equal.
         *
         * @return ture if edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            return ((e.node1()==this->node1() && e.node2()==this->node2())
                    || (e.node1()==this->node2() && e.node2()==this->node1()));
        }



        /** @return this edge's index, a number in the range [0, graph_edge_size). */
        size_type index() const {
            return size_type(uid_);
        }


        /** @return this edge's index, a number in the range [0, graph_edge_size). */
        size_type &index() {
            return uid_;
        }


        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         * @return true if the index of the edge <  e.index
         */
        bool operator<(const Edge& e) const {
            if(this->graph_ == e.graph())
                return (this->index() < e.index());
            else
                return this->graph_ < e.graph();//edges in different graph, compare the graph pointer
        }


        edge_value_type & value() {return graph_->edge_values_[uid_];}
        const edge_value_type & value() const {return graph_->edge_values_[uid_]; }

    private:
        graph_type * graph_;
        size_type uid_;
        size_type n1_, n2_;//node uid!!
        //edge_value_type value_;
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

    };

    /** @return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return edges_.size();
    }

    /** @return  the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        if(i >= num_edges()) {
            std::cout << "In func edge: @pre 0 <= @a i < num_edges()" <<std::endl;
            return Edge();
        }
        else{
            return Edge(this, edges_[i].uid_n1_, edges_[i].uid_n2_, edges_[i].idx_);
        }

    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        for(auto it = a.edge_begin(); it != a.edge_end(); ++it)
            if((*it).node2() == b)   return true;
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


        for(auto it = a.edge_begin(); it != a.edge_end(); ++it)
            if((*it).node2() == b)   return *it;


        size_type num_edges = edges_.size();
        Internal_Edge new_edge = Internal_Edge(a.uid(), b.uid(), num_edges);
        edges_.push_back(new_edge);

        incident_edges_[a.index()].push_back(new_edge);
        incident_edges_[b.index()].push_back(Internal_Edge(b.uid(), a.uid(), num_edges));
        edge_values_.push_back(edge_value_type());

        return Edge(this, a.uid(), b.uid(), num_edges);

    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes_.clear();     //Node vector
        i2n_.clear();       //Node map

        edges_.clear();     //Edge vector
        incident_edges_.clear();  //Node incident edge vector
        edge_values_.clear();     //Edge values
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator: private totally_ordered<NodeIterator>{
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid NodeIterator. */
        NodeIterator() {}

        /** Construct a NodeIterator.
         *  @param[in] graph: graph pointer
         *  @param[in] p: Node id
         */
        NodeIterator(const graph_type * graph, size_type p): graph_(const_cast<graph_type *>(graph)), p_(p) {}
        /** Dereference a NodeIterator.
        *   @return the the node
        */
        Node operator *() const { return Node(graph_, graph_->i2n_[p_]); }
        /** Compare two iterators.
         *  @param[in] constant reference of a Nodeiterator
         *  @return true, if these two iterator are the same
         */
        bool operator ==( const NodeIterator & x ) const { return p_ == x.p_ ; }

        /** Self increment.
         *  @return next iterator
         */
        NodeIterator& operator ++(){ p_++;  return *this; }


    private:
        friend class Graph;
        graph_type* graph_;
        size_type p_ ;

    };

    /**
     *  @return the head of the node iterator
     */
    node_iterator node_begin() const {return NodeIterator(this, 0);}
    /**
     *  @return the tail of the node iterator
     */
    node_iterator node_end() const {return NodeIterator(this, num_nodes());}

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator:private totally_ordered<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {}

        /** Construct an IncidentIterator.
         *  @param[in] graph: graph pointer
         *  @param[in] idx: node index
         *  @param[in] p: Edge id
         */
        IncidentIterator(const graph_type * graph, size_type idx, size_type p):graph_(const_cast<graph_type *>(graph)), idx_(idx), p_(p) {}

        /** Dereference.
         *   @return the constant reference of the edge
         */
        Edge operator*() const{ return Edge(graph_, graph_->incident_edges_[idx_][p_].uid_n1_,
                                            graph_->incident_edges_[idx_][p_].uid_n2_, graph_->incident_edges_[idx_][p_].idx_);}
        /** Dereference a IncidentIterator.
         *   @return the reference of the node
         */
        IncidentIterator& operator++(){ p_++; return *this;}


        /** Compare two iterators.
         *  @param[in] constant reference of a Incidentiterator
         *  @return true, if these two iterator are the same
         */
        bool operator==(const IncidentIterator& x) const{ return p_ == x.p_;}

    private:
        friend class Graph;
        graph_type* graph_;
        size_type idx_;
        size_type p_ ;

    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator: private totally_ordered<EdgeIterator>{
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {}
        /** Construct an EdgeIterator.
         *  @param[in] graph: graph pointer
         *  @param[in] p: Edge id
         */
        EdgeIterator(const graph_type * graph, size_type p): graph_(const_cast<graph_type *>(graph)), p_(p) {}

        /** Dereference an EdgeIterator.
         *   @return the constant reference of the edge
         */

        Edge operator*() const{ return Edge(graph_, graph_->edges_[p_].uid_n1_,
                                            graph_->edges_[p_].uid_n2_, graph_->edges_[p_].idx_);}
        /** Self increment.
         *  @return next iterator
         */
        EdgeIterator& operator++(){p_++; return *this;}

        /** Compare two iterators.
         *  @param[in] constant reference of a Edgeiterator
         *  @return true, if these two iterator are the same
         */
        bool operator==(const EdgeIterator& x) const{return p_ == x.p_;}

    private:
        friend class Graph;
        graph_type* graph_;
        size_type p_ ;

    };


    /**
     *  @return the head of the edge iterator
     */
    edge_iterator edge_begin() const{ return EdgeIterator(this, 0);}
    /**
     *  @return the head of the edge iterator
     */
    edge_iterator edge_end() const{ return EdgeIterator(this, num_edges());}


    /** Remove the node and its related edges from the graph.
     * @param[in] n: node 1
     * @return id: the node index
     *
     * erase all its related edges, erase its incident_edges_ list, see void remove_incident_edge ( const Node & n1, const Node & n2)
     * erase the node from node map, but keep it in the node list for efficiency(for the case node value is large).
     * update node index
     *
     * @post  If the node does not exist, do nothing
     * @post  All  node unique index in i2n_  with  position i < this node's index are not effected.
     *        All  node unique index  with  position i >this node's index are  moved  to  position i-1
     * @post  All  node with  index greater than this node's index are affected, their index is reduced by 1, but
     *        their unique id is remained.
     *
     * @post  All  nodes' _incident_edges_ with node index greater than this node are not effected.
     *        All edges_ with node index greater than the node's index are not effected.
     *
     * @note  Invalidates  all  iterators
     * At most O(2num_nodes + 4num_edges)  operations.
     */
    size_type remove_node ( const Node & n) {

        size_type id = n.index();
        if (n.graph_ != this || id >= num_nodes()) {
            return id;
        }//The node does not belong to this graph




        //remove edges
        while (n.edge_begin() != n.edge_end()) {
            remove_edge(*(n.edge_begin()));
        }


        incident_edges_.erase(incident_edges_.begin() + id);

        //remove nodes
        i2n_.erase(i2n_.begin() + id);
        for (auto it = i2n_.begin() + id; it != i2n_.end(); ++it)
            nodes_[*it].idx_ -= 1;


        return id;
    }

    /** The same as size_type remove_node ( const Node & n)
     * @param[in] n_it: node iterator
     * @return id: the node index
     */
    size_type remove_node ( node_iterator n_it ){
        const Node n = *n_it;
        return remove_node(n);
    }




    /** Remove the edge from the node incident edge vector.
      * @param[in] n1: node 1
      * @param[in] n2: node 2
      *
      * erase edge(n1, n2) from both n1's and n2's  incident_edges_.
      * @pre  n1 and n2 are valid nodes of the graph.
      *
      * @post  If the edge does not exist, do nothing
      * @post  All  edges in the incident_edges_[n1.index()]  with  position i < this edge are not effected.
      *        All  edges  with  position i >this edge are  moved  to  position i-1.
      * @post  All  edges in the incident_edges_[n2.index()]  with  position i < this edge are not effected.
      *        All  edges  with  position i >this edge are  moved  to  position i-1.
      * @note  Invalidates  all  iterators [this edge, end())
      * At most O(2degree())  operations.
      */
    void remove_incident_edge ( const Node & n1, const Node & n2){
        for(auto it = incident_edges_[n1.index()].begin(); it != incident_edges_[n1.index()].end(); ++it)
            if((*it).uid_n2_ == n2.uid()) {
                incident_edges_[n1.index()].erase(it);
                break;
            }

        for(auto it = incident_edges_[n2.index()].begin(); it != incident_edges_[n2.index()].end(); ++it)
            if((*it).uid_n2_ == n1.uid()) {
                incident_edges_[n2.index()].erase(it);
                break;
            }

    }
    /** Remove the edge from the graph.
      * @param[in] n1: node 1
      * @param[in] n2: node 2
      * @return id: the edge index
      *
      * erase edge(n1, n2) from edges_ vector, edge_values_ vector and both n1's and n2's  incident_edges_
      * and update edges index, .
      * @pre  n1 and n2 are valid nodes of the graph.
      *
      * @post  If the edge does not exist, do nothing
      * @post  All  edges in the edges_, and all values in edge_values_ with  position i < this edge are not effected.
      *        All  edges  with  position i >this edge are  moved  to  position i-1, and index is reduced by 1
      * @post  All  edges in the incident_edges_[n1.index()]  with  position i < this edge are not effected.
      *        All  edges  with  position i >this edge are  moved  to  position i-1, and index is reduced by 1
      * @post  All  edges in the incident_edges_[n2.index()]  with  position i < this edge are not effected.
      *        All  edges  with  position i >this edge are  moved  to  position i-1, and index is reduced by 1
      * @note  Invalidates  all  iterators
      * At most O(4num_edges)  operations.
      */
    size_type remove_edge ( const Node & n1, const Node & n2){
        size_type id {};
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            if(((*it).uid_n1_ == n1.uid() && (*it).uid_n2_ == n2.uid())||
               ((*it).uid_n1_ == n2.uid() && (*it).uid_n2_ == n1.uid())) {
                id = (*it).idx_;
                edges_.erase(it);
                edge_values_.erase(edge_values_.begin() + id);
                break;
            }
        //update edge idx
        for (auto it = edges_.begin() + id; it != edges_.end(); ++it) {
            (*it).idx_ = (*it).idx_ - 1;
        }
        remove_incident_edge(n1,n2);

        //update incident edge idx
        for (auto it = incident_edges_.begin(); it != incident_edges_.end(); ++it) {
            for (auto it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
                if ((*it2).idx_ > id)
                    (*it2).idx_ -= 1;
            }
        }

        return id;
    }
    /**
     * The same as size_type remove_edge ( const Node & n1, const Node & n2)
     * @param e: edge reference
     * @return edge index
     * @pre  e is a valid edge.
     */
    size_type remove_edge ( const Edge & e){
        const Node &n1 = e.node1(), &n2 = e.node2();
        return remove_edge(n1,n2);
    }

    /**
     * The same as size_type remove_edge ( const Node & n1, const Node & n2)
     * @param e_it: edge iterator
     * @return edge index
     * @pre  e is a valid edge iterator.
     */
    edge_iterator remove_edge ( edge_iterator e_it ){
        const Node n1 = (*e_it).node1(), n2 = (*e_it).node2();
        return remove_edge(n1,n2);
    }

    struct Internal_Node{
        Point position_;   //node position
        node_value_type value_;//node value
        size_type idx_;    //node idx in i2n_ array
        Internal_Node(Point position, node_value_type value, size_type idx): position_(position), value_(value), idx_(idx){}

    };

    struct Internal_Edge{
        size_type uid_n1_; //edge first node unique id
        size_type uid_n2_; //edge second node unique id
        size_type idx_;    //edge index, in edge_values_ array
        Internal_Edge(size_type uid_n1, size_type uid_n2, size_type idx): uid_n1_(uid_n1), uid_n2_(uid_n2), idx_(idx){}

    };

private:





    std::vector<Internal_Node> nodes_;     //Node vector
    std::vector<size_type> i2n_;           //idx to uid

    std::vector<Internal_Edge> edges_;     //Edge vector
    std::vector<std::vector<Internal_Edge>> incident_edges_;  //Node incident edge vector
    std::vector<edge_value_type> edge_values_;                //edge value vector


};

#endif // CME212_GRAPH_HPP
