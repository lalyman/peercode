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

template <typename V>
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

    /** Type of modifiable value in nodes */
    using node_value_type = V;


    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        // HW0: YOUR CODE HERE
        num_edges_ = 0;
        num_nodes_ = 0;
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
            isvalid_ = false;
        }

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            return graph_->node_elements_[id_];
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return id_;
        }

        /** Return this node's graph as a pointer to cosnt*/
        const graph_type* graph() const {
            return graph_;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Return the number of incident edges of this node */
        size_type degree() const{
            return graph_->incident_edges_[id_].size();
        }

        /** Begin of iterator for incident edges **/
        incident_iterator edge_begin() const{
            return IncidentIterator(graph_,id_,0);
        }

        /** End of iterator for incident edges */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_,id_,graph_->incident_edges_[id_].size());
        }

        /** Return the value of this node */
        node_value_type& value(){
            return value_;
        }

        /** Return the value of this node as constant */
        const node_value_type& value() const {
            return value_;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            return (n.index() == id_ and n.graph() == graph_);
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
            assert(isvalid_);
            return (id_ < n.id_);
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        bool isvalid_;
        graph_type* graph_;
        size_type id_;
        node_value_type value_;
        /** Constructor of a valid Node */
        Node(const graph_type* graph, size_type id, const node_value_type& value)
                : graph_(const_cast<Graph*>(graph)), id_(id), value_(value), isvalid_(true){
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
    Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
        // HW0: YOUR CODE HERE
        node_elements_.push_back(position);
        nodes_.push_back(Node(this, num_nodes_, value));
        incident_edges_.push_back(std::vector<size_type>());
        ++num_nodes_;
        return nodes_[num_nodes_-1];
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        // (void) n;            // Quiet compiler warning
        return (this == n.graph());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node& node(size_type i) {
        // HW0: YOUR CODE HERE
        assert(0 <= i and i < num_nodes_);        assert(0 <= i and i < num_nodes_);
        return nodes_[i];
    }

    Node node(size_type i) const {
        assert(0 <= i and i < num_nodes_);
        return nodes_[i];
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
            isvalid_ = false;
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            assert(isvalid_);
            if (is_major){
                return graph_->node(graph_->edge_elements_[id_][0]);
            }
            else{
                return graph_->node(graph_->edge_elements_[id_][1]);
            }
        }
        Node& node1() {
            // HW0: YOUR CODE HERE
            assert(isvalid_);
            if (is_major){
                return graph_->node(graph_->edge_elements_[id_][0]);
            }
            else{
                return graph_->node(graph_->edge_elements_[id_][1]);
            }
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            if (is_major){
                return graph_->node(graph_->edge_elements_[id_][1]);
            }
            else{
                return graph_->node(graph_->edge_elements_[id_][0]);
            }

        }
        Node& node2() {
            // HW0: YOUR CODE HERE
            if (is_major){
                return graph_->node(graph_->edge_elements_[id_][1]);
            }
            else{
                return graph_->node(graph_->edge_elements_[id_][0]);
            }

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
            return id_;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            return (id_ < e.index());
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        friend class IncidentIterator;
        // HW0: YOUR CODE HERE
        bool isvalid_;
        bool is_major;  // flag if node1 for this edge is of interest
        graph_type* graph_;
        size_type id_;

        /** Constructor of a valid edge */
        Edge(const graph_type* graph, size_type id, bool major=true)
                : graph_(const_cast<Graph*>(graph)), id_(id), isvalid_(true), is_major(major){
        }

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
        return num_edges_;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        // (void) i;             // Quiet compiler warning
        assert(i >= 0 and i < num_edges_);
        return Edge(this, i);        // Invalid Edge
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE

        assert(has_node(a) && has_node(b) && (!(a == b)));
        int i = edge_index(a, b);
        return (i > -1);
    }

    /** Return the index of an edge in this graph, return -1 if not exist*/
    int edge_index(const Node&a, const Node& b) const {
        assert(has_node(a) && has_node(b) && (!(a == b)));
        for (int i = 0; i < num_edges_; i++){
            if (edge_elements_[i][0] == a.index() and edge_elements_[i][1] == b.index()){
                return i;
            }else if (edge_elements_[i][0] == b.index() and edge_elements_[i][1] == a.index()){
                return i;
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
        assert(this->has_node(a) and this->has_node(b) and (!(a == b)));
        if (this->has_edge(a, b)){
            int eid = edge_index(a, b);
            return Edge(this, eid);
        }
        else{
            edge_elements_.push_back(std::vector<size_type> {a.index(), b.index()});
            incident_edges_[a.index()].push_back(num_edges_);
            incident_edges_[b.index()].push_back(num_edges_);
            ++num_edges_;
            return Edge(this, num_edges_-1);

        }
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        nodes_.clear();
        node_elements_.clear();
        edge_elements_.clear();
        num_nodes_ = 0;
        num_edges_ = 0;
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
            isvalid_ = false;
        }

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** dereference operator for NodeIterator */
        Node operator*() const{
            assert(isvalid_);
            return graph_->nodes_[id_];
        }

        Node& operator*() {
            return graph_->nodes_[id_];
        }

        /** increment NodeIterator, point to the next Node in range */
        NodeIterator& operator++(){
            assert(isvalid_);
            ++id_;
            return *this;
        }

        /** return const graph of this node*/
        const Graph* graph() const {
            return graph_;
        }

        /** Return the index of Node the iterator is pointing to */
        const size_type index() const{
            return id_;
        }

        /** Test whether this iterator and iterator ni is the same
        /*  Same NodeIterator point to the same node in the same graph.
        */
        bool operator==(const NodeIterator& ni) const {
            return id_ == ni.index() and ni.graph_ == graph_;
        }

    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        graph_type* graph_;
        size_type id_;
        bool isvalid_;
        /** Constructor for a valid NodeIterator */
        NodeIterator(const graph_type * graph, size_type index): graph_(const_cast<Graph*>(graph)), id_(index), isvalid_(true){}

    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return the begin of node iterator */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }
    /** Return the end of node iterator */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes_);
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

        /** Return the Edge this iterator is pointing to
        /* Node(1) of this edge is the parent node
        /* Node(2) of this edge is other node
        */
        Edge operator*() const{
            if(node_id_ == graph_->edge_elements_[graph_->incident_edges_[node_id_][node_ie_id]][0]){
                return Edge(graph_, graph_->incident_edges_[node_id_][node_ie_id], true);
            }
            else{
                return Edge(graph_, graph_->incident_edges_[node_id_][node_ie_id], false);
            }
        }

        /** increment the iterator to point to the next edge incidents to this node */
        IncidentIterator& operator++(){
            ++node_ie_id;
            return *this;
        }

        /** Return the graph this iterator is poiting to */
        const Graph* graph() const {
            return graph_;
        }

        /** Return the index of the parent node in graph */
        const size_type node_index() const {
            return node_id_;
        }

        /** Return the index of this edge in incident edges*/
        const size_type index() const{
            return node_ie_id;
        }

        /** Test whether is iterator is the same as IncidentIterator iit
        /* Same IncidentIterator points to the same incident edge of the same node
        /* in the same graph
        */
        bool operator==(const IncidentIterator& iit) const {
            return iit.graph_ == graph_ and iit.node_index() == node_id_ and iit.index() == node_ie_id;
        }
    private:
        friend class Graph;
        friend class Node;
        // HW1 #3: YOUR CODE HERE
        graph_type* graph_;
        /** index of the parent node in the graph */
        size_type node_id_;
        /** index of the incident edge of the parent node */
        size_type node_ie_id;
        /** Constructor of a valid IncidentIterator*/
        IncidentIterator(const graph_type * graph, size_type node_index, size_type index)
                : graph_(const_cast<Graph*>(graph)), node_id_(node_index), node_ie_id(index){}
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
        // Supply definitions AND SPECIFICATIONS for:

        /** Return the edge this iterator is pointing to */
        Edge operator*() const {
            return graph_->edge(id_);
        }
         EdgeIterator& operator++() {
             ++id_;
             return *this;
         }

        /** Return the graph of this iterator*/
        const graph_type* graph() const {
            return graph_;
        }

        /** Return the index of edge this iterator is poiting to */
        size_type index() const {
            return id_;
        }

        /** Test if this iterator and EdgeIterator ei are the same
        /*  the same EdgeIterators points to the same edge in the same graph
        */
         bool operator==(const EdgeIterator& ei) const {
            return (ei.graph() == graph_) and (ei.index() == id_);
         }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        graph_type* graph_;
        size_type id_;
        bool isvalid_;
        /** Constructor of a valid EdgeIterator */
        EdgeIterator(const graph_type * graph, size_type index):
                graph_(const_cast<Graph*>(graph)), id_(index), isvalid_(true){}

    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return an iterator pointing to the first edge in this graph */
     edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }
    /** Return an iterator pointing to the last edge in this graph */
     edge_iterator edge_end() const {
         return EdgeIterator(this,num_edges_);
     }

private:
    std::vector<Node> nodes_;
    std::vector<Point> node_elements_;
    std::vector<std::vector<size_type>> incident_edges_;
    std::vector<std::vector<size_type>> edge_elements_;
    size_type num_nodes_;
    size_type num_edges_;

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
