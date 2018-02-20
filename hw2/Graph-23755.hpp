#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <vector>

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

    /** Type of node value. */
    using node_value_type = V;

    /** Type of edge value. */
    using edge_value_type = E;

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

    /** Type of a pair of indexes. */
    using pair_type = std::pair<size_type, size_type>;

    /** A simple hash function for the pair type. */
    struct pair_hash {
        inline std::size_t operator()(pair_type v) const {
            return v.first * 0xffffffff + v.second;
        }
    };

    /** A simple hash function for the size type. */
    struct size_hash {
        inline std::size_t operator()(size_type v) const {
            return v;
        }
    };

    /** Type of edge map, which uses a hash map such that access the pair of
     * (a,b) is O(1). */
    using edge_map_type = std::unordered_map<pair_type, size_type, pair_hash>;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
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
        }

        /** Return this node's position (setter). */
        Point &position() {
            return graph_->node_data[node_id].p;
        }

        /** Return this node's position (getter). */
        const Point &position() const {
            return graph_->node_data[node_id].p;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return graph_->node_data[node_id].idx;
        }

        /** Return this node's unique ID. */
        size_type uid() const {
            return node_id;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have theWhen does one halfspace contain another and the same index.
         */
        bool operator==(const Node &n) const {
            return index() == n.index();
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node &n) const {
            return index() < n.index();
        }


        /** Return this node's value (setter). */
        node_value_type &value() {
            return graph_->node_data[node_id].v;
        }

        /** Return this node's value (getter). */
        const node_value_type &value() const {
            return value();
        }

        /** Return this node's degree. */
        size_type degree() const {
            return graph_->node_data[node_id].ie.size();
        }

        /** Return this node's starting iterator for incident edge. */
        IncidentIterator edge_begin() const {
            return IncidentIterator(graph_, node_id, 0);
        }

        /** Return this node's ending iterator for incident edge. */
        IncidentIterator edge_end() const {
            return IncidentIterator(graph_, node_id, degree());
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects

        // A pointer to graph and the index of node
        graph_type *graph_;
        // Unique ID
        size_type node_id;

        // Private constructor to create a node
        Node(const graph_type *graph, size_type id)
                : graph_(const_cast<graph_type *>(graph)), node_id(id) {
        }

    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return i2u_n.size();
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
    Node add_node(const Point &position,
                  const node_value_type &value = node_value_type()) {

        // Construct and add node information to STL container
        Node_info new_node;
        new_node.p = position;
        new_node.v = value;
        new_node.idx = num_nodes();
        node_data.push_back(new_node);

        // Add index
        i2u_n.push_back(num_nodes());

        // Return node
        return Node(this, num_nodes() - 1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        return (this == n.graph_) && (n.index() < num_nodes());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // Check index is in range
        assert(i < num_nodes());
        // Convert index to unique id and return
        return Node(this, i2u_n[i]);
    }

    /** Return starting iterator for node. */
    NodeIterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /** Return ending iterator for node. */
    NodeIterator node_end() const {
        return NodeIterator(this, num_nodes());
    }

    /** Remove the node @a n from graph.
     * @param[in]   n       @a n's unique ID is @a uid
     *                      @a n's index is @a idx
     * @post        new node index[i] = old node index[i] -1, for all i > uid
     *              new i2u[i] = old i2u[i+1], for all i >index
     *              new num_nodes() = old num_nodes() - 1
     *              all adjacent edges connected to @a n are removed
     *
     * @return 1 if @a n is currently a Node of this Graph
     *         0 if @a n is currently not a Node of this Graph
     *
     * Complexity: O(num_nodes()+num_edges()) assuming sparsity.
     */
    size_type remove_node(const Node &n) {
        if (has_node(n)) {
            // Remove incident edges
            auto v = node_data[n.node_id].ie;
            for (unsigned i = 0; i < v.size(); i++) {
                remove_edge(Edge(this, v[i]));
            }

            // Update index
            for (unsigned i = n.node_id + 1; i < node_data.size(); i++) {
                node_data[i].idx--;
            }

            // Update mapping
            for (unsigned i = n.index(); i < i2u_n.size() - 1; i++) {
                i2u_n[i] = i2u_n[i + 1];
            }
            i2u_n.pop_back();
            return 1;
        }
        return 0;
    }

    /** Remove the node at this node iterator*/
    node_iterator remove_node(node_iterator n_it) {
        remove_node(*n_it);
        return (n_it);
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
        /** Construct an invalid Edge. */
        Edge() {
        }

        /** Return a node of this Edge
         *
         *  Complexity: O(1).
         * */
        Node node1() const {
            // Get the first node index from STL container using pointer to
            // graph and return node
            return Node(graph_, graph_->edge_data[edge_id].pair.first);
        }

        /** Return the other node of this Edge
         *
         *  Complexity: O(1).
         * */
        Node node2() const {
            // Get the second node index from STL container using pointer to
            // graph and return node
            return Node(graph_, graph_->edge_data[edge_id].pair.second);
        }

        /** Return index of this Edge
         *
         *  Complexity: O(1).
         */
        size_type index() const {
            return graph_->edge_data[edge_id].idx;
        }

        /** Return length of this Edge
         *
         *  Complexity: O(1).
         */
        double length() const {
            return norm(node1().position() - node2().position());
        }

        /** Return this edge's value (setter). */
        edge_value_type &value() {
            return graph_->edge_data[edge_id].v;
        }

        /** Return this edge's value (getter). */
        const edge_value_type &value() const {
            return value();
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            return index() == e.index();
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            return index() < e.index();
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects

        // A pointer to graph
        graph_type *graph_;
        // Unique ID
        size_type edge_id;

        // Private constructor to create an edge
        Edge(const graph_type *graph, size_type id)
                : graph_(const_cast<graph_type *>(graph)), edge_id(id) {
        }

    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     *
     * Complexity: O(1).
     */
    size_type num_edges() const {
        return i2u_e.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     *
     * Complexity: O(1).
     */
    Edge edge(size_type i) const {
        // Check the index of edge is in the range
        assert(i < num_edges());
        return Edge(this, i2u_e[i]);

    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     *
     * Complexity: O(1).
     */
    bool has_edge(const Node &a, const Node &b) const {
        // Check nodes are in graph
        assert(has_node(a));
        assert(has_node(b));

        return (edge_map.count(make_pair(a, b)) ||
                edge_map.count(make_pair(b, a)));
    }

    /** Create a pair of indexes given two node */
    pair_type make_pair(const Node &a, const Node &b) const {
        return std::make_pair(a.node_id, b.node_id);
    }

    /** Search for the index of edge given nodes, assuming the edge exists
     * (unique ID actually)
     *
     *  Complexity: O(1).
     */
    size_type get_edge_index(const Node &a, const Node &b) const {
        if (edge_map.count(make_pair(a, b))) {
            return edge_map.find(make_pair(a, b))->second;
        } else {
            return edge_map.find(make_pair(b, a))->second;
        }
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
     *
     * Complexity: O(1).
     */
    Edge add_edge(const Node &a, const Node &b) {
        // Check nodes are in graph
        assert(has_node(a));
        assert(has_node(b));

        // Return the edge if it exists
        if (has_edge(a, b)) {
            return Edge(this, get_edge_index(a, b));
        }

        // Create a index pair of nodes and add to map
        pair_type new_edge_pair = make_pair(a, b);
        edge_map.insert(std::make_pair(new_edge_pair, num_edges()));

        // Add to adjacency list
        this->node_data[a.index()].ie.push_back(num_edges());
        this->node_data[b.index()].ie.push_back(num_edges());

        // Construct edge information and add to STL container
        Edge_info new_edge;
        new_edge.pair = new_edge_pair;
        new_edge.idx = num_edges();
        edge_data.push_back(new_edge);

        // Update index
        i2u_e.push_back(num_edges());

        // Return the new edge
        return Edge(this, num_edges() - 1);
    }

    /** Return starting iterator of edge. */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** Return ending iterator of edge. */
    edge_iterator edge_end() const {
        return EdgeIterator(this, num_edges());
    }

    /** Remove the edge @a e from graph.
     * @param[in]   e       @a e's unique ID is @a uid
     *                      @a e's index is @a idx
     * @post        new edge index[i] = old edge index[i] -1, for all i > uid
     *              new i2u[i] = old i2u[i+1], for all i >index
     *              new num_edges() = old num_edges() - 1
     *
     * @return 1 if @a e is currently a edge of this Graph
     *         0 if @a e is currently not a edge of this Graph
     *
     * Complexity: O(num_edges()).
     */
    size_type remove_edge(const Edge &e) {
        if (has_edge(e.node1(), e.node2())) {
            // Update map
            pair_type pair = make_pair(e.node1(), e.node2());
            edge_map.erase(pair);
            edge_map.erase(pair);

            // Remove from adjcency list
            auto v1 = node_data[e.node1().node_id].ie;
            auto it1 = std::find(v1.begin(), v1.end(), e.edge_id);
            v1.erase(it1);
            node_data[e.node1().node_id].ie = v1;

            auto v2 = node_data[e.node2().node_id].ie;
            auto it2 = std::find(v2.begin(), v2.end(), e.edge_id);
            v2.erase(it2);
            node_data[e.node2().node_id].ie = v2;

            // Update index
            for (unsigned i = e.edge_id + 1; i < edge_data.size(); i++) {
                edge_data[i].idx--;
            }

            // Update unique id
            for (unsigned i = e.index(); i < i2u_e.size() - 1; i++) {
                i2u_e[i] = i2u_e[i + 1];
            }
            i2u_e.pop_back();
            return 1;
        }
        return 0;
    }

    /** Remove the edge connecting node @a a and @a b. */
    size_type remove_edge(const Node &a, const Node &b) {
        if (has_edge(a, b)) {
            return remove_edge(Edge(this, get_edge_index(a, b)));
        }
        return 0;
    }

    /** Remove the edge at this edge iterator. */
    edge_iterator remove_edge(edge_iterator e_it) {
        remove_edge(*e_it);
        return e_it;
    }

    /** Reset the value of nodes. */
    void reset_value(V value) {
        for (auto start = node_data.begin();
             start != node_data.end(); ++start) {
            start->v = value;
        }
    }

    /** Set the position of nodes by unique ID. */
    void set_pos(size_type uid, Point point) {
        node_data[uid].p = point;
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Clear all STL containers
        edge_map.clear();
        node_data.clear();
        edge_data.clear();
        i2u_e.clear();
        i2u_n.clear();
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
        using pointer           = Node *;                    // Pointers to elements
        using reference         = Node &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }

        /** Return current node. */
        Node operator*() const {
            return Node(g, g->i2u_n[i]);
        }

        /** Return iterator to previous node. */
        NodeIterator &operator--() {
            i--;
            return *this;
        }

        /** Return iterator to next node. */
        NodeIterator &operator++() {
            i++;
            return *this;
        }

        /** Compare two iterators of node. */
        bool operator==(const NodeIterator &ni) const {
            // Check graph and node equality
            return (g == ni.g) && (i == ni.i);
        }

    private:
        friend class Graph;

        // Store a pointer to graph and index
        graph_type *g;
        size_type i;

        // Private constructor
        NodeIterator(const graph_type *graph, size_type id)
                : g(const_cast<graph_type *>(graph)), i(id) {}

        friend class Node;
    };

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator : private totally_ordered<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {
        }

        /** Return current edge. */
        Edge operator*() const {
            // Get edge id
            size_type edge_id = g->node_data[node1_id].ie[incid_edge];
            return Edge(g, edge_id);
        }

        /** Return iterator of next incident edge. */
        IncidentIterator &operator++() {
            incid_edge++;
            return *this;
        }

        /** Compare two iterators of incident edge. */
        bool operator==(const IncidentIterator &ii) const {
            // Check graph, current node and incident edge equality
            return (g == ii.g) &&
                   (node1_id == ii.node1_id) &&
                   (incid_edge == ii.incid_edge);
        }

    private:
        friend class Graph;

        // Store pointer to graph, current node id and current incident edge
        // position
        graph_type *g;
        size_type node1_id;
        size_type incid_edge;

        // Private constructor
        IncidentIterator(const graph_type *graph, size_type id, size_type len)
                : g(const_cast<graph_type *>(graph)),
                  node1_id(id), incid_edge(len) {}

        friend class Node;
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
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {
        }

        /** Return current edge. */
        Edge operator*() const {
            return Edge(g, g->i2u_e[i]);
        }

        /** Return iterator of previous edge. */
        EdgeIterator &operator--() {
            i--;
            return *this;
        }

        /** Return iterator of next edge. */
        EdgeIterator &operator++() {
            i++;
            return *this;
        }

        /** Compare two iterators of edge. */
        bool operator==(const EdgeIterator &ei) const {
            // Check graph and edge equality
            return (g == ei.g) && (i == ei.i);
        }

    private:
        friend class Graph;

        // Store pointer to graph and edge index
        graph_type *g;
        size_type i;

        // Private constructor
        EdgeIterator(const graph_type *graph, size_type id)
                : g(const_cast<graph_type *>(graph)), i(id) {}

        friend class Edge;
    };

private:

    // STL container to store position of nodes and edges

    struct Node_info {
        Point p;
        node_value_type v;
        std::vector <size_type> ie;
        size_type idx;
    };

    struct Edge_info {
        pair_type pair;
        edge_value_type v;
        size_type idx;
    };

    edge_map_type edge_map;
    std::vector <Node_info> node_data;
    std::vector <Edge_info> edge_data;

    std::vector <size_type> i2u_n;
    std::vector <size_type> i2u_e;
};

#endif // CME212_GRAPH_HPP