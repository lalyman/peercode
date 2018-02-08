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
template<typename V>
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

    /** Type of incident map. */
    using incid_map_type = std::unordered_map <size_type,
    std::vector<size_type>, size_hash>;

    /** Type of edge vector. */
    using edge_vec_type = std::vector<pair_type>;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        // Initialize the number of node and edge
        n_nodes = 0;
        n_edges = 0;

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

        /** Return this node's position. */
        const Point &position() const {
            return graph_->pos_data[node_id];
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return node_id;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have theWhen does one halfspace contain another and the same index.
         */
        bool operator==(const Node &n) const {
            // Check the graph owning two nodes and then indexes
            return (graph_ == n.graph_) && (index() == n.index());
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
            // First check the graph equality
            assert(graph_ == n.graph_);
            return (index() < n.index());
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;

        /** Return this node's value (setter). */
        node_value_type &value() {
            return graph_->value_data[node_id];
        }

        /** Return this node's value (getter). */
        const node_value_type &value() const {
            return value();
        }

        /** Return this node's degree. */
        size_type degree() const {
            return graph_->incid_map[node_id].size();
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
        return n_nodes;
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

        // Add position and value to STL container
        pos_data.push_back(position);
        value_data.push_back(value);

        // Increment node number
        n_nodes++;

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
        assert(i < num_nodes());
        return Node(this, i);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const

    /** Return starting iterator for node. */
    NodeIterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /** Return ending iterator for node. */
    NodeIterator node_end() const {
        return NodeIterator(this, num_nodes());
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
            return Node(graph_, graph_->edge_vec[edge_id].first);
        }

        /** Return the other node of this Edge
         *
         *  Complexity: O(1).
         * */
        Node node2() const {
            // Get the second node index from STL container using pointer to
            // graph and return node
            return Node(graph_, graph_->edge_vec[edge_id].second);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            // First check the graph equality and then edge index
            return (graph_ == e.graph_) && (index() == e.index());
        }

        /** Return index of this Edge
         *
         *  Complexity: O(1).
         */
        size_type index() const {
            return edge_id;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            // First check graph equality
            assert(graph_ == e.graph_);
            return index() < e.index();
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects

        // A pointer to graph and the index of edge
        graph_type *graph_;
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
        return n_edges;
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
        return Edge(this, i);

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
        return std::make_pair(a.index(), b.index());
    }

    /** Search for the index of edge given nodes, assuming the edge exists
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

        // Create a index pair of nodes and add to STL containers
        pair_type new_edge = make_pair(a, b);
        edge_map.insert(std::make_pair(new_edge, n_edges));
        edge_vec.push_back(new_edge);

        // Update incident map
        incid_map[a.index()].push_back(num_edges());
        incid_map[b.index()].push_back(num_edges());

        // Increment number of edge
        n_edges++;

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

    /** Reset the value of nodes. */
    void reset_value(V v) {
        std::fill(value_data.begin(), value_data.end(), v);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Set number of nodes and edges to be 0 and clear all STL containers
        n_nodes = 0;
        n_edges = 0;
        pos_data.clear();
        value_data.clear();
        edge_map.clear();
        edge_vec.clear();
        incid_map.clear();
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

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Node operator*() const
        // NodeIterator& operator++()
        // bool operator==(const NodeIterator&) const

        /** Return current node. */
        Node operator*() const {
            return Node(g, i);
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

        // HW1 #2: YOUR CODE HERE
        // Store a pointer to graph and node id
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

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const

        /** Return current edge. */
        Edge operator*() const {
            // Get edge id
            size_type edge_id=g->incid_map[node1_id][incid_edge];
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

        // HW1 #3: YOUR CODE HERE
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

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // EdgeIterator& operator++()
        // bool operator==(const EdgeIterator&) const

        /** Return current edge. */
        Edge operator*() const {
            return Edge(g, i);
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

        // HW1 #5: YOUR CODE HERE
        // Store pointer to graph and edge id
        graph_type *g;
        size_type i;

        // Private constructor
        EdgeIterator(const graph_type *graph, size_type id)
                : g(const_cast<graph_type *>(graph)), i(id) {}

        friend class Edge;
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const

private:
    // Track number of nodes and edges
    size_type n_nodes, n_edges;

    // STL container to store position of nodes and edges
    std::vector <Point> pos_data;
    std::vector <node_value_type> value_data;
    edge_map_type edge_map;
    edge_vec_type edge_vec;
    incid_map_type incid_map;
};

#endif // CME212_GRAPH_HPP