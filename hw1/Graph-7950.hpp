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

template<typename V>
class Graph {
private:

    /** Intern node and edge types of this graph. */
    struct intern_element_node;
    struct intern_element_edge;

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph;
    using node_value_type = V;

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

    /** Construct an empty graph. */
    Graph() {
        nodes = new std::unordered_map<size_type, intern_element_node>;
        edges = new std::unordered_map<size_type, intern_element_edge>;
        nodes_neighbors = new std::unordered_map<size_type, std::unordered_map<size_type, intern_element_edge>>;
    }


    /** Default destructor */
    ~Graph() {
        delete nodes;
        delete edges;
        delete nodes_neighbors;
    };

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
    // NODES
    //

    /** @class Graph::Node
     * @brief Class representing the graph's nodes.
     *
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private equality_comparable<Node> {
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
        Node(const Graph *g = nullptr) {
            graph = g;
        }

        /** Return this node's position. */
        const Point &position() const {
            return graph->nodes->at(uid).position;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return uid;
        }

        /** Return this node's degree, number of connection of this node, a number in the range [0, graph_size-1). */
        size_type degree() const {
            if (graph->nodes_neighbors->find(uid) != graph->nodes_neighbors->end())
                return graph->nodes_neighbors->at(uid).size();
            else
                return 0;
        }

        /**Return an iterator to iterate through the connections of the node. */
        IncidentIterator edge_begin() const {
            IncidentIterator inc_it(graph);
            if (graph->nodes_neighbors->find(uid) != graph->nodes_neighbors->end())
                inc_it.intern_it = graph->nodes_neighbors->at(uid).begin();
            else
                throw std::out_of_range("Node not connected to anyone");
            return inc_it;
        }

        /** Return the end iterator of the connections of the node. */
        IncidentIterator edge_end() const {
            IncidentIterator inc_it(graph);
            if (graph->nodes_neighbors->find(uid) != graph->nodes_neighbors->end())
                inc_it.intern_it = graph->nodes_neighbors->at(uid).end();
            else
                throw std::out_of_range("Node not connected to anyone");
            return inc_it;
        }

        /** Equality operator with another node. */
        bool operator==(const Node &n) const {
            return (n.graph == graph) && (n.uid == uid);
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
            return uid < n.uid;
        }

        /** Return the value stored into a node, of node_value_type type. */
        node_value_type &value() {
            return (*graph->nodes)[uid].value;
        }

        /** Return the value stored into a node, of node_value_type type. */
        const node_value_type &value() const {
            return graph->nodes->at(uid).value;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph<node_value_type>;

        size_type uid;
        const Graph *graph;


    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return nodes->size();
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
    Node add_node(const Point &position, const node_value_type &value = node_value_type()) {
        size_type new_uid = nodes->size();

        intern_element_node new_elem;
        new_elem.position = position;
        new_elem.value = value;
        nodes->insert({new_uid, new_elem});

        Node node(this);
        node.uid = new_uid;

        return node;
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        return n.uid < num_nodes();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        Node node(this);
        node.uid = i;

        return node;
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
    class Edge : private equality_comparable<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge(const Graph *g = nullptr) {
            graph = g;
        }

        /** Return a node of this Edge */
        Node node1() const {
            return graph->node(n1);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return graph->node(n2);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            return (e.node1() == node1() && e.node2() == node2()) || (e.node2() == node1() && e.node1() == node2());
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            Node min1, min2, max1, max2;

            if (node1() < node2()) {
                min1 = node1();
                max1 = node2();
            } else {
                min1 = node2();
                max1 = node1();
            }

            if (e.node1() < e.node2()) {
                min2 = e.node1();
                max2 = e.node2();
            } else {
                min2 = e.node2();
                max2 = e.node1();
            }


            return (min1 < min2) || (min1 == min2 && max1 < max2);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph<node_value_type>;

        const Graph *graph;
        size_type n1;
        size_type n2;
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return edges->size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        Edge edge(this);
        intern_element_edge e = edges->at(i);
        edge.n1 = e.n1;
        edge.n2 = e.n2;

        return edge;
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node &a, const Node &b) const {
        try {
            nodes_neighbors->at(a.uid).at(b.uid);
        }
        catch (const std::out_of_range &e) {
            return false;
        }

        return true;
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
    Edge add_edge(const Node &a, const Node &b) {
        intern_element_edge e;
        if (has_edge(a, b))
            e = nodes_neighbors->at(a.uid).at(b.uid);
        else {
            size_type new_edge_id = edges->size();
            e.n1 = a.uid;
            e.n2 = b.uid;

            edges->insert({new_edge_id, e});
            (*nodes_neighbors)[a.uid][b.uid] = (*nodes_neighbors)[b.uid][a.uid] = e;
        }

        Edge edge(this);
        edge.n1 = e.n1;
        edge.n2 = e.n2;
        return edge;
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes->clear();
        edges->clear();
        nodes_neighbors->clear();
    }

    //
    // Node Iterator
    //
    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    typedef typename std::unordered_map<size_type, intern_element_node>::const_iterator intern_node_iterator;

    class NodeIterator : private equality_comparable<NodeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node *;                    // Pointers to elements
        using reference         = Node &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


        /** Construct an invalid NodeIterator. */
        NodeIterator(const Graph *graph = nullptr) :
                graph(graph) {
        }

        /** Deferencing operator of an iterator over the nodes of graph g.
         *
         * @return A Node with id corresponding to the current iteration.
         */
        Node operator*() const {
            Node node(graph);
            node.uid = (*intern_it).first;

            return node;
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return NodeIterator to the new current item.
         */
        NodeIterator &operator++() {
            intern_it++;
            return (*this);
        }

        /**
         * Equality operator.
         * @param other
         * @return Returns true if @a other points to the same item as this iterator; otherwise returns false.
         */
        bool operator==(const NodeIterator &other) const {
            return intern_it == other.intern_it;
        }

    private:
        friend class Graph<node_value_type>;

        const Graph *graph;
        intern_node_iterator intern_it;

    };

    /**
     * Iterator to the first node in the graph.
     * @return A NodeIterator pointing to the first node in the graph.
     */
    NodeIterator node_begin() const {
        NodeIterator node_it(this);
        node_it.intern_it = nodes->begin();
        return node_it;
    }

    /**
     * Iterator to the end of the containers of the nodes.
     * @return A NodeIterator pointing to the end of the container.
     */
    NodeIterator node_end() const {
        NodeIterator node_it(this);
        node_it.intern_it = nodes->end();
        return node_it;
    }

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    typedef typename std::unordered_map<size_type, intern_element_edge>::const_iterator intern_edge_iterator;

    class IncidentIterator : private equality_comparable<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator(const Graph *graph = nullptr) :
                graph(graph) {
        }

        /** Deferencing operator of an iterator over the edges of graph g starting from a specific node,
         * specified during the choice of the underlying container.
         *
         * @return An Edge with id corresponding to the current iteration.
         */
        Edge operator*() const {
            Edge edge(graph);
            edge.n1 = (*intern_it).second.n1;
            edge.n2 = (*intern_it).second.n2;

            return edge;
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return IncidentIterator to the new current item.
         */
        IncidentIterator &operator++() {
            intern_it++;
            return (*this);
        }

        /**
         * Equality operator.
         * @param other
         * @return Returns true if @a other points to the same item as this iterator; otherwise returns false.
         */
        bool operator==(const IncidentIterator &other) const {
            return intern_it == other.intern_it;
        }

    private:
        friend class Graph<node_value_type>;

        const Graph *graph;
        intern_edge_iterator intern_it;

    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private equality_comparable<EdgeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator(const Graph *graph = nullptr) :
                graph(graph) {
        }

        /** Deferencing operator of an iterator over the edges of graph g.
         *
         * @return An Edge with id corresponding to the current iteration.
         */
        Edge operator*() const {
            Edge edge(graph);
            edge.n1 = (*intern_it).second.n1;
            edge.n2 = (*intern_it).second.n2;

            return edge;
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return EdgeIterator to the new current item.
         */
        EdgeIterator &operator++() {
            intern_it++;
            return (*this);
        }

        /**
         * Equality operator.
         * @param other
         * @return Returns true if @a other points to the same item as this iterator; otherwise returns false.
         */
        bool operator==(const EdgeIterator &other) const {
            return intern_it == other.intern_it;
        }


    private:
        friend class Graph<node_value_type>;

        const Graph *graph;
        intern_edge_iterator intern_it;
    };

    /**
     * Iterator to the first edge in the graph.
     * @return A EdgeIterator pointing to the first edge in the graph.
     */
    edge_iterator edge_begin() const {
        EdgeIterator edge_it(this);
        edge_it.intern_it = edges->begin();
        return edge_it;
    }

    /**
     * Iterator to the end of the containers of the edges.
     * @return A EdgeIterator pointing to the end of the container.
     */
    edge_iterator edge_end() const {
        EdgeIterator edge_it(this);
        edge_it.intern_it = edges->end();
        return edge_it;
    }


private:

    /** Intern containers storing nodes and edges. */
    std::unordered_map<size_type, intern_element_node> *nodes;
    std::unordered_map<size_type, intern_element_edge> *edges;
    std::unordered_map<size_type, std::unordered_map<size_type, intern_element_edge>> *nodes_neighbors;


    struct intern_element_node {
        Point position;
        node_value_type value;
    };
    struct intern_element_edge {
        size_type n1;
        size_type n2;
    };

};

#endif // CME212_GRAPH_HPP
