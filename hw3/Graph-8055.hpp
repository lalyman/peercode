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

template<typename V, typename E>
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
    using edge_value_type = E;

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
        nodes_neighbors = new std::unordered_map<size_type, std::unordered_map<size_type, size_type>>;
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
        Node(const Graph *g = nullptr) :
                graph(g) {
        }

        /** Return this node's position. */
        const Point &position() const {
            return graph->nodes->at(uid).position;
        }

        /** Return this node's position, mutable. */
        Point &position() {
            return graph->nodes->at(uid).position;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return uid;
        }

        /** Return this node's degree, number of connection of this node, a number in the range [0, graph_size-1). */
        size_type degree() const {
            return (*graph->nodes_neighbors)[uid].size();
        }

        /**Return an iterator to iterate through the connections of the node. */
        IncidentIterator edge_begin() const {
            IncidentIterator inc_it(uid, graph);
            inc_it.intern_it = (*graph->nodes_neighbors)[uid].begin();

            return inc_it;
        }

        /** Return the end iterator of the connections of the node. */
        IncidentIterator edge_end() const {
            IncidentIterator inc_it(uid, graph);
            inc_it.intern_it = (*graph->nodes_neighbors)[uid].end();

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
        Node(const Graph *g, size_type uid) :
                graph(g), uid(uid) {
        }

        // Allow Graph to access Node's private member data and functions.
        friend class Graph<node_value_type, edge_value_type>;

        const Graph *graph;
        size_type uid;

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

        return Node(this, new_uid);
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
        return Node(this, i);
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
        Edge(const Graph *g = nullptr) :
                graph(g) {
        }

        /** Return a node of this Edge */
        Node node1() const {
            return graph->node(n1);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return graph->node(n2);
        }

        double length() const {
            return norm(node1().position() - node2().position());
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            return ((e.node1() == node1() && e.node2() == node2())
                    || (e.node2() == node1() && e.node1() == node2()))
                   && (graph == e.graph);
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            if (graph != e.graph)
                return graph < e.graph;

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

        /** Return the value stored into an edge, of edge_value_type type. */
        edge_value_type &value() {
            return graph->edges->at(graph->nodes_neighbors->at(node1().uid).at(node2().uid)).value;
        }

        /** Return the value stored into an edge, of edge_value_type type. */
        const edge_value_type &value() const {
            return graph->edges->at(graph->nodes_neighbors->at(node1().uid).at(node2().uid)).value;
        }

    private:
        Edge(const Graph *g, size_type uid, size_type n1, size_type n2) :
                graph(g), uid(uid), n1(n1), n2(n2) {
        }

        // Allow Graph to access Edge's private member data and functions.
        friend class Graph<node_value_type, edge_value_type>;

        const Graph *graph;
        size_type uid;
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
        intern_element_edge e = edges->at(i);
        return Edge(this, i, e.n1, e.n2);
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
        Edge edge(this);
        edge.n1 = a.uid;
        edge.n2 = b.uid;

        if (has_edge(a, b))
            edge.uid = nodes_neighbors->at(a.uid).at(b.uid);
        else {
            intern_element_edge e;

            size_type new_edge_id = edges->size();
            e.n1 = a.uid;
            e.n2 = b.uid;

            edges->insert({new_edge_id, e});
            (*nodes_neighbors)[a.uid][b.uid] = new_edge_id;
            (*nodes_neighbors)[b.uid][a.uid] = new_edge_id;
        }

        return edge;
    }

    /** Remove node n from the graph and all edges linked to it.
     *
     *  The method is very light as it only swaps node n and the last one and removes the last one from the
     *  unordered map. The swap is done using std::swap (which calls std::move). It also removes all the edges between
     *  node n and its neighbors by deleting its adjacency map and removing this node from the adjacency maps of its
     *  neighbors.
     *
     *  Complexity: Due to the relative small number of neighbors compared to the number of nodes and edges, we can
     *  consider that this function is in O(1).
     *
     * @param node The node to be removed
     * @return 0 if the node was not deleted, 1 if it was.
     * @post node(i).index() == i for all i with 0<= i < num_nodes()
     * @post node(n.index()) = n
     */

    size_type remove_node(const Node &node) {
        //switch d'abord les last_node et celui qu'on veut supprimer
        if (node.uid >= size())
            return 0;

        Node last_node = this->node(nodes->size() - 1);
        if (node.uid != last_node.uid) {
            //swap des edges de last vers n
            if (node.degree() > 0) {
                for (auto e = node.edge_begin(); e != node.edge_end(); ++e) {
                    size_type edge_id = (*e).uid;
                    size_type other;
                    if (edges->at(edge_id).n1 == node.uid) {
                        edges->at(edge_id).n1 = last_node.uid;
                        other = edges->at(edge_id).n2;
                    } else {
                        edges->at(edge_id).n2 = last_node.uid;
                        other = edges->at(edge_id).n1;
                    }
                    nodes_neighbors->at(other).erase(node.uid);
                    nodes_neighbors->at(other)[last_node.uid] = edge_id;
                }
            }

            if (last_node.degree() > 0) {
                for (auto e = last_node.edge_begin(); e != last_node.edge_end(); ++e) {
                    size_type edge_id = (*e).uid;
                    size_type other;
                    if (edges->at(edge_id).n1 == last_node.uid) {
                        edges->at(edge_id).n1 = node.uid;
                        other = edges->at(edge_id).n2;
                    } else {
                        edges->at(edge_id).n2 = node.uid;
                        other = edges->at(edge_id).n1;
                    }

                    nodes_neighbors->at(other).erase(last_node.uid);
                    nodes_neighbors->at(other)[node.uid] = edge_id;
                }
            }

            //swap des nodes
            std::swap(nodes->at(node.uid), nodes->at(last_node.uid));
            //swap des neighbors map
            std::swap((*nodes_neighbors)[node.uid], (*nodes_neighbors)[last_node.uid]);
        }
        if (last_node.degree() > 0) {
            // delete edges from this node
            for (auto e = last_node.edge_begin(); e != last_node.edge_end(); ++e) {
                Node neighbor = (*e).node2();

                if ((*e).uid != num_edges() - 1) {

                    intern_element_edge last_edge = (*edges)[num_edges() - 1];

                    //update id of edge swap in NN
                    nodes_neighbors->at(last_edge.n1).at(last_edge.n2) = (*e).uid;
                    nodes_neighbors->at(last_edge.n2).at(last_edge.n1) = (*e).uid;

                    //swap last edge and this one
                    std::swap(edges->at((*e).uid), edges->at(num_edges() - 1));
                }

                //destroy last one in edges
                edges->erase(num_edges() - 1);

                //erase it from node_neighbors of neighbor
                nodes_neighbors->at(neighbor.uid).erase(last_node.uid);
            }
        }

        nodes_neighbors->erase(last_node.uid);
        nodes->erase(last_node.uid);
        return 1;
    }

    /** See remove_node(const Node&) for the specification of the node removal. Calls the former to remove (*n_it).
     *
     * @param n_it Iterator pointing to the node to remove.
     * @return A valid iterator to the node replacing the one just removed.
     */
    node_iterator remove_node(node_iterator n_it) {
        Node node = (*n_it);
        remove_node(node);

        NodeIterator node_it(this);
        node_it.intern_it = nodes->find(node.uid);

        return node_it;
    }

    /** Deletes the edge linking nodes @a n1 and @a n2. To avoid hole in the edges index, we swap it with the last
     * edge before and then delete the last one.
     *
     * @param n1 First node of the edge being removed.
     * @param n2 Second of the edge being removed.
     * @return 1 if the edge was removed successfully, 0 if not (if it does not exist for example).
     */
    size_type remove_edge(const Node &n1, const Node &n2) {
        if (has_edge(n1, n2)) {
            size_type edge_id = nodes_neighbors->at(n1.uid).at(n2.uid);
            intern_element_edge &to_erase_edge = (*edges)[edge_id];
            intern_element_edge &to_update_edge = (*edges)[edges->size() - 1];

            //update id of edge swap in NN
            nodes_neighbors->at(to_update_edge.n2).at(to_update_edge.n1) = edge_id;
            nodes_neighbors->at(to_update_edge.n1).at(to_update_edge.n2) = edge_id;

            nodes_neighbors->at(to_erase_edge.n1).erase(to_erase_edge.n2);
            nodes_neighbors->at(to_erase_edge.n2).erase(to_erase_edge.n1);

            std::swap(to_erase_edge, to_update_edge);

            //destroy last one in edges
            edges->erase(edges->size() - 1);

            return 1;
        } else
            return 0;
    }

    /**
     * Remove the edge @a edge. It calls the function remove_edge(const Node&, const Node&).
     * @param edge The edge to be removed.
     * @return
     */
    size_type remove_edge(const Edge &edge) {
        return remove_edge(edge.node1(), edge.node2());
    }

    /**
     * Remove the edge the edge iterator points to. Returns then an iterator on the element that took the place of the
     * edge just deleted.
     * @param e_it Iterator pointing to the edge we want to delete.
     * @return A valid iterator pointing to the new element at the id (*e_it).uid.
     */
    edge_iterator remove_edge(edge_iterator e_it) {
        Edge edge = (*e_it);
        remove_edge(edge);

        EdgeIterator edge_it(this);
        edge_it.intern_it = edges->find(edge.uid);

        return edge_it;
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
            return Node(graph, (*intern_it).first);
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return NodeIterator to the new current item.
         */
        NodeIterator &operator++() {
            ++intern_it;
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
        friend class Graph<node_value_type, edge_value_type>;

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
    typedef typename std::unordered_map<size_type, size_type>::const_iterator intern_edge_map_iterator;

    class IncidentIterator : private equality_comparable<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator(size_type origin, const Graph *graph = nullptr) :
                graph(graph), origin_node(origin) {
        }

        /** Deferencing operator of an iterator over the edges of graph g starting from a specific node,
         * specified during the choice of the underlying container.
         *
         * @return An Edge with id corresponding to the current iteration.
         */
        Edge operator*() const {
            intern_element_edge e = (*graph->edges)[(*intern_it).second];
            Edge edge(graph);
            if (origin_node == e.n1) {
                edge.n1 = e.n1;
                edge.n2 = e.n2;
            } else {
                edge.n2 = e.n1;
                edge.n1 = e.n2;
            }
            edge.uid = (*intern_it).second;

            return edge;
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return IncidentIterator to the new current item.
         */
        IncidentIterator &operator++() {
            ++intern_it;
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
        friend class Graph<node_value_type, edge_value_type>;

        const Graph *graph;
        intern_edge_map_iterator intern_it;
        size_type origin_node;

    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */

    typedef typename std::unordered_map<size_type, intern_element_edge>::const_iterator intern_edge_vector_iterator;

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
            return Edge(graph, (*intern_it).first, (*intern_it).second.n1, (*intern_it).second.n2);
        }

        /** Prefix ++ operator, advances the iterator to the next item and returns an iterator to the new current item.
         *
         * @return EdgeIterator to the new current item.
         */
        EdgeIterator &operator++() {
            ++intern_it;
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
        friend class Graph<node_value_type, edge_value_type>;

        const Graph *graph;
        intern_edge_vector_iterator intern_it;
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
    std::unordered_map<size_type, std::unordered_map<size_type, size_type>> *nodes_neighbors;


    struct intern_element_node {
        Point position;
        node_value_type value;
    };
    struct intern_element_edge {
        size_type n1;
        size_type n2;
        edge_value_type value;
    };

};

#endif // CME212_GRAPH_HPP
