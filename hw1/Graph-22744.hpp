#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include <unordered_map>
#include <utility>

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

    /** Synonym for the Value of a node */
    using node_value_type = V;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    /** Predeclaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;

    /** Type of node iterators, which iterate over all graph nodes. */
    class NodeIterator;
    /** Synonym for NodeIterator */
    using node_iterator = NodeIterator;

    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;

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
    Graph() = default;

    /** Default destructor */
    ~Graph() = default;




    //
    // NODES
    //
    /** @class Graph::NodeInfo
     * @brief Class representing a node's content (Point and V).
     */
    struct NodeInfo {
        Point p_;
        V v_;
    };

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
            index_ = (size_type) -1;
            graph_ = nullptr;
        }

        /** Return this node's position. */
        const Point &position() const {
            // HW0: YOUR CODE HERE
            return graph_->valueVector[index_].p_;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return index_;
        }

        /** Return this node's graph. */
        const Graph *graph() const {
            // HW0: YOUR CODE HERE
            return graph_;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            // HW0: YOUR CODE HERE
            return (graph() == n.graph()) && (index() == n.index());
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
            // HW0: YOUR CODE HERE
            if (graph_ == n.graph_) {
                return index() < n.index();
            } else {
                return graph_ < n.graph_;
            }
        }

        /**
         *
         * @return const reference to the Node's value.
         */
        const node_value_type &value() const {
            return graph_->valueVector[index_].v_;
        }

        /**
         *
         * @return reference to the Node's value.
         */
        node_value_type &value() {
            return graph_->valueVector[index_].v_;
        }

        /**
         *
         * @return the undirected degree of a node
         */
        size_type degree() const {
            return (size_type) graph_->adjacencyLists[index_].size();
        }

        /**
         *
         * @return an IncidentIterator iterating through a Node's neighbors.
         */
        incident_iterator edge_begin() const {
            return IncidentIterator(this->graph_, *this, this->graph_->adjacencyLists.at(index()).begin());
        }

        /**
         *
         * @return an IncidentIterator pointing to the end of a Node's neighbors.
         */
        incident_iterator edge_end() const {
            return IncidentIterator(this->graph_, *this, this->graph_->adjacencyLists.at(index()).end());
        }


    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // HW0: YOUR CODE HERE
        size_type index_;
        const Graph *graph_;

        Node(const Graph *graph, size_type index) {
            index_ = index;
            graph_ = graph;
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
        return (size_type) valueVector.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position the new node's position
     * @param value the new node's value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point &position, const node_value_type &value = node_value_type()) {
        // HW0: YOUR CODE HERE
        size_type index = size();
        valueVector.push_back({position, value});
        adjacencyLists[index] = std::unordered_map<size_type, size_type>();
        return Node(this, index);
    }

    /**
     * Sets the value of a Node to a provided value.
     * @param node Node whose value has to be changed
     * @param v value to be stored as @a node 's new value
     */
    void set_value(Node &node, const node_value_type &v) {
        valueVector[node.index()].v_ = v;
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     */
    bool has_node(const Node &n) const {
        // HW0: YOUR CODE HERE
        return n.graph_ == this;
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
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
        Edge() {
            // HW0: YOUR CODE HERE
            index_ = (size_type) -1;
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return Node(graph_, myNodeIndex1);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return Node(graph_, myNodeIndex2);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            return (graph() == e.graph()) &&
                   (e.node1().index() == myNodeIndex1) &&
                   (e.node2().index() == myNodeIndex2);
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            // HW0: YOUR CODE HERE
            if (graph_ == e.graph_) {
                return index() < e.index();
            } else {
                return graph_ < e.graph_;
            }
        }

        /** Return this edge's index, a number in the range [0, edge_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return index_;
        }

        /** Return this edge's graph. */
        const Graph *graph() const {
            // HW0: YOUR CODE HERE
            return graph_;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        // HW0: YOUR CODE HERE
        size_type index_;
        const Graph *graph_;
        size_type myNodeIndex1;
        size_type myNodeIndex2;

        Edge(const Graph *graph, size_type index, size_type nodeIndex1, size_type nodeIndex2) {
            graph_ = graph;
            index_ = index;
            myNodeIndex1 = nodeIndex1;
            myNodeIndex2 = nodeIndex2;
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
        return (size_type) edgeMap.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        return *edgeMap.find(i)->second;
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node &a, const Node &b) const {
        return adjacencyLists.at(a.index()).find(b.index()) != adjacencyLists.at(a.index()).end();
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
        // HW0: YOUR CODE HERE
        if (not has_edge(a, b)) {
            size_type newEdgeIndex = num_edges();
            auto *newEdge = new Edge(this, newEdgeIndex, a.index(), b.index());
            edgeMap[newEdgeIndex] = newEdge;
            adjacencyLists[a.index()][b.index()] = newEdgeIndex;
            adjacencyLists[b.index()][a.index()] = newEdgeIndex;
            return *newEdge;
        } else {
            return Edge(this, adjacencyLists[a.index()][b.index()], a.index(), b.index());
        }
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        for (auto edge : edgeMap) {
            delete edge.second;
        }
        edgeMap.clear();
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
        NodeIterator() = default;

        /**
         * Construct a NodeIterator from a pointer to a graph and a start index.
         * @param g pointer to a Graph whose Nodes we want to iterate through
         * @param startIndex the index of the Node starting from which we want to iterate.
         */
        NodeIterator(const Graph *g, size_type startIndex) {
            g_ = g;
            currentIndex = startIndex;
        }

        // HW1 #2: YOUR CODE HERE
        /**
         *
         * @return the Node object currently pointed to by the iterator.
         */
        Node operator*() const {
            return g_->node(currentIndex);
        }

        /**
         *
         * @return an NodeIterator which starts at the next Node in the iterator.
         */
        NodeIterator &operator++() {
            currentIndex++;
            return *this;
        }

        /**
         *
         * @param otherNodeIterator NodeIterator to compare equality with
         * @return true if both NodeIterator objects are pointing to the same Node in the same graph.
         */
        bool operator==(const NodeIterator &otherNodeIterator) const {
            return (g_ == otherNodeIterator.g_) && (currentIndex == otherNodeIterator.currentIndex);
        }

    private:
        friend class Graph;

        // HW1 #2: YOUR CODE HERE
        const Graph *g_;
        size_type currentIndex{};

    };

    // HW1 #2: YOUR CODE HERE
    /**
     *
     * @return the start of a NodeIterator through the nodes of the graph.
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /**
     *
     * @return the end of a NodeIterator through the nodes of the graph.
     */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes());
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
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() = default;


        /**
         * Construct a IncidentIterator from a pointer to a graph, a starting node and an iterator
         * through the starting node's adjacency list.
         * @param g pointer to a graph
         * @param startingNode starting node in the graph @a g
         * @param it iterator through @a startingNode 's adjacency list in @a g
         */
        IncidentIterator(const Graph *g, const Node &startingNode,
                         std::unordered_map<size_type, size_type>::const_iterator it) {
            g_ = g;
            startingNodeId = startingNode.index_;
            currentIt = it;
        }

        // HW1 #3: YOUR CODE HERE
        /**
         *
         * @return the Edge object currently pointed to by the iterator.
         */
        Edge operator*() const {
            return Edge(g_, (*currentIt).second, startingNodeId, (*currentIt).first);
        }

        /**
         *
         * @return an IncidentIterator which starts at the next Edge in the iterator.
         */
        IncidentIterator &operator++() {
            currentIt++;
            return *this;
        }

        /**
         *
         * @param otherIncidentIterator IncidentIterator to compare equality with
         * @return true if both IncidentIterator objects are pointing to the same Edge in the same adjacency list.
         */
        bool operator==(const IncidentIterator &otherIncidentIterator) const {
            return currentIt == otherIncidentIterator.currentIt;
        }


    private:
        friend class Graph;

        // HW1 #3: YOUR CODE HERE
        const Graph *g_;
        size_type startingNodeId{};
        std::unordered_map<size_type, size_type>::const_iterator currentIt;
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
        EdgeIterator() = default;

        /**
         * Constructs an EdgeIterator from an unordered_map whose values are pointers to Edges
         * we will iter through.
         * @param it unordered_map whose values are pointers to Edges
         */
        EdgeIterator(typename std::unordered_map<size_type, Edge *>::const_iterator it) {
            currentIt = it;
        }

        // HW1 #5: YOUR CODE HERE
        /**
         *
         * @return the Edge object currently pointed to by the iterator.
         */
        Edge operator*() const {
            return *((*currentIt).second);
        }

        /**
         *
         * @return an EdgeIterator which starts at the next Edge in the iterator.
         */
        EdgeIterator &operator++() {
            currentIt++;
            return *this;
        }

        /**
         *
         * @param otherEdgeIterator EdgeIterator to compare equality with
         * @return true if both EdgeIterator objects are pointing to the same Edge in the same graph.
         */
        bool operator==(const EdgeIterator &otherEdgeIterator) const {
            return currentIt == otherEdgeIterator.currentIt;
        }

    private:
        friend class Graph;

        // HW1 #5: YOUR CODE HERE
        typename std::unordered_map<size_type, Edge *>::const_iterator currentIt;
    };

    // HW1 #5: YOUR CODE HERE
    /**
     *
     * @return an iterator through the edges of the graph.
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(edgeMap.begin());
    }

    /**
     *
     * @return the end of the iterator through of the edgers of the graph.
     */
    edge_iterator edge_end() const {
        return EdgeIterator(edgeMap.end());
    }

private:
    // HW0: YOUR CODE HERE
    std::vector<NodeInfo> valueVector;
    std::unordered_map<size_type, std::unordered_map<size_type, size_type>> adjacencyLists;

    /** Return a node's position based on its index. */
    const Point &position(size_type index) const {
        // HW0: YOUR CODE HERE
        return valueVector[index]->p_;
    }

    std::unordered_map<size_type, Edge *> edgeMap;
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
