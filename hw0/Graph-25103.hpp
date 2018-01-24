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
class Graph {
private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)
    struct internal_element_node;
    struct internal_element_edge;

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

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        // HW0: YOUR CODE HERE
        begin_edge_index = 0;
        begin_node_index = 0;
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
    class Node {
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
        Node(const Graph* g = nullptr) {
            // HW0: YOUR CODE HERE
            graph = g;
        }

        /** Return this node's position. */
        const Point &position() const {
            // HW0: YOUR CODE HERE
            bool valid = (graph!=nullptr)&&(graph->has_node(*this));
            if(valid){
                return graph->nodes[uid].position;
            }
            else
                return Point();
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            bool valid = (graph!=nullptr)&&(graph->has_node(*this));
            if(valid)
                return uid;
            else
                return size_type(-1);
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            // HW0: YOUR CODE HERE
            return (n.graph ==  graph)&&(n.uid == uid);
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
            return uid<n.uid;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        size_type uid;
        const Graph* graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects


    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        /// HW0: YOUR CODE HERE
        return nodes.size();
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
    Node add_node(const Point &position) {
        // HW0: YOUR CODE HERE
        size_type uid = nodes.size();

        internal_element_node elem;
        elem.position = position;
        nodes.push_back(elem);

        Node node(this);
        node.uid = uid;

        return node;        // Invalid node
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        // HW0: YOUR CODE HERE
        return n.uid<num_nodes() && n.uid >= 0;
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        if(i<num_nodes() && i>=0) {
            Node node(this);
            node.uid = i;

            return node;
        }
        else
            return Node();
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge(const Graph* g=nullptr) {
            graph = g;
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            if(graph!=nullptr)
                return graph->node(n1);
            else
                return Node();
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            if(graph!=nullptr)
                return graph->node(n2);
            else
                return Node();
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            return (e.node1()==node1()&&e.node2()==node2())||(e.node2()==node1()&&e.node1()==node2());
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            Node min1, min2, max1, max2;

            if(node1()<node2()){
                min1 = node1();
                max1 = node2();
            }
            else{
                min1 = node2();
                max1 = node1();
            }

            if(e.node1()<e.node2()){
                min2 = e.node1();
                max2 = e.node2();
            }
            else{
                min2 = e.node2();
                max2 = e.node1();
            }


            return (min1<min2)||(min1==min2&&max1<max2);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        const Graph* graph;
        size_type n1;
        size_type n2;
        size_type uid;
        // HW0: YOUR CODE HERE
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
        return edges.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        if(i<num_edges() && i>=0) {
            Edge edge(this);
            internal_element_edge e = edges[i];
            edge.uid = i;
            edge.n1 = e.n1.uid;
            edge.n2 = e.n2.uid;

            return edge;
        }
        else
            return Edge();
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node &a, const Node &b) const {
        // HW0: YOUR CODE HERE
        for(internal_element_edge e : edges){
            if((e.n1==a && e.n2==b)||(e.n2==a && e.n1==b))
                return true;
        }

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
    Edge add_edge(const Node &a, const Node &b) {
        // HW0: YOUR CODE HERE
//        std::cout<<'Adding';
        for(internal_element_edge e : edges){
            if((e.n1==a && e.n2==b)||(e.n2==a && e.n1==b)){
                Edge edge(this);
                edge.uid = e.uid;
                edge.n1 = e.n1.uid;
                edge.n2 = e.n2.uid;

                return edge;
            }
        }
        internal_element_edge e;
        e.uid = edges.size();
        e.n1 = a;
        e.n2 = b;
        edges.push_back(e);

        Edge edge(this);
        edge.uid = e.uid;
        edge.n1 = e.n1.uid;
        edge.n2 = e.n2.uid;

        return edge;
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        begin_edge_index = edges.size();
        begin_node_index = nodes.size();
        nodes.clear();
        edges.clear();
    }

private:

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    std::vector<internal_element_node> nodes;
    std::vector<internal_element_edge> edges;

    size_type begin_node_index;
    size_type begin_edge_index;

    struct internal_element_node {
        Point position;
        size_type uid;
    };
    struct internal_element_edge {
        Node n1;
        Node n2;
        size_type uid;
    };

};

#endif // CME212_GRAPH_HPP
