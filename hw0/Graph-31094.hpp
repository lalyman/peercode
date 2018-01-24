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
#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

using namespace std;
/** @class Graph
 * @brief A template for 3D undirected graphs.
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
public:
    //
    // PUBLIC TYPE DEFINITIONS
    //
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;


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

private:

    struct edge_struct;
    struct node_struct;
    vector<node_struct> nodes;
    vector<vector<size_type>> adjacency_list;
    vector<edge_struct> edge_list;

public:


    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph(): nodes(),adjacency_list(), edge_list() {
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
        Node() {
            // HW0: YOUR CODE HEREbu
        }

        /** Return this node's position. */
        const Point &position() const {
            return fetch().point;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        // The .idx may be redundant
        size_type index() const {
            return fetch().idx;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            return (n.graph_ == graph_) && (uid_ == n.uid_);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node &n) const {
            if (graph_ == n.graph_) {
                return uid_ < n.uid_;
            } else {
                return graph_ < n.graph_;
            }
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        graph_type *graph_;
        size_type uid_;

        /** Private Constructor **/
        Node(const graph_type *graph, size_type uid)
                : graph_(const_cast<graph_type *>(graph)), uid_(uid) {
        }

        node_struct &fetch() const {
            assert(uid_ < graph_->size());
            return graph_->nodes[uid_];
        }

        friend class Graph;
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return static_cast<size_type>(nodes.size());
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
        DEBUG_MSG("Adding node");
        DEBUG_MSG("Num nodes is " <<  nodes.size());
        node_struct new_node;
        new_node.point = position;
        new_node.idx = static_cast<size_type>(nodes.size());
        nodes.emplace_back(new_node);
        DEBUG_MSG("Adding with index " <<  nodes.size() - 1);
        adjacency_list.emplace_back(vector<size_type>());
        return Node(this, static_cast<size_type>(nodes.size() - 1));
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        // HW0: YOUR CODE HERE
        if (n.graph_ == const_cast<graph_type *>(this)) {
            return static_cast<size_type>(nodes.size()) > n.uid_;
        } else {
            return false;
        }
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        assert(i < size());
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {
        }
        /** Return a node of this Edge */
        Node node1() const {
            DEBUG_MSG("Getting e1 " << n1_id_);
            return graph_->node(n1_id_);

        }
        /** Return the other node of this Edge */
        Node node2() const {
            DEBUG_MSG("Getting e2 " << n2_id_);
            return graph_->node(n2_id_);
        }
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            if (graph_ == e.graph_) {
                return (min(n1_id_, n2_id_) == min(e.n1_id_, e.n2_id_)) &&
                       max(n1_id_, n2_id_) == max(e.n1_id_, e.n2_id_);
            }
            return false;
        }
        /** Test whether this edge is less than @a e in a global order.
         * This ordering function is useful for STL containers such as
         * map<>. It need not have any interpretive meaning.
         */
        /* Guess just have minimum */
        bool operator<(const Edge &e) const {
            if (graph_ == e.graph_) {
                if (min(n1_id_, n2_id_) < min(e.n1_id_, e.n2_id_)) {
                    return true;
                } else if (min(n1_id_, n2_id_) == min(e.n1_id_, e.n2_id_)) {
                    return max(n1_id_, n2_id_) < max(e.n1_id_, e.n2_id_);
                }
            }
            else if (graph_ < e.graph_){
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        graph_type *graph_;
        size_type n1_id_;
        size_type n2_id_;

        Edge(const graph_type *graph, size_type n1_id, size_type n2_id)
                : graph_(const_cast<graph_type *>(graph)),
                  n1_id_(graph->nodes[n1_id].idx),
                  n2_id_(graph->nodes[n2_id].idx) {
        }
    };

    /** Return the total number of edges in the graph.
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return n_edges;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        assert(i < n_edges);
        return Edge(this, edge_list[i].n1_id, edge_list[i].n2_id);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node &a, const Node &b) const {
        for (size_type i = 0; i < adjacency_list[a.uid_].size(); ++i) {
            if (b.uid_ == adjacency_list[a.uid_][i]) {
                return true;
            }
        }
        /**     for (size_type i = 0; i < num_edges(); ++i){
                 if (edge_list[i].n1_id == a.uid_  && edge_list[i].n2_id == b.uid_){
                     return true;
                 }
                 if (edge_list[i].n1_id == b.uid_ && edge_list[i].n2_id == a.uid_ ){
                     return  true;
                 }

             } **/
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

        if (has_edge(a, b)) {
            return Edge(this, a.uid_, b.uid_);
        }
        DEBUG_MSG("adding edge " << a.uid_ << " to " << b.uid_ );
        n_edges += 1;
        adjacency_list[a.uid_].push_back(b.uid_);
        adjacency_list[b.uid_].push_back(a.uid_);
        edge_struct new_edge;
        new_edge.n1_id = a.uid_;
        new_edge.n2_id = b.uid_;
        edge_list.emplace_back(new_edge);
        return Edge(this, a.uid_, b.uid_);        // Invalid Edge
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes.clear();
        adjacency_list.clear();
        edge_list.clear();
    }

private:
    size_type n_edges = 0;
    struct node_struct {
        Point point;
        // THis may be redundant whatever
        size_type idx;
    };
    struct edge_struct {
        size_type n1_id;
        size_type n2_id;
    };
};

#endif // CME212_GRAPH_HPP
