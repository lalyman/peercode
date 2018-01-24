#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <set>

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
        struct node_edge_inds;

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
        Graph()
            : nodes(), edges(), edge_ind(), num_nodes_(0), num_edges_(0) {
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
                 * if (...should pick the.node_ind node...)
                 *   x = graph.node(0);
                 * else
                 *   x = some other node using a complicated calculation
                 * do_something(x);
                 * @endcode
                 */
                Node() {
                }

                /** Return this node's position. */
                const Point& position() const {
                    return graph_->nodes[this->index()];
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return ind_;
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    return n.index() == this->index() && n.graph_ == this->graph_;
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
                    return this->index() < n.index();
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                Graph* graph_;  // Graph object that contains the node.
                size_type ind_; // Unique index identifier

                // Private constructor
                Node(const Graph* graph, size_type ind):
                    graph_(const_cast<Graph*>(graph)), ind_(ind) {
                    }
        };

        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         */
        size_type size() const {
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
        Node add_node(const Point& position) {
            nodes[num_nodes_] = position;
            num_nodes_ += 1;
            edges.push_back(std::vector<node_edge_inds>());

            return Node(this, num_nodes_ - 1);
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            return nodes.find(n.index()) != nodes.end();
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
            assert(0 <= i && i < num_nodes());
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
                    return graph_->node(v1_ind_);
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    return graph_->node(v2_ind_);
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    return (v1_ind_ == e.v1_ind_ && v2_ind_ == e.v2_ind_) ||
                           (v1_ind_ == e.v2_ind_ && v2_ind_ == e.v1_ind_);
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    return ind_ < e.ind_;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                Graph* graph_;
                size_type ind_;
                size_type v1_ind_;
                size_type v2_ind_;

                // Private constructor
                Edge(const Graph* graph, size_type ind, size_type v1_ind, size_type v2_ind):
                    graph_(const_cast<Graph*>(graph)), ind_(ind),
                    v1_ind_(v1_ind), v2_ind_(v2_ind) {
                    }
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: O(1)
         */
        size_type num_edges() const {
            return num_edges_;
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: O(1)
         */
        Edge edge(size_type i) const {
            assert(0 <= i && i < num_edges());
            return edge_ind[i];
        }
        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: O(num_edges())
         */
        bool has_edge(const Node& a, const Node& b) const {
            assert(0 <= a.index() && a.index() < num_nodes());
            assert(0 <= b.index() && b.index() < num_nodes());

            for (const node_edge_inds & it : edges[a.index()]) {
                if (it.node_ind == b.index()) {
                    return true;
                }
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
         * Complexity: O(num_edges())
         */
        Edge add_edge(const Node& a, const Node& b) {
            assert(!(a == b));

            size_type a_ind = a.index();
            size_type b_ind = b.index();

            // Check if edge already exists
            for (const node_edge_inds & it: edges[a_ind]) {
                if (it.node_ind == b_ind) {
                    return Edge(this, it.edge_ind, a_ind, b_ind);
                }
            }

            edges[a_ind].push_back({.node_ind = b_ind, .edge_ind = num_edges_});
            edges[b_ind].push_back({.node_ind = a_ind, .edge_ind = num_edges_});

            num_edges_ += 1;

            assert(has_edge(a, b));
            assert(has_edge(b, a));

            Edge e = Edge(this, num_edges_ - 1, a_ind, b_ind);
            edge_ind.push_back(e);

            return e;
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
            nodes.clear();
            edges.clear();
            edge_ind.clear();
            num_edges_ = num_nodes_ = 0;
        }

    private:
        struct node_edge_inds {
            size_type node_ind;
            size_type edge_ind;
        };

        std::unordered_map<size_type, Point> nodes; // Node_ind -> position of point

        // edges[v1_ind] = [(v2_ind, edge_ind), ...]
        std::vector<std::vector<node_edge_inds>> edges;
        std::vector<Edge> edge_ind; // Vector of edge_ind -> Edge object

        size_type num_nodes_;
        size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
