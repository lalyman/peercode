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

template<typename V>

class Graph {

public:
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    /** Type of this graph. */
    using size_type = unsigned;

private:
    struct edge_struct;
    struct node_struct;
    vector<node_struct> nodes;
    vector<vector<size_type>> adjacency_list;
    vector<edge_struct> edge_list;

public:
    /* Type of graph*/
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


    /** Construct an empty graph. */
    Graph() : nodes(), adjacency_list(), edge_list() {
    }

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
        Node() {
            // HW0: YOUR CODE HEREbu
        }

        /** Return this node's position. */
        const Point &position() const {
            return fetch_node().point;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        // The .idx may be redundant
        size_type index() const {
            return fetch_node().idx;
        }



        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */

        /**
         *value() returns a referece to Node's value so that .value() can be called on node
         * @pre must call on valid node within the graph
         * @return reference to @ val for node of NodeType V
         */
        node_value_type &value() {
            return fetch_node().val;
        }

        /**
         *value() returns a referece to Node's value so that .value() can be called on node when
         * @pre must call on valid node within the graph
         * @return reference to const @val for node of NodeType V
         */
        const node_value_type &value() const {
            return fetch_node().val;
        }

        /**
         *degree() determines the number of edges connected to a Node
         * @pre Must be called on valid node in a graph
         * @return @a deg of type size_type
         */
        size_type degree() const {
            size_type deg = static_cast<size_type>(graph_->adjacency_list[uid_].size());
            return deg;
        }

        /**
         *edge_begin() constructs a @a Iiter Incident_iterator pointing to first neighbor
         * to which the node it is called on is connected
         * @pre must be called on a valid node
         * @return  @a IIter object of type incident_iterator
         */
        incident_iterator edge_begin() const {
            incident_iterator IIter(uid_, 0, const_cast<Graph *>(graph_));
            return IIter;
        }


        /**
         *edge_endn() constructs a @a Iiter Incident_iterator pointing to one past the last
         * to which the node it is called on is connected
         * @pre must be called on a valid node
         * @return  @a IIter object of type incident_iterator
         */
        incident_iterator edge_end() const {
            incident_iterator IIter(uid_, degree(), const_cast<Graph *>(graph_));
            return IIter;
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

        node_struct &fetch_node() const {
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

    Node add_node(const Point &position, const node_value_type &a = node_value_type()) {
/* DEBUG_MSG("Adding node"); DEBUG_MSG("Num nodes is " << nodes.size());
        node_struct new_node;
        new_node.point = position;
        new_node.val = a;
        new_node.idx =  static_cast<size_type>(nodes.size(); */
        node_struct new_node{
                .point = position,
                .val = a,
                .idx =  static_cast<size_type>(nodes.size())
        };
        nodes.emplace_back(new_node);
        /*DEBUG_MSG("Adding with index " << nodes.size() - 1);*/
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
    class Edge : private totally_ordered<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge() {
        }

        /** Return a node of this Edge */
        Node node1() const {
            return graph_->node(n1_id_);

        }

        /** Return the other node of this Edge */
        Node node2() const {
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
            } else if (graph_ < e.graph_) {
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
    // I could use my edge iterator here using *std::next(edge_begin(), i)
    // But this is faster and memory doesnt seem to be an issue
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
        n_edges += 1;
        adjacency_list[a.uid_].push_back(b.uid_);
        adjacency_list[b.uid_].push_back(a.uid_);
        edge_struct new_edge{
                .n1_id = a.uid_,
                .n2_id = b.uid_
        };

//        edge_struct new_edge;
//        new_edge.n1_id = a.uid_;
//        new_edge.n2_id = b.uid_;
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

        /**
         * operator*() returns the Node at the current position of the NodeIterator
         * @pre iterator node at last node of graph (graph.node_end())
         * @return Node at current step of iterator
         */
        Node operator*() const {
            return ngraph_->node(nidx_);
        }

        /**
         *operator++() increments the postition of the iterator by 1
         *@returns reference to iterator after incrementing
         *@pre: iterator is at the last node of the graph
         */
        NodeIterator &operator++() {
            nidx_ = nidx_ + 1;
            return *this;
        }


        /**
         * operator==() compares iterator on which it is called with @a other NodeIterator
         * @param[in] @a other node iterator for comparison
         * @return a boolean if the two iterators are acting on the same graph and at the same position
         */
        bool operator==(const NodeIterator &other) const {
            return other.ngraph_ == ngraph_ && other.nidx_ == nidx_;

        }

    private:
        friend class Graph;

        size_type nidx_;
        graph_type *ngraph_;

        NodeIterator(size_type idx, const graph_type *graph) : nidx_(idx), ngraph_(const_cast<Graph *>(graph)) {
        };


    };


    /**
     * node_begin() returns a pointer to the beggining of the node container of the graph
     * @return a node_iterator that points to the node with id zero
     **/
    node_iterator node_begin() const {
        return NodeIterator(0, this);
    }


    /**
     * node_end() returns a pointer the the end of nodes + 1
     * @return returns a node iterator pointing to end of nodes + 1
     */
    node_iterator node_end() const {
        return NodeIterator(size(), this);
    }


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


        /**
         *operator*() returns the Edge at the current position of the IncidentIterator
         * @pre the iterator is not at node.edge_end()
         * @return An Edge object corresponding to the current position
         */
        Edge operator*() const {
            return Edge(Igraph_, outer_id_, Igraph_->adjacency_list[outer_id_][inner_id_]);
        }


        /**
         *operator++() increments the IncidentIterator by one
         * @pre the iterator is not at node.edge_end()
         * @return reference to the iterator after incrementing
         */
        IncidentIterator &operator++() {
            inner_id_ = inner_id_ + 1;
            return *this;

        }


        /**
         *operator==() compares the given IncidentIterator to @a iit IncidentIterator passed in
         * @return boolean true if the two iterators act on the same graph, the same node, and the same incident edge
          * @param[in] other node iterator @a iit for comparison
         */
        bool operator==(const IncidentIterator &iit) const {
            return iit.Igraph_ == Igraph_ && iit.inner_id_ == inner_id_ && iit.outer_id_ == outer_id_;
        }

    private:
        friend class Graph;
        size_type outer_id_;
        size_type inner_id_;
        graph_type *Igraph_;

        IncidentIterator(size_type outer_id, size_type inner_id, const graph_type *graph) : outer_id_(outer_id),
                                                                                            inner_id_(inner_id),
                                                                                            Igraph_(const_cast<Graph *> (graph)) {

        }
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


        /**
         *operator*() returns the edge at the current position of the EdgeIterator
         * @return Edge object
         * @pre the EdgeIterator != graph.edge_end()
         */
        Edge operator*() const {
            return *curr_inc;
        }

        /**
         *operator++() increments the EdgeIterator by one
         * @pre EdgeIterator != graph.edge_end()
         * @return a reference to the EdgeIterator after incrementing by one
         */
        EdgeIterator &operator++() {
            ++curr_inc;
            fix_edges();
            return *this;

        }

        /**
         *operator==() compares the given EdgeIterator to the supplied EdgeIterator @a eother
         * @param[in] @a eother is the EdgeIterator object to compare to
         * @return boolean type
         */
        bool operator==(const EdgeIterator &eother) const {
            return eother.Egraph_ == Egraph_ && eother.curr_node == curr_node && eother.curr_inc == curr_inc;
        }

    private:
        friend class Graph;

        graph_type *Egraph_;
        NodeIterator curr_node = (*Egraph_).node_begin();
        IncidentIterator curr_inc = (*curr_node).edge_begin();

        EdgeIterator(const graph_type *g, NodeIterator N, IncidentIterator I) :
                Egraph_(const_cast<Graph *>(g)), curr_node(N), curr_inc(I) {
            fix_edges();
        }

        void fix_edges() {
            int break_cond = false;
            // while we arent at the last node
            while (curr_node != Egraph_->node_end()) {
                // not at the last incident edge
                while (curr_inc != (*curr_node).edge_end()) {
                    // specify ordering to prevent duplicates
                    // If the current node is less than the curr_incident node
                    // we can break, onyl seeing repeats
                    if ((*curr_node) < (*curr_inc).node2()) {
                        break_cond = true;
                        break;
                    }
                        // Otherwise increment
                    else {
                        ++curr_inc;
                    }
                }
                if (break_cond) {
                    break;
                }
                ++curr_node;

                if (curr_node == Egraph_->node_end()) {
                    curr_inc = (*(Egraph_->node_begin())).edge_begin();
                    break;
                } else {
                    curr_inc = (*curr_node).edge_begin();
                }

            }
        }

    };

    /**
     *edge_begin() returns an  pointer to an edge at the first edge of the graph
     *
     * @return of type edge_iterator
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, node_begin(), (*node_begin()).edge_begin());
    }

    /**
     *edge_end() returns a pointer to one past the last edge in the graph
     * @return of type edge_iterator
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, node_end(), (*node_begin()).edge_begin());

    }


private:
    size_type n_edges = 0;
    struct node_struct {
        Point point;
        // THis may be redundant whatever
        size_type idx;
        node_value_type val;

    };
    struct edge_struct {
        size_type n1_id;
        size_type n2_id;
    };


};

#endif // CME212_GRAPH_HPP
