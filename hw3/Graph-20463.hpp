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
template <typename V, typename E = double>
class Graph {
    typedef V node_value_type;
    typedef E edge_value_type;
    private:
        struct edge_info;
        struct node_info;

    public:
        using graph_type = Graph;

        class Node;
        using node_type = Node;

        class Edge;
        using edge_type = Edge;

        class NodeIterator;
        using node_iterator = NodeIterator;

        class EdgeIterator;
        using edge_iterator = EdgeIterator;

        class IncidentIterator;
        using incident_iterator = IncidentIterator;

        /** Type of indexes and sizes.
          Return type of Graph::Node::index(), Graph::num_nodes(),
          Graph::num_edges(), and argument type of Graph::node(size_type) */
        using size_type = unsigned;


        /** Construct an empty graph. */
        Graph()
            : nodes(), edges(), num_nodes_(0), num_edges_(0) {
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
        class Node : private totally_ordered<Node>{
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

                /** Return this node's position **/
                Point& position() {
                    size_type index = this->index();
                    assert(index < graph_->nodes.size());
                    assert(index >= 0);
                    return graph_->nodes[index].pos;
                }

                /** Return this node's position. */
                const Point& position() const {
                    size_type index = this->index();
                    assert(index < graph_->nodes.size());
                    assert(index >= 0);
                    return graph_->nodes[index].pos;
                }


                /** Returns node's previous index **/
                size_type prev_ind() const {
                    return graph_->nodes[this->index()].prev_ind;
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return ind_;
                }

                /* @pre Assumes graph has this node
                 * @brief Returns value of this node **/
                node_value_type& value() {
                    return graph_->nodes[this->index()].val;
                }

                /* @pre Assumes graph has this node
                 * @brief Returns value of this node **/
                const node_value_type& value() const {
                    return graph_->nodes[this->index()].val;
                }

                /* @pre Assumes graph has this node
                 * @brief Returns degree of this node **/
                size_type degree() const {
                    return graph_->edges[this->index()].size();
                }

                /* @brief Returns iterator to edges incident to this node */
                incident_iterator edge_begin() const {
                    return (degree() > 0) ?
                        IncidentIterator(graph_, this->index()) : IncidentIterator();
                }

                /* Returns invalid iterator **/
                incident_iterator edge_end() const {
                    return IncidentIterator();
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
                    return graph_ < n.graph_ || index() < n.index();
                }

            private:
                friend class Graph;
                friend class NodeIterator;
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

        /* @brief removes a node and all its incident edges
         * @param n. Node to be removed
         * @post Invalidates node @n and node with index @param n.graph_->num_nodes()
         * @post All other nodes are still valid
         * @post !has_node(n)
         * @post num_nodes() decremeneted by 1, and num_edges() decremented by n.degree()
         * @post Invalidates all edges incident to node @n, and node(n.graph_->num_nodes()-1)
         * @return 1 if node is removed, 0 if node is not present in the graph
         *
         * Runtime O(n.degree())
         */
        size_type remove_node(const Node& n) {
            if (!has_node(n)) {
                return 0;
            }
            size_type ind = n.index();

            // Remove all edges to node ind

            size_type degree = n.degree();
            for (size_type nind = 0; nind < degree; nind++) {
                edge_info ei = edges[ind][nind];
                size_type v2_ind = ei.v2_ind;
                size_type rev_ind = ei.rev_ind;

                if (edges[v2_ind].size() == 1) {
                    edges[v2_ind].pop_back();
                    continue;
                }
                edge_info e_rep = edges[v2_ind].back();
                edges[v2_ind][rev_ind] = e_rep;
                edges[v2_ind].pop_back();

                if (e_rep.v2_ind == ind) {
                    continue;
                }

                // Update reverse edge (v2_ind, rev_ind) index
                edges[e_rep.v2_ind][e_rep.rev_ind].rev_ind = rev_ind;
            }
            num_edges_ -= degree;
            edges[ind].clear();

            // Relabel last node as node ind
            size_type last_ind = num_nodes() - 1;
            nodes[ind] = nodes.back();
            nodes[ind].prev_ind = last_ind;
            nodes.pop_back();

            // Relabel any edges last_ind has as edges from ind;
            degree = edges[last_ind].size();
            for (size_type nind = 0; nind < degree; nind++) {
                // Update reverse edge to point to index
                edge_info ei = edges[last_ind][nind];
                size_type v2_ind = ei.v2_ind;
                size_type rev_ind = ei.rev_ind;

                edges[v2_ind][rev_ind].v2_ind = ind;

                if (ind < ei.v2_ind && ei.v2_ind < last_ind) {
                    // Update values, since they will now be invalid
                    ei.val = edges[ei.v2_ind][ei.rev_ind].val;
                }
                edges[ind].push_back(ei); // Add edge from last_ind to ind
            }
            edges[last_ind].clear();

            num_nodes_--;

            return 1;
        }

        /* @brief removes node *@param n_it
         * @param[n_it] iterator to node to be removed
         * @post Invalidates iterators to any nodes with index larger than *@param n_it
         * @post Iterators to any node with index less than *@param n_it remain valid
         * @return iterator to next node
         *
         * Runtime O(n.degree())
         */

        // Removes node, returns iterator to a next node.
        node_iterator remove_node(node_iterator n_it) {
            remove_node(*n_it);
            if (n_it.ind_ >= n_it.graph_->num_nodes()) {
                n_it = node_end();
            }
            return n_it;
        }

        /** Add a node to the graph, returning the added node.
         * @param[in] position The new node's position
         * @post new num_nodes() == old num_nodes() + 1
         * @post result_node.index() == old num_nodes()
         *
         * Complexity: O(1) amortized operations.
         */
        Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
            nodes.push_back({.pos = position, .val = value, .prev_ind = 0});
            num_nodes_ += 1;
            edges.push_back(std::vector<edge_info>());

            return Node(this, num_nodes_ - 1);
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            return nodes.size() > n.index() && n.graph_ == this;
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

                /* @pre Assumes graph has this edge e = (v1,v2)
                 * @brief Returns value of this edge
                 * Runtime O(degree(v1))*/
                edge_value_type& value() {
                    size_type v1_ind = v1_ind_;
                    size_type v2_ind = v2_ind_;

                    if (v1_ind > v2_ind) {
                        std::swap(v1_ind, v2_ind);
                    }

                    std::vector<edge_info> eis = graph_->edges[v1_ind];
                    for (size_type nind = 0; nind < eis.size(); nind++) {
                        edge_info& ei = eis[nind];
                        if (ei.v2_ind == v2_ind)
                            return graph_->edges[v1_ind][nind].val;
                    }
                    assert(false);
                }

                /* @pre Assumes graph has this edge e = (v1,v2)
                 * @brief Returns value of this edge
                 * Runtime O(degree(v1))*/
                const edge_value_type& value() const {
                    return const_cast<edge_value_type&>(value());
                }

                /** Returns Euclidean distance between node endpoint positions **/
                double length() const {
                    return norm(graph_->node(v1_ind_).position() -
                                graph_->node(v2_ind_).position());
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    return ((v1_ind_ == e.v1_ind_ && v2_ind_ == e.v2_ind_) ||
                           (v1_ind_ == e.v2_ind_ && v2_ind_ == e.v1_ind_)) &&
                           (graph_   == e.graph_);
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    return graph_  <= e.graph_ || v1_ind_ < e.v1_ind_
                                               || v2_ind_ < e.v2_ind_;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                Graph* graph_;
                size_type v1_ind_;
                size_type v2_ind_;

                // Private constructor
                Edge(const Graph* graph, size_type v1_ind, size_type v2_ind):
                    graph_(const_cast<Graph*>(graph)) {
                        v1_ind_ = v1_ind;
                        v2_ind_ = v2_ind;
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
         * Complexity: O(num_edges()))
         */
        Edge edge(size_type i) const {
            assert(0 <= i && i < num_edges());
            size_type ind = -1;

            for(size_type v1_ind = 0; v1_ind < edges.size(); ++v1_ind) {
                std::vector<edge_info> v1_edges = edges[v1_ind];
                for(size_type nind = 0; nind < v1_edges.size(); ++nind) {
                    size_type v2_ind = v1_edges[nind].v2_ind;
                    if (node(v1_ind) < node(v2_ind))
                        ind++;
                    if (ind == i)
                        return Edge(this, v1_ind, v2_ind);
                }
            }
            assert(false);
            return Edge();
        }
        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: O(num_edges())
         */
        bool has_edge(const Node& a, const Node& b) const {
            for (const edge_info & it : edges[a.index()]) {
                if (it.v2_ind == b.index()) {
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Removes edge (@param a,@param b)
         * @param[a] Node endpoint of edge to be removed
         * @param[b] Node endpoint of edge to be removed
         * @post !has_edge(@param a,@param b)
         * @post @param a.degree( and @param b.degree() both decremented by 1
         * @post num_edges() decremented by 1
         * @post invalidates all edges
         * @return 1 if edge(a,b) is removed, 0 if edge was not present in graph
         *
         * Runtime O(num_edges())
         */
        size_type remove_edge(const Node& a, const Node& b) {
            if (!has_edge(a,b))
                return 0;

            std::vector<edge_info>& a_info = edges[a.index()];
            std::vector<edge_info>& b_info = edges[b.index()];
            for (size_type bind = 0; bind < a_info.size(); ++bind) {
                if (a_info[bind].v2_ind == b.index()) {
                    edge_info bback = b_info.back();
                    edge_info aback = a_info.back();

                    size_type revind = a_info[bind].rev_ind;

                    b_info[revind] = bback;
                    edges[bback.v2_ind][bback.rev_ind].rev_ind = revind;

                    a_info[bind] = aback;
                    if (aback.v2_ind != b.index()) {
                        edges[aback.v2_ind][aback.rev_ind].rev_ind = bind;
                    }

                    a_info.pop_back();
                    b_info.pop_back();

                    break;
                }
            }

            num_edges_--;
            return 0;
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
        Edge add_edge(const Node& a, const Node& b,
                      const edge_value_type& edge_value = edge_value_type()) {
            assert(!(a == b));

            Node v1 = a < b ? a : b;
            Node v2 = a < b ? b : a;

            size_type v1_ind = v1.index();
            size_type v2_ind = v2.index();

            // Check if edge already exists
            for (const edge_info & it: edges[v1_ind]) {
                if (it.v2_ind == v2_ind) {
                    return Edge(this, v1_ind, v2_ind);
                }
            }

            edges[v1_ind].push_back({.v2_ind = v2_ind,
                                    .rev_ind = static_cast<size_type>(edges[v2_ind].size()),
                                     .val = edge_value });
            edges[v2_ind].push_back({.v2_ind = v1_ind,
                                    .rev_ind = static_cast<size_type>(edges[v1_ind].size() - 1),
                                    .val = edge_value });

            num_edges_++;

            assert(has_edge(a, b));
            assert(has_edge(b, a));

            Edge e = Edge(this, v1_ind, v2_ind);

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
            num_edges_ = num_nodes_ = 0;
        }

        //

        /** @class Graph::NodeIterator
         * @brief Iterator class for nodes. A forward iterator. */
        class NodeIterator : private totally_ordered<NodeIterator> {
            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Node;                     // Element type
                using pointer           = Node*;                    // Pointers to elements
                using reference         = Node&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                // HW1 #2: YOUR CODE HERE
                /* Construct an invalid NodeIterator.
                 * @post NodeIterator() == node_end() */
                NodeIterator() {
                    invalidate_node_iterator();
                }

                /* Construct an invalid NodeIterator.
                 * @post NodeIterator() == node_end() */
                void invalidate_node_iterator() {
                    graph_ = nullptr;
                    ind_ = -1;
                    //node_ = &Node();
                }

                /** @return node referenced to **/
                Node operator*() const {
                    if (ind_ == size_type(-1)) {
                        std::cout << "returning invalid node " << std::endl;
                        return Node();
                    }
                    return graph_->node(ind_);
                }

                /* @pre Assumes graph_ is valid or node iterator is invalid
                 * @ brief Increments node iterator. Node() is invalidated if
                 * ind_ is invalid or graph_->num_nodes() - 1
                 * @return incremented iterator
                 */
                NodeIterator& operator++() {
                    if (ind_ == size_type(-1)) {
                        return *this;
                    }

                    if (ind_ >= graph_->num_nodes() - 1) {
                        invalidate_node_iterator();
                        return *this;
                    }
                    ind_++;
                    return *this;
                }
                bool operator==(const NodeIterator& it) const {
                    return ind_ == it.ind_ && graph_ == it.graph_;
                }
                size_type get_ind() {
                    return ind_;
                }

            private:
                friend class Graph;

                Graph* graph_;
                size_type ind_; // Denotes our index into the node list

                NodeIterator(const Graph* graph, size_type ind):
                    graph_(const_cast<Graph*>(graph)), ind_(ind) {
                    }
        };

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        /* @returns node_iterator to the first node in the graph
         * @post Node() is valid iff NodeIterator() is valid
         *
         */
        node_iterator node_begin() const {
           return (num_nodes_ > 0) ?
               NodeIterator(this, 0) : NodeIterator();
        }
        /* @returns invalid node iterator
         * @post NodeIterator() == node_end()
         */
        node_iterator node_end() const {
            return NodeIterator();
        }

        //
        // Incident Iterator
        //

        /** @class Graph::IncidentIterator
         * @brief Iterator class for edges incident to a node. A forward iterator. */
        class IncidentIterator : private totally_ordered<IncidentIterator>{
            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator() {
                    invalidate_incident_iterator();
                }

                // HW1 #3: YOUR CODE HERE
                // @return edge object
                Edge operator*() const {
                    return edge_;
                }

                // @brief invalidates incident iterator
                // @post edge_ is invalid
                // @post *this == edge_end()
                void invalidate_incident_iterator() {
                    v1_ind_ = n_ind_ = -1;
                    graph_ = nullptr;
                    edge_ = Edge();
                }

                /* @brief increments incident iterator.
                 * @return incremented incident iteratoer
                 * @post invalid iterator remains invalid */
                IncidentIterator& operator++() {
                    if (n_ind_ >= graph_->node(v1_ind_).degree() - 1) {
                        invalidate_incident_iterator();
                        return *this;
                    }
                    n_ind_++;
                    edge_info e_info = graph_->edges[v1_ind_][n_ind_];
                    edge_ = Edge(graph_, v1_ind_, e_info.v2_ind);

                    return *this;
                }

                /* @brief returns edge info associated with current iterator **/
                edge_info& get_edge_info() {
                    return e_info;
                }

                /* @brief checks equality incident iterators. Incident
                 * iterators are equal if their graphs are the same,
                 * their neighbor indices are the same, and their
                 * node indices are the same.
                 * @return true if iterators are equal, false otherwise */
                bool operator==(const IncidentIterator& it) const {
                    return v1_ind_ == it.v1_ind_ && graph_ == it.graph_
                                                 && n_ind_ == it.n_ind_;
                }

            public:
                edge_info e_info;
            private:
                friend class Graph;

                // HW1 #3: YOUR CODE HERE
                Graph* graph_;
                size_type v1_ind_;
                size_type n_ind_; // Neighbor index
                Edge edge_;

                /* @brief Private constructor
                 * @pre Assumes graph contains node v1_ind_ */
                IncidentIterator(const Graph* graph, size_type v1_ind_):
                    graph_(const_cast<Graph*>(graph)), v1_ind_(v1_ind_) {
                        n_ind_ = 0;
                        e_info = graph_->edges[v1_ind_][n_ind_];
                        size_type v2_ind_ = e_info.v2_ind;
                        edge_ = Edge(graph_, v1_ind_, v2_ind_);
                    }
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

                const size_type invalid_ind = -1;

                /** Construct an invalid EdgeIterator. */
                EdgeIterator() {
                    invalidate_edge_iterator();
                }

                void invalidate_edge_iterator() {
                    v1_ind_ = v2_ind_ = nind_ = invalid_ind;
                    graph_ = nullptr;

                }

                // HW1 #5: YOUR CODE HERE
                /* @brief returns edge object referred to.
                 * @return dereferenced edge
                 * @post returned edge is invalid iff v1_ind_, v2_ind_, or graph_
                 * are invalid
                 **/
                Edge operator*() const {
                    if (v1_ind_ == invalid_ind || v2_ind_ == invalid_ind || !graph_)
                        return Edge();
                    return Edge(graph_, v1_ind_, v2_ind_);
                }

                /* @brief Increments edge iterator. Skips duplicate edges.
                 * Ammortized O(1) time
                 * @return incremented edge iterator
                 **/
                EdgeIterator& operator++() {
                    if (v1_ind_ == graph_->num_nodes() - 1
                       && nind_ == Node(graph_, v1_ind_).degree() - 1) {
                        invalidate_edge_iterator();
                        return *this;
                    }
                    Node v1 = Node(graph_, v1_ind_);
                    std::vector<edge_info> v1_edges = graph_->edges[v1_ind_];
                    while (true) {
                        nind_++;

                        if (nind_ > v1.degree() - 1 || v1.degree() == 0) {
                            v1_ind_++;
                            if (v1_ind_ > graph_->num_nodes() - 1) {
                                invalidate_edge_iterator();
                                return *this;
                            }
                            v1 = Node(graph_, v1_ind_);
                            v1_edges = graph_->edges[v1_ind_];
                            nind_ = 0;
                        }
                        Node v2 = Node(graph_, v1_edges[nind_].v2_ind);
                        v2_ind_ = v1_edges[nind_].v2_ind;

                        if (v1 < v2) { // Don't double-count edges
                            return *this;
                        }
                    }
                    assert(false);
                    return *this;
                }
                bool operator==(const EdgeIterator& ei) const {
                    return v1_ind_ == ei.v1_ind_ && nind_ == ei.nind_
                                                 &&  graph_ == ei.graph_;
                }

            private:
                friend class Graph;
                // HW1 #5: YOUR CODE HERE
                Graph* graph_;
                size_type v1_ind_; // Index of first node in edge
                size_type v2_ind_; // Index of second node in edge
                size_type nind_;   // Neighbor index of second node in neighbors of first node

                /* @brief Constructs edge iterator object
                 * @pre Assumes node with index v1_ind_ has degree at least nind_
                 * @param graph valid graph object
                 * @param v1_ind_ index of first node in edge
                 * @param nind_ neighbor index of second node in neighbors of first node */
                EdgeIterator(const Graph* graph, size_type v1_ind_, size_type nind_):
                    graph_(const_cast<Graph*>(graph)), v1_ind_(v1_ind_), nind_(nind_) {
                        std::vector<edge_info> v1_edges = graph_->edges[v1_ind_];
                        v2_ind_ = v1_edges[nind_].v2_ind;
                    }
        };

        // HW1 #5: YOUR CODE HERE
        /* @brief Creates edge iterator starting at beginning of 0'th vertex
         * @return edge iterator to first edge
         * @post Edge(graph_, 0) == *(@return)
         */
        edge_iterator edge_begin() const {
            return (num_edges() > 0) ?
                EdgeIterator(this, 0, 0) : EdgeIterator();
        }

        /* @brief Creates invalid edge iterator
         * @return invalid edge iterator
         * @post @return == EdgeIterator()
         */
        edge_iterator edge_end() const {
            return EdgeIterator();
        }

        // @brief Prints graph.
        void print() {
            std::cout << "G(|V| = " << num_nodes()
                      << ", |E| = " << num_edges() << ")" << std::endl;
            for(size_type i = 0; i < size(); i++) {
                Node n = node(i);
                Point pos = n.position();
                std::cout << "Node " << i << " at (" << pos.x <<
                    ", " << pos.y << ", " << pos.z << ")" << std::endl;
            }
            for (size_type i = 0; i < size(); i++) {
                Node n = node(i);
                std::cout << i << " --> ";
                for (size_type j = 0 ; j < n.degree(); j++) {
                    auto ei = edges[i][j];
                    std::cout << "(" << ei.v2_ind << "), ";
                }
                std::cout << std::endl;
            }
        }

        /** @brief Checks if edge adjacency list indices are consistent.
         *  @return true if for each edge, the reverse edge is in the list,
         *  and with the correct rev_ind index. Returns false otherwise
         *
         *  Useful for debugging
         */
        bool consistent() {
            for (size_type v1 = 0; v1 != num_nodes(); ++v1) {
                std::vector<edge_info> v1_info = edges[v1];
                for (size_type nind = 0; nind != v1_info.size(); ++nind) {
                    size_type v2_ind = v1_info[nind].v2_ind;
                    size_type rev_ind = v1_info[nind].rev_ind;

                    edge_info rev_edge = edges[v2_ind][rev_ind];
                    if (rev_edge.v2_ind != v1 ||
                        rev_edge.rev_ind != nind) {
                        print();
                        std::cout << "Not consistent at " << v1 << " neighbor "
                                  << nind << std::endl;
                        return false;
                    }
                }
            }
            return true;
        }

    private:
        struct edge_info {
            size_type v2_ind;
            size_type rev_ind; // Index of reverse edge in v2's adjacency list
            edge_value_type val;
        };

        struct node_info {
            Point pos;
            node_value_type val;
            size_type prev_ind; // Previous index, in case has been relabeled
        };

        std::vector<node_info> nodes; // i'th entry is i'th node

        // Adjacency list of edges. Edges stored forward and backward
        std::vector<std::vector<edge_info>> edges;

        size_type num_nodes_;
        size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
