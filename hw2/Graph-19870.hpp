#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
public:

        //
        // PUBLIC TYPE DEFINITIONS
        //

        /** Type of this graph. */
        using graph_type = Graph;

        /** The type of the value in the nodes of this Graph */
        using node_value_type = V;

        /** The type of the value in the edges of this Graph */
        using edge_value_type = E;

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

        /** Type of indexes and sizes.
         Return type of Graph::Node::index(), Graph::num_nodes(),
         Graph::num_edges(), and argument type of Graph::node(size_type)
         */
        using size_type = unsigned;

        //
        // CONSTRUCTORS AND DESTRUCTOR
        //

        /** Construct an empty graph. */
        Graph()
        : graph_nodes_(), graph_edges_() {
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
                        // no code needed
                }

                /** Return this node's position. */
                Point& position() {
                        return(this->node_graph_->graph_nodes_[this->node_idx_].pos);
                }

                /** Return this node's position. */
                const Point& position() const {
                        return(this->node_graph_->graph_nodes_[this->node_idx_].pos);
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                        return(this->node_idx_);
                }

                /** Return the value of the node */
                node_value_type& value() {
                        return(this->node_graph_->graph_nodes_[this->node_idx_].node_val);
                }

                /** Return the value of the node */
                const node_value_type& value() const {
                        return(this->node_graph_->graph_nodes_[this->node_idx_].node_val);
                }

                /** Return the degree of the node.
                 *   The degree of a node is the number of nodes that this node is adjacent to
                 */
                size_type degree() const {
                        return(this->node_graph_->graph_nodes_[this->node_idx_].adj_edges.size());
                }

                /** Return an incident iterator of the begining */
                incident_iterator edge_begin() const {
                        return(IncidentIterator(this->node_graph_, this->node_idx_));
                }

                /** Return an incident iterator of the end */
                incident_iterator edge_end() const {
                        return(IncidentIterator(this->node_graph_, this->node_idx_, this->degree()));
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                        return(n.node_graph_ == this->node_graph_ && n.node_idx_ == this->node_idx_);
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
                        if(n.node_graph_ == this->node_graph_) {
                                return(this->node_idx_ < n.node_idx_);      // no possibility of n.index == this->node_idx_
                        } else {
                                // from different graphs
                                if(n.node_idx_ != this->node_idx_) {
                                        return(this->node_idx_ < n.node_idx_);  // if the indexes are different
                                } else {
                                        // if the index is the same and from different graphs, then ** this ** is less
                                        return(true);
                                }
                        }
                }

        private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;

                // pointer to the graph of the node
                Graph* node_graph_;

                // the id of the node
                size_type node_idx_;

                /** Construct a node object.
                 * @param[in] graph     A pointer to the graph of the node
                 * @param[in] node_idx  The id of the node to construct
                 */
                Node(const Graph* graph, size_type node_idx)
                : node_graph_(const_cast<Graph*>(graph)), node_idx_(node_idx) {
                }
        };

        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         */
        size_type size() const {
                return(this->graph_nodes_.size());
        }

        /** Synonym for size(). */
        size_type num_nodes() const {
                return(size());
        }

        /** Add a node to the graph, returning the added node.
         * @param[in] position The new node's position
         * @post new num_nodes() == old num_nodes() + 1
         * @post result_node.index() == old num_nodes()
         *
         * Complexity: O(1) amortized operations.
         */
        Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
                proxy_node new_node = proxy_node(position, value);
                size_type id = this->graph_nodes_.size();
                this->graph_nodes_.push_back(new_node);
                return(Node(this, id));
        }

        /** Remove node from the graph.
         * @param[in] n 	node to remove from the graph
         * @return True if the node was removed, false otherwise.
         *
         * @post @a n was removed from the graph, including all of its incident edges.
         *
         * ALL GRAPH INVARIENTS ARE PRESERVED DURING OPERATION.
         *
         * Complexity: O(num_nodes())
         */
        size_type remove_node(const Node &n) {
                if(!has_node(n)) return(0);
                std::vector<Edge> edges_to_remove;
                for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
                        edges_to_remove.push_back(*it);
                }
                for(size_type i = 0; i < edges_to_remove.size(); i++) {
                        remove_edge(edges_to_remove[i]);
                }
                this->graph_nodes_.erase(this->graph_nodes_.begin() + n.index());
                size_type counter = 0;
                for(auto it = node_begin(); it != node_end(); ++it) {
                        Node nd = *it;
                        nd.node_idx_ = counter;
                        counter++;
                }
                return(1);
        }

        /** Remove node from the graph.
         * @param[in] n_it	 	iterator to node to remove from the graph
         * @return Iterator to the node that was removed
         *
         * @post node pointed to by @a n_it was removed from the graph, including all of its incident edges.
         *
         * ALL GRAPH INVARIENTS ARE PRESERVED DURING OPERATION.
         *
         * Complexity: O(num_nodes())
         */
        node_iterator remove_node(node_iterator n_it) {
                remove_node(*n_it);
                return(NodeIterator(this, (*n_it).index()));
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
                return(n.node_graph_ == this);
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
                return(Node(this, i));
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
        class Edge : private totally_ordered<Edge>{
        public:
                /** Construct an invalid Edge. */
                Edge() {
                        // no code here
                }

                /** Return a node of this Edge */
                Node node1() const {
                        return(Node(this->edge_graph_, this->edge_node_1_));
                }

                /** Return the other node of this Edge */
                Node node2() const {
                        return(Node(this->edge_graph_, this->edge_node_2_));
                }

                /** Return the value of the edge */
                edge_value_type& value() {
                        return(this->edge_graph_->graph_edges_[this->edge_idx_].edge_val);
                }

                /** Return the value of the edge */
                const edge_value_type& value() const {
                        return(this->edge_graph_->graph_edges_[this->edge_idx_].edge_val);
                }

                /** Return the length of the edge */
                double length() const {
                        return(norm(this->node1().position() - this->node2().position()));
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                        return(e.edge_graph_ == this->edge_graph_ && this->edge_idx_ == e.edge_idx_);
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                        if(e.edge_graph_ == this->edge_graph_) {
                                return(this->edge_idx_ < e.edge_idx_);      // no possibility of e.edge_idx_ == this->edge_idx_
                        } else {
                                // from different graphs
                                if(e.edge_idx_ != this->edge_idx_) {
                                        return(this->edge_idx_ < e.edge_idx_);  // if the indexes are different
                                } else {
                                        // if the index is the same and from different graphs, then ** this ** is less
                                        return(true);
                                }
                        }
                }

        private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                // pointer to the graph of the edge
                Graph* edge_graph_;

                // the id of the edge
                size_type edge_idx_;

                // id of the root node
                size_type edge_node_1_;

                // id of the leaf node
                size_type edge_node_2_;

                /** Construct a edge object.
                 * @param[in] graph     pointer to the graph of the edge
                 * @param[in] edge_idx  id of the edge to construct
                 * @param[in] node_1    id of the root node
                 * @param[in] node_2    id of the leaf node
                 */
                Edge(const Graph* graph, size_type edge_idx, size_type node_1, size_type node_2)
                : edge_graph_(const_cast<Graph*>(graph)), edge_idx_(edge_idx), edge_node_1_(node_1), edge_node_2_(node_2) {
                }
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const {
                return(this->graph_edges_.size());
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const {
                size_type node_a = this->graph_edges_[i].edge_node_a;
                size_type node_b = this->graph_edges_[i].edge_node_b;
                return(Edge(this, i, node_a, node_b));
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(a.degree()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const {
                for(auto it = a.edge_begin(); it != a.edge_end(); ++it) {
                        Edge e = *it;
                        if(e.node1() == a && e.node2() == b) {
                                return(true);
                        }
                }
                return(false);
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
         * Complexity: No more than O(a.degree()), hopefully less
         */
        Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
                if(has_edge(a, b)) {
                        for(auto it = a.edge_begin(); it != a.edge_end(); ++it) {
                                Edge e = *it;
                                if((e.node1() == a && e.node2() == b) ||
                                   (e.node1() == b && e.node2() == a)) {
                                        return(e);
                                }
                        }
                }
                proxy_edge new_proxy_edge = proxy_edge(a.index(), b.index(), value);
                Edge new_edge = Edge(this, this->num_edges(), a.index(), b.index());
                this->graph_nodes_[a.index()].adj_edges.push_back(this->num_edges());       // add edge that is adj to a
                this->graph_nodes_[b.index()].adj_edges.push_back(this->num_edges());       // add edge that is adj to b
                this->graph_edges_.push_back(new_proxy_edge);                               // add edge to vector of edges
                return(new_edge);
        }

        /** Remove edge from the graph.
         * @param[in] a    node1 of the edge to remove
         * @param[in] b    node2 of the edge to remove
         * @return True if the node was removed, false otherwise.
         *
         * @post edge e defined by @a a and @a b was removed from the graph
         *
         * ALL GRAPH INVARIENTS ARE PRESERVED DURING OPERATION.
         *
         * Complexity: O(num_edges())
         */
        size_type remove_edge(const Node &a, const Node& b) {
                for(auto it = a.edge_begin(); it != a.edge_end(); ++it) {
                        Edge e = *it;
                        if((e.node1() == a && e.node2() == b) ||
                           (e.node1() == b && e.node2() == a)) {
                                return(remove_edge(e));
                        }
                }
                return(0);
        }

        /** Remove edge from the graph.
         * @param[in] e    edge to remove
         * @return True if the node was removed, false otherwise.
         *
         * @post edge @a e was removed from the graph
         *
         * ALL GRAPH INVARIENTS ARE PRESERVED DURING OPERATION.
         *
         * Complexity: O(num_edges())
         */
        size_type remove_edge(const Edge& e) {
                if(!has_edge(e.node1(), e.node2())) return(0);
                this->graph_edges_.erase(this->graph_edges_.begin() + e.edge_idx_);
                size_type counter = 0;
                for(auto it = edge_begin(); it != edge_end(); ++it) {
                        Edge ed = *it;
                        ed.edge_idx_ = counter;
                        counter++;
                }
                // removes e from node1 adjacency list
                for(size_type i = 0; i < this->graph_nodes_[e.node1().index()].adj_edges.size(); i++) {
                        if(e.edge_idx_ == this->graph_nodes_[e.node1().index()].adj_edges[i]) {
                                this->graph_nodes_[e.node1().index()].adj_edges.erase(this->graph_nodes_[e.node1().index()].adj_edges.begin() + i);
                                break;
                        }
                }
                // removes e from node2 adjacency list
                for(size_type i = 0; i < this->graph_nodes_[e.node2().index()].adj_edges.size(); i++) {
                        if(e.edge_idx_ == this->graph_nodes_[e.node2().index()].adj_edges[i]) {
                                this->graph_nodes_[e.node2().index()].adj_edges.erase(this->graph_nodes_[e.node2().index()].adj_edges.begin() + i);
                                break;
                        }
                }
                return(1);
        }

        /** Remove edge from the graph.
         * @param[in] e_it  	iterator to edge to remove
         * @return An iterator to the edge that was removed
         *
         * @post edge pointed to by @a e_it was removed from the graph
         *
         * ALL GRAPH INVARIENTS ARE PRESERVED DURING OPERATION.
         *
         * Complexity: O(num_edges())
         */
        edge_iterator remove_edge(edge_iterator e_it) {
                remove_edge(*e_it);
                return(NodeIterator(this, (*e_it).edge_idx_));
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
                this->graph_edges_.clear();
                this->graph_nodes_.clear();
        }

        //
        // Node Iterator
        //

        /** @class Graph::NodeIterator
         * @brief Iterator class for nodes. A forward iterator.
         */
        class NodeIterator : private totally_ordered<NodeIterator> {
        public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Node;                     // Element type
                using pointer           = Node*;                    // Pointers to elements
                using reference         = Node&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid NodeIterator. */
                NodeIterator() {
                        // no code needed
                }

                /** Derefrence the node iterator to get the node at the iterators current location.
                 * @pre @a it >= 0 && it < Number of nodes in the graph
                 * @return The node that is at the current loaction in the graph
                 *
                 * Has a possibility of returning a node that is not in the graph if
                 *    @a it >= Number of nodes in the graph.
                 *
                 * Complexity: O(1)
                 */
                Node operator*() const {
                        return(Node(this->node_iterator_graph_, this->graph_pos_));
                }

                /** Increment the node iterator by one.
                 * @pre @a it >= 0 && it < Number of nodes in the graph
                 * @return A node iterator object that has been incremented by 1
                 * @post @a it = @a it + 1
                 *
                 * Complexity: O(1)
                 */
                node_iterator& operator++() {
                        this->graph_pos_++;
                        return(*this);
                }

                /** Test whether this NodeIterator and @a it are equal.
                 *
                 * Equal node iterators are incident to the same graph and are at the same position in the graph.
                 */
                bool operator==(const node_iterator& it) const {
                        return(it.node_iterator_graph_ == this->node_iterator_graph_ && it.graph_pos_ == this->graph_pos_);
                }

        private:
                // Allow Graph to access NodeIterator's private member data and functions.
                friend class Graph;

                // pointer to the graph of this node iterator
                Graph* node_iterator_graph_;

                // position that the iterator currently has
                size_type graph_pos_;

                /** Construct a node iterator object.
                 * @param[in] graph     pointer to the graph of the iterator
                 * @param[in] position  the position that the iterator starts with (default = 0)
                 */
                NodeIterator(const Graph* graph, const size_type position = 0)
                : node_iterator_graph_(const_cast<Graph*>(graph)), graph_pos_(position) {
                }
        };

        /** Return a node iterator to the start of the graph nodes */
        node_iterator node_begin() const {
                return(NodeIterator(this));
        }

        /** Return a node iterator to the end of the graph nodes */
        node_iterator node_end() const {
                return(NodeIterator(this, this->graph_nodes_.size()));
        }

        //
        // Incident Iterator
        //

        /** @class Graph::IncidentIterator
         * @brief Iterator class for edges incident to a node. A forward iterator.
         */
        class IncidentIterator : private totally_ordered<IncidentIterator> {
        public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator() {
                        // no code needed here
                }

                /** Derefrence the incident iterator to get the edge of the current node.
                 * @pre @a it >= 0 && it < Number of nodes in the graph
                 * @return The edge that is at the current loaction adjacent to the iterators root node
                 *
                 * Has a possibility of returning a edge that is not adjacent to the node if
                 *    @a it >= Number of nodes in the graph.
                 *
                 * Complexity: O(1)
                 */
                Edge operator*() const {
                        size_type edge_id = this->incident_iterator_graph_->graph_nodes_[this->og_node_idx_].adj_edges[this->inc_it_pos_];
                        size_type node_a = this->og_node_idx_;
                        size_type node_b = this->incident_iterator_graph_->graph_edges_[edge_id].edge_node_b;
                        if(node_a == node_b) {
                                node_b = this->incident_iterator_graph_->graph_edges_[edge_id].edge_node_a;
                        }
                        return(Edge(this->incident_iterator_graph_, edge_id, node_a, node_b));
                }

                /** Increment the incident iterator by one.
                 * @pre @a it >= 0 && it < Number of nodes in the graph
                 * @return A incident iterator object that has been incremented by 1
                 * @post @a it = @a it + 1
                 *
                 * Complexity: O(1)
                 */
                IncidentIterator& operator++() {
                        this->inc_it_pos_++;
                        return(*this);
                }

                /** Test whether this IncidentIterator and @a it are equal.
                 *
                 * Equal incident iterators are incident to the same graph,
                 *   have the same root node,
                 *   and are at the same position
                 */
                bool operator==(const IncidentIterator& it) const {
                        return(it.incident_iterator_graph_ == this->incident_iterator_graph_
                                   && it.inc_it_pos_ == this->inc_it_pos_
                                   && this->og_node_idx_ == it.og_node_idx_);
                }

        private:
                // Allow Graph to access IncidentIterator's private member data and functions.
                friend class Graph;

                // pointer to the graph of incident iterator
                graph_type* incident_iterator_graph_;

                // the id of the node (root) we are iterating around
                size_type og_node_idx_;

                // the position that the iterator currently has
                size_type inc_it_pos_;

                /** Construct a incident iterator object.
                 * @param[in] graph     pointer to the graph of the iterator
                 * @param[in] node_idx  id of the root node
                 * @param[in] position  the position that the iterator starts with (default = 0)
                 */
                IncidentIterator(const Graph* graph, const size_type node_idx, const size_type position = 0)
                : incident_iterator_graph_(const_cast<Graph*>(graph)), og_node_idx_(node_idx), inc_it_pos_(position) {
                }
        };

        //
        // Edge Iterator
        //

        /** @class Graph::EdgeIterator
         * @brief Iterator class for edges. A forward iterator.
         */
        class EdgeIterator : private totally_ordered<EdgeIterator> {
        public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid EdgeIterator. */
                EdgeIterator() {
                        // no code needed here
                }

                /** Derefrence the edge iterator to get the edge.
                 * @pre @a it >= 0 && it < Number of edges in the graph
                 * @return The edge that is at the current loaction in the graph
                 *
                 * Has a possibility of returning a edge that is not in the graph if
                 *    @a it >= Number of nodes in the graph.
                 *
                 * Complexity: O(1)
                 */
                Edge operator*() const {
                        size_type node_a = this->edge_iterator_graph_->graph_edges_[this->graph_pos_].edge_node_a;
                        size_type node_b = this->edge_iterator_graph_->graph_edges_[this->graph_pos_].edge_node_b;
                        return(Edge(this->edge_iterator_graph_, this->graph_pos_, node_a, node_b));
                }

                /** Increment the edge iterator by one.
                 * @pre @a it >= 0 && it < Number of edges in the graph
                 * @return A edge iterator object that has been incremented by 1
                 * @post @a it = @a it + 1
                 *
                 * Complexity: O(1)
                 */
                EdgeIterator& operator++() {
                        this->graph_pos_++;
                        return(*this);
                }

                /** Test whether this EdgeIterator and @a it are equal.
                 *
                 * Equal edge iterators are incident to the same graph,
                 *   and are at the same position
                 */
                bool operator==(const EdgeIterator& it) const {
                        return(this->edge_iterator_graph_ == it.edge_iterator_graph_ && this->graph_pos_ == it.graph_pos_);
                }

        private:
                // Allow Graph to access EdgeIterator's private member data and functions.
                friend class Graph;

                // pointer to the graph of the edge iterator
                Graph* edge_iterator_graph_;

                // position in the graph of this edge iterator
                size_type graph_pos_;

                /** Construct a edge iterator object.
                 * @param[in] graph     pointer to the graph of the iterator
                 * @param[in] position  the position that the iterator starts with (default = 0)
                 */
                EdgeIterator(const Graph* graph, const size_type position = 0)
                : edge_iterator_graph_(const_cast<Graph*>(graph)), graph_pos_(position) {
                }
        };

        /** Return an edge iterator to the start of the graph edges */
        edge_iterator edge_begin() const {
                return(EdgeIterator(this));
        }

        /** Return an edge iterator to the end of the graph edges */
        edge_iterator edge_end() const {
                return(EdgeIterator(this, this->graph_edges_.size()));
        }

private:

        /** Struct that represents all of the information of a node */
        struct proxy_node {
                // positon of the node in Euclidean Space
                Point pos;

                // the value that is stored in the nodes
                node_value_type node_val;

                // vector of all the edge ids that are adjacent to the node
                std::vector<size_type> adj_edges;

                /** Construct a proxy node.
                 * @param[in] position  the position of the node in 3D Euclidean Space
                 * @param[in] value     the value that the node stores
                 */
                proxy_node(const Point &position, const node_value_type value)
                : pos(position), node_val(value), adj_edges() {
                }
        };

        /** Struct that represents all of the information of an edge */
        struct proxy_edge {
                // id of the root node of the edge
                size_type edge_node_a;

                // id of the leaf node of the edge
                size_type edge_node_b;

                // value of the node
                edge_value_type edge_val;

                /** Construct a proxy node.
                 * @param[in] a     id of the root node of the edge
                 * @param[in] b     id of the leaf node of the edge
                 */
                proxy_edge(const size_type a, const size_type b, const edge_value_type& val)
                : edge_node_a(a), edge_node_b(b), edge_val(val) {
                }
        };

        // vector that stores all of the proxy nodes of the graph
        std::vector<proxy_node> graph_nodes_;

        // vector that stores all of the proxy edges of the graph
        std::vector<proxy_edge> graph_edges_;

};

#endif // CME212_GRAPH_HPP
