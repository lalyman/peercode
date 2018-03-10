#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

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

		/** Synonyms for graph template types (following STL conventions). */
		typedef V node_value_type;
		typedef E edge_value_type;

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
		Graph::num_edges(), and argument type of Graph::node(size_type) */
		using size_type = unsigned;

		//
		// CONSTRUCTORS AND DESTRUCTOR
		//

		/** Construct an empty graph. */
		Graph() : _num_nodes(0), _num_edges(0) {}

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
				Node() {}

				/** Return this node's position as reference. */
				Point& position() {
					return _graph->_points[_uid];
				}

				/** Return this node's position as const reference. */
				const Point& position() const {
					return _graph->_points[_uid];
				}

				/** Return this node's index, a number in the range [0, graph_size). */
				size_type index() const {
					return _uid;
				}

				/** Return node value as reference. */
				node_value_type& value() {
					return _graph->_values[_uid];
				}

				/** Return node value as const reference. */
				const node_value_type& value() const {
					return _graph->_values[_uid];
				}

				/** Return node degree. */
				size_type degree() const {
					return _graph->_adj_edges[_uid].size();
				}

				/** Return an iterator for the node's first incident edge. */
				IncidentIterator edge_begin() const {
					return IncidentIterator(_graph, _uid, 0);
				}

				/** Return an iterator for the past-the-end element of the node. */
				IncidentIterator edge_end() const {
					return IncidentIterator(_graph, _uid, degree());
				}


				/** Test whether this node and @a n are equal.
				 *
				 * Equal nodes have the same graph and the same index.
				 */
				bool operator==(const Node& n) const {
					return ((n._graph != _graph) || (n._uid != _uid)) ? false : true;
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
					if (n._graph == _graph) {
						return (_uid < n._uid);
					}
					else {
						return (_graph < n._graph);
					}
				}

			private:

				// Allow Graph to access Node's private member data and functions.
				friend class Graph;
    
				/** Private Constructor */
				Node(const Graph* graph, size_type uid) 
					: _graph(const_cast<Graph*>(graph)), _uid(uid) {}

				// Pointer back to the Graph container
				Graph* _graph;
				// This node's unique identification number
				size_type _uid;
		};
		
		/** Return the number of nodes in the graph.
		*
		* Complexity: O(1).
		*/
		size_type size() const {
			return _num_nodes;
		}

		/** Synonym for size(). */
		size_type num_nodes() const {
			return size();
		}
		
		/** Add a node to the graph, returning the added node.
		* @param[in] position The new node's position
		* @param[in] value The new node's value
		* @post new num_nodes() == old num_nodes() + 1
		* @post result_node.index() == old num_nodes()
		*
		* Complexity: O(1) amortized operations.
		*/
		Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
			_points.push_back(position);
			_values.push_back(value);
			_adj_edges[_num_nodes];
			++_num_nodes;
			return Node(this, _num_nodes - 1);
		}

		/** Determine if a Node belongs to this Graph
		* @return True if @a n is currently a Node of this Graph
		*
		* Complexity: O(1).
		*/
		bool has_node(const Node& n) const {
			return ((this == n._graph) && (n._uid < _num_nodes)) ? true : false;
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

		/** Remove the node if it exists.
		* @param[in] n      A valid node of the graph
		* @return           1 if node was removed, 0 otherwise
		*
		* @pre  @a n is a valid node of the graph.
		* @post @a n not a node of the graph.
		*		- new @a _num_nodes = old @a _num_nodes - 1
		*       - new @a _num_edges = old @a _num_edges - degree(n) .
		*
		* Can invalidate node and edge indexes, and node and edge iterators.
		*
		* Complexity : O(max degree) .
		*/
		size_type remove_node(const Node& n) {
			if (has_node(n)) {
				// Get node id
				size_type node_id = n._uid;

				// Swap node with last node
				swap_nodes(node_id, _num_nodes - 1);

				// Pop last node
				pop_node();

				return 1;
			}
			return 0;
		}

		/** Remove the node referenced by the node iterator if it exists.
		* @param[in] n_it   A valid node iterator of the graph
		* @return           A node iterator that points to a node @a n1 such that
		*                   old node index of @a n1 = old @a _num_edges - 1 .
		*
		* @pre  Let @a n be a valid node of the graph to which iterator @a n_it points to.
		* @post @a n not a node of the graph.
		*		- new @a _num_nodes = old @a _num_nodes - 1
		*       - new @a _num_edges = old @a _num_edges - degree(n) .
		*
		* Can invalidate node and edge indexes, and node and edge iterators.
		*
		* Complexity : O(max degree) .
		*/
		NodeIterator remove_node(NodeIterator n_it) {
			remove_node(*n_it);
			return n_it;
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
				Edge() {}

				/** Return a node of this Edge */
				Node node1() const {
					if (_override) {
						return Node(_graph, _graph->_edge_node2[_uid]);
					}
					else {
						return Node(_graph, _graph->_edge_node1[_uid]);
					}
				}

				/** Return the other node of this Edge */
				Node node2() const {
					if (_override) {
						return Node(_graph, _graph->_edge_node1[_uid]);
					}
					else {
						return Node(_graph, _graph->_edge_node2[_uid]);
					}
				}

				/** Return edge value as reference. */
				edge_value_type& value() {
					return _graph->_edge_values[_uid];
				}

				/** Return edge value as const reference. */
				const edge_value_type& value() const {
					return _graph->_edge_values[_uid];
				}

				/** Test whether this edge and @a e are equal.
				 *
				 * Equal edges represent the same undirected edge between two nodes.
				 */
				bool operator==(const Edge& e) const {
					if (e._graph == _graph) {
						if ((e.node1() == node1()) && (e.node2() == node2())) {
							return true;
						}
						if ((e.node1() == node2()) && (e.node2() == node1())) {
							return true;
						}
					}
					return false;
				}
				
				/** Test whether this edge is less than @a e in a global order.
				 *
				 * This ordering function is useful for STL containers such as
				 * std::map<>. It need not have any interpretive meaning.
				 */
				bool operator<(const Edge& e) const {
					if (e._graph == _graph) {
						return (_uid < e._uid);
					}
					else {
						return (_graph < e._graph);
					}
				}

			private:
				// Allow Graph to access Edge's private member data and functions.
				friend class Graph;

				/** Private Constructors */
				Edge(const Graph* graph, size_type uid)
					: _graph(const_cast<Graph*>(graph)), _uid(uid), _override(false) {}

				Edge(const Graph* graph, size_type uid, bool override)
					: _graph(const_cast<Graph*>(graph)), _uid(uid), _override(override) {}

				// Pointer back to the Graph container
				Graph* _graph;
				// This edge's unique identification number
				size_type _uid;
				// Inverts node1 and node2 if true
				bool _override;
		};

		/** Return the total number of edges in the graph.
		*
		* Complexity: O(1) .
		*/
		size_type num_edges() const {
			return _num_edges;
		}

		/** Return the edge with index @a i.
		* @pre 0 <= @a i < num_edges()
		*
		* Complexity: O(1) .
		*/
		Edge edge(size_type i) const {
			return Edge(this, i);
		}

		/** Test whether two nodes are connected by an edge.
		* @param[in] a         A valid node
		* @param[in] b         A valid node
		* @return True if for some @a i, edge(@a i) connects @a a and @a b.
		*
		* @pre @a a and @a b are valid nodes
		* 
		* Complexity: O(max degree)
		*/
		bool has_edge(const Node& a, const Node& b) const {

			if (has_node(a) && has_node(b)) {
				if (a.degree() < b.degree()) {
					for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
						if ((*ei).node2() == b) return true;
					}
				}
				else {
					for (auto ei = b.edge_begin(); ei != b.edge_end(); ++ei) {
						if ((*ei).node2() == a) return true;
					}
				}
			}

			return false;
		}

		/** Add an edge to the graph, or return the current edge if it already exists.
		* @param[in] a         A valid node of the graph
		* @param[in] b         A valid node of the graph
		* @param[in] value     Edge value
		* @return An Edge object e with e.node1() == @a a and e.node2() == @a b
		* 
		* @pre @a a and @a b are distinct valid nodes of this graph
		* @post has_edge(@a a, @a b) == true
		* @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
		*       Else,                        new num_edges() == old num_edges() + 1.
		*
		* Can invalidate edge indexes -- in other words, old edge(@a i) might not
		* equal new edge(@a i). Must not invalidate outstanding Edge objects.
		*
		* Complexity: O(max degree) .
		*/
		Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {

			if (!has_node(a) || !has_node(b)) {
				assert(1 == 2);
			}

			/** If edge exists, return existing edge. 
			Update edge values.
			*/
			if (a.degree() < b.degree()) {
				for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
					Edge e = *ei;
					if (e.node2() == b) {
						_edge_values[e._uid] = value;
						return e;
					}
				}
			}
			else {
				for (auto ei = b.edge_begin(); ei != b.edge_end(); ++ei) {
					Edge e = *ei;
					if (e.node2() == a) {
						_edge_values[e._uid] = value;
						return Edge(this, e._uid, true);
					}
				}
			}

			/** If edge does not exist, create a new edge. */
			_edge_node1.push_back(a.index());
			_edge_node2.push_back(b.index());
			_edge_values.push_back(value);

			_adj_edges[a.index()].push_back(_num_edges);
			_adj_edges[b.index()].push_back(_num_edges);

			++_num_edges;
			return Edge(this, _num_edges - 1);
		}

		/** Remove the edge if it exists.
		* @param[in] a      A valid node of the graph
		* @param[in] b      A valid node of the graph
		* @return           1 if edge was removed, 0 otherwise
		*
		* @pre  @a a and @a b are valid nodes of the graph.
		* @post If an edge @a e existed between @a a and @a b : 
		*       - @a e not an edge of the graph.
		*		- new @a _num_edges = old @a _num_edges - 1 .
		*
		* Can invalidate edge indexes, and edge iterators.
		*
		* Complexity : O(1) .
		*/
		size_type remove_edge(const Node& a, const Node& b) {
			// Test if edge exists and if yes get it
			if (has_edge(a, b)) {
				// Overwriting edge value does not matter as edge will be removed
				Edge e = add_edge(a, b);

				// Get edge id
				size_type edge_id = e._uid;

				// Swap edge with last edge
				swap_edges(edge_id, _num_edges - 1);

				// Pop last edge
				pop_edge();

				return 1;
			}
			return 0;
		}

		/** Remove the edge if it exists.
		* @param[in] e      A valid edge of the graph
		* @return           1 if edge was removed, 0 otherwise
		* 
		* @pre  @a e is a valid edge of the graph.
		* @post @a e not an edge of the graph.
		*		new @a _num_edges = old @a _num_edges - 1 .
		*
		* Can invalidate edge indexes, and edge iterators.
		*
		* Complexity : O(1) .
		*/
		size_type remove_edge(const Edge& e) {
			return remove_edge(e.node1(), e.node2());
		}

		/** Remove the edge referenced by the edge iterator if it exists.
		* @param[in] e_it   A valid edge iterator of the graph
		* @return           An edge iterator that points to an edge @a e1 such that
		*                   old edge index of @a e1 = old @a _num_edges - 1 .
		* 
		* @pre  Let @a e be a valid edge of the graph to which iterator @a e_it points to.
		* @post @a e not an edge of the graph.
		*		new @a _num_edges = old @a _num_edges - 1 .
		*
		* Can invalidate edge indexes, and edge iterators.
		*
		* Complexity : O(1) .
		*/
		EdgeIterator remove_edge(EdgeIterator e_it) {
			remove_edge(*e_it);
			return e_it;
		}

		/** Remove all nodes and edges from this graph.
		* @post num_nodes() == 0 && num_edges() == 0
		*
		* Invalidates all outstanding Node and Edge objects.
		*/
		void clear() {
			/** Clear all the vectors and releases storage. */
			_points.clear();
			_values.clear();
			_adj_edges.clear();
			_edge_node1.clear();
			_edge_node2.clear();
			_edge_values.clear();

			_points.shrink_to_fit();
			_values.shrink_to_fit();
			_edge_node1.shrink_to_fit();
			_edge_node2.shrink_to_fit();
			_edge_values.shrink_to_fit();

			/** Set number of nodes, edges, unique edges to zero. */
			_num_nodes = 0;
			_num_edges = 0;
		}

		//
		// Node Iterator
		//

		/** @class Graph::NodeIterator
		* @brief Iterator class for nodes. A forward iterator. 
		*/
		class NodeIterator : private equality_comparable<NodeIterator> {
			public:
				// These type definitions let us use STL's iterator_traits.
				using value_type = Node;                            // Element type
				using pointer = Node*;                              // Pointers to elements
				using reference = Node&;                            // Reference to elements
				using difference_type = std::ptrdiff_t;             // Signed difference
				using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

				/** Construct an invalid NodeIterator. */
				NodeIterator() {}

				/** Return a node that the iterator refers to. */
				Node operator*() const {
					return _graph->node(_uid);
				}

				/** Increment the iterator. */
				NodeIterator& operator++() {
					++_uid;
					return *this;
				}

				/** Return a boolean indicating whether two iterators are equal. */
				bool operator==(const NodeIterator& iterator) const {
					return ((_graph == iterator._graph) && (_uid == iterator._uid));
				}

			private:
				friend class Graph;
				/** Private Constructor */
				NodeIterator(const Graph* graph, size_type uid)
					: _graph(const_cast<Graph*>(graph)), _uid(uid) {
				}
				// Pointer back to the Graph container
				Graph* _graph;
				// The node unique identification number to which iterator points to
				size_type _uid;
		};

		/** Return an iterator for the first node in the graph. */
		NodeIterator node_begin() const {
			return NodeIterator(this, 0);
		}

		/** Return an iterator for the past-the-end element of the nodes in the graph. */
		NodeIterator node_end() const {
			return NodeIterator(this, _num_nodes);
		}

		//
		// Incident Iterator
		//

		/** @class Graph::IncidentIterator
		* @brief Iterator class for edges incident to a node. A forward iterator. 
		*/
		class IncidentIterator : private equality_comparable<IncidentIterator> {
			public:
				// These type definitions let us use STL's iterator_traits.
				using value_type = Edge;                            // Element type
				using pointer = Edge*;                              // Pointers to elements
				using reference = Edge&;                            // Reference to elements
				using difference_type = std::ptrdiff_t;             // Signed difference
				using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

				/** Construct an invalid IncidentIterator. */
				IncidentIterator() {}

				/** Return an edge that the iterator refers to. */
				Edge operator*() const {
					size_type edge_id = _graph->_adj_edges[_node_id][_adj_edge_id];

					if (_graph->edge(edge_id).node1().index() == _node_id) {
						return Edge(_graph, edge_id);
					}
					else {
						return Edge(_graph, edge_id, true);
					}
				}

				/** Return the next iterator for the incident edge iterator. */
				IncidentIterator& operator++() {
					++_adj_edge_id;
					return *this;
				}

				/** Return a boolean indicating whether two iterators are equal. */
				bool operator==(const IncidentIterator& iterator) const {
					return ((_graph == iterator._graph) && 
						(_node_id == iterator._node_id) && 
						(_adj_edge_id == iterator._adj_edge_id));
				}


			private:
				friend class Graph;
				/** Private Constructor */
				IncidentIterator(const Graph* graph, size_type node_id, size_type adj_edge_id)
					: _graph(const_cast<Graph*>(graph)), _node_id(node_id), _adj_edge_id(adj_edge_id) {
				}
				// Pointer back to the Graph container
				Graph* _graph;
				// The node unique identification number to which iterator points to
				size_type _node_id;
				// The adjacent edge id to which iterator points to
				size_type _adj_edge_id;
		};

		//
		// Edge Iterator
		//

		/** @class Graph::EdgeIterator
		* @brief Iterator class for edges. A forward iterator. 
		*/
		class EdgeIterator : private equality_comparable<EdgeIterator> {
			public:
				// These type definitions let us use STL's iterator_traits.
				using value_type = Edge;                            // Element type
				using pointer = Edge*;                              // Pointers to elements
				using reference = Edge&;                            // Reference to elements
				using difference_type = std::ptrdiff_t;             // Signed difference
				using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

				/** Construct an invalid EdgeIterator. */
				EdgeIterator() {}

				/** Return an edge that the iterator refers to. */
				Edge operator*() const {
					return _graph->edge(_uid);
				}

				/** Increment the iterator. */
				EdgeIterator& operator++() {
					++_uid;
					return *this;
				}

				/** Return a boolean indicating whether two iterators are equal. */
				bool operator==(const EdgeIterator& iterator) const {
					return ((_graph == iterator._graph) && (_uid == iterator._uid));
				}

			private:
				friend class Graph;
				/** Private Constructor */
				EdgeIterator(const Graph* graph, size_type uid)
					: _graph(const_cast<Graph*>(graph)), _uid(uid) {
				}
				// Pointer back to the Graph container
				Graph* _graph;
				// The edge unique identification number to which iterator points to
				size_type _uid;
		};

		/** Return an iterator for the first edge in the graph. */
		EdgeIterator edge_begin() const {
			return EdgeIterator(this, 0);
		}

		/** Return an iterator for the past-the-end element of the edges in the graph. */
		EdgeIterator edge_end() const {
			return EdgeIterator(this, _num_edges);
		}

	private:
		/** Number of nodes, edges in graph. */
		size_type _num_nodes;
		size_type _num_edges;

		/** 
		Vector containing all the points corresponding to the nodes in graph.
		Vector containing all the values corresponding to the nodes in graph.
		Dictionary of adjacent edges incident to a node in graph.
		*/
		std::vector<Point> _points;
		std::vector<node_value_type> _values;
		std::map<size_type, std::vector<size_type>> _adj_edges;

		/** 
		Vectors containing node 1 and node 2 of each unique edge in graph.
		Vector containing all the values corresponding to the unique edges in graph.
		*/
		std::vector<size_type> _edge_node1;
		std::vector<size_type> _edge_node2;
		std::vector<edge_value_type> _edge_values;

	private:

		/** Swap two nodes in the graph.
		* @param[in] n1   Valid node index of graph
		* @param[in] n2   Valid node index of graph
		* 
		* @pre  Let @a v1 and @a v2 be valid old nodes of the graph such that
		*       old @a v1 index = @a n1 , old @a v2 index = @a n2 .
		* @post new @a v1 index = @a n2 , new @a v2 index = @a n1 .
		* 
		* Can invalidate edge and node indexes, and node and edge iterators.
		* 
		* Complexity : O(max degree) assuming graph is sparse.
		*/
		void swap_nodes(size_type n1, size_type n2) {
			if ((n1 < num_nodes()) && (n2 < num_nodes())) {
				// Fix nodes of edges involved
				Node node1 = node(n1);
				for (auto it = node1.edge_begin(); it != node1.edge_end(); ++it) {
					size_type edge_id = (*it)._uid;
					if ((_edge_node1[edge_id] == n1) && (_edge_node2[edge_id] != n2)) {
						_edge_node1[edge_id] = n2;
					}
					if ((_edge_node2[edge_id] == n1) && (_edge_node1[edge_id] != n2)) {
						_edge_node2[edge_id] = n2;
					}
				}
				Node node2 = node(n2);
				for (auto it = node2.edge_begin(); it != node2.edge_end(); ++it) {
					size_type edge_id = (*it)._uid;
					if ((_edge_node1[edge_id] == n2) && (_edge_node2[edge_id] != n1)) {
						_edge_node1[edge_id] = n1;
					}
					if ((_edge_node2[edge_id] == n2) && (_edge_node1[edge_id] != n1)) {
						_edge_node2[edge_id] = n1;
					}
				}

				// Swap points and values
				Point p = _points[n1];
				node_value_type node_val = _values[n1];
				_points[n1] = _points[n2];
				_points[n2] = p;
				_values[n1] = _values[n2];
				_values[n2] = node_val;

				// Swap adjacency lists
				std::vector<size_type> adj_edges = _adj_edges[n1];
				_adj_edges[n1] = _adj_edges[n2];
				_adj_edges[n2] = adj_edges;
			}
		}

		/** Pop last node in the graph. 
		*
		* @pre  Let @a v1 be node of the graph such that
		*       old @a v1 index = @a _num_nodes - 1 .
		* @post v1 not a node of graph.
		*		new @a _num_nodes = old @a _num_nodes - 1 
		*       new @a _num_edges = old @a _num_edges - degree(v1)
		*
		* Can invalidate edge and node indexes, and node and edge iterators.
		*
		* Complexity : O(max degree) assuming graph is sparse.
		*/
		void pop_node() {
			// Pop node and value arrays
			_points.pop_back();
			_values.pop_back();

			// Delete edges incident on the node
			while (_adj_edges[_num_nodes - 1].size() > 0) {
				size_type edge_id = _adj_edges[_num_nodes - 1][0];
				swap_edges(edge_id, _num_edges - 1);
				pop_edge();
			}

			// Delete adjacency list
			_adj_edges.erase(_num_nodes - 1);

			// Update number of nodes
			--_num_nodes;
		}

		/** Swap two edges in the graph. 
		* @param[in] n1   Valid edge index of graph
		* @param[in] n2   Valid edge index of graph
		*
		* @pre  Let @a e1 and @a e2 be valid old edges of the graph such that
		*       old @a e1 edge index = @a n1 , old @a e2 edge index = @a n2 .
		* @post new @a e1 edge index = @a n2 , new @a e2 edge index = @a n1 .
		*
		* Can invalidate edge indexes, and edge iterators.
		*
		* Complexity : O(1).
		*/
		void swap_edges(size_type n1, size_type n2) {
			if ((n1 < num_edges()) && (n2 < num_edges())) {
				// Update adjacency lists
				size_type node1_id = _edge_node1[n1];
				size_type node2_id = _edge_node2[n1];
				if ((node1_id != _edge_node1[n2]) && (node1_id != _edge_node2[n2])) {
					for (auto &it : _adj_edges[node1_id]) {
						if (it == n1) {
							it = n2;
							break;
						}
					}
				}
				if ((node2_id != _edge_node1[n2]) && (node2_id != _edge_node2[n2])) {
					for (auto &it : _adj_edges[node2_id]) {
						if (it == n1) {
							it = n2;
							break;
						}
					}
				}

				node1_id = _edge_node1[n2];
				node2_id = _edge_node2[n2];
				if ((node1_id != _edge_node1[n1]) && (node1_id != _edge_node2[n1])) {
					for (auto &it : _adj_edges[node1_id]) {
						if (it == n2) {
							it = n1;
							break;
						}
					}
				}
				if ((node2_id != _edge_node1[n1]) && (node2_id != _edge_node2[n1])) {
					for (auto &it : _adj_edges[node2_id]) {
						if (it == n2) {
							it = n1;
							break;
						}
					}
				}

				// Swap edge nodes and values
				node1_id = _edge_node1[n1];
				node2_id = _edge_node2[n1];
				edge_value_type edge_val = _edge_values[n1];
				_edge_node1[n1] = _edge_node1[n2];
				_edge_node1[n2] = node1_id;
				_edge_node2[n1] = _edge_node2[n2];
				_edge_node2[n2] = node2_id;
				_edge_values[n1] = _edge_values[n2];
				_edge_values[n2] = edge_val;
			}
		}

		/** Pop last edge in the graph. 
		*
		* @pre  Let @a e1 be node of the graph such that
		*       old @a e1 edge index = @a _num_edges -1 .
		* @post e1 not an edge of graph.
		*		new @a _num_edges = old @a _num_edges - 1 .
		*
		* Can invalidate edge indexes, and edge iterators.
		*
		* Complexity : O(1) .
		*/
		void pop_edge() {
			size_type node1_id = _edge_node1[_num_edges - 1];
			size_type node2_id = _edge_node2[_num_edges - 1];

			// Pop edge node and edge value arrays
			_edge_node1.pop_back();
			_edge_node2.pop_back();
			_edge_values.pop_back();
			
			// Update adjacency lists
			size_type t = 0;
			for (auto &it : _adj_edges[node1_id]) {
				if (it == (_num_edges - 1)) {
					break;
				}
				++t;
			}
			_adj_edges[node1_id][t] = _adj_edges[node1_id][_adj_edges[node1_id].size() - 1];
			_adj_edges[node1_id].pop_back();

			t = 0;
			for (auto &it : _adj_edges[node2_id]) {
				if (it == (_num_edges - 1)) {
					break;
				}
				++t;
			}
			_adj_edges[node2_id][t] = _adj_edges[node2_id][_adj_edges[node2_id].size() - 1];
			_adj_edges[node2_id].pop_back();

			// Update number of edges
			--_num_edges;
		}
};

#endif // CME212_GRAPH_HPP