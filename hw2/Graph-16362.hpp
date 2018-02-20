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
	private:
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
		* Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		*/
		size_type num_edges() const {
			return _num_edges;
		}

		/** Return the edge with index @a i.
		* @pre 0 <= @a i < num_edges()
		*
		* Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		*/
		Edge edge(size_type i) const {
			return Edge(this, i);
		}

		/** Test whether two nodes are connected by an edge.
		* @pre @a a and @a b are valid nodes of this graph
		* @return True if for some @a i, edge(@a i) connects @a a and @a b.
		*
		* Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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
};

#endif // CME212_GRAPH_HPP