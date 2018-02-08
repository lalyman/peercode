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
				Node() {}

				/** Return this node's position. */
				const Point& position() const {
					return _graph->_points[_uid];
				}

				/** Return this node's index, a number in the range [0, graph_size). */
				size_type index() const {
					return _uid;
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
					return (_uid >= n._uid) ? false : true;
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
		* @post new num_nodes() == old num_nodes() + 1
		* @post result_node.index() == old num_nodes()
		*
		* Complexity: O(1) amortized operations.
		*/
		Node add_node(const Point& position) {
			_points.push_back(position);
			++_num_nodes;
			return Node(this, _num_nodes - 1);
		}

		/** Determine if a Node belongs to this Graph
		* @return True if @a n is currently a Node of this Graph
		*
		* Complexity: O(1).
		*/
		bool has_node(const Node& n) const {
			return (n._uid < _num_nodes) ? true : false;
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
		class Edge {
			public:
				/** Construct an invalid Edge. */
				Edge() {}

				/** Return a node of this Edge */
				Node node1() const {
					return Node(_graph, _graph->_edge_node1[_uid]);
				}

				/** Return the other node of this Edge */
				Node node2() const {
					return Node(_graph, _graph->_edge_node2[_uid]);
				}

				/** Test whether this edge and @a e are equal.
				 *
				 * Equal edges represent the same undirected edge between two nodes.
				 */
				bool operator==(const Edge& e) const {
					return ((e._graph != _graph) || (e._uid != _uid)) ? false : true;
				}

				/** Test whether this edge is less than @a e in a global order.
				 *
				 * This ordering function is useful for STL containers such as
				 * std::map<>. It need not have any interpretive meaning.
				 */
				bool operator<(const Edge& e) const {
					return (_uid >= e._uid) ? false : true;
				}

			private:
				// Allow Graph to access Edge's private member data and functions.
				friend class Graph;

				/** Private Constructor */
				Edge(const Graph* graph, size_type uid)
					: _graph(const_cast<Graph*>(graph)), _uid(uid) {}

				// Pointer back to the Graph container
				Graph* _graph;
				// This edge's unique identification number
				size_type _uid;
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
			size_type a_index = a.index();
			size_type b_index = b.index();

			for (size_type i = 0; i < _num_edges; ++i) {
				if ((a_index == _edge_node1[i]) && (b_index == _edge_node2[i])) {
					return true;
				}
				if ((a_index == _edge_node2[i]) && (b_index == _edge_node1[i])) {
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
		Edge add_edge(const Node& a, const Node& b) {
			size_type a_index = a.index();
			size_type b_index = b.index();

			/** If edge exists, return existing edge. */
			for (size_type i = 0; i < _num_edges; ++i) {
				if ((a_index == _edge_node1[i]) && (b_index == _edge_node2[i])) {
					return Edge(this, i);
					break;
				}
				if ((a_index == _edge_node2[i]) && (b_index == _edge_node1[i])) {
					return Edge(this, i);
					break;
				}
			}

			/** If edge does not exist, create a new edge. */
			_edge_node1.push_back(a_index);
			_edge_node2.push_back(b_index);
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
			_edge_node1.clear();
			_edge_node2.clear();
			_points.shrink_to_fit();
			_edge_node1.shrink_to_fit();
			_edge_node2.shrink_to_fit();

			/** Set number of nodes and edges to zero. */
			_num_nodes = 0;
			_num_edges = 0;
		}

	private:
		/** Number of nodes in graph. */
		size_type _num_nodes;
		/** Number of edges in graph. */
		size_type _num_edges;
		/** Vector containing all the points corresponding to the nodes in graph. */
		std::vector<Point> _points;
		/** Vector containing node 1 of each edge in graph. */
		std::vector<size_type> _edge_node1;
		/** Vector containing node 2 of each edge in graph. */
		std::vector<size_type> _edge_node2;
};

#endif // CME212_GRAPH_HPP
