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
using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
private:
    // Def of the adjacency list, a representation of a finite graph as a
    // collection of ordered lists.
    vector< vector< unsigned > > adjacencyList;
    
    struct edgeSt; // Def of an edge structure.
    vector< edgeSt > edgeVec; // Ordered vector of all edges, i.e., its heads and tails.
    
    struct vertexSt; // Def of a vertex structure, i.e, its location.
    vector< vertexSt > verticesVec; // Ordered vector of all vertices.
    
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
    // Initializes an empty adjacency list, vector of edges and vertices.
    Graph() : adjacencyList(), edgeVec(), verticesVec() {}
    
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
        Node() {} // Invalid node.
        
        /** Return this node's position. */
        // Getting position from a vector of vertices structs.
        const Point& position() const {
            return grid -> verticesVec[vertexIdx].location;
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        // Getting idx of a graph's vertex.
        size_type index() const {
            return vertexIdx;
        }
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        // Checking for equal idecies in a graph.
        bool operator==(const Node& n) const {
            bool valReturn = (vertexIdx == n.vertexIdx && grid == n.grid) ? true : false;
            return valReturn;
        }
        
        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        // Checking for equality of graphs using the global order.
        bool operator<(const Node& n) const {
            bool valReturn = ((grid != n.grid && grid < n.grid) ||
                              (grid == n.grid && grid < n.grid)) ? true : false;
            return valReturn;
        }
        
    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        
        // All require 8 bits of memory.
        Graph* grid; // Pointer to a graph.
        size_type vertexIdx; // Index of a vertex in the graph.
        size_type edgeIdx; // Index of an edge in the graph.
    };
    
    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    // Number of unique vertices in a graph.
    size_type size() const {
        return verticesVec.size();
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
        vertexSt newVertex;
        newVertex.location = position; // Def location of a vertex.
        verticesVec.push_back(newVertex); // Adding an extra node to the graph.
        // Pushing a placeholder vector to account for the new vector of vertices.
        adjacencyList.emplace_back(vector <size_type>());
        
        // Def a new vertex with current graph and vertex index.
        Node vertex;
        vertex.grid = const_cast<Graph*>(this);
        vertex.vertexIdx = verticesVec.size() - 1;
        return vertex;
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    // Checks num of nodes and if graphs are the same.
    bool has_node(const Node& n) const {
        bool varReturn;
        varReturn = (const_cast< Graph* >(this) == n.grid &&
                     n.vertexIdx < verticesVec.size()) ? true : false;
        return varReturn;
    }
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        Node vertexIth; // Def vertex at index i that the function returns.
        vertexIth.vertexIdx = i; // Def index of a vertex to be i.
        vertexIth.grid = const_cast< Graph* >(this); // Def graph for a given vertex.
        return vertexIth;
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
        Edge() {} // Invalid edge.
        
        /** Return a node of this Edge */
        // Def head vertex, its index and current graph.
        Node node1() const {
            Node head;
            head.vertexIdx = edgeHead;
            head.grid = const_cast< Graph* >(gridEdge);
            return head;
        }
        
        /** Return the other node of this Edge */
        // Def tail vertex, its index and current graph.
        Node node2() const {
            Node tail;
            tail.vertexIdx = edgeTail;
            tail.grid = const_cast< Graph* >(gridEdge);
            return tail;
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        // Checking equality of edges by testing max and min ind of vertices.
        bool operator==(const Edge& e) const {
            bool varReturn = (max(e.edgeHead, e.edgeTail) == max(edgeHead, edgeTail) &&
                              min(e.edgeHead, e.edgeTail) == min(edgeHead, edgeTail) &&
                              e.gridEdge == gridEdge) ? true : false;
            return varReturn;
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            bool varReturn = false;
            size_type minEdge = min(edgeHead, edgeTail);
            size_type minEEdge = min(e.edgeHead, e.edgeTail);
            
            // Using smaller memory location when edges are in different graphs.
            if( e.gridEdge > gridEdge) varReturn = true;
            // Using an index of a lower degree of a head vertex when edges are in the same
            // graph and index of a lower degree of a tail vertex if degrees are not equal.
            else if((minEdge < minEEdge) || (minEdge == minEEdge &&
                                             max(edgeHead, edgeTail) < max(e.edgeHead, e.edgeTail))) {
                varReturn = true;}
            return varReturn;
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        Graph* gridEdge; // Pointer to the graph acessible from this class.
        size_type edgeHead, edgeTail; // Head and tail edges each with 8 bits of memory.
        size_type edgeIdx; // Idx of an edge.
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    // Returning number of edges in a graph.
    size_type num_edges() const {
        return edgeVec.size();
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        Edge edgeIth; // Def edge at index i, the function returns.
        edgeSt newEdge = edgeVec[i]; // Accessing ith edge.
        // Def head and tail attributes of the edge object.
        edgeIth.edgeHead = newEdge.head;
        edgeIth.edgeTail = newEdge.tail;
        return edgeIth;
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    // Checking if an edge exists in a graph.
    bool has_edge(const Node& a, const Node& b) const {
        vector< size_type > vec = adjacencyList[a.vertexIdx];
        for( size_type j = 0; j < vec.size(); j++) {
            if( vec[j] == b.vertexIdx) return true; }
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
        // Creating an edge between vertices a and b.
        if(!has_edge(a, b)) {
            edgeSt edgeRef;
            edgeRef.tail = b.vertexIdx;
            edgeRef.head = a.vertexIdx;
            edgeVec.push_back(edgeRef);
            
            // Adding edge into the adjList.
            adjacencyList[a.vertexIdx].push_back(b.vertexIdx);
            adjacencyList[b.vertexIdx].push_back(a.vertexIdx);
        }
        
        // Creating an Edge object that will be returned.
        Edge newEdge;
        newEdge.edgeHead = a.vertexIdx; 
        newEdge.edgeTail = b.vertexIdx;
        newEdge.gridEdge = const_cast< Graph* >(this);
        newEdge.edgeIdx = edgeVec.size() - 1;
        
        return newEdge;
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Removes all elements of the vertices and edges vectors and adjacencyList.
        verticesVec.clear();
        edgeVec.clear();
        adjacencyList.clear();
    }
    
private:
    
    // Def of edge's head and tail.
    struct edgeSt {
        size_type head;
        size_type tail;
    };
    
    // Def of a vertex's structure, i.e, its location.
    struct vertexSt {
        Point location;
    };
    
};

#endif // CME212_GRAPH_HPP
