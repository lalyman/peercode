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
template <typename V>
class Graph {
 private:
    
    struct edgeSt; // Def of an edge structure.
    std::vector< edgeSt > edgeVec; // Ordered vector of all edges, i.e., its heads and tails.
    
    struct vertexSt; // Def of a vertex structure, i.e, its location.
    std::vector< vertexSt > verticesVec; // Ordered vector of all vertices.
    
    //track total number of edges
    unsigned total_edges = 0;
    
    // Def of the adjacency list, a representation of a finite graph as a
    // collection of ordered lists.
    std::vector< std::vector< unsigned > > adjacencyList;
    
 public:
    
    //
    // PUBLIC TYPE DEFINITIONS
    //
    
    /** Type of this graph. */
    using graph_type = Graph;
    
    /** Type of graph's input. */
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
    
    /** Type of indexes and sizes.
     Return type of Graph::Node::index(), Graph::num_nodes(),
     Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;
    
    //
    // CONSTRUCTORS AND DESTRUCTOR
    //
    
    /** Construct an empty graph. */
    Graph() : verticesVec(), adjacencyList() {}
    
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
        Node() {} // Invalid node.
        
        /* Getting node's position.
         *
         * position() const accessess vertex's position from the adjacency list.
         * @return @a location by reference to Point.
         *
         * @pre function must be called on a valid node.
         */
        const Point& position() const {
            return grid->verticesVec[vertexIdx].location;
        }
        
        /* Getting idx of a graph's vertex.
         *
         * index() const gets the idx.
         * @return @a vertexIdx of type unsigned. 
         *
         * @pre function must be called on a valid node.
         * @post @a vertexIdx is a number in the range [0, graph_size).
         */
        size_type index() const {
            return vertexIdx;
        }
        
        /* Adds new information to a vertex.
         *
         * value() allows to add additional info regarding the vertex.
         * @return @a val by reference to Node of type V.
         *
         * @pre function must be called on a valid node.
         */
        node_value_type& value() {
            return grid->verticesVec[vertexIdx].val;
        }
        
        /** Add new information to a vertex.
         *
         * value() const allows to add additional info regardin a node the vertex.
         * @return const @a val by reference for node of type V. @a val is of type const.
         *
         * @pre graph pointer must point to a valid memeory location storing graph object.
         */
        const node_value_type& value() const {
            return grid->verticesVec[vertexIdx].val;
        }
        
        /** Gives a number of neighboord a vertex has.
         *
         * degree() const number of adacent vertices.
         * @return @a d of type size_type.
         *
         * @pre function must be called on a valid node.
         * @post @a size() is an integer in the range [0, total_edges].
         */
        size_type degree() const {
            return grid->adjacencyList[vertexIdx].size();
        }
        
        /** Employs incident_iterator constructor to construct an iterator pointing to the
         * first neightboor of a vertex.
         *
         * edge_begin() constructs an indicent_iterator @a I, which points to the first
         * vertex adjacent to the vertex passed in the function call.
         * @return @a I of type incident_iterator.
         *
         * @pre function must be called on a valid node.
         */
        incident_iterator edge_begin() const {
            incident_iterator I(vertexIdx, 0, const_cast<Graph*>(grid));
            return I;
        }
        
        /** Employs incident_iterator constructor to construct an iterator pointing to the
         * last neightboor of a vertex.
         *
         * edge_end() constructs an indicent_iterator @a I, which points to the last
         * vertex adjacent to the vertex passed in the function call.
         * @return @a I of type incident_iterator.
         *
         * @pre function must be called on a valid node.
         */
        incident_iterator edge_end() const {
            incident_iterator I(vertexIdx, degree(), const_cast<Graph*>(grid));
            return I;
        }
        
        /** Checks whether vertexIdx and @a n are equivalent.
         *
         * == estabishes equality between two vertices.
         * @param[in] n vertex to compare.
         * @return True if @a n and vertexIdx have the same graph and index.
         * @return False if @a n and vertexIdx have different indecies or don't have the same
         * graph.
         */
        bool operator==(const Node& n) const {
            bool valReturn = (vertexIdx == n.vertexIdx && grid == n.grid) ? true : false;
            return valReturn;
        }
        
        /** Established if in graph order @n is greather than vertex.
         *
         * < compares global order of two vertices.
         * @param[in] n vertex to compare.
         * @return True if @a n and vertexIdx have the same graph and index of vertex if
         * less than index of @a n, or when graph of the vertex has smaller gloabl ordering 
         * than that of @a n's graph.
         * @return False above conditions are false.
         *
         * @pre ordering relation of vertices must obey trichotomy: given two vertices 
         * x and y, exactly one of the following holds true:
         * x < y, y < x, xor y == x.
         */
        bool operator<(const Node& n) const {
            //if they are in the same graph, order based on node id
            bool valReturn = ((n.grid == grid && vertexIdx < n.vertexIdx) ||
                              (grid <  n.grid)) ? true : false;
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
    
    /** Number of unique vertices in a graph.
     *
     * size() const gives a total number of vertices in a graph.
     * @return the number of nodes in the graph, which are of type unsigned.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return verticesVec.size();
    }
    
    /** Gives a number of unique vertices in a graph.
     *
     * num_nodes() const gives a total number of vertices in a graph.
     * @return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type num_nodes() const {
        return size();
    }
    
    /** Adds a node to the graph.
     *
     * add_node(const Point& position, const node_value_type& a = node_value_type())
     * appends a new node to a grid with location @a position and value @a a.
     * @param[in] position new node's location within a grid.
     * @param[in] a value associated with the node.
     * @return vertex with @a position and value @a a.
     *
     * @post num_nodes() is incremented by 1.
     * @post result_node.index() is assigned the artibutes of the old num_nodes().
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& a = node_value_type()) {
        vertexSt newVertex;
        newVertex.location = position; // Def location of a vertex.
        newVertex.val = a; // Def value of a vertex.
        verticesVec.push_back(newVertex); // Adding an extra node to the graph.
        // Pushing a placeholder vector to account for the new vector of vertices.
        adjacencyList.emplace_back(std::vector <size_type>());
        
        // Def a new vertex with current graph and vertex index.
        Node vertex;
        vertex.grid = const_cast<Graph*>(this);
        vertex.vertexIdx = verticesVec.size() - 1;
        return vertex;
    }
    
    /** Checks num of nodes and if graphs are the same.
     *
     * has_node(const Node& n) const establishes if graph contains a given vertex.
     * @parm[in] n is a vertex to check containment against.
     * @return True when @a n is in this Graph.
     * @return False when @a n is not in this Graph.
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        bool varReturn;
        varReturn = (const_cast< Graph* >(this) == n.grid &&
                     n.vertexIdx < verticesVec.size()) ? true : false;
        return varReturn;
    }
    
    /** Gives the node with the specified index.
     *
     * node(size_type i) const accesses the vertex at @a i.
     * parm[in] i unsigned index.
     * @return vertex stored at index @a i.
     *
     * @pre @a i must be in [0, num_nodes()).
     * @post @a i is assigned to result_node.index().
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
    class Edge : private totally_ordered<Edge> {
     public:
        /** Construct an invalid Edge. */
        Edge() {} // Invalid edge.
        
        /** Def head vertex, its index and current graph.
         *
         * node1() const gets the current head vertex.
         * @return current head vertex.
         *
         * @pre Graph is not empty.
         *
         * Complexity: O(1).
         */
        Node node1() const {
            Node head;
            head.vertexIdx = edgeHead;
            head.grid = const_cast< Graph* >(gridEdge);
            return head;
        }
        
        /** Def tail vertex, its index and current graph.
         *
         * node2() const gets the current head vertex.
         * @return current tail vertex.
         *
         * @pre Graph is not empty.
         *
         * Complexity: O(1).
         */
        Node node2() const {
            Node tail;
            tail.vertexIdx = edgeTail;
            tail.grid = const_cast< Graph* >(gridEdge);
            return tail;
        }
        
        /** Check if this edge and @a e are equal.
         *
         * == compares two edges.
         * @param[in] e edge to compare.
         * @return True if distance between @e's head and tail is the same as
         * the distane between the head and tail of this edge.
         * @return False otherwise.
         *
         * @pre Graph needs to contain another edge to compare @a e against.
         */
        bool operator==(const Edge& e) const {
            bool varReturn = (max(e.edgeHead, e.edgeTail) == max(edgeHead, edgeTail) &&
                              min(e.edgeHead, e.edgeTail) == min(edgeHead, edgeTail) &&
                              e.gridEdge == gridEdge) ? true : false;
            return varReturn;
        }
        
        /** Check if this edge is less than @a e in the global order.
         *
         * < compares two edges.
         * @param[in] e edge to compare.
         * @return True if @e is at a larger memory location than this edges, or when the
         * min, or max distance between head and tail of @a e is less than respective
         * distances of this edge.
         * @return False otherwise.
         *
         * @pre Graph needs to contain another edge to compare @a e against.
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
                                             max(edgeHead, edgeTail) < max(e.edgeHead, 						     e.edgeTail))) {
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
    
    /** Getting a total number of edges in a graph.
     *
     * num_edges() const gets full number of edges.
     * @return total_edges edge at @i.
     *
     * @pre Graph needs to contain edge at the given index @i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return total_edges;
    }
    
    /** Getting edge at the given index.
     *
     * edge(size_type i) const gets an edge at @i.
     * @param[in] i index of an edge to return.
     * @return Edge at @i.
     *
     * @pre Graph needs to contain edge at the given index @i.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        return *std::next(edge_begin(), i);
    }
    
    /** Checking if an edge exists between two nodes.
     *
     * has_edge(const Node& a, const Node& b) const tests if edge between @a and 
     * @b is in Graph.
     * @param[in] a starting node.
     * @param[in] b ending node.
     * @return True if @a a and @a b are connected by an edge.
     * @return False otherwise.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
     */
    bool has_edge(const Node& a, const Node& b) const {
        std::vector< size_type > vec = adjacencyList[a.vertexIdx];
        for( size_type j = 0; j < vec.size(); j++) {
            if( vec[j] == b.vertexIdx) return true; }
        return false;
    }
    
    /** Adds an edge between two nodes.
     *
     * add_edge(const Node& a, const Node& b) adds a new edge between @a a and 
     * @a b to the Graph.
     * @param[in] a starting node.
     * @param[in] b ending node.
     * @return an object newEdge starting @a a and ending @a b.
     *
     * @post checks if has_edge(@a a, @a b) was true and then returns the edge between 
     * the two vertices. If false a new edge is created and returned and number of
     * edges is incremented by 1..
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b) {
        // Creating an edge between vertices a and b.
        if(!has_edge(a, b)) {
            total_edges += 1;
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
    
    /** Remove all vertices and edges from this graph.
     *
     * clear() empties vertices and adjacencyList.
     * @return an object newEdge starting @a a and ending @a b.
     *
     * @post num_nodes() and num)edge() becomes 0.
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Removes all elements of the vertices and adjacencyList.
        verticesVec.clear();
        adjacencyList.clear();
    }
    
    //
    // Node Iterator
    //
    
    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        /** Construct an invalid NodeIterator. */
        NodeIterator() {}
        
        /** Obtains access to the vertex the NodeIterater points to.
         *
         * * gets the Node that NodeIterater points to at the moment.
         * @return Node @a node(niLocation).
         *
         * @pre iterator mustn't point to the last node in the Graph, i.e., iterator != graph.node_end() 
         */
        Node operator*() const {
            return niGrid -> node(niLocation);
        }
        
        /** Incrementing iterator by one.
         *
         * ++ increments the iterator by one position.
         * @return a reference to the incremented NodeIterator.
         *
         * @pre iterator mustn't be at the end of iterator, i.e., iterator != graph.node_end().
         */
        NodeIterator& operator++() {
            this->niLocation = niLocation + 1;
            return *this;
        }
        
        /** Checks if NodeIterator and the one provided are equal.
         *
         * == compares two node iterators.
         * @param[in] @a n node iterator for comparison.
         * @return True if node iterator graph are the same and if their locations were the same.
         * @return False otherwise.
         */
        bool operator==(const NodeIterator& n) const {
            bool valReturn = (n.niGrid == niGrid  && n.niLocation == niLocation) ? true : false;
            return valReturn;
        }
        
        /** Checks if NodeIterator and the one provided are equal.
         *
         * != compares two node iterators, provided and current using == operator.
         * @param[in] @a n node iterator for comparison.
         * @return True if node iterator graph are the same and if their locations were the same.
         * @return False otherwise.
         */
        bool operator!=(const NodeIterator& n) const {
            return !(*this == n);
        }
        
    private:
        
        friend class Graph;
        size_type niLocation;
        Graph* niGrid;
        
        NodeIterator(size_type p1, const Graph* g) : niLocation(p1), niGrid(const_cast<Graph *>(g)) {}
    };
    
    /** Gets the pointer to the first vertex in the vertex vector.
     *
     * node_begin() obtains the pointer to the start of the vertex vector.
     * @return @a NodeIterator, the node_iterator that points to the node id with zero.
     */
    node_iterator node_begin() const {
        return NodeIterator(0, this);
    }
    
    /** Gets the pointer to the last vertex in the vertex vector.
     *
     * node_end() obtains the pointer to the end of the vertex vector.
     * @return @a NodeIterator, the node_iterator that points to the the last stored node.
     */
    node_iterator node_end() const {
        return NodeIterator(size(), this);
    }
    
    //
    // Incident Iterator
    //
    
    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {}
        
        /** Accessess current Edge that IncidentIterater points to.
         *
         * * obtains the Edge that the IncidentIterater is pointing to.
         * @return Edge @a invEdge.
         *
         * @pre Iterator mustn't be at the end of edges vec., e.i., IncidentIterator is not equal to 
         * node.edge_end().
         */
        Edge operator*() const {
            Edge invEdge;
            invEdge.edgeHead = verIdxIter;
            invEdge.edgeTail = iterGrid -> adjacencyList[verIdxIter][innerIter];
            invEdge.gridEdge = iterGrid;
            
            return invEdge;
        }
        
        /** Incrementing iterator by one.
         *
         * ++ increments the iterator by one position.
         * @return a reference to the incremented IncidentIterator.
         *
         * @pre iterator mustn't be at the end of edges of the node, 
         * i.e., iterator != node.edge_end().
         */
        IncidentIterator& operator++() {
            this->innerIter = innerIter + 1;
            return *this;
        }
        
        /** Checks if IncidentIterator and the one provided are equal.
         *
         * == compares two incident iterators.
         * @param[in] @a i input incident iterator for comparison.
         * @return True if both incident iterators are equal 
         * @return False otherwise.
         */
        bool operator==(const IncidentIterator& i) const {
            if (i.iterGrid == iterGrid  && i.innerIter == innerIter && i.verIdxIter == verIdxIter) {
                return true;
            }
            return false;
        }
        
        /** Checks if IncidentIterator and the one provided are equal.
         *
         * != compares two incident iterators.
         * @param[in] @a i input incident iterator for comparison.
         * @return True if both incident iterators are equal using ==.
         * @return False otherwise.
         */
        bool operator!=(const IncidentIterator& i) const {
            return !(*this == i);
        }
        
    private:
        friend class Graph;
            size_type verIdxIter; // Vertex based iterator.
            size_type innerIter; // Iterator going across connectivity of a chosen vertex.
            Graph* iterGrid; // Graph pointer.
        
            IncidentIterator(size_type n1, size_type i1, const Graph* g) : verIdxIter(n1),
                             innerIter(i1), iterGrid(const_cast<Graph *>(g)) {}
    };
    
    //
    // Edge Iterator
    //
    
    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {}
        
        /** Accessess current Edge that EdgeIterater points to.
         *
         * * obtains the Edge that the EdgeIterater is pointing to.
         * @return Edge @a connection.
         *
         * @pre Iterator mustn't be at the edges of the graph, e.i., EdgeIterator is not equal to
         * graph.edge_end().
         */
        Edge operator*() const {
            Edge connection;
            connection.gridEdge = edgeIterGrid;
            connection.edgeHead = (*iterHead).index();
            connection.edgeTail = (*iterTail).node2().index();
            
            return connection;
        }
        
        /** Incrementing iterator by one.
         *
         * ++ increments the iterator by one position.
         * @return a reference to the incremented EdgeIterator.
         *
         * @pre iterator mustn't be at the end of edges of the graph,
         * i.e., iterator != graph.edge_end().
         */
        EdgeIterator& operator++() {
            ++iterTail;
            IncIterFun();
            return *this;
        }
        
        /** Checks if EdgeIterator and the one provided are equal.
         *
         * == compares two edge iterators.
         * @param[in] @a e input edge iterator for comparison.
         * @return True if both edge iterators are equal
         * @return False otherwise.
         */
        bool operator==(const EdgeIterator& e) const {
            return ( e.iterHead == iterHead && e.iterTail ==iterTail &&
                     e.edgeIterGrid == e.edgeIterGrid);
        }
        
        /** Checks if EdgeIterator and the one provided are equal.
         *
         * != compares two edge iterators using ==.
         * @param[in] @a e input edge iterator for comparison.
         * @return True if both edge iterators are equal
         * @return False otherwise.
         */
        bool operator!=(const EdgeIterator& e) const {
            return !( e.iterHead == iterHead && e.iterTail == iterTail &&
                      e.edgeIterGrid == e.edgeIterGrid);
        }
        
    private:
        friend class Graph;
        
        NodeIterator iterHead = (*edgeIterGrid).node_begin();
        IncidentIterator iterTail = (*iterHead).edge_begin();
        Graph* edgeIterGrid;
        
        /* Iterates over incident edges. */
        void IncIterFun() {
            // Executes if iterHead is not poiting to the last node.
            while (iterHead != (*edgeIterGrid).node_end()) {
                
                // Executes if iterTail is not pointing to the end of incidents.
                while (iterTail != (*iterHead).edge_end()) {
                    // Forces global ordering of edges to prevent their duplicates.
                    if ((*iterHead) < (*iterTail).node2()) ++iterTail;
                    else break;
                }
                
                // Updates location.
                if (iterTail == (*iterHead).edge_end()) {
                    ++iterHead;
                    iterTail = (*iterHead).edge_begin();
                }
                else break;
            }
            
            if (iterHead == (*edgeIterGrid).node_end() ) {
                iterHead = (*edgeIterGrid).node_begin();
                iterTail = (*iterHead).edge_begin();
            }
        }
        
        EdgeIterator(NodeIterator n1, IncidentIterator i1, const Graph* g) :
        iterHead(n1), iterTail(i1), edgeIterGrid(const_cast<Graph *>(g)) { IncIterFun();}
    };
    
    /** Obtains a pointer to the first edge in a graph.
     *
     * edge_begin() function to the pointer of the first edge of a graph.
     * @return an edge iterator.
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(node_begin(), (*node_begin()).edge_begin(), this);
    }
    
    /** Obtains a pointer to the last edge in a graph.
     *
     * edge_end() function to the pointer of the last edge of a graph.
     * @return an edge iterator.
     */
    edge_iterator edge_end() const {
        return EdgeIterator(node_end(), (*node_begin()).edge_end(), this);
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
        node_value_type val;
    };
    
};

#endif // CME212_GRAPH_HPP
