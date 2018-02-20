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
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template<typename V, typename E>

class Graph {    
private:
    struct edgeSt;                  //< Def of an edge structure.
    struct vertexSt;                //< Def of a vertex structure, i.e, its location.
    vector< vertexSt > verticesVec; //< Ordered vector of all vertices.
    vector< unsigned > idxVec;      //< Vector of indicies corresponding to each node.
    unsigned edgeNum = 0;           //< Total number of edges.
    
public:
    
    //
    // PUBLIC TYPE DEFINITIONS
    //
    
    /** Type of this graph. */
    using graph_type = Graph;
    
    /** Type of graph's input. */
    using node_value_type = V;
    
    /** Type of edge's input. */
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
     Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //
    
    /** Construct an empty graph. */
    Graph() : verticesVec(), idxVec() {}
    
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
        Node() {} //< Invalid node.
        
        
        /* Getting node's position.
         *
         * position() accessess vertex's position from the adjacency list.
         * @return @a location by reference to Point.
         *
         * @pre function must be called on a valid node.
         */
        Point &position() {
            return grid -> verticesVec[nodeIdx].location;
        }
        
        /* Getting node's position.
         *
         * position() const accessess vertex's position from the adjacency list.
         * @return @a location by reference to Point.
         *
         * @pre function must be called on a valid node.
         */
        const Point &position() const {
            return grid -> verticesVec[nodeIdx].location;
        }
        
        /* Getting idx of a graph's vertex.
         *
         * index() const gets the idx.
         * @return @a vertexIdx of type unsigned.
         *
         * @pre function must be called on a valid node.
         * @post @a vertexIdx is a number in the range [0, verticesVec[nodeIdx].size()).
         */
        size_type index() const {
            return grid -> verticesVec[nodeIdx].vertexIdx;
        }
        
        /* Adds new information to a vertex.
         *
         * value() allows to add additional info regarding the vertex.
         * @return @a val by reference to Node of type V.
         *
         * @pre function must be called on a valid node.
         */
        node_value_type &value() {
            return grid -> verticesVec[nodeIdx].vertexVal;
        }
        
        /** Add new information to a vertex.
         *
         * value() const allows to add additional info regardin a node the vertex.
         * @return const @a val by reference for node of type V. @a val is of type const.
         *
         * @pre graph pointer must point to a valid memeory location storing graph object.
         */
        const node_value_type &value() const {
            return grid -> verticesVec[nodeIdx].vertexVal;
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
            return grid -> verticesVec[nodeIdx].adjacencyList.size();
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
            incident_iterator I(grid -> verticesVec[nodeIdx].vertexIdx, 0, const_cast<Graph *>(grid));
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
            incident_iterator I(grid -> verticesVec[nodeIdx].vertexIdx, degree(),
                                const_cast<Graph *>(grid));
            return I;
        }
        
        /** Checks whether vertexIdx and @a n are equivalent.
         *
         * == estabishes equality between two vertices.
         * @param[in] n vertex to compare.
         * @return True if @a n and nodeIdx have the same graph and index.
         * @return False if @a n and nodeIdx have different indecies or don't have the same
         * graph.
         */
        bool operator==(const Node &n) const {
            return ( n.nodeIdx == nodeIdx && grid == n.grid ) ? true : false;
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
        bool operator<(const Node &n) const {
            return ((grid == n.grid && nodeIdx < n.nodeIdx) ||
                    (grid < n.grid))  ? true : false;
        }
        
    private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
        graph_type *grid;
        size_type nodeIdx;
        /* Constructor */
        Node(const graph_type *g, size_type i): grid(const_cast<graph_type *>(g)), nodeIdx(i) {}

    };
    
    /** Number of unique vertices in a graph.
     *
     * size() const gives a total number of vertices in a graph.
     * @return the number of nodes in the graph, which are of type unsigned.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return idxVec.size();
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
     * @return vertex with @a position, value @a a, vertex index @a idxVec.size(), and
     * adj. list @a adjacencyList equal to an empty vector of edge structs.
     *
     * @post num_nodes() is incremented by 1.
     * @post result_node.index() is assigned the artibutes of the old num_nodes().
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& a = node_value_type()) {
        vertexSt newVertex {
            //< Def different aspects of the vertex struct.
            .vertexVal = a,
            .location = position,
            .vertexIdx = static_cast<size_type>( idxVec.size() ),
            //< Pushing a placeholder vector to account for the new vector of edges.
            .adjacencyList = vector< edgeSt >()
        };
        verticesVec.emplace_back( newVertex );
        auto refVar = verticesVec.size() - 1;
        //< Populates index vector.
        idxVec.push_back( refVar );
        return Node( this, refVar );
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
    bool has_node(const Node &n) const {
        return (const_cast<graph_type *>(this) == n.grid) ? n.nodeIdx < verticesVec.size() : false;
    }
    
    /** Removes vertex from the graph.
     *
     * remove_node(const Node& n) ejects vertex @a n.
     * @param[in] n vertex to be removed from the graph.
     * @return 0 when failed and 1 when successful removal.
     *
     * @pre graph needs to have at least one valid vertex @a n.
     * @pre node @n must have a unique index in a graph.
     * @post graph has one less vertex which number is bounded below by 0.
     * @post indices vector has one less index which number is bounded below by 0.
     *
     * Complexity: O(deg(@a n)).
     */
    size_type remove_node(const Node &n) {
        if (! has_node(n) ) return 0;
        // < Removing vertex.
        auto& vertexStructure = verticesVec[n.nodeIdx];
        for (unsigned long i = 0; i < idxVec.size(); ++i) {
            auto& idxVer = verticesVec[idxVec[i]].vertexIdx;
            if( n.index() < idxVer ) idxVer  -= 1;
        }
        
        vertexStructure.vertexIdx = -1; // Invalidate vertex.
        auto& var = vertexStructure.adjacencyList;
        for( unsigned long i = 0; i < var.size(); ++i ) {
            // Removes vertex for all adjacent edges.
            auto &neigh = verticesVec[var[i].tail].adjacencyList;
            for (unsigned long j = 0; j < neigh.size(); ++j) {
                if( n.nodeIdx == neigh[j].tail ) {
                    neigh[j] = neigh[neigh.size() - 1];
                    neigh.pop_back();
                    --edgeNum;
                    break;
                }
            }
        }
        //< Moves pointer.
        for (auto i = idxVec.begin(); i != idxVec.end(); ++i) {
            if( n.nodeIdx == *i ) {
                idxVec.erase(i);
                break;
            }
        }
        return 1;
    }
    
    /** Gives the node with the specified index.
     *
     * node(size_type i) const accesses the vertex at @a i.
     * parm[in] i unsigned index.
     * @return vertex stored at index @a i.
     *
     * @pre @a i must be in [0, num_nodes()).
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        return Node(this, idxVec[i]);
    }
    
    /** Removes node using an iterator.
     *
     * remove_node(node_iterator n_it) ejects @a nodeIter from graph.
     * parm[in] nodeIter node iterator pointing to the vertex that will be removed.
     * @return iterator pointing to the next vertex.
     *
     * @pre initially @a nodeIter must point to a valid node.
     * @post @a nodeIter points to a new valid vertex.
     *
     * Complexity: O(deg(@a n)).
     */
    node_iterator remove_node(node_iterator nodeIter) {
        remove_node(*nodeIter);
        return nodeIter;
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
        
        /** Def head vertex, its index and current graph.
         *
         * node1() const gets the current head vertex.
         * @return current head vertex.
         *
         * @pre graph is not empty.
         *
         * Complexity: O(1).
         */
        Node node1() const {
            return gridEdge -> node(edgeHead);
        }
        
        /** Def tail vertex, its index and current graph.
         *
         * node2() const gets the current head vertex.
         * @return current tail vertex.
         *
         * @pre graph is not empty.
         *
         * Complexity: O(1).
         */
        Node node2() const {
            return gridEdge -> node(edgeTail);
        }
        
        /** Gets capacity value of an edge.
         *
         * const value() const obtains the value of an edge @a i.
         * @return const reference to the @a capacity of type E.
         *
         * @pre there must be a valid node in the graph.
         */
        edge_value_type &value() {
            for (unsigned long i = 0; i < gridEdge -> verticesVec[node1().nodeIdx].adjacencyList.size(); ++i) {
                if (node2().nodeIdx == gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[i].tail)
                    return gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[i].capacity;
            }
            return gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[0].capacity;
        }
        
        /** Gets capacity value of an edge.
         *
         * value() obtains the value of an edge @a i.
         * @return reference to the @a capacity of type E.
         *
         * @pre there must be a valid node in the graph.
         */
        const edge_value_type &value() const {
            for (auto i = 0; i < gridEdge-> verticesVec[node1().nodeIdx].adjacencyList.size(); ++i) {
                if (node2().nodeIdx == gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[i].tail)
                    return gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[i].capacity;
            }
            return gridEdge -> verticesVec[node1().nodeIdx].adjacencyList[0].capacity;
            
        }
        
        /** Computes length .
         *
         * length() const gets the current head vertex.
         * @return length, which is of type double.
         *
         * @pre graph has at least one edge.
         */
        double length() const {
            return norm( node1().position() - node2().position() );
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
            return ( min(edgeHead, edgeTail) == min(e.edgeHead, e.edgeTail) && max(edgeHead, edgeTail) == max(e.edgeHead, e.edgeTail) && e.gridEdge == gridEdge ) ? true : false;
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
            // Different graphs
            if( e.gridEdge > gridEdge ) return true;
            // Using an index of a lower degree of a head vertex when edges are in the same
            // graph and index of a lower degree of a tail vertex if degrees are not equal.
            else if( gridEdge == e.gridEdge ) {
                if( min(edgeHead, edgeTail) < min(e.edgeHead, e.edgeTail)) return true;
                else if (min(edgeHead, edgeTail) == min(e.edgeHead, e.edgeTail)) {
                    if (max(edgeHead, edgeTail) < max(e.edgeHead, e.edgeTail)) return true;
                }
            }
            return false;
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        
        graph_type *gridEdge; //< Pointer to the graph acessible from this class.
        size_type edgeHead;   //< Head and tail edges each with 8 bits of memory.
        size_type edgeTail;   //< Idx of an edge.
        
        // < Constructor.
        Edge(const graph_type* g, size_type head, size_type tail)
        : gridEdge(const_cast<graph_type *>(g)), edgeHead(g -> verticesVec[head].vertexIdx),
        edgeTail(g -> verticesVec[tail].vertexIdx) {}
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
        return edgeNum;
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
        for( unsigned long j = 0; j < verticesVec[a.nodeIdx].adjacencyList.size(); ++j ) {
            if( b.nodeIdx == verticesVec[a.nodeIdx].adjacencyList[j].tail ) return true;
        }
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
     * edges is incremented by 1.
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node &a, const Node &b) {
        //< Creating an edge between vertices a and b.
        if(!has_edge(a, b)) {
            edgeSt edgeRefHead {
                .tail = b.nodeIdx};
            edgeSt edgeRefTail {
                .tail = a.nodeIdx};
            
            //< Adding edge into the vector of stucts.
            verticesVec[a.nodeIdx].adjacencyList.emplace_back(edgeRefHead);
            verticesVec[b.nodeIdx].adjacencyList.emplace_back(edgeRefTail);
            edgeNum += 1;
            return Edge(this, a.nodeIdx, b.nodeIdx);
        }
        return Edge(this, a.nodeIdx, b.nodeIdx);
    }
    
    /** Removes edge from the graph.
     *
     * remove_edge(const Node& a, const Node& b) ejects an edge
     * connecting @a a and @a b.
     * @param[in] a head vertex.
     * @param[in] b tail vertex.
     * @return 0, a size_type integer.
     *
     * @pre @a a and @a b must be valid vartices in a graph.
     * @post graph has one less edge which number is bounded below by 0.
     * @post @a a and @a b both have one less degree.
     *
     * Complexity: O(max(deg(@a a), deg(@a b))).
     */
    size_type remove_edge(const Node &a, const Node &b) {
        if( !has_edge(a, b) ) return 0;
        
        //< Iterates over all edges of a to remove one connected to b.
        for (size_type i = 0; i < verticesVec[a.nodeIdx].adjacencyList.size(); ++i) {
            if (verticesVec[a.nodeIdx].adjacencyList[i].tail == b.nodeIdx) {
                verticesVec[a.nodeIdx].adjacencyList[i] =
                verticesVec[a.nodeIdx].adjacencyList[verticesVec[a.nodeIdx].adjacencyList.size() - 1];
                verticesVec[a.nodeIdx].adjacencyList.pop_back();
                break;
            }
        }
        edgeNum -= 1; //< Decrements number of edges.
        //< Iterates over all edges of b to remove one connected to a.
        for (size_type i = 0; i < verticesVec[b.nodeIdx].adjacencyList.size(); ++i) {
            if (verticesVec[b.nodeIdx].adjacencyList[i].tail == a.nodeIdx) {
                verticesVec[b.nodeIdx].adjacencyList[i] =
                verticesVec[b.nodeIdx].adjacencyList[verticesVec[b.nodeIdx].adjacencyList.size() - 1];
                verticesVec[b.nodeIdx].adjacencyList.pop_back();
                break;
            }
        }
        return 1;
    }
    
    /** Removes edge from the graph.
     *
     * remove_edge(const Edge& e) ejects an edge connecting @a e.
     * @param[in] e edge to be removed.
     * @return 0, a size_type integer.
     *
     * @pre @a e must be valid edge in a graph.
     * @post graph has one less edge which number is bounded below by 0.
     * @post @a a and @a b both have one less degree.
     *
     * Complexity: O(max(deg(@a a), deg(@a b))).
     */
    size_type remove_edge(const Edge& e) {
        return remove_edge(e.node1(), e.node2()); //< Utilize other remove edge function.
    }
    
    /** Removes edge from the graph using an edge iterator.
     *
     * remove_edge(edge_iterator e_it) ejects an edge.
     * @param[in] e_it edge iterator to be used to remove an edge.
     * @return next edge iterator from @a e_it.
     *
     * @pre @a e_it can't be at the last edges in the vector.
     * @post @a a and @a b, vertices to which edge is connected both have one less degree.
     *
     * Complexity: O(max(deg(@a a), deg(@a b))).
     */
    edge_iterator remove_edge(edge_iterator e_it) {
        remove_edge(*e_it);
        return e_it;
    }
    
    /** Empties the graph.
     *
     * clear() removes everything from the graph.
     *
     * @post num_edges() == 0 && verticesVec.size() == 0 && idxVec.size() == 0.
     */
    void clear() {
        verticesVec.clear();
        idxVec.clear();
        edgeNum = 0;
    }
    
    //
    // Node Iterator
    //
    
    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : private totally_ordered<NodeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     //< Element type
        using pointer           = Node *;                   //< Pointers to elements
        using reference         = Node &;                   //< Reference to elements
        using difference_type   = std::ptrdiff_t;           //< Signed difference
        using iterator_category = std::input_iterator_tag;  //< Weak Category, Proxy
        
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
        
        /** Increments the NodeIterater to the next position.
         *
         * ++ ups the NodeIterater by one.
         * @return reference to incremented NodeIterator.
         *
         * @pre iterator mustn't point to the last node in the Graph, i.e., iterator != graph.node_end()
         */
        NodeIterator &operator++() {
            this -> niLocation += 1;
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
            return (n.niGrid == niGrid  && n.niLocation == niLocation) ? true : false;
            
        }
        
    private:
    friend class Graph;
        size_type niLocation;
        graph_type *niGrid;
        NodeIterator(size_type p1, const graph_type* g) : niLocation(p1), niGrid(const_cast<Graph *>(g)) {}
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
            return Edge(iterGrid, iterGrid -> idxVec[verIdxIter], iterGrid -> verticesVec[ iterGrid -> idxVec[verIdxIter] ].adjacencyList[innerIter].tail);
        }
        
        /** Incrementing iterator by one.
         *
         * ++ increments the iterator by one position.
         * @return a reference to the incremented IncidentIterator.
         *
         * @pre iterator mustn't be at the end of edges of the node,
         * i.e., iterator != node.edge_end().
         */
        IncidentIterator &operator++() {
            this -> innerIter += 1;
            return *this;
            
        }
        
        /** Checks if IncidentIterator and the one provided are equal.
         *
         * == compares two incident iterators.
         * @param[in] @a i input incident iterator for comparison.
         * @return True if both incident iterators are equal
         * @return False otherwise.
         */
        bool operator==(const IncidentIterator &i) const {
            return (i.iterGrid == iterGrid  && i.innerIter == innerIter && i.verIdxIter == verIdxIter);
        }
        
    private:
    friend class Graph;
        
        size_type verIdxIter; //< Vertex based iterator.
        size_type innerIter;  //< Iterator going across connectivity of a chosen vertex.
        graph_type *iterGrid; //< Graph pointer.
        
        IncidentIterator(size_type n1, size_type i1, const graph_type* g) : verIdxIter(n1),
        innerIter(i1), iterGrid(const_cast<Graph *> (g)) {}
    };
    
    //
    // Edge Iterator
    //
    
    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private totally_ordered<EdgeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     //< Element type
        using pointer           = Edge *;                   //< Pointers to elements
        using reference         = Edge &;                   //< Reference to elements
        using difference_type   = std::ptrdiff_t;           //< Signed difference
        using iterator_category = std::input_iterator_tag;  //< Weak Category, Proxy
        
        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {}
        
        /** Accessess current Edge that EdgeIterater points to.
         *
         * * obtains the Edge that the EdgeIterater is pointing to.
         * @return Edge @a *iterTail.
         *
         * @pre Iterator mustn't be at the edges of the graph, e.i., EdgeIterator is not equal to
         * graph.edge_end().
         */
        Edge operator*() const {
            return *iterTail;
        }
        
        /** Incrementing iterator by one.
         *
         * ++ increments the iterator by one position.
         * @return a reference to the incremented EdgeIterator.
         *
         * @pre iterator mustn't be at the end of edges of the graph,
         * i.e., iterator != graph.edge_end().
         */
        EdgeIterator &operator++() {
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
                    e.edgeIterGrid == e.edgeIterGrid );
        }
        
    private:
    friend class Graph;
        //< Initialization of iterator variables.
        graph_type *edgeIterGrid;
        NodeIterator iterHead = (*edgeIterGrid).node_begin();
        IncidentIterator iterTail = (*iterHead).edge_begin();
        
        /* Iterates over incident edges. */
        void IncIterFun() {
            //< Executes if iterHead is not poiting to the last node.
            while (iterHead != (*edgeIterGrid).node_end()) {
                //< Executes if iterTail is not pointing to the end of incidents.
                while (iterTail != (*iterHead).edge_end()) {
                    //< Forces global ordering of edges to prevent their duplicates.
                    if ((*iterHead) < (*iterTail).node2()) ++iterTail;
                    else break;
                }
                
                //< Updates location.
                if (iterTail == (*iterHead).edge_end()) {
                    ++iterHead;
                    iterTail = (*iterHead).edge_begin();
                } else break;
            }
            
            if (iterHead == (*edgeIterGrid).node_end() ) {
                iterHead = (*edgeIterGrid).node_begin();
                iterTail = (*iterHead).edge_begin();
            }
        }
        
        EdgeIterator(const Graph* g, NodeIterator n1, IncidentIterator i1 ) :
        edgeIterGrid( const_cast<Graph *>(g) ), iterHead(n1), iterTail(i1) { IncIterFun(); }
        
    };
    
    
    /** Obtains a pointer to the first edge in a graph.
     *
     * edge_begin() function to the pointer of the first edge of a graph.
     * @return an edge iterator.
     */
    edge_iterator edge_begin() const {
        return EdgeIterator( this, node_begin(), (*node_begin()).edge_begin() );
    }
    
    /** Obtains a pointer to the last edge in a graph.
     *
     * edge_end() function to the pointer of the last edge of a graph.
     * @return an edge iterator.
     */
    edge_iterator edge_end() const {
        return EdgeIterator( this, node_end(), (*node_begin()).edge_end() );
    }
    
private:
    //< Def of edge's attributes.
    struct edgeSt {
        size_type tail;
        edge_value_type capacity;
    };

    //< Def of a vertex's structure.
    struct vertexSt {
        Point location;
        vector< edgeSt > adjacencyList;
        size_type vertexIdx;
        node_value_type vertexVal;
    };
    
    
};

#endif // CME212_GRAPH_HPP
