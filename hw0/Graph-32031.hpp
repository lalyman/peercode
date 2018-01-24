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

  // PUBLIC TYPE DEFINITIONS
    public:

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


        // CONSTRUCTORS AND DESTRUCTOR

        /** Construct an empty graph. */
        /* The constructor initialises two empty vectors which store the information about nodes and edges of the graph.
           It also initializes the count of number of nodes and edges in the graph to be zero. */
        Graph() :
            nodes_(),edges_(), num_nodes_(0),num_edges_(0){ } 

        /** Default destructor */
        ~Graph() = default;

        // NODES

        /** @class Graph::Node
        * @brief Class representing the graph's nodes.
        * Node objects are used to access information about the Graph's nodes.*/

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
                Node() { } //Invalid Node constructor when called by the user directly

                /** Return this node's position. */
                const Point& position() const {
                    return graph_ -> nodes_[index_]; 
                    //Access the position in the nodes vector of the graph containing the node at index_
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return index_; //This returns the index associated with Node object
                }

                /** Test whether this node and @a n are equal.
                *
                * Equal nodes have the same graph and the same index.
                */
                bool operator==(const Node& n) const {
                    if(n.index_ == index_  && n.graph_ == graph_) //Checking if the index and graph pointers of two nodes is the same
                        return true;
                    else
                        return false;
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
                    if(index_ < n.index() && graph_ == n.graph_)
                        return true;
                    else
                        return false;
                }
            // Use this space to declare private data members and methods for Node
            // that will not be visible to users, but may be useful within Graph.
            // i.e. Graph needs a way to construct valid Node objects
            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;

                size_type index_; // Works like a uid for a node
                const Graph* graph_; // Pointer to a Graph to which the node belongs

                // Private Node constructor which takes in index and pointer to a graph
                Node(size_type index, const graph_type* graph) : index_(index),graph_(graph) { }

        };

        /** Return the number of nodes in the graph.
        *
        * Complexity: O(1).
        */
        size_type size() const {
            return num_nodes_; //Private data member of  graph class
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
            nodes_.push_back(position); //nodes_ vector is a vector of Points/position  of nodes 
            num_nodes_++; //The number of nodes is incremented by one
            return Node(num_nodes_-1,this);// The new Node is returned by calling the Node constructor
        }

        /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
        */
        bool has_node(const Node& n) const {
            if(n.graph_ == this) 
            //Check if the graph attribute associated with the node is same as this graph object 
                return true;
            else
                return false;
        }

        /** Return the node with index @a i.
        * @pre 0 <= @a i < num_nodes()
        * @post result_node.index() == i
        *
        * Complexity: O(1).
        */
        Node node(size_type i) const {
            if(i >= size())
                return Node();// Return an Invalid node if index  is greater than the number of nodes in graph
            else
                return Node(i, this); //Return the node by calling the private constructor  of Node
        }

        // EDGES

        /** @class Graph::Edge
        * @brief Class representing the graph's edges.
        *
        * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
        * are considered equal if they connect the same nodes, in either order.
        */
        class Edge {
            public:
        
                /** Construct an invalid Edge. */
                Edge() {    } // An invalid edge constructor
    
                /** Return a node of this Edge */
                Node node1() const {
                    return Node(index1_,graph1_);  //Node constructor called by taking the index and graph pointer of node 
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    return Node(index2_,graph2_);      // call the Node constructor
                }

                /** Test whether this edge and @a e are equal.
                *
                * Equal edges represent the same undirected edge between two nodes.
                */
                bool operator==(const Edge& e) const {
                    // Check if the two nodes of two edges are same by calling the node1() and node2() methods of edge objects
                    if((e.node1() == node1() && e.node2() == node2()) || (e.node1() ==  node2() && e.node2() == node1()))
                        return true;
                    else
                        return false;
                }

                /** Test whether this edge is less than @a e in a global order.
                *
                * This ordering function is useful for STL containers such as
                * std::map<>. It need not have any interpretive meaning.
                */
                bool operator<(const Edge& e) const {
                    // Compare the edge_index of two edges which is a private attribute of edge object
                    if (index_edge_ < e.index_edge_)
                        return true;
                    else
                        return false;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                size_type index1_;// index associated with  node 1 of this edge object
                const Graph* graph1_; // graph pointer associated with node 1 of this edge object
                size_type index2_; //index associated with node 2 of this edge object
                const Graph* graph2_; // graph pointer associated with node 2 of this edge object
                size_type index_edge_; //index associated with this edge object

                // Private edge constructor. Takes in index and graph pointers of two nodes and an edge index 

                Edge(size_type  index1, const Graph* graph1, size_type index2, const Graph* graph2, size_type index_edge): index1_(index1),
                     graph1_(graph1), index2_(index2), graph2_(graph2), index_edge_(index_edge) { }
                // Use this space to declare private data members and methods for Edge
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Edge objects
        };

        /** Return the total number of edges in the graph.
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        size_type num_edges() const {
            return num_edges_;// Private data member of Graph class
        }

        /** Return the edge with index @a i.
        * @pre 0 <= @a i < num_edges()
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        Edge edge(size_type i) const {
            if (i >= num_edges_) // Since i is of unsigned data type, I dont check the condition that it has to be >0
            // If i is greater than the number of edges in the system  return an invalid edge
                return Edge();        // Invalid Edge
            else
                return edges_[i]; //edges_ is a vector containing the "Edges" objects of the graph 
                // I return the Edge object at index i stored in edges_vector of Graph object
        }

        /** Test whether two nodes are connected by an edge.
        * @pre @a a and @a b are valid nodes of this graph
        * @return True if for some @a i, edge(@a i) connects @a a and @a b.
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        bool has_edge(const Node& a, const Node& b) const {
            // call the has_node() method of Graph class to check if the nodes passed are valid nodes of the graph                
            if ( has_node (a)  &&  has_node (b) ){ // if the nodes passed are valid
                for ( size_type i = 0; i < num_edges() ; i++){ //iterate over edges_ vector of graph
                    if( (edges_[i].node1() == a && edges_[i].node2() == b) || // check if the nodes of an edge match the given nodes
                        (edges_[i].node1() == b && edges_[i].node2() == a)) // call the node1() and node2() methods
                        return true;
                }
                return false;// the nodes are present in graph but they don't have an edge
            }
            else
                return false;// the nodes are not valid nodes of the graph
        }

        /** Add an edge to the graph, or return the current edge if it already exists.
        * @pre @a a and @a b are distinct valid nodes of this graph
        * @return an Edge object e with  e.node1() == @a a and e.node2() == @a b
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
            if ( has_node(a) && has_node(b) && !(a==b)){// valid nodes of graph and distinct
                for(size_type i = 0; i< num_edges(); i++){//iterate over edges to find if edge exists_
                    if ((edges_[i].node1() == a && edges_[i].node2() == b ) || 
                        (edges_[i].node1() == b && edges_[i].node2() == a )) {
                        return edges_[i] ;//return the edge object
                    }
                }
                /* Couldn't use has_edge because it doesn't store the edge index and to 
                   return an edge object, I need to have edge index. So iterated over all elements again.*/
                /* if( has_edge(a, b) )
                    return Edge(a.index_,a.graph_,b.index_,b.graph_,*/
                num_edges_++; // if no edge exists then create a new edge and increment the counter
                // edges_ vector holds the edge objects which are created using private constructor
                edges_.push_back(Edge(a.index_,a.graph_,b.index_,b.graph_,num_edges_-1));
                return edges_[num_edges_-1]; //return the last added edge
            }
            else
                return Edge();        //If invalid node then return Invalid Edge
        }

        /** Remove all nodes and edges from this graph.
        * @post num_nodes()  == 0 && num_edges() == 0
        *
        * Invalidates all outstanding Node and Edge objects.
        */
        void clear() {
            edges_.clear(); //clears all the vectors and assigns 0 to all the counters
            nodes_.clear();
            num_nodes_=0;
            num_edges_=0;
        }

    // PRIVATE TYPE DEFINITIONS
    private:
        std::vector<Point> nodes_; //nodes_ vector to hold all the nodes information i.e. Point
        std::vector<Edge> edges_; //edges_ vector to hold all information about edges i.e. Edges
        size_type num_nodes_; // counter for number of nodes in graph
        size_type num_edges_; // counter for number of edges in graph

        // Use this space for your Graph class's internals:
        //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
