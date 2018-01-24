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
class Graph 
{
    private:
        // HW0: YOUR CODE HERE
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
        Graph(): nodes(), edges() 
        {
            // HW0: YOUR CODE HERE
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
        class Node 
        {
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
                Node() 
                {
                    // HW0: YOUR CODE HERE
                }

                /** Return this node's position. */
                const Point& position() const 
                {
                    // HW0: YOUR CODE HERE
                    return graph_->nodes[node_id_];
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const 
                {
                    // HW0: YOUR CODE HERE
                    return node_id_;
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const 
                {
                    // HW0: YOUR CODE HERE
                    if (node_id_ == n.node_id_)
                    {
                        return true;
                    }
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
                bool operator<(const Node& n) const 
                {
                    // HW0: YOUR CODE HERE
                    if (node_id_ < n.node_id_)
                    {
                        return true;
                    }
                    return false;
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                // HW0: YOUR CODE HERE
                // Use this space to declare private data members and methods for Node
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Node objects
                Graph* graph_;
                size_type node_id_;
                Node(const Graph* graph, size_type node_id)
                    : graph_(const_cast<Graph*>(graph)), node_id_(node_id) 
                {
                    }
        };
        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         */
        size_type size() const 
        {
            // HW0: YOUR CODE HERE
            return nodes.size();
        }

        /** Synonym for size(). */
        size_type num_nodes() const 
        {
            return size();
        }

        /** Add a node to the graph, returning the added node.
         * @param[in] position The new node's position
         * @post new num_nodes() == old num_nodes() + 1
         * @post result_node.index() == old num_nodes()
         * Complexity: O(1) amortized operations.
         */
        Node add_node(const Point& position) 
        {
            // HW0: YOUR CODE HERE
            nodes.push_back(position);
            return Node(this, size()-1);        // Invalid node
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const 
        {
            // HW0: YOUR CODE HERE
            if (n.node_id_ < this->nodes.size())
            {
                return true;
            }
            return false;
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const 
        {
            // HW0: YOUR CODE HERE
            return Node(this, i);        // Invalid node
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
        class Edge 
        {
            public:
                /** Construct an invalid Edge. */
                Edge() 
                {
                    // HW0: YOUR CODE HERE
                }

                /** Return a node of this Edge */
                Node node1() const 
                {
                    // HW0: YOUR CODE HERE
                    return Node(graph_,node1_id_);      
                }

                /** Return the other node of this Edge */
                Node node2() const 
                {
                    // HW0: YOUR CODE HERE
                    return Node(graph_, node2_id_);      
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const 
                {
                    if ((node1_id_ == e.node1_id_ and node2_id_ == e.node2_id_) or (node1_id_ == e.node2_id_ and node2_id_ == e.node1_id_))
                    {
                        return true;
                    }
                    return false;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const 
                {
                    if (node_id_ < e.node_id_)
                    {
                        return true;
                    }
                    return false;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // HW0: YOUR CODE HERE
                // Use this space to declare private data members and methods for Edge
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Edge objects
                Graph* graph_;
                size_type node_id_;
                size_type node1_id_;
                size_type node2_id_;

                Edge(const Graph* graph , size_type node_id , size_type node1_id , size_type node2_id)
                    : graph_(const_cast<Graph*>(graph)) , node_id_(node_id) , node1_id_(node1_id) , node2_id_(node2_id) 
                {
                }
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const 
        {
            // HW0: YOUR CODE HERE
            
            // Number of edges is the size of the edges vector
            return edges.size();
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const 
        {
            // HW0: YOUR CODE HERE
            return Edge(this , i , edges[i].first , edges[i].second); 
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const 
        {
            // HW0: YOUR CODE HERE
            
            // Iterate through all the edges.
            // Scenario 1: a == edges.first and b == edges.second
            // Scenario 2: a == edges.second and b == edges.first
            for (unsigned i = 0; i < edges.size(); i++)
            {
                if (this->edges[i].first == a.node_id_ and this->edges[i].second == b.node_id_) 
                {
                    return true;
                }
                else if (this->edges[i].first == b.node_id_ and this->edges[i].second == a.node_id_)
                {
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
            // HW0: YOUR CODE HERE
            
            // Iterate through all the edges. Check to see if there is same edge.
            // Scenario 1: a == edge.first and b == edge.second
            // Scenario 2: a == edge.second and b == edge.first
            for (unsigned i = 0; i < edges.size(); i++)
            {
                if (this->edges[i].first == a.node_id_ and this->edges[i].second == b.node_id_) 
                {
                    return Edge(this, i ,a.node_id_ , b.node_id_);
                }
                else if (this->edges[i].first == b.node_id_ and this->edges[i].second == a.node_id_)
                {
                    return Edge(this, i , a.node_id_ , b.node_id_);
                }
            }
            
            // Push back new edge to the back of the vector
            // if same edge is not found
            edges.push_back(std::make_pair(a.node_id_ , b.node_id_));

            // Return the last edge
            return Edge(this, edges.size()-1 , a.node_id_ , b.node_id_);        
        }

        /* Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() 
        {
            nodes.clear();
            edges.clear();
        }

    private:

        // HW0: YOUR CODE HERE
        // Use this space for your Graph class's internals:
        //   helper functions, data members, and so forth.
        
        /* A vector of points.
         * x, y and z coordinates are stored at each point.
         * The index number on the vector is same as the node number.
         */
        std::vector<Point> nodes;

        /* A vector of edges.
         * In each edge, there is a pair of node stating which nodes are connected.
         */
        std::vector<std::pair<size_type, size_type>> edges;
};

#endif // CME212_GRAPH_HPP
