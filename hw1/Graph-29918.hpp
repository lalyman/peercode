#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/// @file Graph.hpp
/// @brief An undirected graph type
/// 

#include <algorithm>
#include <cassert>
#include <set>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/// @class Graph
/// @brief A template for 3D undirected graphs.
/// 
/// Users can add and retrieve nodes and edges. Edges are unique (there
/// is at most one edge between any pair of distinct nodes).
/// 


template <typename Value>
class Graph {
  
    // Pre-declare structs.
    struct InternalNode;
    struct InternalEdge;

  public:

    /// Type of this graph.
    using graph_type = Graph;

    /// Predeclaration of Node type.
    class Node;
    /// Synonym for Node (following STL conventions).
    using node_type = Node;
    using node_value_type = Value;

    /// Predeclaration of Edge type.
    class Edge;
    /// Synonym for Edge (following STL conventions).
    using edge_type = Edge;

    // PLEASE READ:
    // 
    // After speaking with Slaven and Noam and Andreas about the proxy
    // design pattern, I have done a large refactor of my code from
    // HW0. My proxy objects no longer contain pointers to the graph
    // container. Instead, they contain pointers to their corresponding
    // 'real' objects: InternalNode and InternalEdge.
    //
    // This implementation allows my underlying data structures to be
    // vectors of the proxy objects, which makes my code is fully
    // analagous with Slaven's examples from Lecture 3. In particular,
    // my member variables are now of type std::vector<Node> and
    // std::vector<Edge>. Additionally, the incident_edges_ member of
    // the InternalEdge struct is also a container of proxy objects: an
    // std::set<Edge>.
    //
    // This implementation allowed me to come up with a very sly
    // implementation of NodeIterator, EdgeIterator, and
    // IncidentIterator: I simply type alias'ed them to
    // std::vector<Node>::iterator, std::vector<Edge>::iterator, and
    // std::set<Edge>::iterator respectively.
    //
    // I personally think this is a cleaner implementation. However, in
    // order to assuage any concerns that I haven't actaully
    // implemented the iterator classes, I have provided two additional
    // files: viewer2.cpp and Graph2.hpp. These contain the full
    // implementation of the Iterator classes (but they're really just
    // wrappers on the standard library classes I mentioned above---so
    // I hope that for future assignments you'll let me keep the
    // implementation I have in this file). In order to check that
    // those work as well, simply do a 'make viewer2' instead of 'make
    // viewer'.
    //
    // Best,
    // Brett

    /// Type of node iterators, which iterate over all graph nodes.
    /// Synonym for NodeIterator
    using NodeIterator = typename std::vector<Node>::iterator;
    using node_iterator = NodeIterator;

    /// Type of edge iterators, which iterate over all graph edges.
    /// Synonym for EdgeIterator
    using EdgeIterator = typename std::vector<Edge>::iterator;
    using edge_iterator = EdgeIterator;

    /// Type of incident iterators, which iterate incident edges to a node.
    /// Synonym for IncidentIterator
    using IncidentIterator = typename std::set<Edge>::iterator;
    using incident_iterator = IncidentIterator;

    /// Type of indexes and sizes.
    /// Return type of Graph::Node::index(), Graph::num_nodes(),
    /// Graph::num_edges(), and argument type of Graph::node(size_type)
    using size_type = unsigned;


    // CONSTRUCTORS AND DESTRUCTOR

    /// Construct an empty graph.
    Graph() : nodes_(new std::vector<Node>()), edges_(new std::vector<Edge>()) {
    }

    /// Destroy the graph and free memory.
    ~Graph() {

      // Since edges_ and nodes_ store pointers to structs on the heap,
      // we need to free that memory.
      nodes_->clear();
      edges_->clear();
      delete nodes_;
      delete edges_;
    }


    /// @class Graph::Node
    /// @brief Class representing the graph's nodes.
    /// 
    /// Node objects are used to access info about the Graph's nodes.
    /// 
    class Node : private totally_ordered<Node> {
      
      public:
        
        /// Construct an invalid node.
        /// 
        /// Valid nodes are obtained from the Graph class, but it is
        /// occasionally useful to declare an @i invalid node, and
        /// assign a valid node to it later. For example:
        //
        /// @code
        /// Graph::node_type x;
        /// if (...should pick the first node...)
        ///   x = graph.node(0);
        /// else
        ///   x = some other node using a complicated calculation
        /// do_something(x);
        /// @endcode
        /// 
        Node() {
        }

        /// Return this node's position.
        const Point &position() const {
          return node_->pos_;
        }

        /// Return this node's value.
        const Value &value() const {
          return node_->val_;
        }

        /// Return this node's value for setting.
        Value &value() {
          return node_->val_;
        }

        /// Return this node's index, a number in the range [0, graph_size).
        size_type index() const {
          return index_;
        }

        /// Return the degree of the node by getting the number of
        /// edges incident to the node.
        size_type degree() const {
          return node_->incident_edges_->size();
        }

        /// Return iterator pointing to first incident edge.
        IncidentIterator edge_begin() const {
          // Since an InternalNode already contains a std::set<Edge> of
          // incident Edges, we can just return an iterator to the
          // beginning of that set.
          return node_->incident_edges_->begin();
        }

        /// Return iterator to one past the last incident edge.
        IncidentIterator edge_end() const {
          return node_->incident_edges_->end();
        }

        /// Test whether this node and @a n are equal.
        /// 
        /// Equal nodes have the same graph and the same index.
        /// 
        bool operator==(const Node &n) const {
          // We can just compare the addresses of the 'real' objects.
          return node_ == n.node_;
        }

        /// Test whether this node is less than @a n in a global order.
        /// 
        /// This ordering function is useful for STL containers such as
        /// std::map<>. It need not have any geometric meaning.
        /// 
        /// The node ordering relation must obey trichotomy: For any
        /// two nodes x and y, exactly one of x == y, x < y, and y < x
        /// is true.
        ///
        bool operator<(const Node &n) const {
          return node_ < n.node_;
        }

      private:
        
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Pointer back to the real object as per proxy design pattern.
        InternalNode *node_;

        // Node's unique index in the set of nodes.
        size_type index_;

        /// Private constructor.
        Node(const Point &position, Value &value, size_type ind) 
            : node_(new InternalNode(position, value)), index_(ind) {
        }
    };

    /// Return the number of nodes in the graph.
    /// 
    /// Complexity: O(1).
    /// 
    size_type size() const {
      return nodes_->size();
    }

    /// Synonym for size().
    size_type num_nodes() const {
      return size();
    }

    /// Add a node to the graph, returning the added node.
    /// @param[in] position The new node's position
    /// @post new num_nodes() == old num_nodes() + 1
    /// @post result_node.index() == old num_nodes()
    /// 
    /// Complexity: O(1) amortized operations.
    /// 
    Node &add_node(const Point &position, const Value &val = Value()) {
      // Just add a new internal_node to the end of nodes_.
      nodes_->push_back(Node(position, const_cast<Value &>(val), num_nodes()));
      return (*nodes_)[num_nodes() - 1];
    }

    /// Determine if a Node belongs to this Graph
    /// @return True if @a n is currently a Node of this Graph
    /// 
    /// Complexity: O(1).
    /// 
    bool has_node(const Node& n) const {
      // Checking equality ensures that their Graphs are the same.
      return n == node(n.index());
    }

    /// Return the node with index @a i.
    /// @pre 0 <= @a i < num_nodes()
    /// @post result_node.index() == i
    /// 
    /// Complexity: O(1).
    /// 
    const Node &node(size_type i) const {
      // Simply create a new instance of the proxy.
      return (*nodes_)[i];
    }


    /// @class Graph::Edge
    /// @brief Class representing the graph's edges.
    /// 
    /// Edges are order-insensitive pairs of nodes. Two Edges with the
    /// same nodes are considered equal if they connect the same nodes,
    /// in either order.
    /// 
    class Edge : private totally_ordered<Edge> {
      
      public:
        
        /// Construct an invalid Edge.
        Edge() {
        }

        // There's no need to store the node indices in the proxy
        // object since we can easily get them from the internal_edge
        // object defined below.

        /// Return a node of this Edge.
        Node node1() const {
          return edge_->one_;
        }

        /// Return the other node of this Edge.
        Node node2() const {
          return edge_->two_;
        }

        /// Test whether this edge and @a e are equal.
        /// 
        /// Equal edges are the same undirected edge between two nodes.
        /// 
        bool operator==(const Edge &e) const {
          // Can simply compare the addresses of the 'real' objects.
          return edge_ = e.edge_;
        }

        /// Test whether this edge is less than @a e in a global order.
        /// 
        /// This ordering function is useful for STL container such as
        /// std::map<>. It need not have any interpretive meaning.
        ///
        bool operator<(const Edge &e) const {
          return edge_ < e.edge_;
        }

      private:
      
        // Allow Graph to access Edge's private members.
        friend class Graph;

        // Pointer to the 'real' object as per proxy design pattern.
        InternalEdge *edge_;

        // Edge's unique index in the set of edges.
        size_type index_;

        /// Private constructor.
        Edge(const Node &a, const Node &b, size_type ind)
            : edge_(new InternalEdge(a, b)), index_(ind) {
        }
    };

    /// Return the total number of edges in the graph.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    size_type num_edges() const {
      return edges_->size();
    }

    /// Return the edge with index @a i.
    /// @pre 0 <= @a i < num_edges()
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    const Edge edge(size_type i) const {
      return (*edges_)[i];
    }

    /// Test whether two nodes are connected by an edge.
    /// @pre @a a and @a b are valid nodes of this graph
    /// @return True if for some @a i, edge(@a i) connects @a a and @a b.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    bool has_edge(const Node &a, const Node &b) const {
      return find_edge(a, b) != a.edge_end();
    }

    /// Add an edge to the graph, or return the current edge if it
    /// already exists.
    /// @pre @a a and @a b are distinct valid nodes of this graph
    /// @return an Edge object e with e.node1() == @a a and 
    ///                               e.node2() == @a b
    /// @post has_edge(@a a, @a b) == true
    /// @post If old has_edge(@a a, @a b),
    ///         new num_edges() == old num_edges().
    ///       Else
    ///         new num_edges() == old num_edges() + 1.
    /// 
    /// Can invalidate edge indexes -- in other words, old edge(@a i)
    /// might not equal new edge(@a i). Must not invalidate outstanding
    /// Edge objects.
    /// 
    /// Complexity: No more than O(num_nodes() + num_edges()),
    /// hopefully less
    /// 
    const Edge &add_edge(const Node &a, const Node &b) {

      IncidentIterator it = find_edge(a, b);
      if (it != a.edge_end()) {
        return *it;
      }

      edges_->push_back(Edge(a, b, num_edges()));
      (*nodes_)[a.index()].node_->incident_edges_->insert((*edges_)[num_edges() - 1]);
      (*nodes_)[b.index()].node_->incident_edges_->insert((*edges_)[num_edges() - 1]);
      return (*edges_)[num_edges() - 1];
    }

    /// Remove all nodes and edges from this graph.
    /// @post num_nodes() == 0 && num_edges() == 0
    /// 
    /// Invalidates all outstanding Node and Edge objects.
    /// 
    void clear() {
      nodes_->clear();
      edges_->clear();
    }

    /// Returns iterator pointing to the first node in the graph.
    NodeIterator node_begin() const {
      return nodes_->begin();
    }

    /// Returns iterator pointing to one past the last node.
    NodeIterator node_end() const {
      return nodes_->end();
    }

    /// Returns iterator pointing to first edge in the graph.
    EdgeIterator edge_begin() const {
      return edges_->begin();
    }

    /// Returns iterator pointing to one past the last edge.
    EdgeIterator edge_end() const {
      return edges_->end();
    }

  private:

    // Helper function to get an incident edge if it exists.
    // 
    // Searches through the edges incident to Node a for one whose
    // other endpoint is b. If such a node is found, return iterator to
    // that edge, otherwise return iterator to end.
    //
    IncidentIterator find_edge(const Node &a, const Node &b) const {

      // Search incident edges to find the edge whose other endpoint
      // is b using the std::find_if function with a custom predicate.
      // References for find_if with custom predicates:
      //   http://www.cplusplus.com/reference/algorithm/find_if/
      //   https://stackoverflow.com/q/6679096/902812
      IncidentIterator it = std::find_if(
        a.edge_begin(), a.edge_end(), 

        // Use a lambda function for the custom predicate so that we
        // can access the node b.
        [b](const Edge &candidate_edge) {
          return (
            // Either endpoint can be b, so check both. Since we can
            // assume a and b are from this graph, we just need to
            // check indices.
            candidate_edge.node1() == b or candidate_edge.node2() == b); } );
      
      return it;
    }

    // This graph's internals are based on the description in the last
    // bullet of this section from Wikipedia on adjacency lists:
    // https://en.wikipedia.org/wiki/Adjacency_list#Implementation_details

    // Internal struct to bundle node data (Point) with a set of
    // indices of edges to which it's connected.
    struct InternalNode {

      // Constructor makes pointer to set that's stored on the heap.
      InternalNode(const Point &pos, Value &val)
          : pos_(pos), val_(val), incident_edges_(new std::set<Edge>()) {
      }

      // Free the heap memory.
      ~InternalNode() {
        delete incident_edges_;
      }

      // Store node data.
      const Point pos_;
      Value val_;

      // Store a pointer to the set instead of the entire set in order
      // to keep the amount of memory on the stack small.
      std::set<Edge> *incident_edges_;
    };

    // Internal struct to store an edge.
    struct InternalEdge {
      
      InternalEdge(const Node &a, const Node &b)
          : one_(a), two_(b) {
      }
      
      // Keep pointers to node endpoints (as indices).
      const Node one_;
      const Node two_;
    };
    
    // Positions associated with each node (access by index).
    std::vector<Node> *nodes_;
    std::vector<Edge> *edges_;
};

#endif // CME212_GRAPH_HPP
