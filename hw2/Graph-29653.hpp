#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
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

using node_value_type = V;
using edge_value_type = E;

//
// CONSTRUCTORS AND DESTRUCTOR
//

/** Construct an empty graph. */
Graph() : nodes(new std::unordered_map<size_type, std::pair<Point, node_value_type>>()),
        edges(new std::unordered_map<
                      size_type,
                      std::pair<std::pair<size_type, size_type>, edge_value_type>>()),
        adjacency(new std::unordered_map<size_type, std::unordered_map<size_type, size_type> >()){
}

/** Destructor */
~Graph() {
        delete nodes;
        delete edges;
        delete adjacency;
}

//
// NODES
//

/** @class Graph::Node
 * @brief Class representing the graph's nodes.
 *
 * Node objects are used to access information about the Graph's nodes.
 */
class Node : totally_ordered<Node>{
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
Node(){
}

/** Return this node's position. */
const Point& position() const {
        return graph_->nodes->at(index()).first;
}

Point& position() {
        return graph_->nodes->at(index()).first;
}

/** Return this node's index, a number in the range [0, graph_size). */
size_type index() const {
        return nid_;
}

// Supply definitions AND SPECIFICATIONS for:
// node_value_type& value();
// const node_value_type& value() const;
// size_type degree() const;
// incident_iterator edge_begin() const;
// incident_iterator edge_end() const;

/** Test whether this node and @a n are equal.
 *
 * Equal nodes have the same graph and the same index.
 */
bool operator==(const Node& n) const {
        return index() == n.index() && graph_ == n.graph_;
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
        return index() < n.index();
}

/**
 * Returns value of the node
 * @return value of the node
 */
node_value_type& value(){
        return graph_->nodes->at(index()).second;
}

/**
 * Returns value of the node
 * @return value of the node
 */
const node_value_type& value() const {
        return graph_->nodes->at(index()).second;
}

/**
 * Returns degree of the node (number of incident edges)
 * @return degree of the node
 */
size_type degree() const {
        if(graph_->adjacency->find(index()) == graph_->adjacency->end())
                return 0;
        return graph_->adjacency->at(index()).size();
}

/**
 * Iterator through the edges incident to the given node
 * @return incident_iterator object
 */
incident_iterator edge_begin() const {
        return incident_iterator(
                       graph_->adjacency->at(index()).begin(),
                       nid_,
                       graph_);
}

/**
 * Iterator to the end of the edges incident to the given node
 * @return incident_iterator object
 */
incident_iterator edge_end() const {
        return incident_iterator(
                       graph_->adjacency->at(index()).end(),
                       nid_,
                       graph_);
}

private:
size_type nid_;
graph_type* graph_;

Node(const size_type id, const graph_type* graph) :
        nid_(id), graph_(const_cast<graph_type*>(graph)) {
}

// Allow Graph to access Node's private member data and functions.
friend class Graph;
};

/** Return the number of nodes in the graph.
 *
 * Complexity: O(1).
 */
size_type size() const {
        return nodes->size();
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

node_type add_node(const Point& position, const node_value_type& value = node_value_type()){
        std::pair<size_type, std::pair<Point, node_value_type>> newpair(size(), std::pair<Point, node_value_type>(position, value));
        nodes->insert(newpair);
        return node_type(size()-1, this);
}

/** Determine if a Node belongs to this Graph
 * @return True if @a n is currently a Node of this Graph
 *
 * Complexity: O(1).
 */
bool has_node(const node_type& n) const {
        return nodes->find(n.index()) != nodes->end();
}

/** Return the node with index @a i.
 * @pre 0 <= @a i < num_nodes()
 * @post result_node.index() == i
 *
 * Complexity: O(1).
 */
node_type node(size_type i) const {
        // HW0: DONE
        // if(i >= size()) // DEBUG
        //   std::cout << "i = " << i << ", size = " << size() << '\n';
        assert(i < size());
        return node_type(i, this);
}

/**
 * Sets the value associated with given node with specified value
 * @param i index of the node to modify
 * @param v value to give to this node
 */
void set_value(const size_type i, const node_value_type& v) {
        auto it = nodes->find(i);
        if(it != nodes->end())
                it->second.second = v;
}

/**
 * Removes node and all edges incident to it
 *
 * @param  n node to remove
 * @return   index of the node just after the one removed
 *
 * @post     (old) size() == (new) size()+1
 * @post     (old) node((old) size() - 1) is now indexed at (old) n.index()
 *
 * In order to remove the given node, we have replaced it with the last node of
 * the graph, and rewired the edges and values to make sure nothing changes from
 * its point of view. Then we remove the last element.
 */
size_type remove_node(const Node& n) {
        size_type s(size());
        // need to handle:
        //  - suppressing the incident edges
        //  - suppressing the node
        //      * remove its value
        //      * remove its position
        //  - changing the indices in adjacency matrix and in edges

        // removing incident edges
        while(n.degree() > 0) {
                Edge incident = *(n.edge_begin());
                remove_edge(incident);
        }
        // at this point the adjacency map at index 'index' is empty

        size_type index(n.index());
        size_type replacement(size()-1);
        auto replacement_node = node(replacement);

        // adjacency and edges
        while (replacement_node.degree() > 0) {
                // we iterate over the neighbors of replacement_node
                auto it = replacement_node.edge_begin();
                auto neighbor = (*it).node2();
                auto edge_index = (*it).index();

                // we change the index of replacement_node in the neighbor's map
                std::pair<size_type, size_type> newpair_for_neigh(index, edge_index);
                adjacency->at(neighbor.index()).insert(newpair_for_neigh);
                adjacency->at(neighbor.index()).erase(replacement);

                // we insert the edge to neighbor at index's position in
                // adjacency
                std::pair<size_type, size_type> newpair_for_rep(neighbor.index(), edge_index);
                if (adjacency->find(index) == adjacency->end()){
                    // in case index didn't have any neighbor, we need to
                    // initialize the neighboring map
                    std::unordered_map<size_type, size_type> newemptymap;
                    std::pair<size_type, std::unordered_map<size_type, size_type>> newpair(index, newemptymap);
                    adjacency->insert(newpair);
                }
                adjacency->at(index).insert(newpair_for_rep);

                // we replace 'replacement' by 'index' in the edge data
                if(edges->at(edge_index).first.first == replacement)
                        edges->at(edge_index).first.first = index;
                else
                        edges->at(edge_index).first.second = index;

                // we get rid of the neighbor's entry in replacement's map
                adjacency->at(replacement).erase(neighbor.index());
        }

        // if necessary, we get rid of replacement's entry in adjacency
        if (adjacency->find(replacement) != adjacency->end()) {
                adjacency->at(replacement).clear();
                adjacency->erase(replacement);
        }

        // migrate and erase node data
        nodes->erase(index);
        if(index != replacement) {
                std::pair<size_type, std::pair<Point, node_value_type>> newposition(index, nodes->at(replacement));
                nodes->insert(newposition);
                nodes->erase(replacement);
        }

        // sanity check, necessary for good functioning, but not sufficient
        // of course
        assert(s == size()+1);
        return index+1;
}

/**
 * Removes node and all edges incident to it
 * @param  n_it Iterator pointing to the node to remove
 * @return      Valid iterator to the node that replaced the one removed, ie
 *              at the same index.
 *
 * @post        return value is a valid iterator
 * @post        same as remove_node(const Node&)
 */
node_iterator remove_node(node_iterator n_it) {
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
class Edge : totally_ordered<Edge>{
public:
/** Construct an invalid Edge. */
Edge(){
}

/** Return the index of this edge */
size_type index() const {
        return eid_;
}

/** Return a node of this Edge */
node_type node1() const {
        return node_type(nid_a_, graph_);
}

/** Return the other node of this Edge */
node_type node2() const {
        return node_type(nid_b_, graph_);
}

/**
 * Returns the euclidean length of the edge_end
 */
double length() const {
        return norm(node2().position() - node1().position());
}

/**
 * Returns value of the edge
 * @return value of the edge
 */
edge_value_type& value(){
        return graph_->edges->at(index()).second;
}

/**
 * Returns value of the edge
 * @return value of the edge
 */
const edge_value_type& value() const {
        return graph_->edges->at(index()).second;
}

/** Test whether this edge and @a e are equal.
 *
 * Equal edges represent the same undirected edge between two nodes.
 */
bool operator==(const Edge& e) const {
        return nid_a_ == e.nid_a_ && nid_b_ == e.nid_b_ && graph_ == e.graph_;
}

/** Test whether this edge is less than @a e in a global order.
 *
 * This ordering function is useful for STL containers such as
 * std::map<>. It need not have any interpretive meaning.
 */
bool operator<(const Edge& e) const {
        if (graph_ == e.graph_)
                return index() < e.index();
        return graph_ < e.graph_;
}

private:

size_type eid_, nid_a_, nid_b_;
graph_type* graph_;

Edge(const size_type id, const size_type id1, const size_type id2, const graph_type* graph) :
        eid_(id), nid_a_(id1), nid_b_(id2), graph_(const_cast<graph_type*>(graph)){
}

// Allow Graph to access Edge's private member data and functions.
friend class Graph;
/* HW0: DONE
 * Use this space to declare private data members and methods for Edge
 * that will not be visible to users, but may be useful within Graph.
 * i.e. Graph needs a way to construct valid Edge objects
 */
};

/** Return the total number of edges in the graph.
 *
 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
 */
size_type num_edges() const {
        return edges->size();
}

/** Return the edge with index @a i.
 * @pre 0 <= @a i < num_edges()
 *
 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
 */
Edge edge(size_type i) const {
        assert(i < num_edges());
        std::pair<size_type, size_type> pair_i(edges->at(i).first);
        return Edge(i, pair_i.first, pair_i.second, this);
}

/** Test whether two nodes are connected by an edge.
 * @pre @a a and @a b are valid nodes of this graph
 * @return True if for some @a i, edge(@a i) connects @a a and @a b.
 *
 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
 */
bool has_edge(const node_type& a, const node_type& b) const {
        return adjacency->find(a.index()) != adjacency->end() && adjacency->at(a.index()).find(b.index()) != adjacency->at(a.index()).end();
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
Edge add_edge(const node_type& a, const node_type& b, const edge_value_type& newval = edge_value_type()) {
        size_type index;
        if (has_edge(a, b)) {
                return Edge(adjacency->at(a.index()).at(b.index()), a.index(), b.index(), this);
        }
        else {
                index = num_edges();

                std::pair<std::pair<size_type, size_type>, edge_value_type> newvaluepair(std::pair<size_type, size_type>(a.index(), b.index()), newval);
                std::pair<size_type, std::pair<std::pair<size_type, size_type>, edge_value_type>> newpair(index, newvaluepair);

                edges->insert(newpair);

                std::pair<size_type, size_type> newpaira(b.index(), index);
                std::pair<size_type, size_type> newpairb(a.index(), index);

                if(adjacency->find(a.index()) == adjacency->end()) {
                        std::unordered_map<size_type, size_type> newmapa;
                        std::pair<size_type, std::unordered_map<size_type, size_type> > newmapitema(a.index(), newmapa);
                        adjacency->insert(newmapitema);
                }
                if(adjacency->find(b.index()) == adjacency->end()) {
                        std::unordered_map<size_type, size_type> newmapb;
                        std::pair<size_type, std::unordered_map<size_type, size_type> > newmapitemb(b.index(), newmapb);
                        adjacency->insert(newmapitemb);
                }

                adjacency->at(a.index()).insert(newpaira);
                adjacency->at(b.index()).insert(newpairb);


                return Edge(index, a.index(), b.index(), this);
        }
}

/**
 * Remove edge and edge information from the data
 * @param  a first node defining the edge
 * @param  b second node defining the edge
 * @return   index of the index after the one we removed, or 0 if nothing was
 *           removed
 *
 * @post     (old) num_edges() == (new) num_edges() +1
 * @post     (old) edge((old) num_edges()-1) is now indexed at the index of
 *           the one we removed, as is its data.
 *
 * To remove the edge, we just need to transfer the data from the last edge to
 * the spot of the one to remove, and then erase the last edge from the map.
 */
size_type remove_edge(const Node& a, const Node& b) {
        size_type n(num_edges());
        // need to handle edges, adjacency:
        //  in adjacency,   just two values to remove from maps
        //  in edges,       remove one element and swap the last one with the gap

        // check if anything to remove at all
        if (!has_edge(a, b))
                return 0;
        size_type index(adjacency->at(a.index()).at(b.index()));
        size_type replacement(num_edges()-1);

        // we get rid of the edge's information in adjacency
        adjacency->at(a.index()).erase(b.index());
        adjacency->at(b.index()).erase(a.index());

        // we get rid of the edge's data
        edges->erase(index);
        if(index != replacement) {
                // we replace the erased data with that of the replacement edge
                std::pair<size_type, std::pair<std::pair<size_type, size_type>, edge_value_type> > newedge(index, edges->at(replacement));
                edges->insert(newedge);
                edges->erase(replacement);
        }

        // we make sure the replacement edge's links refer to it with its new
        // index, ie `index` and no longer `replacement`
        if(index != replacement) {
                auto it = adjacency->at(edge(index).node1().index()).find(edge(index).node2().index());
                auto it2 = adjacency->at(edge(index).node2().index()).find(edge(index).node1().index());
                it->second=index;
                it2->second=index;
        }

        // sanity check
        assert(n == num_edges()+1);
        return index+1;
}

/**
 * Remove edge and edge information from the data
 * @param  e edge to remove
 * @return   index after that of the removed edge, or 0 if nothing happenned
 *
 * @post     same as remove_edge(const Node&, const Node&)
 */
size_type remove_edge(const Edge& e) {
        return remove_edge(e.node1(), e.node2());
}

/**
 * Remove edge and edge information from the data
 * @param  e_it iterator to the edge to remove
 * @return      valid iterator to the edge that now replaces the one we removed
 *
 * @post        same as remove_edge(const Node&, const Node&)
 */
edge_iterator remove_edge(edge_iterator e_it) {
        size_type removed(remove_edge(*e_it));
        return edge_begin() + removed;
}

/** Remove all nodes and edges from this graph.
 * @post num_nodes() == 0 && num_edges() == 0
 *
 * Invalidates all outstanding Node and Edge objects.
 */
void clear() {
        nodes->clear();
        edges->clear();
        adjacency->clear();
}

//
// Node Iterator
//

/** @class Graph::NodeIterator
 * @brief Iterator class for nodes. A forward iterator. */
class NodeIterator : totally_ordered<NodeIterator> {
public:
// These type definitions let us use STL's iterator_traits.
using value_type        = Node;                         // Element type
using pointer           = Node*;                        // Pointers to elements
using reference         = Node&;                        // Reference to elements
using difference_type   = std::ptrdiff_t;               // Signed difference
using iterator_category = std::input_iterator_tag;      // Weak Category, Proxy

/** Construct an invalid NodeIterator. */
NodeIterator() {
}

NodeIterator(const size_type x, const graph_type* graph) :
        p(x), graph_(const_cast<graph_type*>(graph)) {
}

// HW1 #2: YOUR CODE HERE
// Supply definitions AND SPECIFICATIONS for:
// Node operator*() const
// NodeIterator& operator++()
// bool operator==(const NodeIterator&) const

/**
 * dereferencing operator
 * @param p index of the node currently pointed to by the iterator
 */
value_type operator*() const {
        return graph_->node(p);
}

/**
 * pre-incrementing operator
 * @return reference to current iterator, that has been incremented
 */
NodeIterator& operator++(){
        ++p;
        return *this;
}

/**
 * equality operator
 * @return whether the iterator is equal to the one passed in argument
 */
bool operator==(const NodeIterator& otherIt) const {
        return p == otherIt.p && graph_ == otherIt.graph_;
}

private:
friend class Graph;
size_type p;
graph_type* graph_;
};

// Supply definitions AND SPECIFICATIONS for:
// node_iterator node_begin() const
// node_iterator node_end() const

/**
 * Returns the iterator that begins the nodes of the graph
 * @return iterator to first node of the graph
 */
node_iterator node_begin() const {
        return node_iterator(0, this);
}

/**
 * Returns the iterator that ends the nodes of the graph
 * @return iterator to last node of the graph
 */
node_iterator node_end() const {
        return node_iterator(size(), this);
}

//
// Incident Iterator
//

/** @class Graph::IncidentIterator
 * @brief Iterator class for edges incident to a node. A forward iterator. */
class IncidentIterator : totally_ordered<IncidentIterator> {
public:
// These type definitions let us use STL's iterator_traits.
using value_type        = Edge;                         // Element type
using pointer           = Edge*;                        // Pointers to elements
using reference         = Edge&;                        // Reference to elements
using difference_type   = std::ptrdiff_t;               // Signed difference
using iterator_category = std::input_iterator_tag;      // Weak Category, Proxy

/** Construct an invalid IncidentIterator. */
IncidentIterator() {
}

IncidentIterator(std::unordered_map<size_type, size_type>::iterator i,
                 const size_type node_index,
                 const graph_type* graph) :
        it(i), source_node_index(node_index),
        graph_(const_cast<graph_type*>(graph)) {
}

// HW1 #3: YOUR CODE HERE
// Supply definitions AND SPECIFICATIONS for:
// Edge operator*() const
// IncidentIterator& operator++()
// bool operator==(const IncidentIterator&) const

/**
 * dereferencing operator
 * @param  it iterator over the edges incident to the source_node
 * @return    incident edge currently pointed to by the iterator
 */
value_type operator*() const {
        auto pair(*it);  //pair represents (second_node, index_of_the_edge)
        return Edge(pair.second,
                    source_node_index,
                    pair.first, graph_);
}

/**
 * pre-incrementing operator
 * @return reference to current iterator, that has been incremented
 */
IncidentIterator& operator++(){
        ++it;
        return *this;
}

/**
 * equality iterator
 * @return whether the current iterator is equal to the one given in argument
 */
bool operator==(const IncidentIterator& otherIncid) const {
        return it == otherIncid.it
               && source_node_index == otherIncid.source_node_index
               && graph_ == otherIncid.graph_;
}


private:
friend class Graph;
std::unordered_map<size_type, size_type>::iterator it;
size_type source_node_index;
graph_type* graph_;
};

//
// Edge Iterator
//

/** @class Graph::EdgeIterator
 * @brief Iterator class for edges. A forward iterator. */
class EdgeIterator : totally_ordered<EdgeIterator>{
public:
// These type definitions let us use STL's iterator_traits.
using value_type        = Edge;                         // Element type
using pointer           = Edge*;                        // Pointers to elements
using reference         = Edge&;                        // Reference to elements
using difference_type   = std::ptrdiff_t;               // Signed difference
using iterator_category = std::input_iterator_tag;      // Weak Category, Proxy

/** Construct an invalid EdgeIterator. */
EdgeIterator() {
}

EdgeIterator(const size_type x, const graph_type* graph) :
        p(x), graph_(const_cast<graph_type*>(graph)) {
}

// HW1 #5: YOUR CODE HERE
// Supply definitions AND SPECIFICATIONS for:
// Edge operator*() const
// EdgeIterator& operator++()
// bool operator==(const EdgeIterator&) const

/**
 * dereferencing operator
 * @param p index of the edge currently pointed to by the iterator
 */
value_type operator*() const {
        return graph_->edge(p);
}

/**
 * pre-incrementing operator
 * @return reference to current iterator, that has been incremented
 */
EdgeIterator& operator++(){
        ++p;
        return *this;
}

/**
 * equality iterator
 * @return whether the current iterator is equal to the one given in argument
 */
bool operator==(const EdgeIterator& otherIt) const {
        return p == otherIt.p && graph_ == otherIt.graph_;
}

private:
friend class Graph;
size_type p;
graph_type* graph_;
};

// HW1 #5: YOUR CODE HERE
// Supply definitions AND SPECIFICATIONS for:
// edge_iterator edge_begin() const
// edge_iterator edge_end() const

/**
 * Returns the iterator that begins the edges of the graph
 * @return iterator to first edge of the graph
 */
edge_iterator edge_begin() const {
        return edge_iterator(0, this);
}

/**
 * Returns the iterator that ends the edges of the graph
 * @return iterator to last edge of the graph
 */
edge_iterator edge_end() const {
        return edge_iterator(num_edges(), this);
}

private:

// maps index to node positions and values
std::unordered_map<size_type, std::pair<Point, node_value_type>>* nodes;

// maps index to edge information and value
std::unordered_map<size_type, std::pair<std::pair<size_type, size_type>, edge_value_type>>* edges;

// maps a node index to a map containing pairs of (index of adjacent node, index of the given edge)
// we trade off storage and data redundancy for constant-time access
// which is eventually probably more important
std::unordered_map<size_type, std::unordered_map<size_type, size_type> >* adjacency;

};

#endif // CME212_GRAPH_HPP
