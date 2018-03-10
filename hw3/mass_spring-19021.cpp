/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>
#include <math.h>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
    double L;
    double K;

    EdgeData() : L(0), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // fixes the two points that hold the graph
        if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // apply constraints to the sheet
    constraint(g, t);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
    /** Return the force applying to @a n at time @a t.
         *
         * For HW2 #1, this is a combination of mass-spring force and gravity,
         * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
         * model that by returning a zero-valued force. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        // HW2 #1: YOUR CODE HERE
        (void) t; // silence compiler warnings
        Point sprg_f;
        Point grav_f = n.value().mass * Point(0, 0, -grav);
        if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
            return(Point(0, 0, 0));
        }
        for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Edge e = *it;
            Point x_i = n.position();
            Point x_j = e.node2().position();
            sprg_f += -e.value().K * (x_i - x_j) / (norm(x_i - x_j)) * (norm(x_i - x_j) - e.value().L);
        }
        return(sprg_f + grav_f);
    }
};

// =======================================================================================================================
// ====================================================== Forces =========================================================
// =======================================================================================================================

/** Force function object for HW2 #1. */
struct GravityForce {

    /** Applies the gravity force to a given node.
         * @tparam[in] n   Node object
         * @param[in]  t   Time at current step
         * @return The gravity force applied to the given node
         *
         * @pre t > 0
         * @post For node n,
         *			n.value().mass * Point(0,0,-grav)
         *
         * Complexity: O(1)
         */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return(n.value().mass * Point(0,0,-grav));
    }
};

/** Defines the mass spring force */
struct MassSpringForce {

    /** Applies the mass spring force to a given node.
         * @tparam[in] n   Node object
         * @param[in]  t   Time at current step
         * @return The mass spring force applied to the given node
         *
         * @pre t > 0
         * @post For edge e, let x_i = e.node1().position() and x_j = e.node2().position(), then
         * 			force = sum(-e.value().K * (x_i - x_j) / (norm(x_i - x_j)) * (norm(x_i - x_j) - e.value().L) for all e
         *
         * Complexity: O(n.degree())
         */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        Point sprg_f;
        for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Edge e = *it;
            Point x_i = n.position();
            Point x_j = e.node2().position();
            sprg_f += -e.value().K * (x_i - x_j) / (norm(x_i - x_j)) * (norm(x_i - x_j) - e.value().L);
        }
        return(sprg_f);
    }
};

/** Defines the damping force */
struct DampingForce {

    /** Applies the damping force to a given node.
         * @tparam[in] n   Node object
         * @param[in]  t   Time at current step
         * @return The damping force applied to the given node
         *
         * @pre t > 0
         * @post For node n,
         * 			n.value().vel = - c * n.value().vel
         *
         * Complexity: O(1)
         */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return(-(n.value().vel) * c);
    }
};

/** Defines a null force, used for defualt value. */
struct null_f {

    /** Returns a zero force.
         * @tparam[in] n    Node object
         * @param[in]  t    Time at current step
         * @return A zero force to a node.
         *
         * Complexity: O(1)
         */
    template <typename NODE>
    Point operator() (NODE n, double t) {
        (void) t; (void) n;
        return(Point(0, 0, 0));
    }
};

/** Combines up to three different forces. */
template <typename f_1, typename f_2, typename f_3>
struct make_combined_force {
    f_1 one_; f_2 two_; f_3 three_;

    // init the combined force struct
    make_combined_force(f_1 one = null_f(), f_2 two = null_f(), f_3 three = null_f()) : one_(one), two_(two), three_(three) {}

    /** Applies forces to one node.
         * @param[in] g		Node object
         * @param[in] t     Time at current step
         *
         * @pre t > 0
         * @post If node n has specified forces applied to it.
         *
         * Complexity: O(1)
         */
    template <typename NODE>
    Point operator() (NODE n, double t) {
        (void) t;
        return(one_(n, t) + two_(n, t) + three_(n, t));
    }
};

// =========================================================================================================================
// ===================================================== Constraints =======================================================
// =========================================================================================================================

/** Defines the plane constraint for the graph. */
struct plane {

    /** Finds all nodes that violate the condition.
         * @tparam[in] g   Graph object
         * @param[in]  t   Time at current step
         *
         * @pre t > 0
         * @post For all nodes n in g s.t n.position.z < -0.75,
         *			n.position.z = -0.75 and n.value.vel.z = 0
         *
         * Complexity: O(g.num_nodes())
         */
    template <typename G>
    void operator()(G &g, double t) const {
        (void) t;
        for(auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if(n.position().z < -0.75) {
                n.position().z = -0.75;
                n.value().vel.z = 0;
            }
        }
    }
};

/** Defines the sphere constraint for the graph. */
struct sphere {

    /** Finds all nodes that violate the condition.
         * @tparam[in] g	Graph object
         * @param[in]  t    Time at current step
         *
         * @pre t > 0
         * @post Given sphere center point c, and radius r. For all nodes n in g s.t norm(n.position() - c) < r,
         *			n.position.z = c + r/norm(n.position() - c) * (n.position() - c)\
         *			n.value().vel = n.value().vel - (inner_prod(n.value().vel, n.position() - c) / norm(n.position() - c))
         * 					* n.position() - c) / norm(n.position() - c);
         *
         * Complexity: O(g.num_nodes())
         */
    template <typename G>
    void operator()(G &g, double t) const {
        (void) t;
        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;
        for(auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if(norm(n.position() - c) < r) {
                n.position() = c + r/norm(n.position() - c) * (n.position() - c);
                Point r_i = (n.position() - c) / norm(n.position() - c);
                n.value().vel = n.value().vel - (inner_prod(n.value().vel, r_i)) * r_i;
            }
        }
    }
};

/** Defines the sphere constraint for the graph to remove nodes. */
struct sphere_remove {

    /** Finds all nodes that violate the condition.
         * @tparam[in] g	Graph object
         * @param[in] t     Time at current step
         *
         * @pre t > 0
         * @post Given sphere center point c, and radius r. For all nodes n in g s.t norm(n.position() - c) < r,
         *			g.remove_node(n)
         *
         * Complexity: O(g.num_nodes())
         */
    template <typename G>
    void operator()(G &g, double t) const {
        (void) t;
        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;
        for(auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if(norm(n.position() - c) < r) {
                g.remove_node(n);
            }
        }
    }
};

/** Defines a null constraint, used for default value. */
struct null_c {

    /** Does nothing.
         * @tparam[in] g    Graph object
         * @param[in]  t    Time at current step
         */
    template <typename G>
    void operator()(G &g, double t) {
        (void) g; (void) t;
    }
};

/** Combines up to three different constraints. */
template <typename c_1, typename c_2>
struct make_combined_constraint {
    c_1 one_; c_2 two_;

    // init the combined constraint struct
    make_combined_constraint(c_1 one = null_c(), c_2 two = null_c()) : one_(one), two_(two) {}

    /** Applies constraints to one node.
         * @tparam[in] g    Graph object
         * @param[in]  t    Time at current step
         *
         * @pre t > 0
         * @post All nodes in graph g have appropriate constraint applied.
         *
         * Complexity: O(1)
         */
    template <typename G>
    void operator()(G &g, double t) {
        one_(g, t); two_(g, t);
    }
};

int main(int argc, char** argv)
{
    // Check arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
        exit(1);
    }

    // Construct an empty graph
    GraphType graph;

    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    Point p;
    std::vector<typename GraphType::node_type> nodes;
    while (CME212::getline_parsed(nodes_file, p))
        nodes.push_back(graph.add_node(p));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
        graph.add_edge(nodes[t[0]], nodes[t[1]]);
        graph.add_edge(nodes[t[0]], nodes[t[2]]);

        // Diagonal edges: include as of HW2 #2
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);

        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }

    // HW2 #1 YOUR CODE HERE
    // Set initial conditions for your nodes, if necessary.

    for(size_t i = 0; i < nodes.size(); i++) {
        nodes[i].value().mass = 1.0 / graph.num_nodes();
        nodes[i].value().vel = Point(0, 0, 0);
    }

    for(auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
        Edge e = *it;
        e.value().L = e.length();
    }

    c = 1.0 / graph.num_nodes();

    // Print out the stats
    std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

    // Launch the Viewer
    CME212::SFML_Viewer viewer;
    auto node_map = viewer.empty_node_map(graph);

    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    viewer.center_view();

    // We want viewer interaction and the simulation at the same time
    // Viewer is thread-safe, so launch the simulation in a child thread
    bool interrupt_sim_thread = false;
    auto sim_thread = std::thread([&](){

        // Begin the mass-spring simulation
        double dt = 0.001;
        double t_start = 0;
        double t_end = 5.0;

        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            //std::cout << "t = " << t << std::endl;
            symp_euler_step(graph, t, dt,
                            make_combined_force<GravityForce, MassSpringForce, DampingForce>(GravityForce(), MassSpringForce(), DampingForce()),
                            make_combined_constraint<plane, sphere_remove>(plane(), sphere_remove()));

            viewer.clear();
            node_map.clear();

            // Update viewer with nodes' new positions
            viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
            viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
            viewer.set_label(t);

            // These lines slow down the animation for small graphs, like grid0_*.
            // Feel free to remove them or tweak the constants.
            if (graph.size() < 100)
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

    });  // simulation thread

    viewer.event_loop();

    // If we return from the event loop, we've killed the window.
    interrupt_sim_thread = true;
    sim_thread.join();

    return 0;
}
