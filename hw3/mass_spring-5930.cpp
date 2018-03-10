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
#include <iostream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;
    double L;
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using NodeType=Node;
using Edge = typename GraphType::edge_type;
using size_type = typename GraphType::size_type;
using IncidentIterator =typename GraphType::incident_iterator;
using NodeIter = typename GraphType::node_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */
NodeIter nearest_node(const GraphType &g, const Point &point) {
    // Initialize iterator of graph nodes
    NodeIter start = g.node_begin();
    NodeIter end = g.node_end();

    // Use std library to search for minimum element
    return std::min_element(start, end,
                            [point](const NodeType &a, const NodeType &b) {
                                return norm(a.position() - point) <
                                       norm(b.position() - point);
                            });
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     cst    Function object defining the constraint per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template<typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C cst) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // Apply constraint
    cst(g, t);

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
    template<typename NODE>
    Point operator()(NODE n, double t) {
        // HW2 #1: YOUR CODE HERE
        // Fix points
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
            return Point(0, 0, 0);
        }

        // Gravitational force
        Point f_grav = n.value().mass * Point(0, 0, -grav);
        // Initialize spring force
        Point f_spri = Point(0, 0, 0);

        // Sum over all edges
        auto start = n.edge_begin();
        auto end = n.edge_end();
        for (; start != end; ++start) {
            Edge e = *start;

            Node na;
            // Get adjacent node
            if (e.node1() == n) {
                na = e.node2();
            } else {
                na = e.node1();
            }

            Point diff = n.position() - na.position();
            f_spri -= e.value().K * diff / norm(diff) *
                    (norm(diff) - e.value().L);
        }
        (void) t;
        return f_grav + f_spri;
    }

    double L;
    double K;
    Problem1Force(double l, double k) : L(l), K(k) {};
};

/** Gravitational force function object. */
struct GravityForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point f_grav = n.value().mass * Point(0, 0, -grav);
        (void) t;
        return f_grav;
    }
};

/** Spring force function object. */
struct MassSpringForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point f_spri = Point(0, 0, 0);

        // Sum over all edges
        auto start = n.edge_begin();
        auto end = n.edge_end();
        for (; start != end; ++start) {
            Edge e = *start;

            Node na;
            // Get adjacent node
            if (e.node1() == n) {
                na = e.node2();
            } else {
                na = e.node1();
            }

            Point diff = n.position() - na.position();
            // Spring force
            f_spri -= e.value().K * diff / norm(diff) *
                    (norm(diff) - e.value().L);
        }
        (void) t;
        return f_spri;
    }
};

/** Damping force function object. */
struct DampingForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point f_damp = -c * n.value().vel;
        (void) t;
        return f_damp;
    }
    double c;
};

/** Variational combined force function object. */
template<typename...>
struct make_combined_force;

/** Two combined force function object. */
template<typename F1, typename F2>
struct make_combined_force<F1, F2> {
    F1 f1;
    F2 f2;

    template<typename NODE>
    Point operator()(NODE n, double t) {
        return f1(n, t) + f2(n, t);
    }

    make_combined_force(F1 force1, F2 force2) {
        f1 = force1;
        f2 = force2;
    }
};

/** Three combined force function object. */
template<typename F1, typename F2, typename F3>
struct make_combined_force<F1, F2, F3> {
    F1 f1;
    F2 f2;
    F3 f3;

    template<typename NODE>
    Point operator()(NODE n, double t) {
        return f1(n, t) + f2(n, t) + f3(n, t);
    }

    make_combined_force(F1 force1, F2 force2, F3 force3) {
        f1 = force1;
        f2 = force2;
        f3 = force3;
    }
};

/** Constant constraint function object. */
struct Const_cst {
    void operator()(GraphType &g, double t) {
        for (unsigned i = 0; i < idx.size(); i++) {
            g.set_pos(idx[i], pt[i]);

        }
        (void) t;
    }
    // All fixed point and corresponding value
    std::vector <size_type> idx;
    std::vector <Point> pt;
};

/** Z plane constraint function object. */
struct Zpln_cst {
    void operator()(GraphType &g, double t) {
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            auto n = *it;
            if (n.position().z < val) {
                n.position().z = val;
                n.value().vel.z = 0;
            }
        }
        (void) t;
    }
    double val;
};

/** Sphere constraint function object. */
struct Sphere_cst {
    void operator()(GraphType &g, double t) {
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            auto n = *it;
            double dist = norm(n.position() - c);
            if (dist < r) {
                g.remove_node(n);
                --it;
            }
        }
        (void) t;
    }
    Point c;
    double r;
};

/** Make combined constraint function object. */
template<typename C1, typename C2>
struct make_combined_cst {
    C1 c1;
    C2 c2;

    void operator()(GraphType &g, double t) {
        c1(g, t);
        c2(g, t);
        (void) t;
    }

    make_combined_cst(C1 cst1, C2 cst2) {
        c1 = cst1;
        c2 = cst2;
    }
};

int main(int argc, char **argv) {
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
    std::array<int, 4> t;
    while (CME212::getline_parsed(tets_file, t)) {
        graph.add_edge(nodes[t[0]], nodes[t[1]]);
        graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
        // Diagonal edges: include as of HW2 #2
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }


    // HW2 #1 YOUR CODE HERE
    // Set initial conditions for your nodes, if necessary.

    // Initialize mass and velocity
    double m = 1. / ((double) graph.num_nodes());
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        auto n = *it;
        n.value().vel = Point(0, 0, 0);
        n.value().mass = m;
    }

    // Initialize stiffness and rest length
    for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
        auto e = *it;
        e.value().K = 100;
        e.value().L = e.length();
    }

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
    auto sim_thread = std::thread([&]() {

        // Begin the mass-spring simulation
        double dt = 0.001;
        double t_start = 0;
        double t_end = 5.0;

        // Get force functions and combine
        GravityForce gf;
        MassSpringForce sf;
        DampingForce df;
        df.c = 1. / ((double) graph.num_nodes());

        using combine_force_type =
        make_combined_force<GravityForce, MassSpringForce, DampingForce>;
        combine_force_type force = combine_force_type(gf, sf, df);

        // Get constraint functions and combine
        Const_cst c0;
        c0.idx.push_back((*nearest_node(graph, Point(0, 0, 0))).uid());
        c0.idx.push_back((*nearest_node(graph, Point(1, 0, 0))).uid());
        c0.pt.push_back(Point(0, 0, 0));
        c0.pt.push_back(Point(1, 0, 0));

        Sphere_cst c1;
        c1.c = Point(0.5, 0.5, -0.5);
        c1.r = 0.15;

        using combine_cst_type = make_combined_cst<Const_cst, Sphere_cst>;
        combine_cst_type cst = combine_cst_type(c0, c1);

        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            //std::cout << "t = " << t << std::endl;
            symp_euler_step(graph, t, dt, force, cst);

            // Clear the viewer’s nodes and edges
            viewer.clear();
            node_map.clear();

            // Update viewer with nodes’ new positions and new edges
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
