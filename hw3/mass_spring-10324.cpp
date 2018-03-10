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
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"
#include "Constraint.hpp"
#include "Force.hpp"


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
    double initialLength;       //< Node velocity
    double springConstant;     //< Node mass
    EdgeData() : initialLength(.01), springConstant(100) {}
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
 * @tparam G::node_value_type has a Point vel field.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template<typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        if (not(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)))
            n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    constraint(g);

    return t + dt;
}


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
#if 1
        // Diagonal edges: include as of HW2 #2
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }

    for (auto n: nodes) {
        n.value().mass = 1. / nodes.size();
    }

    for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
        (*it).value().initialLength = (*it).length();
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

        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            //std::cout << "t = " << t << std::endl;
            symp_euler_step(graph, t, dt,
                            make_combined_force(GravityForce{grav},
                                                MassSpringForce<GraphType>(),
                                                DampingForce(1. / nodes.size())),
                            make_combined_constraint(ConstraintApplication<PlaneConstraint>(PlaneConstraint(Point(0, 0, 1), .75)),
                                                     DeleteApplication<SpherePred>{
                                                             SpherePred(Point(.5, .5, -.5), .15)}));

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
