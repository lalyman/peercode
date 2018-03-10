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

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
double N;

/** Custom structure of data to store with Nodes */
struct NodeData {
    double fixed;
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    NodeData() : fixed(0),vel(0), mass(1) {}
};


struct EdgeData {
    double K; // spring constant
    double L; // length of edge
    EdgeData() : K(100), L(1) {}
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
template<typename G, typename F, typename C>
//template<typename G, typename F>
//, C constraints
double symp_euler_step(G &g, double t, double dt, F force, C constraints) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {

        auto n = *it;
//        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){ continue;}
        n.position() += n.value().vel * dt;
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt

    }
    constraints(g, t);
    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
//        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){ continue;}
        n.value().vel += force(n, t) * (dt / n.value().mass);
        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
//        n.value().vel += force(n, t) * (dt / n.value().mass);
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
        (void)t;
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
            return Point(0, 0, 0);
        } else {

            Point Fg = (n.value().mass) * Point(0, 0, -grav);
            Point Fs = Point(0, 0, 0);
            for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
                auto n2 = (*e).node2();
                DEBUG_MSG("Here in " << n.index() << " and " << n2.index());
                auto dist = n.position() - n2.position();
                Fs = Fs - ((*e).value().K) * (dist) * (norm(dist) - (*e).value().L) / norm(dist);
            }
            return Fg + Fs;
        }


    }
};


// MAKE MASS FORCE
struct MassSpringForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        Point Fs = Point(0, 0, 0);
        for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
            auto n2 = (*e).node2();
            auto dist = n.position() - n2.position();
            Fs -= ((*e).value().K) * (dist) * (norm(dist) - (*e).value().L) / norm(dist);
        }
        return Fs;
    }
};


// MAKE GRAVITY FORCE
struct GravityForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        return (n.value().mass) * Point(0, 0, -grav);
    }
};


// Make Damping Force
struct DampingForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void)t;
        auto c = -1.0 / N;
        return c * (n.value().vel);
    }
};


// Make Paired Force
template<typename FORCE1, typename FORCE2>
struct PairedForce {
    FORCE1 force1;
    FORCE2 force2;
    template<typename NODE>
    Point operator()(NODE n, double t) {
        return force1(n, t) + force2(n, t);
    }
};


//  Paired Force Fucntion
template<typename FORCE1, typename FORCE2>
PairedForce<FORCE1, FORCE2> make_combined_force(FORCE1 force1, FORCE2 force2) {
    return {force1, force2};
};


//  Three force Force Fucntion
template<typename FORCE1, typename FORCE2, typename DRAKE>
PairedForce<FORCE1, PairedForce<FORCE2, DRAKE>> make_combined_force(FORCE1 force1, FORCE2 force2, DRAKE force3) {
    return {force1, {force2, force3}};
};


// Sphere Constraint One
struct SphereConstraint{
    Point center = Point(.5, .5, -.5);
    double radius = .15;
    template<typename G>
    void operator()(G& g, double t){
        (void)t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            auto dist = n.position() - center;
            if (norm(dist)< radius){
                auto R = (dist / norm(dist));
                n.position() = center + radius*R;
                n.value().vel = n.value().vel -  dot(n.value().vel, R)*R;
            }
        }

    }
};


// Plane Constraint
struct PlaneConstraint{
    template<typename G>
    void operator()(G& g, double t){
        (void)t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (n.position().z < -.75){
                n.value().vel.z = 0.0;
                n.position().z = -.75;
            }

        }

    }
};


// Sphere Constraint Two
struct SphereConstraintNew{
    Point center = Point(.5, .5, -.5);
    double radius = .15;
    template<typename G>
    void operator()(G& g, double t){
        (void)t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            auto dist = n.position() - center;
            if (norm(dist)< radius){
                g.remove_node(n);
            }
        }
    }
};


// Fixed Constraint
struct FixedConstraint{
    template<typename G>
    void operator()(G& g, double t){
        (void)t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (n.value().fixed == 1){
                n.position() = Point(0,0,0);
                n.value().vel = Point(0,0,0);
            }
            if (n.value().fixed == 2){
                n.position() = Point(1,0,0);
                n.value().vel = Point(0,0,0);
            }
        }
    }
};


// Paired Constraint
template<typename CONSTRAINT1, typename CONSTRAINT2>
struct PairedConstraint {
    CONSTRAINT1 gods;
    CONSTRAINT2 plan;
    template<typename G>
    void operator()(G& g, double t) {
        gods(g,t);
        plan(g,t);

    }
};


// Paired Constraint Function
template<typename CONSTRAINT1, typename CONSTRAINT2>
PairedConstraint<CONSTRAINT1, CONSTRAINT2> make_combined_constraint(CONSTRAINT1 c1, CONSTRAINT2 c2) {
    return {c1, c2};
};


// Triple Constraint Function
template<typename CONSTRAINT1, typename CONSTRAINT2, typename CONSTRAINT3>
PairedConstraint<CONSTRAINT1, PairedConstraint<CONSTRAINT2, CONSTRAINT3>> make_combined_constraint(CONSTRAINT1 c1, CONSTRAINT2 c2, CONSTRAINT3 c3) {
    return {c1, {c2, c3}};
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

        // Diagonal edges: include as of HW2 #2

        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);

        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }

    // HW2 #1 YOUR CODE HERE
    // Set initial conditions for your nodes, if necessary.

    N = graph.num_nodes();
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        auto n = *it;
        n.value().mass =1.0/N;
        if (n.position() == Point(0,0,0)){
            n.value().fixed = 1;
        }
        if (n.position() == Point(1,0,0)){
            n.value().fixed = 2;
        }
        for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
            (*e).value().L = (*e).length();
        }
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
            auto forces = make_combined_force(GravityForce(), MassSpringForce(), DampingForce());
//            auto constraints = make_combined_constraint(FixedConstraint(), SphereConstraint(), PlaneConstraint());
            auto constraints = make_combined_constraint(PlaneConstraint(),FixedConstraint(), SphereConstraintNew());
            symp_euler_step(graph, t, dt, forces, constraints);

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
