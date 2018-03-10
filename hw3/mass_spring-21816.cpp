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
    double L;       //< Node velocity
    double K;     //< Node mass
    EdgeData() : L(1), K(1) {}
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
template<typename G, typename F>
double symp_euler_step(G &g, double t, double dt, F force) {
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
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
            continue;
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}

template<typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint) {
    // Compute the t+dt position

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // Compute the constraints
    constraint(g, t);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
            continue;
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
        (void) t;
        double mass = n.value().mass;
        Point gravity(0, 0, -mass * grav);

        Point spring(0, 0, 0);
        for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
            NODE n_neighbor = (*e).node2();
            Point diff_pos = n_neighbor.position() - n.position();
            double K = (*e).value().K;
            double L = (*e).value().L;

            spring += K * diff_pos / norm(diff_pos) * (norm(diff_pos) - L);
        }

        return gravity + spring;
    }
};

// Structure representing the sum of two forces f1 and f2.
template<typename F1, typename F2>
struct PairForces{
    PairForces(F1 f1, F2 f2):f1_(f1), f2_(f2){
    }

    template<typename NODE>
    Point operator()(NODE n, double t){
        return f1_(n, t) + f2_(n, t);
    }

    F1 f1_;
    F2 f2_;
};

// Function to help build a summed force f1+f2.
template<typename F1, typename F2>
PairForces<F1, F2> make_combined_forces(F1 f1, F2 f2){
    return PairForces<F1, F2>(f1, f2);
};

// Structure representing the sum of three forces f1, f2 and f3.
template<typename F1, typename F2, typename F3>
struct TripleForces{
    TripleForces(F1 f1, F2 f2, F3 f3):f1_(f1), f2_(f2), f3_(f3){
    }

    template<typename NODE>
    Point operator()(NODE n, double t){
        return f1_(n, t) + f2_(n, t) + f3_(n, t);
    }

    F1 f1_;
    F2 f2_;
    F3 f3_;
};

// Function to help build a summed forces f1+f2+f3.
template<typename F1, typename F2, typename F3>
TripleForces<F1, F2, F3> make_combined_forces(F1 f1, F2 f2, F3 f3){
    return TripleForces<F1, F2, F3>(f1, f2, f3);
};


/** Gravity Force function object. */
struct GravityForce{
    /** Return the gravity force applying to @a n at time @a t.
     *
    */

    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        double mass = n.value().mass;
        return Point(0, 0, -mass * grav);
    }
};


/** Mass Spring Force function object. */
struct MassSpringForce{
    /** Return the mass spring force applying to @a n at time @a t. It is modeled by an elastic link between every
     * connected nodes, with a constant K and an idle length L.
     */

    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point spring(0, 0, 0);
        (void) t;

        for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
            Edge current_edge = (*e);
            Node n_neighbor = current_edge.node2();
            Point diff_pos = n_neighbor.position() - n.position();
            double K = current_edge.value().K;
            double L = current_edge.value().L;
//
            spring += K * diff_pos * (1 - L / norm(diff_pos));
        }

        return spring;
    }
};


/** Damping Force function object. */
struct DampingForce{
    /** Return the dumping force applying to @a n at time @a t.
     * Models friction to the nodes.
     */
    DampingForce(double c):c(c){
    }
    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        Point velocity = n.value().vel;

        Point damping = -c*velocity;
        return damping;
    }

    double c;
};

// Allow to build a single constraint being the combionation of the two constraints f1 and f2.
template<typename C1, typename C2>
struct PairConstraints{
    PairConstraints(C1 c1, C2 c2):c1_(c1), c2_(c2){
    }

    template<typename G>
    void operator()(G& g, double t){
        c1_(g, t);
        c2_(g, t);
    }

    C1 c1_;
    C2 c2_;
};

// Allow to build a constraint being the the total of the two constraints c1 and c2.
template<typename C1, typename C2>
PairConstraints<C1, C2> make_combined_constraints(C1 c1, C2 c2){
    return PairConstraints<C1, C2>(c1, c2);
};

// Constraint representing a solid plan.
struct PlaneConstraint{
    PlaneConstraint(double p): z_plane(p){
    }

    template<typename G>
    void operator()(G& graph, double t){
        (void) t;
        for(auto n = graph.node_begin(); n!= graph.node_end(); ++n){
            if((*n).position().z<z_plane){
                (*n).position().z = z_plane;
                (*n).value().vel.z = 0;
            }
        }
    }

    double z_plane;
};

// Constraint representing a solid sphere.
struct SphereConstraint{
    SphereConstraint(Point c, double r): center(c), radius(r){
    }

    template<typename G>
    void operator()(G& graph, double t){
        (void) t;
        for(auto n = graph.node_begin(); n!= graph.node_end(); ++n){
            if(norm((*n).position()-center)<radius){
                Point diff_to_radius = ((*n).position()-center);
                Point R = (diff_to_radius)/norm(diff_to_radius);
                (*n).position() = R*radius+center;
                (*n).value().vel -= (inner_prod((*n).value().vel, R))*R;
            }
        }
    }

    Point center;
    double radius;
};

// Constraint representing a plan which removes touched nodes.
struct PlaneRemovalConstraint{
    PlaneRemovalConstraint(double p): z_plane(p){
    }

    template<typename G>
    void operator()(G& graph, double t){
        (void) t;
        for(auto n = graph.node_begin(); n!= graph.node_end(); ){
            if((*n).position().z<z_plane){
                n = graph.remove_node(n);
            }
            else{
                ++n;
            }
        }
    }

    double z_plane;
};

// Constraint representing a sphere which removes touched nodes.
struct SphereRemovalConstraint{
    SphereRemovalConstraint(Point c, double r): center(c), radius(r){
    }

    template<typename G>
    void operator()(G& graph, double t){
        (void) t;
        for(auto n = graph.node_begin(); n!= graph.node_end(); ){
            if(norm((*n).position()-center)<radius){
                n = graph.remove_node(n);
            }
            else{
                ++n;
            }
        }
    }

    Point center;
    double radius;
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

    int number_nodes = graph.size();

    // Set initial conditions for your nodes, if necessary.
    for(auto n = graph.node_begin(); n!= graph.node_end(); ++n){
        (*n).value().mass = 1/float(number_nodes);
    }

    for(auto e = graph.edge_begin(); e!= graph.edge_end(); ++e){
        (*e).value().K = 100;
        (*e).value().L = (*e).length();
    }

    // Print out the stats
    std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

    // Launch the Viewer
    CME212::SFML_Viewer viewer;
    auto node_map = viewer.empty_node_map(graph);

    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    viewer.center_view();

    GravityForce gravity_force;
    MassSpringForce spring_force;
    DampingForce damping_force(1/float(number_nodes));
    auto total_force = make_combined_forces(gravity_force, spring_force, damping_force);

    SphereConstraint sphere_constraint(Point(0.5, 0.5, -0.5), 0.15);
    SphereRemovalConstraint sphere_removal_constraint(Point(0.5, 0.5, -0.5), 0.15);
    PlaneConstraint plane_constraint(-0.75);
    PlaneRemovalConstraint plane_removal_constraint(-0.75);
    auto total_constraint = make_combined_constraints(plane_constraint, sphere_removal_constraint);


    // We want viewer interaction and the simulation at the same time
    // Viewer is thread-safe, so launch the simulation in a child thread
    bool interrupt_sim_thread = false;
    auto sim_thread = std::thread([&]() {

        // Begin the mass-spring simulation
        double dt = 0.001;
        double t_start = 0;
        double t_end = 5.0;

        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            symp_euler_step(graph, t, dt, total_force, total_constraint);

            viewer.clear();
            node_map.clear();

            // Update viewer with nodes' new positions
            viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
            viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
            viewer.set_label(t);

            // These lines slow down the animation for small graphs, like grid0_*.
            // Feel free to remove them or tweak the constants.
//            if (graph.size() < 100)
//                std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

    });  // simulation thread

    viewer.event_loop();

    // If we return from the event loop, we've killed the window.
    interrupt_sim_thread = true;
    sim_thread.join();

    return 0;
}
