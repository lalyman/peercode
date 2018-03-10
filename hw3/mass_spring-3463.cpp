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
    Point vel;          // Node velocity
    double mass;        // Node mass
    NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double rest_length; // rest length
    double K;           //spring constant
    EdgeData() : rest_length(-1), K(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g          Graph
 * @param[in]     t          The current time (useful for time-dependent forces)
 * @param[in]     dt         The time step
 * @param[in]     force      Function object defining the force per node
 * @param[in]     constraint Function object defining the constraints per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports .vel and
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C& constraint) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }
    constraint.reset_violations(g, t);

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
        Point operator()(NODE n, double t) const {
            (void)t;
            // HW2 #1: YOUR CODE HERE
            if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
                return Point(0,0,0);

            Point force(0, 0, -grav * n.value().mass);
            for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
                Edge e = *ii;
                double L = e.value().rest_length;
                double K = e.value().K;
                Point delta = n.position() - e.node2().position();
                Point f_spring = -K * delta + K*L * delta / e.length();

                force += f_spring;
            }
            return force;
        }
};

// Gravitational force object
struct GravityForce {
    // Return the gravitational force applying to @param n at time @param t
    template <typename NODE>
        Point operator()(NODE n, double t) const {
            (void)t;
            return Point(0, 0, -grav * n.value().mass);
        }
};

// Mass spring force object
struct MassSpringForce {
    // Return the cumulative spring force applying to @param n at time @param t
    template <typename NODE>
        Point operator()(NODE n, double t) const {
            (void)t;

            Point force(0,0,0);
            for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
                Edge e = *ii;
                double L = e.value().rest_length;
                double K = e.value().K;

                Point delta = n.position() - e.node2().position();
                Point f_spring = -K * delta + K * L * delta / e.length();
                force += f_spring;
            }
            return force;
        }
};

// Damping force object
struct DampingForce {
    // Return the damping force applying to @param n at time @param t
    template <typename NODE>
        Point operator()(NODE n, double t) const {
            (void)t;
            return -c_* n.value().vel;
        }

    DampingForce(double c = .01) : c_(c) {};
    double c_; // Damping constant
};

// Constant no-force object
struct NoForce {
    // Returns Point representing no force
    template <typename NODE>
        Point operator()(NODE n, double t) const {
            (void)t;
            (void)n;
            return Point(0,0,0);
        };
};

// Force object supporting up to three types of forces
template <typename F, typename G, typename H = NoForce>
struct Force {
    template <typename NODE>
        Point operator()(NODE n, double t) const {
            (void) t;
            return f1_(n,t) + f2_(n,t) + f3_(n,t);
        };

    Force(F f1, G f2, H f3) : f1_(f1), f2_(f2), f3_(f3) {};

    Force(F f1, G f2) : f1_(f1), f2_(f2), f3_(NoForce()) {};

    F f1_;
    G f2_;
    H f3_;
};

/* @brief Creates combined force of @a f1 and @a f2
 * @param[f1] functor which supports ()(NODE n, double t)
 * @param[f2] functor which supports ()(NODE n, double t)
 * @tparam[F,G] Functor suppporting ()(NODE n, double t) for same type NODE
 * @return force object representing the sum of the forces f1 and f2
 */
template <typename F, typename G>
Force<F,G> make_combined_force(F f1, G f2) {
    Force<F,G> f(f1, f2);
    return f;

};

/* @brief Creates combined force of @a f1, @a f2, and @a f3
 * @param[f1] functor which supports ()(NODE n, double t)
 * @param[f2] functor which supports ()(NODE n, double t)
 * @param[f3] functor which supports ()(NODE n, double t)
 * @tparam[F,G,H] Functor suppporting ()(NODE n, double t) for same type NODE
 * @return force object representing the sum of the forces f1, f2, and f3
 */
template <typename F, typename G, typename H>
Force<F,G,H> make_combined_force(F f1, G f2, H f3) {
    Force<F,G, H> f(f1, f2, f3);
    return f;
};

// Constraint on the z-plane.
class ZPlaneConstraint {

    public:
        ZPlaneConstraint(double w) : w_(w) {};

    // @brief Returns whether or node @a n violates the z-plane constraint
    // Complexity O(1)
    template <typename NODE>
        bool violates_constraint(NODE n) {
            return n.position().z < w_;
        };

    // @brief Sets @a n z-position to constraint value w_, and z-velocity to 0
    // Complexity O(1)
    template <typename NODE>
        void fix_constraint(NODE n) {
            n.position().z = w_;
            n.value().vel.z = 0;
        };

    /* @brief Fixes any nodes violating the constraint
     * @param[g] Graph object
     * @param[t] Time at which constraints imposted
     * @post !violates_constraint(n) for any node n in g
     *
     * Complexity O(g.num_nodes())
     */
    template <typename G>
        void reset_violations(G& g, double t) {
            (void) t;
            for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
                auto n = *ni;
                if (violates_constraint(n)) {
                    fix_constraint(n);
                }
            }
        };

    double w_; // z-value representing constraint
};

// Constraint in the shape of a sphere
// When a node enters the sphere, it is reset by being placed
// on the closest boundary point of the sphere, with 0 normal
// velocity to the sphere surface.
struct SphereConstraintAvoid {
    template <typename NODE>
        // Returns true if node @a n is within distance r_ of c_
        // Complexity O(1)
        bool violates_constraint(NODE n) {
            return norm(n.position() - c_) < r_;
        }

    /* @param[n] Node violating constraint
     * @brief Set the position to the nearest point on the surface of the sphere
     * Set the component of the velocity that is normal to the sphere’s surface to zero:
     *
     * Complexity O(1)
     */
    template <typename NODE>
        void fix_constraint(NODE n) {
            Point delta = n.position() - c_;
            n.position() = c_ + (delta) * r_ / norm(delta);
            Point Ri = delta / norm(delta);

            Point ip = inner_prod(n.value().vel, Ri) * Ri;
            n.value().vel = n.value().vel - ip;
        }

    /* @brief Fixes any nodes violating the constraint
     * @param[g] Graph object
     * @param[t] Time at which constraints imposted
     * @post !violates_constraint(n) for any node n in g
     *
     * Complexity O(g.num_nodes())
     */
    template <typename G>
        void reset_violations(G& g, double t) {
            (void) t;
            for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
                auto n = *ni;
                if (violates_constraint(n)) {
                    fix_constraint(n);
                }
            }
        };

    SphereConstraintAvoid(Point c, double r) : c_(c), r_(r) {};

    Point c_;  // Center of sphere
    double r_; // Radius of spehre

};

// Constraint in the shape of a sphere
// When a node enters the sphere, it is reset by being removed
struct SphereConstraintRemove {
    // Returns true if node @a n is within distance r_ of c_
    // Complexity O(1)
    template <typename NODE>
        bool violates_constraint(NODE n) {
            return norm(n.position() - c_) < r_;
        }

    /* @brief Removes node n
     * @param[g] Graph object containing node n
     * @param[n_it] Node iterator at a node violating the  constraint
     * @post !g.num_nodes() is decremented by 1
     * @post n_it has been incremented
     *
     * Removes node n, and increments n_it
     *
     * Complexity O(g.num_nodes())
     */
    template <typename G, typename NODE_IT>
        void fix_constraint(G& g, NODE_IT& n_it) {
            n_it = g.remove_node(n_it);
        }

    /* @brief Fixes any nodes violating the constrain
     * @param[g] Graph object
     * @param[t] Time at which constraints imposted
     * @post !violates_constraint(n) for any node n in g
     *
     * Complexity O(g.num_nodes()**2)
     */
    template <typename G>
        void reset_violations(G& g, double t) {
            (void) t;
            for (auto ni = g.node_begin(); ni != g.node_end(); ) {
                auto n = *ni;
                if (violates_constraint(n)) {
                    fix_constraint(g, ni);
                } else {
                    ++ni;
                }
            }
            // Prune any loose nodes that may have occurred, due to numerical issues
            // or just bad luck
            for (auto ni = g.node_begin(); ni != g.node_end(); ) {
                auto n = *ni;
                if (n.degree() <= 1) {
                   g.remove_node(ni);
                } else {
                    ++ni;
                }
            }
        };

    SphereConstraintRemove(Point c, double r) : c_(c), r_(r) {};

    Point c_;  // Center of sphere
    double r_; // Radius of sphere
};

// Fixed node constraint. Nodes are never moved from their original positions
template <typename G>
struct FixedNodeConstraint {
    using Node = typename G::node_type;

    /* @brief Returns fixed nodes to their stationary positions
     * @param[g] Graph object
     * @param[t] Time at which constraints imposted
     * @post !violates_constraint(n) for any node n in g
     *
     * Complexity O(1) ammortized time
     */
    void reset_violations(G& g, double t) {
        (void) t; // To quiet compiler warning

        for (typename G::size_type i = 0; i < inds_.size(); ++i) {
            auto ind = inds_[i].first;
            Point point = inds_[i].second;
            if (ind > g.num_nodes() - 1) { // Index invalid due to node removal
                for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
                    Node n = *ni;
                    if (n.prev_ind() == ind) {
                        inds_[i].first = n.index();
                        break;
                    }
                }
            }
            g.node(inds_[i].first).position() = point;
        }
    }

    /* @param[g] Graph object
     * @param[points] vector of points. Any node with position in
     * points should remain stationary
     *
     * Complexity O(g.num_nodes())
     */
    FixedNodeConstraint(G& g, std::vector<Point> points) {
        for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
            Node n = *ni;

            for (Point& point : points) {
                if (point == n.position()) {
                    inds_.push_back(std::make_pair(n.index(), point));
                }
            }
        }
    }
    std::vector<std::pair<typename G::size_type, Point>> inds_; // (index, position) of fixed nodes
};

// Constraint object representing no constraint
struct NoConstraint {
    template <typename G>
        /* @brief Does nothing
         * @param[g] Graph object
         * @param[t] Time at which to reset violations
         *
         * Complexity O(1)
         */
        void reset_violations(G& g, double t) {
            (void) g;
            (void) t; // To quiet compiler warning
        };

    NoConstraint() {};
};

// @brief Object representing up to 3 constraints
template <typename C, typename D, typename E = NoConstraint>
struct Constraint {
    /* @brief Applies constraints to each of up to three constraints in order
     * @param[g] Graph object
     * @param[t] Time at which constraints imposted
     *
     * The constraints are applied in order. In particular, the order of
     * resetting violations is not commutative.
     */
    template <typename G>
        void reset_violations(G& g, double t) {
            (void) t;

            c1_.reset_violations(g, t);
            c2_.reset_violations(g, t);
            c3_.reset_violations(g, t);
        };

    Constraint(C c1, D c2, E c3) : c1_(c1), c2_(c2), c3_(c3) {};

    Constraint(C c1, D c2) : c1_(c1), c2_(c2), c3_(NoConstraint()) {};

    C c1_; // First constraint to be enforced
    D c2_; // Second constraint to be enforced
    E c3_; // Third consrtaint to be enforced
};

/* @brief Constructs a constraint object from two constraints
 * @param[c1] First constraint
 * @param[c2] Second constraint
 * @tparam[C,D] Objects supporting reset_violations(G, t) function
 * @return Constraint object representing @a c1 followed by @a c2
 */
template <typename C, typename D>
Constraint<C,D> make_combined_constraint(C c1, D c2) {
    Constraint<C,D> c(c1, c2);
    return c;
};

/* @brief Constructs a constraint object from three constraints
 * @param[c1] First constraint
 * @param[c2] Second constraint
 * @param[c3] Second constraint
 * @tparam[C,D,E] Objects supporting reset_violations(G, t) function
 * @return Constraint object representing @a c1 followed by @a c2 and then @a c3
 */
template <typename C, typename D, typename E>
Constraint<C,D,E> make_combined_constraint(C c1, D c2, E c3) {
    Constraint<C,D,E> c(c1, c2, c3);
    return c;
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
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
        Node node = *ni;
        node.value().vel = Point(0,0,0);
        node.value().mass = 1. / graph.num_nodes();

        for(auto ii = node.edge_begin(); ii != node.edge_end(); ++ii) {
            Edge e = *ii;
            e.value().rest_length = e.length();
            e.value().K = 100;
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

    auto comb_force = make_combined_force(GravityForce(), MassSpringForce(),
                                          DampingForce(1./graph.num_nodes()));

    auto sphere_constraint_avoid = SphereConstraintAvoid(Point(0.5,0.5,-0.5),0.15);
    auto sphere_constraint_remove = SphereConstraintRemove(Point(0.5,0.5,-0.5),0.15);
    (void)sphere_constraint_remove; // To avoid compiler warning
    (void)sphere_constraint_avoid;  // To avoid compiler warning

    auto z_plane_constraint = ZPlaneConstraint(-0.75);

    std::vector<Point> fixed_points = {Point(0,0,0), Point(1,0,0)};

    auto fixed_node_constraint = FixedNodeConstraint<GraphType>(graph, fixed_points);

    //auto comb_constraint = make_combined_constraint(fixed_node_constraint, z_plane_constraint);
    auto comb_constraint = make_combined_constraint(fixed_node_constraint,
                                                    sphere_constraint_remove,
                                                    z_plane_constraint);

    // auto comb_constraint = fixed_node_constraint;

    // We want viewer interaction and the simulation at the same time
    // Viewer is thread-safe, so launch the simulation in a child thread
    bool interrupt_sim_thread = false;
    auto sim_thread = std::thread([&](){
            // Begin the mass-spring simulation
            double dt = 0.001;
            double t_start = 0;
            double t_end = 5.0;

            for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
                viewer.clear();
                node_map.clear();
                symp_euler_step(graph, t, dt, comb_force, comb_constraint);

                // Update viewer with nodes' new positions
                viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
                viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
                viewer.set_label(t);

                // These lines slow down the animation for small graphs, like grid0_*.
                // Feel free to remove them or tweak the constants.
                //if (graph.size() < 100)
                std::this_thread::sleep_for(std::chrono::milliseconds(15));
            }
    });  // simulation thread

    viewer.event_loop();

    // If we return from the event loop, we've killed the window.
    interrupt_sim_thread = true;
    sim_thread.join();

    return 0;
}
