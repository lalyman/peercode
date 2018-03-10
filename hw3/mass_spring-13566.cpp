/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
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

//< Sphere dimensions.
const double R = 0.15;
const Point REFONE = Point( 0, 0, 0 );
const Point REFTWO = Point( 0.5, 0.5, -0.5 );

//< Gravity in meters/sec^2
static constexpr double grav = 9.81;
double N;

/** Custom structure of data to store with Edges */
struct EdgeData {
    double L;                        //< Length of an edge
    double K;                        //< Spring const. associated with an edge
    EdgeData() : L(1.0), K(100.0) {} //< Constructor
};

/** Custom structure of data to store with Nodes */
struct NodeData {
    int pin;                                //< Fixed position of a vertex.
    Point vel;                              //< Velocity associated with the vertex.
    double mass;                            //< Mass associated with the vertex.
    NodeData() : pin(0), vel(0), mass(1) {} //< Constructor
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
template< typename G, typename F, typename Constraint >
double symp_euler_step( G& g, double t, double dt, F force, Constraint constraint ) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }
    
    // Enforcing constraints before updating velocity.
    constraint( g );
    
    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        
        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n) * (dt / n.value().mass);
    }
    return t + dt;
}

//
// FORCES
//
/* Gravitational force on @a vertex.
 */
struct GravityForce {
    template <typename NODE>
    Point operator()( NODE vertex ) { return Point( 0, 0, -grav ) * vertex.value().mass; }
};

/* Spring force on @a vertex.
 */
struct SpringForce {
    template <typename NODE>
    Point operator()( NODE vertex ) {
        Point forceSpring = REFONE;
        for(auto i = vertex.edge_begin(); i != vertex.edge_end(); ++i) {
            auto tail = (*i).node2();
            Point tailDistance = vertex.position() - tail.position();
            
            //< Updating the value of the force.
            forceSpring -= ( norm( tailDistance ) - ( *i ).value().L ) / norm( tailDistance )
            * tailDistance * ( *i ).value().K;
        }
        return forceSpring;
    }
};

/* Damping force on @a vertex.
 */
struct DampingForce {
    template <typename NODE>
    Point operator()( NODE vertex ) { return (-1) * vertex.value().vel / N; }
};

//< Two forces stuct.
template< typename ForceOne, typename ForceTwo >
struct TwoCombinedForces {
    ForceOne fOne;
    ForceTwo fTwo;
    
    template< typename Node >
    Point operator()( Node vertex ) { return fOne( vertex ) + fTwo( vertex ); }
};

//<
//< Handles combined forces.
//<
template< typename ForceOne, typename ForceTwo >
TwoCombinedForces< ForceOne, ForceTwo > make_combined_forces( ForceOne fOne, ForceTwo fTwo ) {
    return { fOne, fTwo };
}

template< typename ForceOne, typename ForceTwo, typename ForceThree >
TwoCombinedForces< TwoCombinedForces < ForceOne, ForceTwo >, ForceThree >
make_combined_forces( ForceOne fOne, ForceTwo fTwo, ForceThree fThree ) {
    return { { fOne, fTwo }, fThree };
}

//< Multiple constraints.
template< typename ConstraintOne, typename ConstraintTwo >
struct MultipleConstraints {
    ConstraintOne cOne;
    ConstraintTwo cTwo;
    template< typename G >
    void operator()( G& g ) {
        cOne( g );
        cTwo( g );
    }
};

//< Two constraints functions.
template< typename ConstraintOne, typename ConstraintTwo >
MultipleConstraints< ConstraintOne, ConstraintTwo >
make_combined_constraints( ConstraintOne cOne, ConstraintTwo cTwo ) { return { cOne, cTwo }; }

//< Three constraints functions.
template< typename ConstraintOne, typename ConstraintTwo, typename ConstraintThree >
MultipleConstraints< MultipleConstraints< ConstraintOne, ConstraintTwo >, ConstraintThree >
make_combined_constraints( ConstraintOne cOne, ConstraintTwo cTwo, ConstraintThree cThree) {
    return { { cOne, cTwo }, cThree };
}

//< Four constraints functions.
template< typename ConstraintOne, typename ConstraintTwo, typename ConstraintThree,
typename ConstraintFour > MultipleConstraints < MultipleConstraints < ConstraintOne,
ConstraintTwo >, MultipleConstraints< ConstraintThree, ConstraintFour > >
make_combined_constraints(ConstraintOne cOne, ConstraintTwo cTwo, ConstraintThree cThree,
                          ConstraintFour cFour) {
    return { { cOne, cTwo }, { cThree, cFour } };
}

//<
//< CONSTRAINTS
//<

/* Fixed poisiton constraint @a vertex.
 */
struct FixedPositionConstraint {
    template< typename G >
    void operator()( G& g ) {
        for( auto i = g.node_begin(); i != g.node_end(); ++i ) {
            Node vertex = *i;
            
            //< Executes when vertex is constrained.
            switch( vertex.value().pin ){
                case ( -1 ):
                    vertex.position() = REFONE;
                    vertex.value().vel = REFONE;
                    break;
                case 1:
                    vertex.position() = Point(1,0,0);
                    vertex.value().vel = REFONE;
                    break;
            }
        }
    }
};

/* Plane constraint @a vertex.
 */
struct PlaneConstraint {
    template< typename G >
    void operator()( G& g ) {
        for( auto i = g.node_begin(); i != g.node_end(); ++i ) {
            Node vertex = *i;
            if( vertex.position().z < - 0.75 ) {
                vertex.value().vel = REFONE;
                vertex.position().z = - 0.75;
            }
        }
    }
};

/* Transparent sphere constraint @a vertex.
 */
struct SphereConstraintTwo {
    template< typename G >
    void operator()( G& g ) {
        for( auto idx = g.node_begin(); idx != g.node_end(); ++idx ) {
            Node vertex = *idx;
            if( norm( vertex.position() - REFTWO ) < R ) g.remove_node( vertex );
        }
    }
};

/* Stiff sphere constraint @a vertex.
 */
struct SphereConstraint {
    template< typename G >
    void operator()( G& g ) {
        for( auto i = g.node_begin(); i != g.node_end(); ++i ) {
            Node vertex = *i;
            Point unit = ( vertex.position() - REFTWO ) / norm( vertex.position() - REFTWO );
            if( norm( vertex.position() - REFTWO ) < R ) {
                vertex.position() = REFTWO + R * unit;
                vertex.value().vel = vertex.value().vel - unit * ( vertex.value().vel * unit );
            }
        }
    }
};

int main(int argc, char** argv) {
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
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);
        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
    
    //
    // INITIALIZATION
    //
    N = graph.num_nodes(); //< Number of vertices in a graph.
    for (auto idx = graph.node_begin(); idx != graph.node_end(); ++idx) {
        (*idx).value().mass = 1.0 / N; //< Adds mass by dereferencing a pointer.
        if( (*idx).position() == Point( 1, 0, 0 ) ) (*idx).value().pin = 1;
        else if( (*idx).position() == REFONE ) (*idx).value().pin = -1;
        for( auto i = (*idx).edge_begin(); i != (*idx).edge_end(); ++i ) {
            (*i).value().L = (*i).length();
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
    auto sim_thread = std::thread([&](){
        
        // Begin the mass-spring simulation
        double dt = 0.001;
        double t_start = 0;
        double t_end = 5.0;
        
        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            //std::cout << "t = " << t << std::endl;
            
            auto forcesVar = make_combined_forces( GravityForce(), SpringForce(), DampingForce() );
            auto constraintsVar = make_combined_constraints( FixedPositionConstraint(),  PlaneConstraint(), SphereConstraintTwo(), SphereConstraint() );
            
            symp_euler_step(graph, t, dt, forcesVar, constraintsVar );
            viewer.clear();   //< Clears the viewer.
            node_map.clear(); //< Clears nodes & edges.
            
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
