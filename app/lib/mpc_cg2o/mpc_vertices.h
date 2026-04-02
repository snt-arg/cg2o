#ifndef  G2O_GRAPH_VERTEXES_MPC_H_
#define  G2O_GRAPH_VERTEXES_MPC_H_

#include <g2o/core/base_vertex.h>
#include <Eigen/Core>

namespace cg2o::mpc {


// Traits for defining the Estimate type
template <int D>
struct EstimateTraits {
    using Type = Eigen::Matrix<double, D, 1>;
};

// Specialization for D = 1
template <>
struct EstimateTraits<1> {
    using Type = double;
};

// VertexMPC class template declaration
template <int D>
class VertexMPC : public g2o::BaseVertex<D, typename EstimateTraits<D>::Type> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Base = g2o::BaseVertex<D, typename EstimateTraits<D>::Type>;
    using EstimateType = typename EstimateTraits<D>::Type;

    // Set the vertex state to the origin
    virtual void setToOriginImpl() override;

    // Apply an increment (delta) to the state
    virtual void oplusImpl(const double* update) override;

    // Write the vertex state to an output stream
    virtual bool write(std::ostream& os) const override;

    // Read the vertex state from an input stream
    virtual bool read(std::istream& is) override;
};

}  // namespace cg2o::mpc

# include "mpc_vertices.hpp"


// List the need vertices 
class Vertex_scalar : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for longitudinal velocity

class Vertex_v_h : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for longitudinal velocity
class Vertex_d_h : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for longitudinal distance
class Vertex_f_t : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for traction force
class Vertex_f_b : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for braking force
class Vertex_slack_1 : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for slack variable 1
class Vertex_slack_2 : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for slack variable 2
class Vertex_slack_3 : public cg2o::mpc::VertexMPC<1> { };   // 1D vertex for slack variable 3

 
#endif  // G2O_GRAPH_VERTEXES_MPC_H_

