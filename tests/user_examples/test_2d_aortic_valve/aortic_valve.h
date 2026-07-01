/**
 * @file aortic_valve.h
 * @brief Case setup for the 2D aortic valve tutorial example.
 */

#ifndef AORTIC_VALVE_CASE_H
#define AORTIC_VALVE_CASE_H

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
// Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.16;                           /**< Channel length. */
Real DH = 0.02;                           /**< Channel height. */
Real R = 0.02;                            /**< Sinus cavity radius. */
Real Leaflet_L = 0.024;                   /**< Leaflet length. */
Real Leaflet_T = 2.4e-3;                  /**< Leaflet thickness. */
int sinus_resolution = 180;               /**< Point count for the sinus cavity arc. */
Vec2d insert_circle_center(0.5 * DL, 0.5 * DH);
Real resolution_ref = Leaflet_T / 4.0;    /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0;   /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;           /**< Boundary width. */
Real Leaflet_base_length = BW;            /**< Length of constrained leaflet. */

//----------------------------------------------------------------------
// Material properties.
//----------------------------------------------------------------------
Real rho0_f = 1.0e3;                      /**< Fluid density. */
Real U_f = 0.13889;                       /**< Characteristic fluid velocity. */
Real c_f = 10.0 * U_f;                    /**< Speed of sound. */
Real mu_f = 4.3e-3;                       /**< Dynamic viscosity. */
Real rho0_s = 1.0e3;                      /**< Solid reference density. */
Real poisson = 0.49;                      /**< Poisson ratio. */
Real Youngs_modulus = 1.5e5;              /**< Young's modulus. */

//----------------------------------------------------------------------
// Geometry helpers.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.5 * DH));
    for (int i = 0; i < sinus_resolution + 1; ++i)
    {
        water_block_shape.push_back(Vecd(insert_circle_center[0] - R * cos(i * Pi / sinus_resolution),
                                         insert_circle_center[1] + R * sin(i * Pi / sinus_resolution)));
    }
    water_block_shape.push_back(Vecd(DL, 0.5 * DH));
    water_block_shape.push_back(Vecd(DL, -0.5 * DH));
    water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));

    return water_block_shape;
}

std::vector<Vecd> createOuterWallShape1()
{
    std::vector<Vecd> outer_wall_shape1;
    outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
    outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, 0.5 * DH + BW));
    outer_wall_shape1.push_back(Vecd(DL + BW, 0.5 * DH + BW));
    outer_wall_shape1.push_back(Vecd(DL + BW, -0.5 * DH - BW));
    outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

    return outer_wall_shape1;
}

std::vector<Vecd> createOuterWallShape2()
{
    std::vector<Vecd> outer_wall_shape2;
    for (int i = 0; i < sinus_resolution + 1; ++i)
    {
        outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - (R + BW) * cos(i * Pi / sinus_resolution),
                                         insert_circle_center[1] + (R + BW) * sin(i * Pi / sinus_resolution)));
    }
    outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - R - BW, insert_circle_center[1]));

    return outer_wall_shape2;
}

std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
    inner_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH));
    inner_wall_shape.push_back(Vecd(DL + BW, 0.5 * DH));
    inner_wall_shape.push_back(Vecd(DL + BW, -0.5 * DH));
    inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

    return inner_wall_shape;
}

Vec2d BLB(0.5 * DL - R - (Leaflet_base_length + 0.5 * Leaflet_T) * cos(0.25 * Pi),
          0.5 * DH + (Leaflet_base_length - 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BLT(0.5 * DL - R - (Leaflet_base_length - 0.5 * Leaflet_T) * cos(0.25 * Pi),
          0.5 * DH + (Leaflet_base_length + 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BRT(0.5 * DL - R + (Leaflet_L + 0.5 * Leaflet_T) * cos(0.25 * Pi),
          0.5 * DH - (Leaflet_L - 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BRB(0.5 * DL - R + (Leaflet_L - 0.5 * Leaflet_T) * cos(0.25 * Pi),
          0.5 * DH - (Leaflet_L + 0.5 * Leaflet_T) * sin(0.25 * Pi));

std::vector<Vecd> createLeafletShape()
{
    std::vector<Vecd> leaflet_shape;
    leaflet_shape.push_back(BLB);
    leaflet_shape.push_back(BLT);
    leaflet_shape.push_back(BRT);
    leaflet_shape.push_back(BRB);
    leaflet_shape.push_back(BLB);

    return leaflet_shape;
}

Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, -0.5 * DH) + buffer_halfsize;

namespace SPH
{
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(createWaterBlockShape(), GeometricOps::add);
        multi_polygon_.addPolygon(createLeafletShape(), GeometricOps::sub);
    }
};

class InsertedBody : public MultiPolygonShape
{
  public:
    explicit InsertedBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(createOuterWallShape1(), GeometricOps::add);
        multi_polygon_.addPolygon(createOuterWallShape2(), GeometricOps::add);
        multi_polygon_.addPolygon(createInnerWallShape(), GeometricOps::sub);
        multi_polygon_.addPolygon(createWaterBlockShape(), GeometricOps::sub);
        multi_polygon_.addPolygon(createLeafletShape(), GeometricOps::add);
    }
};

MultiPolygon createLeafletBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addPolygon(createOuterWallShape1(), GeometricOps::add);
    multi_polygon.addPolygon(createOuterWallShape2(), GeometricOps::add);
    multi_polygon.addPolygon(createInnerWallShape(), GeometricOps::sub);
    multi_polygon.addPolygon(createWaterBlockShape(), GeometricOps::sub);
    multi_polygon.addPolygon(createLeafletShape(), GeometricOps::sub);
    return multi_polygon;
}

struct InflowVelocity
{
    Real u_ref_, t_ref_;
    OrientedBox &oriented_box_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(1.0), oriented_box_(boundary_condition.getOrientedBox()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = u_ref_ * 0.5 * (1.0 + sin(Pi * current_time / t_ref_ - 0.5 * Pi));
        if (oriented_box_.checkInBounds(position))
        {
            target_velocity[0] = (-6.0 * position[1] * position[1] / DH / DH + 1.5) * u_ave;
            target_velocity[1] = 0.0;
        }
        return target_velocity;
    }
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    size_t number_observation_points = 21;
    Real range_of_measure = DH - resolution_ref * 4.0;
    Real start_of_measure = resolution_ref * 2.0 - 0.5 * DH;
    for (size_t i = 0; i < number_observation_points; ++i)
    {
        Vec2d point_coordinate(0.0, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
}
} // namespace SPH

#endif // AORTIC_VALVE_CASE_H
