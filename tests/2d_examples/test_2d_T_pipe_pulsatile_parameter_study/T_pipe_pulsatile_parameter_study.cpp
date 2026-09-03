/**
 * @file    T_pipe_pulsatile_parameter_study.cpp
 * @brief   Pulsatile flow in a 2D T-shaped bifurcation - dimensionless parameter study.
 *
 * @details Derived from T_pipe_pulsatile_poiseuille_study.cpp (which itself extends the
 *          T_shaped_pipe benchmark). That file remains the single-case version; this one
 *          is the scriptable study. Differences:
 *
 *            1. The case is defined by three dimensionless groups (Re, alpha, A) supplied
 *               on the command line. omega is DERIVED from alpha, never prescribed, because
 *               Re and alpha are coupled through nu.  See case_params.h.
 *            2. Run length is derived from alpha: the start-up transient decays on the
 *               viscous time tau_visc = R^2/nu, so cycles-to-periodicity ~ alpha^2/(2 pi).
 *               A fixed "5 cycles" is far too short for alpha ~ 10.
 *            3. Observer sampling is decoupled from VTP snapshot writing.
 *            4. Volume flow rate through three cross-sections is reduced directly on the
 *               particles, so Q1 = Q2 + Q3 can actually be evaluated.
 *            5. Wall shear stress is obtained from the SPH velocity-gradient tensor at
 *               near-wall probe lines at two offsets, then linearly extrapolated to the
 *               wall in post-processing -> TAWSS, OSI, RRT.
 *            6. Two geometric extensions: asymmetric lower branch (--branch_ratio) and a
 *               cosine stenosis in the upper branch (--stenosis).
 *
 *          Lecture connections:
 *            - Poiseuille parabolic profile:              Lecture 1
 *            - Mass conservation Q1 = Q2 + Q3:            Lecture 1
 *            - Pulsatility index PI:                      Lecture 3
 *            - Womersley number alpha, Stokes layer:      Lecture 4
 *
 *          Geometry (global coordinates, all lengths relative to DH = 3):
 *            inlet channel   x in [-DL_sponge, DL1], y in [0, DH]      (flow in +x)
 *            vertical duct   x in [DL1, DL],  y in [-b*DH, 2*DH]
 *            upper outlet    y = 2*DH   (flow out in +y)
 *            lower outlet    y = -b*DH  (flow out in -y)
 *          with b = --branch_ratio. Outlet duct width DL - DL1 = 1.5 each, so the two
 *          outlet areas sum exactly to the inlet area DH = 3.
 *
 * @author  Based on T_shaped_pipe by Xiangyu Hu and Shuoguo Zhang, and on
 *          T_pipe_pulsatile_poiseuille_study.cpp. Extended into a parameter study.
 */

#include "sphinxsys.h"

#include "case_params.h"

using namespace SPH;
using TPipeStudy::CaseParams;

//----------------------------------------------------------------------
//  Global case parameters.
//
//  This has to be a file-scope object: InflowVelocity is constructed by template
//  deduction inside fluid_dynamics::InflowVelocityCondition and therefore cannot
//  be given constructor arguments of its own.
//----------------------------------------------------------------------
CaseParams g_params;

//----------------------------------------------------------------------
//  Measurement stations, expressed as fractions of the branch length so that they
//  stay inside the branch for every --branch_ratio. The upper branch spans
//  y in [DH, 2 DH] and the lower y in [-b DH, 0]; each has an outflow disposer
//  occupying the last BW of it, which the stations must clear.
//
//  Both stations sit the SAME fraction of a branch away from the junction, so the
//  volume of fluid between the junction and each measuring station is equal when
//  branch_ratio = 1. That symmetry matters: any fluid stored between the junction
//  and a station appears as a Q1 != Q2 + Q3 imbalance, and unequal stations would
//  bias that imbalance systematically rather than just adding noise to it.
inline Real stationFraction() { return 0.50; }
inline Real upperStation(const CaseParams &p)
{
    return p.DH + stationFraction() * p.DH;
}
inline Real lowerStation(const CaseParams &p)
{
    return -stationFraction() * p.branch_ratio * p.DH;
}
/** How far the outflow disposer band reaches into a branch, plus a safety margin.
 *  The disposer boxes are BW deep in the flow direction (halfsize 0.5*BW, shifted
 *  so they sit flush with the outlet), so probes must stay BW + a margin clear. */
inline Real disposerDepth(const CaseParams &p) { return p.BW + 2.0 * p.dp; }

//----------------------------------------------------------------------
//  Stenosis in the upper branch (extension E).
//
//  A cosine constriction acting symmetrically from both side walls of the upper
//  branch. `stenosis` is the fraction of the branch width removed at the throat,
//  i.e. the 2D equivalent of an area reduction.
//----------------------------------------------------------------------
//  Kept strictly below upperStation() so the flow-rate slab never cuts the throat.
inline Real stenosisY0(const CaseParams &p) { return p.DH + 0.05 * p.DH; }
inline Real stenosisY1(const CaseParams &p) { return p.DH + 0.35 * p.DH; }

/** Depth by which ONE side wall intrudes into the branch at height y. */
inline Real stenosisIntrusion(const CaseParams &p, Real y)
{
    if (p.stenosis <= 0.0)
        return 0.0;
    const Real y0 = stenosisY0(p), y1 = stenosisY1(p);
    if (y <= y0 || y >= y1)
        return 0.0;
    const Real L = y1 - y0, yc = 0.5 * (y0 + y1);
    const Real W = p.DL - p.DL1;
    return 0.5 * p.stenosis * W * 0.5 * (1.0 + cos(2.0 * Pi * (y - yc) / L));
}

/** d(intrusion)/dy, needed for the local wall normal on the stenosed wall. */
inline Real stenosisIntrusionSlope(const CaseParams &p, Real y)
{
    if (p.stenosis <= 0.0)
        return 0.0;
    const Real y0 = stenosisY0(p), y1 = stenosisY1(p);
    if (y <= y0 || y >= y1)
        return 0.0;
    const Real L = y1 - y0, yc = 0.5 * (y0 + y1);
    const Real W = p.DL - p.DL1;
    return -0.5 * p.stenosis * W * (Pi / L) * sin(2.0 * Pi * (y - yc) / L);
}

/**
 * @brief Closed polygon for one side of the stenosis.
 * @param left true for the x = DL1 wall, false for the x = DL wall.
 *
 * The polygon runs along the constricted contour and closes outside the fluid
 * (BW beyond the nominal wall) so that it merges cleanly with the existing wall
 * when added to the solid body and carves cleanly when subtracted from the fluid.
 */
std::vector<Vecd> stenosisPolygon(const CaseParams &p, bool left, int n = 41)
{
    std::vector<Vecd> pts;
    const Real y0 = stenosisY0(p), y1 = stenosisY1(p);
    const Real x_wall = left ? p.DL1 : p.DL;
    const Real x_back = left ? p.DL1 - p.BW : p.DL + p.BW;
    const Real dir = left ? 1.0 : -1.0; // into the fluid

    pts.emplace_back(Vecd(x_back, y0));
    for (int i = 0; i < n; ++i)
    {
        const Real y = y0 + (y1 - y0) * Real(i) / Real(n - 1);
        pts.emplace_back(Vecd(x_wall + dir * stenosisIntrusion(p, y), y));
    }
    pts.emplace_back(Vecd(x_back, y1));
    pts.emplace_back(Vecd(x_back, y0));
    return pts;
}

//----------------------------------------------------------------------
//  Geometry polygons. Same topology as the baseline T_shaped_pipe; the lower
//  branch extent -DH is replaced by -branch_ratio * DH.
//----------------------------------------------------------------------
std::vector<Vecd> waterBlockShape(const CaseParams &p)
{
    const Real Ls = p.DL_sponge, DH = p.DH, DL = p.DL, DL1 = p.DL1;
    const Real yb = -p.branch_ratio * DH;
    return {Vecd(-Ls, 0.0), Vecd(-Ls, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH),
            Vecd(DL, 2.0 * DH), Vecd(DL, yb), Vecd(DL1, yb), Vecd(DL1, 0.0),
            Vecd(-Ls, 0.0)};
}

std::vector<Vecd> outerWallShape(const CaseParams &p)
{
    const Real Ls = p.DL_sponge, DH = p.DH, DL = p.DL, DL1 = p.DL1, BW = p.BW;
    const Real yb = -p.branch_ratio * DH;
    return {Vecd(-Ls - BW, -BW), Vecd(-Ls - BW, DH + BW),
            Vecd(DL1 - BW, DH + BW), Vecd(DL1 - BW, 2.0 * DH + BW),
            Vecd(DL + BW, 2.0 * DH + BW), Vecd(DL + BW, yb - BW),
            Vecd(DL1 - BW, yb - BW), Vecd(DL1 - BW, -BW),
            Vecd(-Ls - BW, -BW)};
}

std::vector<Vecd> innerWallShape(const CaseParams &p)
{
    const Real Ls = p.DL_sponge, DH = p.DH, DL = p.DL, DL1 = p.DL1, BW = p.BW;
    const Real yb = -p.branch_ratio * DH;
    return {Vecd(-Ls - BW, 0.0), Vecd(-Ls - BW, DH), Vecd(DL1, DH),
            Vecd(DL1, 2.0 * DH + BW), Vecd(DL, 2.0 * DH + BW), Vecd(DL, yb - BW),
            Vecd(DL1, yb - BW), Vecd(DL1, 0.0), Vecd(-Ls - BW, 0.0)};
}

//----------------------------------------------------------------------
//  Body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(waterBlockShape(g_params), GeometricOps::add);
        if (g_params.stenosis > 0.0)
        {
            multi_polygon_.addPolygon(stenosisPolygon(g_params, true), GeometricOps::sub);
            multi_polygon_.addPolygon(stenosisPolygon(g_params, false), GeometricOps::sub);
        }
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(outerWallShape(g_params), GeometricOps::add);
        multi_polygon_.addPolygon(innerWallShape(g_params), GeometricOps::sub);
        if (g_params.stenosis > 0.0)
        {
            multi_polygon_.addPolygon(stenosisPolygon(g_params, true), GeometricOps::add);
            multi_polygon_.addPolygon(stenosisPolygon(g_params, false), GeometricOps::add);
        }
    }
};

//----------------------------------------------------------------------
//  Pulsatile inflow velocity condition.
//
//  Phase 1 (t < t_ramp): cosine ramp 0 -> U_f, avoids an impulsive start.
//  Phase 2 (t >= t_ramp): u_ave = U_f * (1 + A * sin(omega * (t - t_ramp)))
//
//  Subtracting t_ramp means the pulsation always starts from its mean value with
//  a rising slope, so the phase reference is identical in every run - essential
//  for comparing phase lags across the study.
//
//  Spatial shape: Poiseuille parabola across the channel (Lecture 1).
//  Time variation: sinusoidal -> Womersley-type pulsatile flow (Lecture 4).
//
//  NOTE: `position` is in the inlet-box LOCAL frame (InflowVelocityCondition
//  applies transform_.shiftBaseStationToFrame first), so y = 0 is the channel
//  centreline and halfsize_[1] is the half-height R.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_, omega_, A_pulse_;
    OrientedBox &oriented_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(g_params.U_f), t_ref_(g_params.t_ramp),
          omega_(g_params.omega), A_pulse_(g_params.A),
          oriented_box_(boundary_condition.getOrientedBox()),
          halfsize_(oriented_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;

        Real u_ave;
        if (current_time < t_ref_)
        {
            u_ave = 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_));
        }
        else
        {
            u_ave = u_ref_ * (1.0 + A_pulse_ * sin(omega_ * (current_time - t_ref_)));
        }

        target_velocity[0] = 1.5 * u_ave *
                             SMAX(Real(0.0), Real(1.0) - position[1] * position[1] /
                                                             (halfsize_[1] * halfsize_[1]));
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//  Volume flow rate through an axis-aligned cross-section slab.
//
//  Two sums are reduced at once:
//      component 0 :  sum_i sign * (u_i . n) * Vol_i
//      component 1 :  sum_i Vol_i
//  and post-processing forms
//      Q = duct_width * component0 / component1
//
//  Dividing by the ACCUMULATED volume rather than by the nominal slab thickness is
//  deliberate. A slab of a few particle spacings contains a whole number of lattice
//  layers, so normalising by the nominal thickness carries a bias of up to one layer
//  (tens of percent). The volume-weighted mean velocity has no such bias: it is
//  exact for a uniform profile and second-order accurate otherwise.
//
//  A particle-level reduction is used rather than integrating an interpolated
//  profile: it needs no quadrature rule. Note that BodyPartByCell cannot be used
//  here - it loops over whole CELLS, so it would include particles outside the slab.
//----------------------------------------------------------------------
class SectionFlowRate : public LocalDynamicsReduce<ReduceSum<Vecd>>
{
  public:
    /**
     * @param flow_axis      0 for an x-normal section, 1 for a y-normal section
     * @param station        coordinate of the slab centre along flow_axis
     * @param half_thickness half slab thickness (1.5 particle spacings works well)
     * @param span_lo/hi     bounds along the transverse axis
     * @param sign           +1/-1 so that flow LEAVING the domain is positive
     */
    SectionFlowRate(SPHBody &body, const std::string &name, int flow_axis, Real station,
                    Real half_thickness, Real span_lo, Real span_hi, Real sign)
        : LocalDynamicsReduce<ReduceSum<Vecd>>(body),
          flow_axis_(flow_axis), transverse_axis_(1 - flow_axis),
          station_(station), half_thickness_(half_thickness),
          span_lo_(span_lo), span_hi_(span_hi), sign_(sign),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "FluxAndVolume_" + name;
    }
    virtual ~SectionFlowRate() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        if (std::abs(pos_[index_i][flow_axis_] - station_) > half_thickness_)
            return Vecd::Zero();
        const Real t = pos_[index_i][transverse_axis_];
        if (t < span_lo_ || t > span_hi_)
            return Vecd::Zero();
        return Vecd(sign_ * vel_[index_i][flow_axis_] * Vol_[index_i], Vol_[index_i]);
    }

  protected:
    int flow_axis_, transverse_axis_;
    Real station_, half_thickness_, span_lo_, span_hi_, sign_;
    Vecd *pos_, *vel_;
    Real *Vol_;
};

//----------------------------------------------------------------------
//  max |rho - rho0| / rho0 - the weak-compressibility validity check.
//  Should stay below ~1 % for a WCSPH result to be trustworthy.
//----------------------------------------------------------------------
class MaxDensityDeviation : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit MaxDensityDeviation(SPHBody &body)
        : LocalDynamicsReduce<ReduceMax>(body),
          rho0_(g_params.rho0_f),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          indicator_(particles_->getVariableDataByName<int>("Indicator"))
    {
        quantity_name_ = "MaxDensityDeviation";
    }
    virtual ~MaxDensityDeviation() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        // Skip free-surface / near-boundary particles (Indicator != 0): their
        // kernel support is truncated, so DensitySummationFreeStream reports a
        // deficit there that is a discretisation artefact, not compressibility.
        if (indicator_[index_i] != 0)
            return 0.0;
        return std::abs(rho_[index_i] - rho0_) / rho0_;
    }

  protected:
    Real rho0_;
    Real *rho_;
    int *indicator_;
};

//----------------------------------------------------------------------
//  Maximum speed among BULK particles only.
//
//  The library's MaximumSpeed reduces over every particle, including those in the
//  inflow buffer and in the outflow disposer bands, where the velocity is either
//  prescribed or belongs to a particle about to be deleted. Those are not part of
//  the solution, so including them overstates the Mach number that the weak-
//  compressibility assumption actually has to satisfy. Both are recorded so the
//  report can show the difference.
//----------------------------------------------------------------------
class MaxBulkSpeed : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit MaxBulkSpeed(SPHBody &body)
        : LocalDynamicsReduce<ReduceMax>(body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          indicator_(particles_->getVariableDataByName<int>("Indicator")),
          x_buffer_end_(0.0), y_up_(2.0 * g_params.DH - g_params.BW),
          y_lo_(-g_params.branch_ratio * g_params.DH + g_params.BW)
    {
        quantity_name_ = "MaxBulkSpeed";
    }
    virtual ~MaxBulkSpeed() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i] != 0)
            return 0.0;                                  // free-surface particle
        if (pos_[index_i][0] < x_buffer_end_)
            return 0.0;                                  // inside the inflow buffer
        if (pos_[index_i][1] > y_up_ || pos_[index_i][1] < y_lo_)
            return 0.0;                                  // inside a disposer band
        return vel_[index_i].norm();
    }

  protected:
    Vecd *pos_, *vel_;
    int *indicator_;
    Real x_buffer_end_, y_up_, y_lo_;
};

//----------------------------------------------------------------------
//  Observer point construction.
//----------------------------------------------------------------------
namespace probes
{
/** One straight (or stenosed) wall that WSS is sampled along. */
struct Wall
{
    std::string name;
    int normal_axis;  /**< 0 if the wall is x = const, 1 if y = const */
    Real station;     /**< the constant coordinate of the nominal wall */
    Real inward_sign; /**< +1 if the fluid lies at larger coordinate, else -1 */
    Real span_lo, span_hi;
    bool stenosed; /**< follow the stenosis contour along this wall */
};

inline std::vector<Wall> wallList(const CaseParams &p)
{
    const Real DH = p.DH, DL = p.DL, DL1 = p.DL1;
    const Real yb = -p.branch_ratio * DH;
    const Real d = disposerDepth(p) + 0.1; // keep clear of the outflow disposers
    const bool sten = p.stenosis > 0.0;

    return {
        // inlet channel, starting downstream of the inflow buffer (x > 0) where the
        // velocity is no longer prescribed and WSS is therefore meaningful
        {"inlet_bottom", 1, 0.0, +1.0, 0.2, DL1 - 0.2, false},
        {"inlet_top", 1, DH, -1.0, 0.2, DL1 - 0.2, false},
        // the wall the inlet jet impinges on: carries the stagnation point and the
        // outer wall of both branches
        {"divider", 0, DL, -1.0, yb + d, 2.0 * DH - d, sten},
        // inner (flow-divider side) walls of the two branches
        {"upper_inner", 0, DL1, +1.0, DH + 0.1, 2.0 * DH - d, sten},
        {"lower_inner", 0, DL1, +1.0, yb + d, -0.1, false},
    };
}

/** Flattened wall-probe cloud: for each wall, for each offset, for each station. */
struct ProbeCloud
{
    StdVec<Vecd> points;
    std::vector<int> wall_id, offset_id;
    std::vector<Real> arclength;
    std::vector<Real> nx, ny; /**< inward unit normal at the wall foot point */
    std::vector<Real> offsets;
};

inline ProbeCloud buildWallProbes(const CaseParams &p, int n_per_wall = 40)
{
    ProbeCloud c;
    c.offsets = {1.0 * p.dp, 2.0 * p.dp};
    const auto walls = wallList(p);

    for (size_t io = 0; io < c.offsets.size(); ++io)
    {
        const Real d = c.offsets[io];
        for (size_t iw = 0; iw < walls.size(); ++iw)
        {
            const Wall &w = walls[iw];
            for (int i = 0; i < n_per_wall; ++i)
            {
                const Real s = w.span_lo + (w.span_hi - w.span_lo) *
                                               Real(i) / Real(n_per_wall - 1);
                Real wall_coord = w.station;
                Vecd n_in = Vecd::Zero();
                n_in[w.normal_axis] = w.inward_sign;

                // On a stenosed wall the surface is displaced into the branch and the
                // normal tilts with the contour slope.
                if (w.stenosed && w.normal_axis == 0)
                {
                    const Real f = stenosisIntrusion(p, s);
                    const Real fp = stenosisIntrusionSlope(p, s);
                    wall_coord = w.station + w.inward_sign * f;
                    n_in = Vecd(w.inward_sign, -fp) / std::sqrt(1.0 + fp * fp);
                }

                Vecd foot = Vecd::Zero();
                foot[w.normal_axis] = wall_coord;
                foot[1 - w.normal_axis] = s;

                c.points.emplace_back(foot + d * n_in);
                c.wall_id.push_back(int(iw));
                c.offset_id.push_back(int(io));
                c.arclength.push_back(s);
                c.nx.push_back(n_in[0]);
                c.ny.push_back(n_in[1]);
            }
        }
    }
    return c;
}

/** Cross-section velocity profiles: inlet, upper branch, lower branch. */
struct ProfileCloud
{
    StdVec<Vecd> points;
    std::vector<std::string> names;
    std::vector<int> begin, count;
};

inline ProfileCloud buildProfiles(const CaseParams &p, int n = 21)
{
    ProfileCloud c;
    // Span 5 %..95 % of the duct so that every point stays inside the fluid at every
    // resolution; the fractions are resolution-independent so profiles from different
    // dp are sampled at identical stations (needed for the convergence study).
    auto add = [&](const std::string &name, int axis, Real station, Real lo, Real hi) {
        c.names.push_back(name);
        c.begin.push_back(int(c.points.size()));
        for (int i = 0; i < n; ++i)
        {
            const Real t = lo + (hi - lo) * Real(i) / Real(n - 1);
            Vecd pt = Vecd::Zero();
            pt[axis] = station;
            pt[1 - axis] = t;
            c.points.emplace_back(pt);
        }
        c.count.push_back(n);
    };

    const Real DH = p.DH, DL = p.DL, DL1 = p.DL1;
    add("inlet", 0, 0.5 * DL1, 0.05 * DH, 0.95 * DH);
    add("upper", 1, upperStation(p), DL1 + 0.05 * (DL - DL1), DL - 0.05 * (DL - DL1));
    add("lower", 1, lowerStation(p), DL1 + 0.05 * (DL - DL1), DL - 0.05 * (DL - DL1));
    return c;
}
} // namespace probes

//----------------------------------------------------------------------
//  Tiny JSON helpers (no external dependency).
//----------------------------------------------------------------------
namespace json
{
template <typename T>
std::string arr(const std::vector<T> &v)
{
    std::ostringstream o;
    o << std::setprecision(10) << "[";
    for (size_t i = 0; i < v.size(); ++i)
        o << (i ? ", " : "") << v[i];
    o << "]";
    return o.str();
}

inline std::string strArr(const std::vector<std::string> &v)
{
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < v.size(); ++i)
        o << (i ? ", " : "") << "\"" << v[i] << "\"";
    o << "]";
    return o.str();
}

inline std::string pointArr(const StdVec<Vecd> &v)
{
    std::ostringstream o;
    o << std::setprecision(10) << "[";
    for (size_t i = 0; i < v.size(); ++i)
        o << (i ? ", " : "") << "[" << v[i][0] << ", " << v[i][1] << "]";
    o << "]";
    return o.str();
}
} // namespace json

//----------------------------------------------------------------------
//  Main program.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //  Parse our own case options and hand SPHinXsys a filtered argv.
    //----------------------------------------------------------------------
    std::vector<char *> filtered_argv;
    TPipeStudy::parseCaseOptions(ac, av, g_params, filtered_argv);
    const CaseParams &p = g_params;
    p.print();

    //----------------------------------------------------------------------
    //  SPH system. The output folder is switched immediately so that every run of
    //  the sweep keeps its own results.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vec2d(-p.DL_sponge - p.BW, -p.branch_ratio * p.DH - p.BW),
        Vec2d(p.DL + p.BW, 2.0 * p.DH + p.BW));
    SPHSystem sph_system(system_domain_bounds, p.dp);
    IO::getEnvironment().resetOutputFolder("./output_" + p.case_name);
    sph_system.handleCommandlineOptions(int(filtered_argv.size()), filtered_argv.data());

    //----------------------------------------------------------------------
    //  Fluid and wall bodies.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(p.rho0_f, p.c_f);
    water_block.addMaterialProperty<Viscosity>(p.mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //  Observer bodies.
    //
    //  centreline observers - one point each, sampled over the WHOLE run so that
    //      cycle-to-cycle periodicity can be assessed (velocity and pressure).
    //  profile observer     - 3 x 21 cross-section points, analysis window only.
    //  wall probe observer  - near-wall probes at two offsets on five walls,
    //      analysis window only; carries VelocityGradient for WSS.
    //----------------------------------------------------------------------
    ObserverBody inlet_cl_observer(sph_system, "InletCentrelineObserver");
    inlet_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * p.DL1, 0.5 * p.DH)});

    ObserverBody upper_cl_observer(sph_system, "UpperOutletObserver");
    upper_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * (p.DL1 + p.DL), upperStation(p))});

    ObserverBody lower_cl_observer(sph_system, "LowerOutletObserver");
    lower_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * (p.DL1 + p.DL), lowerStation(p))});

    const probes::ProfileCloud profile_cloud = probes::buildProfiles(p);
    ObserverBody profile_observer(sph_system, "ProfileObserver");
    profile_observer.generateParticles<ObserverParticles>(profile_cloud.points);

    const probes::ProbeCloud wall_cloud = probes::buildWallProbes(p);
    ObserverBody wall_probe_observer(sph_system, "WallProbeObserver");
    wall_probe_observer.generateParticles<ObserverParticles>(wall_cloud.points);

    std::cout << "observer points: profile = " << profile_cloud.points.size()
              << ", wall probes = " << wall_cloud.points.size() << "\n";

    //----------------------------------------------------------------------
    //  Body relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    ContactRelation inlet_cl_contact(inlet_cl_observer, {&water_block});
    ContactRelation upper_cl_contact(upper_cl_observer, {&water_block});
    ContactRelation lower_cl_contact(lower_cl_observer, {&water_block});
    ContactRelation profile_contact(profile_observer, {&water_block});
    ContactRelation wall_probe_contact(wall_probe_observer, {&water_block});

    //----------------------------------------------------------------------
    //  Physics methods (identical to the baseline apart from the velocity gradient).
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann>
        pressure_relaxation(water_block_inner, water_wall_contact);
    // The Riemann variant of the continuity step, where the baseline T_shaped_pipe
    // uses NoRiemann. The baseline is validated only over a short run (t ~ 33) with
    // a steady inflow; this study runs to t ~ 150 with a pulsating inflow, i.e. an
    // order of magnitude more exposure to the sharp re-entrant corners at the
    // junction. With NoRiemann the run is flaky - repeated identical runs either
    // finish cleanly or collapse to dt ~ 5e-5 and exhaust the particle buffer,
    // depending only on thread scheduling. The Riemann flux adds the dissipation
    // needed to make the long pulsatile runs reproducible.
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann>
        density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>
        viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex>
        update_density_by_summation(water_block_inner, water_wall_contact);

    // Velocity gradient tensor -> wall shear stress. DistanceFromWall must exist and
    // be up to date before the gradient's wall contribution is evaluated.
    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<NoKernelCorrection>>
        update_velocity_gradient(water_block_inner, water_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep>
        get_fluid_advection_time_step_size(water_block, p.U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep>
        get_fluid_time_step_size(water_block);

    //----------------------------------------------------------------------
    //  Inlet emitter, outlet disposers, particle sorting.
    //----------------------------------------------------------------------
    Vec2d emitter_halfsize = Vec2d(0.5 * p.BW, 0.5 * p.DH);
    Vec2d emitter_translation = Vec2d(-p.DL_sponge, 0.0) + emitter_halfsize;
    OrientedBoxByParticle emitter(water_block,
                                  OrientedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection>
        emitter_inflow_injection(emitter, inlet_particle_buffer);

    Vec2d inlet_buffer_halfsize = Vec2d(0.5 * p.DL_sponge, 0.5 * p.DH);
    Vec2d inlet_buffer_translation = Vec2d(-p.DL_sponge, 0.0) + inlet_buffer_halfsize;
    OrientedBoxByCell inlet_flow_buffer(water_block,
                                        OrientedBox(xAxis, Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>>
        inflow_condition(inlet_flow_buffer);

    Vec2d disposer_up_halfsize = Vec2d(0.3 * p.DH, 0.5 * p.BW);
    Vec2d disposer_up_translation = Vec2d(p.DL + 0.05 * p.DH, 2.0 * p.DH) - disposer_up_halfsize;
    OrientedBoxByCell disposer_up(water_block,
                                  OrientedBox(yAxis, Transform(Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up);

    Vec2d disposer_down_halfsize = disposer_up_halfsize;
    Vec2d disposer_down_translation =
        Vec2d(p.DL1 - 0.05 * p.DH, -p.branch_ratio * p.DH) + disposer_down_halfsize;
    OrientedBoxByCell disposer_down(water_block,
                                    OrientedBox(yAxis, Transform(Rotation2d(Pi), Vec2d(disposer_down_translation)), disposer_down_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_down_outflow_deletion(disposer_down);

    ParticleSorting particle_sorting(water_block);

    //----------------------------------------------------------------------
    //  I/O.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure");
    write_body_states.addToWrite<Real>(water_block, "Density");
    write_body_states.addToWrite<int>(water_block, "Indicator");
    write_body_states.addToWrite<Vecd>(water_block, "Velocity");
    write_body_states.addToWrite<Matd>(water_block, "VelocityGradient");

    // Global diagnostics, whole run: periodicity indicator + WCSPH validity checks.
    ReducedQuantityRecording<TotalKineticEnergy> write_kinetic_energy(water_block);
    ReducedQuantityRecording<MaximumSpeed> write_max_speed(water_block);
    ReducedQuantityRecording<MaxBulkSpeed> write_max_bulk_speed(water_block);
    ReducedQuantityRecording<MaxDensityDeviation> write_max_density_deviation(water_block);

    // Volume flow rate through three cross-sections; outflow positive.
    // Q_inlet = Q_upper + Q_lower is mass conservation (Lecture 1).
    const Real span_pad = 2.0 * p.BW;
    const Real slab_half = 1.5 * p.dp;
    ReducedQuantityRecording<SectionFlowRate> write_Q_inlet(
        water_block, std::string("inlet"), 0, Real(0.5) * p.DL1, slab_half,
        -span_pad, p.DH + span_pad, Real(1.0));
    ReducedQuantityRecording<SectionFlowRate> write_Q_upper(
        water_block, std::string("upper"), 1, upperStation(p), slab_half,
        p.DL1 - span_pad, p.DL + span_pad, Real(1.0));
    ReducedQuantityRecording<SectionFlowRate> write_Q_lower(
        water_block, std::string("lower"), 1, lowerStation(p), slab_half,
        p.DL1 - span_pad, p.DL + span_pad, Real(-1.0));

    // Centreline point histories, whole run.
    ObservedQuantityRecording<Vecd> write_inlet_velocity("Velocity", inlet_cl_contact);
    ObservedQuantityRecording<Vecd> write_upper_velocity("Velocity", upper_cl_contact);
    ObservedQuantityRecording<Vecd> write_lower_velocity("Velocity", lower_cl_contact);
    ObservedQuantityRecording<Real> write_inlet_pressure("Pressure", inlet_cl_contact);
    ObservedQuantityRecording<Real> write_upper_pressure("Pressure", upper_cl_contact);
    ObservedQuantityRecording<Real> write_lower_pressure("Pressure", lower_cl_contact);

    // Cross-section profiles and wall probes, analysis window only.
    ObservedQuantityRecording<Vecd> write_profile_velocity("Velocity", profile_contact);
    ObservedQuantityRecording<Vecd> write_probe_velocity("Velocity", wall_probe_contact);
    ObservedQuantityRecording<Real> write_probe_pressure("Pressure", wall_probe_contact);
    ObservedQuantityRecording<Matd> write_probe_velocity_gradient("VelocityGradient", wall_probe_contact);

    //----------------------------------------------------------------------
    //  Initialise system.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    distance_to_wall.exec();
    update_velocity_gradient.exec();

    const int n_fluid_particles = int(water_block.getBaseParticles().TotalRealParticles());
    std::cout << "fluid particles = " << n_fluid_particles
              << ", wall particles = " << wall_boundary.getBaseParticles().TotalRealParticles()
              << "\n\n";

    //----------------------------------------------------------------------
    //  Time-stepping control. All cadences are derived from the pulse period so
    //  that every case is sampled at the same number of points per cycle.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    const int screen_output_interval = 200;

    Real next_sample = 0.0;   // centreline + flow rate, whole run
    Real next_ke = 0.0;       // global diagnostics, whole run
    Real next_window = p.analysis_start; // profiles + wall probes, analysis window
    Real next_vtp = p.analysis_start;    // VTP snapshots, analysis window
    const Real window_interval = p.T / 50.0;
    Real dt = 0.0;

    auto writeWholeRunSamples = [&]() {
        inlet_cl_contact.updateConfiguration();
        upper_cl_contact.updateConfiguration();
        lower_cl_contact.updateConfiguration();
        write_inlet_velocity.writeToFile(number_of_iterations);
        write_upper_velocity.writeToFile(number_of_iterations);
        write_lower_velocity.writeToFile(number_of_iterations);
        write_inlet_pressure.writeToFile(number_of_iterations);
        write_upper_pressure.writeToFile(number_of_iterations);
        write_lower_pressure.writeToFile(number_of_iterations);
        write_Q_inlet.writeToFile(number_of_iterations);
        write_Q_upper.writeToFile(number_of_iterations);
        write_Q_lower.writeToFile(number_of_iterations);
    };

    auto writeWindowSamples = [&]() {
        profile_contact.updateConfiguration();
        wall_probe_contact.updateConfiguration();
        write_profile_velocity.writeToFile(number_of_iterations);
        write_probe_velocity.writeToFile(number_of_iterations);
        write_probe_pressure.writeToFile(number_of_iterations);
        write_probe_velocity_gradient.writeToFile(number_of_iterations);
    };

    TickCount t1 = TickCount::now();
    TimeInterval output_interval_accum;

    writeWholeRunSamples();
    next_sample += p.sample_interval;
    write_kinetic_energy.writeToFile(number_of_iterations);
    write_max_speed.writeToFile(number_of_iterations);
    write_max_bulk_speed.writeToFile(number_of_iterations);
    write_max_density_deviation.writeToFile(number_of_iterations);
    next_ke += p.ke_interval;
    write_body_states.writeToFile(); // t = 0 snapshot for the geometry check

    //----------------------------------------------------------------------
    //  Main time loop. One pass = one advection (Dt) step; all output is driven by
    //  the cadence counters above rather than by the loop structure, which is what
    //  decouples observer sampling from VTP writing.
    //----------------------------------------------------------------------
    while (physical_time < p.end_time)
    {
        Real Dt = get_fluid_advection_time_step_size.exec();
        inlet_outlet_surface_particle_indicator.exec();
        update_density_by_summation.exec();
        viscous_force.exec();
        transport_velocity_correction.exec();

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            dt = SMIN(p.cfl_scale * get_fluid_time_step_size.exec(), Dt - relaxation_time);
            pressure_relaxation.exec(dt);
            inflow_condition.exec();
            density_relaxation.exec(dt);

            relaxation_time += dt;
            physical_time += dt;
        }

        if (number_of_iterations % screen_output_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6)
                      << "N=" << number_of_iterations
                      << "  t=" << physical_time
                      << "  t/T=" << (physical_time - p.t_ramp) / p.T
                      << "  Dt=" << Dt << "  dt=" << dt << "\n";
        }
        number_of_iterations++;

        emitter_inflow_injection.exec();
        disposer_up_outflow_deletion.exec();
        disposer_down_outflow_deletion.exec();

        if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            particle_sorting.exec();

        water_block.updateCellLinkedList();
        water_block_complex.updateConfiguration();

        // The velocity gradient is a post-processing quantity here, so it only has to
        // be current when the wall probes are about to be sampled.
        const bool in_window = physical_time >= p.analysis_start;
        if (in_window && physical_time >= next_window)
        {
            distance_to_wall.exec();
            update_velocity_gradient.exec();
        }

        TickCount t2 = TickCount::now();

        if (physical_time >= next_sample)
        {
            writeWholeRunSamples();
            while (next_sample <= physical_time)
                next_sample += p.sample_interval;
        }

        if (physical_time >= next_ke)
        {
            write_kinetic_energy.writeToFile(number_of_iterations);
            write_max_speed.writeToFile(number_of_iterations);
            write_max_bulk_speed.writeToFile(number_of_iterations);
            write_max_density_deviation.writeToFile(number_of_iterations);
            while (next_ke <= physical_time)
                next_ke += p.ke_interval;
        }

        if (in_window && physical_time >= next_window)
        {
            writeWindowSamples();
            while (next_window <= physical_time)
                next_window += window_interval;
        }

        if (in_window && physical_time >= next_vtp)
        {
            write_body_states.writeToFile();
            while (next_vtp <= physical_time)
                next_vtp += p.vtp_interval;
        }

        output_interval_accum += TickCount::now() - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - output_interval_accum;
    std::cout << "\nTotal wall time for computation: " << tt.seconds() << " s\n";

    //----------------------------------------------------------------------
    //  Case metadata for the Python analysis: every derived parameter plus the
    //  exact observer point tables, so nothing has to be re-derived downstream.
    //----------------------------------------------------------------------
    std::ostringstream obs;
    obs << "  \"observers\": {\n";
    obs << "    \"profile\": {\n";
    obs << "      \"names\": " << json::strArr(profile_cloud.names) << ",\n";
    obs << "      \"begin\": " << json::arr(profile_cloud.begin) << ",\n";
    obs << "      \"count\": " << json::arr(profile_cloud.count) << ",\n";
    obs << "      \"points\": " << json::pointArr(profile_cloud.points) << "\n";
    obs << "    },\n";

    const auto walls = probes::wallList(p);
    std::vector<std::string> wall_names;
    for (const auto &w : walls)
        wall_names.push_back(w.name);
    obs << "    \"wall_probes\": {\n";
    obs << "      \"wall_names\": " << json::strArr(wall_names) << ",\n";
    obs << "      \"offsets\": " << json::arr(wall_cloud.offsets) << ",\n";
    obs << "      \"wall_id\": " << json::arr(wall_cloud.wall_id) << ",\n";
    obs << "      \"offset_id\": " << json::arr(wall_cloud.offset_id) << ",\n";
    obs << "      \"arclength\": " << json::arr(wall_cloud.arclength) << ",\n";
    obs << "      \"normal_x\": " << json::arr(wall_cloud.nx) << ",\n";
    obs << "      \"normal_y\": " << json::arr(wall_cloud.ny) << ",\n";
    obs << "      \"points\": " << json::pointArr(wall_cloud.points) << "\n";
    obs << "    },\n";

    obs << "    \"centrelines\": {\n";
    obs << "      \"names\": [\"inlet\", \"upper\", \"lower\"],\n";
    obs << "      \"points\": [[" << 0.5 * p.DL1 << ", " << 0.5 * p.DH << "], ["
        << 0.5 * (p.DL1 + p.DL) << ", " << upperStation(p) << "], ["
        << 0.5 * (p.DL1 + p.DL) << ", " << lowerStation(p) << "]]\n";
    obs << "    },\n";

    obs << "    \"sections\": {\n";
    obs << "      \"names\": [\"inlet\", \"upper\", \"lower\"],\n";
    obs << "      \"axis\": [0, 1, 1],\n";
    obs << "      \"station\": [" << 0.5 * p.DL1 << ", " << upperStation(p) << ", "
        << lowerStation(p) << "],\n";
    obs << "      \"sign\": [1, 1, -1],\n";
    obs << "      \"half_thickness\": " << slab_half << ",\n";
    obs << "      \"width\": [" << p.DH << ", " << p.DL - p.DL1 << ", " << p.DL - p.DL1 << "]\n";
    obs << "    }\n";
    obs << "  },\n";

    p.writeJson("./output_" + p.case_name + "/case_params.json",
                n_fluid_particles, tt.seconds(), obs.str());

    std::cout << "wrote ./output_" << p.case_name << "/case_params.json\n";
    return 0;
}
