/**
 * @file    T_pipe_pulsatile_poiseuille_study.cpp
 * @brief   Pulsatile Poiseuille flow in a 2D T-shaped bifurcation — parameter study.
 *
 * @details Extends the T_shaped_pipe benchmark with:
 *            1. Sinusoidal pulsatile inflow (amplitude A_pulse around mean U_f).
 *            2. Parabolic Poiseuille spatial profile at all times.
 *            3. Observer points for post-processing velocity time series.
 *            4. Cross-section profile observer for Poiseuille parabola comparison.
 *
 *          Lecture connections:
 *            - Poiseuille parabolic profile:      Lecture 1
 *            - Mass conservation Q1 = Q2 + Q3:   Lecture 1
 *            - Womersley number alpha:            Lecture 4
 *            - Pulsatility index PI:              Lecture 3
 *
 *          Parameter study matrix (edit Re and omega_pulse below, adjust end_time):
 *          -----------------------------------------------------------------------
 *          Run | Re  | omega_pulse | alpha (Womersley) | Period T | end_time
 *           1  |  50 |  0.50       |  ~2.7             |  12.6    |  65.0
 *           2  | 100 |  0.50       |  ~3.9             |  12.6    |  65.0
 *           3  | 100 |  1.00       |  ~8.7             |   6.3    |  34.0
 *           4  | 200 |  1.00       |  ~12.2            |   6.3    |  34.0
 *           5  | 100 |  2.00       |  ~12.2            |   3.1    |  18.0
 *          -----------------------------------------------------------------------
 *          alpha = (DH/2) * sqrt(omega_pulse * rho0_f / mu_f)
 *          Low alpha  (~1-3): quasi-steady, parabolic profile at every instant.
 *          High alpha (~10+): oscillatory boundary layer, nearly flat core profile.
 *
 * @author  Based on T_shaped_pipe by Xiangyu Hu and Shuoguo Zhang.
 *          Extended for pulsatile parameter study.
 */

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//  Geometry parameters.
//----------------------------------------------------------------------
Real DL             = 5.0;                          /**< Domain reference length. */
Real DH             = 3.0;                          /**< Main channel height. */
Real DL1            = 0.7 * DL;                     /**< Length of main inlet channel. */
Real global_resolution = 0.15;                     /**< Particle spacing. */
Real BW             = global_resolution * 4;        /**< Wall boundary width. */
Real DL_sponge      = global_resolution * 20;       /**< Emitter buffer length. */

//----------------------------------------------------------------------
//  Fluid material parameters.
//  To run a different Re: change Re here; mu_f updates automatically.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                          /**< Reference density. */
Real U_f    = 1.0;                          /**< Mean inflow velocity. */
Real c_f    = 10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * (DL - DL1)));
Real Re     = 100.0;                        /**< Reynolds number — edit for parameter study. */
Real mu_f   = rho0_f * U_f * DH / Re;      /**< Dynamic viscosity. */

//----------------------------------------------------------------------
//  Pulsatile inflow parameters.
//  omega_pulse: angular frequency of oscillation [rad/s].
//  A_pulse:     amplitude as a fraction of U_f (0.5 = ±50% around mean).
//  To run a different Womersley number: change omega_pulse here.
//----------------------------------------------------------------------
Real omega_pulse = 1.0;     /**< Pulsation angular frequency [rad/s] — edit for parameter study. */
Real A_pulse     = 0.5;     /**< Pulsation amplitude (fraction of mean velocity). */

//----------------------------------------------------------------------
//  Geometry polygons (unchanged from baseline T_shaped_pipe).
//----------------------------------------------------------------------
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH),
    Vecd(DL, 2.0 * DH),    Vecd(DL, -DH),        Vecd(DL1, -DH), Vecd(DL1, 0.0),
    Vecd(-DL_sponge, 0.0)};

std::vector<Vecd> outer_wall_shape{
    Vecd(-DL_sponge - BW, -BW),         Vecd(-DL_sponge - BW, DH + BW),
    Vecd(DL1 - BW, DH + BW),            Vecd(DL1 - BW, 2.0 * DH + BW),
    Vecd(DL + BW,  2.0 * DH + BW),      Vecd(DL + BW,  -DH - BW),
    Vecd(DL1 - BW, -DH - BW),           Vecd(DL1 - BW, -BW),
    Vecd(-DL_sponge - BW, -BW)};

std::vector<Vecd> inner_wall_shape{
    Vecd(-DL_sponge - BW, 0.0),         Vecd(-DL_sponge - BW, DH),
    Vecd(DL1, DH),                       Vecd(DL1, 2.0 * DH + BW),
    Vecd(DL,  2.0 * DH + BW),           Vecd(DL,  -DH - BW),
    Vecd(DL1, -DH - BW),                 Vecd(DL1, 0.0),
    Vecd(-DL_sponge - BW, 0.0)};

//----------------------------------------------------------------------
//  Body shape classes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(water_block_shape, GeometricOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(outer_wall_shape, GeometricOps::add);
        multi_polygon_.addPolygon(inner_wall_shape, GeometricOps::sub);
    }
};

//----------------------------------------------------------------------
//  Pulsatile inflow velocity condition.
//
//  Two phases:
//    Phase 1 (t < t_ref): cosine ramp-up 0 → U_f  (avoids impulsive start).
//    Phase 2 (t >= t_ref): U_f * (1 + A_pulse * sin(omega * t))
//
//  Spatial shape: Poiseuille parabola across channel height (Lecture 1).
//  Time variation: sinusoidal → Womersley-type pulsatile flow (Lecture 4).
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_, omega_, A_pulse_;
    OrientedBox &oriented_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          omega_(omega_pulse), A_pulse_(A_pulse),
          oriented_box_(boundary_condition.getOrientedBox()),
          halfsize_(oriented_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;

        Real u_ave;
        if (current_time < t_ref_)
        {
            // smooth ramp-up — avoids numerical shock at startup
            u_ave = 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_));
        }
        else
        {
            // pulsatile mean flow with sinusoidal oscillation
            u_ave = u_ref_ * (1.0 + A_pulse_ * sin(omega_ * current_time));
        }

        // Poiseuille parabolic profile: u(y) = 1.5 * u_ave * (1 - y^2/R^2)
        // position[1] is measured from the channel centreline (halfsize_[1] = R)
        target_velocity[0] = 1.5 * u_ave *
            SMAX(0.0, 1.0 - position[1] * position[1] / (halfsize_[1] * halfsize_[1]));
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//  Helper: 21-point cross-section profile in the inlet channel.
//  Points span y = 0 to DH at x = DL1/2 (midway along inlet channel).
//  Used to compare simulated profile against analytical Poiseuille parabola.
//----------------------------------------------------------------------
StdVec<Vecd> createProfileObserverPoints()
{
    StdVec<Vecd> points;
    Real x = 0.5 * DL1;                  // midway along the main inlet channel
    int  n = 21;
    for (int i = 0; i < n; ++i)
    {
        Real y = DH * Real(i) / Real(n - 1);  // y from 0 to DH
        points.emplace_back(Vecd(x, y));
    }
    return points;
}

//----------------------------------------------------------------------
//  Main program.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //  Print run parameters and derived quantities.
    //----------------------------------------------------------------------
    Real nu      = mu_f / rho0_f;
    Real R       = 0.5 * DH;
    Real alpha_w = R * std::sqrt(omega_pulse / nu);
    Real T_pulse = 2.0 * Pi / omega_pulse;

    std::cout << "\n===== Pulsatile T-pipe parameter study run =====\n"
              << "  Re               = " << Re          << "\n"
              << "  omega_pulse      = " << omega_pulse << " rad/s\n"
              << "  A_pulse          = " << A_pulse     << "\n"
              << "  Womersley alpha  = " << alpha_w     << "\n"
              << "  Pulse period  T  = " << T_pulse     << "\n"
              << "  Kinematic visc   = " << nu          << "\n"
              << "=================================================\n\n";

    //----------------------------------------------------------------------
    //  SPH system.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH), Vec2d(DL, 2.0 * DH));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);

    //----------------------------------------------------------------------
    //  Fluid and wall bodies.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //  Observer bodies.
    //
    //  Geometry reference:
    //    Main inlet channel:  x in [-DL_sponge, DL1], y in [0, DH]
    //    Upper outlet branch: x in [DL1, DL],         y in [DH, 2*DH]
    //    Lower outlet branch: x in [DL1, DL],         y in [-DH, 0]
    //
    //  inlet_cl_observer  — centreline of inlet channel (midway along x).
    //                       Measures incoming velocity waveform.
    //  upper_cl_observer  — centreline of upper outlet branch.
    //                       Combined with lower: verify Q1 = Q2 + Q3 (Lecture 1).
    //                       Compute pulsatility index PI = (Vmax-Vmin)/Vmean (Lecture 3).
    //  lower_cl_observer  — centreline of lower outlet branch.
    //  profile_observer   — 21 points across inlet channel height at x = DL1/2.
    //                       Compare simulated profile against Poiseuille parabola.
    //----------------------------------------------------------------------
    ObserverBody inlet_cl_observer(sph_system, "InletCentrelineObserver");
    inlet_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * DL1, 0.5 * DH)});

    ObserverBody upper_cl_observer(sph_system, "UpperOutletObserver");
    upper_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * (DL1 + DL), 1.5 * DH)});

    ObserverBody lower_cl_observer(sph_system, "LowerOutletObserver");
    lower_cl_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{Vecd(0.5 * (DL1 + DL), -0.5 * DH)});

    ObserverBody profile_observer(sph_system, "ProfileObserver");
    profile_observer.generateParticles<ObserverParticles>(createProfileObserverPoints());

    //----------------------------------------------------------------------
    //  Body relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    ContactRelation inlet_cl_contact(inlet_cl_observer,  {&water_block});
    ContactRelation upper_cl_contact(upper_cl_observer,  {&water_block});
    ContactRelation lower_cl_contact(lower_cl_observer,  {&water_block});
    ContactRelation profile_contact(profile_observer,    {&water_block});

    //----------------------------------------------------------------------
    //  Physics methods (unchanged from baseline).
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann>
        pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann>
        density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>
        viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex>
        update_density_by_summation(water_block_inner, water_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep>
        get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep>
        get_fluid_time_step_size(water_block);

    //----------------------------------------------------------------------
    //  Inlet emitter, outlet disposers, particle sorting (unchanged).
    //----------------------------------------------------------------------
    Vec2d emitter_halfsize       = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d emitter_translation    = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
    OrientedBoxByParticle emitter(water_block,
        OrientedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection>
        emitter_inflow_injection(emitter, inlet_particle_buffer);

    Vec2d inlet_buffer_halfsize      = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d inlet_buffer_translation   = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;
    OrientedBoxByCell inlet_flow_buffer(water_block,
        OrientedBox(xAxis, Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>>
        inflow_condition(inlet_flow_buffer);

    Vec2d disposer_up_halfsize      = Vec2d(0.3 * DH, 0.5 * BW);
    Vec2d disposer_up_translation   = Vec2d(DL + 0.05 * DH, 2.0 * DH) - disposer_up_halfsize;
    OrientedBoxByCell disposer_up(water_block,
        OrientedBox(yAxis, Transform(Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up);

    Vec2d disposer_down_halfsize    = disposer_up_halfsize;
    Vec2d disposer_down_translation = Vec2d(DL1 - 0.05 * DH, -DH) + disposer_down_halfsize;
    OrientedBoxByCell disposer_down(water_block,
        OrientedBox(yAxis, Transform(Rotation2d(Pi), Vec2d(disposer_down_translation)), disposer_down_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_down_outflow_deletion(disposer_down);

    ParticleSorting particle_sorting(water_block);

    //----------------------------------------------------------------------
    //  I/O: VTP body states + observer time series.
    //  Velocity is written for all particles (needed for ParaView visualisation).
    //  Each observer writes a separate dat file — load in Python/MATLAB to plot.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure");
    write_body_states.addToWrite<int>(water_block,  "Indicator");
    write_body_states.addToWrite<Vecd>(water_block, "Velocity");

    ReducedQuantityRecording<TotalKineticEnergy>
        write_water_kinetic_energy(water_block);

    ObservedQuantityRecording<Vecd> write_inlet_velocity("Velocity",  inlet_cl_contact);
    ObservedQuantityRecording<Vecd> write_upper_velocity("Velocity",  upper_cl_contact);
    ObservedQuantityRecording<Vecd> write_lower_velocity("Velocity",  lower_cl_contact);
    ObservedQuantityRecording<Vecd> write_profile_velocity("Velocity", profile_contact);

    //----------------------------------------------------------------------
    //  Initialise system.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();

    //----------------------------------------------------------------------
    //  Time-stepping control.
    //  end_time      = ramp-up time + 5 complete pulse cycles.
    //  output_interval = T_pulse / 20  (20 snapshots per cycle for smooth waveforms).
    //----------------------------------------------------------------------
    Real &physical_time     = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval  = 100;
    Real t_ramp          = 2.0;
    Real end_time        = t_ramp + 5.0 * T_pulse;
    Real output_interval = T_pulse / 20.0;
    Real dt              = 0.0;

    std::cout << "end_time       = " << end_time        << "\n"
              << "output_interval = " << output_interval << "\n\n";

    //----------------------------------------------------------------------
    //  First output (t = 0).
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    write_body_states.writeToFile();
    write_inlet_velocity.writeToFile(number_of_iterations);
    write_upper_velocity.writeToFile(number_of_iterations);
    write_lower_velocity.writeToFile(number_of_iterations);
    write_profile_velocity.writeToFile(number_of_iterations);

    //----------------------------------------------------------------------
    //  Main time loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time   += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "  Time=" << physical_time
                          << "  Dt=" << Dt << "  dt=" << dt << "\n";
                write_water_kinetic_energy.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            emitter_inflow_injection.exec();
            disposer_up_outflow_deletion.exec();
            disposer_down_outflow_deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                particle_sorting.exec();

            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();

        // Write VTP snapshot
        write_body_states.writeToFile();

        // Update observer contacts and write velocity time series
        inlet_cl_contact.updateConfiguration();
        upper_cl_contact.updateConfiguration();
        lower_cl_contact.updateConfiguration();
        profile_contact.updateConfiguration();

        write_inlet_velocity.writeToFile(number_of_iterations);
        write_upper_velocity.writeToFile(number_of_iterations);
        write_lower_velocity.writeToFile(number_of_iterations);
        write_profile_velocity.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "\nTotal wall time for computation: " << tt.seconds() << " s\n";

    return 0;
}
