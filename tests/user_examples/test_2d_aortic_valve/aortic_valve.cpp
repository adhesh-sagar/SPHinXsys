/**
 * @file aortic_valve.cpp
 * @brief 2D aortic valve example ported from the SPHinXsys tutorial.
 */

#include "aortic_valve.h"
#include "sphinxsys.h"

using namespace SPH;

int main(int ac, char *av[])
{
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5 * DH - BW),
                                      Vec2d(DL + BW, 0.5 * DH + R + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody inserted_body(sph_system, makeShared<InsertedBody>("InsertedBody"));
    inserted_body.defineAdaptationRatios(1.15, 1.0);
    inserted_body.defineBodyLevelSetShape().writeLevelSet();
    inserted_body.defineMatterMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? inserted_body.generateParticles<BaseParticles, Reload>(inserted_body.Name())
        : inserted_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody leaflet_observer(sph_system, "LeafletObserver");
    StdVec<Vecd> leaflet_observation_location = {0.5 * (BRT + BRB)};
    leaflet_observer.generateParticles<ObserverParticles>(leaflet_observation_location);
    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(createObservationPoints());

    if (sph_system.RunParticleRelaxation())
    {
        InnerRelation inserted_body_inner(inserted_body);

        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(inserted_body);
        RelaxationStepInner relaxation_step_inner(inserted_body_inner);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(inserted_body);
        ReloadParticleIO write_particle_reload_files(inserted_body);

        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }

    InnerRelation water_block_inner(water_block);
    InnerRelation inserted_body_inner(inserted_body);
    ContactRelation water_block_contact(water_block, {&inserted_body});
    ContactRelation inserted_body_contact(inserted_body, {&water_block});
    ContactRelation leaflet_observer_contact(leaflet_observer, {&inserted_body});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    SimpleDynamics<NormalDirectionFromBodyShape> inserted_body_normal_direction(inserted_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> inserted_body_corrected_configuration(inserted_body_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> inserted_body_stress_relaxation_first_half(inserted_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> inserted_body_stress_relaxation_second_half(inserted_body_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> inserted_body_computing_time_step_size(inserted_body);
    BodyRegionByParticle leaflet_base(inserted_body, makeShared<MultiPolygonShape>(createLeafletBaseShape()));
    SimpleDynamics<FixBodyPartConstraint> constrain_leaflet_base(leaflet_base);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    OrientedBoxByCell inflow_buffer(water_block, OrientedBox(xAxis, Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);

    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);

    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(inserted_body);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> inserted_body_update_normal(inserted_body);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(inserted_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(inserted_body_contact);

    ParticleSorting particle_sorting(water_block);

    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force_from_fluid(inserted_body, "ViscousForceFromFluid");
    ObservedQuantityRecording<Vecd> write_leaflet_tip_displacement("Position", leaflet_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);

    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    inserted_body_normal_direction.exec();
    inserted_body_corrected_configuration.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 4.0;
    Real output_interval = end_time / 100.0;

    TickCount t1 = TickCount::now();
    TimeInterval interval;

    write_real_body_states.writeToFile();
    write_leaflet_tip_displacement.writeToFile(number_of_iterations);
    write_fluid_velocity.writeToFile(number_of_iterations);

    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_correction.exec();

            viscous_force_from_fluid.exec();
            inserted_body_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                density_relaxation.exec(dt);

                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(inserted_body_computing_time_step_size.exec(), dt - dt_s_sum);
                    inserted_body_stress_relaxation_first_half.exec(dt_s);
                    constrain_leaflet_base.exec();
                    inserted_body_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                parabolic_inflow.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "	Time = " << physical_time
                          << "	Dt = " << Dt
                          << "	Dt / dt = " << inner_ite_dt
                          << "	dt / dt_s = " << inner_ite_dt_s << "\n";
                write_leaflet_tip_displacement.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();
            inserted_body.updateCellLinkedList();
            inserted_body_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
