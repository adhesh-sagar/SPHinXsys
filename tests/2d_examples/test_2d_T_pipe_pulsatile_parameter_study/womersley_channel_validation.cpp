/**
 * @file    womersley_channel_validation.cpp
 * @brief   Textbook Womersley benchmark: oscillatory flow in a plane channel.
 *
 * @details WHY THIS CASE EXISTS, SEPARATELY FROM THE T-PIPE
 *
 *          The T-pipe cannot validate the Womersley solution, and it is worth being
 *          explicit about why rather than quietly comparing anyway:
 *
 *            - its inlet velocity profile is PRESCRIBED as a parabola by the inflow
 *              buffer, so the profile is imposed rather than computed;
 *            - the inlet channel is only 3.5 long, while the viscous entrance length
 *              is ~0.04 Re D = 12 and the oscillatory development length is
 *              U/omega = 3.0. The flow therefore never develops away from the
 *              imposed parabola before it reaches the junction.
 *
 *          Measured on the T-pipe, the profile "peakedness" comes out at 1.41
 *          against a Womersley prediction of 1.21 - i.e. it is still essentially
 *          the imposed parabola (1.5), exactly as the length scales say it must be.
 *
 *          This case removes both problems: a straight channel, periodic in x, with
 *          NO velocity boundary condition at all. The flow is driven by an
 *          oscillatory body force, which is precisely the configuration Womersley's
 *          analytical solution describes:
 *
 *              rho du/dt = G cos(omega t) + mu d2u/dy2,   u(+-R) = 0.
 *
 *          Everything in the interior is then genuinely computed, so comparing the
 *          simulated profile against the analytical one at several phases is a real
 *          validation of the SPH scheme's unsteady viscous behaviour, and it is the
 *          figure that belongs in the report's validation section.
 *
 *          The case is defined by alpha alone (plus a resolution). Following the
 *          same convention as the parameter study, omega is DERIVED:
 *
 *              nu    = U_ref * (2R) / Re
 *              omega = nu * (alpha / R)^2
 *
 *          and the forcing amplitude G is chosen so that the peak cross-sectional
 *          mean velocity is U_ref, which keeps the Mach number under control at
 *          every alpha.
 *
 * @author  Written for the pulsatile T-pipe study; geometry follows the SPHinXsys
 *          test_2d_poiseuille_flow benchmark.
 */

#include "sphinxsys.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;

//----------------------------------------------------------------------
//  Case parameters (command line, same style as the parameter study).
//----------------------------------------------------------------------
struct WomersleyCase
{
    Real alpha{5.0};
    Real Re{100.0};
    Real dp_over_DH{1.0 / 40.0}; /**< particle spacing as a fraction of the channel height */
    int cycles{8};
    int analysis_cycles{3};
    std::string case_name{"womersley"};

    Real DH{1.0};    /**< channel height */
    Real DL{2.0};    /**< periodic length; only needs to be a few dp wide */
    Real rho0_f{1.0};
    Real U_ref{1.0}; /**< target peak cross-sectional mean velocity */

    // derived
    Real R{0.0}, dp{0.0}, BW{0.0}, nu{0.0}, mu_f{0.0}, omega{0.0}, T{0.0};
    Real G{0.0}, c_f{0.0}, delta{0.0}, end_time{0.0}, analysis_start{0.0};
    Real sample_interval{0.0};

    void derive()
    {
        R = 0.5 * DH;
        dp = DH * dp_over_DH;
        BW = 4.0 * dp;
        nu = U_ref * DH / Re;
        mu_f = rho0_f * nu;
        omega = nu * (alpha / R) * (alpha / R);
        T = 2.0 * Pi / omega;
        delta = std::sqrt(2.0 * nu / omega);

        // Choose G so that the peak sectional-mean velocity is U_ref. The complex
        // flow-rate response per unit forcing is
        //     Qhat/G = (-i / (rho w)) 2R [1 - tanh(lam R)/(lam R)],  lam = sqrt(i w / nu)
        // so |U_mean| = |Qhat| / (2R) and G = U_ref / |Qhat/G / (2R)|.
        std::complex<Real> lam = std::sqrt(std::complex<Real>(0.0, omega / nu));
        std::complex<Real> lamR = lam * R;
        std::complex<Real> resp =
            (std::complex<Real>(0.0, -1.0) / (rho0_f * omega)) *
            (Real(2.0) * R) * (std::complex<Real>(1.0, 0.0) - std::tanh(lamR) / lamR);
        G = U_ref / (std::abs(resp) / (2.0 * R));

        c_f = 10.0 * 1.5 * U_ref;
        end_time = Real(cycles) * T;
        analysis_start = end_time - Real(analysis_cycles) * T;
        sample_interval = T / 100.0;
    }

    /** Analytical complex velocity amplitude at height y (measured from centreline). */
    std::complex<Real> uHat(Real y) const
    {
        std::complex<Real> lam = std::sqrt(std::complex<Real>(0.0, omega / nu));
        return (std::complex<Real>(0.0, -1.0) * G / (rho0_f * omega)) *
               (std::complex<Real>(1.0, 0.0) - std::cosh(lam * y) / std::cosh(lam * R));
    }
};

WomersleyCase g_case;

//----------------------------------------------------------------------
//  Oscillatory driving force. GravityForce already passes the physical time to
//  InducedAcceleration, so an oscillatory "gravity" is all that is needed.
//----------------------------------------------------------------------
class OscillatoryGravity : public Gravity
{
    Real g0_, omega_;

  public:
    OscillatoryGravity(Real g0, Real omega)
        : Gravity(Vecd(g0, 0.0)), g0_(g0), omega_(omega) {}

    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        // Purely oscillatory, with no mean: the analytical solution above is the
        // response to G cos(omega t), and starting at a maximum avoids an
        // impulsive start in the velocity (which starts from rest at a phase where
        // the flow rate is zero).
        return Vecd(g0_ * cos(omega_ * physical_time), 0.0);
    }
};

//----------------------------------------------------------------------
//  Geometry: straight channel, periodic in x.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &name) : MultiPolygonShape(name)
    {
        const Real DL = g_case.DL, DH = g_case.DH;
        multi_polygon_.addPolygon({Vecd(0.0, 0.0), Vecd(0.0, DH), Vecd(DL, DH),
                                   Vecd(DL, 0.0), Vecd(0.0, 0.0)},
                                  GeometricOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &name) : MultiPolygonShape(name)
    {
        const Real DL = g_case.DL, DH = g_case.DH, BW = g_case.BW;
        multi_polygon_.addPolygon({Vecd(-BW, -BW), Vecd(-BW, DH + BW),
                                   Vecd(DL + BW, DH + BW), Vecd(DL + BW, -BW),
                                   Vecd(-BW, -BW)},
                                  GeometricOps::add);
        multi_polygon_.addPolygon({Vecd(-2.0 * BW, 0.0), Vecd(-2.0 * BW, DH),
                                   Vecd(DL + 2.0 * BW, DH), Vecd(DL + 2.0 * BW, 0.0),
                                   Vecd(-2.0 * BW, 0.0)},
                                  GeometricOps::sub);
    }
};

//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // -- minimal option parsing, mirroring case_params.h ------------------
    std::vector<char *> filtered{av[0]};
    for (int i = 1; i < ac; ++i)
    {
        std::string a(av[i]);
        auto eq = a.find('=');
        std::string key = a.rfind("--", 0) == 0
                              ? a.substr(2, eq == std::string::npos ? std::string::npos : eq - 2)
                              : "";
        std::string val = eq == std::string::npos ? "" : a.substr(eq + 1);
        if (key == "alpha") g_case.alpha = std::stod(val);
        else if (key == "Re") g_case.Re = std::stod(val);
        else if (key == "n_across") g_case.dp_over_DH = 1.0 / std::stod(val);
        else if (key == "cycles") g_case.cycles = std::stoi(val);
        else if (key == "case") g_case.case_name = val;
        else filtered.push_back(av[i]);
    }
    g_case.derive();
    const WomersleyCase &p = g_case;

    std::cout << "\n===== Womersley channel validation =====\n"
              << std::fixed << std::setprecision(5)
              << "  alpha        = " << p.alpha << "\n"
              << "  Re           = " << p.Re << "\n"
              << "  nu           = " << p.nu << "\n"
              << "  omega        = " << p.omega << "   T = " << p.T << "\n"
              << "  forcing G    = " << p.G << "\n"
              << "  dp           = " << p.dp << "   (" << 1.0 / p.dp_over_DH << " across DH)\n"
              << "  Stokes delta = " << p.delta << "   delta/dp = " << p.delta / p.dp << "\n"
              << "  cycles       = " << p.cycles << "   end_time = " << p.end_time << "\n"
              << "========================================\n\n";

    BoundingBoxd bounds(Vec2d(-p.BW, -p.BW), Vec2d(p.DL + p.BW, p.DH + p.BW));
    SPHSystem sph_system(bounds, p.dp);
    IO::getEnvironment().resetOutputFolder("./output_" + p.case_name);
    sph_system.handleCommandlineOptions(int(filtered.size()), filtered.data());

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(p.rho0_f, p.c_f);
    water_block.addMaterialProperty<Viscosity>(p.mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    // A column of observers spanning the channel at mid-length, inset by half a
    // particle spacing from each wall.
    //
    // The inset is not cosmetic. An observer placed exactly ON the wall has only
    // fluid particles on one side of it, so the Shepard-normalised interpolation
    // averages a one-sided stencil and cannot return the no-slip value: measured
    // that way the wall point reports |u| = 0.17 where the analytical answer is 0,
    // which alone inflates the L-infinity error of an otherwise 0.3 %-accurate
    // profile to 14 %. Sampling from the first particle position inwards is the
    // standard remedy and keeps every probe fully supported.
    const int n_obs = 41;
    StdVec<Vecd> obs_points;
    const Real y_lo = 0.5 * p.dp, y_hi = p.DH - 0.5 * p.dp;
    for (int i = 0; i < n_obs; ++i)
        obs_points.emplace_back(
            Vecd(0.5 * p.DL, y_lo + (y_hi - y_lo) * Real(i) / Real(n_obs - 1)));
    ObserverBody profile_observer(sph_system, "ProfileObserver");
    profile_observer.generateParticles<ObserverParticles>(obs_points);

    InnerRelation water_inner(water_block);
    ContactRelation water_wall(water_block, {&wall_boundary});
    ComplexRelation water_complex(water_inner, water_wall);
    ContactRelation profile_contact(profile_observer, {&water_block});

    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal(wall_boundary);
    OscillatoryGravity oscillatory_gravity(p.G / p.rho0_f, p.omega);
    SimpleDynamics<GravityForce<OscillatoryGravity>> apply_forcing(water_block, oscillatory_gravity);

    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);

    // Numerics copied from the validated SPHinXsys test_2d_poiseuille_flow
    // benchmark, which is the same periodic-channel configuration. In particular
    // the AllParticles variant of the transport correction is required: the
    // BulkParticles variant needs the free-surface Indicator, and this domain is
    // fully enclosed, so no free-surface indicator is ever computed.
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>
        kernel_correction(water_inner, water_wall);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann>
        pressure_relaxation(water_inner, water_wall);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann>
        density_relaxation(water_inner, water_wall);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWallCorrection>
        viscous_force(water_inner, water_wall);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionCorrectedComplex<AllParticles>>
        transport_correction(water_inner, water_wall);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density(water_inner, water_wall);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_advection_dt(water_block, p.U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_acoustic_dt(water_block);

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(water_block, "Velocity");
    ObservedQuantityRecording<Vecd> write_profile("Velocity", profile_contact);

    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_normal.exec();
    kernel_correction.exec();
    apply_forcing.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t iterations = 0;
    Real next_sample = 0.0, dt = 0.0;

    write_profile.writeToFile(iterations);
    write_states.writeToFile();

    while (physical_time < p.end_time)
    {
        Real Dt = get_advection_dt.exec();
        update_density.exec();
        kernel_correction.exec();
        // The forcing is TIME DEPENDENT, so it must be re-evaluated every step.
        // (test_2d_poiseuille_flow calls its gravity once before the loop because
        // there the force is constant; copying that pattern here freezes the
        // forcing at G cos(0) = G and produces a steady flow with no oscillation
        // at all - the first harmonic collapses to noise.)
        apply_forcing.exec();
        viscous_force.exec();
        transport_correction.exec();

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            dt = SMIN(get_acoustic_dt.exec(), Dt - relaxation_time);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            relaxation_time += dt;
            physical_time += dt;
        }
        iterations++;
        if (iterations % 500 == 0)
            std::cout << "N=" << iterations << "  t/T=" << physical_time / p.T << "\n";

        periodic_condition.bounding_.exec();
        water_block.updateCellLinkedList();
        periodic_condition.update_cell_linked_list_.exec();
        water_complex.updateConfiguration();

        if (physical_time >= next_sample)
        {
            profile_contact.updateConfiguration();
            write_profile.writeToFile(iterations);
            while (next_sample <= physical_time)
                next_sample += p.sample_interval;
        }
    }
    write_states.writeToFile();

    // -- metadata + the analytical solution on the observer points ---------
    std::ofstream f("./output_" + p.case_name + "/case_params.json");
    f << std::setprecision(12) << "{\n";
    f << "  \"kind\": \"womersley_channel\",\n";
    f << "  \"case_name\": \"" << p.case_name << "\",\n";
    f << "  \"alpha\": " << p.alpha << ",\n";
    f << "  \"Re\": " << p.Re << ",\n";
    f << "  \"nu\": " << p.nu << ",\n";
    f << "  \"mu_f\": " << p.mu_f << ",\n";
    f << "  \"rho0_f\": " << p.rho0_f << ",\n";
    f << "  \"omega\": " << p.omega << ",\n";
    f << "  \"T\": " << p.T << ",\n";
    f << "  \"G\": " << p.G << ",\n";
    f << "  \"R\": " << p.R << ",\n";
    f << "  \"DH\": " << p.DH << ",\n";
    f << "  \"dp\": " << p.dp << ",\n";
    f << "  \"delta_stokes\": " << p.delta << ",\n";
    f << "  \"delta_over_dp\": " << p.delta / p.dp << ",\n";
    f << "  \"U_ref\": " << p.U_ref << ",\n";
    f << "  \"c_f\": " << p.c_f << ",\n";
    f << "  \"t_ramp\": 0,\n";
    f << "  \"n_cycles\": " << p.cycles << ",\n";
    f << "  \"analysis_cycles\": " << p.analysis_cycles << ",\n";
    f << "  \"analysis_start\": " << p.analysis_start << ",\n";
    f << "  \"end_time\": " << p.end_time << ",\n";
    f << "  \"y\": [";
    for (int i = 0; i < n_obs; ++i)
        f << (i ? ", " : "") << obs_points[i][1];
    f << "],\n";
    // analytical first harmonic at each observer point, as (real, imag) pairs
    f << "  \"u_hat_analytic\": [";
    for (int i = 0; i < n_obs; ++i)
    {
        std::complex<Real> uh = p.uHat(obs_points[i][1] - p.R);
        f << (i ? ", " : "") << "[" << uh.real() << ", " << uh.imag() << "]";
    }
    f << "]\n}\n";
    f.close();

    std::cout << "\nwrote ./output_" << p.case_name << "/case_params.json\n";
    return 0;
}
