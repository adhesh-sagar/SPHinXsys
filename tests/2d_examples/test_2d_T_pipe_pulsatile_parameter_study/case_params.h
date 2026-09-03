/**
 * @file    case_params.h
 * @brief   Case definition, command-line parsing and metadata output for the
 *          pulsatile T-pipe parameter study.
 *
 * @details The case is defined by THREE dimensionless groups and nothing else:
 *
 *            Re    = rho * U_f * DH / mu          (inertia / viscous)
 *            alpha = R * sqrt(omega / nu)         (Womersley: unsteady inertia / viscous)
 *            A     = dU / U_mean                  (pulsatility amplitude)
 *
 *          Re and alpha are COUPLED through the kinematic viscosity, so they must
 *          not be varied by editing omega directly. Instead the definitions are
 *          inverted, and omega is a DERIVED quantity:
 *
 *            nu    = U_f * DH / Re
 *            omega = nu * (alpha / R)^2 ,   R = DH / 2
 *            T     = 2 * pi / omega
 *
 *          Everything in the study is then reported in dimensionless form
 *          (u/U_f, t/T, p/(rho U_f^2), tau_w/(mu U_f / R)).
 *
 * @note    SPHSystem::handleCommandlineOptions() builds a FIXED Boost
 *          options_description and calls exit(1) on any option it does not know.
 *          We therefore parse our own options here and hand SPHinXsys a FILTERED
 *          argv containing only the flags it recognises. The library is untouched.
 */

#ifndef CASE_PARAMS_H
#define CASE_PARAMS_H

#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace TPipeStudy
{
using SPH::Real;

/** Options understood by SPHinXsys itself; anything else that starts with "--"
 *  and is not one of ours is passed through untouched so that new library
 *  options keep working without editing this list. */
inline bool isOurOption(const std::string &key)
{
    static const std::vector<std::string> ours{
        "Re", "alpha", "A", "dp", "cycles", "case", "branch_ratio",
        "stenosis", "cfl_scale", "c_f", "carreau", "samples_per_cycle",
        "vtp_per_cycle", "analysis_cycles", "help_case"};
    for (const auto &o : ours)
        if (key == o)
            return true;
    return false;
}

/**
 * @struct CaseParams
 * @brief  One run of the study: dimensionless inputs plus everything derived.
 */
struct CaseParams
{
    //---- dimensionless inputs (command line) --------------------------------
    Real Re{100.0};           /**< Reynolds number based on DH and U_f. */
    Real alpha{5.0};          /**< Womersley number; omega is derived from it. */
    Real A{0.5};              /**< Pulsatility amplitude as a fraction of U_f. */
    Real dp{0.15};            /**< Particle spacing (resolution study). */
    int cycles{0};            /**< Number of analysed+transient cycles; 0 = auto. */
    int analysis_cycles{3};   /**< Trailing cycles treated as the analysis window. */
    Real branch_ratio{1.0};   /**< Lower-branch length ratio (extension D). */
    Real stenosis{0.0};       /**< Upper-branch area reduction, 0..0.9 (extension E). */
    Real cfl_scale{1.0};      /**< Multiplies the acoustic time step (CFL study). */
    Real c_f_override{0.0};   /**< >0 pins the artificial sound speed (regression). */
    bool carreau{false};      /**< Reserved: shear-thinning blood rheology. */
    int samples_per_cycle{200};
    int vtp_per_cycle{20};
    std::string case_name{"default"};

    /** Allowance for flow acceleration around the junction corners, used when
     *  sizing the artificial sound speed. Calibrated against MaxBulkSpeed (see
     *  derive()): the peak bulk speed is ~1.2x the prescribed inlet peak. */
    static constexpr Real CORNER_FACTOR = 1.25;

    //---- fixed reference quantities ----------------------------------------
    Real DL{5.0};     /**< Domain reference length. */
    Real DH{3.0};     /**< Inlet channel height (the reference diameter). */
    Real U_f{1.0};    /**< Mean inflow velocity. */
    Real rho0_f{1.0}; /**< Reference density. */

    //---- derived (filled by derive()) --------------------------------------
    Real DL1{0.0};       /**< Length of the horizontal inlet channel. */
    Real BW{0.0};        /**< Wall boundary width. */
    Real DL_sponge{0.0}; /**< Emitter buffer length. */
    Real R{0.0};         /**< Channel half-height, DH/2. */
    Real nu{0.0};        /**< Kinematic viscosity. */
    Real mu_f{0.0};      /**< Dynamic viscosity. */
    Real c_f{0.0};       /**< Artificial sound speed. */
    Real omega{0.0};     /**< Pulsation angular frequency (DERIVED from alpha). */
    Real T{0.0};         /**< Pulse period. */
    Real delta_stokes{0.0}; /**< Stokes layer thickness sqrt(2 nu / omega). */
    Real tau_visc{0.0};  /**< Viscous diffusion time R^2 / nu. */
    Real u_peak{0.0};    /**< Peak inlet centreline speed 1.5 U_f (1 + A). */
    Real mach{0.0};      /**< u_peak / c_f, must stay <= 0.1 for WCSPH. */
    Real t_ramp{0.0};    /**< Cosine ramp-up duration. */
    Real end_time{0.0};
    Real analysis_start{0.0}; /**< Physical time at which the analysis window opens. */
    Real sample_interval{0.0};
    Real vtp_interval{0.0};
    Real ke_interval{0.0};
    int n_cycles{0};
    int decay_cycles{0};        /**< Cycles estimated to kill the start-up transient. */
    bool periodicity_ok{false}; /**< n_cycles * T >= 2.5 * tau_visc. */
    bool stokes_resolved{false};/**< delta_stokes >= 4 * dp. */

    //----------------------------------------------------------------------
    /** Invert the dimensionless definitions and fill every derived member. */
    void derive()
    {
        DL1 = 0.7 * DL;
        BW = 4.0 * dp;
        DL_sponge = 20.0 * dp;
        R = 0.5 * DH;

        nu = U_f * DH / Re;
        mu_f = rho0_f * nu;

        // omega is DERIVED from the target Womersley number, never prescribed.
        omega = nu * (alpha / R) * (alpha / R);
        T = 2.0 * SPH::Pi / omega;
        delta_stokes = std::sqrt(2.0 * nu / omega);
        tau_visc = R * R / nu;

        // Peak BULK speed, which is what sets the sound speed.
        //
        //   1.5      parabolic profile peak over the section mean
        //   (1 + A)  pulsatile peak over the cycle mean
        //   1.25     measured allowance for acceleration at the junction corners
        //
        // The last factor was calibrated against MaxBulkSpeed, not the library's
        // MaximumSpeed. MaximumSpeed reduces over every particle including those in
        // the inflow buffer and the outflow disposers, where it reaches absurd
        // values (70-400 in a flow of order 1) that reflect a handful of particles
        // about to be deleted rather than the solution. Sizing c_f against that
        // number inflates it by 2x, which measurably DESTABILISES the run: at
        // c_f = 45 two runs in three failed, at c_f = 22.5 none did.
        u_peak = CORNER_FACTOR * 1.5 * U_f * (1.0 + A);
        c_f = c_f_override > 0.0
                  ? c_f_override
                  : 10.0 * U_f * SPH::SMAX(u_peak, DH / (Real(2.0) * (DL - DL1)));
        mach = u_peak / c_f;

        // Number of cycles needed before the flow becomes periodic.
        //
        // The slowest-decaying start-up mode of a plane channel decays at
        // lambda_1 = nu (pi / 2R)^2 = (pi^2 / 4) / tau_visc, i.e. with e-folding time
        // 0.405 * tau_visc. Reaching a 0.2 % residual therefore needs
        //     t_decay ~ 2.5 * tau_visc,
        // and since tau_visc / T = alpha^2 / (2 pi),
        //     decay_cycles = ceil(2.5 * alpha^2 / (2 pi)) = ceil(0.398 * alpha^2).
        //
        // This is the whole reason a fixed cycle count cannot work: alpha = 2 is
        // periodic after 2 cycles while alpha = 10 needs 40.
        decay_cycles = int(std::ceil(2.5 * alpha * alpha / (2.0 * SPH::Pi)));
        if (cycles > 0)
        {
            n_cycles = cycles;
        }
        else
        {
            n_cycles = SPH::SMAX(analysis_cycles + 2, decay_cycles + analysis_cycles);
        }

        t_ramp = 2.0 * T; // scale the ramp with the period, never a fixed 2.0 s
        end_time = t_ramp + Real(n_cycles) * T;
        analysis_start = end_time - Real(analysis_cycles) * T;

        sample_interval = T / Real(samples_per_cycle);
        vtp_interval = T / Real(vtp_per_cycle);
        ke_interval = T / 50.0;

        periodicity_ok = Real(n_cycles) * T >= 2.5 * tau_visc;

        // The oscillatory (Stokes) boundary layer is the smallest physical length in
        // the problem: delta = sqrt(2 nu / omega) = R * sqrt(2) / alpha, which is
        // INDEPENDENT of Re. It therefore sets the production resolution:
        //   >= 4 particles across delta  <=>  dp <= 0.53 * R / alpha.
        stokes_resolved = delta_stokes >= 4.0 * dp;
    }

    /** Largest dp that resolves the Stokes layer with 4 particles - used in messages. */
    Real dpForStokes() const { return delta_stokes / 4.0; }

    //----------------------------------------------------------------------
    /** Human-readable banner, including the two checks that can invalidate a run. */
    void print() const
    {
        std::cout << "\n===== Pulsatile T-pipe parameter study =====\n"
                  << std::fixed << std::setprecision(5)
                  << "  case              = " << case_name << "\n"
                  << "  --- dimensionless inputs ---\n"
                  << "  Re                = " << Re << "\n"
                  << "  alpha (Womersley) = " << alpha << "\n"
                  << "  A (pulsatility)   = " << A << "\n"
                  << "  dp (resolution)   = " << dp << "\n"
                  << "  branch_ratio      = " << branch_ratio << "\n"
                  << "  stenosis          = " << stenosis << "\n"
                  << "  --- derived ---\n"
                  << "  nu                = " << nu << "\n"
                  << "  mu_f              = " << mu_f << "\n"
                  << "  omega             = " << omega << " rad/s\n"
                  << "  T (period)        = " << T << "\n"
                  << "  Stokes layer d    = " << delta_stokes
                  << "   (d/dp = " << delta_stokes / dp << ")\n"
                  << "  tau_visc = R^2/nu = " << tau_visc << "\n"
                  << "  c_f               = " << c_f << "\n"
                  << "  --- run control ---\n"
                  << "  n_cycles          = " << n_cycles
                  << "   (transient decay needs " << decay_cycles << ")\n"
                  << "  t_ramp            = " << t_ramp << "\n"
                  << "  end_time          = " << end_time << "\n"
                  << "  analysis window   = [" << analysis_start << ", " << end_time << "]"
                  << "  (" << analysis_cycles << " cycles)\n"
                  << "  sample_interval   = " << sample_interval
                  << "   (" << samples_per_cycle << " samples/cycle)\n";

        std::cout << "  --- validity checks ---\n";
        std::cout << "  Mach = u_peak/c_f  = " << mach
                  << (mach <= 0.1 ? "   [OK, <= 0.1]\n"
                                  : "   [WARNING: > 0.1, weak compressibility violated]\n");
        std::cout << "  n_cycles*T >= 2.5*tau_visc? " << Real(n_cycles) * T
                  << " vs " << 2.5 * tau_visc
                  << (periodicity_ok ? "   [OK]\n"
                                     : "   [WARNING: transient may not have decayed]\n");
        // The lattice generator fills the channel with whole particle layers. If
        // DH/dp is not an integer the leftover partial layer biases every
        // volume-weighted average, and the inlet flow-rate error can flip sign
        // between resolutions - which destroys the monotonicity a convergence
        // study depends on.
        const Real layers = DH / dp;
        const Real layer_defect = std::abs(layers - std::round(layers));
        std::cout << "  DH / dp (layers)   = " << layers
                  << (layer_defect < 1e-6
                          ? "   [OK, integer]\n"
                          : "   [WARNING: not an integer - partial near-wall layer will\n"
                            "                          bias volume-weighted averages; prefer a dp\n"
                            "                          that divides DH exactly]\n");

        std::cout << "  delta_stokes / dp  = " << delta_stokes / dp;
        if (stokes_resolved)
            std::cout << "   [OK, >= 4]\n";
        else
            std::cout << "   [WARNING: < 4, use dp <= " << dpForStokes()
                      << " to resolve the oscillatory boundary layer]\n";
        std::cout << "============================================\n\n";
    }

    //----------------------------------------------------------------------
    /** Machine-readable metadata so the Python side never re-derives anything. */
    void writeJson(const std::string &path, int n_fluid_particles, Real wall_time_s,
                   const std::string &probe_json) const
    {
        std::ofstream f(path.c_str());
        f << std::setprecision(12);
        f << "{\n";
        f << "  \"kind\": \"t_pipe_study\",\n";
        f << "  \"case_name\": \"" << case_name << "\",\n";
        f << "  \"Re\": " << Re << ",\n";
        f << "  \"alpha\": " << alpha << ",\n";
        f << "  \"A\": " << A << ",\n";
        f << "  \"dp\": " << dp << ",\n";
        f << "  \"branch_ratio\": " << branch_ratio << ",\n";
        f << "  \"stenosis\": " << stenosis << ",\n";
        f << "  \"cfl_scale\": " << cfl_scale << ",\n";
        f << "  \"carreau\": " << (carreau ? "true" : "false") << ",\n";
        f << "  \"DL\": " << DL << ",\n";
        f << "  \"DH\": " << DH << ",\n";
        f << "  \"DL1\": " << DL1 << ",\n";
        f << "  \"DL_sponge\": " << DL_sponge << ",\n";
        f << "  \"BW\": " << BW << ",\n";
        f << "  \"R\": " << R << ",\n";
        f << "  \"U_f\": " << U_f << ",\n";
        f << "  \"rho0_f\": " << rho0_f << ",\n";
        f << "  \"nu\": " << nu << ",\n";
        f << "  \"mu_f\": " << mu_f << ",\n";
        f << "  \"c_f\": " << c_f << ",\n";
        f << "  \"omega\": " << omega << ",\n";
        f << "  \"T\": " << T << ",\n";
        f << "  \"delta_stokes\": " << delta_stokes << ",\n";
        f << "  \"tau_visc\": " << tau_visc << ",\n";
        f << "  \"u_peak\": " << u_peak << ",\n";
        f << "  \"mach\": " << mach << ",\n";
        f << "  \"t_ramp\": " << t_ramp << ",\n";
        f << "  \"end_time\": " << end_time << ",\n";
        f << "  \"analysis_start\": " << analysis_start << ",\n";
        f << "  \"analysis_cycles\": " << analysis_cycles << ",\n";
        f << "  \"n_cycles\": " << n_cycles << ",\n";
        f << "  \"decay_cycles\": " << decay_cycles << ",\n";
        f << "  \"stokes_resolved\": " << (stokes_resolved ? "true" : "false") << ",\n";
        f << "  \"delta_stokes_over_dp\": " << delta_stokes / dp << ",\n";
        f << "  \"layers_across_DH\": " << DH / dp << ",\n";
        f << "  \"integer_layers\": "
          << (std::abs(DH / dp - std::round(DH / dp)) < 1e-6 ? "true" : "false") << ",\n";
        f << "  \"samples_per_cycle\": " << samples_per_cycle << ",\n";
        f << "  \"sample_interval\": " << sample_interval << ",\n";
        f << "  \"periodicity_ok\": " << (periodicity_ok ? "true" : "false") << ",\n";
        f << "  \"n_fluid_particles\": " << n_fluid_particles << ",\n";
        f << "  \"wall_time_s\": " << wall_time_s << ",\n";
        f << probe_json; // observer point tables, already indented + comma-terminated
        f << "  \"schema\": 1\n";
        f << "}\n";
        f.close();
    }
};

//----------------------------------------------------------------------
/** Usage text for our own options. */
inline void printCaseUsage()
{
    std::cout <<
        R"(Pulsatile T-pipe parameter study - case options
  --Re=<Real>                Reynolds number rho U DH / mu           [100]
  --alpha=<Real>             Womersley number; omega is derived      [5]
  --A=<Real>                 pulsatility amplitude dU/U              [0.5]
  --dp=<Real>                particle spacing                        [0.15]
  --cycles=<int>             cycles to run; 0 = auto from alpha      [0]
  --analysis_cycles=<int>    trailing cycles used for analysis       [3]
  --branch_ratio=<Real>      lower branch length ratio               [1.0]
  --stenosis=<Real>          upper branch area reduction 0..0.9      [0.0]
  --cfl_scale=<Real>         multiplies the acoustic time step       [1.0]
  --c_f=<Real>               pin the sound speed (regression only)   [auto]
  --samples_per_cycle=<int>  observer samples per cycle              [200]
  --vtp_per_cycle=<int>      VTP snapshots per cycle (last window)   [20]
  --carreau=<0|1>            shear-thinning blood rheology           [0]
  --case=<name>              output goes to ./output_<name>          [default]
  --help_case                print this and exit
All other options (--relax, --reload, --regression, --state_recording,
--restart_step, --log_level, --help) are forwarded to SPHinXsys unchanged.
)";
}

//----------------------------------------------------------------------
/**
 * @brief Split argv into our case options and the ones SPHinXsys must see.
 * @param[out] params        filled from our options
 * @param[out] filtered_argv argv containing argv[0] plus only foreign options
 *
 * Accepts both "--key=value" and "--key value".
 */
inline void parseCaseOptions(int ac, char *av[], CaseParams &params,
                             std::vector<char *> &filtered_argv)
{
    filtered_argv.clear();
    filtered_argv.push_back(av[0]);

    auto as_real = [](const std::string &s) { return Real(std::stod(s)); };
    auto as_int = [](const std::string &s) { return std::stoi(s); };
    auto as_bool = [](const std::string &s) { return !(s == "0" || s == "false" || s == "off"); };

    for (int i = 1; i < ac; ++i)
    {
        std::string arg(av[i]);
        if (arg.rfind("--", 0) != 0)
        {
            filtered_argv.push_back(av[i]); // positional, not ours
            continue;
        }

        std::string body = arg.substr(2);
        std::string key = body;
        std::string value;
        bool has_inline_value = false;
        const auto eq = body.find('=');
        if (eq != std::string::npos)
        {
            key = body.substr(0, eq);
            value = body.substr(eq + 1);
            has_inline_value = true;
        }

        if (!isOurOption(key))
        {
            filtered_argv.push_back(av[i]); // let Boost handle it
            continue;
        }

        if (key == "help_case")
        {
            printCaseUsage();
            std::exit(0);
        }

        // "--key value" form: consume the following token as the value.
        if (!has_inline_value)
        {
            if (i + 1 < ac && std::string(av[i + 1]).rfind("--", 0) != 0)
            {
                value = av[++i];
            }
            else
            {
                std::cerr << "error: case option --" << key << " requires a value\n";
                std::exit(1);
            }
        }

        try
        {
            if (key == "Re")
                params.Re = as_real(value);
            else if (key == "alpha")
                params.alpha = as_real(value);
            else if (key == "A")
                params.A = as_real(value);
            else if (key == "dp")
                params.dp = as_real(value);
            else if (key == "cycles")
                params.cycles = as_int(value);
            else if (key == "analysis_cycles")
                params.analysis_cycles = as_int(value);
            else if (key == "branch_ratio")
                params.branch_ratio = as_real(value);
            else if (key == "stenosis")
                params.stenosis = as_real(value);
            else if (key == "cfl_scale")
                params.cfl_scale = as_real(value);
            else if (key == "c_f")
                params.c_f_override = as_real(value);
            else if (key == "samples_per_cycle")
                params.samples_per_cycle = as_int(value);
            else if (key == "vtp_per_cycle")
                params.vtp_per_cycle = as_int(value);
            else if (key == "carreau")
                params.carreau = as_bool(value);
            else if (key == "case")
                params.case_name = value;
        }
        catch (const std::exception &e)
        {
            std::cerr << "error: could not parse --" << key << "=" << value
                      << " (" << e.what() << ")\n";
            std::exit(1);
        }
    }

    // Guard rails: these ranges are what the geometry and the WCSPH model support.
    if (params.Re <= 0.0 || params.alpha <= 0.0 || params.dp <= 0.0)
    {
        std::cerr << "error: Re, alpha and dp must all be positive\n";
        std::exit(1);
    }
    if (params.A < 0.0)
    {
        std::cerr << "error: A must be >= 0 (A = 0 is the steady validation case)\n";
        std::exit(1);
    }
    if (params.A > 1.0)
    {
        // At A > 1 the prescribed inlet velocity goes negative over part of the
        // cycle. The emitter + inflow buffer inherited from T_shaped_pipe is
        // unidirectional and cannot accept backflow: tested at A = 1.5 the run
        // collapses to Dt ~ 3e-5 and stalls, while A <= 1.0 runs cleanly.
        // Genuine flow reversal needs a bidirectional buffer (see
        // tests/extra_source_and_tests/extra_src/shared/pressure_boundary/).
        std::cerr << "error: A = " << params.A << " > 1 requires a bidirectional inlet.\n"
                  << "       This inlet treatment is unidirectional; the run will stall.\n"
                  << "       Use A <= 1.0, or implement a bidirectional buffer.\n";
        std::exit(1);
    }
    if (params.stenosis < 0.0 || params.stenosis > 0.9)
    {
        std::cerr << "error: stenosis must be in [0, 0.9]\n";
        std::exit(1);
    }
    if (params.branch_ratio <= 0.2 || params.branch_ratio > 1.0)
    {
        std::cerr << "error: branch_ratio must be in (0.2, 1.0]\n";
        std::exit(1);
    }
    if (params.analysis_cycles < 1)
    {
        std::cerr << "error: analysis_cycles must be >= 1\n";
        std::exit(1);
    }

    params.derive();

    if (params.cycles > 0 && params.cycles <= params.analysis_cycles)
    {
        std::cerr << "error: cycles (" << params.cycles << ") must exceed analysis_cycles ("
                  << params.analysis_cycles << ")\n";
        std::exit(1);
    }
}

} // namespace TPipeStudy

#endif // CASE_PARAMS_H
