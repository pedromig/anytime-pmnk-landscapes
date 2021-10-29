/**
 * @file app.cpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Simple app to test the apmnkl library implementation of some search
 *        heuristics. This app serves as an example program showing the usage
 *        of the library, while also simplifing the gathering of anytime data
 *        relevant to the study of the performance of the implemented algorithms
 *        in the context of the pmnk-landscapes problem.
 * @version 0.2.0
 * @date 29-10-2021
 *
 * @copyright Copyright (c) 2021
 *
 */

// CLI11 Command Line Parser
#include <CLI/CLI.hpp>

// anytime pmnk-landscapes (apmnkl) library includes
#include <apmnkl/gsemo.hpp>
#include <apmnkl/ibea.hpp>
#include <apmnkl/operators.hpp>
#include <apmnkl/pls.hpp>

// Standard Includes
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

// Helper IBEA Subcommand CLI Enums

/**
 * @brief Enum representing which indicator to use as the IBEA algorithm indicator operator.
 *        Possible values are:
 *          - eps (Additive Epsilon Indicator)
 *          - ihd (Hypervolume Based Indicator)
 */
enum class indicator { eps, ihd };

/**
 * @brief Enum representing which crossover method to use as the IBEA algorithm
 *        crossover operator.
 *        Possible values are:
 *          - npc (N-Point Crossover)
 *          - uc (Uniform Crossover)
 */
enum class crossover { npc, uc };

/**
 * @brief Enum representing which mutation operator to use as the IBEA algorithm
 *        mutation operator.
 *        Possible values are:
 *          - um (Uniform Mutation)
 */
enum class mutation { um };

/**
 * @brief Enum representing which selection method to use as the IBEA algorithm
 *        selection operator.
 *        Possible values are:
 *          - kwt (K-Way Tournament Selection)
 */
enum class selection { kwt };

// CLI Options

/**
 * @brief Set the CLI positional arguments options/flags
 *
 * @param app CLI::App object that will hold all the options/flags of the positional arguments *
 * (below).
 * @param instance The path for the "rmnk" instance (.dat) file to be used.
 *        These files can be generated using the rmnkGenerator.R script
 */
inline void set_positional_arguments(CLI::App &app, std::string &instance) {
  app.add_option("instance", instance,
                 "= pmnk-landscapes instance file path\n(instances can be generated using the "
                 "rmnkGenerator.R script).")
      ->type_name(".dat")
      ->check(CLI::ExistingFile)
      ->required()
      ->group("Positionals");
}

/**
 * @brief Set the CLI general options/flags
 *
 * @param app CLI::App object that will hold all the global options/flags (below).
 * @param maxeval The maximum number of evaluations performed by the algorithms (stopping criterion)
 * @param seed The seed used by the pseudo random number generator used in these algorithms
 * @param output The name of the output file where the standard output stream should be redirected
 * @param ref The reference point considered by hypervolume indicator whilst running the algorithms
 * (anytime measure)
 */
inline void set_general_options(CLI::App &app, std::size_t &maxeval, unsigned int &seed,
                                std::string &output, apmnkl::objective_vector &ref) {
  app.add_option("-m,--maxeval", maxeval,
                 "= maximum number of evaluations to be performed (stopping criterion).")
      ->needs(app.get_option("instance"))
      ->required()
      ->check(CLI::NonNegativeNumber)
      ->group("Options");

  app.add_option("-s,--seed", seed, "= pseudo random generator seed used by the search heuristics.")
      ->needs(app.get_option("instance"))
      ->check(CLI::NonNegativeNumber)
      ->group("Options");

  app.add_option("-o,--output", output,
                 "= specifies the file to which the output stream should be redirected.")
      ->needs(app.get_option("instance"))
      ->group("Options");

  app.add_option("-r,--hvref", ref, "= reference point considered in the hypervolume calculation.")
      ->needs(app.get_option("instance"))
      ->group("Options");

  app.set_help_flag("-h,--help", "= print this help message and exit.");
  app.set_help_all_flag("-H,--help-all", "= expand all help.");
}

/**
 * @brief Set the PLS algorithm options/flags
 *
 * @param app CLI::App object that will hold all the PLS options/flags (below).
 * @param pac The PLS algorithm solution acceptance criterion
 * @param pne The PLS algorithm neighboorhood solutions exploration criterion
 */
inline void set_pls_options(CLI::App &app, apmnkl::pls::pac &pac, apmnkl::pls::pne &pne) {
  std::map<std::string, apmnkl::pls::pac> acceptance_opts{
      {"NON_DOMINATING", apmnkl::pls::pac::non_dominating},
      {"DOMINATING", apmnkl::pls::pac::dominating},
      {"BOTH", apmnkl::pls::pac::both}};

  app.add_option(
         "-a,--pls-acceptance-criterion", pac,
         "= acceptance criterion considered whilst running pls.\n  => (NON_DOMINATING): accept "
         "every non-dominated neighbor. \n  => (DOMINATING): accept only neighbors that "
         "dominate current solution.\n  => (BOTH): first try to accept only neighbors that "
         "dominate\n      the current solution, if none exist accept non-dominated solutions.")
      ->transform(CLI::CheckedTransformer(acceptance_opts, CLI::ignore_case));

  std::map<std::string, apmnkl::pls::pne> exploration_opts{
      {"BEST_IMPROVEMENT", apmnkl::pls::pne::best_improvement},
      {"FIRST_IMPROVEMENT", apmnkl::pls::pne::first_improvement},
      {"BOTH", apmnkl::pls::pne::both}};

  app.add_option(
         "-e,--pls-neighborhood-exploration", pne,
         "= neighborhood exploration criterion considered whilst running pls.\n "
         "=> (BEST_IMPROVEMENT): explore every acceptable neighboor. \n => (FIRST_IMPROVEMENT):"
         " stop once on neighbor is accepted. \n => (BOTH): use FIRST_IMPROVEMENT until PLS stops,"
         " afterwards use BEST_IMPROVEMENT")
      ->transform(CLI::CheckedTransformer(exploration_opts, CLI::ignore_case));
}

/**
 * @brief Set the IBEA algorithm options/flags
 *
 * @param app CLI::App object that will hold all the IBEA options/flags (below).
 * @param population_size The maximum size the population
 * @param generations The maximum number of generations (stopping criterion)
 * @param scaling_factor The IBEA scaling factor
 * @param adaptive A boolean indicative of the version of the algorithm to be used.
 *                   True for B-IBEA (Basic IBEA) and false for A-IBEA (Adaptive IBEA)
 */
inline void set_ibea_options(CLI::App &app, std::size_t &population_size, std::size_t &generations,
                             double &scaling_factor, bool &adaptive) {
  app.add_option("-p,--pop-size", population_size, "= max population size.")
      ->required()
      ->check(CLI::NonNegativeNumber);

  app.add_option("-g,--generations", generations, "= number of generations. (stopping criterion)")
      ->required()
      ->check(CLI::NonNegativeNumber);

  app.add_option("-k,--scaling-factor", scaling_factor, "= scaling factor.")
      ->required()
      ->check(CLI::NonNegativeNumber);

  app.add_flag("-a,--adaptive", adaptive, "= use the adaptive version of the algorithm");
}

/**
 * @brief Set the IBEA algorithm mutation operators options/flags
 *
 * @param app CLI::App object that will hold all the IBEA mutation operator options/flags (below).
 * @param mutation_probability The probability of occurrence of mutations in the individual's
 * genotype
 */
inline void set_ibea_mutation_options(CLI::App &app, double &mutation_probability) {
  app.add_option("-p,--mutation-probability", mutation_probability,
                 "= probability of occurrence of a mutation in the individual's genotype")
      ->required()
      ->check(CLI::Range(0.0, 1.0));
}

/**
 * @brief Set the IBEA algorithm crossover operators options/flags
 *
 * @param app CLI::App object that will hold all the IBEA crossover operator options/flags (below).
 * @param crossover_probability The probability of occurrence of crossover between two individual's
 * genotypes
 * @param npoints If the crossover operator is the n-point crossover this option is set. This option
 * represents the number of crossover points to be considered when swapping both individual's
 * genetic content
 */
inline void set_ibea_crossover_options(CLI::App &app, double &crossover_probability,
                                       std::size_t &npoints) {
  app.add_option("-p,--crossover_probability", crossover_probability,
                 "= probability of occurrence of a mutation in the individual's genotype")
      ->required()
      ->check(CLI::Range(0.0, 1.0));

  if (app.get_name() == "NPointCrossover") {
    app.add_option("-n,--n-points", npoints, "= number of randomly picked crossover points")
        ->required()
        ->check(CLI::NonNegativeNumber);
  }
}

/**
 * @brief Set the ibea selection options object
 *
 * @param app CLI::App object that will hold all the IBEA selection operator options/flags (below).
 * @param matting_pool_size The size of the matting pool.
 * @param tournament_size The size of the tournament used for selection.
 */
inline void set_ibea_selection_options(CLI::App &app, std::size_t &matting_pool_size,
                                       std::size_t &tournament_size) {
  app.add_option("-s,--matting-pool-size", matting_pool_size,
                 "= target size of the matting pool to be obtained from the selection step")
      ->required()
      ->check(CLI::NonNegativeNumber);

  app.add_option("-t,--tournament-size", tournament_size,
                 "= size of the tournament used for individual selection")
      ->required()
      ->check(CLI::NonNegativeNumber);
}

// Write CSV ouput file (Utils)

/**
 * @brief Unpack tuple containing the elements of a row and dump
 *        them into the output stream in csv format. (Worker function)
 *
 * @tparam I the index of the current element being printed
 * @tparam Ts The types of the elements contained in the tuple
 * @param os The output stream where the data (in csv format)
 *           should be redirected to
 * @param row  The current row of data
 */
template <std::size_t I, typename... Ts>
auto csv_row(std::ostream &os, std::tuple<Ts...> const &row) {
  if constexpr (I == sizeof...(Ts)) {
    os << "\n";
  } else {
    os << std::setprecision(12) << std::get<I>(row);
    if constexpr (I + 1 != sizeof...(Ts)) {
      os << ",";
    }
    csv_row<I + 1>(os, row);
  }
}

/**
 * @brief Unpack tuple containing the elements of a row and dump
 *        them into the output stream in csv format. (Wrapper function)
 *
 * @tparam Ts The types of the elements contained in the tuple
 * @param os The output stream where the data (in csv format)
 *           should be redirected to
 * @param row  The current row of data
 */
template <typename... Ts>
auto write_csv_row(std::ostream &os, std::tuple<Ts...> const &row) {
  csv_row<0>(os, row);
}

/**
 * @brief Write anytime data gathered by the search heuristics and
 *        dump them into a stream in the csv format with delimiter=","
 *
 * @tparam R The type of for the container having the row data
 * @param os The output stream where the data (in csv format)
 *           should be redirected to
 * @param header The header row of the csv
 * @param data The payload (rows) with the anytime data
 */
template <typename R>
void to_csv(std::ostream &os, std::string const &header, std::vector<R> const &data) {
  os << header << '\n';
  for (auto &row : data) {
    write_csv_row(os, row);
  }
}

// Algorithm Callbacks

/**
 * @brief CLI::App callback for the gsemo algorithm
 *
 * @param instance The path for the "rmnk" instance (.dat) file to be used.
 *        These files can be generated using the rmnkGenerator.R script
 * @param maxeval The maximum number of evaluations performed by the algorithms (stopping criterion)
 * @param seed The seed used by the pseudo random number generator used in these algorithms
 * @param os The name of the output file where the standard output stream should be redirected
 * @param ref The reference point considered by hypervolume indicator whilst running the algorithms
 * (anytime measure)
 */
inline void gsemo(std::string const &instance, std::size_t const maxeval, unsigned int const seed,
                  std::ostream &os, apmnkl::objective_vector &ref) {
  if (ref.empty()) {
    apmnkl::gsemo gsemo(instance, seed);
    gsemo.run(maxeval);
    to_csv(os, "evaluation,hypervolume", gsemo.anytime());
  } else {
    apmnkl::gsemo gsemo(instance, seed, ref);
    gsemo.run(maxeval);
    to_csv(os, "evaluation,hypervolume", gsemo.anytime());
  }
}

/**
 * @brief CLI::App callback for the pls algorithm
 *
 * @param instance The path for the "rmnk" instance (.dat) file to be used.
 *        These files can be generated using the rmnkGenerator.R script
 * @param maxeval The maximum number of evaluations performed by the algorithms (stopping criterion)
 * @param seed The seed used by the pseudo random number generator used in these algorithms
 * @param pac The PLS algorithm solution acceptance criterion
 * @param pne The PLS algorithm neighboorhood solutions exploration criterion
 * @param os The name of the output file where the standard output stream should be redirected
 * @param ref The reference point considered by hypervolume indicator whilst running the algorithms
 * (anytime measure)
 */
inline void pls(std::string const &instance, std::size_t maxeval, unsigned int seed,
                apmnkl::pls::pac const pac, apmnkl::pls::pne const pne, std::ostream &os,
                apmnkl::objective_vector &ref) {
  if (ref.empty()) {
    apmnkl::pls pls(instance, seed);
    pls.run(maxeval, pac, pne);
    to_csv(os, "evaluation,hypervolume", pls.anytime());
  } else {
    apmnkl::pls pls(instance, seed, ref);
    pls.run(maxeval, pac, pne);
    to_csv(os, "evaluation,hypervolume", pls.anytime());
  }
}

/**
 * @brief CLI::App callback for the ibea algorithm
 *
 * @param instance The path for the "rmnk" instance (.dat) file to be used.
 *        These files can be generated using the rmnkGenerator.R script
 * @param maxeval The maximum number of evaluations performed by the algorithms (stopping criterion)
 * @param seed The seed used by the pseudo random number generator used in these algorithms
 * @param ps The maximum population size
 * @param gen The maximum number of generations
 * @param k The scaling factor
 * @param mp The mutation probability considered by the mutation operator
 * @param cp The crossover probability considered by the crossover operator
 * @param npts The number of crossover points considered by the NPointCrossover operator if used
 * @param mps The matting pool size
 * @param ts The size of the tournament considered by the KWayTournament operator
 * @param indicator The indicator to be used by the IBEA indicator operator
 * @param crossover The crossover method considered by the IBEA crossover operator
 * @param mutation The mutation method considered by the IBEA mutation operator
 * @param selection The selection method considered by the IBEA selection operator
 * @param adaptive boolean indicative of version of IBEA to be used.
 *                   If true use adaptive version of (A-IBEA) else use (B-IBEA)
 * @param os The name of the output file where the standard output stream should be redirected
 * @param ref The reference point considered by hypervolume indicator whilst running the algorithms
 * (anytime measure)
 */
inline void ibea(std::string const &instance, std::size_t const maxeval, unsigned int const seed,
                 std::size_t const ps, std::size_t const gen, double const k, double const mp,
                 double const cp, std::size_t npts, std::size_t const mps, std::size_t const ts,
                 indicator const indicator, crossover const crossover, mutation const mutation,
                 selection const selection, bool adaptive, std::ostream &os,
                 apmnkl::objective_vector &ref) {
  std::random_device dev;
  std::mt19937 rng(dev());

/**
 * @brief Helper define to avoid the use of runtime polymorphism methods to distinguish
 * between selection operators that ibea is going to use during its execution
 */
#define SELECTION(MAXEVAL, POP, GEN, FACTOR, I, C, M, S, ADAPT)                                    \
  switch (S) {                                                                                     \
    case selection::kwt:                                                                           \
      if (ref.empty()) {                                                                           \
        apmnkl::ibea ibea(instance, seed);                                                         \
        ibea.run(MAXEVAL, POP, GEN, FACTOR, I, C, M, apmnkl::selection::kwt(ts, mps, rng), ADAPT); \
        to_csv(os, "evaluation,generation,hypervolume", ibea.anytime());                           \
      } else {                                                                                     \
        apmnkl::ibea ibea(instance, seed, ref);                                                    \
        ibea.run(MAXEVAL, POP, GEN, FACTOR, I, C, M, apmnkl::selection::kwt(ts, mps, rng), ADAPT); \
        to_csv(os, "evaluation,generation,hypervolume", ibea.anytime());                           \
      }                                                                                            \
      break;                                                                                       \
    default:                                                                                       \
      throw("Unknown selection method!");                                                          \
  }

/**
 * @brief Helper define to avoid the use of runtime polymorphism methods to distinguish
 * between crossover operators that ibea is going to use during its execution
 *
 */
#define MUTATION(MAXEVAL, POP, GEN, FACTOR, I, C, M, S, ADAPT)                            \
  switch (M) {                                                                            \
    case mutation::um:                                                                    \
      SELECTION(MAXEVAL, POP, GEN, FACTOR, I, C, apmnkl::mutation::um(mp, rng), S, ADAPT) \
      break;                                                                              \
    default:                                                                              \
      throw("Unknown mutation operator!\n");                                              \
  }

/**
 * @brief Helper define to avoid the use of runtime polymorphism methods to distinguish
 * between crossover operators that ibea is going to use during its execution.
 *
 */
#define CROSSOVER(MAXEVAL, POP, GEN, FACTOR, I, C, M, S, ADAPT)                                  \
  switch (C) {                                                                                   \
    case crossover::npc:                                                                         \
      MUTATION(MAXEVAL, POP, GEN, FACTOR, I, apmnkl::crossover::npc(npts, cp, rng), M, S, ADAPT) \
      break;                                                                                     \
    case crossover::uc:                                                                          \
      MUTATION(MAXEVAL, POP, GEN, FACTOR, I, apmnkl::crossover::uc(cp, rng), M, S, ADAPT)        \
      break;                                                                                     \
    default:                                                                                     \
      throw("Unknown crossover operator!\n");                                                    \
  }

/**\
 * @brief Helper define to avoid the use of runtime polymorphism methods to distinguish
 * between indicator operators that ibea is going to use during its execution.
 *
 */
#define INDICATOR(MAXEVAL, POP, GEN, FACTOR, I, C, M, S, ADAPT)                                 \
  switch (I) {                                                                                  \
    case indicator::eps:                                                                        \
      CROSSOVER(MAXEVAL, POP, GEN, FACTOR, apmnkl::indicator::eps(), C, M, S, ADAPT)            \
      break;                                                                                    \
    case indicator::ihd:                                                                        \
      CROSSOVER(MAXEVAL, POP, GEN, FACTOR,                                                      \
                apmnkl::indicator::ihd(apmnkl::objective_vector(ibea.eval.getM(), 0)), C, M, S, \
                ADAPT)                                                                          \
      break;                                                                                    \
    default:                                                                                    \
      throw("Unknown indicator!\n");                                                            \
  }

/// Helper define that forwards all the information to the others
#define RUN_IBEA_LOOP(...) INDICATOR(__VA_ARGS__)
  RUN_IBEA_LOOP(maxeval, ps, gen, k, indicator, crossover, mutation, selection, adaptive);
}

int main(int argc, char **argv) {
  // App Global Settings
  CLI::App app(
      "Driver app to test the implementation of some search heuristics and "
      "gather\ndata relevant to the study of their performance from an anytime "
      "perspective\nin the context of pmnk-landscapes problem.\n",
      "anytime-pmnk-landscapes");

  // Setup a custom formatter for this app
  auto fmt = std::make_shared<CLI::Formatter>();
  fmt->column_width(40);
  fmt->label("SUBCOMMAND", "ALGORITHM");
  app.formatter(fmt);

  // Positional Arguments
  std::string instance;
  set_positional_arguments(app, instance);

  // General Settings
  std::size_t maxeval;
  unsigned int seed = std::random_device()();
  std::string outfile;
  apmnkl::objective_vector ref;
  set_general_options(app, maxeval, seed, outfile, ref);

  // App Parse Complete Callback (DEBUG)
  app.parse_complete_callback([&]() {
    // Required
    std::cerr << "Instance: " << instance << "\n";
    std::cerr << "Maxeval: " << maxeval << "\n";
    std::cerr << "Seed: " << seed << "\n";

    // Optionals
    if (!outfile.empty()) {
      std::cerr << "Output File:  " << outfile << "\n";
    }
    if (!ref.empty()) {
      std::cerr << "HV Reference: (" << ref[0];
      for (std::size_t i = 1; i < ref.size(); ++i) {
        std::cerr << ", " << ref[i];
      }
      std::cerr << ")\n";
    }
  });

  // Subcommands (Algorithms)
  app.require_subcommand(1);

  // GSEMO Subcommand
  auto gsemo_subcommand =
      app.add_subcommand("GSEMO",
                         "Run the global simple evolutionary multiobjective optimizer "
                         "algorithm\non the instance.")
          ->ignore_case()
          ->group("Algorithms");

  // GSEMO DEBUG
  gsemo_subcommand->callback([&]() { std::cerr << "Algorithm: GSEMO\n"; });

  // PLS Subcommand
  auto pls_subcommand =
      app.add_subcommand("PLS", "Run the pareto local search algorithm on the instance.")
          ->ignore_case()
          ->group("Algorithms");

  // PLS Subcommand Options
  auto pls_acceptance_criterion = apmnkl::pls::pac::non_dominating;
  auto pls_neighborhood_exploration = apmnkl::pls::pne::best_improvement;
  set_pls_options(*pls_subcommand, pls_acceptance_criterion, pls_neighborhood_exploration);

  // PLS Callback (DEBUG)
  pls_subcommand->callback([&pls_acceptance_criterion, &pls_neighborhood_exploration]() {
    std::cerr << "Algorithm: PLS\n";
    std::cerr << "Acceptance Criterion: " << static_cast<int>(pls_acceptance_criterion) << "\n";
    std::cerr << "Neighboorhood Exploration: " << static_cast<int>(pls_neighborhood_exploration)
              << "\n";
  });

  // IBEA Subcommand
  auto ibea_subcommand =
      app.add_subcommand("IBEA", "Run the indicator-based evolutionary algorithm on the instance.")
          ->ignore_case()
          ->group("Algorithms");

  // Setup a custom formatter for the this ibea_subcommand
  auto ibea_fmt = std::make_shared<CLI::Formatter>();
  ibea_fmt->column_width(40);
  ibea_fmt->label("SUBCOMMAND", "INDICATOR MUTATION CROSSOVER SELECTION");
  ibea_subcommand->formatter(ibea_fmt);

  // IBEA Subcommand Options
  std::size_t pop;
  std::size_t gen;
  double factor;
  bool adaptive = false;
  set_ibea_options(*ibea_subcommand, pop, gen, factor, adaptive);

  // IBEA Subcommands
  ibea_subcommand->require_subcommand(4);

  // IBEA Indicators (IBEA Sub-Subcommand)
  auto ihd = ibea_subcommand->add_subcommand("IHD", "Run using the hypervolume indicator")
                 ->ignore_case()
                 ->group("Indicators");

  // IHD Callback (DEBUG)
  ihd->callback([]() { std::cerr << "Indicator: IHD\n"; });

  auto eps = ibea_subcommand->add_subcommand("EPS", "Run using the epsilon (+) indicator")
                 ->ignore_case()
                 ->excludes(ihd)
                 ->group("Indicators");

  // EPS(+) Callback (DEBUG)
  eps->callback([]() { std::cerr << "Indicator: EPS(+)\n"; });

  // IBEA Mutation Operators (IBEA Sub-Subcommand)
  auto um =
      ibea_subcommand->add_subcommand("UniformMutation", "Run using a uniform mutation operator")
          ->ignore_case()
          ->alias("UM")
          ->group("Mutation Operators");

  // IBEA Mutation Operator Options
  double mutation_probability;
  set_ibea_mutation_options(*um, mutation_probability);

  // Uniform Mutation Callback (DEBUG)
  um->callback([&]() {
    std::cerr << "Mutation Operator: UniformMutation\n";
    std::cerr << "Mutation probability: " << mutation_probability << "\n";
  });

  // IBEA Crossover Operators (IBEA Sub-Subcommand)
  auto npc =
      ibea_subcommand->add_subcommand("NPointCrossover", "Run using a n-point crossover operator")
          ->ignore_case()
          ->alias("NPC")
          ->group("Crossover Operators");

  auto uc =
      ibea_subcommand->add_subcommand("UniformCrossover", "Run using a uniform crossover operator")
          ->ignore_case()
          ->excludes(npc)
          ->alias("UC")
          ->group("Crossover Operators");

  // IBEA Crossover Operator Options
  std::size_t npoints;
  double crossover_probability;
  set_ibea_crossover_options(*npc, crossover_probability, npoints);
  set_ibea_crossover_options(*uc, crossover_probability, npoints);

  // N-Point Crossover Callback (DEBUG)
  npc->callback([&]() {
    std::cerr << "Crossover Operator: N-Point Crossover\n";
    std::cerr << "Crossover Probability: " << crossover_probability << "\n";
    std::cerr << "Number of Crossover Points: " << npoints << "\n";
  });

  // Uniform Crossover Callback (DEBUG)
  uc->callback([&]() {
    std::cerr << "Crossover Operator: Uniform Crossover\n";
    std::cerr << "Crossover Probability: " << crossover_probability << "\n";
  });

  // IBEA Selection Operators (IBEA Sub-Subcommand)
  auto kwt =
      ibea_subcommand
          ->add_subcommand("KWayTournament", "Run using a k-way tournament selection operator")
          ->ignore_case()
          ->alias("KWT")
          ->group("Selection Operators");

  // IBEA Selection Operator Options
  std::size_t matting_pool_size;
  std::size_t tournament_size;
  set_ibea_selection_options(*kwt, matting_pool_size, tournament_size);

  // K-Way Tournament Selection Callback (DEBUG)
  kwt->callback([&]() {
    std::cerr << "Selection Operator: K-Way Tournament\n";
    std::cerr << "Matting Pool Size: " << matting_pool_size << "\n";
    std::cerr << "Tournament Size: " << tournament_size << "\n";
  });

  // IBEA Callback (DEBUG)
  ibea_subcommand->callback([&]() {
    std::cerr << "Algorithm: IBEA\n";
    std::cerr << "Population Size: " << pop << "\n";
    std::cerr << "Generations: " << gen << "\n";
    std::cerr << "Scaling Factor: " << factor << "\n";
    std::cerr << "Adaptive: " << std::boolalpha << adaptive << "\n";
  });

  // Main App
  app.callback([&]() {
    std::ofstream of;
    auto buf = !outfile.empty() ? (of.open(outfile), of.rdbuf()) : std::cout.rdbuf();
    std::ostream os(buf);

    if (app.got_subcommand("GSEMO")) {
      gsemo(instance, maxeval, seed, os, ref);

    } else if (app.got_subcommand("PLS")) {
      pls(instance, maxeval, seed, pls_acceptance_criterion, pls_neighborhood_exploration, os, ref);

    } else if (app.got_subcommand("IBEA")) {
      indicator ind = ibea_subcommand->got_subcommand("IHD") ? indicator::ihd : indicator::eps;
      mutation mut =
          ibea_subcommand->got_subcommand("UM") ? mutation::um : static_cast<mutation>(-1);
      crossover cross = ibea_subcommand->got_subcommand("NPC") ? crossover::npc : crossover::uc;
      selection sel =
          ibea_subcommand->got_subcommand("KWT") ? selection::kwt : static_cast<selection>(-1);
      ibea(instance, maxeval, seed, pop, gen, factor, mutation_probability, crossover_probability,
           npoints, matting_pool_size, tournament_size, ind, cross, mut, sel, adaptive, os, ref);
    }
  });

  CLI11_PARSE(app, argc, argv);
  return EXIT_SUCCESS;
}
