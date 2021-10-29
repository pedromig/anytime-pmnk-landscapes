/**
 * @file pls.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief PLS (Pareto Local Search) algorithm implementation
 * @version 0.1.0
 * @date 13-09-2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef PLS_HPP
#define PLS_HPP

#include <csignal>
#include <iomanip>

#include "utils/solution.hpp"
#include "utils/utils.hpp"
#include "utils/wfg.hpp"

namespace apmnkl {

/// Wrapper class for the PLS algorithm
class pls {
  using objv_type = typename apmnkl::objective_vector;
  using decv_type = typename apmnkl::decision_vector;
  using hv_type = typename objv_type::value_type;
  using solution_type = typename priv::solution;

 public:
  priv::RMNKEval eval;

 private:
  std::mt19937 m_generator;

  priv::hvobj<hv_type> m_hvo;
  std::vector<std::tuple<std::size_t, hv_type>> m_anytime;

  std::vector<solution_type> m_solutions;
  std::vector<solution_type> m_non_visited_solutions;

 public:
  /** Acceptance criterion:
   *  - 0 -> accept every non-dominated neighbor (non_dominating).
   *  - 1 -> accept only neighbors that dominate current solution (dominating).
   *  - 2 -> first try to accept only neighbors that dominate the (both)
   *         current solution, if none exist accept non-dominated.
   */
  enum class pls_acceptance_criterion { non_dominating, dominating, both };

  /// Helper `using` to avoid typing too much
  using pac = pls_acceptance_criterion;

  /** Neighborhood exploration:
   *  - 0 -> explore every acceptable neighbor (best_improvement).
   *  - 1 -> stop once one neighbor is accepted (first_improvement).
   *  - 2 -> use 1 until PLS stops, afterwards restart and use 0 (best_improvement).
   */
  enum class pls_neighborhood_exploration { best_improvement, first_improvement, both };

  /// Helper `using` to avoid typing too much
  using pne = pls_neighborhood_exploration;

  /**
   * @brief Construct a new pls object.
   *
   * @tparam Str the type used to store the instance path.
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script.
   * @param seed The seed used by the pseudo random number generator used in PLS
   * @param ref The reference point considered by hypervolume indicator whilst running the
   * algorithms (anytime measure)
   */

  template <typename Str = std::string, typename Ref = objv_type>
  pls(Str &&instance, unsigned int seed, Ref &&ref)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new pls object.
   *
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in PLS
   */
  template <typename Str = std::string>
  pls(Str &&instance, unsigned int seed)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(objv_type(eval.getM(), 0.0)) {}

  /**
   * @brief Construct a new pls object
   *
   * @tparam Str the type used to store the instance path.
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script.
   * @param ref The reference point considered by hypervolume indicator whilst running the
   * algorithms (anytime measure).
   */
  template <typename Str = std::string, typename Ref = objv_type>
  pls(Str &&instance, Ref &&ref)
      : pls(std::forward<Str>(instance), std::random_device()(), std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new pls object
   *
   * @tparam Str the type used to store the instance path.
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script.
   */
  template <typename Str = std::string>
  pls(Str &&instance)
      : pls(std::forward<Str>(instance), std::random_device()()) {}

  /**
   * @brief Getter for the vector of solutions found by this algorithm.
   *
   * @return auto const& Read-Only reference to a vector of
   *         visited solutions produced by the PLS algorithm.
   */
  auto const &solutions() const {
    return m_solutions;
  }

  /**
   * @brief Getter for the vector of solutions found by this algorithm
   *        that were not visited in the process of local search.
   *
   * @return auto const& Read-Only reference to a vector of
   *         non visited solutions produced by PLS algorithm.
   */
  auto const &non_visited_solutions() const {
    return m_non_visited_solutions;
  }

  /**
   * @brief Getter for the anytime data produced by this algorithm.
   *
   * @return auto const& Read-Only reference to a vector of pairs <evaluation,
   *         hypervolume> obtained in the run of the PLS algorithm.
   */
  auto const &anytime() const {
    return m_anytime;
  }

  /**
   * @brief PLS implementation runner. This effectively starts the algorithm and runs it
   * until the maximum number of evaluations has been reached.
   *
   * @param maxeval The maximum number of evaluations performed by PLS (stopping criterion)
   * @param acceptance_criterion The PLS algorithm solution acceptance criterion
   * @param neighborhood_exploration The PLS algorithm solution exploration criterion
   */
  void run(std::size_t maxeval, pac const acceptance_criterion,
           pne const neighborhood_exploration) {
    auto rand_solution = solution_type::random_solution(eval, m_generator);
    m_hvo.insert(rand_solution.objective_vector());

    add_non_dominated(m_non_visited_solutions, std::move(rand_solution));
    m_solutions = m_non_visited_solutions;

    std::size_t evaluation = 0;
    m_anytime.push_back({evaluation, m_hvo.value()});

#define RUNLOOP(FIRSTIMPROV)                                         \
  switch (acceptance_criterion) {                                    \
    case pac::non_dominating:                                        \
      m_loop<FIRSTIMPROV, pac::non_dominating>(evaluation, maxeval); \
      break;                                                         \
    case pac::dominating:                                            \
      m_loop<FIRSTIMPROV, pac::dominating>(evaluation, maxeval);     \
      break;                                                         \
    case pac::both:                                                  \
      m_loop<FIRSTIMPROV, pac::both>(evaluation, maxeval);           \
      break;                                                         \
    default:                                                         \
      throw("Unknown value for acceptance criterion");               \
  }

    if (neighborhood_exploration == pne::best_improvement) {
      RUNLOOP(false);
    } else if (neighborhood_exploration == pne::first_improvement) {
      RUNLOOP(true);
    } else if (neighborhood_exploration == pne::both) {
      RUNLOOP(true);
      RUNLOOP(false);
    } else {
      throw("Unknown value for neighborhood exploration");
    }
  }

 private:
  /**
   * @brief Helper function that provides the implementation of
   *        multiple acceptance/exploration criterion exploration methods.
   *
   * @tparam FirstImprov Boolean template parameter indicating if the first
   *         improvement technique should be used
   * @tparam Acceptance Neighborhood exporation criterion type indication the
   *         exploration method to be used.
   * @param evaluation The current evaluation number.
   * @param maxeval The maximum number of evaluations.
   */
  template <bool FirstImprov, pac Acceptance>
  void m_loop(size_t &evaluation, size_t maxeval) {
    while (evaluation < maxeval && !m_non_visited_solutions.empty()) {
      std::uniform_int_distribution<std::size_t> distrib(0, m_non_visited_solutions.size() - 1);
      std::size_t index = distrib(m_generator);

      auto original = std::move(m_non_visited_solutions[index]);
      m_non_visited_solutions[index] = std::move(m_non_visited_solutions.back());
      m_non_visited_solutions.pop_back();

      if constexpr (Acceptance == pac::non_dominating) {
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          auto decv = original.decision_vector();
          decv[i] = !decv[i];
          auto solution = solution_type(eval, std::move(decv));
          ++evaluation;
          if (add_non_dominated(m_solutions, solution)) {
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_anytime.push_back({evaluation, m_hvo.value()});
            if constexpr (FirstImprov) {
              break;
            }
          }
        }
      } else if constexpr (Acceptance == pac::dominating) {
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          auto decv = original.decision_vector();
          decv[i] = !decv[i];
          auto solution = solution_type(eval, std::move(decv));
          ++evaluation;
          if (solution.dominance(original) == priv::dominance_type::dominates &&
              add_non_dominated(m_solutions, solution)) {
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_anytime.push_back({evaluation, m_hvo.value()});
            if constexpr (FirstImprov) {
              break;
            }
          }
        }
      } else if constexpr (Acceptance == pac::both) {
        std::vector<std::pair<solution_type, size_t>> remaining;
        remaining.reserve(original.decision_vector().size());
        bool use_remaining = true;
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          auto decv = original.decision_vector();
          decv[i] = !decv[i];
          auto solution = solution_type(eval, std::move(decv));
          ++evaluation;
          if (solution.dominance(original) == priv::dominance_type::dominates &&
              add_non_dominated(m_solutions, solution)) {
            use_remaining = false;
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_anytime.push_back({evaluation, m_hvo.value()});
            if constexpr (FirstImprov) {
              break;
            }
          } else if (use_remaining) {
            remaining.emplace_back(std::move(solution), evaluation);
          }
        }
        if (use_remaining) {
          for (auto &&[solution, iteration] : remaining) {
            if (add_non_dominated(m_solutions, solution)) {
              m_hvo.insert(solution.objective_vector());
              add_non_dominated(m_non_visited_solutions, std::move(solution));
              m_anytime.push_back({evaluation, m_hvo.value()});
              if constexpr (FirstImprov) {
                break;
              }
            }
          }
        }
      }
    }
  }
};
}  // namespace apmnkl
#endif  // PLS_HPP
