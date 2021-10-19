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

#include "../Utils/solution.hpp"
#include "../Utils/utils.hpp"
#include "../Utils/wfg.hpp"

namespace pmnk {

// Acceptance criterion:
//   - 0 -> accept every non-dominated neighbor (NON_DOMINATING)
//   - 1 -> accept only neighbors that dominate current solution (DOMINATING)
//   - 2 -> first try to accept only neighbors that dominate the (BOTH)
//          current solution, if none exist accept non-dominated

enum class PLSAcceptanceCriterion { NON_DOMINATING, DOMINATING, BOTH };

// Neighborhood exploration:
//   - 0 -> explore every acceptable neighbor (BEST_IMPROVEMENT)
//   - 1 -> stop once one neighbor is accepted (FIRST_IMPROVEMENT)
//   - 2 -> use 1 until PLS stops, afterwards restart and use 0 (BEST_IMPROVEMENT)
enum class PLSExplorationCriterion { BEST_IMPROVEMENT, FIRST_IMPROVEMENT, BOTH };

/// Wrapper class for PLS
class PLS {
 public:
  RMNKEval eval;

 private:
  std::mt19937 m_generator;
  std::ostream &m_os;

  hvobj<typename ObjectiveVector::value_type> m_hvo;
  std::vector<Solution> m_solutions;
  std::vector<Solution> m_non_visited_solutions;

 public:

  /**
  * @brief Construct a new PLS object
  * 
  * @tparam Str the type used to store the instance path
  * @tparam Ref the type used to store the reference point of the hvobj obj 
  * @param instance The path for the "rmnk" instance (.dat) file to be used. 
  *        These files can be generated using the rmnkGenerator.R script 
  * @param seed The seed used by the pseudo random number generator used in PLS 
  * @param os The name of the output file where the standard output stream should be redirected 
  * @param ref The reference point considered by hypervolume indicator whilst running the algorithms (anytime measure)
  */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  PLS(Str &&instance, unsigned int seed, std::ostream &os, Ref &&ref)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_os(os)
      , m_hvo(std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new PLS object
   * 
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used. 
   *        These files can be generated using the rmnkGenerator.R script 
   * @param seed The seed used by the pseudo random number generator used in PLS 
   * @param os The name of the output file where the standard output stream should be redirected 
   */
  template <typename Str = std::string>
  PLS(Str &&instance, unsigned int seed, std::ostream &os)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_os(os)
      , m_hvo(ObjectiveVector(eval.getM(), 0.0)) {}

  /**
  * @brief Construct a new PLS object
  * 
  * @tparam Str the type used to store the instance path
  * @tparam Ref the type used to store the reference point of the hvobj obj 
  * @param instance The path for the "rmnk" instance (.dat) file to be used. 
  *        These files can be generated using the rmnkGenerator.R script 
  * @param seed The seed used by the pseudo random number generator used in PLS 
  * @param ref The reference point considered by hypervolume indicator whilst running the algorithms (anytime measure)
  */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  PLS(Str &&instance, unsigned int seed, Ref &&ref)
      : PLS(std::forward<Str>(instance), seed, std::cout, std::forward<Ref>(ref)) {}

  /**
  * @brief Construct a new PLS object
  * 
  * @tparam Str the type used to store the instance path
  * @param instance The path for the "rmnk" instance (.dat) file to be used. 
  *        These files can be generated using the rmnkGenerator.R script 
  * @param seed The seed used by the pseudo random number generator used in PLS 
  */
  template <typename Str = std::string>
  PLS(Str &&instance, unsigned int seed)
      : PLS(std::forward<Str>(instance), seed, std::cout) {}

  /**
  * @brief Construct a new PLS object
  * 
  * @tparam Str the type used to store the instance path
  * @tparam Ref the type used to store the reference point of the hvobj obj 
  * @param instance The path for the "rmnk" instance (.dat) file to be used. 
  *        These files can be generated using the rmnkGenerator.R script 
  * @param ref The reference point considered by hypervolume indicator whilst running the algorithms (anytime measure)
  */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  PLS(Str &&instance, Ref &&ref)
      : PLS(std::forward<Str>(instance), std::random_device()(), std::cout,
            std::forward<Ref>(ref)) {}

  /**
  * @brief Construct a new PLS object
  * 
  * @tparam Str the type used to store the instance path
  * @param instance The path for the "rmnk" instance (.dat) file to be used. 
  *        These files can be generated using the rmnkGenerator.R script 
  */
  template <typename Str = std::string>
  PLS(Str &&instance)
      : PLS(std::forward<Str>(instance), std::random_device()(), std::cout) {}

  /**
   * @brief Getter for the vector of solutions found by this algorithm.
   * 
   * @return std::vector<Solution> const& Read-Only reference to a vector of 
   *         visited solutions produced by PLS.
   */
  std::vector<Solution> const solutions() const {
    return m_solutions;
  }

  /**
   * @brief Getter for the vector of solutions found by this algorithm
   *        that were not visited in the process of local search.
   * 
   * @return std::vector<Solution> const& Read-Only reference to a vector of 
   *         non visited solutions produced by PLS.
   */
  std::vector<Solution> const non_visited_solutions() const {
    return m_non_visited_solutions;
  }

  /**
   * @brief PLS implementation runner. This effectively starts the algorithm and runs it
   * until the maximum number of evaluations has been reached.
   * 
   * @param maxeval The maximum number of evaluations performed by PLS (stopping criterion)
   * @param acceptance_criterion The PLS algorithm solution acceptance criterion
   * @param neighborhood_exploration The PLS algorithm solution exploration criterion
   */
  void run(std::size_t maxeval, PLSAcceptanceCriterion const acceptance_criterion,
           PLSExplorationCriterion const neighborhood_exploration) {
    auto rand_solution = Solution::random_solution(eval, m_generator);
    m_hvo.insert(rand_solution.objective_vector());

    add_non_dominated(m_non_visited_solutions, std::move(rand_solution));
    m_solutions = m_non_visited_solutions;

    std::size_t eval = 0;
    m_os << "evaluation,hypervolume\n";
    m_os << std::setprecision(12) << eval << "," << m_hvo.value() << "\n";

#define RUNLOOP(FIRSTIMPROV)                                                      \
  switch (acceptance_criterion) {                                                 \
    case PLSAcceptanceCriterion::NON_DOMINATING:                                  \
      m_loop<FIRSTIMPROV, PLSAcceptanceCriterion::NON_DOMINATING>(eval, maxeval); \
      break;                                                                      \
    case PLSAcceptanceCriterion::DOMINATING:                                      \
      m_loop<FIRSTIMPROV, PLSAcceptanceCriterion::DOMINATING>(eval, maxeval);     \
      break;                                                                      \
    case PLSAcceptanceCriterion::BOTH:                                            \
      m_loop<FIRSTIMPROV, PLSAcceptanceCriterion::BOTH>(eval, maxeval);           \
      break;                                                                      \
    default:                                                                      \
      throw("Unknown value for acceptance criterion");                            \
  }

    if (neighborhood_exploration == PLSExplorationCriterion::BEST_IMPROVEMENT) {
      RUNLOOP(false);
    } else if (neighborhood_exploration == PLSExplorationCriterion::FIRST_IMPROVEMENT) {
      RUNLOOP(true);
    } else if (neighborhood_exploration == PLSExplorationCriterion::BOTH) {
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
  template <bool FirstImprov, PLSAcceptanceCriterion Acceptance>
  void m_loop(size_t &evaluation, size_t maxeval) {
    while (evaluation < maxeval && !m_non_visited_solutions.empty()) {
      std::uniform_int_distribution<std::size_t> distrib(0, m_non_visited_solutions.size() - 1);
      std::size_t index = distrib(m_generator);

      auto original = std::move(m_non_visited_solutions[index]);
      m_non_visited_solutions[index] = std::move(m_non_visited_solutions.back());
      m_non_visited_solutions.pop_back();

      if constexpr (Acceptance == PLSAcceptanceCriterion::NON_DOMINATING) {
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          DecisionVector decision_vector = original.decision_vector();
          decision_vector[i] = !decision_vector[i];
          auto solution = Solution(eval, std::move(decision_vector));
          ++evaluation;
          if (add_non_dominated(m_solutions, solution)) {
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_os << std::setprecision(12) << evaluation << "," << m_hvo.value() << "\n";
            if constexpr (FirstImprov) {
              break;
            }
          }
        }
      } else if constexpr (Acceptance == PLSAcceptanceCriterion::DOMINATING) {
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          DecisionVector decision_vector = original.decision_vector();
          decision_vector[i] = !decision_vector[i];
          auto solution = Solution(eval, std::move(decision_vector));
          ++evaluation;
          if (solution.dominance(original) == DominanceType::DOMINATES &&
              add_non_dominated(m_solutions, solution)) {
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_os << std::setprecision(12) << evaluation << "," << m_hvo.value() << "\n";
            if constexpr (FirstImprov) {
              break;
            }
          }
        }
      } else if constexpr (Acceptance == PLSAcceptanceCriterion::BOTH) {
        std::vector<std::pair<Solution, size_t>> remaining;
        remaining.reserve(original.decision_vector().size());
        bool use_remaining = true;
        for (size_t i = 0; i < original.decision_vector().size() && evaluation < maxeval; ++i) {
          DecisionVector decision_vector = original.decision_vector();
          decision_vector[i] = !decision_vector[i];
          auto solution = Solution(eval, std::move(decision_vector));
          ++evaluation;
          if (solution.dominance(original) == DominanceType::DOMINATES &&
              add_non_dominated(m_solutions, solution)) {
            use_remaining = false;
            m_hvo.insert(solution.objective_vector());
            add_non_dominated(m_non_visited_solutions, std::move(solution));
            m_os << std::setprecision(12) << evaluation << "," << m_hvo.value() << "\n";
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
              m_os << std::setprecision(12) << iteration << "," << m_hvo.value() << "\n";
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
}  // namespace pmnk
#endif  // PLS_HPP
