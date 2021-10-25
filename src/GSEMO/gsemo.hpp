/**
 * @file gsemo.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief GSEMO (Global Simple Multi-objective Optimizer) algorithm
 * implementation
 * @version 0.1.0
 * @date 13-09-2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef GSEMO_HPP
#define GSEMO_HPP

#include <algorithm>
#include <csignal>
#include <iomanip>
#include <random>

#include "../Utils/solution.hpp"
#include "../Utils/utils.hpp"
#include "../Utils/wfg.hpp"

namespace pmnk {

/// Wrapper class for GSEMO.
class GSEMO {
  using hv_t = typename ObjectiveVector::value_type;

 public:
  RMNKEval eval;

 private:
  std::mt19937 m_generator;

  hvobj<hv_t> m_hvo;
  std::vector<Solution> m_solutions;
  std::vector<std::tuple<std::size_t, hv_t>> m_anytime;

 public:
  /**
   * @brief Construct a new GSEMO object
   *
   * @tparam Str the type used to store the instance path
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in
   * GSEMO
   * should be redirected
   * @param ref The reference point considered by hypervolume indicator whilst
   * running the algorithms (anytime measure)
   */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  GSEMO(Str &&instance, unsigned int const seed, Ref &&ref)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new GSEMO object
   *
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in
   * GSEMO
   */
  template <typename Str = std::string>
  GSEMO(Str &&instance, unsigned int const seed)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(ObjectiveVector(eval.getM(), 0)) {}

  /**
   * @brief Construct a new GSEMO object
   *
   * @tparam Str the type used to store the instance path
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param ref The reference point considered by hypervolume indicator whilst
   * running the algorithms (anytime measure)
   */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  explicit GSEMO(Str &&instance, Ref &&ref)
      : GSEMO(std::forward<Str>(instance), std::random_device()(), std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new GSEMO object
   *
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   */
  template <typename Str = std::string>
  explicit GSEMO(Str &&instance)
      : GSEMO(std::forward<Str>(instance), std::random_device()()) {}

  /**
   * @brief Getter for the vector of solutions found by this algorithm.
   *
   * @return std::vector<Solution> const& Read-Only reference to a vector of
   * solutions found by the GSEMO algorithm.
   */
  std::vector<Solution> const &solutions() const {
    return m_solutions;
  }

  /**
   * @brief Getter for the anytime data produced by this algorithm.
   *
   * @return std::vecto<std::tuple<std::size_t, hv_t>> const&
   *         Read-Only reference to a vector of pairs <evaluation,
   *         hypervolume> obtained in the run of the GSEMO algorithm.
   */
  auto const &anytime() const {
    return m_anytime;
  }

  /**
   * @brief GSEMO implementation runner. This effectively starts the algorithm
   * and runs it until the maximum number of evaluations has been reached.
   *
   * @param maxeval The maximum number of evaluations performed by GSEMO
   * (stopping criterion)
   */
  void run(std::size_t maxeval) {
    auto rand_solution = Solution::random_solution(eval, m_generator);
    m_hvo.insert(rand_solution.objective_vector());

    add_non_dominated(m_solutions, std::move(rand_solution));
    m_anytime.push_back({0, m_hvo.value()});

    for (std::size_t i = 0; i < maxeval; ++i) {
      std::uniform_int_distribution<std::size_t> randint(0, m_solutions.size() - 1);

      std::size_t index = randint(m_generator);
      auto solution = Solution::uniform_bit_flip_solution(eval, m_generator, m_solutions[index]);

      auto sov = solution.objective_vector();
      if (add_non_dominated(m_solutions, std::move(solution))) {
        m_hvo.insert(sov);
        m_anytime.push_back({i + 1, m_hvo.value()});
      }
    }
  }
};
}  // namespace pmnk
#endif  // GSEMO_HPP
