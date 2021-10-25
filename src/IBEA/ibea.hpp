/**
 * @file ibea.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief IBEA (Indicator Based Evolutionary Algorithm) implementation
 * @version 0.1.0
 * @date 13-09-2021
 *
 * @copyright Copyright (c) 2021
 */

#ifndef IBEA_HPP
#define IBEA_HPP

#include <csignal>
#include <iomanip>
#include <random>

#include "../Utils/solution.hpp"
#include "../Utils/utils.hpp"
#include "../Utils/wfg.hpp"
#include "functor.hpp"

namespace pmnk {

/// Wrapper class for IBEA
class IBEA {
  using hv_t = typename ObjectiveVector::value_type;

 public:
  RMNKEval eval;

 private:
  std::mt19937 m_generator;

  hvobj<ObjectiveVector::value_type> m_hvo;
  std::vector<std::tuple<std::size_t, std::size_t, hv_t>> m_anytime;

  std::vector<GASolution> m_solutions;

 public:
  /**
   * @brief Construct a new IBEA object
   *
   * @tparam Str the type used to store the instance path
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in IBEA
   * @param ref The reference point considered by hypervolume indicator whilst running the
   * algorithms (anytime measure)
   */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  IBEA(Str &&instance, unsigned int const seed, Ref &&ref)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new IBEA object
   *
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in IBEA
   */
  template <typename Str = std::string>
  IBEA(Str &&instance, unsigned int const seed)
      : eval(std::forward<Str>(instance).c_str())
      , m_generator(seed)
      , m_hvo(ObjectiveVector(eval.getM(), 0.0)) {}

  /**
   * @brief Construct a new IBEA object
   *
   * @tparam Str the type used to store the instance path
   * @tparam Ref the type used to store the reference point of the hvobj obj
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   * @param seed The seed used by the pseudo random number generator used in IBEA
   * @param ref The reference point considered by hypervolume indicator whilst running the
   * algorithms (anytime measure)
   */
  template <typename Str = std::string, typename Ref = ObjectiveVector>
  explicit IBEA(Str &&instance, Ref &&ref)
      : IBEA(std::forward<Str>(instance), std::random_device()(), std::forward<Ref>(ref)) {}

  /**
   * @brief Construct a new IBEA object
   *
   * @tparam Str the type used to store the instance path
   * @param instance The path for the "rmnk" instance (.dat) file to be used.
   *        These files can be generated using the rmnkGenerator.R script
   */
  template <typename Str = std::string>
  explicit IBEA(Str &&instance)
      : IBEA(std::forward<Str>(instance), std::random_device()()) {}

  /**
   * @brief Getter for the vector of solutions found by this algorithm.
   *
   * @return std::vector<Solution> const& Read-Only reference to a vector of solutions
   *         found by the IBEA.
   */
  std::vector<GASolution> const &solutions() const {
    return m_solutions;
  }

  /**
   * @brief Getter for the anytime data produced by this algorithm.
   *
   * @return std::vecto<std::tuple<std::size_t, hv_t>> const&
   *         Read-Only reference to a vector of pairs <evaluation,
   *         hypervolume> obtained in the run of IBEA.
   */
  auto const &anytime() const {
    return m_anytime;
  }

  /**
   * @brief IBEA implementation runner. This effectively starts the algorithm and runs it
   * until the maximum number of evaluations has been reached.
   *
   * @tparam I The type used to store an IBEA indicator
   * @tparam S The type used to store an IBEA selection operator
   * @tparam M The type used to store and IBEA mutation operator
   * @tparam C The type used to store and IBEA crossover operator
   * @param maxeval The maximum number of evaluations performed by the algorithms (stopping
   * criterion)
   * @param population_max_size The maximum population size
   * @param max_generations The maximum number of generations
   * @param scaling_factor The scaling factor
   * @param indicator  The indicator to be used by the IBEA indicator operator
   * @param crossover_method The crossover method considered by the IBEA mutation operator
   * @param mutation_method The mutation method considered by the IBEA mutation operator
   * @param selection_method The selection method considered by the IBEA selection operator
   * @param adaptive boolean indicative of version of IBEA to be used.
   *                   If true use adaptive version of (A-IBEA) else use (B-IBEA)
   */
  template <typename I, typename S, typename M, typename C>
  void run(std::size_t const maxeval, std::size_t const population_max_size,
           std::size_t const max_generations, double const scaling_factor, I &&indicator,
           C &&crossover_method, M &&mutation_method, S &&selection_method, bool adaptive) {
    if (adaptive) {
      m_run<true>(maxeval, population_max_size, max_generations, scaling_factor, indicator,
                  crossover_method, mutation_method, selection_method);
    } else {
      m_run<false>(maxeval, population_max_size, max_generations, scaling_factor, indicator,
                   crossover_method, mutation_method, selection_method);
    }
  }

 private:
  /**
   * @brief IBEA implementation runner. This effectively starts the algorithm and runs it
   * until the maximum number of evaluations has been reached.
   *
   * @tparam Adaptive boolean template type indicative of version of IBEA to be used.
   *                    If true use adaptive version of (A-IBEA) else use (B-IBEA)
   * @tparam I The type used to store an IBEA indicator
   * @tparam S The type used to store an IBEA selection operator
   * @tparam M The type used to store and IBEA mutation operator
   * @tparam C The type used to store and IBEA crossover operator
   * @param maxeval The maximum number of evaluations performed by the algorithms (stopping
   * criterion)
   * @param pop_max The maximum population size
   * @param max_generations The maximum number of generations
   * @param scaling_factor The scaling factor
   * @param indicator  The indicator to be used by the IBEA indicator operator
   * @param crossover_method The crossover method considered by the IBEA mutation operator
   * @param mutation_method The mutation method considered by the IBEA mutation operator
   * @param selection_method The selection method considered by the IBEA selection operator
   */
  template <bool Adaptive, typename I, typename S, typename M, typename C>
  void m_run(std::size_t const maxeval, std::size_t const pop_max,
             std::size_t const max_generations, double const scaling_factor, I &&indicator,
             C &&crossover_method, M &&mutation_method, S &&selection_method) {
    std::size_t evaluation = 0, gen = 0;
    double c = 1;

    std::vector<GASolution> population;
    population.reserve(pop_max);

    for (std::size_t i = 0; i < pop_max && evaluation < maxeval; ++i) {
      auto sol = GASolution(Solution::random_solution(eval, m_generator));
      if (add_non_dominated(m_solutions, sol)) {
        m_hvo.insert(sol.objective_vector());
        m_anytime.push_back({evaluation, gen, m_hvo.value()});
      }
      population.push_back(std::move(sol));
      ++evaluation;
    }

    if (evaluation < maxeval) {
      if constexpr (Adaptive) {
        c = m_adaptive_factor(population, indicator);
      }
      m_fitness_assignment(population, scaling_factor * c, indicator);
    }

    for (; evaluation < maxeval && gen < max_generations; ++gen) {
      auto matting_pool = selection_method(population);

      for (std::size_t i = 0; i < matting_pool.size() - 1; i += 2) {
        crossover_method(matting_pool[i], matting_pool[i + 1]);
      }

      for (auto &individual : matting_pool) {
        mutation_method(individual);
        individual.eval(eval);
      }

      if constexpr (Adaptive) {
        c = m_adaptive_factor(population, indicator);
      }
      m_fitness_assignment(population, scaling_factor * c, indicator);

      for (auto &individual : matting_pool) {
        if (add_non_dominated(m_solutions, individual)) {
          m_hvo.insert(individual.objective_vector());
          m_anytime.push_back({evaluation, gen, m_hvo.value()});
        }
        population.push_back(std::move(individual));
        ++evaluation;
      }
      m_environmental_selection(population, scaling_factor * c, pop_max, indicator);
    }
    m_anytime.push_back({evaluation, gen, m_hvo.value()});
  }

  /**
   * @brief Calculate objective functions lower and upper bounds
   *        for scaling
   *
   * @tparam S The type for the genetic algorithm's solution
   * @param population The IBEA population to be scaled
   * @return auto A pair containing the upper and lower objective value bounds
   */
  template <typename S = GASolution>
  auto m_objective_bounds(std::vector<S> const &population) {
    auto ub = std::numeric_limits<ObjectiveVector::value_type>::min();
    auto lb = std::numeric_limits<ObjectiveVector::value_type>::max();
    for (auto const &individual : population) {
      for (auto const v : individual.objective_vector()) {
        lb = std::min(lb, v);
        ub = std::max(ub, v);
      }
    }
    return std::make_pair(lb, ub);
  }

  /**
   * @brief Scale population objective vectors values using the bounds
   *        previously calculated for the population.
   *
   * @tparam S The type for the genetic algorithm's solution
   * @param population The IBEA population to be scaled
   * @param lb The population objective vectors lower bound
   * @param ub the population objective vectors upper bound
   * @return auto The population with the objective vectors scaled
   */
  template <typename S = GASolution>
  auto m_scale_objective_vectors(std::vector<S> const &population, double const lb,
                                 double const ub) {
    std::vector<S> s = population;
    for (auto &individual : s) {
      auto ov = individual.objective_vector();
      for (auto &i : ov) {
        i = (i - ub) / (ub - lb);
      }
      individual.set_objv(std::move(ov));
    }
    return s;
  }

  /**
   * @brief Calculate the adaptive factor for IBEA
   *
   * @tparam I The type for the IBEA indicator
   * @tparam S The type for the genetic algorithm's solution
   * @param population The IBEA population to be scaled
   * @param indicator The indicator used by IBEA
   * @return auto The adaptive factor.
   */
  template <typename I, typename S = GASolution>
  auto m_adaptive_factor(std::vector<S> const &population, I &&indicator) {
    auto &&[lb, ub] = m_objective_bounds(population);
    auto s = m_scale_objective_vectors(population, lb, ub);
    auto c = std::numeric_limits<ObjectiveVector::value_type>::min();
    for (std::size_t i = 0; i < s.size(); ++i) {
      for (std::size_t j = 0; j < s.size(); ++j) {
        if (i != j) {
          c = std::max(c, std::abs(indicator(s[i], s[j])));
        }
      }
    }
    return c;
  }

  /**
   * @brief Calculate the fitness values for the population individuals.
   *
   * @tparam I The type for the IBEA indicator.
   * @tparam S The type for the genetic algorithm's solution.
   * @param population The IBEA population to be scaled.
   * @param k IBEA scaling factor.
   * @param indicator The indicator used by IBEA.
   */
  template <typename I, typename S = GASolution>
  void m_fitness_assignment(std::vector<S> &population, double const k, I &&indicator) const {
    for (std::size_t i = 0; i < population.size(); ++i) {
      population[i].set_fitness(0);
      for (std::size_t j = 0; j < population.size(); ++j) {
        if (i != j) {
          population[i].set_fitness(population[i].fitness() -
                                    std::exp(-indicator(population[j], population[i]) / k));
        }
      }
    }
  }

  /**
   * @brief Do the environmental selection step on the IBEA population.
   *        Reduce the number of individuals (that was increased after the matting process)
   *        to the maximum size allowed.
   *
   * @tparam I The type for the IBEA indicator.
   * @tparam S The type for the genetic algorithm's solution.
   * @param population The IBEA population from which the individuals will be selected.
   * @param k IBEA scaling factor.
   * @param population_max_size The max population size.
   * @param indicator The indicator used by IBEA.
   */
  template <typename I, typename S = GASolution>
  void m_environmental_selection(std::vector<S> &population, double const k,
                                 std::size_t population_max_size, I &&indicator) {
    while (population.size() > population_max_size) {
      std::size_t worst = 0;
      for (std::size_t i = 0; i < population.size(); ++i) {
        worst = population[i].fitness() < population[worst].fitness() ? i : worst;
      }

      std::swap(population[worst], population.back());
      for (std::size_t i = 0; i < population.size() - 1; ++i) {
        population[i].set_fitness(population[i].fitness() +
                                  std::exp(-indicator(population.back(), population[i]) / k));
      }
      population.pop_back();
    }
  }
};
}  // namespace pmnk
#endif  // IBEA_HPP
