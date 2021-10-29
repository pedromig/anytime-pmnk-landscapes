/**
 * @file operators.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Implementation of IBEA operators (using functors)
 * @version 0.2.0
 * @date 29-10-2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IBEA_OPERATORS_HPP
#define IBEA_OPERATORS_HPP

#include <iostream>
#include <random>

#include "utils/solution.hpp"
#include "utils/wfg.hpp"

namespace apmnkl {

namespace indicator {
/**
 * @brief Hypervolume Based IBEA indicator
 *
 * @tparam R The IBEA reference point type.
 */
template <typename R = objective_vector>
struct ihd {
  R m_ref;

  /**
   * @brief Construct a new ihd object
   *
   * @param ref The reference point used for indicator calculation.
   */
  constexpr explicit ihd(R &&ref)
      : m_ref(std::forward<R>(ref)) {}

  /**
   * @brief Function call operator overload. Implements the
   *        indicator functionality.
   *
   * @tparam S The type used to store an genetic algorithm (IBEA) solution
   * @param s1 A solution to be evaluated.
   * @param s2 A solution to be evaluated.
   * @return double The value for the indicator.
   */
  template <typename S = priv::gasolution>
  [[nodiscard]] double operator()(S const &s1, S const &s2) const {
    auto const &o1 = s1.objective_vector();
    auto const &o2 = s2.objective_vector();

    if (priv::weakly_dominates(o1, o2)) {
      return priv::point_hv(o2, m_ref) - priv::point_hv(o1, m_ref);
    } else {
      priv::hvobj<typename R::value_type> hvo(m_ref);
      hvo.insert(o1), hvo.insert(o2);
      return hvo.value() - priv::point_hv(o1, m_ref);
    }
  }
};

/**
 * @brief Additive epsilon indicator.
 *
 */
struct eps {
  /**
   * @brief Construct a new eps object.
   *
   */
  constexpr explicit eps() = default;

  /**
   * @brief Function call operator overload. Implements the
   *        indicator functionality.
   *
   * @tparam S The type used to store an genetic algorithm (IBEA) solution
   * @param s1 A solution to be evaluated.
   * @param s2 A solution to be evaluated.
   * @return double The value for the indicator.
   */
  template <typename S = priv::gasolution>
  [[nodiscard]] constexpr double operator()(S const &s1, S const &s2) const noexcept {
    double indicator = std::numeric_limits<double>::min();
    for (std::size_t i = 0; i < s1.objective_vector().size(); ++i) {
      double const o1 = s1.objective_vector()[i];
      double const o2 = s2.objective_vector()[i];
      indicator = std::max(indicator, o2 - o1);
    }
    return indicator;
  }
};
}  // namespace indicator

namespace crossover {

/**
 * @brief N-Point crossover operator.
 *
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct n_point_crossover {
  std::size_t m_crossover_points;
  double m_crossover_probability;
  RNG m_rng;
  std::uniform_real_distribution<double> m_distrib;

  /**
   * @brief Construct a new n_point_crosover object.
   *
   * @param crossover_points Number of crossover points considered by this operator.
   * @param crossover_probability Crossover probability considered by this operator.
   * @param rng The random number generator object instance.
   */
  constexpr n_point_crossover(std::size_t const crossover_points,
                              double const crossover_probability, RNG &rng)
      : m_crossover_points(crossover_points)
      , m_crossover_probability(crossover_probability)
      , m_rng(rng)
      , m_distrib(0.0, 1.0) {}

  /**
   * @brief Function call operator overload. Implements the
   *        crossover operator functionality.
   *
   * @param s1 A solution to be recombinated.
   * @param s2 A solution to be recombinated.
   */
  template <typename S = priv::gasolution>
  void operator()(S &s1, S &s2) noexcept {
    if (m_distrib(m_rng) < m_crossover_probability) {
      std::size_t p1 = 0, p2 = 0;
      for (std::size_t i = 0; i < m_crossover_points; ++i, p1 = p2) {
        std::uniform_int_distribution<std::size_t> randint(p1, s1.size() - 1);
        p2 = randint(m_rng);
        for (std::size_t j = p1; j < p2; ++j) {
          std::swap(s1[j], s2[j]);
        };
      }
    }
  }
};

/// Helper `using` to avoid typing full class name.
template <typename RNG>
using npc = apmnkl::crossover::n_point_crossover<RNG>;

/**
 * @brief Uniform Crossover operator.
 *
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct uniform_crossover {
  double m_crossover_probability;
  RNG m_rng;
  std::bernoulli_distribution m_distrib;

  /**
   * @brief Construct a new uniform_crossover object.
   *
   * @param crossover_probability Crossover probability considered in this operator.
   * @param rng The random number generator object instance.
   */
  constexpr uniform_crossover(double crossover_probability, RNG &rng)
      : m_crossover_probability(crossover_probability)
      , m_rng(rng)
      , m_distrib() {}

  /**
   * @brief Function call operator overload. Implements the
   *        crossover operator functionality.
   *
   * @tparam S The type used to store an genetic algorithm (IBEA) solution
   * @param s1 A solution to be recombinated.
   * @param s2 A solution to be recombinated.
   */
  template <typename S = priv::gasolution>
  constexpr void operator()(S &s1, S &s2) noexcept {
    for (std::size_t i = 0; i < s1.size(); ++i) {
      if (m_distrib(m_rng)) {
        std::swap(s1[i], s2[i]);
      }
    }
  }
};

/// Helper `using` to avoid typing full class name.
template <typename RNG>
using uc = apmnkl::crossover::uniform_crossover<RNG>;

}  // namespace crossover

namespace mutation {

/**
 * @brief Uniform Mutation operator.
 *
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct uniform_mutation {
  double m_mutation_probability;
  RNG m_rng;
  std::uniform_real_distribution<double> m_distrib;

  /**
   * @brief Construct a new uniform_mutation object.
   *
   * @param mutation_probability Mutation probability considered in this operator.
   * @param rng The random number generator object instance.
   */
  constexpr uniform_mutation(double const mutation_probability, RNG &rng)
      : m_mutation_probability(mutation_probability)
      , m_rng(rng)
      , m_distrib(0.0, 1.0) {}

  /**
   * @brief Function call operator overload. Implements the
   *        crossover operator functionality.
   *
   * @tparam S The type used to store an genetic algorithm (IBEA) solution
   * @param s A solution to be mutated.
   */
  template <typename S = priv::gasolution>
  constexpr void operator()(S &s) noexcept {
    for (std::size_t i = 0; i < s.size(); ++i) {
      if (m_distrib(m_rng) < m_mutation_probability) {
        s[i] = !s[i];
      }
    }
  }
};

/// Helper `using` to avoid typing full class name.
template <typename RNG>
using um = apmnkl::mutation::uniform_mutation<RNG>;

}  // namespace mutation

namespace selection {
/**
 * @brief K-Way Tournament Selection operator.
 *
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct k_way_tournament {
  std::size_t m_tournament_size;
  std::size_t m_matting_pool_size;
  RNG m_rng;

  /**
   * @brief Construct a new k_way_tournament object.
   *
   * @param tournament_size The size of the tournament (K)
   * @param matting_pool_size The the maximum number of individuals allowed in the matting pool.
   * @param rng The random number generator objects.
   */
  constexpr k_way_tournament(std::size_t const tournament_size, std::size_t const matting_pool_size,
                             RNG &rng)
      : m_tournament_size(tournament_size)
      , m_matting_pool_size(matting_pool_size)
      , m_rng(rng){};

  /**
   * @brief Function call operator overload. Implements the
   *        selection operator functionality.
   *
   * @tparam S The type used to store an genetic algorithm (IBEA) solution
   * @param population The population from which the individuals will be selected.
   * @return std::vector<S> The matting pool obtain as a result from the selection
   *                        of the individuals of the population.
   */
  template <typename S = priv::gasolution>
  [[nodiscard]] std::vector<S> operator()(std::vector<S> const &population) noexcept {
    std::vector<S> matting_pool;
    matting_pool.reserve(m_matting_pool_size);
    std::uniform_int_distribution<std::size_t> distrib(0, population.size() - 1);

    for (std::size_t i = 0; i < m_matting_pool_size; ++i) {
      auto best = distrib(m_rng);
      for (std::size_t j = 0; j < m_tournament_size - 1; ++j) {
        auto other = distrib(m_rng);
        best = population[other].fitness() > population[best].fitness() ? other : best;
      }
      matting_pool.push_back(population[best]);
    }
    return matting_pool;
  }
};

/// Helper `using` to avoid typing full class name.
template <typename RNG>
using kwt = apmnkl::selection::k_way_tournament<RNG>;

}  // namespace selection
}  // namespace apmnkl
#endif  // IBEA_OPERATORS_HPP
