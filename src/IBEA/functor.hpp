/**
 * @file functor.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Implementation of IBEA operators (using functors)
 * @version 0.1
 * @date 13-09-2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef IBEA_FUNCTOR_HPP
#define IBEA_FUNCTOR_HPP

#include <iostream>
#include <random>

#include "../Utils/solution.hpp"
#include "../Utils/wfg.hpp"

namespace pmnk {

/**
 * @brief Hypervolume Based IBEA indicator
 * 
 * @tparam R The IBEA reference point type.
 */
template <typename R = ObjectiveVector>
struct IHD {
  R m_ref;

  /**
   * @brief Construct a new IHD object
   * 
   * @param ref The reference point used for indicator calculation.
   */
  constexpr explicit IHD(R &&ref)
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
  template <typename S = GASolution>
  [[nodiscard]] double operator()(S const &s1, S const &s2) const {
    auto const &o1 = s1.objective_vector();
    auto const &o2 = s2.objective_vector();

    if (weakly_dominates(o1, o2)) {
      return point_hv(o2, m_ref) - point_hv(o1, m_ref);
    } else {
      hvobj<typename R::value_type> hvo(m_ref);
      hvo.insert(o1), hvo.insert(o2);
      return hvo.value() - point_hv(o1, m_ref);
    }
  }
};


/**
 * @brief Additive epsilon indicator.
 * 
 */
struct EPS {

  /**
   * @brief Construct a new EPS object.
   * 
   */
  constexpr explicit EPS() = default;

  /**
   * @brief Function call operator overload. Implements the 
   *        indicator functionality.
   * 
   * @tparam S The type used to store an genetic algorithm (IBEA) solution 
   * @param s1 A solution to be evaluated.
   * @param s2 A solution to be evaluated.
   * @return double The value for the indicator.
   */
  template <typename S = GASolution>
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


/**
 * @brief N-Point crossover operator.
 * 
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct NPointCrossover {
  std::size_t m_crossover_points;
  double m_crossover_probability;
  RNG m_rng;
  std::uniform_real_distribution<double> m_distrib;

  /**
   * @brief Construct a new NPointCrossover.
   * 
   * @param crossover_points Number of crossover points considered by this operator.
   * @param crossover_probability Crossover probability considered by this operator.
   * @param rng The random number generator object instance.
   */
  constexpr NPointCrossover(std::size_t const crossover_points, double const crossover_probability,
                            RNG &rng)
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
  template<typename S = GASolution>
  void operator()(S &s1, S &s2) noexcept {
    if (m_distrib(m_rng) < m_crossover_probability) {
      std::size_t p1 = 0, p2 = 0;
      for (std::size_t i = 0; i < m_crossover_points; ++i, p1 = p2) {
        std::uniform_int_distribution<std::size_t> randint(p1, s1.size() - 1);
        p2 = randint(m_rng);
        for (auto i = p1; i < p2; ++i) {
          std::swap(s1[i], s2[i]);
        };
      }
    }
  }
};

/**
 * @brief Uniform Crossover operator.
 * 
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct UniformCrossover {
  double m_crossover_probability;
  RNG m_rng;
  std::bernoulli_distribution m_distrib;

  /**
   * @brief Construct a new UniformCrossover object.
   * 
   * @param crossover_probability Crossover probability considered in this operator.
   * @param rng The random number generator object instance.
   */
  constexpr UniformCrossover(double crossover_probability, RNG &rng)
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
  template <typename S = GASolution>
  constexpr void operator()(S &s1, S &s2) noexcept {
    for (std::size_t i = 0; i < s1.size(); ++i) {
      if (m_distrib(m_rng)) {
        std::swap(s1[i], s2[i]);
      }
    }
  }
};

/**
 * @brief Uniform Mutation operator.
 * 
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct UniformMutation {
  double m_mutation_probability;
  RNG m_rng;
  std::uniform_real_distribution<double> m_distrib;

  /**
   * @brief Construct a new UniformMutation object.
   * 
   * @param mutation_probability Mutation probability considered in this operator.
   * @param rng The random number generator object instance.
   */
  constexpr UniformMutation(double const mutation_probability, RNG &rng)
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
  template <typename S = GASolution>
  constexpr void operator()(S &s) noexcept {
    for (std::size_t i = 0; i < s.size(); ++i) {
      if (m_distrib(m_rng) < m_mutation_probability) {
        s[i] = !s[i];
      }
    }
  }
};

/**
 * @brief KWayTournamentSelection operator. 
 * 
 * @tparam RNG The type for the random number generator object.
 */
template <typename RNG>
struct KWayTournamentSelection {
  std::size_t m_tournament_size;
  std::size_t m_matting_pool_size;
  RNG m_rng;

  /**
   * @brief Construct a new KWayTournamentSelection object
   *
   * @param tournament_size  The size of the tournament (K)
   * @param matting_pool_size The the maximum number of individuals allowed in the matting pool.
   * @param rng The random number generator objects.
   */
  constexpr KWayTournamentSelection(std::size_t const tournament_size,
                                    std::size_t const matting_pool_size, RNG &rng)
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
  template <typename S = GASolution>
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

}  // namespace pmnk
#endif  // IBEA_FUNCTOR_HPP
