/**
 * @file solution.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Implementation of a solution class.
 * @version 0.2.0
 * @date 29-10-2021
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <cassert>
#include <random>

#include "rMNKEval.hpp"

namespace apmnkl {

using decision_vector = std::vector<bool>;
using objective_vector = std::vector<double>;

namespace priv {

enum class dominance_type { dominates, equal, dominated, incomparable };

/// Standard solution class
class solution {
  using objv_type = apmnkl::objective_vector;
  using decv_type = apmnkl::decision_vector;

 protected:
  decv_type m_decision;
  objv_type m_objective;

 public:
  explicit solution() = default;

  /**
   * @brief Construct a new solution object (copy constructor)
   *
   * @param other A const lvalue reference to a solution to be coppied
   */
  solution(solution const &other) = default;

  /**
   * @brief Construct a new solution object (move constructor)
   *
   * @param other A rvalue reference to a solution to be moved
   */
  solution(solution &&other) = default;

  /**
   * @brief Defaulted copy assignment operator
   *
   * @param other A const lvalue reference to a solution to be coppied.
   * @return solution& A lvalue reference to the solution to be assigned.
   */
  solution &operator=(solution const &other) = default;

  /**
   * @brief Defaulted move assignment operator
   *
   * @param other A rvalue reference to solution to be coppied
   * @return solution& A reference to the solution to be assigned.
   */
  solution &operator=(solution &&other) = default;

  /**
   * @brief Construct a new solution object
   *
   * @param rmnk A lvalue reference to the RMNK instance evaluator.
   * @param decision The solution's decision vector.
   */
  solution(RMNKEval &rmnk, decision_vector const &decision)
      : m_decision(decision) {
    eval(rmnk);
  }

  /**
   * @brief Construct a new solution object
   *
   * @param rmnk A lvalue reference to the RMNK instance evaluator.
   * @param decision The solution's decision vector.
   */
  solution(RMNKEval &rmnk, decision_vector &&decision)
      : m_decision(std::move(decision)) {
    eval(rmnk);
  }

  /**
   * @brief Getter for the solution's decision vector.
   *
   * @return decv_type const& A Read-Only reference
   *         to the solution's decision vector.
   */
  decv_type const &decision_vector() const {
    return m_decision;
  }

  /**
   * @brief Getter for the solution's objective vector.
   *
   * @return objv_type const& A Read-Only reference to
   *         the solution's objective vector.
   */
  objv_type const &objective_vector() const {
    return m_objective;
  }

  /**
   * @brief Getter for the solution's decision/objective vector size.
   *
   * @return std::size_t The size of this solution's decision/objective vector.
   */
  std::size_t size() const noexcept {
    return m_decision.size();
  };

  /**
   * @brief Calculate the objective dominance type of this solution
   *        with respect to another.
   *
   * @param solution A solution whose dominance type of this will be tested against.
   * @return dominance_type The solution's dominance type
   */
  dominance_type dominance(solution const &s) const {
    assert(m_decision.size() == s.m_decision.size());

    auto res = dominance_type::equal;
    for (decltype(m_objective.size()) i = 0; i < m_objective.size(); ++i) {
      if (m_objective[i] < s.m_objective[i]) {
        if (res == dominance_type::dominates) {
          return dominance_type::incomparable;
        }
        res = dominance_type::dominated;
      } else if (m_objective[i] > s.m_objective[i]) {
        if (res == dominance_type::dominated) {
          return dominance_type::incomparable;
        }
        res = dominance_type::dominates;
      }
    }
    return res;
  }

  /**
   * @brief Extraction operator overload for this object.
   *
   * @param os The output stream to where the solution string representation
   *           will be redirected.
   * @param solution The solution whose representation will be inserted in
   *                 the output stream
   * @return std::ostream& A output stream object with the solution's data appended to it.
   */
  friend std::ostream &operator<<(std::ostream &os, solution const &solution) {
    for (auto const &i : solution.objective_vector()) {
      os << i << " ";
    }
    return os;
  }

  /**
   * @brief Array Indexing operator overload. This operator
   *        provides a way to access the i-th bit in the solution's
   *        decision vector implementation.
   *
   * @param i The index of the element to be accessed.
   * @return decv_type::reference A lvalue reference to value contained in the
   *                              accessed index.
   */
  decv_type::reference operator[](std::size_t const i) {
    return m_decision[i];
  }

  /**
   * @brief Wrapper method that calls the RMNK solution evaluator suplied
   *        on the solution's decision vector, calculating the respective
   *        objective vector.
   *
   * @param rmnk The RMNK instance evaluator instance that provides the method used
   *             for solution evaluation.
   */
  void eval(RMNKEval &rmnk) {
    rmnk.eval(m_decision, m_objective);
  }

  /**
   * @brief Build and evaluate a new random solution object.
   *
   * @tparam RNG The type for the random number generator object.
   * @param eval  The instance evaluator object.
   * @param generator The random number generator object
   * @return solution A new Solution object containing a random solution.
   */
  template <typename RNG>
  static solution random_solution(RMNKEval &eval, RNG &generator) {
    std::uniform_int_distribution<int> distrib(0, 1);

    decv_type decision_vector(eval.getN(), 0);
    for (std::size_t i = 0; i < decision_vector.size(); ++i)
      decision_vector[i] = distrib(generator);

    return solution(eval, std::move(decision_vector));
  }

  /**
   * @brief Build and evaluate a new Random solution object
   *        that results from a mutation in another solution's
   *        decision vector bit representation.
   *
   *
   * @tparam RNG The type for the random number generator object.
   * @param eval  The instance evaluator object.
   * @param generator The random number generator object
   * @param original The parent solution to be mutated.
   * @return solution A new olution object containing a random solution.
   */
  template <typename RNG>
  static solution uniform_bit_flip_solution(RMNKEval &eval, RNG &generator,
                                            solution const &original) {
    std::bernoulli_distribution distrib(1 / static_cast<double>(original.decision_vector().size()));

    decv_type flipped = original.decision_vector();
    for (decltype(flipped.size()) i = 0; i < flipped.size(); ++i) {
      if (distrib(generator)) {
        flipped[i] = !flipped[i];
      }
    }
    return solution(eval, std::move(flipped));
  }

  /**
   * @brief Calculate all the neighboor solutions of the current one.
   *
   * @param eval  The instance evaluator object.
   * @param generator The random number generator object
   * @return std::vector<solution> A vector of solutions containing the current
   *                               solution neighboor solutions.
   *
   */
  static std::vector<solution> neighborhood_solutions(RMNKEval &eval, solution const &original) {
    std::vector<solution> neighborhood;
    neighborhood.reserve(original.decision_vector().size());

    for (size_t i = 0; i < original.decision_vector().size(); ++i) {
      auto decv = original.decision_vector();
      decv[i] = !decv[i];

      neighborhood.emplace_back(eval, std::move(decv));
    }

    for (size_t i = 0; i < original.decision_vector().size(); ++i) {
      for (size_t j = i + 1; j < original.decision_vector().size(); ++j) {
        if (original.decision_vector()[i] == original.decision_vector()[j]) {
          continue;
        }

        auto decv = original.decision_vector();
        swap(decv[i], decv[j]);
        neighborhood.emplace_back(eval, std::move(decv));
      }
    }
    return neighborhood;
  }
};

/// Genetic algorithm solution wrapper (adds a fitness attribute to the Solution class)
class gasolution : public solution {
  using objv_type = apmnkl::objective_vector;
  using decv_type = apmnkl::decision_vector;

  double m_fitness;

 public:
  /**
   * @brief Construct a new gasolution object
   *
   */
  explicit gasolution() = default;

  /**
   * @brief Construct a new gasolution object (move)
   *
   * @param sol A rvalue reference to a solution.
   * @param fitness The fitness value of the solution.
   */
  gasolution(solution &&sol, double fitness)
      : solution(std::forward<solution>(sol))
      , m_fitness(fitness) {}

  /**
   * @brief Construct a new gasolution object
   *
   * @param sol A rvalue reference to a solution.
   *            (fitness value defaults to 0)
   */
  gasolution(solution &&sol)
      : gasolution(std::forward<solution>(sol), 0) {}

  /**
   * @brief Getter method for the gasolution's fitness value.
   *
   * @return constexpr double const& The fitness value.
   */
  [[nodiscard]] constexpr double const &fitness() const {
    return m_fitness;
  }

  /**
   * @brief Set the objv object of the current solution.
   *
   * @param objv A rvalue reference to the objective vector.
   */
  void set_objv(objv_type &&objv) {
    m_objective = std::move(objv);
  }

  /**
   * @brief Setter method for the fitness value
   *
   * @param fitness The fitness value.
   */
  constexpr void set_fitness(double const fitness) {
    m_fitness = fitness;
  }
};
}  // namespace priv
}  // namespace apmnkl
#endif  // SOLUTION_HPP
