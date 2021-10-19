/**
 * @file solution.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Implementation of a solution class.
 * @version 0.1.0
 * @date 13-09-2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <cassert>
#include <random>

#include "rMNKEval.hpp"

namespace pmnk {

using DecisionVector = std::vector<bool>;
using ObjectiveVector = std::vector<double>;

enum class DominanceType { DOMINATES, EQUAL, DOMINATED, INCOMPARABLE };

/// Standard solution class
class Solution {
 protected:
  DecisionVector m_decision;
  ObjectiveVector m_objective;

 public:
  explicit Solution() = default;

  /**
   * @brief Construct a new Solution object (copy constructor)
   * 
   * @param other A const lvalue reference to a solution to be coppied
   */
  Solution(Solution const &other) = default;

  /**
   * @brief Construct a new Solution object (move constructor)
   * 
   * @param other A rvalue reference to a solution to be moved
   */
  Solution(Solution &&other) = default;

  /**
   * @brief Defaulted copy assignment operator
   * 
   * @param other A const lvalue reference to a solution to be coppied.
   * @return Solution& A lvalue reference to the solution to be assigned.
   */
  Solution &operator=(Solution const &other) = default;

  /**
   * @brief Defaulted move assignment operator 
   * 
   * @param other A rvalue reference to solution to be coppied 
   * @return Solution& A refernce to the solution to be assigned. 
   */
  Solution &operator=(Solution &&other) = default;

  /**
   * @brief Construct a new Solution object
   * 
   * @param rmnk A lvalue reference to the RMNK instance evaluator.
   * @param decision The solution's decision vector.
   */
  Solution(RMNKEval &rmnk, DecisionVector const &decision)
      : m_decision(decision) {
    eval(rmnk);
  }

  /**
   * @brief Construct a new Solution object
   * 
   * @param rmnk A lvalue reference to the RMNK instance evaluator.
   * @param decision The solution's decision vector.
   */
  Solution(RMNKEval &rmnk, DecisionVector &&decision)
      : m_decision(std::move(decision)) {
    eval(rmnk);
  }

  /**
   * @brief Getter for the solution's decision vector.
   * 
   * @return DecisionVector const& A Read-Only reference 
   *         to the solution's decision vector.
   */
  DecisionVector const &decision_vector() const {
    return m_decision;
  }

  /**
   * @brief Getter for the solution's objective vector.
   * 
   * @return ObjectiveVector const& A Read-Only reference to 
   *         the solution's objective vector.
   */
  ObjectiveVector const &objective_vector() const {
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
   * @return DominanceType The solution's dominance type 
   */
  DominanceType dominance(Solution const &solution) const {
    assert(m_decision.size() == solution.m_decision.size());

    DominanceType res = DominanceType::EQUAL;
    for (decltype(m_objective.size()) i = 0; i < m_objective.size(); ++i) {
      if (m_objective[i] < solution.m_objective[i]) {
        if (res == DominanceType::DOMINATES) {
          return DominanceType::INCOMPARABLE;
        }
        res = DominanceType::DOMINATED;
      } else if (m_objective[i] > solution.m_objective[i]) {
        if (res == DominanceType::DOMINATED) {
          return DominanceType::INCOMPARABLE;
        }
        res = DominanceType::DOMINATES;
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
  friend std::ostream &operator<<(std::ostream &os, Solution const &solution) {
    for (auto const &i : solution.objective_vector()) {
      os << i << " ";
    }
    os << "\n";
    return os;
  }

  /**
   * @brief Array Indexing operator overload. This operator 
   *        provides a way to access the i-th bit in the solution's 
   *        decision vector implementation. 
   * 
   * @param i The index of the element to be accessed.
   * @return DecisionVector::reference  A lvalue reference to value contained in the 
   *                                    accessed index.
   */
  DecisionVector::reference operator[](std::size_t const i) {
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
   * @return Solution A new Solution object containing a random solution.
   */
  template <typename RNG>
  static Solution random_solution(RMNKEval &eval, RNG &generator) {
    std::uniform_int_distribution<int> distrib(0, 1);

    DecisionVector decision_vector(eval.getN(), 0);
    for (std::size_t i = 0; i < decision_vector.size(); ++i)
      decision_vector[i] = distrib(generator);

    return Solution(eval, std::move(decision_vector));
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
   * @return Solution A new Solution object containing a random solution.
   */
  template <typename RNG>
  static Solution uniform_bit_flip_solution(RMNKEval &eval, RNG &generator,
                                            Solution const &original) {
    std::bernoulli_distribution distrib(1 / (double)original.decision_vector().size());

    DecisionVector flipped = original.decision_vector();
    for (decltype(flipped.size()) i = 0; i < flipped.size(); ++i) {
      if (distrib(generator)) {
        flipped[i] = !flipped[i];
      }
    }
    return Solution(eval, std::move(flipped));
  }

  /**
   * @brief Calculate all the neighboor solutions of the current one.
   * 
   * @param eval  The instance evaluator object.
   * @param generator The random number generator object 
   * @return std::vector<Solution> A vector of solutions containing the current
   *                               solution neighboor solutions.
   *                            
   */
  static std::vector<Solution> neighborhood_solutions(RMNKEval &eval, Solution const &original) {
    std::vector<Solution> neighborhood;
    neighborhood.reserve(original.decision_vector().size());

    for (size_t i = 0; i < original.decision_vector().size(); ++i) {
      auto decision_vector = original.decision_vector();
      decision_vector[i] = !decision_vector[i];

      neighborhood.emplace_back(eval, std::move(decision_vector));
    }

    for (size_t i = 0; i < original.decision_vector().size(); ++i) {
      for (size_t j = i + 1; j < original.decision_vector().size(); ++j) {
        if (original.decision_vector()[i] == original.decision_vector()[j]) {
          continue;
        }

        auto decision_vector = original.decision_vector();
        swap(decision_vector[i], decision_vector[j]);
        neighborhood.emplace_back(eval, std::move(decision_vector));
      }
    }
    return neighborhood;
  }
};

/// Genetic algorithm solution wrapper (adds a fitness attribute to the Solution class)
class GASolution : public Solution {
  double m_fitness;

 public:

  /**
  * @brief Construct a new GASolution object
  * 
  */
  explicit GASolution() = default;

  /**
   * @brief Construct a new GASolution object (move)
   * 
   * @param sol A rvalue reference to a solution.
   * @param fitness The fitness value of the solution.
   */
  GASolution(Solution &&sol, double fitness)
      : Solution(std::forward<Solution>(sol))
      , m_fitness(fitness) {}

  /**
   * @brief Construct a new GASolution object
   * 
   * @param sol A rvalue reference to a solution. 
   *            (fitness value defaults to 0)
   */
  GASolution(Solution &&sol)
      : GASolution(std::forward<Solution>(sol), 0) {}

  /**
   * @brief Getter method for the GASolution's fitness value.
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
  void set_objv(ObjectiveVector &&objv) {
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
}  // namespace pmnk
#endif  // SOLUTION_HPP
