/**
 * @file utils.hpp
 * @author Pedro Rodrigues (pedror@student.dei.uc.pt)
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Project Utility functions.
 * @version 0.1.0
 * @date 13-09-2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef UTILS_HPP
#define UTILS_HPP

#include "solution.hpp"

namespace apmnkl {

namespace priv {

/**
 * @brief Utility function resposible for maintaining a
 *        container of non-dominated solutions
 *
 * @tparam Vec The type for the container holding the solutions.
 * @tparam S The type for a solution to be added to the solution container.
 * @param solutions The container for the non-dominated solutions.
 * @param solution The solution to be added to the container.
 * @return true If the solution was successfully added to the container
 * @return false If the solution to be added is dominated by other solution already
 *               present and failed to be inserted.
 */
template <typename Vec, typename S>
bool add_non_dominated(Vec &solutions, S &&solution) {
  for (std::size_t i = 0; i < solutions.size();) {
    auto d = solution.dominance(solutions[i]);
    if (d == dominance_type::equal) {
      if (solution.decision_vector() == solutions[i].decision_vector()) {
        return false;
      } else {
        for (++i; i < solutions.size(); ++i) {
          if (solution.decision_vector() == solutions[i].decision_vector()) {
            return false;
          }
        }
        break;
      }
    } else if (d == dominance_type::dominates) {
      solutions[i] = std::move(solutions.back());
      solutions.pop_back();
    } else if (d == dominance_type::dominated) {
      return false;
    } else {
      ++i;
    }
  }
  solutions.push_back(std::forward<S>(solution));
  return true;
}
}  // namespace priv
}  // namespace apmnkl
#endif  // UTILS_HPP
