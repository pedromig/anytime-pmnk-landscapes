/**
 * @file wfg.hpp
 * @author Alexandre Jesus (ajesus@dei.uc.pt)
 * @brief Implementation of the wfg algorithm to calculate hypervolumes
 * @version 0.2.0
 * @date 29-10-2021
 * @copyright Copyright (c) 2021
 */

#ifndef WFG_H
#define WFG_H

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>
#include <vector>

// This code assumes maximizing objective functions

namespace apmnkl {

namespace priv {

/**
 * @brief Check pareto weakly dominance between two points.
 *
 * @tparam T  The type for the point being supplied as a parameter.
 * @param lhs The first point (left hand side point)
 * @param lhs The second point (left hand side point)
 * @return auto A boolean with the value true if the left hand side
 *         point (lhs) weakly dominated the right hand side point (rhs)
 */
template <typename T>
auto weakly_dominates(T const& lhs, T const& rhs) {
  for (decltype(lhs.size()) i = 0; i < lhs.size(); ++i) {
    if (lhs[i] < rhs[i]) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Insert a point (solution) into a non dominated (solution)
 *        set of points
 *
 * @tparam T  The type for the point being supplied as a parameter.
 * @param sol The point (solution) to be added to the set
 * @param set The (solution) set of non dominated points (solutions)
 */
template <typename T>
void insert_non_dominated(T&& sol, std::vector<T>& set) {
  for (auto it = set.begin(); it != set.end(); ++it) {
    if (weakly_dominates(*it, sol)) {
      return;
    } else if (weakly_dominates(sol, *it)) {
      *it = std::move(set.back());
      set.pop_back();
      set.erase(
          std::remove_if(it, set.end(), [&sol](auto const& s) { return weakly_dominates(sol, s); }),
          set.end());
      break;
    }
  }
  set.push_back(std::move(sol));
}

/**
 * @brief Limit the set of points by replacing them with the points
 * whose value in each objective is limited to be no better than the
 * contributing point
 *
 * @tparam Iter The type for and iterator of a set of points
 * @tparam T The type for the contributing point being supplied as a parameter
 * @param begin An iterator for the begining of the set
 * @param end An iterator for the end of the set
 * @param sol The contributing solution
 * @return auto The limited set
 */
template <typename Iter, typename T>
auto limit_set(Iter begin, Iter end, T const& sol) {
  std::vector<T> res;
  res.reserve(static_cast<typename std::vector<T>::size_type>(std::distance(begin, end)));
  for (; begin != end; ++begin) {
    auto aux = *begin;
    for (size_t i = 0; i < sol.size(); ++i) {
      aux[i] = std::min(aux[i], sol[i]);
    }
    insert_non_dominated(std::move(aux), res);
  }
  return res;
}

/**
 * @brief Compute hypervolume of a point in relation to a reference.
 *
 * @tparam T The type for the point being supplied as a parameter
 * @tparam R The type for the reference being supplied as a parameter
 * @param p The point whose hypervolume value is to be computed
 * @param r The reference point.
 * @return auto  The hypervolume value.
 */
template <typename T, typename R>
auto point_hv(T const& p, R const& r) {
  auto res = p[0] - r[0];
  for (size_t i = 1; i < p.size(); ++i) {
    res *= p[i] - r[i];
  }
  return res;
}

/**
 * @brief Compute a set hypervolume value given a reference point
 *        using the wfg algorithm. (Worker function)
 *
 * @tparam S The type for the set supplied as a parameter
 * @tparam T The type for the point being supplied as a parameter
 * @param s The solution set
 * @param ref  The reference point.
 * @return std::common_type_t<typename S::value_type::value_type, typename T::value_type>
 * The resulting set hypervolume value.
 */
template <typename S, typename T>
std::common_type_t<typename S::value_type::value_type, typename T::value_type> set_hv_wfg(
    S const& s, T const& ref) {
  using result_t = std::common_type_t<typename S::value_type::value_type, typename T::value_type>;
  auto res = result_t(0);
  for (auto it = s.begin(); it != s.end(); ++it) {
    res += point_hv(*it, ref) - set_hv_wfg(limit_set(std::next(it), s.end(), *it), ref);
  }
  return res;
}

/**
 * @brief Compute a set hypervolume value given a reference point
 *        using the wfg algorithm. (Wrapper function)
 *
 * @tparam S The type for the set supplied as a parameter
 * @tparam T The type for the point being supplied as a parameter
 * @param s The solution set
 * @param ref  The reference point.
 * @return std::common_type_t<typename S::value_type::value_type, typename T::value_type>
 * The resulting set hypervolume value.
 */
template <typename S, typename T>
auto set_hv(S const& s, T const& ref) {
  std::vector<T> v;
  v.reserve(s.size());
  for (auto const& sol : s) {
    v.push_back(sol.objective_vector());
  }
  sort(v.begin(), v.end(), [](auto const& a, auto const& b) { return a[0] < b[0]; });
  return set_hv_wfg(v, ref);
}

/**
 * @brief Calculate a point hypervolume contribution to a set.
 *
 * @tparam S The type for the set supplied as a parameter.
 * @tparam T The type for the point being supplied as a parameter.
 * @tparam R The type for the reference point supplied as a parameter.
 * @param p The point whose contribution is going to calculated.
 * @param s The solution set.
 * @param ref  The reference point.
 * @return auto The point hypervolume contribution.
 */
template <typename T, typename S, typename R>
auto point_hvc(T const& p, S const& s, R const& ref) {
  std::vector<T> v;
  v.reserve(s.size());
  for (auto const& sol : s) {
    v.push_back(sol.objective_vector());
  }
  sort(v.begin(), v.end(), [](auto const& a, auto const& b) { return a[0] < b[0]; });
  return point_hv(p, ref) - set_hv_wfg(limit_set(v.begin(), v.end(), p), ref);
}

/// Implementation of an API that supports among others, set/point hypervolume calculations (using
/// the WFG algorithm)
template <typename T>
class [[nodiscard]] hvobj {
 public:
  using hv_type = T;
  using ovec_type = std::vector<hv_type>;
  using set_type = std::vector<ovec_type>;

  constexpr explicit hvobj(ovec_type const& r)
      : m_hv(0)
      , m_set()
      , m_ref(r) {}

  constexpr hvobj(hvobj const& other) = default;
  constexpr hvobj(hvobj&& other) noexcept = default;

  /// Get the current hypervolume value
  [[nodiscard]] constexpr auto value() const {
    return m_hv;
  }

  /// Get the contribution of a new vector w.r.t. to the current set
  template <typename V>
  [[nodiscard]] constexpr auto contribution(V const& v) const {
    return m_point_hv(v, m_ref) - m_set_hv(m_limit_set(m_set, v), m_ref);
  }

  /// Inserts a new objective vector and returns its contribution
  template <typename V>
  constexpr auto insert(V&& v) {
    auto hvc = contribution(v);
    if (hvc != 0) {
      m_insert_non_dominated(std::forward<V>(v), m_set);
      m_hv += hvc;
    }
    return hvc;
  }

  /// Removes a objective vector and returns its contribution (i.e. the lost hv) or -1.0 if no
  /// objective vector was found.
  template <typename V>
  constexpr auto remove(V const& v) {
    auto it = std::find(m_set.begin(), m_set.end(), v);
    if (it == m_set.end())
      return -1.0;
    m_set.erase(it);
    auto hvc = contribution(v);
    m_hv -= hvc;
    return hvc;
  }

 private:
  /// Check if a[1..] >= b[1..]
  template <typename V>
  [[nodiscard]] constexpr auto m_weakly_dominates(V const& a, V const& b) const {
    for (size_t i = 1; i < a.size(); ++i) {
      if (a[i] < b[i]) {
        return false;
      }
    }
    return true;
  }

  /// insert non dominated point into a set implementation
  template <typename V, typename C>
  void m_insert_non_dominated(V&& v, C& set) const {
    auto it = set.begin();

    for (; it != set.end() && (*it)[0] > v[0]; ++it) {
      if (m_weakly_dominates(*it, v)) {
        return;
      }
    }

    for (; it != set.end() && (*it)[0] == v[0]; ++it) {
      if (m_weakly_dominates(*it, v)) {
        return;
      } else if (m_weakly_dominates(v, *it)) {
        *it = std::forward<V>(v);
        set.erase(std::remove_if(std::next(it), set.end(),
                                 [this, it](auto const& a) { return m_weakly_dominates(*it, a); }),
                  set.end());
        return;
      }
    }

    if (it == set.end()) {
      set.push_back(std::forward<V>(v));
    } else {
      auto aux = std::forward<V>(v);
      std::swap(aux, *it);
      for (auto jt = std::next(it); jt != set.end(); ++jt) {
        if (m_weakly_dominates(*it, aux)) {
          set.erase(
              std::remove_if(jt, set.end(),
                             [this, it](auto const& a) { return m_weakly_dominates(*it, a); }),
              set.end());
          return;
        } else {
          std::swap(aux, *jt);
        }
      }
      if (!m_weakly_dominates(*it, aux)) {
        set.push_back(std::move(aux));
      }
    }
  }

  // limit_set (implementation)
  template <typename S, typename V>
  auto m_limit_set(S const& s, V const& v) const {
    S res;
    res.reserve(s.size());
    for (auto const& p : s) {
      auto aux = p;
      for (size_t i = 0; i < aux.size(); ++i) {
        aux[i] = std::min(aux[i], v[i]);
      }
      m_insert_non_dominated(std::move(aux), res);
    }
    return res;
  }

  /// hypervolume of a point (implementation)
  template <typename V, typename R>
  auto m_point_hv(V const& v, R const& r) const {
    auto res = v[0] - r[0];
    for (size_t i = 1; i < v.size(); ++i) {
      res *= v[i] - r[i];
    }
    return res;
  }

  /// 3d hypervolume of a set (implementation)
  template <typename S, typename R>
  auto m_set_hv3d(S const& s, R const& r) const {
    using array2_t = std::array<hv_type, 2>;

    auto aux = std::vector<array2_t>{{r[1], std::numeric_limits<hv_type>::max()},
                                     {std::numeric_limits<hv_type>::max(), r[2]}};

    hv_type v = 0;
    hv_type a = 0;
    hv_type z = 0;

    for (auto const& p : s) {
      v += a * (z - p[0]);
      z = p[0];

      auto tmp = array2_t{p[1], p[2]};
      auto it = std::lower_bound(aux.begin(), aux.end(), tmp,
                                 [](auto const& x, auto const& y) { return x[1] > y[1]; });
      auto jt = it;

      auto r0 = (*std::prev(it))[0];
      auto r1 = tmp[1];
      for (; (*it)[0] <= tmp[0]; ++it) {
        a += (tmp[0] - r0) * (r1 - (*it)[1]);
        r0 = (*it)[0];
        r1 = (*it)[1];
      }
      a += (tmp[0] - r0) * (r1 - (*it)[1]);
      if (jt != it) {
        *jt = tmp;
        aux.erase(++jt, it);
      } else {
        aux.insert(it, tmp);
      }
    }
    v += a * (z - r[0]);
    return v;
  }

  /// Assumes set is sorted by increasing i
  template <typename S, typename R>
  auto m_set_hv(S const& s, R const& r, hv_type c = 1) const -> hv_type {
    hv_type v = 0;

    if (s.size() == 0) {
      return 0;
    }

    if (s.begin()->size() == 2) {
      hv_type r1 = r[1];
      for (auto const& p : s) {
        v += (p[1] - r1) * (p[0] - r[0]);
        r1 = p[1];
      }
      v *= c;
    } else if (s.begin()->size() == 3) {
      v = c * m_set_hv3d(s, r);
    } else {
      auto newr = std::vector<hv_type>();
      for (size_t i = 1; i < r.size(); ++i) {
        newr.push_back(r[i]);
      }
      auto newl = std::vector<std::vector<hv_type>>();
      newl.reserve(s.size());
      for (auto const& p : s) {
        auto newc = c * (p[0] - r[0]);
        auto newp = std::vector<hv_type>();
        for (size_t i = 1; i < p.size(); ++i) {
          newp.push_back(p[i]);
        }

        v += newc * m_point_hv(newp, newr) - m_set_hv(m_limit_set(newl, newp), newr, newc);

        m_insert_non_dominated(std::move(newp), newl);
      }
    }
    return v;
  }

  hv_type m_hv;
  set_type m_set;
  ovec_type m_ref;
};

}  // namespace priv
}  // namespace apmnkl
#endif  // WFG_H
