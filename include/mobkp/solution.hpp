#ifndef MOBKP_SOLUTION_HPP_
#define MOBKP_SOLUTION_HPP_

#include <algorithm>
#include <functional>
#include <utility>

namespace mobkp {

template <typename Problem, typename DVec, typename CVec, typename OVec>
class solution {
 public:
  using problem_type = Problem;
  using decision_vector_type = DVec;
  using decision_vector_size_type = typename decision_vector_type::size_type;
  using objective_vector_type = CVec;
  using objective_vector_size_type = typename objective_vector_type::size_type;
  using constraint_vector_type = OVec;
  using constraint_vector_size_type = typename constraint_vector_type::size_type;

  solution(problem_type const &problem, decision_vector_type &&dvec, objective_vector_type &&ovec,
           constraint_vector_type &&cvec)
      : m_problem(problem)
      , m_dvec(std::move(dvec))
      , m_ovec(std::move(ovec))
      , m_cvec(std::move(cvec)) {}

  /**
   * Generates an empty solution, i.e., a solution with no items
   * selected for the knapsack.
   */
  static auto empty(problem_type const &problem) {
    auto dvec = decision_vector_type(problem.num_items(), false);
    auto ovec = objective_vector_type(problem.num_objectives(), 0);
    auto cvec = constraint_vector_type(problem.num_constraints(), 0);
    return solution(problem, std::move(dvec), std::move(ovec), std::move(cvec));
  }

  /**
   * Returns a reference to the decision vector of this solution, such
   * that each decision variable is given by a boolean that denotes
   * whether or not an item is in the knapsack.
   */
  [[nodiscard]] constexpr auto decision_vector() const -> decision_vector_type const & {
    return m_dvec;
  }

  /**
   * Returns a reference to the objective vector of this solution, that
   * denotes the sums of the values of the selected items for each
   * objective.
   */
  [[nodiscard]] constexpr auto objective_vector() const -> objective_vector_type const & {
    return m_ovec;
  }

  /**
   * Returns a reference to the objective vector of this solution, that
   * denotes the sums of the weights of the selected items for each.
   */
  [[nodiscard]] constexpr auto constraint_vector() const -> constraint_vector_type const & {
    return m_cvec;
  }

  /**
   * Flips a decision variable to 1 (i.e. to true), which effectively
   * adds this item to the knapsack. If the decision variable is already
   * 1, then this function does nothing.
   */
  constexpr auto flip_to_one(decision_vector_size_type i) {
    if (!m_dvec[i]) {
      this->flip_to_one_unchecked(i);
    }
  }

  /**
   * Flips a decision variable to 1 (i.e. to true), which effectively
   * adds this item to the knapsack. This variant does not check whether
   * the variable is already 1, instead, it assumes that the value of
   * the decision variable was 0.
   *
   * Undefined behavior if the decision variable value is already 1.
   */
  constexpr auto flip_to_one_unchecked(decision_vector_size_type i) {
    m_dvec[i] = 1;
    std::transform(m_ovec.begin(), m_ovec.end(), m_problem.get().item_values(i).begin(),
                   m_ovec.begin(), std::plus<>{});
    std::transform(m_cvec.begin(), m_cvec.end(), m_problem.get().item_weights(i).begin(),
                   m_cvec.begin(), std::plus<>{});
  }

  /**
   * Checks whether or not flipping a decision variable to 1 leads to a
   * feasible solution. If the decision variable is already 1, returns
   * whether the current solution is feasible.
   */
  [[nodiscard]] constexpr auto flip_to_one_feasible(decision_vector_size_type i) const {
    if (m_dvec[i]) {
      return this->feasible();
    } else {
      return this->flip_to_one_feasible_unchecked(i);
    }
  };

  /**
   * Checks whether or not flipping a decision variable to 1 leads to a
   * feasible solution. This variant does not check if the decision
   * variable is currently 0 or 1, instead, it assumes that it is 0.
   *
   * Undefined behavior if the decision variable is already 1.
   */
  [[nodiscard]] constexpr auto flip_to_one_feasible_unchecked(decision_vector_size_type i) const {
    auto iter1 = m_cvec.begin();
    auto last1 = m_cvec.end();
    auto iter2 = m_problem.get().item_weights(i).begin();
    auto iter3 = m_problem.get().weight_capacities().begin();
    while (iter1 != last1) {
      if (*(iter1++) + *(iter2++) > *(iter3++)) {
        return false;
      }
    }
    return true;
  };

  /**
   * Flips a decision variable to 0, which effectively removes this item
   * from the knapsack. If the decision variable is already 0, then this
   * function does nothing.
   */
  constexpr auto flip_to_zero(decision_vector_size_type i) {
    if (!m_dvec[i]) {
      this->flip_to_zero_unchecked(i);
    }
  }

  /**
   * Flips a decision variable to 0, which effectively removes this item
   * from the knapsack. This variant does not check whether the current
   * variable is already 0, instead it assumes that the value is 1.
   *
   * Undefined behavior if the decision variable value is already 0.
   */
  constexpr auto flip_to_zero_unchecked(decision_vector_size_type i) {
    m_dvec[i] = 0;
    std::transform(m_ovec.begin(), m_ovec.end(), m_problem.get().item_values(i).begin(),
                   m_ovec.begin(), std::minus<>{});
    std::transform(m_cvec.begin(), m_cvec.end(), m_problem.get().item_weights(i).begin(),
                   m_cvec.begin(), std::minus<>{});
  }

  /**
   * Checks whether or not flipping a decision variable to 0 leads to a
   * feasible solution. If the decision variable is already 0, returns
   * whether the current solution is feasible.
   */
  [[nodiscard]] constexpr auto flip_to_zero_feasible(decision_vector_size_type i) const {
    if (m_dvec[i]) {
      return this->feasible();
    } else {
      return this->flip_to_zero_feasible_unchecked(i);
    }
  };

  /**
   * Checks whether or not flipping a decision variable to 0 leads to a
   * feasible solution. This variant does not check if the decision
   * variable is currently 0 or 1, instead, it assumes that it is 1.
   *
   * Undefined behavior if the decision variable is already 0.
   */
  [[nodiscard]] constexpr auto flip_to_zero_feasible_unchecked(decision_vector_size_type i) const {
    auto iter1 = m_cvec.begin();
    auto last1 = m_cvec.end();
    auto iter2 = m_problem.get().item_weights(i).begin();
    auto iter3 = m_problem.get().weight_capacities().begin();
    while (iter1 != last1) {
      if (*(iter1++) - *(iter2++) > *(iter3++)) {
        return false;
      }
    }
    return true;
  };

  /**
   * Flips the value of a variable.
   *
   * @param i index of the variable
   */
  constexpr auto flip(decision_vector_size_type i) {
    if (this->m_dvec[i]) {
      this->flip_to_zero_unchecked(i);
    } else {
      this->flip_to_one_unchecked(i);
    }
  }

  /**
   * Check whether flipping the value of a variable leads to a feasible
   * solution or not.
   *
   * @param i index of the variable.
   * @return true if the resulting solution is feasible, false otherwise.
   */
  [[nodiscard]] constexpr auto flip_feasible(decision_vector_size_type i) const {
    if (this->m_dvec[i]) {
      return this->flip_to_one_feasible_unchecked(i);
    } else {
      return this->flip_to_zero_feasible_unchecked(i);
    }
  };

  /**
   * Exchange the values of two decision variables in the solution to
   * the knapsack problem. If the two values are the same no operation
   * is performed.
   *
   * @param i index to the first decision variable
   * @param j index to the second decision variable
   */
  constexpr auto exchange(decision_vector_size_type i, decision_vector_size_type j) {
    if (m_dvec[i] != m_dvec[j]) {
      this->exchange_unchecked(i, j);
    }
  }

  /**
   * Exchange the values of two decision variables in the solution to
   * the knapsack problem without checking if the associated values are
   * the same. Undefined behavior if the two values are the same.
   *
   * @param i index to the first decision variable
   * @param j index to the second decision variable
   */
  constexpr auto exchange_unchecked(decision_vector_size_type i, decision_vector_size_type j) {
    if (m_dvec[i]) {
      this->flip_to_zero_unchecked(i);
      this->flip_to_one_unchecked(j);
    } else {
      this->flip_to_one_unchecked(i);
      this->flip_to_zero_unchecked(j);
    }
  }

  /**
   * Returns whether or not exchanging two variables leads to a feasible
   * solution. This function assumes the values of the two variables are
   * distinct. Otherwise, it is undefined behavior.
   *
   * Time complexity: linear in the number of constraints.
   *
   * Memory complexity: constant.
   *
   * @param i index of the first variable to be exchanged
   * @param j index of the second variable to be exchanged
   * @return true if resulting solution is feasible, false otherwise
   */
  [[nodiscard]] constexpr auto exchange_feasible_unchecked(decision_vector_size_type i,
                                                           decision_vector_size_type j) const {
    auto iter1 = m_cvec.begin();
    auto last1 = m_cvec.end();
    auto iter2 = m_problem.get().item_weights(i).begin();
    auto iter3 = m_problem.get().item_weights(j).begin();
    auto iter4 = m_problem.get().weight_capacities().begin();
    if (m_dvec[i]) {
      while (iter1 != last1) {
        if (*(iter1++) - *(iter2++) + *(iter3++) > *(iter4++)) {
          return false;
        }
      }
      return true;
    } else {
      while (iter1 != last1) {
        if (*(iter1++) + *(iter2++) - *(iter3++) > *(iter4++)) {
          return false;
        }
      }
      return true;
    }
  }

  /**
   * Check if the current solutions is feasible or not.
   *
   * Time complexity: linear in the number of constraints.
   *
   * Memory complexity: constant.
   *
   * @return true if solution is feasible, false otherwise.
   */
  [[nodiscard]] constexpr auto feasible() const {
    auto iter1 = m_cvec.begin();
    auto last1 = m_cvec.end();
    auto iter2 = m_problem.get().weight_capacities().begin();
    while (iter1 != last1) {
      if (*(iter1++) > *(iter2++)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] constexpr auto problem() const {
    return m_problem;
  }

  [[nodiscard]] constexpr auto operator==(solution const &other) const {
    return decision_vector() == other.decision_vector();
  }

 private:
  std::reference_wrapper<const problem_type> m_problem;
  decision_vector_type m_dvec;
  objective_vector_type m_ovec;
  constraint_vector_type m_cvec;
};

}  // namespace mobkp

#endif
