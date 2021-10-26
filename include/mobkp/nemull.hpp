#pragma once

#include <mobkp/problem.hpp>

#include <mooutils/orders.hpp>

#include <chrono>

namespace mobkp {

template <typename Problem, typename AnytimeTrace>
auto nemull(Problem const& problem, AnytimeTrace& anytime_trace, double const timeout) {
  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  using data_type = typename Problem::data_type;
  using dvec_type = std::vector<bool>;
  using ovec_type = std::vector<data_type>;
  using cvec_type = std::vector<data_type>;
  using solution_type = mobkp::solution<Problem, dvec_type, ovec_type, cvec_type>;

  auto sols = std::vector<solution_type>{};
  // TODO add initial solution

  size_t iteration = 0;
  for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
    auto aux = std::vector<solution_type>{};

    auto s = queue.pop();

    if (std::any_of(solutions.begin(), solutions.end(),
                    [&s](auto const& sol) { return mooutils::dominates(sol, s); })) {
      continue;
    }
    if (mooutils::strictly_dominates(solutions, s)) {
      continue;
    }

    // 1-flip neighborhood
    bool flip = false;
    for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
      ++iteration;

      if (s.decision_vector()[i]) {
        continue;
      }

      if (!s.flip_to_one_feasible_unchecked(i)) {
        continue;
      }

      auto new_sol = s;
      new_sol.flip_to_one_unchecked(i);
      auto nst = solutions.insert(std::move(new_sol));
      if (nst != solutions.end()) {
        flip = true;
        queue.push(*nst);
        anytime_trace.add_solution(iteration, *nst);
      }
    }

    if (flip) {
      continue;
    }

    // Exchange neighborhood
    for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
      for (size_t j = i + 1; j < problem.num_items() && elapsed_sec() < timeout; ++j) {
        iteration += 1;

        if (s.decision_vector()[i] == s.decision_vector()[j]) {
          continue;
        }

        if (!s.exchange_feasible_unchecked(i, j)) {
          continue;
        }

        auto new_sol = s;
        new_sol.exchange(i, j);
        auto nst = solutions.insert(std::move(new_sol));
        if (nst != solutions.end()) {
          queue.push(*nst);
          anytime_trace.add_solution(iteration, *nst);
          goto loopend;
        }
      }
    }
  loopend:;
  }
}

}  // namespace mobkp
