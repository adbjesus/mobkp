#ifndef MOBKP_BB_HPP_
#define MOBKP_BB_HPP_

#include <mooutils/sets.hpp>

#include <chrono>
#include <utility>

namespace mobkp {

namespace bbdetails {

template <typename Solution>
[[nodiscard]] auto branches(Solution&& s) {
  throw("Unimplemented!");
}

template <typename Solution>
[[nodiscard]] auto lower_bound(Solution const& s) {
  throw("Unimplemented!");
}

template <typename Solution>
[[nodiscard]] auto upper_bound(Solution const& s) {
  throw("Unimplemented!");
}

}  // namespace bbdetails

template <typename Solution, typename Problem, typename Queue, typename AnytimeTrace>
[[nodiscard]] constexpr auto eager_branch_and_bound(Problem const& problem, Queue& queue, AnytimeTrace& anytime_trace,
                                                    double timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto lb = mooutils::set<solution_type>();
  auto iteration = static_cast<size_t>(0);

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration](auto& sol) {
    for (auto& s : bbdetails::lower_bound(sol)) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  auto empty = solution_type::empty(problem);
  update_lb_and_trace(empty);
  queue.push(std::move(empty));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    auto node = queue.pop();
    for (auto& branch : bbdetails::branches(std::move(node))) {
      if (!mooutils::strictly_dominates(lb, bbdetails::upper_bound(branch))) {
        update_lb_and_trace(branch);
        queue.push(std::move(branch));
      }
    }
  }

  return lb;
}

template <typename Solution, typename Problem, typename Queue, typename AnytimeTrace>
[[nodiscard]] constexpr auto lazy_branch_and_bound(Problem const& problem, Queue& queue, AnytimeTrace& anytime_trace,
                                                   double timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto lb = mooutils::set<solution_type>();
  auto iteration = static_cast<size_t>(0);

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration](auto& node) {
    for (auto& sol : bbdetails::lower_bound(node)) {
      auto it = lb.insert(sol);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, sol);
      }
    }
  };

  queue.push(solution_type::empty(problem));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    auto node = queue.pop();
    if (!mooutils::strictly_dominates(lb, bbdetails::upper_bound(node))) {
      update_lb_and_trace(node);
      for (auto& branch : bbdetails::branches(std::move(node))) {
        queue.push(std::move(branch));
      }
    }
  }

  return lb;
}

}  // namespace mobkp

#endif
