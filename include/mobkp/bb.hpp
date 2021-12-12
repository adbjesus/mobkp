#ifndef MOBKP_BB_HPP_
#define MOBKP_BB_HPP_

#include "bounds.hpp"
#include "orders.hpp"

#include <mooutils/sets.hpp>

#include <chrono>
#include <utility>

namespace mobkp {

namespace bbdetails {

template <typename Solution>
struct node {
 public:
  using solution_type = Solution;

  node(size_t index, solution_type&& solution)
      : m_index(index)
      , m_solution(std::move(solution)) {}

  node(node const&) = delete;
  node(node&&) = default;
  node& operator=(node const&) = delete;
  node& operator=(node&&) = default;

  [[nodiscard]] auto solution() const -> solution_type const& {
    return m_solution;
  }

  [[nodiscard]] auto index() const {
    return m_index;
  }

  template <typename ObjectivesOrdersList, typename RatioSumOrderList>
  [[nodiscard]] auto lower_bound(ObjectivesOrdersList const& ool, RatioSumOrderList const& rsol) const {
    auto lb = std::vector<solution_type>{};
    if (m_index >= m_solution.problem().get().num_items()) {
      lb.emplace_back(m_solution);
    } else {
      auto no = m_solution.problem().get().num_objectives();
      lb.reserve(no + 1);
      for (size_t i = 0; i < no; ++i) {
        lb.emplace_back(mobkp::lower_bound(m_solution, ool[m_index][i]));
      }
      lb.emplace_back(mobkp::lower_bound(m_solution, rsol[m_index]));
    }
    return lb;
  }

  template <typename ObjectivesOrdersList>
  [[nodiscard]] auto upper_bound(ObjectivesOrdersList const& ool) const {
    if (m_index >= m_solution.problem().get().num_items()) {
      return m_solution.objective_vector();
    } else {
      return mobkp::upper_bound(m_solution, ool[m_index]);
    }
  }

  [[nodiscard]] auto branches() && {
    auto res = std::vector<node>{};
    if (m_index >= m_solution.problem().get().num_items()) {
      return res;
    } else if (m_solution.flip_to_one_feasible(m_index)) {
      auto sol1 = m_solution;
      auto sol2 = std::move(m_solution);
      sol2.flip_to_one_unchecked(m_index);
      res.reserve(2);
      res.emplace_back(node(m_index + 1, std::move(sol1)));
      res.emplace_back(node(m_index + 1, std::move(sol2)));
    } else {
      res.reserve(1);
      res.emplace_back(node(m_index + 1, std::move(m_solution)));
    }
    return res;
  }

 private:
  size_t m_index;
  solution_type m_solution;
};

template <typename Solution>
[[nodiscard]] auto branches(Solution&& s) {
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

  auto lb = mooutils::flat_set<solution_type>();
  auto iteration = static_cast<size_t>(0);

  auto objectives_orders = mobkp::objectives_orders(problem);
  auto ratio_sum_order = mobkp::ratio_sum_order(problem);

  auto objectives_orders_list = std::vector<decltype(objectives_orders)>();
  objectives_orders_list.reserve(problem.num_items());

  auto ratio_sum_order_list = std::vector<decltype(ratio_sum_order)>();
  ratio_sum_order_list.reserve(problem.num_items());

  auto remove_from_order = [](auto& order, size_t i) {
    std::erase(order, i);
  };

  auto remove_from_orders = [](auto& orders, size_t i) {
    for (auto& ord : orders) {
      std::erase(ord, i);
    }
  };

  for (size_t i = 0; i < problem.num_items(); ++i) {
    objectives_orders_list.emplace_back(objectives_orders);
    ratio_sum_order_list.emplace_back(ratio_sum_order);
    remove_from_orders(objectives_orders, i);
    remove_from_order(ratio_sum_order, i);
  }

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration, &objectives_orders_list,
                              &ratio_sum_order_list](auto& node) {
    for (auto& s : node.lower_bound(objectives_orders_list, ratio_sum_order_list)) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  auto empty = bbdetails::node<solution_type>(0, solution_type::empty(problem));
  update_lb_and_trace(empty);
  queue.push(std::move(empty));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    // if (iteration % 1000 == 0) {
    //   std::cerr << iteration << " " << queue.size() << " " << lb.size() << "\n";
    // }
    auto node = queue.pop();
    for (auto& branch : std::move(node).branches()) {
      if (!mooutils::strictly_dominates(lb, branch.upper_bound(objectives_orders_list))) {
        update_lb_and_trace(branch);
        queue.push(std::move(branch));
      }
    }
  }

  auto res = mooutils::unordered_set<solution_type>{};
  for (auto& s : lb) {
    res.insert_unchecked(std::move(s));
  }
  return res;
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

  auto objectives_orders = mobkp::objectives_orders(problem);
  auto ratio_sum_order = mobkp::ratio_sum_order(problem);

  auto objectives_orders_list = std::vector<decltype(objectives_orders)>();
  objectives_orders_list.reserve(problem.num_items());

  auto ratio_sum_order_list = std::vector<decltype(ratio_sum_order)>();
  ratio_sum_order_list.reserve(problem.num_items());

  auto remove_from_order = [](auto& order, size_t i) {
    std::erase(order, i);
  };

  auto remove_from_orders = [](auto& orders, size_t i) {
    for (auto& ord : orders) {
      std::erase(ord, i);
    }
  };

  for (size_t i = 0; i < problem.num_items(); ++i) {
    objectives_orders_list.emplace_back(objectives_orders);
    ratio_sum_order_list.emplace_back(ratio_sum_order);
    remove_from_orders(objectives_orders, i);
    remove_from_order(ratio_sum_order, i);
  }

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration, &objectives_orders_list,
                              &ratio_sum_order_list](auto& node) {
    for (auto& s : node.lower_bound(objectives_orders_list, ratio_sum_order_list)) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  queue.push(bbdetails::node(0, solution_type::empty(problem)));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    // if (iteration % 1000 == 0) {
    //   std::cerr << iteration << " " << queue.size() << " " << lb.size() << "\n";
    // }
    auto node = queue.pop();
    if (!mooutils::strictly_dominates(lb, node.upper_bound(objectives_orders_list))) {
      update_lb_and_trace(node);
      for (auto& branch : std::move(node).branches()) {
        queue.push(std::move(branch));
      }
    }
  }

  auto res = mooutils::unordered_set<solution_type>{};
  for (auto& s : lb) {
    res.insert_unchecked(std::move(s));
  }
  return res;
}

}  // namespace mobkp

#endif
