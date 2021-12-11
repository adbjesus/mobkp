#ifndef MOBKP_ORDERS_HPP_
#define MOBKP_ORDERS_HPP_

#include <vector>

namespace mobkp {

template <typename Problem>
[[nodiscard]] constexpr auto objectives_orders(Problem const& problem) {
  std::vector<std::vector<size_t>> orders(problem.num_objectives());
  for (size_t i = 0; i < problem.num_objectives(); ++i) {
    for (size_t j = 0; j < problem.num_items(); ++j) {
      orders[i].push_back(j);
    }
    std::sort(orders[i].begin(), orders[i].end(), [&problem, i](auto const& lhs, auto const& rhs) {
      auto r1 = problem.item_value(lhs, i) * problem.item_weight(rhs, 0);
      auto r2 = problem.item_value(rhs, i) * problem.item_weight(lhs, 0);
      if (r1 > r2) {
        return true;
      } else if (r1 == r2 && problem.item_value(lhs, i) < problem.item_value(rhs, i)) {
        return true;
      }
      return false;
    });
  }
  return orders;
}

template <typename Problem, typename ObjectivesOrders>
[[nodiscard]] constexpr auto rank_max_order(Problem const& problem, ObjectivesOrders const& orders) {
  std::vector<size_t> rank_max(problem.num_items(), 0);
  std::vector<size_t> rank_sum(problem.num_items(), 0);
  for (size_t i = 0; i < problem.num_objectives(); ++i) {
    for (size_t j = 0; j < problem.num_items(); ++j) {
      rank_sum[orders[i][j]] += j + 1;
      rank_max[orders[i][j]] = std::max(rank_max[orders[i][j]], j + 1);
    }
  }
  std::vector<size_t> order_rank_max(problem.num_items());
  for (size_t j = 0; j < problem.num_items(); ++j) {
    order_rank_max[j] = j;
  }
  std::sort(order_rank_max.begin(), order_rank_max.end(), [&rank_max, &rank_sum](auto const& lhs, auto const& rhs) {
    return rank_max[lhs] < rank_max[rhs] || (rank_max[lhs] == rank_max[rhs] && rank_sum[lhs] < rank_sum[rhs]);
  });
  return order_rank_max;
}

template <typename Problem, typename ObjectivesOrders>
[[nodiscard]] constexpr auto rank_min_order(Problem const& problem, ObjectivesOrders const& orders) {
  std::vector<size_t> rank_min(problem.num_items(), 0);
  std::vector<size_t> rank_sum(problem.num_items(), 0);
  for (size_t i = 0; i < problem.num_objectives(); ++i) {
    for (size_t j = 0; j < problem.num_items(); ++j) {
      rank_sum[orders[i][j]] += j + 1;
      rank_min[orders[i][j]] = std::min(rank_min[orders[i][j]], j + 1);
    }
  }
  std::vector<size_t> order_rank_min(problem.num_items());
  for (size_t j = 0; j < problem.num_items(); ++j) {
    order_rank_min[j] = j;
  }
  std::sort(order_rank_min.begin(), order_rank_min.end(), [&rank_min, &rank_sum](auto const& lhs, auto const& rhs) {
    return rank_min[lhs] < rank_min[rhs] || (rank_min[lhs] == rank_min[rhs] && rank_sum[lhs] < rank_sum[rhs]);
  });
  return order_rank_min;
}

template <typename Problem, typename ObjectivesOrders>
[[nodiscard]] constexpr auto rank_sum_order(Problem const& problem, ObjectivesOrders const& orders) {
  std::vector<size_t> rank_sum(problem.num_items(), 0);
  for (size_t i = 0; i < problem.num_objectives(); ++i) {
    for (size_t j = 0; j < problem.num_items(); ++j) {
      rank_sum[orders[i][j]] += j + 1;
    }
  }
  std::vector<size_t> order_rank_sum(problem.num_items());
  for (size_t j = 0; j < problem.num_items(); ++j) {
    order_rank_sum[j] = j;
  }
  std::sort(order_rank_sum.begin(), order_rank_sum.end(),
            [&rank_sum](auto const& lhs, auto const& rhs) { return rank_sum[lhs] < rank_sum[rhs]; });
  return order_rank_sum;
}

}  // namespace mobkp

#endif  // MOBKP_ORDERS_HPP_
