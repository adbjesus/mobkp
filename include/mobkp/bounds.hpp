#ifndef MOBKP_BOUNDS_HPP_
#define MOBKP_BOUNDS_HPP_

#include <utility>
#include <vector>

namespace mobkp {

template <typename Solution, typename Order>
constexpr auto lower_bound(Solution&& s, Order const& order) {
  auto lb = std::forward<Solution>(s);
  for (auto i : order) {
    if (lb.flip_to_one_feasible(i)) {
      lb.flip_to_one_unchecked(i);
    } else {
      break;
    }
  }
  return lb;
}

// This gives the Martello and Toth bound or Dantzig's as a fallback
template <typename Solution, typename ObjectivesOrders>
[[nodiscard]] constexpr auto upper_bound(Solution const& s, ObjectivesOrders const& orders) {
  using data_type = typename std::remove_cvref_t<decltype(s.objective_vector())>::value_type;
  auto ub = std::vector<data_type>(s.objective_vector().size(), 0);
  for (size_t i = 0; i < ub.size(); ++i) {
    auto tmp = s;
    auto rem = data_type(0);
    for (size_t j = 0; j < orders[i].size(); ++j) {
      auto k = orders[i][j];
      if (tmp.flip_to_one_feasible(k)) {
        tmp.flip_to_one_unchecked(k);
      } else {
        // k is the breaking item
        auto p = tmp.problem();
        auto vk = p.get().item_value(k, i);
        auto wk = p.get().item_weight(k, 0);
        auto remw = p.get().weight_capacity(0) - tmp.constraint_vector()[0];
        if (j > 0 && j + 1 < orders[i].size()) {
          // This gives Martello and Toth remainder
          auto lasti = orders[i][j - 1];
          auto nexti = orders[i][j + 1];
          auto b1 = (remw * p.get().item_value(nexti, i)) / p.get().item_weight(nexti, 0);
          auto b2 = vk + (remw * p.get().item_value(lasti, i)) / p.get().item_weight(lasti, 0);
          rem = std::max(b1, b2);
        } else {
          // This gives Dantzig's remainder
          rem = (remw * vk) / wk + 1;
        }
        break;
      }
    }
    ub[i] = tmp.objective_vector()[i] + rem;
  }
  return ub;
}

}  // namespace mobkp

#endif  // MOBKP_BOUNDS_HPP_
