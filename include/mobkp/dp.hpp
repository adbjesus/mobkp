#ifndef MOBKP_DP_HPP_
#define MOBKP_DP_HPP_

#include "bounds.hpp"
#include "orders.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <mooutils/orders.hpp>
#include <mooutils/sets.hpp>

#include <algorithm>
#include <chrono>
#include <type_traits>
#include <utility>
#include <vector>

namespace mobkp {

namespace dpdetails {

template <typename SolutionSet>
constexpr auto build_aux(SolutionSet const& sols, size_t i) {
  using solution_type = typename SolutionSet::value_type;
  auto aux = std::vector<solution_type>{};
  aux.reserve(sols.size());
  for (auto const& s : sols) {
    if (s.flip_to_one_feasible(i)) {
      auto tmp = s;
      tmp.flip_to_one_unchecked(i);
      aux.push_back(std::move(tmp));
    }
  }
  return aux;
}

constexpr auto dominates_with_weight(auto const& s1, auto const& s2) {
  int eq = 1;
  auto w1iter = s1.constraint_vector().begin();
  auto w1last = s1.constraint_vector().end();
  auto w2iter = s2.constraint_vector().begin();
  for (; w1iter != w1last; ++w1iter, ++w2iter) {
    if (*w1iter > *w2iter) {
      return 0;
    }
    if (*w1iter < *w2iter) {
      eq = 2;
    }
  }
  auto o1iter = s1.objective_vector().begin();
  auto o1last = s1.objective_vector().end();
  auto o2iter = s2.objective_vector().begin();
  for (; o1iter != o1last; ++o1iter, ++o2iter) {
    if (*o1iter < *o2iter) {
      return 0;
    }
    if (*o1iter > *o2iter) {
      eq = 2;
    }
  }
  return eq;
}

constexpr auto merge(auto& sols, auto&& aux, auto& anytime_trace, auto& iteration) {
  // Returns:
  //  - 0 if s1 does not (weakly) dominate s2
  //  - 1 if s1 is equivalent to s2
  //  - 2 if s1 dominates s2

  sols.reserve(sols.size() + aux.size());
  for (size_t j = 0, l = sols.size(); j < aux.size(); ++j) {
    ++iteration;
    size_t k = 0;
    for (; k < l; ++k) {
      if (dominates_with_weight(sols[k], aux[j]) == 2) {
        break;
      }

      auto c = dominates_with_weight(aux[j], sols[k]);
      if (c == 1) {
        break;
      } else if (c == 2) {
        if (k != l - 1) {
          sols[k] = std::move(sols[l - 1]);
        }
        if (l != sols.size()) {
          sols[l - 1] = std::move(sols.back());
        }
        sols.pop_back();
        --l;
        for (; k < l;) {
          if (dominates_with_weight(aux[j], sols[k]) == 2) {
            if (k != l - 1) {
              sols[k] = std::move(sols[l - 1]);
            }
            if (l != sols.size()) {
              sols[l - 1] = std::move(sols.back());
            }
            sols.pop_back();
            --l;
          } else {
            ++k;
          }
        }
        break;
      }
    }
    if (k == l) {
      anytime_trace.add_solution(iteration, aux[j]);
      sols.push_back(std::move(aux[j]));
    }
  }
}

template <typename SolutionSet>
[[nodiscard]] constexpr auto filter_last(SolutionSet&& sols) {
  using solution_type = SolutionSet::value_type;
  auto res = mooutils::unordered_set<solution_type>{};
  for (auto&& s : sols) {
    res.insert(std::move(s));
  }
  return res;
}

[[nodiscard]] constexpr auto lex_ge(auto const& s1, auto const& s2) {
  if (s1.constraint_vector()[0] < s2.constraint_vector()[0]) {
    return true;
  }
  if (s1.constraint_vector()[0] == s2.constraint_vector()[0]) {
    return !mooutils::lexicographically_less(s1, s2);
  }
  return false;
}

[[nodiscard]] constexpr auto lex_wdom(auto const& s1, auto const& s2) {
  return !mooutils::lexicographically_less(s1, s2);
}

template <typename S, typename M, typename C>
auto mantain_non_dominated(S&& s, M& m, C& c) {
  if (m.insert(s) != m.end()) {
    c.push_back(std::forward<S>(s));
  }
}

}  // namespace dpdetails

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] constexpr auto nemull_dp(Problem const& problem, AnytimeTrace& anytime_trace, double const timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto sols = std::vector<solution_type>{solution_type::empty(problem)};

  size_t iteration = 0;
  for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
    // std::cerr << i << " " << sols.size() << "\n";
    dpdetails::merge(sols, dpdetails::build_aux(sols, i), anytime_trace, iteration);
  }

  return dpdetails::filter_last(std::move(sols));
}

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] constexpr auto bhv_dp(Problem const& problem, AnytimeTrace& anytime_trace, double const timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto remove_from_order = [](auto& order, size_t i) {
    std::erase(order, i);
  };

  auto remove_from_orders = [](auto& orders, size_t i) {
    for (auto& ord : orders) {
      std::erase(ord, i);
    }
  };

  auto orders = mobkp::objectives_orders(problem);
  auto order_max = mobkp::rank_max_order(problem, orders);
  auto order_sum = mobkp::rank_sum_order(problem, orders);

  auto sols = std::vector<solution_type>{solution_type::empty(problem)};
  anytime_trace.add_solution(0, *sols.begin());

  auto const& pc = problem.weight_capacities();
  auto rw = std::vector<typename Problem::data_type>(problem.num_constraints(), 0);
  for (size_t i = 0; i < problem.num_items(); ++i) {
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter += *iciter;
    }
  }

  for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
    // std::cerr << ind << " " << sols.size() << "\n";
    remove_from_orders(orders, i);
    remove_from_order(order_sum, i);
    remove_from_order(order_max, i);

    auto aux = std::vector<solution_type>{};
    aux.reserve(sols.size() * 2);
    auto M = mooutils::set<solution_type>{};

    size_t j = 0;
    size_t k = 0;
    for (; k < sols.size() && sols[k].constraint_vector()[0] + rw[0] <= pc[0]; ++k)
      ;
    for (; j < sols.size() && sols[j].flip_to_one_feasible(i); ++j) {
      auto tmp = sols[j];
      tmp.flip_to_one_unchecked(i);
      for (; k < sols.size() && dpdetails::lex_ge(sols[k], tmp); ++k) {
        if (j < k) {
          dpdetails::mantain_non_dominated(sols[k], M, aux);
        } else {
          dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
        }
      }
      anytime_trace.add_solution(0, tmp);
      dpdetails::mantain_non_dominated(std::move(tmp), M, aux);
    }
    for (; k < sols.size(); ++k) {
      dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
    }

    auto lb = mooutils::set<solution_type>();
    for (auto&& s : M) {
      lb.insert(mobkp::lower_bound(s, order_sum));
      lb.insert(mobkp::lower_bound(std::move(s), order_max));
    }

    // Dom3
    for (auto it = aux.begin(); it != aux.end(); ++it) {
      if (!mooutils::strictly_dominates(lb, mobkp::upper_bound(*it, orders))) {
        sols = decltype(sols)(std::make_move_iterator(it), std::make_move_iterator(aux.end()));
        break;
      }
    }

    // Update rw
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter -= *iciter;
    }
  }

  return dpdetails::filter_last(std::move(sols));
  // return lb;
}

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] constexpr auto fpsv_dp(Problem const& problem, AnytimeTrace& anytime_trace, double const timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto remove_from_order = [](auto& order, size_t i) {
    std::erase(order, i);
  };

  auto remove_from_orders = [](auto& orders, size_t i) {
    for (auto& ord : orders) {
      std::erase(ord, i);
    }
  };

  if (problem.num_objectives() > 2) {
    throw("FPSV-DP not implemented for more than 2 objectives");
  }

  auto lb = mooutils::set<solution_type>();
  // auto fk_at = fake_anytime_trace();
  for (auto&& s : dws<solution_type>(problem, anytime_trace, timeout)) {
    lb.insert_unchecked(std::move(s));
  }

  auto orders = mobkp::objectives_orders(problem);
  auto order_max = mobkp::rank_max_order(problem, orders);
  auto order_sum = mobkp::rank_sum_order(problem, orders);

  auto sols = std::vector<solution_type>{solution_type::empty(problem)};
  if (lb.insert(*sols.begin()) != lb.end())
    anytime_trace.add_solution(0, *sols.begin());

  auto const& pc = problem.weight_capacities();
  auto rw = std::vector<typename Problem::data_type>(problem.num_constraints(), 0);
  for (size_t i = 0; i < problem.num_items(); ++i) {
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter += *iciter;
    }
  }

  for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
    // std::cerr << ind << " " << sols.size() << "\n";
    remove_from_orders(orders, i);
    remove_from_order(order_sum, i);
    remove_from_order(order_max, i);

    auto aux = std::vector<solution_type>{};
    aux.reserve(sols.size() * 2);
    auto M = mooutils::set<solution_type>{};

    size_t j = 0;
    size_t k = 0;
    for (; k < sols.size() && sols[k].constraint_vector()[0] + rw[0] <= pc[0]; ++k)
      ;
    for (; j < sols.size() && sols[j].flip_to_one_feasible(i); ++j) {
      auto tmp = sols[j];
      tmp.flip_to_one_unchecked(i);
      for (; k < sols.size() && dpdetails::lex_ge(sols[k], tmp); ++k) {
        if (j < k) {
          dpdetails::mantain_non_dominated(sols[k], M, aux);
        } else {
          dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
        }
      }
      if (lb.insert(tmp) != lb.end())
        anytime_trace.add_solution(0, tmp);
      dpdetails::mantain_non_dominated(std::move(tmp), M, aux);
    }
    for (; k < sols.size(); ++k) {
      dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
    }

    // Dom3 - BDP-1 variant
    for (auto it = aux.begin(); it != aux.end(); ++it) {
      if (!mooutils::strictly_dominates(lb, mobkp::upper_bound(*it, orders))) {
        sols = decltype(sols)(std::make_move_iterator(it), std::make_move_iterator(aux.end()));
        break;
      }
    }

    // Update rw
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter -= *iciter;
    }
  }

  auto res = mooutils::unordered_set<solution_type>();
  for (auto&& s : lb)
    res.insert_unchecked(std::move(s));
  return res;
  // return dpdetails::filter_last(std::move(sols));
  // return lb;
}

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] constexpr auto anytime_dp(Problem const& problem, AnytimeTrace& anytime_trace, double const timeout) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto remove_from_order = [](auto& order, size_t i) {
    std::erase(order, i);
  };

  auto remove_from_orders = [](auto& orders, size_t i) {
    for (auto& ord : orders) {
      std::erase(ord, i);
    }
  };

  auto lb = mooutils::set<solution_type>();

  if (problem.num_objectives() == 2) {
    for (auto&& s : dws<solution_type>(problem, anytime_trace, timeout)) {
      lb.insert_unchecked(std::move(s));
    }
  }

  auto orders = mobkp::objectives_orders(problem);
  auto order_max = mobkp::rank_max_order(problem, orders);
  auto order_sum = mobkp::rank_sum_order(problem, orders);

  auto sols = std::vector<solution_type>{solution_type::empty(problem)};
  anytime_trace.add_solution(0, *sols.begin());

  auto add_to_lb = [&lb, &orders, &order_sum, &order_max, &anytime_trace](auto const& s) {
    for (auto ord : orders) {
      auto it = lb.insert(mobkp::lower_bound(s, ord));
      if (it != lb.end()) {
        anytime_trace.add_solution(0, *it);
      }
    }
    {
      auto it = lb.insert(mobkp::lower_bound(s, order_sum));
      if (it != lb.end()) {
        anytime_trace.add_solution(0, *it);
      }
    }
    {
      auto it = lb.insert(mobkp::lower_bound(s, order_max));
      if (it != lb.end()) {
        anytime_trace.add_solution(0, *it);
      }
    }
  };

  add_to_lb(*sols.begin());

  auto const& pc = problem.weight_capacities();
  auto rw = std::vector<typename Problem::data_type>(problem.num_constraints(), 0);
  for (size_t i = 0; i < problem.num_items(); ++i) {
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter += *iciter;
    }
  }

  for (size_t i = 0; i < problem.num_items() && elapsed_sec() < timeout; ++i) {
    // std::cerr << ind << " " << sols.size() << "\n";
    remove_from_orders(orders, i);
    remove_from_order(order_sum, i);
    remove_from_order(order_max, i);

    auto aux = std::vector<solution_type>{};
    aux.reserve(sols.size() * 2);
    auto M = mooutils::set<solution_type>{};

    size_t j = 0;
    size_t k = 0;
    for (; k < sols.size() && sols[k].constraint_vector()[0] + rw[0] <= pc[0]; ++k)
      ;
    for (; j < sols.size() && sols[j].flip_to_one_feasible(i); ++j) {
      auto tmp = sols[j];
      tmp.flip_to_one_unchecked(i);
      for (; k < sols.size() && dpdetails::lex_ge(sols[k], tmp); ++k) {
        if (j < k) {
          dpdetails::mantain_non_dominated(sols[k], M, aux);
        } else {
          dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
        }
      }
      add_to_lb(tmp);
      dpdetails::mantain_non_dominated(std::move(tmp), M, aux);
    }
    for (; k < sols.size(); ++k) {
      dpdetails::mantain_non_dominated(std::move(sols[k]), M, aux);
    }

    // Dom3
    auto rend = std::remove_if(aux.begin(), aux.end(), [&lb, &orders](auto const& s) {
      return mooutils::strictly_dominates(lb, mobkp::upper_bound(s, orders));
    });
    aux.erase(rend, aux.end());
    sols = std::move(aux);

    // Update rw
    auto rwiter = rw.begin();
    auto rwlast = rw.end();
    auto iciter = problem.item_weights(i).begin();
    for (; rwiter != rwlast; ++rwiter, ++iciter) {
      *rwiter -= *iciter;
    }
  }

  auto res = mooutils::unordered_set<solution_type>();
  for (auto& s : lb)
    res.insert_unchecked(s);
  return res;
}

}  // namespace mobkp

#endif
