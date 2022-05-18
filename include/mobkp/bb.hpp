#ifndef MOBKP_BB_HPP_
#define MOBKP_BB_HPP_

#include "bounds.hpp"
#include "orders.hpp"

#include <mooutils/queues.hpp>
#include <mooutils/sets.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <ranges>
#include <utility>
#include <variant>

namespace mobkp {

namespace bbdetails {

template <typename Solution, typename ObjectiveOrdersList = std::vector<std::vector<std::vector<size_t>>>,
          typename RatioSumOrderList = std::vector<std::vector<size_t>>>
struct node {
 public:
  using solution_type = Solution;
  using objective_orders_list = ObjectiveOrdersList;
  using ratio_sum_order_list = RatioSumOrderList;

  node(size_t index, solution_type&& solution, objective_orders_list const& ool, RatioSumOrderList const& rsol)
      : m_index(index)
      , m_solution(std::move(solution))
      , m_ool(ool)
      , m_rsol(rsol) {}

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

  [[nodiscard]] auto lower_bound() const {
    auto lb = std::vector<solution_type>{};
    if (m_index >= m_solution.problem().get().num_items()) {
      lb.emplace_back(m_solution);
    } else {
      auto no = m_solution.problem().get().num_objectives();
      lb.reserve(no + 1);
      for (size_t i = 0; i < no; ++i) {
        lb.emplace_back(mobkp::lower_bound(m_solution, m_ool.get()[m_index][i]));
      }
      lb.emplace_back(mobkp::lower_bound(m_solution, m_rsol.get()[m_index]));
    }
    return lb;
  }

  [[nodiscard]] auto upper_bound() const {
    if (m_index >= m_solution.problem().get().num_items()) {
      return m_solution.objective_vector();
    } else {
      return mobkp::upper_bound(m_solution, m_ool.get()[m_index]);
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
      res.emplace_back(node(m_index + 1, std::move(sol1), m_ool.get(), m_rsol.get()));
      res.emplace_back(node(m_index + 1, std::move(sol2), m_ool.get(), m_rsol.get()));
    } else {
      res.reserve(1);
      res.emplace_back(node(m_index + 1, std::move(m_solution), m_ool.get(), m_rsol.get()));
    }
    return res;
  }

 private:
  size_t m_index;
  solution_type m_solution;
  std::reference_wrapper<const objective_orders_list> m_ool;
  std::reference_wrapper<const ratio_sum_order_list> m_rsol;
};

template <typename Node, typename Hypervolume,
          typename Container = std::vector<std::tuple<Node, typename Hypervolume::value_type, std::size_t>>>
class hypervolume_befs_queue
    : public mooutils::base_queue<hypervolume_befs_queue<Node, Hypervolume, Container>, Node, Container> {
 private:
  struct HeapCmp {
    auto operator()(auto const& a, auto const& b) const -> bool {
      return std::get<1>(a) < std::get<1>(b);
    }
  };

 public:
  using node_type = Node;
  using hypervolume_type = Hypervolume;
  using hypervolume_value_type = typename hypervolume_type::value_type;

  template <typename... Args>
  explicit hypervolume_befs_queue(Args&&... args, hypervolume_type&& hv_object)
      : base_class_type(std::forward<Args>(args)...)
      , m_hv_object(std::move(hv_object)) {}

  [[nodiscard]] constexpr auto nodes() {
    return std::views::all(this->c) | std::views::transform([](auto& v) -> node_type& { return std::get<0>(v); });
  }

 private:
  using base_class_type = mooutils::base_queue<hypervolume_befs_queue<Node, Hypervolume, Container>, Node, Container>;
  using typename base_class_type::size_type;
  using typename base_class_type::value_type;

  friend base_class_type;

  // Push a new solution to the queue
  template <typename T>
  constexpr auto push_impl(T&& node) {
    m_no_update_counter += 1;
    for (auto const& v : node.lower_bound()) {
      if (m_hv_object.insert(v.objective_vector()) > 0) {
        m_no_update_counter = 0;
      }
    }

    if (m_no_update_counter == 0) {
      ++m_counter;
      m_container_is_heap = false;
      auto hvc = m_hv_object.contribution(node.upper_bound());
      this->c.emplace_back(std::forward<T>(node), hvc, m_counter);
    } else if (m_no_update_counter == 10) {
      m_container_is_heap = true;
      for (auto& [node, hvc, lastupdated] : this->c) {
        if (lastupdated < m_counter) {
          hvc = m_hv_object.contribution(node.upper_bound());
          lastupdated = m_counter;
        }
      }
      auto hvc = m_hv_object.contribution(node.upper_bound());
      this->c.emplace_back(std::forward<T>(node), hvc, m_counter);
      std::ranges::make_heap(this->c, HeapCmp());
    } else if (m_no_update_counter < 10) {
      auto hvc = m_hv_object.contribution(node.upper_bound());
      this->c.emplace_back(std::forward<T>(node), hvc, m_counter);
    } else {  // > 10
      auto hvc = m_hv_object.contribution(node.upper_bound());
      this->c.emplace_back(std::forward<T>(node), hvc, m_counter);
      std::ranges::push_heap(this->c, HeapCmp());
    }
  }

  constexpr auto pop_impl() -> value_type {
    // Find largest
    auto best_it = this->c.begin();
    if (!m_container_is_heap) {
      for (auto it = std::next(this->c.begin()); it != this->c.end(); ++it) {
        if (std::get<1>(*it) > std::get<1>(*best_it)) {
          best_it = it;
        }
      }
    }

    if (m_counter > std::get<2>(*best_it)) {
      // Try to only update the first value
      auto& [nodef, hvcf, lastupdatedf] = *best_it;
      hvcf = m_hv_object.contribution(nodef.upper_bound());
      lastupdatedf = m_counter;

      auto best_hvc = hvcf;
      for (auto it = this->c.begin(); it != this->c.end(); ++it) {
        if (m_counter > std::get<2>(*it) && std::get<1>(*it) > best_hvc) {
          auto& [node, hvc, lastupdated] = *it;
          hvc = m_hv_object.contribution(node.upper_bound());
          lastupdated = m_counter;
          if (hvc > best_hvc) {
            best_it = it;
            best_hvc = hvc;
          }
        }
      }
    }

    if (m_container_is_heap) {
      std::ranges::pop_heap(this->c, HeapCmp());
    } else {
      std::swap(*best_it, this->c.back());
    }

    auto res = std::move(std::get<0>(this->c.back()));
    this->c.pop_back();

    return res;
  }

  hypervolume_type m_hv_object;
  bool m_container_is_heap = false;
  std::size_t m_counter = 0;
  std::size_t m_no_update_counter = 0;
};

template <typename Node, typename Hypervolume, typename Container = std::vector<Node>>
class hypervolume_bedfs_queue
    : public mooutils::base_queue<hypervolume_bedfs_queue<Node, Hypervolume, Container>, Node, Container> {
 public:
  using node_type = Node;
  using hypervolume_type = Hypervolume;
  using hypervolume_value_type = typename hypervolume_type::value_type;

  template <typename... Args>
  explicit hypervolume_bedfs_queue(Args&&... args, hypervolume_type&& hv_object)
      : base_class_type(std::forward<Args>(args)...)
      , m_hv_object(std::move(hv_object)) {}

  [[nodiscard]] constexpr auto nodes() {
    return std::views::all(this->c);
  }

  [[nodiscard]] constexpr auto hv_object() -> hypervolume_type& {
    return m_hv_object;
  }

 private:
  using base_class_type = mooutils::base_queue<hypervolume_bedfs_queue<Node, Hypervolume, Container>, Node, Container>;
  using typename base_class_type::size_type;
  using typename base_class_type::value_type;

  friend base_class_type;

  // Push a new solution to the queue
  template <typename T>
  constexpr auto push_impl(T&& node) {
    for (auto const& v : node.lower_bound()) {
      m_hv_object.insert(v.objective_vector());
    }
    this->c.push_back(std::forward<T>(node));
  }

  constexpr auto pop_impl() -> value_type {
    auto best_it = this->c.rbegin();

    if (std::next(best_it) != this->c.rend() && std::next(best_it)->index() == best_it->index()) {
      auto best_hvc = m_hv_object.contribution(best_it->upper_bound());

      for (auto it = std::next(best_it); it != this->c.rend() && it->index() == best_it->index(); ++it) {
        auto hvc = m_hv_object.contribution(it->upper_bound());
        if (hvc > best_hvc) {
          best_it = it;
          best_hvc = hvc;
        }
      }
    }

    auto res = std::move(*best_it);
    *best_it = std::move(this->c.back());
    this->c.pop_back();

    return res;
  }

  hypervolume_type m_hv_object;
};

template <typename Node, typename Container = std::vector<std::pair<Node, double>>>
class epsilon_befs_queue : public mooutils::base_queue<epsilon_befs_queue<Node, Container>, Node, Container> {
 private:
  using objective_value_type = typename Node::solution_type::objective_vector_type::value_type;

  struct HeapCmp {
    auto operator()(auto const& a, auto const& b) const -> bool {
      return std::get<1>(a) > std::get<1>(b);
    }
  };

 public:
  using node_type = Node;

  template <typename... Args>
  explicit epsilon_befs_queue(Args&&... args, size_t no)
      : base_class_type(std::forward<Args>(args)...)
      , m_support(no, static_cast<objective_value_type>(0)) {}

  template <typename... Args>
  explicit epsilon_befs_queue(Args&&... args, std::vector<objective_value_type>&& eps_point)
      : base_class_type(std::forward<Args>(args)...)
      , m_support(std::move(eps_point)) {}

  [[nodiscard]] constexpr auto nodes() {
    return std::views::all(this->c) | std::views::transform([](auto& v) -> node_type& { return std::get<0>(v); });
  }

 private:
  using base_class_type = mooutils::base_queue<epsilon_befs_queue<Node, Container>, Node, Container>;
  using typename base_class_type::size_type;
  using typename base_class_type::value_type;

  friend base_class_type;

  // Push a new solution to the queue
  template <typename T>
  constexpr auto push_impl(T&& node) {
    for (auto const& s : node.lower_bound()) {
      for (std::size_t i = 0; i < m_support.size(); ++i) {
        if (s.objective_vector()[i] > m_support[i]) {
          m_support[i] = s.objective_vector()[i];
          m_support_changed = true;
        }
      }
    }

    auto eps = m_eps(node.upper_bound());
    if (m_support_changed) {
      this->c.emplace_back(std::forward<T>(node), eps);
      for (auto& [node, eps] : this->c) {
        eps = m_eps(node.upper_bound());
      }
      std::ranges::make_heap(this->c, HeapCmp());
      m_support_changed = false;
    } else {
      auto eps = m_eps(node.upper_bound());
      this->c.emplace_back(std::forward<T>(node), eps);
      std::ranges::push_heap(this->c, HeapCmp());
    }
  }

  constexpr auto pop_impl() -> value_type {
    std::ranges::pop_heap(this->c, HeapCmp());
    auto res = std::move(this->c.back().first);
    this->c.pop_back();
    return res;
  }

  [[nodiscard]] constexpr auto m_eps(auto const& v) const {
    double res = 0;
    for (std::size_t i = 0; i < m_support.size(); ++i) {
      if (m_support[i] == 0) {
        return std::numeric_limits<double>::max();
      } else {
        res = std::max(res, static_cast<double>(m_support[i]) / static_cast<double>(v[i]));
      }
    }
    return res;
  };

  std::vector<objective_value_type> m_support;
  bool m_support_changed = false;
};

template <typename Node, typename Container = std::vector<Node>>
class epsilon_bedfs_queue : public mooutils::base_queue<epsilon_bedfs_queue<Node, Container>, Node, Container> {
 private:
  using objective_value_type = typename Node::solution_type::objective_vector_type::value_type;

 public:
  using node_type = Node;

  template <typename... Args>
  explicit epsilon_bedfs_queue(Args&&... args, size_t no)
      : base_class_type(std::forward<Args>(args)...)
      , m_support(no, static_cast<objective_value_type>(0)) {}

  template <typename... Args>
  explicit epsilon_bedfs_queue(Args&&... args, std::vector<objective_value_type>&& eps_point)
      : base_class_type(std::forward<Args>(args)...)
      , m_support(std::move(eps_point)) {}

  [[nodiscard]] constexpr auto nodes() {
    return std::views::all(this->c);
  }

 private:
  using base_class_type = mooutils::base_queue<epsilon_bedfs_queue<Node, Container>, Node, Container>;
  using typename base_class_type::size_type;
  using typename base_class_type::value_type;

  friend base_class_type;

  // Push a new solution to the queue
  template <typename T>
  constexpr auto push_impl(T&& node) {
    for (auto const& s : node.lower_bound()) {
      for (std::size_t i = 0; i < m_support.size(); ++i) {
        if (s.objective_vector()[i] > m_support[i]) {
          m_support[i] = s.objective_vector()[i];
        }
      }
    }

    this->c.push_back(std::forward<T>(node));
  }

  constexpr auto pop_impl() -> value_type {
    auto best_it = this->c.rbegin();
    auto best_eps = m_eps(best_it->upper_bound());
    for (auto it = std::next(best_it); it != this->c.rend() && it->index() == best_it->index(); ++it) {
      auto eps = m_eps(it->upper_bound());
      if (eps < best_eps) {
        best_it = it;
        best_eps = eps;
      }
    }

    auto res = std::move(*best_it);
    *best_it = std::move(this->c.back());
    this->c.pop_back();

    return res;
  };

  [[nodiscard]] constexpr auto m_eps(auto const& v) const {
    double res = 0;
    for (std::size_t i = 0; i < m_support.size(); ++i) {
      if (m_support[i] == 0) {
        return std::numeric_limits<double>::max();
      } else {
        res = std::max(res, static_cast<double>(m_support[i]) / static_cast<double>(v[i]));
      }
    }
    return res;
  };

  std::vector<objective_value_type> m_support;
};

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

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration](auto& node) {
    for (auto& s : node.lower_bound()) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  auto empty =
      bbdetails::node<solution_type>(0, solution_type::empty(problem), objectives_orders_list, ratio_sum_order_list);
  update_lb_and_trace(empty);
  queue.push(std::move(empty));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    // if (iteration % 1000 == 0) {
    //   std::cerr << iteration << " " << queue.size() << " " << lb.size() << "\n";
    // }
    auto node = queue.pop();
    for (auto& branch : std::move(node).branches()) {
      if (!mooutils::strictly_dominates(lb, branch.upper_bound())) {
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

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration](auto& node) {
    for (auto& s : node.lower_bound()) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  queue.push(bbdetails::node(0, solution_type::empty(problem), objectives_orders_list, ratio_sum_order_list));

  for (iteration = 1; !queue.empty() && elapsed_sec() < timeout; ++iteration) {
    // if (iteration % 1000 == 0) {
    //   std::cerr << iteration << " " << queue.size() << " " << lb.size() << "\n";
    // }
    auto node = queue.pop();
    if (!mooutils::strictly_dominates(lb, node.upper_bound())) {
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

template <typename Solution, typename Problem, typename HVObject, typename AnytimeTrace>
[[nodiscard]] constexpr auto anytime_branch_and_bound(Problem const& problem, AnytimeTrace& anytime_trace,
                                                      HVObject&& hvobj, double timeout, size_t k, double delta,
                                                      size_t l) {
  using solution_type = Solution;

  auto elapsed_sec = [&anytime_trace]() {
    return std::chrono::duration<double>(anytime_trace.elapsed()).count();
  };

  auto no = problem.num_objectives();
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

  auto update_lb_and_trace = [&lb, &anytime_trace, &iteration](auto& node) {
    for (auto& s : node.lower_bound()) {
      auto it = lb.insert(s);
      if (it != lb.end()) {
        anytime_trace.add_solution(iteration, s);
      }
    }
  };

  auto empty =
      bbdetails::node<solution_type>(0, solution_type::empty(problem), objectives_orders_list, ratio_sum_order_list);
  update_lb_and_trace(empty);

  using node_type = bbdetails::node<solution_type>;
  using hv_object_type = std::remove_cvref_t<decltype(hvobj)>;
  using queue1_t = bbdetails::hypervolume_bedfs_queue<node_type, hv_object_type>;
  using queue2_t = bbdetails::hypervolume_befs_queue<node_type, hv_object_type>;
  using queue3_t = bbdetails::epsilon_bedfs_queue<bbdetails::node<solution_type>>;
  using queue4_t = bbdetails::epsilon_befs_queue<bbdetails::node<solution_type>>;
  std::variant<queue1_t, queue2_t, queue3_t, queue4_t> queue_v(std::in_place_index<0>, std::move(hvobj));

  auto queue_v_empty = [&queue_v]() -> bool {
    return std::visit([](auto&& arg) -> bool { return arg.empty(); }, queue_v);
  };

  auto queue_v_push = [&queue_v](node_type&& node) {
    std::visit([&node](auto&& arg) { return arg.push(std::move(node)); }, queue_v);
  };

  auto queue_v_pop = [&queue_v]() -> node_type {
    return std::visit([](auto&& arg) -> node_type { return arg.pop(); }, queue_v);
  };

  queue_v_push(std::move(empty));

  int phase = 0;
  std::vector<double> qualities(k, 0.0);
  for (iteration = 1; !queue_v_empty() && elapsed_sec() < timeout; ++iteration) {
    // if (iteration % 1000 == 0) {
    //   std::cerr << iteration << " " << queue.size() << " " << lb.size() << "\n";
    // }
    auto node = queue_v_pop();
    for (auto& branch : std::move(node).branches()) {
      if (!mooutils::strictly_dominates(lb, branch.upper_bound())) {
        update_lb_and_trace(branch);
        queue_v_push(std::move(branch));
      }
    }
    if (phase == 0) {
      double quality = static_cast<double>(anytime_trace.indicator_value());
      double relative_quality = (quality - qualities[iteration % k]) / quality;
      qualities[iteration % k] = quality;
      if (relative_quality < delta) {
        auto nodes = std::vector<node_type>{};
        for (auto& node : std::get<0>(queue_v).nodes()) {
          nodes.emplace_back(std::move(node));
        }
        if (no < 4) {
          auto hvobj = std::move(std::get<0>(queue_v).hv_object());
          queue_v.template emplace<1>(std::move(hvobj));
          for (auto& node : nodes) {
            std::get<1>(queue_v).push(std::move(node));
          }
        } else {
          using objective_value_type = typename solution_type::objective_vector_type::value_type;
          auto eps_point = std::vector<objective_value_type>(no, 0);
          for (auto const& s : lb) {
            for (size_t i = 0; i < no; ++i) {
              eps_point[i] = std::max(eps_point[i], s.objective_vector()[i]);
            }
          }
          queue_v.template emplace<3>(std::move(eps_point));
          for (auto& node : nodes) {
            std::get<3>(queue_v).push(std::move(node));
          }
        }
        phase = 1;
      }
    } else if (phase == 1) {
      if (std::visit([](auto&& arg) -> size_t { return arg.size(); }, queue_v) >= l) {
        auto nodes = std::vector<node_type>{};
        if (no < 4) {
          for (auto& node : std::get<1>(queue_v).nodes()) {
            nodes.emplace_back(std::move(node));
          }
        } else {
          for (auto& node : std::get<3>(queue_v).nodes()) {
            nodes.emplace_back(std::move(node));
          }
        }
        using objective_value_type = typename solution_type::objective_vector_type::value_type;
        auto eps_point = std::vector<objective_value_type>(no, 0);
        for (auto const& s : lb) {
          for (size_t i = 0; i < no; ++i) {
            eps_point[i] = std::max(eps_point[i], s.objective_vector()[i]);
          }
        }
        queue_v.template emplace<2>(std::move(eps_point));
        for (auto& node : nodes) {
          std::get<2>(queue_v).push(std::move(node));
        }
        phase = 2;
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
