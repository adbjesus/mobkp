#pragma once

#include <mooutils/indicators.hpp>

#include <chrono>
#include <memory>
#include <vector>

namespace mobkp {

template <typename IncrementalIndicator>
class anytime_trace {
 public:
  using clock_type = std::chrono::high_resolution_clock;
  using time_point_type = std::chrono::time_point<clock_type>;
  using duration_type = std::chrono::nanoseconds;
  using iteration_type = size_t;
  using incremental_indicator_type = IncrementalIndicator;
  using incremental_indicator_value_type = typename incremental_indicator_type::value_type;
  using measures_container_type =
      std::vector<std::tuple<duration_type, iteration_type, incremental_indicator_value_type>>;

  anytime_trace(incremental_indicator_type&& qi)
      : m_start(clock_type::now())
      , m_qi(std::move(qi)) {}

  template <typename Solution>
  auto add_solution(iteration_type const& iteration, Solution&& solution) {
    if (m_qi.insert(solution.objective_vector()) > 0) {
      m_measures.emplace_back(elapsed(), iteration, m_qi.value());
    }
  }

  auto elapsed() -> duration_type {
    return clock_type::now() - m_start;
  }

  auto begin() {
    return m_measures.begin();
  }

  auto end() {
    return m_measures.end();
  }

 private:
  time_point_type m_start;
  incremental_indicator_type m_qi;
  measures_container_type m_measures;
};

}  // namespace mobkp
