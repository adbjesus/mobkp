#ifndef MOBKP_SCALARIZATION_HPP_
#define MOBKP_SCALARIZATION_HPP_

#include "minknap.hpp"

#include <apm/apm.hpp>
#include <boost/icl/interval_set.hpp>
#include <glpk.h>
#include <mooutils/sets.hpp>

#include <array>
#include <optional>
#include <queue>
#include <vector>

namespace mobkp {

template <typename Problem>
struct epsilon_constraint {
  glp_prob* m_lp;
  glp_iocp* m_iocp;

  std::vector<int> m_ia;
  std::vector<int> m_ja;
  std::vector<double> m_ca;

  std::reference_wrapper<const Problem> m_problem;

  epsilon_constraint(Problem const& problem, int o)
      : m_problem(problem) {
    int const ni = static_cast<int>(m_problem.get().num_items());
    int const no = static_cast<int>(m_problem.get().num_objectives());
    int const nc = static_cast<int>(m_problem.get().num_constraints());

    assert(o >= 0 && o < no);

    m_lp = glp_create_prob();

    glp_set_obj_dir(m_lp, GLP_MAX);
    glp_add_rows(m_lp, no - 1 + nc);
    glp_add_cols(m_lp, ni);

    for (int i = 0; i < ni; ++i) {
      glp_set_col_kind(m_lp, i + 1, GLP_BV);
      glp_set_obj_coef(m_lp, i + 1, double(m_problem.get().item_value(i, o)));
    }

    m_ia.push_back(0);
    m_ja.push_back(0);
    m_ca.push_back(0.0);

    int k = 1;
    for (int i = 0; i < no; ++i) {
      if (i == o)
        continue;
      for (int j = 0; j < ni; ++j) {
        m_ia.push_back(k);
        m_ja.push_back(j + 1);
        m_ca.push_back(static_cast<double>(m_problem.get().item_value(j, i)));
      }
      ++k;
    }

    for (int i = 0; i < nc; ++i) {
      auto bnd = static_cast<double>(m_problem.get().weight_capacity(i));
      glp_set_row_bnds(m_lp, k, GLP_UP, 0.0, bnd);
      for (int j = 0; j < ni; ++j) {
        m_ia.push_back(k);
        m_ja.push_back(j + 1);
        m_ca.push_back(static_cast<double>(m_problem.get().item_weight(j, i)));
      }
      ++k;
    }

    glp_load_matrix(m_lp, ni * (no - 1 + nc), m_ia.data(), m_ja.data(), m_ca.data());

    m_iocp = new glp_iocp;
    glp_init_iocp(m_iocp);
    m_iocp->presolve = GLP_ON;
    m_iocp->msg_lev = GLP_MSG_ERR;
    m_iocp->tol_int = 1e-9;
    m_iocp->tol_obj = 1e-9;
  }

  ~epsilon_constraint() {
    glp_free(m_lp);
    m_lp = NULL;

    delete m_iocp;
    m_iocp = NULL;
  }

  template <typename Solution, typename R>
  auto solve(R const& obounds) -> std::optional<Solution> {
    assert(obounds.size() == m_problem.get().num_objectives() - 1);

    for (int i = 0; i < static_cast<int>(obounds.size()); ++i) {
      auto bnd = static_cast<double>(obounds[i]);
      glp_set_row_bnds(m_lp, i + 1, GLP_LO, bnd, 0.0);
    }

    if (glp_intopt(m_lp, m_iocp) != 0 || glp_mip_status(m_lp) != GLP_OPT) {
      return {};
    }

    auto s = Solution::empty(m_problem.get());

    int ni = static_cast<int>(m_problem.get().num_items());
    for (int i = 0; i < ni; ++i) {
      if (glp_mip_col_val(m_lp, i + 1) >= 0.5) {
        s.flip_to_one_unchecked(i);
      }
    }

    return s;
  }
};

template <typename Problem>
struct weighted_sum {
  glp_prob* m_lp;
  glp_iocp* m_iocp;

  std::vector<int> m_ia;
  std::vector<int> m_ja;
  std::vector<double> m_ca;

  std::reference_wrapper<const Problem> m_problem;

  explicit weighted_sum(Problem const& problem)
      : m_problem(problem) {
    int const ni = static_cast<int>(m_problem.get().num_items());
    int const nc = static_cast<int>(m_problem.get().num_constraints());

    m_lp = glp_create_prob();

    glp_set_obj_dir(m_lp, GLP_MAX);
    glp_add_rows(m_lp, nc);
    glp_add_cols(m_lp, ni);

    for (int i = 0; i < ni; ++i) {
      glp_set_col_kind(m_lp, i + 1, GLP_BV);
    }

    m_ia.push_back(0);
    m_ja.push_back(0);
    m_ca.push_back(0.0);

    for (int i = 0; i < nc; ++i) {
      auto bnd = static_cast<double>(m_problem.get().weight_capacity(i));
      glp_set_row_bnds(m_lp, i + 1, GLP_UP, 0.0, bnd);
      for (int j = 0; j < ni; ++j) {
        m_ia.push_back(i + 1);
        m_ja.push_back(j + 1);
        m_ca.push_back(static_cast<double>(m_problem.get().item_weight(j, i)));
      }
    }

    glp_load_matrix(m_lp, ni * nc, m_ia.data(), m_ja.data(), m_ca.data());

    m_iocp = new glp_iocp;
    glp_init_iocp(m_iocp);
    m_iocp->msg_lev = GLP_MSG_ERR;
    m_iocp->presolve = GLP_ON;
    m_iocp->tol_int = 1e-9;
    m_iocp->tol_obj = 1e-9;
  }

  ~weighted_sum() {
    glp_free(m_lp);
    m_lp = NULL;

    delete m_iocp;
    m_iocp = NULL;
  }

  template <typename Solution, typename Weights>
  auto solve(Weights const& weights) -> std::optional<Solution> {
    assert(weights.size() == m_problem.get().num_objectives());

    int const ni = static_cast<int>(m_problem.get().num_items());
    int const no = static_cast<int>(m_problem.get().num_objectives());

    for (int i = 0; i < ni; ++i) {
      auto const values = m_problem.get().item_values(i);
      double coef = 0.0;
      for (int j = 0; j < no; ++j)
        coef += static_cast<double>(weights[j] * values[j]);
      glp_set_obj_coef(m_lp, i + 1, coef);
    }

    if (glp_intopt(m_lp, m_iocp) != 0 || glp_mip_status(m_lp) != GLP_OPT) {
      return {};
    }

    auto s = Solution::empty(m_problem.get());

    for (int i = 0; i < ni; ++i) {
      if (glp_mip_col_val(m_lp, i + 1) >= 0.5) {
        s.flip_to_one_unchecked(i);
      }
    }

    return s;
  }
};

template <typename Problem>
struct minknap_solver {
  std::reference_wrapper<const Problem> m_problem;
  std::vector<std::int64_t> w;

  explicit minknap_solver(Problem const& problem)
      : m_problem(problem) {
    int const ni = static_cast<int>(m_problem.get().num_items());
    w.reserve(ni);
    for (int i = 0; i < ni; ++i) {
      w.push_back(m_problem.get().item_weight(i, 0));
    }
  }

  template <typename Solution, typename Weights>
  auto solve(Weights const& weights) -> std::optional<Solution> {
    assert(weights.size() == m_problem.get().num_objectives());

    int const ni = static_cast<int>(m_problem.get().num_items());
    int const no = static_cast<int>(m_problem.get().num_objectives());

    std::vector<int> x(ni, 0);
    std::vector<std::int64_t> p;
    p.reserve(ni);

    for (int i = 0; i < ni; ++i) {
      auto const values = m_problem.get().item_values(i);
      std::int64_t coef = 0.0;
      for (int j = 0; j < no; ++j)
        coef += static_cast<std::int64_t>(weights[j] * values[j]);
      p.push_back(coef);
    }

    minknap::minknap(ni, p.data(), w.data(), x.data(), m_problem.get().weight_capacity(0));

    auto s = Solution::empty(m_problem.get());
    for (int i = 0; i < ni; ++i) {
      if (x[i] != 0) {
        s.flip_to_one_unchecked(i);
      }
    }

    return s;
  }
};

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] auto dws(Problem const& problem, AnytimeTrace& anytime_trace, double timeout) {
  using solution_type = Solution;

  if (problem.num_objectives() > 2) {
    throw("DWS algorithm is only implemented for 2 objectives");
  }

  auto sols = mooutils::unordered_set<solution_type>();

  // auto solver = weighted_sum(problem);
  auto solver = minknap_solver(problem);
  auto solver0 = epsilon_constraint(problem, 0);
  auto solver1 = epsilon_constraint(problem, 1);

  // lexmax0
  auto aux0 = solver.template solve<solution_type>(std::array{1, 0});
  if (!aux0.has_value())
    return sols;
  // aux0 = solver1.template solve<solution_type>(std::array{aux0.value().objective_vector()[0]});
  auto ov0 = std::array{aux0.value().objective_vector()[0], aux0.value().objective_vector()[1]};
  anytime_trace.add_solution(1, aux0.value());
  sols.insert_unchecked(std::move(aux0.value()));

  // lexmax1
  auto aux1 = solver.template solve<solution_type>(std::array{0, 1});
  if (!aux1.has_value())
    return sols;
  // aux1 = solver0.template solve<solution_type>(std::array{aux1.value().objective_vector()[1]});
  auto ov1 = std::array{aux1.value().objective_vector()[0], aux1.value().objective_vector()[1]};
  anytime_trace.add_solution(2, aux1.value());
  sols.insert_unchecked(std::move(aux1.value()));

  if (sols.size() < 2)
    return sols;

  std::queue<std::array<decltype(ov1), 2>> q;
  q.push(std::array{ov0, ov1});

  for (size_t i = 3; !q.empty() && anytime_trace.elapsed_sec() < timeout; ++i) {
    auto v = std::move(q.front());
    q.pop();

    auto l1 = v[1][1] - v[0][1];
    auto l2 = v[0][0] - v[1][0];

    auto aux = solver.template solve<solution_type>(std::array{l1, l2});
    if (!aux.has_value())
      continue;

    auto ov = std::array{aux.value().objective_vector()[0], aux.value().objective_vector()[1]};
    if (ov == v[0] || ov == v[1])
      continue;

    anytime_trace.add_solution(i, aux.value());
    sols.insert_unchecked(std::move(aux.value()));

    q.push(std::array{v[0], ov});
    q.push(std::array{ov, v[1]});
  }

  return sols;
}

template <typename Solution, typename Problem, typename AnytimeTrace>
[[nodiscard]] auto aeps(Problem const& problem, AnytimeTrace& anytime_trace, double timeout) {
  using solution_type = Solution;

  if (problem.num_objectives() > 2) {
    throw("AEPS algorithm is only implemented for 2 objectives");
  }

  auto sols = mooutils::unordered_set<solution_type>();

  auto solver0 = epsilon_constraint(problem, 0);
  auto solver1 = epsilon_constraint(problem, 1);

  auto aux = solver0.template solve<solution_type>(std::array{0});

  for (size_t i = 1; aux.has_value() && anytime_trace.elapsed_sec() < timeout; ++i) {
    aux = solver1.template solve<solution_type>(std::array{aux.value().objective_vector()[0]});
    anytime_trace.add_solution(i, aux.value());
    sols.insert_unchecked(aux.value());
    aux = solver0.template solve<solution_type>(std::array{aux.value().objective_vector()[1] + 1});
  }

  return sols;
}

template <typename Solution, typename Problem, typename HvRef, typename AnytimeTrace>
[[nodiscard]] auto anytime_eps(Problem const& problem, size_t l, HvRef const& hvref, AnytimeTrace& anytime_trace,
                               double timeout) {
  using solution_type = Solution;
  using value_type = solution_type::objective_vector_type::value_type;

  if (problem.num_objectives() > 2) {
    throw("Anytime-EPS algorithm is only implemented for 2 objectives");
  }

  auto sols = mooutils::unordered_set<solution_type>();

  // auto solver = weighted_sum(problem);
  auto solver = minknap_solver(problem);
  auto solver0 = epsilon_constraint(problem, 0);
  auto solver1 = epsilon_constraint(problem, 1);

  // lexmax0
  auto aux0 = solver.template solve<solution_type>(std::array{1, 0});
  if (!aux0.has_value())
    return sols;
  auto ov0 = std::array{aux0.value().objective_vector()[0], aux0.value().objective_vector()[1]};
  anytime_trace.add_solution(1, aux0.value());
  sols.insert_unchecked(std::move(aux0.value()));

  // lexmax1
  auto aux1 = solver.template solve<solution_type>(std::array{0, 1});
  if (!aux1.has_value())
    return sols;
  auto ov1 = std::array{aux1.value().objective_vector()[0], aux1.value().objective_vector()[1]};
  anytime_trace.add_solution(2, aux1.value());
  sols.insert_unchecked(std::move(aux1.value()));

  if (sols.size() < 2)
    return sols;

  // Find middle (according to weighted sum) point
  auto l1 = ov1[1] - ov0[1];
  auto l2 = ov0[0] - ov1[0];

  auto aux = solver.template solve<solution_type>(std::array{l1, l2});
  if (!aux.has_value())
    return sols;

  auto ov = std::array{aux.value().objective_vector()[0], aux.value().objective_vector()[1]};
  if (ov == ov0 || ov == ov1)
    return sols;

  anytime_trace.add_solution(3, aux.value());
  sols.insert_unchecked(std::move(aux.value()));

  // Compute p, d, and reference point r (assuming scalarization to [0-1]^2)
  auto denom0 = ov0[0] - ov1[0];
  auto denom1 = ov1[1] - ov0[1];
  auto p0 = static_cast<double>(ov[0] - ov1[0]) / static_cast<double>(denom0);
  auto p1 = static_cast<double>(ov[1] - ov0[1]) / static_cast<double>(denom1);
  auto d = std::log(0.5) / std::log((p0 + p1) / 2.0);
  auto r0 = static_cast<double>(hvref[0] - ov1[0]) / static_cast<double>(denom0);
  auto r1 = static_cast<double>(hvref[1] - ov0[1]) / static_cast<double>(denom1);

  // Get segments
  auto segments = apm::piecewise_segments(l, d);

  // Get model
  auto model = apm::greedy_model(segments, {r0, r1});
  auto iter = model.begin();
  auto last = model.end();

  // Interval set
  auto interval_set = boost::icl::interval_set<value_type>();

  for (size_t i = 4; iter != last && anytime_trace.elapsed_sec() < timeout; ++iter, ++i) {
    auto eps = static_cast<value_type>((*iter).first.y * static_cast<double>(denom1)) + ov0[1] + 1;
    if (boost::icl::contains(interval_set, eps))
      continue;
    aux = solver0.template solve<solution_type>(std::array{eps});
    if (!aux.has_value())
      continue;
    aux = solver1.template solve<solution_type>(std::array{aux->objective_vector()[0]});
    if (!aux.has_value())
      continue;
    interval_set.add(boost::icl::interval<value_type>::closed(eps, aux->objective_vector()[1]));
    anytime_trace.add_solution(i, *aux);
    sols.insert(std::move(*aux));
  }

  return sols;
}

}  // namespace mobkp

#endif
