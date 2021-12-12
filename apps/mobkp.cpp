#include <mobkp/anytime_trace.hpp>
#include <mobkp/bb.hpp>
#include <mobkp/dp.hpp>
#include <mobkp/pls.hpp>
#include <mobkp/problem.hpp>
#include <mobkp/scalarization.hpp>
#include <mobkp/solution.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <mooutils/indicators.hpp>
#include <mooutils/queues.hpp>
#include <mooutils/sets.hpp>

#include <bitset>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

enum class queue_type {
  FIFO,
  LIFO,
  Random,
};

enum class problem_order {
  Default,
  Random,
  RankMin,
  RankMax,
  RankSum,
};

int main(int argc, char** argv) {
  CLI::App app{"Multi-Objective Binary Knapsack Solver"};

  std::string algorithm;
  app.add_option("-a,--algorithm", algorithm, "Algorithm to use")->required();

  int seed;
  app.add_option("-s,--seed", seed, "Seed for random generator")->required();

  double timeout;
  app.add_option("-t,--timeout", timeout, "Timeout")->required();

  std::filesystem::path inputfile;
  app.add_option("-i,--input-file", inputfile, "Input file")->required();

  std::filesystem::path outdir;
  app.add_option("-o,--output-directory", outdir, "Output directory")->required();

  int anytime_eps_l{100};
  app.add_option("--anytime-eps-l", anytime_eps_l, "Parameter 'l' for anytime-eps algorithm", true);

  queue_type queue_t{queue_type::FIFO};
  app.add_option("--queue-type", queue_t, "Queue type for algorithms that require a queue (e.g. pls and b&b)", true)
      ->transform(CLI::CheckedTransformer(
          std::map<std::string, queue_type>{
              {"fifo", queue_type::FIFO}, {"lifo", queue_type::LIFO}, {"random", queue_type::Random}},
          CLI::ignore_case));

  problem_order order{problem_order::Default};
  app.add_option("--order", order, "Order for problem items", true)
      ->transform(CLI::CheckedTransformer(std::map<std::string, problem_order>{{"default", problem_order::Default},
                                                                               {"random", problem_order::Random},
                                                                               {"rank_min", problem_order::RankMin},
                                                                               {"rank_max", problem_order::RankMax},
                                                                               {"rank_sum", problem_order::RankSum}},
                                          CLI::ignore_case));

  CLI11_PARSE(app, argc, argv);

  using data_type = int_fast32_t;
  using dvec_type = std::vector<bool>;
  using ovec_type = std::vector<data_type>;
  using cvec_type = std::vector<data_type>;

  using hv_data_type = boost::multiprecision::int256_t;

  using problem_type = mobkp::ordered_problem<mobkp::problem<data_type>>;
  using solution_type = mobkp::solution<problem_type, dvec_type, ovec_type, cvec_type>;

  auto const orig_problem = mobkp::problem<data_type>::from_stream(std::ifstream(inputfile));
  auto const ni = orig_problem.num_items();
  auto const no = orig_problem.num_objectives();

  auto seed_seq = std::seed_seq{seed};
  auto hvref = ovec_type(no, -1);

  auto anytime_trace = mobkp::anytime_trace(mooutils::incremental_hv<hv_data_type, ovec_type>(hvref));

  std::vector<size_t> index_order;
  switch (order) {
    case problem_order::Default:
      index_order.reserve(ni);
      for (size_t i = 0; i < ni; ++i) {
        index_order.emplace_back(i);
      }
      break;
    case problem_order::Random:
      index_order.reserve(ni);
      for (size_t i = 0; i < ni; ++i) {
        index_order.emplace_back(i);
      }
      {
        auto rng = std::mt19937_64(seed_seq);
        std::shuffle(index_order.begin(), index_order.end(), rng);
      }
      break;
    case problem_order::RankMin:
      index_order = mobkp::rank_min_order(orig_problem, mobkp::objectives_orders(orig_problem));
      break;
    case problem_order::RankMax:
      index_order = mobkp::rank_max_order(orig_problem, mobkp::objectives_orders(orig_problem));
      break;
    case problem_order::RankSum:
      index_order = mobkp::rank_sum_order(orig_problem, mobkp::objectives_orders(orig_problem));
      break;
  }

  auto const problem = problem_type(orig_problem, index_order);

  auto solutions = mooutils::unordered_set<solution_type>();

  if (algorithm == "pls") {
    auto initial_solution = solution_type::empty(problem);
    anytime_trace.add_solution(0, initial_solution);
    solutions.insert_unchecked(std::move(initial_solution));
    switch (queue_t) {
      case queue_type::FIFO: {
        auto queue = mooutils::fifo_queue<solution_type>();
        queue.push(*solutions.begin());
        mobkp::flip_exchange_pls(problem, solutions, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::LIFO: {
        auto queue = mooutils::lifo_queue<solution_type>();
        queue.push(*solutions.begin());
        mobkp::flip_exchange_pls(problem, solutions, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::Random: {
        auto rng = std::mt19937_64(seed_seq);
        auto queue = mooutils::random_queue<solution_type, decltype(rng)>(std::move(rng));
        queue.push(*solutions.begin());
        mobkp::flip_exchange_pls(problem, solutions, queue, anytime_trace, timeout);
        break;
      }
    }
  } else if (algorithm == "nemull-dp") {
    solutions = mobkp::nemull_dp<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "bhv-dp") {
    solutions = mobkp::bhv_dp<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "fpsv-dp") {
    solutions = mobkp::fpsv_dp<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "anytime-dp") {
    solutions = mobkp::anytime_dp<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "dws") {
    solutions = mobkp::dws<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "aeps") {
    solutions = mobkp::aeps<solution_type>(problem, anytime_trace, timeout);
  } else if (algorithm == "anytime-eps") {
    solutions = mobkp::anytime_eps<solution_type>(problem, anytime_eps_l, hvref, anytime_trace, timeout);
  } else if (algorithm == "eager-bb") {
    switch (queue_t) {
      case queue_type::FIFO: {
        auto queue = mooutils::fifo_queue<mobkp::bbdetails::node<solution_type>>();
        solutions = mobkp::eager_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::LIFO: {
        auto queue = mooutils::lifo_queue<mobkp::bbdetails::node<solution_type>>();
        solutions = mobkp::eager_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::Random: {
        auto rng = std::mt19937_64(seed_seq);
        auto queue = mooutils::random_queue<mobkp::bbdetails::node<solution_type>, decltype(rng)>(std::move(rng));
        solutions = mobkp::eager_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
    }
  } else if (algorithm == "lazy-bb") {
    switch (queue_t) {
      case queue_type::FIFO: {
        auto queue = mooutils::fifo_queue<mobkp::bbdetails::node<solution_type>>();
        solutions = mobkp::lazy_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::LIFO: {
        auto queue = mooutils::lifo_queue<mobkp::bbdetails::node<solution_type>>();
        solutions = mobkp::lazy_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
      case queue_type::Random: {
        auto rng = std::mt19937_64(seed_seq);
        auto queue = mooutils::random_queue<mobkp::bbdetails::node<solution_type>, decltype(rng)>(std::move(rng));
        solutions = mobkp::lazy_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
        break;
      }
    }
  } else {
    fmt::print("Error: unknown algorithm.\n");
    return EXIT_FAILURE;
  }

  if (!std::filesystem::exists(outdir)) {
    std::filesystem::create_directory(outdir);
  }

  auto anytime_stream = std::ofstream(outdir / "anytime.dat");
  fmt::print(anytime_stream, "time,iteration,hv\n");
  for (auto const& m : anytime_trace) {
    fmt::print(anytime_stream, "{},{},{}\n",  // noformat
               std::get<0>(m).count(),        // noformat
               std::get<1>(m),                // noformat
               std::get<2>(m));
  }
  anytime_stream.close();

  // Before printing the solutions we sort them by objective vector/constraint
  // vector/decision_vector
  std::sort(solutions.begin(), solutions.end(), [](auto const& lhs, auto const& rhs) {
    if (lhs.objective_vector() < rhs.objective_vector()) {
      return true;
    }
    if (lhs.objective_vector() > rhs.objective_vector()) {
      return false;
    }
    if (lhs.constraint_vector() < rhs.constraint_vector()) {
      return true;
    }
    if (lhs.constraint_vector() > rhs.constraint_vector()) {
      return false;
    }
    if (lhs.decision_vector() < rhs.decision_vector()) {
      return true;
    }
    if (lhs.decision_vector() > rhs.decision_vector()) {
      return false;
    }
    return false;
  });

  auto solution_stream = std::ofstream(outdir / "solutions.dat");
  for (auto const& s : solutions) {
    fmt::print(solution_stream, "{} ", fmt::join(s.objective_vector(), " "));
    fmt::print(solution_stream, "{} ", fmt::join(s.constraint_vector(), " "));
    std::vector<bool> dvec(ni, false);
    for (size_t i = 0; i < ni; ++i) {
      dvec[index_order[i]] = s.decision_vector()[i];
    }
    fmt::print(solution_stream, "{:d}\n", fmt::join(dvec, " "));
  }
  solution_stream.close();

  // Print problem.dat useful to validate solutions
  auto problem_stream = std::ofstream(outdir / "problem.dat");
  fmt::print(problem_stream, "{}\n", orig_problem.num_items());
  fmt::print(problem_stream, "{}\n", orig_problem.num_objectives());
  fmt::print(problem_stream, "{}\n", orig_problem.num_constraints());
  fmt::print(problem_stream, "{}\n", fmt::join(orig_problem.weight_capacities(), " "));
  for (size_t i = 0; i < orig_problem.num_items(); ++i) {
    fmt::print(problem_stream, "{} ", fmt::join(orig_problem.item_values(i), " "));
    fmt::print(problem_stream, "{}\n", fmt::join(orig_problem.item_weights(i), " "));
  }
  problem_stream.close();

  return EXIT_SUCCESS;
}
