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
#include <string>

int main(int argc, char** argv) {
  using rng_type = std::mt19937_64;
  using rng_result_type = typename rng_type::result_type;

  CLI::App app{"Multi-Objective Binary Knapsack Solver"};

  std::string algorithm;
  app.add_option("-a,--algorithm", algorithm, "Algorithm to use")->required();

  rng_result_type seed;
  app.add_option("-s,--seed", seed, "Seed for random generator")->required();

  double timeout;
  app.add_option("-t,--timeout", timeout, "Timeout")->required();

  std::filesystem::path inputfile;
  app.add_option("-i,--input-file", inputfile, "Input file")->required();

  std::filesystem::path outdir;
  app.add_option("-o,--output-directory", outdir, "Output directory")->required();

  CLI11_PARSE(app, argc, argv);

  // if (argc < 6) {
  //   fmt::print("Error: wrong number of arguments.\n");
  //   fmt::print("Usage: {} algorithm timeout seed inputfile outdir.\n", argv[0]);
  //   return EXIT_FAILURE;
  // }

  // auto algorithm = std::string(argv[1]);
  // auto timeout = std::stod(argv[2]);
  // auto seed = std::stoull(argv[3]);
  // auto inputfile = argv[4];
  // auto outdir = std::filesystem::path(argv[5]);

  using data_type = int_fast32_t;
  using dvec_type = std::vector<bool>;
  using ovec_type = std::vector<data_type>;
  using cvec_type = std::vector<data_type>;

  using hv_data_type = boost::multiprecision::int256_t;

  using problem_type = mobkp::problem<data_type>;
  using solution_type = mobkp::solution<problem_type, dvec_type, ovec_type, cvec_type>;

  auto const problem = problem_type::from_stream(std::ifstream(inputfile));
  // auto const ni = problem.num_items();
  auto const no = problem.num_objectives();
  // auto const nc = problem.num_constraints();

  auto hvref = ovec_type(no, -1);
  auto anytime_trace = mobkp::anytime_trace(mooutils::incremental_hv<hv_data_type, ovec_type>(hvref));

  auto solutions = mooutils::unordered_set<solution_type>();

  if (algorithm == "pls") {
    auto initial_solution = solution_type::empty(problem);
    anytime_trace.add_solution(0, initial_solution);

    solutions.insert_unchecked(std::move(initial_solution));

    auto rng = std::mt19937_64(seed);  // std::random_device()());
    auto queue = mooutils::random_queue<solution_type, decltype(rng)>(std::move(rng));
    queue.push(*solutions.begin());

    mobkp::flip_exchange_pls(problem, solutions, queue, anytime_trace, timeout);
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
    // TODO Make this a parameter
    size_t l = 100;
    solutions = mobkp::anytime_eps<solution_type>(problem, l, hvref, anytime_trace, timeout);
  } else if (algorithm == "eager-branch-and-bound-dfs") {
    // auto queue = mooutils::lifo_queue<solution_type>();
    // auto sols = mobkp::eager_branch_and_bound<solution_type>(problem, queue, anytime_trace, timeout);
    // for (auto&& s : sols) {
    //   solutions.insert_unchecked(std::move(s));
    // }
  } else if (algorithm == "eager-branch-and-bound-bfs") {
  } else if (algorithm == "lazy-branch-and-bound") {
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
    fmt::print(solution_stream, "{:d}\n", fmt::join(s.decision_vector(), " "));
  }
  solution_stream.close();

  // Print problem.dat useful to validate solutions
  auto problem_stream = std::ofstream(outdir / "problem.dat");
  fmt::print(problem_stream, "{}\n", problem.num_items());
  fmt::print(problem_stream, "{}\n", problem.num_objectives());
  fmt::print(problem_stream, "{}\n", problem.num_constraints());
  fmt::print(problem_stream, "{}\n", fmt::join(problem.weight_capacities(), " "));
  for (size_t i = 0; i < problem.num_items(); ++i) {
    fmt::print(problem_stream, "{} ", fmt::join(problem.item_values(i), " "));
    fmt::print(problem_stream, "{}\n", fmt::join(problem.item_weights(i), " "));
  }
  problem_stream.close();

  return EXIT_SUCCESS;
}
