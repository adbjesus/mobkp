#include <mobkp/anytime_trace.hpp>
#include <mobkp/pls.hpp>
#include <mobkp/problem.hpp>
#include <mobkp/solution.hpp>

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
  if (argc < 6) {
    fmt::print("Error: wrong number of arguments.\n");
    fmt::print("Usage: {} algorithm timeout seed inputfile outdir.\n", argv[0]);
    return EXIT_FAILURE;
  }

  auto algorithm = std::string(argv[1]);
  auto timeout = std::stod(argv[2]);
  auto seed = std::stoi(argv[3]);
  auto inputfile = argv[4];
  auto outdir = std::filesystem::path(argv[5]);

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

  auto anytime_trace =
      mobkp::anytime_trace(mooutils::incremental_hv<hv_data_type, ovec_type>(ovec_type(no, -1)));

  auto solutions = mooutils::unordered_set<solution_type>();

  if (algorithm == "pls") {
    auto initial_solution = solution_type::empty(problem);
    anytime_trace.add_solution(0, initial_solution);

    solutions.insert_unchecked(std::move(initial_solution));

    auto rng = std::mt19937_64(seed);  // std::random_device()());
    auto queue = mooutils::random_queue<solution_type, decltype(rng)>(std::move(rng));
    queue.push(*solutions.begin());

    mobkp::flip_exchange_pls(problem, solutions, queue, anytime_trace, timeout);
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
