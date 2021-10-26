#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: missing arguments. Usage: validate_solutions.R DIR", call. = FALSE)
}

dir <- args[1]
problem_file <- sprintf("%s/problem.dat", dir)
solutions_file <- sprintf("%s/solutions.dat", dir)

aux <- as.integer(readLines(problem_file, n = 4))
ni <- aux[1]
no <- aux[2]
nc <- aux[3]
capacities <- aux[seq(4,3+nc)]

dat <- read.table(problem_file, skip=4)

sol <- read.table(solutions_file)

for (i in seq(nrow(sol))) {
  aux <- as.matrix(sol[i,seq(no+nc+1,no+nc+ni)]) %*% as.matrix(dat)
  if (!all(aux == sol[i,seq(1,no+nc)])) {
    tmp <- as.matrix(sol[i,seq(1,no+nc)])
    cat(sprintf("Error on solution %d.", i))
    cat(aux)
    cat("  !=  ")
    cat(tmp)
    cat("\n")
  } else if (!all(aux[seq(no+1,no+nc)] <= capacities)) {
    tmp <- aux[no+1:no+nc]
    cat(sprintf("Error on solution %d", i))
    cat(tmp)
    cat("  !<=  ")
    cat(capacities)
    cat("\n")
  }
}
