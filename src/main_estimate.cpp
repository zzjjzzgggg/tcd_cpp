/**
 * Copyright (C) by J.Z. (06/30/2017 23:27)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include "ITS_sampler.h"
#include "ITSC_sampler.h"
#include "SGS_sampler.h"

#include "em_alg.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph full path");
DEFINE_string(theta, "", "theta groundtruth file name");
DEFINE_string(output, "", "output file name");

DEFINE_int32(mx_tc, 2000, "maximum triadic cardinality");
DEFINE_int32(mx_i, 128, "maximum triadic cardinality");
DEFINE_int32(mx_iter_theta, 1000, "maximum iterations for estimating theta");
DEFINE_int32(mx_iter_alpha, 10, "maximum iterations for estimating alpha");
DEFINE_int32(trials, 10, "trials per core");
DEFINE_int32(cores, std::thread::hardware_concurrency(), "cores");

DEFINE_double(p_node, 0.1, "node sampling rate");
DEFINE_double(p_edge, 0.1, "edge sampling rate");
DEFINE_double(eps_theta, 1e-3, "threshold of estimating theta");
DEFINE_double(eps_alpha, 1e-3, "threshold of estimating alpha");

std::mutex print_mutex;
void echo(const int n_suc, vector<int>& states) {
    std::lock_guard<std::mutex> guard(print_mutex);
    for (int core = 0; core < FLAGS_cores; core++)
        printf("  [%d] %2d", core, states[core]);
    printf("  suc: %3d\r", n_suc);
    std::fflush(stdout);
}

std::mutex avg_est_mutex;
int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("usage:");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    Sampler::Config confs;
    confs.graph_fnm = FLAGS_graph;
    confs.mx_tc = FLAGS_mx_tc;
    confs.p_edge = FLAGS_p_edge;
    confs.p_node = FLAGS_p_node;
    confs.echo();

#ifdef S_ITS
    ITSSampler sampler(&confs);
#elif S_ITSC
    ITSCSampler sampler(&confs);
#elif S_SGS
    SGSSampler sampler(&confs);
#else
    printf("Warning: Sampler is not specified! Use ISSampler by default!\n");
    ITSSampler sampler(&confs);
#endif

    // sampler.check();
    sampler.info();

    EM::Config conf;
    conf.mx_tc = FLAGS_mx_tc;
    conf.mx_i = FLAGS_mx_i;
    conf.mx_iter_theta = FLAGS_mx_iter_theta;
    conf.mx_iter_alpha = FLAGS_mx_iter_alpha;
    conf.eps_theta = FLAGS_eps_theta;
    conf.eps_alpha = FLAGS_eps_alpha;
    conf.echo();

    printf("trials per processor: %d\n\n", FLAGS_trials);

    auto truth = ioutils::loadPrVec<int, double>(
        strutils::subFilename(FLAGS_graph, FLAGS_theta));

    int n_suc = 0;
    double alpha = 0;

    vector<double> theta_hat(truth.size(), 0), err(truth.size(), 0);
    vector<int> states(FLAGS_cores, 0);

    // define the task
    auto fun_task = [&](int core) {
        EM em{&conf, &sampler};
        for (int trial = 0; trial < FLAGS_trials; trial++) {
            em.init();
            if (em.exec()) {
                auto rst = em.get();
                std::lock_guard<std::mutex> guard(avg_est_mutex);
                n_suc++;
                alpha += rst[0];
                int pos = 0;
                for (auto&& pr : truth) {
                    theta_hat[pos] += rst[pr.first + 1];
                    err[pos] += std::pow(rst[pr.first + 1] - pr.second, 2);
                    pos++;
                }
            }
            states[core]++;
            echo(n_suc, states);
        }
    };

    echo(n_suc, states);
    vector<std::future<void>> futures;
    futures.reserve(FLAGS_cores);
    for (int core = 0; core < FLAGS_cores; core++)
        futures.push_back(std::async(std::launch::async, fun_task, core));
    for (auto& f : futures) f.get();

    printf("\n\n%d out of %d succed.\n", n_suc, FLAGS_trials * FLAGS_cores);

    // saving ...
    if (n_suc > 0) {
        vector<std::tuple<int, double, double>> rst;
        rst.reserve(truth.size() + 1);

        alpha /= n_suc;
        rst.emplace_back(-10, alpha, 0);
        printf("alpha = %.6e\n", alpha);

        int pos = 0;
        for (auto&& pr : truth) {
            rst.emplace_back(pr.first, theta_hat[pos] / n_suc,
                             std::sqrt(err[pos] / n_suc) / pr.second);
            pos++;
        }
        ioutils::saveTupleVec(rst,
                              strutils::subFilename(FLAGS_graph, FLAGS_output),
                              true, "{}\t{:.6e}\t{:.6e}\n");
    }

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
