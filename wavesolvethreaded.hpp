#pragma once
#include <cstring>
#include <charconv>
#include <vector>
#include <array>
#include <ranges>
#include <thread>
#include <semaphore>
#include <atomic>
#include <barrier>
#include "WaveOrthotope.hpp"



namespace {
    // Create a vector of threads, each of which will run F(thread_id) in a while loop until F() returns false
    auto looping_threadpool(auto thread_count, auto F) {
        std::vector<std::jthread> threads;
        threads.reserve(thread_count);
        for (size_t tid=0; tid<thread_count; tid++) {
            threads.emplace_back([F, tid]{
                while (F(tid));
            });
        }
        return threads;
    }
};



class wavesolvethreaded: public WaveOrthotope {
    // Threading-related members
    bool continue_iteration; // used to tell the looping threadpool to terminate at the end of the simulation
    const size_type nthreads;
    std::barrier<> ds_barrier, step_barrier;
    std::vector<std::jthread> ds_workers, step_workers;
    std::atomic<value_type> ds_aggregator; // used to reduce dsteepness from each thread
    value_type iter_dt; // Used to distribute dt to each thread



    // Determine which rows a certain thread is in charge of
    auto this_thread_cell_range(auto tid) {
        return mr::split_range(rows, tid, nthreads);
    }




    // Help message to show that SOLVER_NUM_THREADS controls thread counts
    inline static const std::string help_message =
            "Set the environment variable SOLVER_NUM_THREADS to a positive integer to set thread count (default 1).";



    // Run base constructor, then build threading infrastructure
    wavesolvethreaded(auto &&...args): WaveOrthotope(args...), // https://tinyurl.com/byusc-parpack
            continue_iteration{true},
            nthreads{[]{ // https://tinyurl.com/byusc-lambdai
                size_type nthreads = 1;
                auto nthreads_str = std::getenv("SOLVER_NUM_THREADS");
                if (nthreads_str != nullptr) std::from_chars(nthreads_str, nthreads_str+std::strlen(nthreads_str), nthreads);
                return nthreads;
            }()},
            ds_barrier(nthreads+1),   // worker threads plus main thread
            step_barrier(nthreads+1), // worker threads plus main thread
            ds_workers(looping_threadpool(nthreads, [this](auto tid){ // https://tinyurl.com/byusc-lambda
                ds_barrier.arrive_and_wait();
                if (!continue_iteration) return false;
                auto [first, last] = this_thread_cell_range(tid);
                auto gfirst = tid==0 ? 1 : first;
                auto glast  = tid==nthreads-1 ? last-1 : last;
                value_type E = 0;
                for (size_t i = gfirst; i<glast; ++i)
                {
                for (size_t j = 1; j<cols-1; ++j)
                {
                E += (velocity(i,j) * velocity(i,j) / 2.0);
                }
                }

                for (size_t i = first; i < glast; ++i) {
                for (size_t j = 1; j < cols - 1; ++j) {
            // Vertical differences
                double n = displacement(i, j) - displacement(i + 1, j);
                E += n * n / 4.0;
                }
                }
            for (size_t i = 1; i < rows -1; ++i) {
            for (size_t j = first; j < glast; ++j) {
            // Horizontal differences

                double n = displacement(i, j) - displacement(i, j + 1);
                E += n * n / 4.0;

            }
            }
                // calculations

                this.ds_aggregator += E;
                this.ds_barrier.arrive_and_wait();
                return true;
            })),
            step_workers(looping_threadpool(nthreads, [this](auto tid){ // https://tinyurl.com/byusc-lambda
                step_barrier.arrive_and_wait();
                if (!continue_iteration) return false;
                auto [first, last] = this_thread_cell_range(tid);
                auto gfirst = tid==0 ? 1 : first;
                auto glast  = tid==nthreads-1 ? last-1 : last;
                for (size_t i = gfirst; i<glast; ++i)
                {
                for (size_t j = 1; j<cols-1; ++j)
                    {
                        L = (displacement(i-1,j) + displacement(i+1,j) + displacement(i,j-1) + displacement(i,j+1)) / 2.0 - 2.0 * displacement(i,j);
                        velocity(i,j) = (1 - dt * c) * velocity(i,j) + dt * L;
                    }
                }
                step_barrier.arrive_and_wait(); // h has to be completely updated before g update can start
                for (size_t i = gfirst; i<glast; ++i)
                {
                    for (size_t j = 1; j<cols-1; ++j)
                    {
                        displacement(i,j) += velocity(i,j) * dt;
                    }
                }
                step_barrier.arrive_and_wait();
                return true;
            }))



    // Destructor just tells threads to exit
    ~wavesolvethreaded() {
        continue_iteration = false;
        ds_barrier.arrive_and_wait();   // signal ds_workers to exit
        step_barrier.arrive_and_wait(); // signal step_workers to exit
    }



    // Steepness derivative
    value_type energy() override {
        // Reset reduction destination
        ds_aggregator = 0;

        // Launch workers
        ds_barrier.arrive_and_wait();

        // Wait for workers to finish this iteration
        ds_barrier.arrive_and_wait();

        return ds_aggregator;
    }

    // Iterate from t to t+dt in one step
    value_type step() override {
        // Let threads know what the time step this iteration is
        iter_dt = dt;

        // Signal workers to update h
        step_barrier.arrive_and_wait();

        // Signal workers to update g
        step_barrier.arrive_and_wait();

        // Wait for workers to finish this iteration
        step_barrier.arrive_and_wait();

        // Enforce boundary condition
        //g[0] = g[1];
        //g[rows-1] = g[rows-2];

        // Increment and return dt
        t += dt;
        return t;
    }
};
