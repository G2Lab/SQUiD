#include <iostream>
#include <helib/helib.h>
#include <chrono>
#include <thread>
#include <atomic>

#include "../server.hpp"
#include "../globals.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


using namespace std;

std::atomic<int> counter(0); // Shared counter
const int numThreads = 30;

const int target1 = 70000;
const int target2 = 1072820;
const int target3 = 8200000;
const int target4 = 84000000;

bool done1 = false, done2 = false, done3 = false, done4 = false;// Flags to track task completion

Server *server;

void task(std::chrono::steady_clock::time_point startTime) {
    std::vector<unsigned long> data = std::vector<unsigned long>(1000);
    std::cout << "Thread " << std::this_thread::get_id() << ": Starting...\n";
    helib::Ctxt r = server->Encrypt(data);
    counter++;


    while (!done4) {
        // Perform the task (increment the counter)
        r = server->Encrypt(data);
        counter++;

        // Check for task completion and print time
        if (counter == target1 && !done1) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();
            std::cout << "Thread " << std::this_thread::get_id() << ": Reached " << target1 << " in " << counter << " iterations. Elapsed Time: " << elapsed << " milliseconds.\n";
            done1 = true;
        }

        if (counter == target2 && !done2) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();
            std::cout << "Thread " << std::this_thread::get_id() << ": Reached " << target2 << " in " << counter << " iterations. Elapsed Time: " << elapsed << " milliseconds.\n";
            done2 = true;
        }

        if (counter == target3 && !done3) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();
            std::cout << "Thread " << std::this_thread::get_id() << ": Reached " << target3 << " in " << counter << " iterations. Elapsed Time: " << elapsed << " milliseconds.\n";
            done3 = true;
        }
        if (counter == target4 && !done4) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();
            std::cout << "Thread " << std::this_thread::get_id() << ": Reached " << target3 << " in " << counter << " iterations. Elapsed Time: " << elapsed << " milliseconds.\n";
            done4 = true;
        }
    }
}

int main()
{
    server = new Server(constants::BenchParams, false);
    std::thread threads[numThreads];

    auto startTime = std::chrono::steady_clock::now();

    // Launch threads
    for (int i = 0; i < numThreads; ++i) {
        threads[i] = std::thread(task, startTime);
    }

    // Join threads
    for (int i = 0; i < numThreads; ++i) {
        threads[i].join();
    }

    std::cout << "Counter reached " << counter << '\n';
    delete server;

    return 0;
}
