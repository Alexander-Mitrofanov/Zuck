#include "core/cli.h"
#include "core/pipeline.h"

#include <chrono>
#include <iomanip>
#include <iostream>

using namespace std;
int main(int argc, char** argv) {
    try {
        Settings settings = parse_arguments(argc, argv);
        auto t_start = chrono::steady_clock::now();
        auto step_elapsed = [](chrono::steady_clock::time_point t0) {
            return chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        };

        if (!settings.cycles_input.empty()) {
            auto t0 = chrono::steady_clock::now();
            cout << "[1/1] Post-processing retained cycles...\n";
            step_post_process(settings);
            cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n"
                 << "Output: " << settings.output_folder << "\n\n";
            return 0;
        }

        SDBG sdbg;
        auto t0 = chrono::steady_clock::now();
        if (settings.graph_input.empty()) {
            cout << "[1/4] Building de Bruijn graph...\n";
            step_build_graph(settings);
        } else {
            cout << "[1/4] Using pre-built graph\n";
        }
        step_load_graph(settings, sdbg);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        t0 = chrono::steady_clock::now();
        cout << "[2/4] Finding CRISPR cycles...\n";
        step_find_cycles(settings);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        t0 = chrono::steady_clock::now();
        cout << "[3/4] Post-processing arrays...\n";
        step_post_process(settings);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        cout << "[4/4] Cleaning up...\n";
        step_cleanup(settings);

        double total = step_elapsed(t_start);
        cout << "\nDone in " << fixed << setprecision(1) << total << " s\n"
             << "Output: " << settings.output_folder << "\n\n";

        return 0;
    } catch (const std::exception& error) {
        cerr << error.what() << "\n";
        return 2;
    }
}
