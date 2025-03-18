#include <cstdio>

#include "check.hpp"

int main() {
    Suite suite;
    struct Entry {
        const char* name;
        void (*run)(Suite&);
    };
    const Entry entries[] = {{"root", test_root},           {"linsolve", test_linsolve},   {"interp", test_interp},           {"curve", test_curve},         {"swap", test_swap},
                             {"bootstrap", test_bootstrap}, {"bond", test_bond},           {"nelson_siegel", test_nelson_siegel}, {"black", test_black},     {"shortrate", test_shortrate},
                             {"hullwhite", test_hullwhite}, {"hwtree", test_hwtree},       {"calibration", test_calibration}, {"quotes", test_quotes}};
    for (const Entry& e : entries) {
        const int before = suite.failures, checks = suite.checks;
        e.run(suite);
        std::printf("%-14s %3d checks  %s\n", e.name, suite.checks - checks, suite.failures == before ? "ok" : "FAILED");
    }
    std::printf("%d checks, %d failures\n", suite.checks, suite.failures);
    return suite.failures == 0 ? 0 : 1;
}
