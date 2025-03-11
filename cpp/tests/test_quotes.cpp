#include <fstream>

#include "check.hpp"
#include "curvekit/quotes.hpp"

using namespace ck;

static void write(const std::string& path, const std::string& text) {
    std::ofstream f(path);
    f << text;
}

void test_quotes(Suite& suite) {
    const std::string path = "curvekit_test_quotes.txt";
    write(path, "# comment\n\ndeposit 0 0.5 0.03 1\nfra 0.5 1 0.034 1\nswap 0 5 0.039 2\n");
    const std::vector<Quote> q = read_quotes(path);
    CHECK(q.size() == 3);
    CHECK(q[0].kind == Quote::Kind::Deposit && q[1].kind == Quote::Kind::Fra && q[2].kind == Quote::Kind::Swap);
    CHECK(q[1].start == 0.5 && q[2].maturity == 5.0 && q[2].frequency == 2 && q[2].rate == 0.039);
    CHECK(bootstrap(q, Interpolation::LogLinear).max_error < 1e-13);
    write(path, "bond 0 5 0.04 2\n");
    CHECK_THROWS(read_quotes(path));
    write(path, "swap 0 five 0.04 2\n");
    CHECK_THROWS(read_quotes(path));
    write(path, "swap 3 2 0.04 2\n");
    CHECK_THROWS(read_quotes(path));
    std::remove(path.c_str());
    CHECK_THROWS(read_quotes("no_such_file.txt"));
}
