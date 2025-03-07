#include "curvekit/quotes.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ck {

std::vector<Quote> read_quotes(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<Quote> quotes;
    std::string line;
    int number = 0;
    while (std::getline(in, line)) {
        ++number;
        if (line.empty() || line[0] == '#' || line == "\r") continue;
        std::istringstream ss(line);
        std::string kind;
        Quote q;
        if (!(ss >> kind >> q.start >> q.maturity >> q.rate >> q.frequency)) throw std::runtime_error("line " + std::to_string(number) + " of " + path + " is malformed");
        if (kind == "deposit") q.kind = Quote::Kind::Deposit;
        else if (kind == "fra") q.kind = Quote::Kind::Fra;
        else if (kind == "swap") q.kind = Quote::Kind::Swap;
        else throw std::runtime_error("line " + std::to_string(number) + ": unknown instrument " + kind);
        if (!(q.maturity > q.start) || q.frequency <= 0) throw std::runtime_error("line " + std::to_string(number) + ": invalid times or frequency");
        quotes.push_back(q);
    }
    return quotes;
}

}  // namespace ck
