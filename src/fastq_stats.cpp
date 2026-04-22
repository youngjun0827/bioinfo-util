#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <chrono>
#include <sys/resource.h>

#include <zlib.h>

namespace {

struct Config {
    int threads = 2;
    int phred_offset = 33;
    uint64_t genome_size = 3100000000ULL;
    bool print_header = true;
    std::vector<uint64_t> cov_cutoffs;
    std::vector<std::string> input_files;
};

struct Stats {
    std::string label;
    uint64_t total_bases = 0;
    uint64_t reads = 0;
    uint64_t gc_bases = 0;
    uint64_t acgt_bases = 0;
    uint64_t q20_bases = 0;
    uint64_t q30_bases = 0;
    uint64_t min_len = std::numeric_limits<uint64_t>::max();
    uint64_t max_len = 0;
    std::map<uint32_t, uint64_t> length_counts;
    std::vector<uint64_t> cutoff_bases;

    explicit Stats(std::string name = "") : label(std::move(name)) {}
};

[[noreturn]] void die(const std::string& msg) {
    throw std::runtime_error(msg);
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    return s.substr(start, end - start);
}

uint64_t parse_scaled_uint(const std::string& raw, const std::string& flag_name) {
    std::string s = trim(raw);
    if (s.empty()) {
        die("Empty value for " + flag_name);
    }

    char suffix = '\0';
    if (std::isalpha(static_cast<unsigned char>(s.back()))) {
        suffix = static_cast<char>(std::tolower(static_cast<unsigned char>(s.back())));
        s.pop_back();
        s = trim(s);
    }

    if (s.empty()) {
        die("Invalid numeric value for " + flag_name + ": " + raw);
    }

    char* endptr = nullptr;
    errno = 0;
    const double value = std::strtod(s.c_str(), &endptr);
    if (errno != 0 || endptr == s.c_str() || *endptr != '\0' || !std::isfinite(value) || value < 0) {
        die("Invalid numeric value for " + flag_name + ": " + raw);
    }

    double multiplier = 1.0;
    if (suffix == 'k') {
        multiplier = 1e3;
    } else if (suffix == 'm') {
        multiplier = 1e6;
    } else if (suffix == 'g') {
        multiplier = 1e9;
    } else if (suffix != '\0') {
        die("Unsupported suffix for " + flag_name + ": " + raw + " (use K, M, or G)");
    }

    const double scaled = value * multiplier;
    if (scaled > static_cast<double>(std::numeric_limits<uint64_t>::max())) {
        die("Value too large for " + flag_name + ": " + raw);
    }

    return static_cast<uint64_t>(std::llround(scaled));
}

std::string human_bp(uint64_t bp) {
    std::ostringstream oss;
    if (bp % 1000000000ULL == 0) {
        oss << (bp / 1000000000ULL) << "g";
    } else if (bp % 1000000ULL == 0) {
        oss << (bp / 1000000ULL) << "m";
    } else if (bp % 1000ULL == 0) {
        oss << (bp / 1000ULL) << "k";
    } else {
        oss << bp;
    }
    return oss.str();
}

void print_help() {
    std::cerr
        << "Usage: fastq_stats [options] file1.fastq[.gz] file2.fastq[.gz] ...\n\n"
        << "Options:\n"
        << "  -t, --threads INT         Number of files to process in parallel [default: 2]\n"
        << "      --phred INT           Quality score offset: 33 or 64 [default: 33]\n"
        << "  -g, --genome SIZE         Genome size for coverage, supports K/M/G [default: 3.1G]\n"
        << "      --cov-cutoff SIZE     Additional coverage cutoff, repeatable, supports K/M/G\n"
        << "      --no-header           Do not print header\n"
        << "  -h, --help                Show this help message\n";
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            std::exit(0);
        }
        if (arg == "-t" || arg == "--threads") {
            if (i + 1 >= argc) {
                die("Missing value for " + arg);
            }
            cfg.threads = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--phred") {
            if (i + 1 >= argc) {
                die("Missing value for --phred");
            }
            cfg.phred_offset = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "-g" || arg == "--genome") {
            if (i + 1 >= argc) {
                die("Missing value for " + arg);
            }
            cfg.genome_size = parse_scaled_uint(argv[++i], arg);
            continue;
        }
        if (arg == "--cov-cutoff") {
            if (i + 1 >= argc) {
                die("Missing value for --cov-cutoff");
            }
            cfg.cov_cutoffs.push_back(parse_scaled_uint(argv[++i], "--cov-cutoff"));
            continue;
        }
        if (arg == "--no-header") {
            cfg.print_header = false;
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            die("Unknown option: " + arg);
        }
        cfg.input_files.push_back(arg);
    }

    if (cfg.input_files.empty()) {
        print_help();
        std::exit(1);
    }
    if (cfg.threads < 1) {
        die("--threads must be >= 1");
    }
    if (cfg.phred_offset != 33 && cfg.phred_offset != 64) {
        die("--phred must be 33 or 64");
    }
    if (cfg.genome_size == 0) {
        die("--genome must be > 0");
    }

    std::sort(cfg.cov_cutoffs.begin(), cfg.cov_cutoffs.end());
    cfg.cov_cutoffs.erase(std::unique(cfg.cov_cutoffs.begin(), cfg.cov_cutoffs.end()), cfg.cov_cutoffs.end());

    return cfg;
}

bool is_fastq_header(const std::string& line) {
    return !line.empty() && line[0] == '@';
}

void update_stats_with_record(
    const std::string& seq,
    const std::string& qual,
    int phred_offset,
    const std::vector<uint64_t>& cov_cutoffs,
    Stats& stats
) {
    const uint64_t len = static_cast<uint64_t>(seq.size());
    if (len != qual.size()) {
        die("Sequence and quality lengths do not match in " + stats.label);
    }

    ++stats.reads;
    stats.total_bases += len;
    stats.min_len = std::min(stats.min_len, len);
    stats.max_len = std::max(stats.max_len, len);
    ++stats.length_counts[static_cast<uint32_t>(len)];

    for (size_t i = 0; i < cov_cutoffs.size(); ++i) {
        if (len >= cov_cutoffs[i]) {
            stats.cutoff_bases[i] += len;
        }
    }

    for (char c : seq) {
        switch (std::toupper(static_cast<unsigned char>(c))) {
            case 'G':
            case 'C':
                ++stats.gc_bases;
                ++stats.acgt_bases;
                break;
            case 'A':
            case 'T':
                ++stats.acgt_bases;
                break;
            default:
                break;
        }
    }

    for (char c : qual) {
        const int q = static_cast<int>(static_cast<unsigned char>(c)) - phred_offset;
        if (q >= 20) {
            ++stats.q20_bases;
        }
        if (q >= 30) {
            ++stats.q30_bases;
        }
    }
}

Stats process_plain_fastq(const std::string& path, int phred_offset, const std::vector<uint64_t>& cov_cutoffs) {
    std::ifstream in(path);
    if (!in) {
        die("Failed to open file: " + path);
    }

    Stats stats(path);
    stats.cutoff_bases.assign(cov_cutoffs.size(), 0);

    std::string h;
    std::string seq;
    std::string plus;
    std::string qual;

    while (std::getline(in, h)) {
        if (!std::getline(in, seq) || !std::getline(in, plus) || !std::getline(in, qual)) {
            die("Truncated FASTQ record in " + path);
        }
        if (!is_fastq_header(h)) {
            die("Invalid FASTQ header in " + path + ": " + h);
        }
        if (plus.empty() || plus[0] != '+') {
            die("Invalid FASTQ separator line in " + path);
        }
        update_stats_with_record(seq, qual, phred_offset, cov_cutoffs, stats);
    }

    return stats;
}

Stats process_gzip_fastq(const std::string& path, int phred_offset, const std::vector<uint64_t>& cov_cutoffs) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) {
        die("Failed to open gzip file: " + path);
    }

    constexpr int kBufferSize = 1 << 20;
    std::vector<char> buffer(kBufferSize);
    std::string pending;
    std::vector<std::string> lines;
    lines.reserve(4);

    Stats stats(path);
    stats.cutoff_bases.assign(cov_cutoffs.size(), 0);

    auto consume_line = [&](std::string&& line) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        lines.push_back(std::move(line));
        if (lines.size() == 4) {
            if (!is_fastq_header(lines[0])) {
                gzclose(gz);
                die("Invalid FASTQ header in " + path + ": " + lines[0]);
            }
            if (lines[2].empty() || lines[2][0] != '+') {
                gzclose(gz);
                die("Invalid FASTQ separator line in " + path);
            }
            update_stats_with_record(lines[1], lines[3], phred_offset, cov_cutoffs, stats);
            lines.clear();
        }
    };

    while (true) {
        char* result = gzgets(gz, buffer.data(), kBufferSize);
        if (result == nullptr) {
            break;
        }
        std::string chunk(result);
        if (!chunk.empty() && chunk.back() == '\n') {
            pending += chunk.substr(0, chunk.size() - 1);
            consume_line(std::move(pending));
            pending.clear();
        } else {
            pending += chunk;
        }
    }

    if (!pending.empty()) {
        consume_line(std::move(pending));
        pending.clear();
    }
    if (!lines.empty()) {
        gzclose(gz);
        die("Truncated FASTQ record in " + path);
    }

    gzclose(gz);
    return stats;
}

Stats process_file(const std::string& path, int phred_offset, const std::vector<uint64_t>& cov_cutoffs) {
    if (ends_with(path, ".gz")) {
        return process_gzip_fastq(path, phred_offset, cov_cutoffs);
    }
    return process_plain_fastq(path, phred_offset, cov_cutoffs);
}

void merge_stats(Stats& total, Stats&& part) {
    total.total_bases += part.total_bases;
    total.reads += part.reads;
    total.gc_bases += part.gc_bases;
    total.acgt_bases += part.acgt_bases;
    total.q20_bases += part.q20_bases;
    total.q30_bases += part.q30_bases;
    total.max_len = std::max(total.max_len, part.max_len);
    if (part.reads > 0) {
        total.min_len = std::min(total.min_len, part.min_len);
    }

    if (total.cutoff_bases.empty()) {
        total.cutoff_bases.assign(part.cutoff_bases.size(), 0);
    }
    for (size_t i = 0; i < part.cutoff_bases.size(); ++i) {
        total.cutoff_bases[i] += part.cutoff_bases[i];
    }

    for (const auto& kv : part.length_counts) {
        total.length_counts[kv.first] += kv.second;
    }
}

uint64_t compute_median(const std::map<uint32_t, uint64_t>& length_counts, uint64_t total_reads) {
    if (total_reads == 0 || length_counts.empty()) {
        return 0;
    }

    const uint64_t lower_rank = (total_reads - 1) / 2;
    const uint64_t upper_rank = total_reads / 2;
    uint64_t seen = 0;
    uint64_t lower_value = 0;
    uint64_t upper_value = 0;
    bool lower_found = false;

    for (const auto& kv : length_counts) {
        const uint64_t next_seen = seen + kv.second;
        if (!lower_found && lower_rank < next_seen) {
            lower_value = kv.first;
            lower_found = true;
        }
        if (upper_rank < next_seen) {
            upper_value = kv.first;
            break;
        }
        seen = next_seen;
    }

    return (lower_value + upper_value) / 2;
}

uint64_t compute_n50(const std::map<uint32_t, uint64_t>& length_counts, uint64_t total_bases) {
    if (total_bases == 0 || length_counts.empty()) {
        return 0;
    }

    const uint64_t threshold = (total_bases + 1) / 2;
    uint64_t cumulative = 0;
    for (auto it = length_counts.rbegin(); it != length_counts.rend(); ++it) {
        cumulative += static_cast<uint64_t>(it->first) * it->second;
        if (cumulative >= threshold) {
            return it->first;
        }
    }
    return 0;
}

std::string format_double(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << value;
    return oss.str();
}

void print_header(const std::vector<uint64_t>& cov_cutoffs) {
    std::cout
        << "file\tyields\treads\tgc_contents\tratio_Q20\tratio_Q30"
        << "\tmin_len\tavg_len\tmed_len\tmax_len\tn50\tcoverage";
    for (uint64_t cutoff : cov_cutoffs) {
        std::cout << "\tcoverage_" << human_bp(cutoff) << "+";
    }
    std::cout << '\n';
}

void print_stats_line(const Stats& stats, uint64_t genome_size) {
    const double gc_pct = stats.acgt_bases > 0
        ? 100.0 * static_cast<double>(stats.gc_bases) / static_cast<double>(stats.acgt_bases)
        : 0.0;
    const double q20_pct = stats.total_bases > 0
        ? 100.0 * static_cast<double>(stats.q20_bases) / static_cast<double>(stats.total_bases)
        : 0.0;
    const double q30_pct = stats.total_bases > 0
        ? 100.0 * static_cast<double>(stats.q30_bases) / static_cast<double>(stats.total_bases)
        : 0.0;
    const double avg_len = stats.reads > 0
        ? static_cast<double>(stats.total_bases) / static_cast<double>(stats.reads)
        : 0.0;
    const uint64_t med_len = compute_median(stats.length_counts, stats.reads);
    const uint64_t n50 = compute_n50(stats.length_counts, stats.total_bases);
    const double coverage = static_cast<double>(stats.total_bases) / static_cast<double>(genome_size);

    std::cout
        << stats.label << '\t'
        << stats.total_bases << '\t'
        << stats.reads << '\t'
        << format_double(gc_pct) << '\t'
        << format_double(q20_pct) << '\t'
        << format_double(q30_pct) << '\t'
        << (stats.reads > 0 ? stats.min_len : 0) << '\t'
        << format_double(avg_len) << '\t'
        << med_len << '\t'
        << stats.max_len << '\t'
        << n50 << '\t'
        << format_double(coverage);

    for (uint64_t cutoff_bases : stats.cutoff_bases) {
        const double cutoff_cov = static_cast<double>(cutoff_bases) / static_cast<double>(genome_size);
        std::cout << '\t' << format_double(cutoff_cov);
    }
    std::cout << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const auto t_start = std::chrono::steady_clock::now();
        const Config cfg = parse_args(argc, argv);

        std::vector<Stats> per_file(cfg.input_files.size());
        std::atomic<size_t> next_index{0};
        const size_t worker_count = std::min(static_cast<size_t>(cfg.threads), cfg.input_files.size());
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        std::mutex error_mutex;
        std::string first_error;

        for (size_t t = 0; t < worker_count; ++t) {
            workers.emplace_back([&]() {
                while (true) {
                    const size_t idx = next_index.fetch_add(1);
                    if (idx >= cfg.input_files.size()) {
                        return;
                    }
                    try {
                        per_file[idx] = process_file(cfg.input_files[idx], cfg.phred_offset, cfg.cov_cutoffs);
                    } catch (const std::exception& e) {
                        std::lock_guard<std::mutex> lock(error_mutex);
                        if (first_error.empty()) {
                            first_error = e.what();
                        }
                        return;
                    }
                }
            });
        }

        for (auto& worker : workers) {
            worker.join();
        }

        if (!first_error.empty()) {
            die(first_error);
        }

        if (cfg.print_header) {
            print_header(cfg.cov_cutoffs);
        }

        Stats total("total");
        total.cutoff_bases.assign(cfg.cov_cutoffs.size(), 0);
        for (const Stats& stats : per_file) {
            print_stats_line(stats, cfg.genome_size);
        }
        if (per_file.size() > 1) {
            for (Stats& stats : per_file) {
                merge_stats(total, std::move(stats));
            }
            print_stats_line(total, cfg.genome_size);
        }

        {
            const auto t_end = std::chrono::steady_clock::now();
            const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();
            struct rusage ru;
            getrusage(RUSAGE_SELF, &ru);
            std::cerr << std::fixed << std::setprecision(2)
                      << "elapsed: " << elapsed_s << " s"
                      << "  peak_rss: " << ru.ru_maxrss << " KB\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}

