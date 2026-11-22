#include "io/io_ops.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
    require(argc == 2, "Provide a temporary JSON path");
    const std::string path(argv[1]);
    const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> expected{
        {0, {{0, 4, UINT64_MAX}, {9}}}, {8, {}}};
    io_ops::write_cycles(path, expected);
    require(io_ops::read_cycles(path) == expected, "Cycle JSON must round trip full-width IDs");
    io_ops::write_cycles(path, {});
    require(io_ops::read_cycles(path).empty(), "Empty cycle JSON must be an object");
    for (const auto* bad : {"[]", "null", "{\"-1\":[]}", "{\"4junk\":[]}",
         "{\"1\":{}}", "{\"1\":[4]}", "{\"1\":[[-1]]}", "{\"1\":[[1.0]]}",
         "{\"1\":[[\"4\"]]}", "{\"1\":[],\"01\":[]}", "{} trailing"}) {
        { std::ofstream file(path); file << bad; }
        bool rejected = false;
        try { io_ops::read_cycles(path); } catch (const std::exception&) { rejected = true; }
        require(rejected, "Malformed cycle JSON must fail");
    }
    std::filesystem::remove(path);
}
