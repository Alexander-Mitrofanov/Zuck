#ifndef ZUCK_CLI_H
#define ZUCK_CLI_H

struct Settings;

// Resolve file defaults, CLI overrides, input validation and output paths.
Settings parse_arguments(int argc, char* argv[]);

#endif
