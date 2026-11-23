cmake_minimum_required(VERSION 3.14)
project(Zuck VERSION 6.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(STATIC_BUILD "Build a dependency-free static executable" OFF)
option(PORTABLE_BUILD "Avoid host-specific CPU instructions" OFF)
option(SANITIZER "Enable address/UB sanitizers" OFF)
option(ENABLE_DEBUG "Enable DEBUG diagnostics" OFF)
option(BUILD_TESTING "Build Zuck tests" OFF)

# Library selection must be configured before dependency discovery. Otherwise
# STATIC_BUILD can silently retain an absolute path to a shared libz.
if(STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(ZLIB_USE_STATIC_LIBS ON)
endif()

# ── Dependencies ──────────────────────────────────────────────────────────────

find_package(ZLIB)
if(NOT ZLIB_FOUND)
    include(FetchContent)
    FetchContent_Declare(zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1)
    FetchContent_MakeAvailable(zlib)
    set(ZLIB_LIBRARIES zlibstatic)
    set(ZLIB_INCLUDE_DIRS ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()

include(FetchContent)

FetchContent_Declare(phmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG v1.3.11)
FetchContent_MakeAvailable(phmap)

FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1)
FetchContent_MakeAvailable(fmt)

find_package(OpenMP REQUIRED)

# ── Options ───────────────────────────────────────────────────────────────────

if(ENABLE_DEBUG)
    add_definitions(-DDEBUG)
    include(FetchContent)
    FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3)
    FetchContent_MakeAvailable(json)
endif()

# ── Compiler flags ────────────────────────────────────────────────────────────

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(PORTABLE_BUILD)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
endif()
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb -O1 -D_GLIBCXX_DEBUG")

if(SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined")
endif()

# ── Sub-libraries ─────────────────────────────────────────────────────────────

set(MEGAHIT_SRC "${CMAKE_SOURCE_DIR}/libs/megahit/src")

# The vendored kseqpp snapshot does not include its optional test directory.
# Keep its tests disabled while retaining the top-level BUILD_TESTING choice.
set(ZUCK_SAVED_BUILD_TESTING ${BUILD_TESTING})
set(BUILD_TESTING OFF)
add_subdirectory(libs/kseqpp)
set(BUILD_TESTING ${ZUCK_SAVED_BUILD_TESTING})

set(spoa_build_exe   OFF CACHE BOOL "" FORCE)
set(spoa_build_tests OFF CACHE BOOL "" FORCE)
set(spoa_use_cereal  OFF CACHE BOOL "" FORCE)
set(spoa_use_simde   OFF CACHE BOOL "" FORCE)
if(PORTABLE_BUILD)
    # SPOA otherwise appends its own -march=native after the top-level release
    # flags have been made portable.
    set(spoa_optimize_for_native OFF CACHE BOOL "" FORCE)
endif()
add_subdirectory(libs/spoa)

# ── Include paths ─────────────────────────────────────────────────────────────

include_directories(include)
include_directories(${MEGAHIT_SRC})
include_directories(libs/cft/src)
include_directories(libs/kseqpp/include)
include_directories(libs/spoa/include)
include_directories(${phmap_SOURCE_DIR})

# ── Source globs ──────────────────────────────────────────────────────────────

file(GLOB_RECURSE SOURCES_WITHOUT_MAIN CONFIGURE_DEPENDS "src/*.cpp")
list(REMOVE_ITEM SOURCES_WITHOUT_MAIN "${CMAKE_SOURCE_DIR}/src/main.cpp")

file(GLOB_RECURSE SDBG_SOURCE    "${MEGAHIT_SRC}/sdbg/*.cpp")
file(GLOB_RECURSE ASMBL_SOURCE   "${MEGAHIT_SRC}/assembly/*.cpp")
file(GLOB_RECURSE IDBA_SOURCE    "${MEGAHIT_SRC}/idba/*.cpp")
file(GLOB_RECURSE SEQ_SOURCE     "${MEGAHIT_SRC}/sequence/*.cpp")
file(GLOB_RECURSE CX1_SOURCE     "${MEGAHIT_SRC}/sorting/*.cpp")
file(GLOB_RECURSE TOOLKIT_SOURCE "${MEGAHIT_SRC}/tools/*.cpp")

set(MEGAHIT_EXTRAS
    ${MEGAHIT_SRC}/main_assemble.cpp
    ${MEGAHIT_SRC}/main_buildlib.cpp
    ${MEGAHIT_SRC}/main_iterate.cpp
    ${MEGAHIT_SRC}/main_sdbg_build.cpp
    ${MEGAHIT_SRC}/utils/options_description.cpp
)

set(MEGAHIT_SOURCES
    ${SDBG_SOURCE}
    ${ASMBL_SOURCE}
    ${IDBA_SOURCE}
    ${SEQ_SOURCE}
    ${CX1_SOURCE}
    ${TOOLKIT_SOURCE}
    ${MEGAHIT_EXTRAS}
)

# ── Executable: Zuck ──────────────────────────────────────────────────────

# Compile shared production code once; helper tests link the same objects.
add_library(zuck_core STATIC ${SOURCES_WITHOUT_MAIN} ${MEGAHIT_SOURCES})
add_executable(Zuck src/main.cpp)
target_link_libraries(Zuck PRIVATE zuck_core)

# Identify the compiled source even when the checkout contains local edits.
# Reconfigure after source edits so the identifier cannot silently go stale.
file(GLOB_RECURSE IDENTITY_SOURCES CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/*.h" "src/*.cpp" "libs/*.h" "libs/*.hpp" "libs/*.c" "libs/*.cpp"
    "libs/*CMakeLists.txt" "cmake/*.in")
list(APPEND IDENTITY_SOURCES "CMakeLists.txt")
list(SORT IDENTITY_SOURCES)
set(ZUCK_SOURCE_IDENTITY "")
foreach(source IN LISTS IDENTITY_SOURCES)
    file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/${source}" source_hash)
    string(APPEND ZUCK_SOURCE_IDENTITY "${source}:${source_hash}\n")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
endforeach()
string(SHA256 ZUCK_SOURCE_SHA256 "${ZUCK_SOURCE_IDENTITY}")
configure_file(cmake/build_identity.h.in generated/Zuck_build_identity.h @ONLY)
target_include_directories(zuck_core PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/generated" ${ZLIB_INCLUDE_DIRS})

target_link_libraries(zuck_core PUBLIC
    ${ZLIB_LIBRARIES}
    spoa::spoa
    fmt::fmt
    OpenMP::OpenMP_CXX
)

include(CTest)
if(BUILD_TESTING)
    # Each test uses production headers; C++ helper tests can also link zuck_core.
    function(zuck_cpp_test name)
        add_executable(test_${name} tests/test_${name}.cpp)
        target_link_libraries(test_${name} PRIVATE zuck_core)
        add_test(NAME ${name} COMMAND test_${name})
    endfunction()
    foreach(name IN ITEMS array_database database_reads database_groups
            element_confidence supported_calls read_recruitment boundary_candidates bounded_cycle_search alternative_repeat_view
            local_boundary_graph motif_cliff_boundary directed_haplotype_paths cross_variant_scaffolds
            short_spacer_filter strand_equivalent_groups terminal_variant_rescue
            read_linked_groups algorithm_helpers spacer_ordering)
        zuck_cpp_test(${name})
    endforeach()

    if(ENABLE_DEBUG)
        add_executable(test_debug_io tests/test_debug_io.cpp)
        target_link_libraries(test_debug_io PRIVATE zuck_core)
        add_test(NAME debug_io COMMAND test_debug_io
            ${CMAKE_CURRENT_BINARY_DIR}/debug-cycle-test.json)
    endif()

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        function(zuck_python_test name script)
            add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/${script} ${ARGN})
        endfunction()
        foreach(name IN ITEMS annotations build_array_database cycle_enumeration search_partitions experimental_paths benchmark_runner)
            zuck_python_test(${name} test_${name}.py)
        endforeach()
        foreach(name IN ITEMS array_database_cli element_confidence_cli supported_calls_cli
                family_repeat_context_cli read_recruitment_cli
                repeat_mixture_refinement_cli
                two_sided_boundaries graph_settings_cli candidate_evidence_cli
                length_constraints_cli result_bundle strand_coalescing_cli
                terminal_variant_rescue_cli read_linked_groups_cli settings_cli)
            zuck_python_test(${name} test_${name}.py $<TARGET_FILE:Zuck>)
        endforeach()
        zuck_python_test(cli_io_errors test_cli_io.py $<TARGET_FILE:Zuck>)
        zuck_python_test(short_spacer_filter_cli test_short_spacer_filter.py $<TARGET_FILE:Zuck>)
        if(STATIC_BUILD)
            zuck_python_test(standalone_binary test_standalone_binary.py $<TARGET_FILE:Zuck>)
        endif()
    endif()
endif()

if(ENABLE_DEBUG)
    target_link_libraries(zuck_core PUBLIC nlohmann_json::nlohmann_json)
    target_compile_definitions(Zuck PRIVATE DEBUG)
endif()

if(STATIC_BUILD)
    set_target_properties(Zuck PROPERTIES
        LINK_SEARCH_START_STATIC ON
        LINK_SEARCH_END_STATIC ON)
    target_link_options(Zuck PRIVATE -static)
endif()

