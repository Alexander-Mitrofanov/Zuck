#ifndef ZUCK_CORE_BUILD_IDENTITY_H
#define ZUCK_CORE_BUILD_IDENTITY_H

#if __has_include("Zuck_build_identity.h")
#include "Zuck_build_identity.h"
#else
#define ZUCK_VERSION "unknown (built without CMake identity)"
#define ZUCK_SOURCE_SHA256 "unknown (built without CMake identity)"
#endif

#endif
