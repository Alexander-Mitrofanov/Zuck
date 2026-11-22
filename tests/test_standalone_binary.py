#!/usr/bin/env python3
"""Verify that the standalone executable has no userspace dependencies."""

from __future__ import annotations

import os
from pathlib import Path
import struct
import subprocess
import sys


PT_DYNAMIC = 2
PT_INTERP = 3


def program_header_types(executable: Path) -> set[int]:
    """Read ELF program-header types without relying on host binutils."""
    data = executable.read_bytes()
    assert data[:4] == b"\x7fELF", "standalone output is not an ELF executable"
    assert data[4] == 2, "standalone output is not 64-bit ELF"
    byte_order = "<" if data[5] == 1 else ">" if data[5] == 2 else None
    assert byte_order is not None, "standalone output has an invalid ELF byte order"

    header = struct.unpack_from(f"{byte_order}16sHHIQQQIHHHHHH", data)
    program_offset, entry_size, entry_count = header[5], header[9], header[10]
    assert entry_size >= 4 and entry_count, "standalone output has no program headers"
    assert program_offset + entry_size * entry_count <= len(data), (
        "standalone output has truncated program headers"
    )
    return {
        struct.unpack_from(f"{byte_order}I", data, program_offset + index * entry_size)[0]
        for index in range(entry_count)
    }


def run(executable: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    environment = {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "OMP_NUM_THREADS": "1",
        "PATH": "/nonexistent",
    }
    return subprocess.run(
        [str(executable), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
        timeout=30,
    )


def main() -> None:
    assert len(sys.argv) == 2, "usage: test_standalone_binary.py /path/to/Zuck"
    executable = Path(sys.argv[1]).resolve()
    assert executable.is_file() and os.access(executable, os.X_OK)

    header_types = program_header_types(executable)
    assert PT_INTERP not in header_types, "standalone binary requests a dynamic loader"
    assert PT_DYNAMIC not in header_types, "standalone binary contains a dynamic section"

    help_result = run(executable, "--help")
    assert help_result.returncode == 0, help_result.stderr
    assert "Usage:" in help_result.stdout
    assert "--input-files" in help_result.stdout

    error_result = run(executable, "--definitely-not-a-zuck-option")
    assert error_result.returncode == 2
    assert "Unknown option" in error_result.stderr


if __name__ == "__main__":
    main()
