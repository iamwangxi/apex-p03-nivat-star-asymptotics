#!/usr/bin/env python3
"""Re-run every evidence computation and compare exact deterministic JSON."""
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
# Coverage proof for area_and_anchor_audit.py and three_line_enumeration.py:
# [Complete line enumeration](../proof/remarks-and-scope.md#11-complete-line-enumeration)
# Coverage proof for the quadrant cells in finite_window_audit.py:
# [Complete quadrant enumeration](../proof/remarks-and-scope.md#12-complete-quadrant-enumeration)
NAMES = ["area_and_anchor_audit", "end_profile_audit", "three_line_enumeration", "finite_window_audit"]


def main():
    for name in NAMES:
        output = HERE / (name + ".json")
        expected = json.loads(output.read_text(encoding="utf-8"))
        subprocess.run([sys.executable, "-B", str(HERE / (name + ".py"))], check=True)
        actual = json.loads(output.read_text(encoding="utf-8"))
        if actual != expected:
            raise AssertionError("reproduction mismatch: " + output.name)
        print("MATCH " + output.name, flush=True)
    print("All four deterministic evidence outputs match.")


if __name__ == "__main__":
    main()
