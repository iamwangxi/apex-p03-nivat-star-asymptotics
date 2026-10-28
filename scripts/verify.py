#!/usr/bin/env python3
"""Re-run every evidence computation and compare exact deterministic JSON."""
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
# Coverage proof for area_and_anchor_audit.py and three_line_enumeration.py:
# [Complete line enumeration](../proof/remarks-and-scope.md#11-complete-line-enumeration)
# Coverage proof for the quadrant cells in finite_window_audit.py:
# [Complete quadrant enumeration](../proof/remarks-and-scope.md#12-complete-quadrant-enumeration)
NAMES = ["area_and_anchor_audit", "end_profile_audit", "three_line_enumeration", "finite_window_audit", "full_support_lines_audit", "beta_geometry_audit"]


def main():
    for name in NAMES:
        output = HERE / (name + ".json")
        expected = json.loads(output.read_text(encoding="utf-8"))
        subprocess.run([sys.executable, "-B", str(HERE / (name + ".py"))], check=True)
        actual = json.loads(output.read_text(encoding="utf-8"))
        if actual != expected:
            raise AssertionError("reproduction mismatch: " + output.name)
        print("MATCH " + output.name, flush=True)
    output = HERE / "presburger_examples_audit.txt"
    expected = output.read_bytes()
    subprocess.run([sys.executable, "-B", str(HERE / "presburger_examples_audit.py")], check=True)
    if output.read_bytes() != expected:
        raise AssertionError("reproduction mismatch: " + output.name)
    print("MATCH " + output.name, flush=True)
    output = HERE / "minimal_period_audit.json"
    expected = json.loads(output.read_text(encoding="utf-8"))
    proc = subprocess.run([sys.executable, "-B", str(HERE / "minimal_period_audit.py")], check=True,
                          capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED="0"))
    if json.loads(proc.stdout) != expected:
        raise AssertionError("reproduction mismatch: " + output.name)
    print("MATCH " + output.name, flush=True)
    print("All eight deterministic evidence outputs match.")
    print("The explicit-coefficient supplement and the minimal-period page each passed one fresh-context GPT-6.1 Sol adversarial review (2026-10-25 and 2026-10-26); no fatal or proof-repair issue was found. AI review is not human certification or proof-assistant verification.")


if __name__ == "__main__":
    main()
