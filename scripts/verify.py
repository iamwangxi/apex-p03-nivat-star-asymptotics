#!/usr/bin/env python3
"""Re-run twelve diagnostics and compare their complete deterministic baselines."""
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
# Complete coverage arguments are linked from README.md and the proof pages.
FILE_JSON = [
    "area_and_anchor_audit", "end_profile_audit", "three_line_enumeration",
    "finite_window_audit", "full_support_lines_audit", "beta_geometry_audit",
    "nonconvex_window_audit", "two_line_product_audit", "integral_sparse_audit",
    "integral_boundary_audit",
]
NEW_JSON = {"nonconvex_window_audit", "two_line_product_audit",
            "integral_sparse_audit", "integral_boundary_audit"}


def run(name):
    proc = subprocess.run([sys.executable, "-B", str(HERE / (name + ".py"))],
                          capture_output=True, text=True,
                          env=dict(os.environ, PYTHONHASHSEED="0"))
    if proc.returncode:
        sys.stderr.write(proc.stdout[-4000:] + proc.stderr)
        raise subprocess.CalledProcessError(proc.returncode, proc.args)
    return proc.stdout


def main():
    for name in FILE_JSON:
        output = HERE / (name + ".json")
        expected = json.loads(output.read_text(encoding="utf-8"))
        stdout = run(name)
        actual = json.loads(output.read_text(encoding="utf-8"))
        if actual != expected:
            raise AssertionError("reproduction mismatch: " + output.name)
        if name in NEW_JSON and json.loads(stdout) != actual:
            raise AssertionError("stdout/file mismatch: " + output.name)
        print("MATCH " + output.name, flush=True)
    output = HERE / "presburger_examples_audit.txt"
    expected = output.read_bytes()
    run("presburger_examples_audit")
    if output.read_bytes() != expected:
        raise AssertionError("reproduction mismatch: " + output.name)
    print("MATCH " + output.name, flush=True)
    output = HERE / "minimal_period_audit.json"
    expected = json.loads(output.read_text(encoding="utf-8"))
    if json.loads(run("minimal_period_audit")) != expected:
        raise AssertionError("reproduction mismatch: " + output.name)
    print("MATCH " + output.name, flush=True)
    print("All twelve deterministic evidence outputs match.")
    print("Review scopes and source/report fingerprints: sources/review-status.md. "
          "AI adversarial review is not human certification or proof-assistant verification.")


if __name__ == "__main__":
    main()
