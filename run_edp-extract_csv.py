"""Run a FreeFem++ input file and export its iteration log to CSV."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[+-]?(?:nan|inf(?:inity)?)"
CONTINUATION_PATTERN = re.compile(
    rf"Continuation\s*:\s*n\s*=\s*({NUMBER})", re.IGNORECASE
)
NEWTON_ITERATION_PATTERN = re.compile(
    rf"Newton iteration\s+(\d+)\s*:\s*\|\|h\|\|\s*=\s*({NUMBER})",
    re.IGNORECASE,
)
FIXED_POINT_ITERATION_PATTERN = re.compile(
    rf"Iteration\s+(\d+)\s*:\s*errW\s*=\s*({NUMBER})"
    rf"\s*,\s*errMu\s*=\s*({NUMBER})",
    re.IGNORECASE,
)


def find_freefem(executable: str | None) -> str:
    if executable:
        candidate = Path(executable).expanduser()
        if candidate.is_file():
            return str(candidate.resolve())
        found = shutil.which(executable)
        if found:
            return found
        raise FileNotFoundError(
            f"FreeFem++ executable not found: {executable!r}. "
            "Pass its path with --freefem."
        )

    for name in ("FreeFem++", "FreeFem++-nw"):
        found = shutil.which(name)
        if found:
            return found

    if os.name == "nt":
        program_files = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"))
        candidate = program_files / "FreeFem++" / "FreeFem++.exe"
        if candidate.is_file():
            return str(candidate)

    raise FileNotFoundError(
        "FreeFem++ was not found on PATH. Pass its path with --freefem."
    )


def parse_iteration_log(output: str) -> list[dict[str, int | float | str]]:
    rows: list[dict[str, int | float | str]] = []
    stage: int | str = ""
    current_n: float | str = ""

    for line in output.splitlines():
        continuation = CONTINUATION_PATTERN.search(line)
        if continuation:
            stage = 0 if stage == "" else int(stage) + 1
            current_n = float(continuation.group(1))
            continue

        newton_iteration = NEWTON_ITERATION_PATTERN.search(line)
        if newton_iteration:
            rows.append(
                {
                    "method": "Newton",
                    "stage": stage,
                    "n": current_n,
                    "iteration": int(newton_iteration.group(1)),
                    "errH": float(newton_iteration.group(2)),
                    "errW": "",
                    "errMu": "",
                }
            )
            continue

        fixed_point_iteration = FIXED_POINT_ITERATION_PATTERN.search(line)
        if fixed_point_iteration:
            rows.append(
                {
                    "method": "FixedPoint",
                    "stage": stage,
                    "n": current_n,
                    "iteration": int(fixed_point_iteration.group(1)),
                    "errH": "",
                    "errW": float(fixed_point_iteration.group(2)),
                    "errMu": float(fixed_point_iteration.group(3)),
                }
            )

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run a FreeFem++ .edp file and save its iteration log as CSV."
    )
    default_edp = Path(__file__).with_name("methode_newton.edp")
    parser.add_argument(
        "--edp",
        type=Path,
        default=default_edp,
        help="FreeFem++ input file (default: methode_newton.edp beside this script)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help="CSV output path (default: iteration_history.csv beside the .edp file)",
    )
    parser.add_argument(
        "--freefem",
        help="FreeFem++ executable path (default: search PATH and standard Windows location)",
    )
    parser.add_argument(
        "--log",
        type=Path,
        help="Optional path for saving the complete FreeFem++ console output",
    )
    args = parser.parse_args()

    edp_path = args.edp.expanduser().resolve()
    if not edp_path.is_file():
        parser.error(f"FreeFem++ input file does not exist: {edp_path}")

    try:
        executable = find_freefem(args.freefem)
    except FileNotFoundError as error:
        parser.error(str(error))

    csv_path = (
        args.csv.expanduser().resolve()
        if args.csv
        else edp_path.with_name("iteration_history.csv")
    )
    command = [executable, "-nw", edp_path.name]
    result = subprocess.run(
        command,
        cwd=edp_path.parent,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )

    if args.log:
        log_path = args.log.expanduser().resolve()
        log_path.write_text(result.stdout, encoding="utf-8")

    rows = parse_iteration_log(result.stdout)
    if rows:
        with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(
                csv_file,
                fieldnames=(
                    "method",
                    "stage",
                    "n",
                    "iteration",
                    "errH",
                    "errW",
                    "errMu",
                ),
            )
            writer.writeheader()
            writer.writerows(rows)

        print(f"Saved {len(rows)} iteration records to {csv_path}")

    if args.log:
        print(f"Saved console output to {args.log.expanduser().resolve()}")

    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-20:])
        if tail:
            sys.stderr.write(f"\n{tail}\n")
        print(
            f"FreeFem++ exited with status {result.returncode}.",
            file=sys.stderr,
        )
        return result.returncode

    if not rows:
        sys.stderr.write(result.stdout)
        print(
            "No Newton iteration records were found in the FreeFem++ output.",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
