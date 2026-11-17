#!/usr/bin/env python3
"""Print the failing/erroring tests of testthat JUnit files as markdown.

    python3 .github/scripts/junit-failures.py junit.xml [more.xml ...]

Used for a shard's job summary and for the CI failure issue on main.
Prints nothing when every test passed.
"""
import sys
import xml.etree.ElementTree as ET

rows = []
for path in sys.argv[1:]:
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as e:
        rows.append(f"- `{path}`: could not read JUnit results ({e})")
        continue
    for suite in root.iter("testsuite"):
        for case in suite.iter("testcase"):
            for kind in ("failure", "error"):
                bad = case.find(kind)
                if bad is None:
                    continue
                msg = (bad.get("message") or bad.text or "").strip().splitlines()
                msg = msg[0][:200] if msg else ""
                rows.append(f"- `test-{suite.get('name')}.R` {case.get('name')} ({kind}): {msg}")

if rows:
    print(f"### Failing tests ({len(rows)})\n")
    print("\n".join(rows[:100]))
    if len(rows) > 100:
        print(f"\n... and {len(rows) - 100} more")
