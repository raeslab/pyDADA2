"""K7 primer matching vs Biostrings::matchPattern on 300 random R cases (match_pattern.json).

Each case is checked in mismatch mode and in indel mode (600 checks); match positions must be
identical (R start is 1-based inclusive, end inclusive; pydada2 returns 0-based start, end
exclusive).
"""

import pytest

import pydada2._core as core

from rparity import load, scalar

pytestmark = pytest.mark.parity

CASES = load("match_pattern")


@pytest.mark.parametrize("mode", ["mismatch", "indels"])
@pytest.mark.parametrize("i", range(len(CASES)))
def test_match_pattern_parity(i, mode):
    case = CASES[i]
    exp = case[mode]
    assert exp is not None, "R errored on this case; regenerate the fixture"
    expected = [(s - 1, e) for s, e in zip(exp["start"], exp["end"], strict=True)]
    got = core.match_pattern(
        scalar(case["subject"]),
        scalar(case["pattern"]),
        scalar(case["max_mismatch"]),
        mode == "indels",
        scalar(case["fixed"]),
    )
    assert [tuple(x) for x in got] == expected
