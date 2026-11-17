# -*- coding: utf-8 -*-
"""Tests for the UTC timestamp helpers used by the AI layer."""

import datetime
import re

from metharct.ai.clock import (
    ISO_FORMAT,
    format_iso,
    parse_iso,
    timestamp_to_iso,
    utc_now,
    utc_now_iso,
)

ISO_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


def test_utc_now_is_timezone_aware():
    moment = utc_now()
    assert moment.tzinfo is not None
    assert moment.utcoffset() == datetime.timedelta(0)


def test_utc_now_iso_matches_the_documented_format():
    text = utc_now_iso()
    assert ISO_PATTERN.match(text), text
    assert parse_iso(text) is not None


def test_format_iso_accepts_naive_datetimes_as_utc():
    naive = datetime.datetime(2026, 11, 17, 8, 30, 0)
    assert format_iso(naive) == "2026-11-17T08:30:00Z"


def test_parse_iso_round_trip():
    moment = datetime.datetime(2026, 11, 17, 8, 30, 0, tzinfo=datetime.timezone.utc)
    assert format_iso(moment) == "2026-11-17T08:30:00Z"
    parsed = parse_iso("2026-11-17T08:30:00Z")
    assert parsed == moment
    assert parsed.strftime(ISO_FORMAT) == moment.strftime(ISO_FORMAT)


def test_parse_iso_rejects_invalid_input():
    assert parse_iso(None) is None
    assert parse_iso("") is None
    assert parse_iso("not-a-timestamp") is None


def test_timestamp_to_iso_handles_valid_and_invalid_values():
    assert timestamp_to_iso(0) == "1970-01-01T00:00:00Z"
    assert timestamp_to_iso(None) is None
    assert timestamp_to_iso("nonsense") is None
