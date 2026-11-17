# -*- coding: utf-8 -*-
"""UTC timestamp helpers.

``datetime.utcnow()`` and ``datetime.utcfromtimestamp()`` are deprecated from
Python 3.12, so every timestamp in the AI layer is produced here with an
explicit UTC timezone and serialised as ISO-8601 with a ``Z`` suffix.
"""

import datetime
from typing import Optional

#: Format used for every timestamp the AI layer writes.
ISO_FORMAT = "%Y-%m-%dT%H:%M:%S"


def utc_now() -> datetime.datetime:
    """Return the current time as a timezone-aware UTC datetime."""
    return datetime.datetime.now(datetime.timezone.utc)


def format_iso(moment: datetime.datetime) -> str:
    """Serialise *moment* as ``YYYY-MM-DDTHH:MM:SSZ`` (UTC, second precision)."""
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=datetime.timezone.utc)
    return moment.astimezone(datetime.timezone.utc).strftime(ISO_FORMAT) + "Z"


def utc_now_iso() -> str:
    """Return the current time as an ISO-8601 UTC string with ``Z`` suffix."""
    return format_iso(utc_now())


def timestamp_to_iso(epoch_seconds: float) -> Optional[str]:
    """Convert a POSIX timestamp to an ISO-8601 UTC string."""
    try:
        moment = datetime.datetime.fromtimestamp(float(epoch_seconds), datetime.timezone.utc)
    except (TypeError, ValueError, OSError, OverflowError):  # pragma: no cover - defensive
        return None
    return format_iso(moment)


def parse_iso(text: Optional[str]) -> Optional[datetime.datetime]:
    """Parse a timestamp produced by :func:`format_iso` (``None`` when invalid)."""
    if not text:
        return None
    cleaned = str(text).strip()
    if cleaned.endswith("Z"):
        cleaned = cleaned[:-1]
    try:
        return datetime.datetime.strptime(cleaned, ISO_FORMAT).replace(tzinfo=datetime.timezone.utc)
    except ValueError:
        return None
