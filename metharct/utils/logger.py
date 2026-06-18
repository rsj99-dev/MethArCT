#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logging utilities for MethArCT

Unified logging using loguru. Provides a thin compatibility layer
for code that still calls standard ``logging.getLogger(name)``.
"""

import sys
import logging
from pathlib import Path
from typing import Optional
from loguru import logger


def setup_logger(
    name: str = "metharct",
    level: str = "INFO",
    log_file: Optional[str] = None,
    console: bool = True,
) -> logging.Logger:
    """Configure the global loguru logger and return a stdlib-compatible logger.

    Parameters
    ----------
    name : str
        Logger name (used as the stdlib logger name).
    level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    log_file : str or None
        Path to a rotating log file. ``None`` disables file logging.
    console : bool
        Whether to emit coloured output to stderr.

    Returns
    -------
    logging.Logger
        A standard-library logger that delegates to loguru.
    """
    # Remove default loguru handler
    logger.remove()

    # Console handler
    if console:
        logger.add(
            sys.stderr,
            level=level,
            format=(
                "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                "<level>{level: <8}</level> | "
                "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
                "<level>{message}</level>"
            ),
            colorize=True,
        )

    # File handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        logger.add(
            log_file,
            level=level,
            format=(
                "{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | "
                "{name}:{function}:{line} - {message}"
            ),
            rotation="10 MB",
            retention="30 days",
            compression="zip",
            encoding="utf-8",
        )

    # Build a stdlib logger that forwards to loguru
    class _LoguruBridge(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            try:
                lvl = logger.level(record.levelname).name
            except ValueError:
                lvl = record.levelno
            frame, depth = logging.currentframe(), 2
            while frame.f_code.co_filename == logging.__file__:
                frame = frame.f_back
                depth += 1
            logger.opt(depth=depth, exception=record.exc_info).log(
                lvl, record.getMessage()
            )

    std_logger = logging.getLogger(name)
    std_logger.handlers = [_LoguruBridge()]
    std_logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    return std_logger


def get_logger(name: str = "metharct") -> logging.Logger:
    """Return an existing stdlib logger (delegates to loguru via the bridge).

    Parameters
    ----------
    name : str
        Logger name.

    Returns
    -------
    logging.Logger
    """
    return logging.getLogger(name)
