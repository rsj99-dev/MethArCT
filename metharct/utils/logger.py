#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logging utilities for MethArCT
"""

import os
import sys
import logging
from pathlib import Path
from typing import Optional
from loguru import logger

def setup_logger(name: str = "metharct", 
                level: str = "INFO",
                log_file: Optional[str] = None,
                console: bool = True) -> logging.Logger:
    """
    Setup logger with both file and console output
    
    Args:
        name: Logger name
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Path to log file (optional)
        console: Whether to output to console
    
    Returns:
        Configured logger instance
    """
    # Remove default loguru handler
    logger.remove()
    
    # Add console handler if requested
    if console:
        logger.add(
            sys.stderr,
            level=level,
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
            colorize=True
        )
    
    # Add file handler if log file is specified
    if log_file:
        # Create log directory if it doesn't exist
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        logger.add(
            log_file,
            level=level,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
            rotation="10 MB",
            retention="30 days",
            compression="zip",
            encoding="utf-8"
        )
    
    # Create standard logging adapter for compatibility
    class LoguruHandler(logging.Handler):
        def emit(self, record):
            try:
                level = logger.level(record.levelname).name
            except ValueError:
                level = record.levelno
            
            frame, depth = logging.currentframe(), 2
            while frame.f_code.co_filename == logging.__file__:
                frame = frame.f_back
                depth += 1
            
            logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())
    
    # Setup standard logging to work with loguru
    std_logger = logging.getLogger(name)
    std_logger.handlers = [LoguruHandler()]
    std_logger.setLevel(getattr(logging, level.upper()))
    
    return std_logger

def get_logger(name: str = "metharct") -> logging.Logger:
    """
    Get existing logger or create a new one
    
    Args:
        name: Logger name
    
    Returns:
        Logger instance
    """
    return logging.getLogger(name)