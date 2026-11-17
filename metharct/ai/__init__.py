# -*- coding: utf-8 -*-
"""MethArCT AI-ready adaptation layer.

The layer is intentionally *additive*: it exposes a stable Python API, a
machine-readable CLI mode, a job manager and an optional MCP server without
changing any scientific algorithm, threshold, model file or legacy output.

Typical use::

    from metharct.ai import MethArCTClient

    client = MethArCTClient(workspace="./metharct_workspace")
    job = client.run_analysis("MZD129.faa", modules=["temperature", "ph"])
    print(job["status"], job["job_id"])

Heavy submodules are imported lazily so ``import metharct.ai`` stays cheap.
"""

from .constants import AI_API_VERSION, SCHEMA_VERSION
from .errors import ALL_ERROR_CODES, ErrorCode, MethArCTError

_LAZY_IMPORTS = {
    "MethArCTClient": ".api",
    "run_doctor": ".doctor",
    "get_capabilities": ".capabilities",
    "MODULE_REGISTRY": ".capabilities",
    "MODULE_NAMES": ".capabilities",
    "inspect_fasta": ".validators",
    "JobManager": ".jobs",
    "ArtifactRegistry": ".artifacts",
    "run_job": ".worker",
}

__all__ = [
    "AI_API_VERSION",
    "SCHEMA_VERSION",
    "ALL_ERROR_CODES",
    "ErrorCode",
    "MethArCTError",
    "MethArCTClient",
    "get_capabilities",
    "run_doctor",
    "MODULE_REGISTRY",
    "MODULE_NAMES",
    "inspect_fasta",
    "JobManager",
    "ArtifactRegistry",
    "run_job",
]


def __getattr__(name):
    """Resolve the public names on first access (PEP 562)."""
    module_path = _LAZY_IMPORTS.get(name)
    if module_path is None:
        raise AttributeError("module {0!r} has no attribute {1!r}".format(__name__, name))
    import importlib

    module = importlib.import_module(module_path, __package__)
    return getattr(module, name)
