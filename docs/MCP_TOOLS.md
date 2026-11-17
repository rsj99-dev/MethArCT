# MethArCT MCP tools

`metharct-mcp` is an **optional** MCP server (local `stdio` transport, no
network port). It exposes the analysis pipeline to agent hosts; the host owns
model inference, MethArCT only provides tools.

```bash
pip install ".[mcp]"                                     # requires Python >= 3.10
metharct-mcp --workspace ./metharct_workspace            # start (stdio)
metharct-mcp --list-tools                                # print the manifest
METHARCT_ALLOWED_ROOTS="D:/genomes:E:/projects" metharct-mcp
```

Supported SDK range is `mcp>=1.2`. Both generations are accepted: 1.x
(`mcp.server.fastmcp.FastMCP`) and 2.x
(`mcp.server.mcpserver.MCPServer`, where FastMCP was renamed). On Python 3.8/3.9
the `mcp` extra is skipped by a version marker; `metharct doctor` reports why
and `metharct-mcp` exits with `DEPENDENCY_MISSING`.

Every tool delegates to `MethArCTClient`; none touches analyzer internals, and
none accepts a command, shell string or Python source.

## Envelope

All tools return one JSON document:

```json
{
  "ok": true,
  "request_id": "mcp_9f1c...",
  "data": {},
  "warnings": [],
  "error": null
}
```

```json
{
  "ok": false,
  "request_id": "mcp_9f1c...",
  "data": null,
  "warnings": [],
  "error": {
    "code": "INPUT_FILE_NOT_FOUND",
    "message": "Input file does not exist.",
    "retryable": false,
    "details": {"path": "..."}
  }
}
```

`message` stays short and actionable; full tracebacks go to
`jobs/<job_id>/logs/stderr.log`. Host absolute paths are never returned —
artifacts are addressed by `artifact_id` and relative path.

### Argument handling

* Over the MCP transport the SDK validates each call against the tool's input
  schema, so an undeclared field (for example `command`) is dropped before
  dispatch and can never reach MethArCT. This is verified end to end by
  `tests/ai/test_mcp_integration.py`.
* The Python tool layer (`metharct.ai.mcp_tools.call_tool`) additionally
  rejects undeclared fields with `INVALID_ARGUMENT` for in-process callers
  (`tests/ai/test_mcp_tools.py`).
* Wrong types, missing required fields and out-of-range values are rejected in
  both paths.

---

## Tool 1 — `metharct_capabilities`

Discover the installed version, modules, options and real dependency status.
Call this first: only modules with `enabled: true` can run on this machine.

```json
{"include_details": true}
```

Returns the capability manifest (see `metharct-capabilities.schema.json`).

## Tool 2 — `metharct_doctor`

Check the local environment. Never installs anything.

```json
{"check_databases": true}
```

Returns `ok`, `blocking_issues`, per-check `status` and `modules[].enabled`.

## Tool 3 — `metharct_inspect_input`

Validate a FASTA file before spending compute on it.

```json
{"path": "MZD129.faa", "expected_type": "auto"}
```

Returns sequence type + rule-level confidence, counts, length statistics,
duplicate ids, invalid characters, SHA-256 and compatible modules.

## Tool 4 — `metharct_run_analysis`

Submit a job. Returns immediately unless `wait` is true.

```json
{
  "input_path": "MZD129.faa",
  "input_type": "protein_fasta",
  "modules": ["pathway", "temperature", "salinity", "ph"],
  "options": {"pathway": {"threads": 8}},
  "wait": false,
  "timeout_seconds": 3600,
  "idempotency_key": "MZD129-v1"
}
```

Returns `job_id`, `status`, `requested_modules`, `unavailable_modules` and the
`status_tool` / `result_tool` hints. It never reports a success that did not
happen.

## Tool 5 — `metharct_get_job_status`

```json
{"job_id": "job_20261117T074626_d3f1ebc6"}
```

Returns the real status, phase, progress counters, completed / failed / skipped
modules and the error summary.

## Tool 6 — `metharct_get_result_summary`

```json
{"job_id": "...", "modules": ["pathway", "temperature"]}
```

Returns the compact report: input block, one envelope per module, warnings,
errors, artifact list and provenance. It never returns whole result tables.

## Tool 7 — `metharct_read_artifact`

Omit `artifact_id` to list artifacts:

```json
{"job_id": "..."}
```

Read one bounded page:

```json
{
  "job_id": "...",
  "artifact_id": "a001",
  "offset": 0,
  "limit": 50,
  "filter": {"column": "Detection_Status", "value": "Detected", "contains": false}
}
```

Returns `columns`, `rows`, `total_rows`, `matched_rows`, `returned_rows`,
`next_offset` and `truncated`. Only manifest-registered CSV/TSV/JSON/text files
can be read; XLSX and other binary media types are refused with
`ARTIFACT_ACCESS_DENIED`.

## Tool 8 — `metharct_cancel_job`

```json
{"job_id": "..."}
```

Returns `cancellation_supported`, `accepted`, `worker_terminated` and the
resulting job status. A worker is stopped between modules: cancellation is
explicitly best effort and never claims to have interrupted a DIAMOND process
mid-run.

---

## Resources

| URI | Content |
|---|---|
| `metharct://capabilities` | live capability manifest |
| `metharct://docs/analysis-modules` | module descriptions and scientific limitations |
| `metharct://docs/result-schema` | the packaged report JSON schema |
| `metharct://docs/limitations` | interpretation limits that must accompany results |

Resources expose official, verified information only — never arbitrary files
from the host.

---

## Host configuration example

```json
{
  "mcpServers": {
    "metharct": {
      "command": "metharct-mcp",
      "args": ["--workspace", "D:/metharct_workspace"],
      "env": {
        "METHARCT_ALLOWED_ROOTS": "D:/genomes;D:/metharct_workspace"
      }
    }
  }
}
```

On Windows use `;` between allowed roots; on Linux/macOS use `:`.

## Recommended agent flow

1. `metharct_capabilities` → confirm which modules are available.
2. `metharct_doctor` (optional) → confirm DIAMOND and databases.
3. `metharct_inspect_input` → confirm the file is usable.
4. `metharct_run_analysis` (`wait: false`) → get `job_id`.
5. `metharct_get_job_status` → poll until terminal.
6. `metharct_get_result_summary` → read predictions and limitations.
7. `metharct_read_artifact` → page through a table only if needed.

Never state a biological conclusion from `not_detected` alone, and never
present a prediction as a measurement.
