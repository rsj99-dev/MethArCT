# MethArCT AI layer — security and execution boundaries

## 1. What an agent can reach

The AI layer exposes four things only:

1. read-only environment inspection (`capabilities`, `doctor`);
2. reading of a **user-supplied input path**;
3. analysis jobs that call the existing analyzers;
4. reading of files that a job itself produced (artifact manifest).

An agent can never:

* execute a shell command, Python snippet or arbitrary program;
* smuggle a `command`, `shell` or `python_code` field into an analysis — over
  the MCP transport the SDK's schema validation drops undeclared fields before
  dispatch, and the Python tool layer rejects them with `INVALID_ARGUMENT`;
* read an arbitrary host path through the artifact tools;
* change database paths, thresholds or model files at runtime beyond the
  documented module options.

## 2. Path handling

| Rule | Implementation |
|---|---|
| All user paths are made absolute before comparison | `paths.normalize` |
| Containment checks are case-insensitive on Windows | `paths.is_within` (`os.path.normcase` + `realpath`) |
| Optional read whitelist | `allowed_roots` (constructor) or `METHARCT_ALLOWED_ROOTS` (`;` on Windows, `:` on POSIX) |
| `..` escapes rejected | `paths.safe_join` |
| Artifacts addressable only by manifest id | `artifacts.ArtifactRegistry.resolve` |
| Absolute / `..` artifact paths rejected | `is_within(candidate, job_dir)` on every resolve |
| Symlink escapes rejected | `realpath` is used by the containment check |
| Binary files not sent as text | media-type allowlist (`text/csv`, `text/tab-separated-values`, `application/json`, `text/plain`, `text/x-fasta`) |
| Large files paged | `offset` / `limit` (`<= 1000`), `MAX_ARTIFACT_BYTES` = 2 MiB |

`allowed_roots` defaults to *unset* for the plain Python/CLI use (the caller
already has filesystem access) and is **set explicitly** for the MCP server.

## 3. Process execution

* The worker is started with an argument array:
  `[sys.executable, "-m", "metharct.ai.worker", "--workspace", ..., "--job-id", ...]`
  and `shell=False`.
* Cancellation uses `signal.SIGTERM` (POSIX) or
  `taskkill /PID <pid> /T /F` as a fixed argument array with `shell=False`
  (Windows).
* Every external tool call in `metharct.core` uses `subprocess.run(..., args)
  with `shell=False`; the AI layer adds none.
* `eval`, `exec` and `shell=True` do not appear anywhere in `metharct/ai/`.

Verify with:

```bash
python -m pytest tests/ai -q
grep -rn "shell=True\|eval(\|exec(" metharct/ai || echo "clean"
```

## 4. Resource limits

| Limit | Default | Where |
|---|---|---|
| Input file size | 500 MiB | `max_file_size_mb` / `DEFAULT_MAX_FILE_SIZE_MB` |
| Concurrent jobs | 1 | `max_parallel_jobs` |
| Job timeout | none | `timeout_seconds` (enforced between modules + a hard parent timeout) |
| Artifact page size | 50 rows | `--limit` / `MCP limit` (`<= 1000`) |
| Artifact response bytes | 2 MiB | `MAX_ARTIFACT_BYTES` |
| Logs | per job | `logs/stdout.log`, `logs/stderr.log` |

A crashed or killed worker cannot leave a "success" behind: the parent marks the
job `failed` with `phase: worker_lost` after three liveness checks.

## 5. Data written to disk

* Everything a job produces lives under `<workspace>/jobs/<job_id>/`.
* No file is written outside the workspace and the job directory.
* Atomic writes (temp file + `os.replace`, with retries on Windows sharing
  violations) prevent half-written `job.json` / report files.
* `doctor` never writes configuration; it only probes a temporary file in the
  output / temp directories to check writability.

## 6. What is never logged or returned

* No credentials, tokens or environment dumps.
* Host absolute paths are stripped from job views, module results, artifact
  pages and MCP responses; only `input_filename` and workspace-relative paths
  are returned.
* Stack traces stay in `logs/stderr.log`; the agent sees a stable error code
  and a short message.

### Deliberate exception: `doctor`

`metharct_doctor` reports the *configuration* it found — for example the path
that `tools.diamond.path` resolves to, or the CheckM2 database directory. That
is the purpose of an environment self-check and it is enabled by design. A host
that must not disclose local paths should expose `metharct_capabilities`
(which returns module status without free-form paths in `reasons`) instead of
`metharct_doctor`, or run the MCP server in a container whose paths are
non-sensitive.

## 7. Reporting a problem

Because error codes are part of the contract, an unexpected `INTERNAL_ERROR`
should be reported together with:

* the job id,
* `logs/stdout.log` and `logs/stderr.log`,
* the output of `metharct doctor --format json`,
* the package version (`metharct --version`).
