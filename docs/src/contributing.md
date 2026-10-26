# Contributing to AstroNbodySim.jl

Thanks for your interest in making this package better! Bug reports,
feature requests, documentation fixes, and pull requests are all
welcome.

## Issues

Before opening a new issue please check whether the problem already
exists among:

1. These docs (`docs/src/…`) and any inline error messages.
2. The source-level docstrings (use `help?>` in the Julia REPL).
3. The change log / release notes.
4. The existing open and closed issue tracker on GitHub.

When filing a **bug**, please provide:

- A *minimal, fully reproducible* example (MWE).
- Your Julia version (`versioninfo()` in the REPL) and the
  AstroNbodySim.jl commit / version (`]status`).
- The full error output (stacktrace, including any `@info` or
  `@warn` lines emitted by the package) — copy-paste as text, not a
  screenshot.
- If the bug is GPU-related, the `nvidia-smi` driver version and the
  CUDA toolkit version.

For **feature requests**, describe:

- The use case (a small program snippet that *would* use the feature
  if it existed).
- The desired API surface (function signature, struct fields).
- Any links to papers, prior implementations, or related projects.

## Pull requests

Every PR should ship with:

- A short, imperative summary line
  (`"Fix syncflag inversion in check_sync_point_local"`).
- A test (`test/`) demonstrating the bug, or showing the new behaviour.
- A documentation change:
  - New public symbols → `@docs` block + a `tutorial/…` or `manual/…`
    page if the change is user-visible.
  - Internal refactors → note in the source comment, no doc change.
- For complex algorithmic changes, a brief comment near the
  implementation explaining the method (and a paper reference if
  applicable).

Before pushing, please run the local test suite to confirm nothing
regressed:

```julia
]test
```

or, from the project root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Coding conventions

- Follow [Julia's style guide](https://docs.julialang.org/en/v1/manual/style-guide/).
- Prefer multiple `using` statements over giant imports.
- Use `Base.@kwdef`-style keyword structs for configuration objects.
- Document every public function with a `# Brief title`-style
  docstring and, where helpful, a usage example in ``jldoctest``
  format.
- Compatibility: aim for the oldest Julia release still listed in
  `Project.toml`'s `[compat]` — currently `1.12`.

## Release process

1. Update `Project.toml` `version`.
2. Update `CHANGELOG.md`.
3. Tag `vMAJOR.MINOR.PATCH` and push. The CI workflow creates the
   Zenodo record automatically.
