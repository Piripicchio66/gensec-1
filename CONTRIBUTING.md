# Contributing to GenSec

Thank you for considering a contribution to GenSec!  This document explains
the process and requirements.

---

## Contributor License Agreement (CLA)

**All contributors must sign the CLA before any pull request can be merged.**

When you open your first PR, the [CLA Assistant](https://cla-assistant.io/)
bot will automatically comment with a link to review and sign the
[CLA](CLA.md).  This is a one-time step — once signed, it covers all your
future contributions.

The CLA exists because GenSec uses a **dual-licensing model**:

- The public repository is licensed under **AGPL-3.0**.
- A separate **commercial license** is available for organizations that cannot
  comply with the AGPL.

The CLA ensures the Author can continue offering both options.  Your
Contribution remains credited to you in the `NOTICE` file and in the git
history.

---

## What You Can Contribute

- Bug fixes and regression tests
- New material models (with validation against published data)
- New section geometries or primitives
- Performance improvements
- Documentation improvements, translations, and examples
- Test cases (especially analytical benchmarks)

---

## Coding Standards

GenSec follows strict documentation and style conventions:

1. **Language**: all code, comments, docstrings, and commit messages in
   **English**.

2. **Docstrings**: [numpydoc](https://numpydoc.readthedocs.io/) style,
   Sphinx-compatible.  Every public function, class, and method must have a
   complete docstring.

3. **Math**: any formula referenced in docstrings or documentation must be
   written in LaTeX, using Sphinx-compatible ``:math:`` roles or ``.. math::``
   directives.

4. **Type hints**: all function signatures must include type annotations.

5. **Tests**: every new feature or bug fix must include corresponding tests.
   Analytical/benchmark test cases are strongly preferred over purely
   numerical regression tests.

6. **Commits**: use [Conventional Commits](https://www.conventionalcommits.org/)
   format:
   ```
   feat(integrator): add Gauss-Lobatto quadrature option
   fix(concrete): correct ultimate strain for high-strength classes
   docs(readme): add biaxial bending example
   ```

---

## Pull Request Process

1. **Fork** the repository and create a feature branch from `main`.
2. Ensure all existing tests pass: `pytest`.
3. Add tests for your changes.
4. Update documentation if needed.
5. Open a PR with a clear description of what and why.
6. Sign the CLA when prompted by the bot.
7. Address review comments.

---

## Attribution

All accepted contributions are recorded in the [`NOTICE`](NOTICE) file.  When
your PR is merged, the maintainer will add your entry.  The format is:

```
Full Name <optional-url-or-email>
  — Brief description of contribution.
  — Date or version.
```

---

## Code of Conduct

Be respectful, constructive, and professional.  Technical disagreements are
welcome; personal attacks are not.
