# Contributing to exposomeSC

Bug reports and pull requests are welcome at
<https://github.com/CuiweiG/exposomeSC/issues>. Please open pull requests
against `main`.

The rest of this page records how changes move between GitHub and
Bioconductor. The two sides use different names for the development
branch, which is the easiest place to make a mistake.

## Which branch a change goes to

exposomeSC is in Bioconductor devel (3.25) and has no release branch yet.

| Change | GitHub branch | Bioconductor branch | Version |
|---|---|---|---|
| Any change before the Bioconductor 3.25 release | `main` | `devel` | 0.99.z |
| After the 3.25 release: new features, changed defaults, refactoring | `main` | `devel` | 1.1.z (odd y) |
| After the 3.25 release: bug fixes and documentation corrections | `RELEASE_3_25`, after `main` | `RELEASE_3_25`, after `devel` | 1.0.z (even y) |

- GitHub `main` is Bioconductor `devel`. Bioconductor only builds
  `devel` and `RELEASE_*` branches, so every push to Bioconductor names
  the target explicitly: `main:devel`.
- The review is over: changes no longer go through
  `BiocStaging/exposomeSC` or the `bioc-revision-2026-11` branch, but
  straight to Bioconductor's own repository.
- Once `RELEASE_3_25` exists, a bug fix goes into `main` first and is
  then cherry-picked to `RELEASE_3_25`. Each branch gets its own version
  bump.

## Version numbers

Bioconductor uses x.y.z version numbers
(<https://contributions.bioconductor.org/versionnum.html>):

- Until the 3.25 release the version stays at 0.99.z. At the release the
  Bioconductor team sets the release branch to 1.0.0 and `devel` to
  1.1.0; after that, y is even on the release branch and odd on `devel`.
  Do not change x or y yourself.
- Bump z in every commit that should reach users (0.99.4 to 0.99.5). A
  pushed commit without a version bump is built, but not propagated:
  `BiocManager::install()` keeps serving the previous version.
- Commits that only touch files left out of the package by
  `.Rbuildignore` (for example `.github/`) do not need a bump.
- Builds run once a day and take about 24 hours, so a push appears on
  the build report (<https://bioconductor.org/checkResults/>) one or two
  days later. It reaches users once a build with the new version
  finishes without errors.

## Remotes

This is the setup Bioconductor recommends
(<https://contributions.bioconductor.org/git-version-control.html>):

```sh
git remote add upstream git@git.bioconductor.org:packages/exposomeSC.git
git remote -v
# origin    https://github.com/CuiweiG/exposomeSC.git (fetch)
# origin    https://github.com/CuiweiG/exposomeSC.git (push)
# upstream  git@git.bioconductor.org:packages/exposomeSC.git (fetch)
# upstream  git@git.bioconductor.org:packages/exposomeSC.git (push)
```

`ssh -T git@git.bioconductor.org info` lists what the registered SSH key
may read (R) and write (W); `packages/exposomeSC` should show `R W`.

## Pushing

Before committing, bring `main` up to date with Bioconductor, which can
commit to `devel` too (for example the version bumps at each release):

```sh
git checkout main
git fetch upstream
git merge --ff-only upstream/devel
```

Then push the same commit to both sides:

```sh
git push upstream main:devel
git push origin main
```

To push exactly the commit that was checked, push by hash:

```sh
git push upstream <hash>:refs/heads/devel
git push origin <hash>:refs/heads/main
```

Never force-push: Bioconductor does not accept rewritten history.

## Continuous integration

The workflows in `.github/workflows/` run on `main` with R-devel, which
Bioconductor devel uses: `R CMD check` on Linux, macOS and Windows,
`BiocCheckGitClone()` and `BiocCheck()`, test coverage, and a check of
`citation()`, links and spelling. The pkgdown site is built from `main`.
When Bioconductor 3.25 is released and `RELEASE_3_25` appears, check the
R versions in the workflow files.
