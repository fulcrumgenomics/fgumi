# Developing fgumi

## Getting Started

If you already have the Rust toolchain installed, make sure you have at least version `1.87.0` (e.g., with `rustup show` or `rustc --version`).
If not, you can run `rustup update` to install the latest version.

If you do not have the Rust toolchain installed, the best way to install it is with [`rustup`](https://rustup.rs/) using the following command (which assumes that you have [`curl`](https://everything.curl.dev/get) installed):

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

It is also recommended to install the following additional tools.
These tools are used in the CI/CD workflows, and there are [aliases](#code-conventions) defined that enable you to use them in the same way locally to simulate the checks that are performed during CI.

* [`cargo-nextest`](https://nexte.st/): A replacement for the Cargo test runner. It is recommended to use this for running tests.
  * Installation: Download a [pre-built binary](https://nexte.st/book/pre-built-binaries.html) or run `cargo install cargo-nextest --locked`.
  * Usage: `cargo nextest run`, or use the alias to emulate how it is run during CI: `cargo ci-test`.
* [`cargo-audit`](https://github.com/RustSec/rustsec/tree/main/cargo-audit) to perform a security audit of dependencies.
  * Installation: `cargo install cargo-audit`
  * Usage: `cargo audit`
## Project structure

This project uses the [conventional Cargo package layout](https://doc.rust-lang.org/cargo/guide/project-layout.html).

* `.cargo/config.toml`: Cargo configuration. This is primarily used to define [aliases](#code-conventions) for commonly used Cargo commands.
* `.github/ISSUE_TEMPLATE`: [Issue form templates](https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests/configuring-issue-templates-for-your-repository#creating-issue-forms) for creating new GitHub issues.
* `.github/workflows`: Workflows for continuous integration (CI) and deployment (CD).
* `.gitignore`: Files to ignore when committing changes with `git`, i.e., files that should not be version controlled.
* `benches/`: Optional directory with benchmarks using [`criterion`](https://github.com/bheisler/criterion.rs).
* `build.rs`: Build script that is executed when building the project.
* `Cargo.toml`: The project manifest. This contains all the necessary information for building and publishing the project crate(s).
* `cliff.toml`: Configuration for [`git-cliff`](https://github.com/orhun/git-cliff), which automatically generates the changelog based on commit messages.
* `clippy.toml`: Configuration for Clippy, the Rust linter that is used to check the project's code against best practices.
* `docs/`: All project documentation.
* `LICENSE`: The project license.
* `README.md`: The project description. This is displayed on GitHub, on the [crates.io](https://crates.io/) page for the project, and on the [docs.rs](https://docs.rs/) documentation for the project.
* `release-plz.toml`: Configuration for [release-plz](https://release-plz.ieni.dev/), which is used to automate releases.
* `rustfmt.toml`: Configuration for [rustfmt](https://github.com/rust-lang/rustfmt), which is used to automatically format the project's code, and to check code formatting in the CI workflow.
* `src/`: All Rust source code, including unit tests.

## Git etiquette

The use of source control - specifically, `git` via GitHub - is central to the development of this project.
As such, there are some guidelines for minimizing friction, especially when multiple developers are working on the project at the same time.

* Do not push to [someone else’s branch](#git-branches) without first getting their approval; this includes committing suggestions in a PR review. It is always acceptable to branch off and make a PR into someone else's branch.
* The person who opens and owns a PR is the only person who may merge the PR (once all checks are passed); anyone else should ask the PR owner for permission before merging the PR.
* Use `rebase` rather than `merge`, to maintain a clean, linear history. Learning to use [`--update-refs`](https://andrewlock.net/working-with-stacked-branches-in-git-is-easier-with-update-refs/) can make this much easier.
* If you’re moving/renaming a source file, use `git mv` to make sure git tracks it properly.
* Be careful and purposeful about what you commit to the repository. Do not commit IDE-specific files, or your own personal test files and notes. Use `.gitignore` to avoid commiting local-only files (although please use `~/.gitignore` for any files specific to your personal development process). Any test data larger than a few kb should, not be committed to the repository; instead, generate it dynamically or pull it from an object store (e.g., AWS S3) during setup of your integration tests.

## Development lifecycle

### Issues

Generally, we follow the practice of [issue-driven development](https://simonwillison.net/2022/Nov/26/productivity/): any development should begin as a GitHub issue, whether from a project developer (such as yourself) or a user.
This project has [issue templates](../.github/ISSUE_TEMPLATE/) to help ensure that each issue has the minimum necessary information to fix the bug, add the feature, or otherwise complete the task.

If you are one of the [project admins](../.github/CODEOWNERS), please add/change issue labels, projects, and milestones as necessary.
If you plan to work on the issue shortly, please assign it to yourself if you have appropriate permissions, otherwise add a comment expressing your interest.

### Git branches

The `main` branch is protected and can only be pushed to from a pull request (PR).
This means you'll need to do your development in a separate git branch and open a PR to `main` when you're done.

As with all aspects of development, you should follow any Fulcrum best practices that apply.
In this case, Fulcrum [git etiquette](#git-etiquette) dictates breaking up larger issues into multiple smaller PRs and maintaining a clean commit history through the use of rebasing or squash-merging.
Each PR should use a separate branch, and branches are automatically deleted after they are merged.
Branches can be independent or cascading, but keep in mind that cascading branches can become unwieldy when rebasing is required (although [`git rebase --update-refs`](https://andrewlock.net/working-with-stacked-branches-in-git-is-easier-with-update-refs/) can make this much easier).

Branches should follow the standard naming convention `<issue-number>/<user>/<type>-<description>`:

* `issue-number` is the GitHub issue number.
* `user` is your GitHub username.
* `type` is the type of change being made in the branch. This is the same type that you will use in the [commit message](#commit-messages) when squash-merging the PR.
* `description` is a short, hyphenated description of the issue.

For example, `42/jdidion/fix-fibonacci-calculation` is a good branch name for a change (by user `jdidion`) that fixes a bug (reported in issue #42) in a function that calculates Fibonacci numbers.

### Code conventions

When writing code, you should follow the Fulcrum best practices, which include writing thorough documentation comments and using Cargo to format (using `rustfmt`) and lint (using `clippy`) your code.

It is recommended to use an IDE that has support for Rust formatting and linting, or to configure your editor to integrate with Cargo and/or [`rust-analyzer`](https://rust-analyzer.github.io/).
For example, [VS Code](https://code.visualstudio.com/) has excellent support for Rust via the [rust-analyzer plugin](https://code.visualstudio.com/docs/languages/rust).

You can use the following Cargo aliases (defined in the project [Cargo configuration](../.cargo/config.toml)) to run the same checks that are run in the CI workflow before you commit your code:

* `ci-test`: Execute the same command used in CI to run tests.
* `ci-fmt`: Execute the same command used in CI to check code formatting.
* `ci-lint`: Execute the same command used in CI to check code style.

For example, the following command will run all CI checks:

```
cargo ci-test && cargo ci-fmt && cargo ci-lint
```

When a particular topic is not covered by Fulcrum best practices, you can [start a Discussion](https://github.com/fulcrumgenomics/fgumi/discussions) or search online (e.g., [Idiomatic Rust](https://github.com/mre/idiomatic-rust)).

#### Documentation

Please read the [rustdoc book](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html) to learn how to write good documentation comments.
At a minimum, all of the functions and data structures that are a part of the project's public API should be documented completely.

You are encouraged to include [examples](https://doc.rust-lang.org/rust-by-example/testing/doc_testing.html) in the function documentation, as these are also run as test cases.


#### Testing

In general, every addition of (non-trivial) code to the project must be accompanied by unit tests that, at a minimum, exercise all the "happy paths" (i.e., non-error cases) through the code.
It is strongly recommended to also test at least the most commonly expected error cases to make sure they are handled correctly.
Most importantly, every bugfix *must* be accompanied by at least one regression test.

Please read the [documentation](https://doc.rust-lang.org/rust-by-example/testing/unit_testing.html) on unit testing in Rust if you haven't already done so.
We use [`rstest`](https://docs.rs/rstest/latest/rstest/) for writing tests, and we use [`nextest`](https://nexte.st/) for running tests in CI and recommend you use it locally.
The Cargo `ci-test` alias will run the project tests in the same way they are run in CI.


### Commit messages

We follow the practice of using [conventional commit messages](https://www.conventionalcommits.org/en/v1.0.0/), which enables automation of versioning and the changelog.
Conventional commit messages are only required for commits that are merged into the `main` branch (e.g., the squashed commit message when merging a PR), but you are encouraged to always use them as it can help you when you need to write the final commit message, as well as others working in your branch or reviewing your PR.

A commit message is structured as follows:

```
<type>[(scope)][!]: <description>

[optional body]

[optional footer(s)]
```

The commit `type` indicates the impact of the change on users of the package.
When a new version of the package is released, the `type`s of the commits included in the release are used to determine how the version is incremented, according to the rules of [semantic versioning](https://semver.org).

* `fix`: patches a bug in your codebase. This is generally associated with issues labeled `bug`. When there are only `fix` commits included in a release, it causes an increment in the `PATCH` version.
* `feat`: introduces a new feature to the codebase. This is generally associated with issues labeled `enhancement`. When there are only `feat` and `fix` commits included in a release, it causes an increment in the `MINOR` version.
* Commits that only change dependencies (e.g., update versions) should have type `build`.
* Types other than `fix` and `feat` are allowed, but typically have no impact on the version. Currently, the following other types are allowed. If you identify a need for a new type, please add it to this list:
  * `build`: changes to how the package is compiled. Typically, such changes involve `build.rs` or `Cargo.toml`.
  * `chore`: catch-all for maintainance activities that don't fall under any of the other types.
  * `ci`: changes to the CI/CD workflows or other files in the `.github` folder.
  * `config`: changes to project configuration that does not affect compilation. Typically, these changes involve `*.toml` files.
  * `docs`: changes to code comments or other documentation (`README.md` or files in `docs/`).
  * `example`: addition of or changes to example code in the documentation comments.
  * `perf`: changes that only affects the performance of existing code. This should be accompanied by regression tests.
  * `refactor`: changes to the structure of the code that does not change any functionality. This should be accompanied by regression tests.
  * `style`: changes to code formatting only. Note that there should not be any `style` commits to `main` (since only properly styled code should pass the CI checks) unless accompanied by a `config` change to the formatting rules.
  * `test`: addition of or changes to test code only.

The commit's `scope` (if any) provides additional context, such as the component of the project that is being changed.
For example, `feat(parser)` implies the addition of features to the project's parser component.

Breaking changes (i.e., changes that break backward compatibility) must be annotated with an exclamation point (`!`) following the type/scope, and must also be accompanied by a footer beginning with `BREAKING CHANGE:` that describes what is broken.
A breaking change causes an increment in the `MAJOR` version of the package.

Additional footers may be added that follow the [git trailer](https://git-scm.com/docs/git-interpret-trailers) format and conventions.

### Pull requests

When you are ready to share your changes with others, please open a pull request using the provided [template](../.github/PULL_REQUEST_TEMPLATE.md).
You are encouraged to create a draft pull request early on (e.g. as soon as you push your first commit) to solicit feedback.
Please provide all the required information and any of the optional information that applies.
If you omit any required information, please provide justification.
Please assign yourself to the PR.

When you prepare a PR, your primary consideration should be "how am I going to make this easy for the reviewer to review?".
Ask yourself:

* How do I organize commits to see more complicated changes?
* How do I enable multiple rounds of review?
* How do I tell a story of the work I am doing?

A PR should be focused on one thing.
Prefer to break up larger multi-logical-commit changes into multiple cascading PRs, e.g., `feature_b => feature_a => main`.
Never bundle unrelated changes together in a PR - especially broad formatting changes or large scale movements of code.

#### Continuous integration

Several checks will be run automatically when you open a PR and each time you push changes to a branch associated with an open PR:

* Run all unit tests
* Check code formatting
* Check code style


The [development lifecycle](#development-lifecycle) section describes how you can run these same checks locally; you should make sure to do so before opening/updating a PR.
If any of these checks fail, you must resolve them before requesting a review.
Please do not globally disable any checks in your branch, as such changes will be rejected.
If you think that any of these checks should be changed, please open a separate PR (type `config` or `ci`) and resolve that first.

#### Review process

When you have completed development and are ready to have your changes reviewed, mark your PR as "Ready for review".
If you have the appropriate permissions, you can request reviewers and add any labels, projects, and milestones as necessary.
Otherwise, add a comment to the corresponding issue with a link to your PR and state that it is ready for review.

While a PR is under review:

* Do not make any further pushes to the associated branch.
* Do not close the PR.
* Do not open another PR from the same branch.

If you believe there are extenuating circumstances, please communicate with and get approval from all active and potential reviewers before taking any of the above actions.

When a round of review has been completed, you should address all the issues raised by the reviewers - either make the suggested change or propose an alternative change and discuss until there is agreement.
After each conversation is resolved (i.e., the change has been pushed or it has been agreed that no change is necessary), click the associated button.
When all conversations have been resolved, click the button(s) by the name of each reviewer who provided a review on the previous round to request a re-review.

Once all reviewers sign off, all checks are successful, and the branch has been rebased to the target branch, click the "Squash and merge" button.
You must provide a [commit message with the correct format](#commit-messages).
Make sure that the commit message type, optional scope, and optional breaking change signifier (`!`) appear in the commit message title.
Once the PR is merged, the branch will be deleted automatically.
You should rebase any branches that were pointing to the merged branch to now point to `main`.

## Continuous deployment

The release process is largely managed by automated workflows.
These workflows require some [initial configuration](https://github.com/fulcrumgenomics/rust-template) that must be completed before the first PR is opened.
This process was inspired by [this excellent blog post](https://blog.orhun.dev/automated-rust-releases/), which can be referenced for more details.

Each time a pull request is merged into `main`, the [publish](.github/workflows/release.yml) workflow uses [`release-plz`](https://release-plz.ieni.dev/docs/github) to open or update a release PR.
The release PR tracks all changes that have been made since the last release, including an automatically updated changelog and version.
The use of [properly formatted commit messages](#commit-messages) is essential to ensuring the changelog and version are updated correctly.

### Making a release

Making a release simply requires merging the current release PR.
When the `publish` workflow runs, it recognizes that the PR being merged is a release PR, and it does the following:

* Creates a new tag with the release version.
* Publishes any updated artifacts to [`crates.io`](https://crates.io/).


The creation of the version tag triggers the [`release`](.github/workflows/release.yml) workflow, which does the following:

* Creates a new [GitHub release](https://github.com/fulcrumgenomics/fgumi/releases) with changes from the (automatically updated) changelog.
* Builds binaries for the most commonly used platforms and attaches them to the GitHub release.

