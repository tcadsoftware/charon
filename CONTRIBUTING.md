# Contributing to Charon

The general workflow for everyone interacting with Charon is the following:

## Contents

1. [Creating Issues](#creating-issues)
   1. [Markdown](#markdown)
   1. [Issue Templates](#issue-templates)
   1. [Related Issues](#related-issues)
   1. [Labels](#labels)
      1. [Blocked Issues](#blocked-issues)
1. [Working Issues](#working-issues)
   1. [Breaking Issues Down](#breaking-issues-down)
   1. [When Work Begins](#when-work-begins)
      1. [If Working in Multiple Repositories](#if-working-in-multiple-repositories)
   1. [As Work Continues](#as-work-continues)
      1. [Commit Messages](#commit-messages)
      1. [Doxygen](#doxygen)
   1. [When Work is Complete](#when-work-is-complete)
   1. [Closing Old Issues](#closing-old-issues)
1. [Merge Requests](#merge-requests)
   1. [Reviewers](#reviewers)
   1. [Work-in-Progress](#work-in-progress)
   1. [Related Merge Requests](#related-merge-requests)
   1. [Merging](#merging)



## Creating Issues

Create issues in GitLab for any work that needs to be done.  Issues should be
created in the project in which the work will be completed.  The possible
locations include:
*  [tcad-charon](https://cee-gitlab.sandia.gov/Charon/tcad-charon/issues):  the
   scripts and CMake that tie the whole Charon codebase together,
*  [src](https://cee-gitlab.sandia.gov/Charon/src/issues):  the source code,
*  [nightlyTests](https://cee-gitlab.sandia.gov/Charon/nightlyTests/issues):
   any relatively short tests that are UUR,
*  [nightlyTestsOUO](https://cee-gitlab.sandia.gov/Charon/nightlyTestsOUO/issues):
   any relatively short tests that are OUO,
*  [heavyTests](https://cee-gitlab.sandia.gov/Charon/heavyTests/issues):  any
   relatively long tests that are UUR,
*  [heavyTestsOUO](https://cee-gitlab.sandia.gov/Charon/heavyTestsOUO/issues):
   any relatively long tests that are OUO, and
*  [docs](https://cee-gitlab.sandia.gov/Charon/docs/issues):  our Doxygen
   configuration and user manual.

Any issues that pertain to the whole of Charon (e.g., ~"Question"s,
~"Discussion"s, etc.) belong in
[tcad-charon](https://cee-gitlab.sandia.gov/Charon/tcad-charon/issues).  If
there is a task that will involve work in multiple repositories, the
recommendation would be to create a ~"Multi-Step Issue" in
[tcad-charon](https://cee-gitlab.sandia.gov/Charon/tcad-charon/issues) and then
create sub-issues in the appropriate projects, linking them back to the
~"Multi-Step Issue".

> **Note:**  If an issue is created in one project and it is later realized
> that it belongs in a different one, you can move the issue from one project
> to another via a button at the bottom of the right-hand panel on the issue
> page.

Newly-created issues will automatically go in the **Open** column on the
[Kanban board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards).

[↑ Contents](#contents)

### Markdown

[Markdown](https://en.wikipedia.org/wiki/Markdown) is a lightweight markup
language with plain text formatting syntax.  GitLab uses a form of it for
rendering issue and merge request descriptions and comments, wiki pages, and
any files in your repositories with a `.md` extension (such as this one).  For
more details on what's possible with GitLab-flavored Markdown, [see GitLab's
documentation on it](https://docs.gitlab.com/ee/user/markdown.html).

[↑ Contents](#contents)

### Issue Templates

In the top-left of the page when creating an issue, use the drop-down menu to
select a template (~"Bug Report", ~Discussion, ~Documentation,
~"Feature Request", ~Task, or ~Question).  Follow the instructions that appear
in the template.

[↑ Contents](#contents)

### Related Issues

Once an issue has been created, a **Related issues** box will appear below the
issue Description.  You can use this feature to indicate which issues are
connected to each other.

[↑ Contents](#contents)

### Labels

The
[Charon > Labels](https://cee-gitlab.sandia.gov/groups/Charon/-/labels)
page shows you all the labels we use, along with their descriptions.  Labels
corresponding to the issue templates will be applied automatically when an
issue is created.  Other labels (e.g., ~"Multi-Step Issue") can be added as
appropriate.  When creating new issues, please choose between
~"Priority:  Low", ~"Priority:  Medium", and ~"Priority:  High".

[↑ Contents](#contents)

#### Blocked Issues

When forward progress on an issue becomes blocked, please add the
~"Stage:  Blocked" label, and then indicate in a comment what issues are
blocking and why.

[↑ Contents](#contents)



## Working Issues

### Breaking Issues Down

Before work can begin on an issue, it must be dragged into ~"Stage:  Breakdown"
on the [Kanban board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards).
During this stage, the issue must be fully defined such that work can actually
start on it.  Make sure you include enough detail in the description such that
another team member could attempt the work without having to ask too many
questions.

The issue templates include a table for calculating an issue's weight. 
Experience suggests that a weight of six to eight is a reasonable amount of
work for a two to three week period.  Larger issues should likely be broken
down by creating sub-issues that are more bite-size.  Once an issue is small
and well-defined enough, assign it to the appropriate team member and move it
into the ~"Stage:  Ready to Work" column.

If you wind up breaking an issue down into sub-issues, there are two
possibilities for what can happen to the initial issue.  First, if it is useful
to keep the initial issue open for the sake of tracking its separate pieces,
label it a ~"Multi-Step Issue" and remove its weight.  When one of its
sub-issues moves into ~"Stage:  Ready to Work", ~"Stage:  In Development", or
~"Stage:  Blocked", move the ~"Multi-Step Issue" into the same column on the
[Kanban board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards).  Once all
the sub-issues have been closed, close the ~"Multi-Step Issue".

Alternatively, if the sub-issues that were created during ~"Stage:  Breakdown"
are self-sufficient, and there is no reason to keep the parent issue around,
close the parent issue and mention in a comment the sub-issues it was broken
into.

[↑ Contents](#contents)

### When Work Begins

First move the issue to ~"Stage:  In Development" on the [Kanban
board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards).  Next make sure
your local `develop` branch is up-to-date with
```bash
git checkout develop
git pull --ff-only
```

> **Note:**  You should never be making commits on your `develop` branch.  The
> `--ff-only` flag ensures you only update your local `develop` branch if it
> can be fast-forwarded.

Once `develop` is updated, you then create a feature branch off of it with `git
checkout -b <branchName>`.  The recommended branch naming convention is to use
the issue number, following by a hyphen, followed by the issue title, all
lowercase, omitting special characters, and replacing spaces with hyphens.  For
instance, if issue number 123 has "Implement Awesome New Feature" as the title,
the corresponding branch name would be `123-implement-awesome-new-feature`.

[↑ Contents](#contents)

#### If Working in Multiple Repositories

If the work you're doing will require changes to more than one of the Charon
repositories at the same time, the procedure is modified slightly.  First move
all the appropriate issues to ~"Stage:  In Development" on the [Kanban
board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards).  Then you can
[use gitdist](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/Using-gitdist)
to simplify your workflow a bit.
```bash
gitdist checkout develop
gitdist pull --ff-only
gitdist checkout -b <branchName>
```

Since your work spans multiple issues with separate issue numbers and titles,
choose a descriptive branch name that encompases all the work being done.

[↑ Contents](#contents)

### As Work Continues

Do whatever work is necessary to address the issue you're tackling.  Break your
work into logical, compilable commits.  Feel free to commit small chunks early
and often in your local repository and then use `git rebase -i` to reorganize
your commits before sharing.

[↑ Contents](#contents)

#### Commit Messages

Make sure your commit messages reference the appropriate issue numbers using
the `#<issueNumber>` syntax.  The first line of the commit message should be a
descriptive title, limited to 50 characters.  This is then followed by a blank
line, and then the rest of the commit message is a description of the changes,
particularly why they were made, limited to 72 characters wide.

[↑ Contents](#contents)

#### Doxygen

Charon uses [Doxygen](http://www.doxygen.nl) to generate documentation from
annotated source code.  Please see [this wiki
page](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/Doxygen) for our
Doxygen guidelines.

[↑ Contents](#contents)

### When Work is Complete

While working on your feature in your local `<branchName>` branch, other
commits will likely make it into the remote `develop` branch.  There are a
variety of ways to merge these changes into your local feature branch.  One
possibility is
```bash
git checkout develop
git pull --ff-only
git checkout <branchName>
git rebase develop
```

though there are others that are equally valid.

> **Note:**  If you're modifying multiple repositories, you can substitute
> `gitdist` for `git` in the commands above.

To ensure your changes haven't broken anything, you'll want to run `ctest` in
your Charon build directory.  You may want to use the command line options
`--output-on-failure --output-to-root-rank-only=-1` so you can see output in
your terminal if any tests happen to fail.

Once all is well, [create a merge request](#merge-requests) (see below).

[↑ Contents](#contents)

### Closing Old Issues

If at any point you encounter an issue that will not be worked in the
foreseeable future, it is worthwhile to close the issue such that we can
maintain a reasonable backlog of upcoming work.  Do be sure to include in the
comments some explanation as to why the issue won't be addressed.

[↑ Contents](#contents)



## Merge Requests

The only way changes get into `develop` is through merge requests.  When you've
completed work on an issue, push your branch to the remote with `git push -u
<remoteName> <branchName>`, and then create a merge request, selecting a
template corresponding to the issue you've worked on.  On the [Kanban
board](https://cee-gitlab.sandia.gov/groups/Charon/-/boards), drag your issue
into ~"Stage:  Under Review".

[↑ Contents](#contents)

### Reviewers

We recommend having your merge request reviewed by at least two other team
members.  The first should be someone who is knowledgable about the code that
you're changing&mdash;this is to make sure you don't accidentally do something
foolish.  The second should be someone who knows little about the code you're
touching&mdash;this is to spread the knowledge of how the code works throughout
the team.  Work with your reviewers to get your changes into an acceptable
state.

[↑ Contents](#contents)

### Work-in-Progress

You may wish to have your changes reviewed by colleagues before they are ready
to be merged into `develop`.  To do so, create a merge request as usual, but
insert "WIP:" at the beginning of the Title.  GitLab will not allow you to
merge a WIP request.

[↑ Contents](#contents)

### Related Merge Requests

If addressing a task has involved work in multiple repositories, such that
multiple merge requests have been created, be sure to mention these
corresponding merge requests in the description.  If they need to be merged in
a certain order (one is needed before the other will work), or if they must be
merged "simultaneously" (all are needed for the code to work again), be sure to
indicate that and make sure to follow those instructions when merging.

[↑ Contents](#contents)

### Merging

When the review is finished and changes are ready to be merged into `develop`:
1. Rebase your feature branch on top of the latest `develop`.
1. Squash your feature branch down to a single commit.
1. Merge the request.
1. Return to the issue the merge request addressed and provide some evidence in
   a comment that the **Done Criteria** have been met.

> **Note:**  The motivation here is we want the code to build and tests to pass
> for every commit that makes it into `develop`, and we'd like a history that
> is as linear as possible.  This makes finding problems with `git bisect`
> significantly easier.  However, there may be situations in which you don't
> want to squash down to a single commit.  In such a case, squash down to the
> smallest number of commits that makes sense, ensuring the code builds and
> tests pass for each commit.

[↑ Contents](#contents)
