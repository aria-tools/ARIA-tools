# Contributing Guidelines #

This document is inspired by similar instructions from PySAR, ISCE, gdal and jupyterhub.
These are several ways to contribute to the ARIA-tools framework:

*   Submitting bug reports and feature requests in ARIA-tools
*   Writing tutorials or jupyter-notebooks in ARIA-tools-docs
*   Fixing typos, code and improving documentation
*   Writing code for everyone to use

If you get stuck at any point you can create an [issue on GitHub](https://github.com/aria-tools/aria-tools/issues).

For more information on contributing to open source projects, [GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point if you are new to version control.

## Git workflows ##

### Setting up the development environment ###

Fork aria-tools/aria-tools from GitHub UI, and then

```.tcsh
git clone https://github.com/aria-tools/aria-tools.git
cd aria-tools
git remote add my_user_name https://github.com/my_user_name/aria-tools.git
```

### Setting up the documentation environment ###

Fork aria-tools/aria-tools-docs from GitHub UI, and then

```.tcsh
git clone https://github.com/aria-tools/aria-tools-docs.git
cd aria-tools-docs
git remote add my_user_name https://github.com/my_user_name/aria-tools-docs.git
```

### Updating your local master against upstream master ###

```.tcsh
git checkout master
git fetch origin
# Be careful: this will loose all local changes you might have done now
git reset --hard origin/master
```

### Working with a feature branch ###

[Here](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history) is a great tutorial if you are new to rewriting history with git.

```.tcsh
git checkout master
(potentially update your local master against upstream, as described above)
git checkout -b my_new_feature_branch

# do work. For example:
git add my_new_file
git add my_modifid_message
git rm old_file
git commit -a

# you may need to resynchronize against master if you need some bugfix
# or new capability that has been added to master since you created your
# branch
git fetch origin
git rebase origin/master

# At end of your work, make sure history is reasonable by folding non
# significant commits into a consistent set
git rebase -i master (use 'fixup' for example to merge several commits together,
and 'reword' to modify commit messages)

# or alternatively, in case there is a big number of commits and marking
# all them as 'fixup' is tedious
git fetch origin
git rebase origin/master
git reset --soft origin/master
git commit -a -m "Put here the synthetic commit message"

# push your branch
git push my_user_name my_new_feature_branch
```

### Issue a pull request from GitHub UI ###
commit locally and push. To get a reasonable history, you may need to

```.tcsh
git rebase -i master
```

, in which case you will have to force-push your branch with

```.tcsh
git push -f origin my_new_feature_branch
```

Once a pull request is issued it will be reviewed by multiple members before it will be approved and integrated into the main.

### Things you should **NOT** do
(For anyone with push rights to [https://github.com/aria-tools/ARIA-tools] or [https://github.com/aria-tools/ARIA-tools-docs]) Never modify a commit or the history of anything that has been committed to [https://github.com/aria-tools/ARIA-tools] and [https://github.com/aria-tools/ARIA-tools-docs].
