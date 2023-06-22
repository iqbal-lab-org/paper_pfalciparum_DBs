from gramtools import version

_, version_dict = version.report()
GMTOOLS_COMMIT = version_dict.get("last_git_commit_hash")
if GMTOOLS_COMMIT is None:
    raise ValueError(
        "Could not get gramtools commit hash from gramtools version. Gramtools is probably not compiled."
    )
else:
    GMTOOLS_COMMIT = GMTOOLS_COMMIT[:8]

print(GMTOOLS_COMMIT)
