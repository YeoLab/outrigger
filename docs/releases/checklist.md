# Release checklist

- 1. Check that version numbers in `outrigger/outrigger/__init__.py` and `outrigger/setup.py` match
- 2. Add release notes in `outrigger/docs/releases`
- 3. Copy release notes to `outrigger/HISTORY.rst`
- 4. Convert `README.md` to RST: `pandoc --from=markdown_github --to=rst README.md > README.rst`
- 5. Create an annotated tag for git and push it:

```
git tag -a v0.2.1 -m "v0.2.1 - Release *with* requirements.txt"
git push origin v0.2.1
```

- 6. If you made changes to the tag, remove the remote tag and force re-add the local tag:

```
git push origin :refs/tags/v0.2.1
git tag -fa v0.2.1
git push origin master --tags
```


- 7. Do a test run of uploading to PyPI test:
```
python setup.py register -r pypitest
python setup.py sdist upload -r pypitest
```
