# Release checklist

- [ ] 1. Check that version numbers in `outrigger/outrigger/__init__.py` and `outrigger/setup.py` match
- [ ] 2. Add release notes in `outrigger/docs/releases`
- [ ] 3. Copy release notes to `outrigger/HISTORY.rst`
- [ ] 4. Convert `README.md` to RST: `pandoc --from=markdown_github --to=rst README.md > README.rst`
- [ ] 5. Create an annotated tag for git and push it:

```
git tag -a v0.2.1 -m "v0.2.1 - Release *with* requirements.txt"
git push origin v0.2.1
```

- [ ] 6. If you messed up and made the tag too early and later made changes that should be added to the tag, remove the remote tag and force re-add the local tag:

```
git push origin :refs/tags/v0.2.1
git tag -fa v0.2.1
git push origin master --tags
```


- [ ] 7. Do a test run of uploading to PyPI test:
```
python setup.py register -r pypitest
python setup.py sdist upload -r pypitest
```

- [ ] 8. If there's an error in either the `HISTORY.rst` or `README.rst` file, you will get an error. To narrow down this error, look for the corresponding line in `outrigger.egg-info`. It won't be exactly the same line because of the header but it will be closer since `PKG-INFO` concatenates your `README.rst` and `HISTORY.rst` files.

```
warning: unexpected indent on line 615 blah blah I'm pypi RST and I'm way too strict
```

- [ ] 10. Check the [PyPI test server](https://testpypi.python.org/pypi) and make sure everything is there and the RST is rendered correctly.
- [ ] 11. Do a test installation in a `conda` environment using the PyPI test server

```
conda create --yes -n outrigger_pypi_test_v3 --file conda_requirements.txt
source activate outrigger_pypi_test_v3
pip install --index-url https://testpypi.python.org/pypi outrigger --extra-index-url https://pypi.python.org/simple
```

- [ ] 12. Check that the installation was successful. `outrigger -h` should have the following output:

```
$ outrigger -h
usage: outrigger [-h] {index,validate,psi} ...

Calculate "percent-spliced in" (Psi) scores of alternative splicing on a *de
novo*, custom-built splicing index

positional arguments:
  {index,validate,psi}  Sub-commands
    index               Build an index of splicing events using a graph
                        database on your junction reads and an annotation
    validate            Ensure that the splicing events found all have the
                        correct splice sites
    psi                 Calculate "percent spliced-in" (Psi) values using the
                        splicing event index built with "outrigger index"

optional arguments:
  -h, --help            show this help message and exit
```

- [ ] 13. Upload to PyPI:

```
python setup.py register -r pypi
python setup.py sdist upload -r pypi
```

- [ ] 14. Upload to Anaconda.org:

```
anaconda login
anaconda upload dist/*.tar.gz
```
