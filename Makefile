.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "dist - package"
	@echo "install - install the package to the active Python's site-packages"
	@echo "tasic2016 - Run outrigger on mouse brain single-cell dataset"
	@echo "treutlein2014 - Run outrigger on mouse lung single-cell dataset"
	@echo "arabdopsis - Run outrigger to check for ENSEMBL compatability"

clean: clean-build clean-pyc clean-test clean-output

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

clean-output:
	rm -rf *outrigger_output

lint:
	flake8 --exclude outrigger/external,doc outrigger

test: clean-pyc
	py.test outrigger

coverage: clean-pyc
	coverage run --source outrigger --omit="*/test*" --module py.test
	coverage report --show-missing

docs:
	rm -f docs/outrigger.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ outrigger
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html

release: clean
	python setup.py sdist upload
	python setup.py bdist_wheel upload

dist: clean
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean
	python setup.py install

N_JOBS = -1

tasic2016:
	rm -rf tasic2016_outrigger_output
	outrigger index --sj-out-tab outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/*SJ.out.tab --gtf outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.subset.gtf --output tasic2016_outrigger_output --n-jobs $(N_JOBS)
	outrigger validate --genome ~/genomes/mm10/mm10.chrom.sizes --fasta ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa --output tasic2016_outrigger_output
	outrigger psi --output tasic2016_outrigger_output --n-jobs $(N_JOBS)

tasic2016_bigger:
	OUTPUT=tasic2016_bigger_outrigger_output
	rm -rf $$OUTPUT
	outrigger index --sj-out-tab outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/originals/CAV_LP_Ipsi_tdTpos_cell1*SJ.out.tab --gtf outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.subset.gtf --output $$OUTPUT --n-jobs $(N_JOBS)
	outrigger validate --genome ~/genomes/mm10/mm10.chrom.sizes --fasta ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa --output $$OUTPUT
	outrigger psi --output $$DIR --n-jobs $(N_JOBS)

treutlein2014: clean-output
	rm -rf treutlein2014
	outrigger index \
		-j outrigger/tests/data/io/star/treutlein2014/sj_out_tab/* \
		-g outrigger/tests/data/io/gtf/treutlein2014/gencode.vM2.annotation.fgfr2.gtf \
		-o treutlein2014

arabdopsis: clean-output
	outrigger index \
		--sj-out-tab outrigger/tests/data/arabdopsis/unprocessed/rna.chr4.subset.SJ.out.tab \
		--gtf outrigger/tests/data/arabdopsis/unprocessed/Arabidopsis_thaliana.TAIR10.31.chr4.subset.gtf \
		--min-reads 1 --n-jobs $(N_JOBS) \
		--output arabdopsis_outrigger_output
	outrigger psi --n-jobs $(N_JOBS) --output arabdopsis_outrigger_output
