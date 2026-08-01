## Run fast unit tests for scE2G's own Python/R scripts against the sc_e2g
## conda env (the same env those scripts run in), instead of building a
## second named env just to run pytest/testthat.
rule test_python_unit:
	output:
		touch("tests/unit/.pytest.done")
	conda:
		"../envs/sc_e2g.yml"
	shell:
		# pytest exits 5 ("no tests collected") once all unit tests are added
		# incrementally; tolerate only that code, not real failures.
		"pytest tests/unit/python -v || [ $? -eq 5 ]"

rule test_r_unit:
	output:
		touch("tests/unit/.testthat.done")
	conda:
		"../envs/sc_e2g.yml"
	shell:
		# testthat::test_dir() errors "No test files found" once all unit
		# tests are added incrementally; tolerate only that, not real failures.
		"""Rscript -e 'tryCatch(testthat::test_dir("tests/unit/r"), error = function(e) if (conditionMessage(e) == "No test files found") invisible(NULL) else stop(e))'"""
