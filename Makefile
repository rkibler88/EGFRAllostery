.PHONY: all figures verify pathways knockout manuscript clean
PY := python

all: verify figures

verify:                      ## recompute every manuscript number from data/processed
	$(PY) analysis/verify_numbers.py

figures:                     ## regenerate all data-driven figure panels into results/figures
	$(PY) analysis/make_figure1.py
	$(PY) analysis/make_figure2.py
	$(PY) analysis/make_figure3.py
	$(PY) analysis/make_figure4.py
	$(PY) analysis/make_figure5.py
	$(PY) analysis/make_figure6.py
	$(PY) analysis/make_figure7_sectors.py
	$(PY) analysis/make_figure8_beta_sheet.py
	$(PY) analysis/make_figureS1_replicates.py

pathways:                    ## shortest-path analysis (L858R, alphaG-AH loop to tail)
	$(PY) analysis/compute_pathways.py

knockout:                    ## global-efficiency knockout table
	$(PY) analysis/network_knockout.py

manuscript:                  ## build the PDF (requires a LaTeX installation)
	cd manuscript && latexmk -pdf frontiers.tex

clean:
	rm -f manuscript/*.aux manuscript/*.bbl manuscript/*.blg manuscript/*.log manuscript/*.out manuscript/*.fls manuscript/*.fdb_latexmk
