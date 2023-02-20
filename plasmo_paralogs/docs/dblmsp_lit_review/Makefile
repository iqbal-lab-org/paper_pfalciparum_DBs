build_dir := built
fname := report
target := $(fname).pdf
build_target := $(build_dir)/$(fname).pdf

latex_build := mkdir -p $(build_dir) && pdflatex -output-directory=$(build_dir) -jobname=$(fname) $(fname).tex

all: $(target)

supmat:
	pdflatex -output-directory=supplementary supplementary/supplementary.tex


lite:
	$(latex_build)
	cp $(build_target) $(target)

$(build_dir)/$(fname).bbl: $(fname).tex $(fname).bib
	$(latex_build)
	biber $(build_dir)/$(fname)

$(target): $(fname).tex $(build_dir)/$(fname).bbl
	$(latex_build)
	$(latex_build)
	cp $(build_target) $(target)

clean: 
	rm $(build_dir)/*

.PHONY: all clean
