PDF=main.pdf

all: $(PDF)

$(PDF): main.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C

check: $(PDF)
	@grep -n "LaTeX Warning" main.log || true
	@grep -n "Undefined" main.log || true
