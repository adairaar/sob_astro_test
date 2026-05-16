# .latexmkrc — force latexmk to use the full TeX Live 2026 xelatex
$xelatex = '/usr/local/texlive/2026/bin/universal-darwin/xelatex -synctex=1 -interaction=nonstopmode -file-line-error %O %S';
$pdf_mode = 5;          # 5 = use xelatex to produce PDF
$postscript_mode = 0;
$dvi_mode = 0;

# Cross-document dependency: sob_paper.tex uses \externaldocument{sob_si}
# (via the xr package) so sob_si.aux must exist before sob_paper is compiled.
# latexmk handles this automatically when both files are tracked together.
# If compiling sob_paper.tex in isolation, run sob_si.tex first.
