# Pumas Optimal Design Workshop

This is a workshop for using Pumas' `OptimalDesign.jl` package.

## Pure Minimalistic Theme

This template is based on the [Pure Minimalistic Theme](https://github.com/kai-tub/latex-beamer-pure-minimalistic)

## Required Packages

```bash
beamer textpos babel biblatex inputenc csquotes xpatch tikz pgfplots silence appendixnumberbeamer fira fontaxes mwe noto
```

## Plots

If you are plotting stuff and want to use the dark theme you'll probably want to add this to the preamble:

```latex
\usepackage{pgfplots}
\pgfplotsset{height=7cm,  % only if needed
             width=12cm,  % only if needed
             compat=1.18, % only if needed
             legend style = {fill = black, draw = white}}
```

 
## License
The Beamer template software is released under the GNU GPL v3.0 [license](LICENSE).
Pumas' logos and software are copyrighted and all rights reserved.
The content of the slides and code are licensed under a Creative Commons license.
