% \vignette{topic} references a PDF vignette in this package via simple substitution, 
% using relative paths compatible with static & dynamic HTML help:
\newcommand{\vignette}{\ifelse{html}{\href{../doc/#1.pdf}{\code{vignette("#1")}}}{\code{vignette("#1")}}}
