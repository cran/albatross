% In older R, macros can't take multi-line arguments. \[d]eqn3 is
% expected to take multi-line arguments, so force R CMD build to cache
% the parse tree for the whole man page. \bibentry doesn't need
% \forcebuild because it consists of \Sexpr[stage=build] itself.
% Macros still can't be multi-line themselves, hence no line breaks.
\newcommand{\forcebuild}{\Sexpr[results=hide,stage=build]{}}
% Note that we still need to cache the build-time parse tree despite
% \Sexpr[stage=install] present because install-time R could be too old
% to understand our use of \eqn3. \ifelse can take a \Sexpr as a
% condition, but not a macro that expands to a \Sexpr, hence the
% duplication of the condition.
\newcommand{\eqn3}{\forcebuild\ifelse{\Sexpr[results=rd,stage=render]{albatross:::.eqn3legacy()}}{\out{<i>#3</i>}}{\eqn{#1}{#2}}}
\newcommand{\deqn3}{\forcebuild\ifelse{\Sexpr[results=rd,stage=render]{albatross:::.eqn3legacy()}}{\out{<p style="text-align: center;"><i>#3</i></p>}}{\deqn{#1}{#2}}}
% This could be enhanced to both create and insert a PNG/SVG/PDF figure.
\newcommand{\makefigure}{\Sexpr[results=rd,stage=build]{if (!dir.exists('man/figures')) dir.create('man/figures'); {#1}; ' '}}
\newcommand{\cite}{\Sexpr[results=rd,stage=build]{albatross:::.Rdcite('#1')}}
\newcommand{\bibliography}{\Sexpr[results=rd,stage=build]{albatross:::.Rdbibliography()}}
\newcommand{\reference}{\Sexpr[results=rd,stage=build]{albatross:::.Rdreference('#1')}}
