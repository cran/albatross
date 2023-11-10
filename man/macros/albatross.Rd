% In older R, macros can't take multi-line arguments. \[d]eqn3 is
% expected to take multi-line arguments, so force R CMD build to cache
% the parse tree for the whole man page. \bibentry doesn't need
% \forcebuild because it consists of \Sexpr[stage=build] itself.
% Macros still can't be multi-line themselves, hence no line breaks.
\newcommand{\forcebuild}{\Sexpr[results=hide,stage=build]{}}
% Note that we still need to cache the build-time parse tree despite
% \Sexpr[stage=render] present because install-time R could be too old
% to understand our use of \eqn3. \ifelse can take a \Sexpr as a
% condition, but not a macro that expands to a \Sexpr, hence the
% duplication of the condition.
\newcommand{\eqn3}{\forcebuild\ifelse{\Sexpr[results=rd,stage=render]{albatross:::.eqn3legacy()}}{\out{<i>#3</i>}}{\eqn{#1}{#2}}}
\newcommand{\deqn3}{\forcebuild\ifelse{\Sexpr[results=rd,stage=render]{albatross:::.eqn3legacy()}}{\out{<p style="text-align: center;"><i>#3</i></p>}}{\deqn{#1}{#2}}}
% This should be changed to create either PNG or PDF+SVG and then insert
% it using conditional syntax if needed. Try to provide alt text if
% possible. Note the workaround for \Sexpr processing in older R
% versions: it has to be a space (not an empty string) with results=rd.
\newcommand{\makefigure}{\Sexpr[results=rd,stage=build]{if (!dir.exists('man/figures')) dir.create('man/figures'); {#1}; ' '}}
% Insert a short reference here, remember the key for \bibliography
\newcommand{\cite}{\Sexpr[results=rd,stage=build]{albatross:::.Rdcite('#1')}}
% Insert full references for every \cited item
\newcommand{\bibliography}{\Sexpr[results=rd,stage=build]{albatross:::.Rdbibliography()}}
% Insert a full reference here
\newcommand{\reference}{\Sexpr[results=rd,stage=build]{albatross:::.Rdreference('#1')}}
% Post R-4.3, \var{foo} renders as <foo> in plain text output
\newcommand{\v}{\ifelse{\Sexpr[results=rd,stage=install]{if (getRversion() > '4.3') 'FALSE' else 'text'}}{<#1>}{\var{#1}}}
