\newcommand{\forcebuild}{\Sexpr[results=hide,stage=build]{}}
\newcommand{\eqn3}{\forcebuild\ifelse{html}{\out{<i>#3</i>}}{\eqn{#1}{#2}}}
\newcommand{\deqn3}{\forcebuild\ifelse{html}{\out{<p style="text-align: center;"><i>#3</i></p>}}{\deqn{#1}{#2}}}
\newcommand{\bibentry}{\Sexpr[results=rd,stage=build]{tools::toRd(bibentry(#1))}}
