%
%  A class for ESSCC documentation
%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{esysdoc}
\LoadClass[a4paper, twoside]{book}

\usepackage[a4paper, left=2.5cm, right=2.5cm]{geometry}
\usepackage{times}
% \usepackage{fancyvrb}	% Hopefully this will give us our code environments
\usepackage{listings}

%\usepackage{avant}
\renewcommand{\familydefault}{\sfdefault}	% Switch to san seriff
\usepackage{url}
\usepackage{xspace}
\usepackage{graphics}
\usepackage{amsmath}	%For equation*
\usepackage{amsfonts,amssymb} % see http://www.ams.org/tex/amslatex.html
\usepackage[nohyphen, strings]{underscore}
\usepackage[colorlinks=false,pdfborder={0 0 0},pdfdisplaydoctitle=true,linktocpage=true]{hyperref}

%Ensures that latex doesn't have an error if we don't specify the version
\providecommand{\RepVersion}{Unknown\xspace}



%This is FAR from perfect for one thing it doesn't allow underscores to be used in filenames

%Now we mess with the underscore char
% {
% \catcode`|=8   % Since we will be redefining _ we need another way to indicate subscript in math mode
%                % This change only holds for the definition since it disappears at the close brace
% \catcode`_=\active	% define _ as a char which triggers a macro (which is what active means)
% \gdef\_{%		% A global definition (persists after } ) for _ when it is a macro
% \ifmmode |		% If we are in math mode, then output subscript char
% \else
% \textunderscore		% else output underscore symbol
% \fi
% }
% }
% \catcode`_=\active	% Set _ to active at the top level


% \newcommand*{\module}[1]{\texttt{#1}}
% \newcommand*{\class}[1]{\texttt{#1}}
% \newcommand*{\var}[1]{\textsl{#1}}
% \newcommand*{\method}[1]{\texttt{#1}}


%Use the \Var form if you need to include % in the URL
%Please note that the \Var form does not work in maths mode
% \DeclareUrlCommand\Var{}
%Same here with \function except \Function is already in use for a FunctionSpace
% \DeclareUrlCommand\FunctionName{}

% \DeclareUrlCommand\umodule{}
% \DeclareUrlCommand\class{}
% \DeclareUrlCommand\method{}
\DeclareUrlCommand\email{}
% \DeclareUrlCommand\program{}
% \DeclareUrlCommand\env{}


% \DeclareUrlCommand\member{}
% \DeclareUrlCommand\file{}

% \DeclareUrlCommand\kbd{}

% \DeclareUrlCommand\usafe{}	%intended as a wrapper to make stuff underscrore safe
% \DeclareUrlCommand\constant{}

%To safely deal with a module that has _ in its name do this:
% \urldef{\rMATPLOTLIB}\umodule{matplotlib}
% \newcommand{\MATPLOTLIB}{\rMATPLOTLIB\index{visualization!matplotlib}\index{matplotlib}\xspace}

\newcommand{\module}[1]{\mbox{\texttt{#1}}}

\newcommand*{\optional}[1]{[ #1 ]}
%This should be converted into a two arg form so we can have different fonts for different parts
\newcommand*{\programopt}[1]{#1}

\newenvironment{classdesc}[3]{\begin{description}\item[class #1]\mbox{}\\}{\end{description}}
\newenvironment{methoddesc}[3][{}]{\begin{description}\item[\textbf{#2}(#3)]\mbox{}\\}{\end{description}}
\newenvironment{funcdesc}[2]{\begin{description}\item[\textbf{#1}(#2)]\mbox{}\\}{\end{description}}
\newenvironment{memberdesc}[2][{}]{\begin{description}\item[\textbf{#2}]\mbox{}\\}{\end{description}}
\newenvironment{datadesc}[1]{\begin{description}\item[\textbf{#1}]\mbox{}\\}{\end{description}}

\newcommand{\file}[1]{\texttt{#1}\xspace}
\newcommand{\program}[1]{\textbf{\texttt{#1}}\xspace}
% \newcommand{\env}[1]{\texttt{#1}\xspace}
\newcommand{\var}[1]{\texttt{#1}\xspace}
\newcommand{\kbd}[1]{\texttt{#1}}
\newcommand{\function}[1]{\texttt{#1}}
\newcommand{\class}[1]{\texttt{#1}}
\newcommand{\constant}[1]{\textsl{#1}}
\newcommand{\method}[1]{\texttt{#1}}
\newcommand{\member}[1]{\texttt{#1}}
\newcommand{\code}[1]{\mbox{\texttt{#1}}}   % Because people try to put code in maths mode


\newenvironment{abstract}{\newpage\mbox{}\vfill\begin{center}\textbf{Abstract}\end{center}}{\vfill\newpage}

\newcommand{\authoraddress}[1]{\gdef \@authoraddress {#1}}

\newcommand{\release}[1]{\gdef \@ESreleasename {#1}}
\release{nightly}	%Just so we have a default

\newcommand*{\subtitleLineOne}{\textsl{Release - \@ESreleasename}}
\newcommand*{\subtitleLineTwo}{\textsl{(r\RepVersion)}}

\renewcommand{\maketitle}{
\begin{titlepage}
\begin{raggedleft}
\Huge
\hrule
\mbox{}\\
\mbox{}\\
\@title
\mbox{}\\
\LARGE
\mbox{}\\
\subtitleLineOne\\
\subtitleLineTwo\\
\vfill
\@author
\vfill
\vfill
\large
\@date\\
\vfill
\@authoraddress\mbox{}\\
\end{raggedleft}

\end{titlepage}}


%Please do not use this one
\newcommand*{\strong}{}

\usepackage{fancyhdr}
\renewcommand{\footrulewidth}{0.4pt}
\renewcommand{\headrulewidth}{0pt}

\fancypagestyle{plain}{%
\fancyhead{}
\cfoot{}
\fancyfoot[LE,RO]{\thepage}	%left on even pages right on odd pages
\fancyfoot[RE]{\nouppercase{\rightmark}}
\fancyfoot[LO]{\nouppercase{\leftmark}}
}
\pagestyle{plain}

\usepackage[Bjarne]{fncychap}

