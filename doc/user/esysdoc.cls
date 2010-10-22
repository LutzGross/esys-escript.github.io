%
%  A class for ESSCC documentation
%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{esysdoc}
\LoadClass[a4paper, twoside]{book}

\usepackage[a4paper, left=2.5cm, right=2.5cm, nohead, vcentering, scale=0.8]{geometry}
\usepackage{times}
% \usepackage{fancyvrb}	% Hopefully this will give us our code environments
\usepackage{listings}
\usepackage{ifpdf}

%\usepackage{avant}
%\renewcommand{\familydefault}{\sfdefault}	% Switch to san seriff
\usepackage{url}
\usepackage{xspace}
\usepackage{graphicx,color}
\usepackage{amsmath}	%For equation*
\usepackage{amsfonts,amssymb} % see http://www.ams.org/tex/amslatex.html
\usepackage[nohyphen, strings]{underscore}

%Blue was chosen here because it looks ok when printed
\usepackage[pdfborder={0 0 0},pdfdisplaydoctitle=true,linktocpage=false, colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}

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
% \ChRuleWidth{1pt}
\ChNameVar{\raggedleft\large\sffamily}	% The word chapter
\ChNumVar{\raggedleft\bfseries\Large}	% The number of the chapter
\ChTitleVar{\raggedleft \Huge\sffamily} % title of the chapter

%Ensures that latex doesn't have an error if we don't specify the version
\providecommand{\RepVersion}{Unknown\xspace}

\DeclareUrlCommand\email{}

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
\textsf{
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
}
\end{titlepage}}


%Please do not use this one
\newcommand*{\strong}{}



