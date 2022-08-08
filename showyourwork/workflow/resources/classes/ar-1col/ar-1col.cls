%% Modified 2022-08-08 by @dfm to rename "marginnote" environment to "armarginnote" for showyourwork compatibility
%%
%% ar.cls - version 1.1 
%% Aptara Inc., dated 07 Mar. 2014
%%
%% Version 1.1 (History)
%% ---------------------
%% 1) Implemented the color.sty to avoide the  
%%    color related problems.
%% 2) Used the ifpdf.sty for providing the PDF Paper
%%    Size in case of pdflatex.
%% 3) Seperation between author and year removed for Harvard style reference (unnumbered)
%% 4) Removed the clashing of marginnote environment with amsmath.sty 
%% 5) Introduced an option to change format of Equation Number and a command ``\firstpagenote''to appear text on first page
%% 
%% For AR journals
%%
%% Steps to compile: latex latex latex
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{ar}[2013/03/07 v1.1 Standard LaTeX document class]
%
\newif\if@restonecol
\newif\if@DotinEqNum\global\@DotinEqNumtrue
%
% Global Variable Declaration
%
\def\doi#1{\gdef\@doi{#1}}\doi{}%
\def\fstpage#1{\gdef\@fstpage{#1}}\fstpage{}%
\def\endpage#1{\gdef\@endpage{#1}}\endpage{}%
\def\jvol#1{\gdef\@jvol{#1}}\jvol{00}%
\def\jyear#1{\gdef\@jyear{#1}}\jyear{0000}%
\def\jname#1{\gdef\@jname{#1}}\jname{xxxxxx}%
%
% Font size declaration
%
\def\rhfont{\fontsize{7.5}{7.5}\itshape\selectfont}
\def\titlefont{\fontsize{20}{24}\selectfont\bfseries\leftskip14.25pc\rightskip0\p@ plus1fill}
\def\authorfont{\fontsize{12}{15}\selectfont\bfseries\leftskip14.25pc\rightskip0\p@ plus1fill\mathversion{bold}}
\def\abstractfont{\fontsize{9}{12}\selectfont}
\def\keywordsfont{\fontsize{9}{12}\selectfont\rightskip0\p@ plus1fill\mathversion{normal}}
\def\affilfont{\fontsize{7}{9}\selectfont\leftskip14.25pc\rightskip0\p@ plus1fil{}}
\def\jinfofont{\fontsize{7}{9}\selectfont\rightskip12\p@ plus1fil}
\def\foliofont{\fontsize{7.5}{7.5}\itshape\selectfont}
\def\figcaptionfont{\fontsize{8}{10}\selectfont\mathversion{normal}}
\def\figcaptionnumfont{\sffamily\fontsize{8}{10}\selectfont\bfseries\mathversion{bold}}
\def\tablecaptionfont{\fontsize{8}{12}\bfseries\selectfont\mathversion{normal}}% 
\def\tablecaptionnumfont{\fontsize{8}{12}\bfseries\selectfont}
\def\tablefont{\fontsize{8}{11}\selectfont\mathversion{normal}}%
\def\tabnotefont{\fontsize{7}{10}\selectfont}%
\def\textboxfont{\normalfont\normalsize\leftskip12\p@\rightskip12\p@}
\def\textboxheadfont{\fontsize{10}{12}\selectfont\bfseries\leftskip12\p@\rightskip12\p@ plus1fill\mathversion{normal}}
%
\def\sectionfont{\sffamily\fontsize{10}{12}\selectfont\bfseries\raggedright\mathversion{bold}}
\def\subsectionfont{\sffamily\fontsize{10}{12}\selectfont\bfseries\raggedright\mathversion{bold}}
\def\subsubsectionfont{\sffamily\fontsize{9}{12}\bfseries\selectfont\mathversion{bold}}
\def\paragraphfont{\bfseries\itshape}
\def\subparagraphfont{\itshape}
\def\extractfont{\fontsize{8}{12}\selectfont\mathversion{normal}}
\def\armarginnotefont{\fontsize{8}{10}\selectfont\rightskip0\p@ plus1fill\mathversion{normal}}
\def\bibarmarginnotefont{\fontsize{7}{10}\selectfont\bfseries\rightskip0\p@ plus1fill\mathversion{normal}}
\def\bibliofont{\fontsize{8}{11}\selectfont\mathversion{normal}}
%
%
%% \if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{ChEqNum}{\global\@DotinEqNumtrue}
\ExecuteOptions{letterpaper}
\ProcessOptions
%
\@twosidetrue\@mparswitchtrue
%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@ixpt\@xiipt
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\normalsize
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{9.5}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\setlength\parindent{15\p@}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\headheight{\z@}%
\setlength\headsep{\z@}%
\setlength\topskip{6.9\p@}
\setlength\footskip{24\p@}
\setlength\maxdepth{.5\topskip}
%
\newdimen\typeheight\typeheight47pc%
\newdimen\typewidth\typewidth38pc%
\newlength\extramargin
\setlength\extramargin{7.5pc}
\setlength\textheight{47pc}%
\setlength\textwidth{30.5pc}%
\setlength\marginparsep {10\p@}
\setlength\marginparpush{5\p@}
\setlength\marginparwidth {2.0cm}
\setlength\oddsidemargin{46.5\p@}%
\setlength\evensidemargin{28.5\p@}%
%
\setlength\topmargin{40.5\p@}%
\setlength\footnotesep{7.35\p@}
\setlength{\skip\footins}{19.5\p@ \@plus 4\p@ \@minus 2\p@}
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@ \@plus 1fil}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\setlength\partopsep{0\p@ \@plus 1\p@ \@minus 1\p@}
\def\@listi{\leftmargin\leftmargini
            \parsep 4\p@ \@plus2\p@ \@minus\p@
            \topsep 8\p@ \@plus2\p@ \@minus4\p@
            \itemsep4\p@ \@plus2\p@ \@minus\p@}
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    4\p@ \@plus2\p@ \@minus\p@
              \parsep    2\p@ \@plus\p@  \@minus\p@
              \itemsep   \parsep}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    2\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep   \topsep}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\def\@listI{\leftmargin\leftmargini
            \parsep 0\p@% \@plus2\p@ \@minus\p@
            \topsep 6\p@ \@plus2\p@% \@minus2\p@
            \itemsep0\p@}% \@plus2\p@ \@minus\p@}
%
\newenvironment{unnumlist}{\list{}{\leftmargin15\p@\itemindent-15\p@}}{\endlist}%
%
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@}% \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{4}
\renewcommand\topfraction{1}
\setcounter{bottomnumber}{4}
\renewcommand\bottomfraction{1}
\setcounter{totalnumber}{8}
\renewcommand\textfraction{0.0001}
\renewcommand\floatpagefraction{.91}
\setcounter{dbltopnumber}{4}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.91}%
%
\font\cms=cmsy10 at 5.7\p@
\def\arblcirc{\lower-.6\p@\hbox{\cms\char'17}}%
%
  \def\ps@headings{%
      \let\@mkboth\@gobbletwo
      \def\@oddfoot{\hbox to\typewidth{\hfill\rhfont www.annualreviews.org\kern3.5\p@{\arblcirc}\kern3.5\p@\rightmark\hbox to \extramargin{\hskip1pc\foliofont\thepage\hfill}}}%
      \def\@evenfoot{\hbox to\typewidth{\hspace*{-\extramargin}\hbox to \extramargin{\hfill\foliofont\thepage\hskip1pc}\rhfont\leftmark\hss}}
      \def\@evenhead{}%
      \def\@oddhead{}%
    \def\sectionmark##1{}%
    \def\subsectionmark##1{}}
%
\def\ps@plain{%
    \let\@mkboth\@gobbletwo%
    \def\@oddhead{}%
    \let\@evenhead\@oddhead%
    \def\@oddfoot{\hbox to\typewidth{\hfill\hbox to \extramargin{\hskip2pc\hfill\foliofont\thepage\hfill}}}%
    \def\@evenfoot{\hbox to\typewidth{\hspace*{-\extramargin}\hbox to \extramargin{\hfill\foliofont\thepage\hfill\hskip2pc}}}%
  }%
%
\newcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{}%{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\sffamily\raggedright\noindent
            %\hb@xt@1.8em{\hss\@textsuperscript{\normalfont\@thefnmark}}
                   ##1\vskip3\p@}%
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle%
     \thispagestyle{plain}\@thanks\clearpage
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
%
\usepackage{color,ifpdf}

\definecolor{titlecolor}{cmyk}{0.4,0.4,0.4,0}
\definecolor{headcolor}{cmyk}{0,1.0,1.0,0.30}
\definecolor{shadecolor}{cmyk}{0.03,0.03,0.12,0.0}
\definecolor{fignumcolor}{cmyk}{1.0,0.4,0,0}
\definecolor{marginrulecolor@val}{cmyk}{0,0,0,0.30}
\definecolor{textboxcolor@val}{cmyk}{0.12,0.04,0.08,0.0}

\def\title@color#1{\textcolor{titlecolor}{#1}}
\def\head@color#1{\textcolor{headcolor}{#1}}
\def\shade@color#1{\textcolor{shadecolor}{#1}}
\def\fignum@color#1{\textcolor{fignumcolor}{#1}}
\def\marginrulecolor#1{\textcolor{marginrulecolor@val}{#1}}
\def\textboxcolor#1{\textcolor{textboxcolor@val}{#1}}

%\def\title@color#1{\special{color push cmyk 0.4 0.4 0.4 0}#1\special{color pop}}
%\def\head@color#1{\special{color push cmyk 0 1.0 1.0 0.30}#1\special{color pop}}%
%\def\shade@color#1{\special{color push cmyk 0.03 0.03 0.12 0.0}#1\special{color pop}}%
%\def\fignum@color#1{\special{color push cmyk 1.0 0.4 0 0}#1\special{color pop}}
%\def\marginrulecolor#1{\special{color push cmyk 0 0 0 0.30}#1\special{color pop}}
%\def\textboxcolor#1{\special{color push cmyk 0.12 0.04 0.08 0.0}#1\special{color pop}}

\def\@maketitle{\newpage%
  \null%
  \begingroup\hsize\typewidth\parindent0\p@%
    \let\footnote\thanks
    \nointerlineskip%
    \vskip-6\p@%
    \noindent{\reset@font\titlefont\title@color{\@title}\par}
    \vskip15\p@
    \noindent{\authorfont\@author\par\vskip5\p@}%
    \vfill%
    \global\setbox\z@\vbox{\hsize\typewidth%
      \begin{tabular*}{\typewidth}{@{}p{13.5pc}@{\hskip.75pc}p{22.25pc}@{}}
      \nointerlineskip{\vskip5.4\p@\vbox{\journalinfo}\vrule height0\p@ depth6\p@ width\z@}%
      &\ifvoid\@kwdbox%
         \ifvoid\@absbbox\else{\box\@absbbox}\fi%
      \else%
         {\vskip5.4\p@\box\@kwdbox}\endgraf\nointerlineskip%
         \ifvoid\@absbbox\else\vskip15\p@{\box\@absbbox}\fi%
      \fi\\%
    \end{tabular*}}
    \@tempdima\typewidth\advance\@tempdima77\p@%
    \vbox{\shade@color{\hskip-46\p@\vrule height\ht\z@ width\@tempdima depth\dp\z@}\vskip-\ht\z@\vskip-\dp\z@%
    {\box\z@}}%
  \endgroup}
%
% First page element and layout declaration
%
\def\firstpagenote#1{\gdef\@firstpagenote{#1}}\firstpagenote{}%
\def\journalinfo{{\jinfofont%
    \@jname\ \@jyear. \@jvol:\@fstpage--\@endpage\endgraf%
    %\ifx\@doi\empty\else\vskip6\p@ This article's doi:\break \@doi\endgraf\vskip6\p@\fi%
    \ifx\@doi\empty\else\vskip6\p@ https://doi.org/\@doi\endgraf\vskip6\p@\fi%
   \ifnum\@jyear>2022%
      Copyright\ \copyright\ \@jyear\ by the  author(s).\endgraf All rights reserved\endgraf%
   \else%
      Copyright\ \copyright\ \@jyear\ by Annual Reviews.\endgraf All rights reserved\endgraf%
   \fi%
    \ifx\@firstpagenote\@empty\else\vskip6\p@\@firstpagenote\endgraf\fi%
}}
%
\def\chk@key{no}%
\def\YES{yes}%
\def\affil#1{\par\ifx\chk@key\YES\else\vskip6\p@\fi{\reset@font\affilfont#1\par}\def\chk@key{yes}}
%
\newbox\@absbbox
\newenvironment{abstract}{\parindent0\p@%
      \global\setbox\@absbbox\vbox\bgroup\hsize23.75pc\abstractfont\noindent\ignorespaces%
      \ifx\abstractname\@empty\else{\head@color{\sffamily\bfseries\abstractname}\endgraf\vskip9\p@}\fi%
      }
      {\par\egroup}

%
\newcommand\keywordsname{Keywords}
\newbox\@kwdbox
\newenvironment{keywords}{\parindent0\p@%
      \global\setbox\@kwdbox\vbox\bgroup\hsize23.75pc\keywordsfont\noindent\ignorespaces%
      \ifx\keywordsname\@empty\else{\head@color{\sffamily\bfseries\keywordsname}\endgraf\vskip4\p@}\fi%
      }{\par\egroup}%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand \thepart {\@Roman\c@part}
\renewcommand \thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\newcommand\part{%
   \if@noskipsec \leavevmode \fi
   \par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname\nobreakspace\thepart
       \par\nobreak
     \fi
     \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
%
\def\@seccntformat#1{\head@color{\csname the#1\endcsname.}\hskip3\p@}
%
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-20\p@ \@plus -3\p@ \@minus -1\p@}%
                                   {4\p@}%
                                   {\sectionfont}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-16\p@ \@plus -1\p@ \@minus -.1\p@}%
                                     {4\p@}%
                                     {\subsectionfont}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-12\p@ \@plus -1\p@ \@minus -0.1\p@}%
                                     {-0.5em}%
                                     {\subsubsectionfont}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {13\p@ \@plus 3.25\p@ \@minus 1\p@}%
                                    {-0.5em}%
                                    {\paragraphfont}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {13\p@ \@plus 3.25\p@ \@minus 1\p@}%
                                       {-0.5em}%
                                      {\subparagraphfont}}
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M \head@color{#8}\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}%
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec \head@color{#8.}}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}
%
\def\@ssect#1#2#3#4#5{%
  \@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M \head@color{#5}\@@par}%
    \endgroup
  \else
    \def\@svsechd{#4{\hskip #1\relax #5}}%
  \fi
  \@xsect{#3}}
%
\setlength\leftmargini  {2.65em}
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\newenvironment{extract}
 {\list{}{\leftmargin1pc\topsep12\p@ plus2\p@\rightmargin\leftmargin}\item[]\extractfont}
 {\endlist}
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \gdef\thesection{\@Alph\c@section}}
\setlength\arraycolsep{2\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand \theequation {\@arabic\c@equation}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename\nobreakspace\thefigure}
\newcount\chkfigcnt
\newenvironment{figure}
               {\global\advance\chkfigcnt\@ne\@float{figure}}
               {\label{chk@figure@\the\chkfigcnt}\end@float}
\newenvironment{figure*}
               {\global\advance\chkfigcnt\@ne\@dblfloat{figure}}
               {\label{chk@figure@\the\chkfigcnt}\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{-6\p@}
\setlength\belowcaptionskip{0\p@}
%
\def\FigName{figure}
\long\def\@makecaption#1#2{%
    \ifx\FigName\@captype
      \vskip\abovecaptionskip
      \@makefigurecaption{#1}{#2}%
    \else
      \@maketablecaption{#1}{#2}%
      \vskip\belowcaptionskip
    \fi}
%
\newdimen\chk@fig@width\chk@fig@width\z@
%
\def\@makefigurecaption#1#2{\edef\chk@fig@temp{\getpagerefnumber{chk@figure@\the\chkfigcnt}}%
  \ifdim\chk@fig@width>\typewidth%
    \noindent{\fignum@color{\figcaptionnumfont#1}}\endgraf\vskip-8.5\p@%
    \normalcolor\noindent\hbox to \hsize{\hrulefill}\endgraf\vskip-.9\p@%
    {\figcaptionfont\raggedright#2\par}%
  \else%
    \ifdim\chk@fig@width>\textwidth%
      \vbox{\par\vskip5\p@\hsize\typewidth\ifodd\chk@fig@temp\raggedright\else\leftskip-7.5pc\rightskip7.5pc plus1fil\fi
      \noindent{\fignum@color{\figcaptionnumfont#1}}\endgraf\vskip-8.5\p@%
      \normalcolor\noindent\hbox to \hsize{\hrulefill}\endgraf\vskip-.9\p@%
      {\figcaptionfont#2\par}\par}%
    \else%
      \noindent{\fignum@color{\figcaptionnumfont#1}}\endgraf\vskip-8.5\p@%
      \normalcolor\noindent\hbox to \hsize{\hrulefill}\endgraf\vskip-.9\p@%
      {\figcaptionfont\rightskip0pt plus1fill#2\par}%
    \fi
  \fi%
}
%
\def\@maketablecaption#1#2{%
    \noindent\tablecaptionfont%
    \fignum@color{{\tablecaptionnumfont#1\quad}{\tablecaptionfont#2}}\par
}%
\input epsf.sty
%
\def\tabular{\tablefont\let\@halignto\@empty\@tabular}
\newenvironment{tabnote}{\par\tabnotefont}{\par}
%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.75em}
\newcommand\@tocrmarg{6.5em}
\newcommand\@dotsep{1}
\setcounter{tocdepth}{2}
\newcommand\tableofcontents{\par\vspace*{-15.75\p@}%
\definecolor{shadecolor}{cmyk}{0.03,0.03,0.12,0}
\fboxsep12\p@\fboxrule0\p@\ifodd\c@page\moveright3.75pc\else\moveright-3.75pc\fi\vbox\bgroup\begin{shaded}\@nobreaktrue\hsize36pc
    \section*{\contentsname}
%%         \@mkboth{%
%%            \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
\end{shaded}\egroup}
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \setlength\@tempdima{3em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section{\@dottedtocline{1}{0em}{1.2em}}
\newcommand*\l@subsection{\@dottedtocline{2}{1.2em}{2em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.2em}{2.8em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
%
\def\numberline#1{\hb@xt@\@tempdima{#1.\hfil}}
\def\@dottedtocline#1#2#3#4#5{%
  \ifnum #1>\c@tocdepth \else
    \vskip \z@ \@plus.2\p@
    {\sffamily\small\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     \@tempdima #3\relax
     \advance\leftskip \@tempdima \null\nobreak\hskip -\leftskip
     {#4}\nobreak
      \leaders\hbox{$\m@th
         \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
         mu$}
     \hfill
     \nobreak
     \hb@xt@\@pnumwidth{\hfil#5}%
     \par}%
  \fi}
%
\newcommand\listoffigures{%
    \section*{\listfigurename}%
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename}%
      \@mkboth{%
          \MakeUppercase\listtablename}%
         {\MakeUppercase\listtablename}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\def\@biblabel#1{#1.}
\newenvironment{thebibliography}[1]
     {\section*{\refname}%
      \bibliofont%
      \def\@tempa{#1}%
      %\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
    \ifx\@tempa\@empty
      \list{}%
           {\labelwidth0\p@\labelsep0\p@%
            \leftmargin12\p@\itemindent-12\p@%
            \itemsep\z@
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}%
           }%
    \else%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \itemsep\z@
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}%
            \def\@biblabel##1{\@ifundefined{bibnote@\romannumeral\theenumiv}{{\reset@font\normalfont\bibliofont##1.}}{{\bfseries##1.}}}
           }%
    \fi%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \columnseprule \z@
                \columnsep 35\p@
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width88.5\p@ height.5\p@
  \kern3\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{LITERATURE CITED}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{12\p@}
\setlength\columnseprule{0\p@}
%
\pagestyle{headings}
\sloppy
\voffset-1pc
\hoffset-2.75pc
%
\pagenumbering{arabic}
%
%%% Body environment declaration
%
\newcount\armarginnotecnt
\newbox\armarginnotebox%
\newenvironment{armarginnote}[1][\relax]
  {\def\@tempa{#1}\global\advance\armarginnotecnt\@ne\let\entry\armarginnoteentry%
   \edef\armarginnote@page{\getpagerefnumber{armarginnote@\the\armarginnotecnt}}
   \setbox\armarginnotebox\vbox\bgroup\hsize6.5pc\armarginnotefont\parindent\z@\if\@tempa\relax\else\vskip\@tempa\fi%
   \noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}\par%
  }%
  {\par\vskip-2\p@\noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}\egroup%
   \ifodd\armarginnote@page
     \marginpar{\hspace*{2.5\p@}\box\armarginnotebox}%
   \else
     \marginpar{\hspace*{-23\p@}\box\armarginnotebox}%
   \fi\label{armarginnote@\the\armarginnotecnt}}
% 
\def\armarginnoteentry#1#2{\par\addvspace{4\p@}\noindent{\fignum@color{\sffamily\bfseries#1:}}\ #2\par}
%
\newcount\bibnotecnt
\def\bibnote{\@ifnextchar[{\@bibnote}{\@bibnote[\relax]}}
\def\@bibnote[#1]#2{\def\chk@tempa{#1}%
  \global\advance\bibnotecnt\@ne%
  \edef\bibnote@page{\getpagerefnumber{bibnote@\the\bibnotecnt}}%
  \ifodd\bibnote@page
    \marginpar{\hbox{\reset@font\hskip2.5\p@\parbox{6.5pc}{\vbox to \z@{\if\chk@tempa\relax\else\vskip\chk@tempa\fi\noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}\par{\bibarmarginnotefont\ifx\chk@natbib@nonumref\YES\else\ifx\@tempa\@empty\else\theenumiv.\ \fi\fi#2\par}\par\vskip-6.6\p@%
    \noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}}}}}%
  \else%
    \marginpar{\hbox{\reset@font\hskip-23\p@\parbox{6.5pc}{\vbox to \z@{\if\chk@tempa\relax\else\vskip\chk@tempa\fi\noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}\par{\reset@font\bibarmarginnotefont\ifx\chk@natbib@nonumref\YES\else\ifx\@tempa\@empty\else\theenumiv.\ \fi\fi#2\par}\vskip-6.6\p@%
    \noindent\hbox to 6.5pc{\marginrulecolor{\hrulefill}}}}}}%
\fi%
\label{bibnote@\the\bibnotecnt}
\immediate\write\@mainaux{\string\gdef\string\bibnote@\romannumeral\theenumiv{}}\bfseries\mathversion{bold}}
%
\def\@bibitem#1{\item\if@filesw \immediate\write\@auxout
       {\string\bibcite{#1}{\the\value{\@listctr}}}\fi\ignorespaces\reset@font\normalfont\bibliofont}
%
\def\summaryhead#1{\noindent{\title@color{\sffamily\bfseries#1}\par\vskip4.2\p@}}
\newenvironment{summary}[1][\relax]%
  {\par\overfullrule0\p@\let\centerline\leftline%
   \definecolor{shadecolor}{cmyk}{0.03,0.03,0.12,0}
  \fboxsep10\p@\advance\hsize-2\fboxsep\begin{shaded}\leftmargini12\p@%
  \if#1\relax\else\summaryhead{#1}\fi}
  {\end{shaded}}
%
\def\issueshead#1{\noindent{\title@color{\sffamily\bfseries#1}\par\vskip4.2\p@}}
\newenvironment{issues}[1][\relax]%
  {\par\overfullrule0\p@\let\centerline\leftline%
   \definecolor{shadecolor}{cmyk}{0.06,0.02,0.04,0}
  \fboxsep10\p@\advance\hsize-2\fboxsep\begin{shaded}\leftmargini12\p@%
  \if#1\relax\else\issueshead{#1}\fi}
  {\end{shaded}}
%
\gdef\fps@textbox{t}
\def\ftype@textbox{3}
\def\textboxhere{h}
\def\textboxbottom{b}
\newcommand\textboxhead[1]{\noindent{\textboxheadfont\fignum@color{\sffamily#1}\par\vskip6\p@}\noindent\ignorespaces}
\newcommand\textboxsubhead[1]{\par\addvspace{12\p@}\noindent{\textboxheadfont\fignum@color{\sffamily#1}\par\vskip6\p@}\noindent\ignorespaces}
\newcommand\textboxsubsubhead[1]{\par\addvspace{12\p@}\noindent{\textboxheadfont\fignum@color{\sffamily#1.}\hskip6\p@}\noindent\ignorespaces}
\newcount\textboxcnt
\newenvironment{textbox}[1][\relax]%
  {\gdef\textboxpos{#1}\global\advance\textboxcnt\@ne%
   \edef\textbox@page{\getpagerefnumber{textbox@\the\textboxcnt}} 
   \setbox\z@\vbox\bgroup\hsize38pc\textboxfont\vskip12\p@\noindent\ignorespaces%
   \let\section\textboxhead%
   \let\subsection\textboxsubhead%
   \let\subsubsection\textboxsubsubhead%
  }
  {\vskip11\p@\egroup%
   \ifx\textboxpos\textboxbottom\@float{textbox}[b]\else\ifx\textboxpos\textboxhere\@float{textbox}[h]\else\@float{textbox}[t]\fi\fi%
   \ifodd\textbox@page\else\hskip-7.5pc\fi%
   \textboxcolor{\vrule height\ht\z@ width\wd\z@ depth\dp\z@}%
   \llap{\box\z@}
   \label{textbox@\the\textboxcnt}
   \end@float}
%
\AtEndDocument{\immediate\write\@mainaux{\string\endpage{\thepage}}}%
%
\usepackage{graphicx,multicol,refcount,framed}%colortbl,array,
%
\def\Gin@setfile#1#2#3{%
  \ifx\\#2\\\Gread@false\fi
  \ifGin@bbox\else
    \ifGread@
      \csname Gread@%
         \expandafter\ifx\csname Gread@#1\endcsname\relax
           eps%
         \else
           #1%
         \fi
      \endcsname{\Gin@base#2}%
    \else
      \Gin@nosize{#3}%
    \fi
  \fi
  \Gin@viewport@code
  \Gin@nat@height\Gin@ury bp%
  \advance\Gin@nat@height-\Gin@lly bp%
  \Gin@nat@width\Gin@urx bp%
  \advance\Gin@nat@width-\Gin@llx bp%
  \Gin@req@sizes
  \expandafter\ifx\csname Ginclude@#1\endcsname\relax
    \Gin@drafttrue
    \expandafter\ifx\csname Gread@#1\endcsname\relax
      \@latex@error{Can not include graphics of type: #1}\@ehc
      \global\expandafter\let\csname Gread@#1\endcsname\@empty
    \fi
  \fi
  \leavevmode
  \ifGin@draft
      \hb@xt@\Gin@req@width{%
        \vrule\hss
        \vbox to \Gin@req@height{%
           \hrule \@width \Gin@req@width
           \vss
           \edef\@tempa{#3}%
           \rlap{ \ttfamily\expandafter\strip@prefix\meaning\@tempa}%
           \vss
           \hrule}%
        \hss\vrule}%
  \else
    \@addtofilelist{#3}%
    \ProvidesFile{#3}[Graphic file (type #1)]%
    \setbox\z@\hbox{\csname Ginclude@#1\endcsname{#3}}%
    \dp\z@\z@
    \ht\z@\Gin@req@height
    \wd\z@\Gin@req@width
    \global\chk@fig@width\Gin@req@width
    \edef\chk@fig@temp{\getpagerefnumber{chk@figure@\the\chkfigcnt}}
    \ifdim\Gin@req@width>\typewidth%
      \centerline{\box\z@}
    \else
      \ifdim\Gin@req@width>\textwidth%
        \@tempdima\typewidth\advance\@tempdima-\Gin@req@width\divide\@tempdima by 2%
        \leftline{\ifodd\chk@fig@temp\hbox{\hspace*{\@tempdima}\box\z@}\else\hbox{\hspace*{-7.5pc}\hspace*{\@tempdima}\box\z@}\fi}
      \else%
        \centerline{\box\z@}
      \fi%  
    \fi%  
  \fi}
%
%
\def\chk@natbib@nonumref{no}%
\AtBeginDocument{\immediate\write\@mainaux{\string\fstpage{\thepage}}%
  \@ifpackageloaded{natbib}{%
    \global\let\theenumiv\theNAT@ctr%
    \ifNAT@numbers%
      \renewcommand\NAT@open{(}%
      \renewcommand\NAT@close{)}%
      \gdef\chk@natbib@nonumref{no}%
      \gdef\bibnumfmt#1{\@ifundefined{bibnote@\romannumeral\theenumiv}{{\reset@font\normalfont\bibliofont#1.}}{{\bfseries#1.}}}%
    \else%
      \gdef\chk@natbib@nonumref{yes}%
      \gdef\NAT@aysep{}%
    \fi%
    \let\bibfont\bibliofont%
    \global\bibsep\z@%
    \def\bibitem@fin{\@ifxundefined\@bibstop{}{\csname bibitem@\@bibstop\endcsname}\reset@font\normalfont\bibliofont}%
  }{\def\@cite#1#2{({#1\if@tempswa , #2\fi})}}%
%
\if@DotinEqNum%
  \@ifpackageloaded{amsmath}{%
    %\def\tagform@#1{\maketag@@@{(\ignorespaces#1\unskip\@@italiccorr)}}
    \def\tagform@#1{\maketag@@@{\ignorespaces#1\unskip.}}
  }{%
    %\def\@eqnnum{{\normalfont \normalcolor (\theequation)}}
    \def\@eqnnum{{\normalfont \normalcolor \theequation.}}
  }%
\fi% 
%%
\@ifpackageloaded{hyperref}{%
  \def\@linkcolor{green}%
  \def\termkey#1{\hyperlink{key#1}{#1}}%
  \def\armarginnoteentry#1#2{\par\addvspace{4\p@}\noindent{\fignum@color{\sffamily\bfseries#1:\hypertarget{key#1}{}}}\ #2\par}%
}{%
  \def\termkey#1{#1}%
}
%
}
%
\def\@outputpage{%
%\special{color push  cmyk 0 0 0 1.0}%
\begingroup           % the \endgroup is put in by \aftergroup
\let\firstmark\botmark
  \let \protect \noexpand
  \@resetactivechars
  \@parboxrestore
  \shipout \vbox{%
    \set@typeset@protect
    \aftergroup \endgroup
    \aftergroup \set@typeset@protect
                                % correct? or just restore by ending
                                % the group?
  \if@specialpage
    \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
  \fi
  \if@twoside
    \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
         \let\@themargin\oddsidemargin
    \else \let\@thehead\@evenhead
       \let\@thefoot\@evenfoot \let\@themargin\evensidemargin
    \fi
  \fi
  \reset@font
  \normalsize
  \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
    \@begindvi\trimmarks
    \vskip \topmargin
    \ifodd\c@page\else\advance\@themargin\extramargin\fi%
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@hbox
          \normalcolor
          \hb@xt@\textwidth {%
            \let \label \@gobble
            \let \index \@gobble
            \let \glossary \@gobble %% 21 Jun 91
            \@thehead
            }%
        \color@endbox
        }%                        %% 22 Feb 87
      \dp\@tempboxa \z@
%      \box\@tempboxa
%      \vskip \headsep
      \box\@outputbox
      \baselineskip \footskip
      \color@hbox
        \normalcolor
        \hb@xt@\textwidth{%
          \let \label \@gobble
          \let \index \@gobble      %% 22 Feb 87
          \let \glossary \@gobble   %% 21 Jun 91
          \@thefoot
          }%
      \color@endbox
      }%
    }%
\global \@colht \textheight
\stepcounter{page}%
%\special{color pop}
}
\def\trimmarks{}%
%
\ifpdf%
  \setlength{\pdfpagewidth}{8.5in}%
  \setlength{\pdfpageheight}{11in}%
\fi%
%%
\endinput
%%
%% End of file `ar.cls'.
