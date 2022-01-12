%%
%% This is file `aa.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% aa.dtx  (with options: `aa')
%% 
%% aa.dtx
%% Copyright EDP Sciences, E-press (2016), J.-C. Charpentier (2010), all rights reserved.
%% Class `aa' to use with LaTeX 2e
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
%%
\def\filedate{2016/09/01}
\let\aafiledate\filedate
\def\fileversion{9.1}
\let\aafileversion\fileversion
\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{aa}
         [\filedate{} v\fileversion{} Class for AA journal]
\LoadClass[fleqn]{article}
\newif\ifaa@letter
\DeclareOption{letter}{\aa@lettertrue}
\newif\ifaa@onlletter
\DeclareOption{onlletter}{\aa@onllettertrue}
\newif\ifaa@online
\DeclareOption{online}{\aa@onlinetrue}
\newif\ifaa@referee
\DeclareOption{referee}{\aa@refereetrue}
\newif\ifaa@contents
\DeclareOption{contents}{\aa@contentstrue}
\newif\ifaa@contentsAut
\DeclareOption{contentsAut}{\aa@contentsAuttrue}
\newlength\aa@fullwidth
\newlength\aa@fullheight
\DeclareOption{a4paper}{%
  \renewcommand*\aa@paper{a4paper}%
  \setlength\aa@fullwidth{184mm}%
  \renewcommand*\aa@leftmarginpaper{(21cm-\aa@fullwidth)/2}%
  \setlength\aa@fullheight{705pt}%
}
\DeclareOption{letterpaper}{%
  \renewcommand*\aa@paper{letterpaper}
  \setlength\aa@fullwidth{184mm}%
  \renewcommand*\aa@leftmarginpaper{(8.5in-\aa@fullwidth)/2}%
  \setlength\aa@fullheight{705pt-1.76cm}%
}
\DeclareOption{legalpaper}{%
  \renewcommand*\aa@paper{legalpaper}%
  \setlength\aa@fullwidth{184mm}%
  \renewcommand*\aa@leftmarginpaper{(8.5in-\aa@fullwidth)/2}%
  \setlength\aa@fullheight{705pt}%
}
\newcommand*\aa@paper{}
\newcommand*\aa@leftmarginpaper{}
\newif\ifaa@twocolumn
\DeclareOption{onecolumn}{%
  \aa@twocolumnfalse
  \renewcommand*\aa@column{onecolumn}%
}
\DeclareOption{twocolumn}{%
  \aa@twocolumntrue
  \renewcommand*\aa@column{twocolumn}%
}
\newcommand*\aa@column{}
\newif\if@wideboxfn
\@wideboxfnfalse
\DeclareOption{wideboxfn}{\@wideboxfntrue}
%
\newif\ifaa@longfn
\DeclareOption{longfn}{\aa@longfntrue}
\newif\ifaa@longauth
\DeclareOption{longauth}{\aa@longauthtrue}
\let\if@myobjects\iffalse
\DeclareOption{ascii}{%
  \def\aa@inputenc{\RequirePackage[ascii]{inputenc}}%
}
\DeclareOption{latin1}{%
  \def\aa@inputenc{\RequirePackage[latin1]{inputenc}}%
}
\DeclareOption{latin9}{%
  \def\aa@inputenc{\RequirePackage[latin9]{inputenc}}%
}
\DeclareOption{ansinews}{%
  \def\aa@inputenc{\RequirePackage[ansinews]{inputenc}}%
}
\DeclareOption{applemac}{%
  \def\aa@inputenc{\RequirePackage[applemac]{inputenc}}%
}
\DeclareOption{utf8}{%
  \def\aa@inputenc{\RequirePackage[utf8]{inputenc}}%
}
\newcommand*\aa@font{}
\DeclareOption{cm}{%
  \def\aa@font{%
    \RequirePackage{textcomp}%
  }%
}
\DeclareOption{showoverfull}{\setlength\overfullrule{5\p@}}
\DeclareOption{hideoverfull}{\setlength\overfullrule{\z@}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
\newif\ifaa@bibnumber
\DeclareOption{bibnumber}{\aa@bibnumbertrue%
  \def\aa@biblio{%
    \RequirePackage[sort&compress]{natbib}
    \bibpunct{[}{]}{,}{n}{}{,}%
    \def\@biblabel##1{##1.}%
    \bibliographystyle{plainnat}%
    \def\bibfont{\aa@bibliographyfont}%
    \setlength{\bibsep}{\z@}%
  }
}
\DeclareOption{bibauthoryear}{%
  \def\aa@biblio{%
    \RequirePackage{natbib}
    \bibpunct{(}{)}{;}{a}{}{,}%
    \def\bibfont{\aa@bibliographyfont}%
    \setlength{\bibsep}{\z@}%
  }%
}
\DeclareOption{bibyear}{%
  \def\aa@biblio{%
    \def\@biblabel##1{}%
    \def\@cite##1##2{{##1\if@tempswa , ##2\fi}}%
  }
}
\ExecuteOptions{a4paper,twocolumn,utf8,hideoverfull,bibauthoryear}
\ProcessOptions
\RequirePackage{calc}
\RequirePackage{etex}
\RequirePackage[T1]{fontenc}
\aa@inputenc
\aa@font
\ifaa@referee
\RequirePackage{xcolor}
\RequirePackage{marginnote}
\AtEndOfClass{
\setlength{\marginparwidth}{50mm}
\setlength{\marginparsep}{5mm}
}
%
\newcounter{lenote}
\AtBeginDocument{%
\@ifpackageloaded{hyperref}%
{\newcommand{\LEm}[1]{\texorpdfstring{\protect\stepcounter{lenote}{\normalfont\textbf{\color{red}[Note \arabic{lenote}] }}
\protect\marginnote{\baselineskip3ex{\normalsize{\normalfont\textbf{\color{red}[Note \arabic{lenote}: #1]}}}}}{}}
\newcommand{\LEt}[1]{\texorpdfstring{\protect\stepcounter{lenote}\normalfont\textbf{\color{red}[Note \arabic{lenote}: #1] }}{}}}%
{\newcommand{\LEm}[1]{\protect\stepcounter{lenote}{\normalfont\textbf{\color{red}[Note \arabic{lenote}] }}
\protect\marginnote{\baselineskip3ex{\normalsize{\normalfont\textbf{\color{red}[Note \arabic{lenote}: #1]}}}}}
\newcommand{\LEt}[1]{\protect\stepcounter{lenote}\normalfont\textbf{\color{red}[Note \arabic{lenote}: #1] }}}%
}
  \renewcommand*\aa@column{onecolumn}
  \aa@twocolumnfalse
  \setlength{\aa@fullwidth}{0.75\aa@fullwidth}
  \renewcommand*\aa@leftmarginpaper{15mm}
  \if!\csname appendix\endcsname!\AtEndDocument{\listofobjects}\fi
\fi
\ifx\pdfoutput\undefined
  \newcount\pdfoutput
\fi
\ifnum\pdfoutput=\z@
  \newcommand\aa@driver{dvips}
\else
  \newcommand\aa@driver{pdftex}
\fi
\RequirePackage[\aa@driver, \aa@paper, left=\aa@leftmarginpaper,
                \ifaa@twocolumn\aa@column\fi, \ifaa@referee\else twoside,\fi
                textheight=\aa@fullheight,headheight=12pt,headsep=12pt,
                footskip=20pt,
                textwidth=\aa@fullwidth, columnsep=4mm,
                marginpar=0in, marginparsep=0in]{geometry}
\RequirePackage{setspace}
\ifaa@referee
  \doublespacing
\fi
\RequirePackage{fancyhdr}
\RequirePackage{cuted}
\RequirePackage[tbtags,fleqn]{amsmath}
\RequirePackage{amssymb}
\setlength{\mathindent}{\z@}
\ifaa@letter
  \RequirePackage{xcolor}
\fi
\RequirePackage{url, enumerate}
\RequirePackage{longtable}
\renewcommand*\LT@makecaption[3]{%
  \LT@mcol\LT@cols
  c{\hbox to\z@{%
      \hss\parbox[t]\LTcapwidth{%
        \aa@longtablecaptionfont
        \sbox\@tempboxa{%
          #1{{\aa@longtableheadcaptionfont#2.} }#3%
        }%
        \ifdim\wd\@tempboxa>\hsize
          #1{{\aa@longtableheadcaptionfont#2.} }#3%
        \else
          \hbox to\hsize{\box\@tempboxa\hfil}%
        \fi
        \endgraf\vskip\baselineskip
      }%
      \hss
    }%
  }%
}
\aa@biblio
\RequirePackage{xspace}
\RequirePackage[english]{babel}
\def\aatolm#1{\gdef\@aatolm{#1}}
\def\@aatolm{}
\renewenvironment{appendix}{\par\clearpage\aatolm{1}
\def\@hangfrom##1{\setbox\@tempboxa\hbox{{##1}}%
      \hangindent \parindent\noindent\box\@tempboxa}
\def\sectcounterend{:}
\def\@seccntformat##1{\csname theap##1\endcsname%
                      \hskip\betweenumberspace}
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}%
  \setcounter{table}{0}%
  \newcommand\theapsection{\appendixname\ \@Alph\c@section\sectcounterend}%
  \newcommand\theapsubsection{\@Alph\c@section.\@arabic\c@subsection.}%
  \newcommand\theapsubsubsection{\theapsubsection\@arabic\c@subsubsection.}%
  \renewcommand\thesection{\@Alph\c@section}%
  \renewcommand\thesubsection{\@Alph\c@section.\@arabic\c@subsection}%
  \renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}%
  \renewcommand\theequation{\@Alph\c@section.\@arabic\c@equation}%
  \renewcommand\thefigure{\@Alph\c@section.\@arabic\c@figure}%
  \renewcommand\thetable{\@Alph\c@section.\@arabic\c@table}%
  \@addtoreset{equation}{section}%
  \@addtoreset{figure}{section}%
  \@addtoreset{table}{section}}
{
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}%
  \setcounter{table}{0}%
  \renewcommand\thesection{\@arabic\c@section}%
  \renewcommand\thesubsection{\@arabic\c@subsection}%
  \renewcommand\thesubsubsection{\@arabic\c@subsubsection}%
  \renewcommand\theequation{\@arabic\c@equation}%
  \renewcommand\thefigure{\@arabic\c@figure}%
  \renewcommand\thetable{\@arabic\c@table}%
}
\RequirePackage[figuresright]{rotating}
\newcount\aa@cnta
\newcount\aa@cntb
\newcount\aa@nbwords
\newcount\aa@nbauthors
\newcount\aa@nbinstitutes
\newcommand*\aa@traditabstractmax{300}
\newcommand*\aa@structabstractmax{300}
\newcommand*\aa@keywordsmax{6}
\newcounter{aa@institutecnt}
\newcount\aa@savecounter
\newlength{\aa@dima}
\newlength{\aa@dimb}
\newlength{\aa@dimc}
\newcommand*\aa@removefootnoteindentspace{-1.5em}
\newcommand*\aa@removefootnotewbindentspace{-2.25em}
\newcommand*\aa@logoruleindent{8\p@}
\newcommand*\aa@logoruleraise{4\p@}
\newcommand*\aa@logorulelength{78\p@}
\newcommand*\aa@logorulethick{1\p@}
\newcommand*\aa@logoraise{25\p@}
\newcommand*\aa@logoraiserest{14.5\p@}
\newlength{\VSpaceBeforeTheorem}
\setlength{\VSpaceBeforeTheorem}{8\p@ \@plus 2\p@ \@minus 2\p@}
\newlength{\VSpaceAfterTheorem}
\setlength{\VSpaceAfterTheorem}{8\p@ \@plus 2\p@ \@minus 2\p@}
\newlength{\HSpaceIndentTheorem}
\setlength{\HSpaceIndentTheorem}{\z@}
\newlength{\HVSpaceSepTheorem}
\setlength{\HVSpaceSepTheorem}{0.25em}
\newlength{\aa@securityheadmargin}
\setlength{\aa@securityheadmargin}{2cm}
\newlength{\aa@bleedmargin}
\setlength{\aa@bleedmargin}{5mm}
\newlength{\aa@voffsetlettername}
\setlength{\aa@voffsetlettername}{2.55cm}
\newlength{\aa@hoffsetlettername}
\setlength{\aa@hoffsetlettername}{2mm}
\newlength{\aa@hoffsetgray}
\setlength{\aa@hoffsetgray}{3mm}
\newlength{\aa@aftercontentsnameskip}
\setlength{\aa@aftercontentsnameskip}{12\p@}
\newlength{\aa@aftercontentsAutnameskip}
\setlength{\aa@aftercontentsAutnameskip}{12\p@}
\newlength{\aa@beforeheadboxskip}
\setlength{\aa@beforeheadboxskip}{-20\p@}
\newlength{\aa@afterheadboxskip}
\setlength{\aa@afterheadboxskip}{49\p@}
\newlength{\aa@beforetopicsskip}
\setlength{\aa@beforetopicsskip}{18\p@}
\newlength{\aa@aftertopicsskip}
\setlength{\aa@aftertopicsskip}{32\p@}
\newlength{\aa@afterletternameskip}
\setlength{\aa@afterletternameskip}{6mm}
\newlength{\aa@beforeletternameskip}
\setlength{\aa@beforeletternameskip}{-4.8mm}
\newlength{\aa@beforeheadnoteskip}
\setlength{\aa@beforeheadnoteskip}{-12\p@}
\newlength{\aa@afterheadnoteskip}
\setlength{\aa@afterheadnoteskip}{3.5mm}
\newlength{\aa@aftertitleskip}
\setlength{\aa@aftertitleskip}{10\p@}
\newlength{\aa@beforesubtitleskip}
\setlength{\aa@beforesubtitleskip}{7\p@}
\newlength{\aa@aftersubtitleskip}
\setlength{\aa@aftersubtitleskip}{10\p@}
\newlength{\aa@afterauthorskip}
\setlength{\aa@afterauthorskip}{16.5\p@}
\newlength{\aa@beforedefaultheadnoteskip}
\setlength{\aa@beforedefaultheadnoteskip}{2pt}
\newlength{\aa@headingleftmargin}
\setlength{\aa@headingleftmargin}{8mm}
\newlength{\aa@headingrightmargin}
\setlength{\aa@headingrightmargin}{8mm}
\newlength{\aa@afterinstituteskip}
\setlength{\aa@afterinstituteskip}{8\p@}
\newlength{\aa@instituteindent}
\newlength{\aa@afterdateskip}
\setlength{\aa@afterdateskip}{16.5\p@}
\newlength{\aa@afterdedicationskip}
\setlength{\aa@afterdedicationskip}{16.5\p@}
\newlength{\aa@afterabstractnameskip}
\setlength{\aa@afterabstractnameskip}{9\p@}
\newlength{\aa@afterabstractskip}
\setlength{\aa@afterabstractskip}{6\p@}
\newlength{\aa@aftersubabstractskip}
\setlength{\aa@aftersubabstractskip}{\z@}
\newlength{\aa@afterkeywordsskip}
\setlength{\aa@afterkeywordsskip}{24\p@}
\newlength{\aa@beforeacknowledgementsskip}
\setlength{\aa@beforeacknowledgementsskip}{6pt}
\newlength{\aa@afteracknowledgementsskip}
\setlength{\aa@afteracknowledgementsskip}{6pt}
\newlength{\aa@figgap}
\setlength{\aa@figgap}{12pt}
\newlength{\aa@sidecaptionminwidth}
\setlength{\aa@sidecaptionminwidth}{5cm}
\newlength{\aa@sidecaptionwidth}
\newlength{\VSpaceBeforeFigCaption}
\setlength{\VSpaceBeforeFigCaption}{5pt}
\newlength{\VSpaceAfterTabCaption}
\setlength{\VSpaceAfterTabCaption}{10pt}
\newlength{\VSpaceBeforeTabBib}
\setlength{\VSpaceBeforeTabBib}{2ex}
\newlength{\VSpaceBeforeTabFoot}
\setlength{\VSpaceBeforeTabFoot}{2ex}
\newlength{\aftertext}
\setlength{\aftertext}{5pt}
\newlength{\aa@entryboxwidth}
\setlength{\aa@entryboxwidth}{8mm}
\newlength{\betweenumberspace}
\setlength{\betweenumberspace}{3.33pt}
\newif\ifaa@structabstract
\newif\ifaa@numberedtheorem
\newif\ifaa@firstonline
\aa@firstonlinetrue
\newcommand*\journalname[1]{\renewcommand*\aa@journalname{#1}}
\newcommand*\aa@journalname{Astronomy \& Astrophysics}
\newcommand*\aa@manuscriptname{%
  manuscript no. \jobname
  \hspace{\stretch{1}}%
  \copyright ESO \the\year
}
\def\Num{1}
\def\Vol{1}
\def\Month{JANUARY~I}
\newcommand*\AANum[1]{\renewcommand*\aa@numarticle{#1}%
\gdef\@AANum{#1}}
\newcommand*\@AANum{}
\newcommand*\aa@numarticle{Article number}
\newcommand*\doi[1]{%
  \renewcommand*\aa@doi{DOI: #1}%
  \renewcommand*\aa@doifig{#1}%
}
\newcommand*\aa@doi{(DOI: will be inserted by hand later)}
\newcommand*\aa@doifig{(to add)}
\newcommand*\AALogo{%
  \settowidth{\aa@dima}{\aa@logo@astrophysicsfont Astrophysics}%
  \setbox0=\hbox{%
    \begin{minipage}{\aa@dima}
      \aa@logo@astronomyfont Astronomy\\
      \raisebox{\aa@logoraise}{%
        \parbox[c][0pt]{\z@}{%
          \rlap{%
            \hspace*{\aa@logoruleindent}%
            \aa@logo@andfont \&%
            \rule[\aa@logoruleraise]{\aa@logorulelength}{\aa@logorulethick}%
          }%
        }%
      }%
      \raisebox{\aa@logoraiserest}{\aa@logo@astrophysicsfont Astrophysics}
    \end{minipage}%
  }%
  \copy0%
}
\newcommand*\aa@pageof{, page \thepage{} of \pageref{LastPage}}
\newcommand*{\aa@publishlink}{\hfil%
 Article published by EDP Sciences%
 \hfil}
\newcommand*\aa@topicsname{Special feature}
\newcommand*\aa@offprintsname{Send offprint requests to\/}
\newcommand*\listobjectname{List of Objects}
\addto\captionsenglish{\renewcommand*\figurename{Fig.}}
\newcommand*\aa@figtype{figure}
\newcommand*\aa@tabtype{table}
\newcommand*\lettername{Letter to the Editor}
\newcommand*\noteaddname{Note added in proof.}
\def\noteadd{\par\addvspace{17pt}\rm
\trivlist\item[\hskip\labelsep
{\it\noteaddname}]}
\def\endnoteadd{\endtrivlist\addvspace{6pt}}
\newcommand*\andname{\unskip, }
\newcommand*\andlastname{\unskip, and }
\newcommand*\andtwoname{\unskip\ and }
\newcommand*\etalname{et al.}
\newcommand*\aa@contentsname{Contents sorted by section}
\newcommand*\aa@contentsAutname{Contents sorted by first author}
\newcommand*\aa@pagename{page}
\newcommand*\aa@AAname{A\&A}
\newcommand*\aa@copyrightname{\copyright~ESO}
\newcommand*\aa@subabstractiname{Context}
\newcommand*\aa@subabstractiiname{Aims}
\newcommand*\aa@subabstractiiiname{Methods}
\newcommand*\aa@subabstractivname{Results}
\newcommand*\aa@subabstractvname{Conclusions}
\newcommand*\aa@keywordname{Key words}
\newcommand*\aa@acknowledgementsname{Acknowledgements}
\newcommand*\aa@floatcounterend{.~}
\newcommand*\aa@tablebibname{References}
\newcommand*\aa@tablefootname{Notes}
\renewcommand\normalsize{%
  \@setfontsize\normalsize\@xpt{11}%
  \abovedisplayskip 8.55\p@ \@plus4\p@ \@minus5\p@
  \abovedisplayshortskip \z@ \@plus3\p@
  \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
  \belowdisplayskip \abovedisplayskip
  \let\@listi\@listI}
\normalsize
\renewcommand\small{%
  \@setfontsize\small\@ixpt{10}%
  \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
  \abovedisplayshortskip \z@ \@plus2\p@
  \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
  \def\@listi{\leftmargin\leftmargini
    \parsep 0\p@ \@plus1\p@ \@minus1\p@
    \topsep 8\p@ \@plus2\p@ \@minus4\p@
    \itemsep0\p@}%
  \belowdisplayskip \abovedisplayskip
}
\let\footnotesize\small
\renewcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\renewcommand\tiny{\@setfontsize\tiny\@viiipt\@ixpt}
\renewcommand\large{\@setfontsize\large\@xipt{13}}
\renewcommand\Large{\@setfontsize\Large\@xivpt{16}}
\renewcommand\LARGE{\@setfontsize\LARGE\@xviipt{20}}
\renewcommand\huge{\@setfontsize\huge\@xxpt{24}}
\renewcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\newcommand*\aa@headfont{\small}
\newcommand*\aa@footfont{\small}
\newcommand*\aa@subabstractfont{\small\sffamily\itshape}
\newcommand*\aa@abstractfont{\small}
\newcommand*\aa@offprintsnamefont{\itshape}
\newcommand*\aa@logo@astronomyfont{\fontsize{18}{21}\sffamily\bfseries}
\newcommand*\aa@logo@andfont{\fontsize{16}{12}\sffamily\mdseries}
\newcommand*\aa@logo@astrophysicsfont{\fontsize{18}{21}\sffamily\bfseries}
\newcommand*\aa@letternamefont{\sffamily\fontsize{19.2\p@}{21\p@}\bfseries}
\newcommand*\aa@headletternamefont{\large\scshape}
\newcommand*\aa@contentsnamefont{\Large\sffamily\bfseries}
\newcommand*\aa@contentsAutnamefont{\Large\sffamily\bfseries}
\newcommand*\aa@journalnamefont{\slshape}
\newcommand*\aa@manuscriptnamefont{}
\newcommand*\aa@todayheadboxfont{}
\newcommand*\aa@msnrheadboxfont{\LARGE\bfseries\smash}
\newcommand*\aa@textidlineempty{{\slshape A\&A proofs:}\ manuscript no.~\jobname}
\newcommand*\aa@idlinefont{}
\newcommand*\aa@doifont{}
\newcommand*\aa@copyrightfont{\small}
\ifaa@letter
  \definecolor{aa@colorletter}{gray}{0.6}
\fi
\newcommand*\aa@topicsfont{\fontsize{12}{14}\sffamily\itshape}
\newcommand*\aa@topicsnamefont{\fontsize{17}{20}\sffamily}
\newcommand*\aa@headnotefont{\LARGE\itshape}
\newcommand*\aa@defaultheadnotefont{\LARGE\sffamily\bfseries\itshape}
\newcommand*\aa@titlefont{\LARGE \sffamily\bfseries\boldmath}
\newcommand*\aa@subtitlefont{\Large \sffamily\bfseries\boldmath}
\newcommand*\aa@authorfont{\large}
\newcommand*\aa@institutefont{\small}
\newcommand*\aa@datefont{\small}
\newcommand*\aa@dedicationfont{\itshape}
\newcommand*\aa@abstractnamefont{\small\sffamily\bfseries\MakeUppercase}
\newcommand*\aa@keywordnamefont{\small\sffamily\bfseries}
\newcommand*\aa@keywordfont{\small}
\newcommand*\aa@longauthfont{\centering\small\itshape}
\newcommand*\aa@longauthname{Affiliations can be found after the references}
\newcommand*\aa@longauthopen{(}
\newcommand*\aa@longauthclose{)}
\newcommand*\aa@acknowledgementsfont{\tiny}
\newcommand*\aa@acknowledgementsnamefont{\itshape}
\newcommand*\aa@bibliographyfont{\tiny}
\newcommand*\aa@captionfont{\normalfont\small}
\newcommand*\aa@headcaptionfont{\bfseries}
\newcommand*\aa@tablebibfont{\small}
\newcommand*\aa@tablebibnamefont{\small\bfseries}
\newcommand*\aa@tablefootfont{\small}
\newcommand*\aa@tablefootnamefont{\small\bfseries}
\newcommand*\aa@descriptionlabelfont{}
\newcommand\aa@longtablecaptionfont{\small}
\newcommand\aa@longtableheadcaptionfont{\small\bfseries}
\newcommand*\aa@emailfont{\ttfamily}
\newcommand*\aa@entryifont{\itshape}
\newcommand*\aa@entryiifont{\bfseries}
\newtoks\authorlink
\newtoks\yearCop
\newbox\aa@stripbox
\newbox\aa@outputbox
\newenvironment{aa@strip}
{%
  \onecolumn\parindent0mm
  \setbox\aa@stripbox=\vbox\bgroup
}
{%
  \egroup
  \aa@@strip
}
\def\aa@@strip{%
  \ifdim\ht\aa@stripbox>\textheight
    \global\setbox\aa@outputbox=\vsplit\aa@stripbox to \textheight
    \dp\aa@outputbox \z@
    \box\aa@outputbox\clearpage
    \expandafter\aa@@strip
  \else
    \dp\aa@stripbox \z@
    \twocolumn[\box\aa@stripbox]
  \fi
}
\newcommand{\ErratumOf}[3][]{}
\newcommand{\NoteFor}[3][]{}
\renewcommand*\@fnsymbol[1]{%
  \ensuremath{%
    \ifcase#1
    \or\star
    \or{\star}{\star}
    \or{\star}{\star}{\star}%
    \or{\star}{\star}{\star}{\star}
    \or \dagger
    \or \ddagger
    \or \mathsection
    \or \mathparagraph
    \or \|
    \or **
    \or \dagger\dagger
    \or \ddagger\ddagger
    \else \@ctrerr
    \fi
  }%
}
\newcommand\longempty{}
\newcommand*\aa@cntwords[1]{%
  \protected@edef\aa@cntwords@arg{#1\space\noexpand\aa@nil\space}%
  \aa@nbwords \z@
  \let\@@protect\protect
  \let\protect\@unexpandable@protect
  \expandafter\aa@@cntwords\aa@cntwords@arg
  \restore@protect
}
\newcommand*\aa@nnil{\aa@nil}
\long\def\aa@@cntwords#1 {%
  \def\aa@arg{#1}%
  \ifx\aa@arg\aa@nnil
  \else
    \advance\aa@nbwords \@ne
    \expandafter\aa@@cntwords
  \fi
}
\newcommand*\aa@roheadletter{%
  \setlength{\aa@dima}
            {1in + \topmargin + \headheight - \dp\strutbox +
              2\aa@bleedmargin}%
  \setlength{\aa@dimb}
            {\paperwidth - (1in + \oddsidemargin + \textwidth) +
              \aa@bleedmargin - \aa@hoffsetgray}%
  \setbox0=\hbox{\aa@letternamefont \lettername}%
  \setlength{\aa@dimc}
            {\paperheight - \aa@dima - \wd0 + 2\aa@bleedmargin}%
  \addtolength{\aa@voffsetlettername}{0.5cm}
  \leavevmode
  \hbox to0pt{%
    \hspace*{\aa@hoffsetgray}%
    \vbox to0pt{%
      \hsize=\aa@dimb
      \vskip-\aa@dima
      \colorbox{aa@colorletter}{%
        \vtop{%
          \vspace*{\aa@voffsetlettername}%
          \hspace*{\aa@hoffsetlettername}%
          \rotatebox{90}{\aa@letternamefont \lettername}%
          \hfill\null
          \vspace*{\aa@dimc}%
        }%
      }%
      \vss
    }%
    \hss
  }%
}
\newcommand*\aa@leheadletter{%
  \setlength{\aa@dima}
            {1in + \topmargin + \headheight - \dp\strutbox +
              2\aa@bleedmargin}%
  \setlength{\aa@dimb}
            {1in + \oddsidemargin + \aa@bleedmargin - \aa@hoffsetgray}%
  \setbox0=\hbox{\aa@letternamefont \lettername}%
  \setlength{\aa@dimc}
            {\paperheight-\aa@dima-\wd0+2\aa@bleedmargin}%
  \addtolength{\aa@voffsetlettername}{0.5cm}
  \leavevmode
  \hspace*{-2\fboxsep}%
  \hbox to0pt{%
    \hss
    \vbox to0pt{%
      \hsize=\aa@dimb
      \vskip-\aa@dima
      \colorbox{aa@colorletter}{%
        \vtop{%
          \vspace*{\aa@voffsetlettername}%
          \hfill
          \rotatebox{90}{\aa@letternamefont \lettername}%
          \hspace*{\aa@hoffsetlettername}%
          \vspace*{\aa@dimc}%
        }%
      }%
      \vss
    }%
    \hspace*{\aa@hoffsetgray}%
  }%
}
\ifaa@online
  \fancypagestyle{firstpage}{%
    \fancyhf{}%
    \renewcommand*{\headrulewidth}{\z@}%
    \renewcommand*{\footrulewidth}{\z@}%
    \fancyfoot[RO]{\aa@footfont\aa@publishlink\aa@numarticle\aa@pageof}%
    \fancyfoot[LE]{\aa@footfont\aa@numarticle\aa@pageof\qquad\aa@publishlink}%
  }
\else
  \fancypagestyle{firstpage}{%
    \fancyhf{}%
    \renewcommand*{\headrulewidth}{\z@}%
    \renewcommand*{\footrulewidth}{\z@}%
    \fancyfoot[RO]{\aa@footfont \aa@numarticle\aa@pageof}%
    \fancyfoot[LE]{\aa@footfont \aa@numarticle\aa@pageof}%
  }
\fi
\ifaa@contents
  \fancypagestyle{otherpage}{%
    \fancyhf{}%
    \fancyhead[C]{\aa@AAname}%
    \fancyfoot[C]{\aa@contentsname, \aa@pagename\space \thepage}%
    \renewcommand*{\headrulewidth}{\z@}%
    \renewcommand*{\footrulewidth}{\z@}%
  }
\else\ifaa@contentsAut
  \fancypagestyle{otherpage}{%
    \fancyhf{}%
    \fancyhead[C]{\aa@AAname}%
    \fancyfoot[C]{\aa@contentsAutname, \aa@pagename\space \thepage}%
    \renewcommand*{\headrulewidth}{\z@}%
    \renewcommand*{\footrulewidth}{\z@}%
  }
\else
  \fancypagestyle{otherpage}{%
    \fancyhf{}%
    \fancyhead[CO]{\aa@headfont \aa@headings}%
    \fancyhead[CE]{\aa@headfont \ifx\aa@idline\empty\aa@textidlineempty\else\aa@idline\fi}%
    \fancyfoot[RO]{\aa@footfont \aa@numarticle\aa@pageof}%
    \fancyfoot[LE]{\aa@footfont \aa@numarticle\aa@pageof}%
    \renewcommand*{\headrulewidth}{\z@}%
    \renewcommand*{\footrulewidth}{\z@}%
  }
\fi\fi
\pagestyle{otherpage}
\newcommand*\aa@headings{}
\renewcommand*\date[1]{\renewcommand*\aa@date{#1}}%
\newcommand*\aa@date{\today}
\renewcommand*\title[1]{\renewcommand*\aa@title{#1}}%
\newcommand*\aa@title{}
\newcommand*\subtitle[1]{\renewcommand*\aa@subtitle{#1}}
\newcommand*\aa@subtitle{}
\renewcommand*\author[1]{%
  \renewcommand*\aa@author{#1}%
  \begingroup
    \let\and\relax
    \renewcommand*\inst[1]{}%
    \renewcommand*\thanks[1]{}%
    \renewcommand*\email[1]{}%
    \let\@@protect\protect
    \let\protect\@unexpandable@protect
    \global\aa@nbauthors \z@
    \expandafter\aa@cntauthors\aa@author\and\aa@nil\and
    \restore@protect
  \endgroup
}%
\newcommand*\aa@author{}
\def\aa@cntauthors#1\and{%
  \def\aa@arg{#1}%
  \ifx\aa@arg\aa@nnil
  \else
    \global\advance\aa@nbauthors \@ne
    \expandafter\aa@cntauthors
  \fi
}
\newcommand*\institute[1]{%
  \renewcommand*\aa@institute{#1}
  \begingroup
    \let\and\relax
    \renewcommand*\inst[1]{}%
    \renewcommand*\thanks[1]{}%
    \renewcommand*\email[1]{}%
    \let\@@protect\protect
    \let\protect\@unexpandable@protect
    \global\aa@nbinstitutes \z@
    \expandafter\aa@cntinstitutes\aa@institute\and\aa@nil\and
    \restore@protect
  \endgroup
}%
\newcommand*\aa@institute{}
\def\aa@cntinstitutes#1\and{%
  \def\aa@arg{#1}%
  \ifx\aa@arg\aa@nnil
  \else
    \global\advance\aa@nbinstitutes \@ne
    \expandafter\aa@cntinstitutes
  \fi
}
\newcommand*\inst[1]{\unskip\hbox{\@textsuperscript{\normalfont$#1$}}}
\newcommand*\fnmsep{\unskip\hbox{\@textsuperscript{\normalfont,}}}
\newcommand*\email[1]{e-mail: {\aa@emailfont #1}}
\let\endabstract\relax
\renewcommand\abstract[1]{%
  \renewcommand\aa@abstract@i{#1}%
  \@ifnextchar\bgroup{\aa@abstract@next}{\aa@abstract@one}
}
\newcommand\aa@abstract@next[4]{%
  \aa@structabstracttrue
  \renewcommand\aa@abstract@ii{#1}%
  \renewcommand\aa@abstract@iii{#2}%
  \renewcommand\aa@abstract@iv{#3}%
  \renewcommand\aa@abstract@v{#4}%
  \begingroup
    \renewcommand\thanks[1]{}%
    \protected@xdef\aa@abstractwords{%
      \aa@abstract@i \space
      \aa@abstract@ii \space
      \aa@abstract@iii \space
      \aa@abstract@iv \space
      \aa@abstract@v \space
    }%
  \endgroup
  \aa@cnta \z@
  \aa@cntb \aa@structabstractmax
  \aa@abstract@cnt
}
\newcommand*\aa@abstract@one{%
  \protected@xdef\aa@abstractwords{\aa@abstract@i}%
  \aa@cnta \z@
  \aa@cntb \aa@traditabstractmax
  \aa@abstract@cnt
}
\def\aa@abstract@cnt{%
  \aa@cntwords\aa@abstractwords
  \advance\aa@nbwords \aa@cnta\relax
  \ifnum\aa@nbwords>\aa@cntb
  \typeout{^^J^^JIMPORTANT}
    \ClassWarning{aa}
      {Around \the\aa@nbwords\space words in the abstract instead
        of \the\aa@cntb\space (max). You should rewrite your abstract}
  \typeout{IMPORTANT
  ^^J}
  \fi
}
\newcommand\aa@abstract@i{}
\newcommand\aa@abstract@ii{}
\newcommand\aa@abstract@iii{}
\newcommand\aa@abstract@iv{}
\newcommand\aa@abstract@v{}
\newcommand*\keywords[1]{%
  \renewcommand*\aa@keywords{#1}%
  \aa@keywords@cnt
  \ifnum\aa@nbwords>\aa@keywordsmax
    \ClassWarning{aa}
      {\the\aa@nbwords\space keywords instead of a maximum of
        \aa@keywordsmax\space}
  \fi
}
\newcommand*\aa@keywords{}
\newcommand*\aa@doublehyphen{--}
\newcommand*\aa@keywords@cnt{%
  \protected@edef\aa@@keywords{\aa@keywords--\noexpand\aa@nil--}%
  \aa@nbwords \z@
  \let\@@protect\protect
  \let\protect\@unexpandable@protect
  \expandafter\aa@@keywords@cnt\aa@@keywords
  \restore@protect
}
\def\aa@@keywords@cnt#1--{%
  \def\aa@arg{#1}%
  \ifx\aa@arg\aa@nnil
  \else
    \advance\aa@nbwords \@ne
    \expandafter\aa@@keywords@cnt
  \fi
}
\newcommand*\topics[1]{\renewcommand*\aa@topics{#1}}
\newcommand*\aa@topics{}
\newcommand*\titlerunning[1]{%
  \renewcommand*\aa@titlerunning{#1}
}
\newcommand*\aa@titlerunning{}
\newcommand*\authorrunning[1]{%
  \renewcommand*\aa@authorrunning{#1}
}
\newcommand*\aa@authorrunning{}
\newcommand*\offprints[1]{%
  \protected@edef\@thanks{%
    \@thanks
    \protect\footnotetext[0]{%
      \unskip%
      \if@wideboxfn%
        \hspace{\aa@removefootnotewbindentspace}%
      \else%
        \hspace{\aa@removefootnoteindentspace}%
      \fi
      {%
        \aa@offprintsnamefont
        \aa@offprintsname
      }: \ignorespaces#1%
    }%
  }%
}
\def\aainstitutename{\par
  \begingroup
  \renewcommand*\and{%
      \par\refstepcounter{aa@institutecnt}%
      \hangindent\aa@instituteindent
      \noindent
      \hbox to\aa@instituteindent{%
        \hss
        \hbox{\@textsuperscript{\normalfont\theaa@institutecnt}}%
        \enspace
      }%
      \ignorespaces
    }%
    \ifnum\aa@nbinstitutes>99
      \setbox0=\hbox{\hbox{\@textsuperscript{\normalfont888}}\enspace}%
    \else\ifnum\aa@nbinstitutes>9
      \setbox0=\hbox{\hbox{\@textsuperscript{\normalfont88}}\enspace}%
    \else
      \setbox0=\hbox{\hbox{\@textsuperscript{\normalfont8}}\enspace}%
    \fi\fi
    \aa@instituteindent \wd0
    {%
      \aa@institutefont
      \par\refstepcounter{aa@institutecnt}%
      \ifnum\aa@nbinstitutes=1\relax
      \noindent\ignorespaces
      \else
        \hangindent\aa@instituteindent
        \noindent
        \hbox to\aa@instituteindent{%
          \hss
          \hbox{\@textsuperscript{\normalfont\theaa@institutecnt}}%
          \enspace
        }\ignorespaces%
      \fi%
      \aa@institute\par
    }
  \endgroup}
\newcommand*\mail[1]{\renewcommand*\aa@mail{#1}}
\newcommand*\aa@mail{}
\newcommand*\msnr[1]{\renewcommand*\aa@msnr{#1}}
\newcommand*\aa@msnr{}
\newcommand*\idline[1]{\renewcommand*\aa@idline{#1}}
\newcommand*\aa@idline{}
\newcommand*\headnote[1]{\renewcommand*\aa@headnote{#1}}
\newcommand*\aa@headnote{}
\newcommand*\aa@defaultheadnote{Corrigendum}
\newcommand*\dedication[1]{\renewcommand*\aa@dedication{#1}}
\newcommand*\aa@dedication{}
\renewcommand*\maketitle{%
  \thispagestyle{firstpage}
\begingroup
    \if@wideboxfn
    \setlength\bibindent{1.4\parindent}
    \else
    \setlength\bibindent{\parindent}
    \fi
    \renewcommand*\thefootnote{\@fnsymbol\c@footnote}%
    \renewcommand\@makefntext[1]{%
    \ifaa@longfn\hsize\textwidth\fi
    \noindent
    \hb@xt@\bibindent{\hss\@makefnmark\enspace}##1}
  \ifaa@twocolumn
  \begingroup
    \twocolumn[\aa@maketitle]
    \@thanks
  \endgroup
  \else
    \begingroup
      \let\thanks\footnote
      \aa@maketitle
    \endgroup
  \fi
\endgroup
  \setcounter{footnote}{0}%
}
\newcommand*\aa@maketitle{%
  \ifaa@contents
    {\aa@contentsnamefont \aa@contentsname\par}
    \addvspace{\aa@aftercontentsnameskip}
  \else\ifaa@contentsAut
    {\aa@contentsAutnamefont \aa@contentsAutname\par}
    \addvspace{\aa@aftercontentsAutnameskip}
  \else
    \vbox to0pt{\vskip\aa@beforeheadboxskip
    \nointerlineskip
  \ifx\aa@idline\empty
    \noindent
    \fbox{%
      \begin{minipage}{\linewidth-2\fboxsep-2\fboxrule}
        \strut
        {\aa@journalnamefont \aa@journalname}
        {\aa@manuscriptnamefont \aa@manuscriptname}\\
        \ifx\aa@msnr\empty
          {\strut\aa@todayheadboxfont \today}%
        \else
          \makebox[0.75\linewidth][r]{%
            \strut
            {\aa@msnrheadboxfont \aa@msnr}%
          }%
        \fi
      \end{minipage}%
    }%
  \else
    \begin{tabular}[b]{@{}l@{}}
      \aa@idlinefont \aa@idline\\
      \aa@doifont \aa@doi\\
      \aa@copyrightfont \aa@copyrightname~\if!\the\yearCop!\the\year\else\the\yearCop\fi
    \end{tabular}
    \hfill\AALogo
  \fi\vss}
  \par
  \ifx\aa@topics\empty%
    \addvspace{\aa@afterheadboxskip}
  \else
    \addvspace{\aa@beforetopicsskip}
    {\noindent\aa@topicsfont \aa@topics}%
    \hfill
    {\aa@topicsnamefont \aa@topicsname\par}
    \addvspace{\aa@aftertopicsskip}
  \fi
  \begingroup
    \centering
    \ifaa@letter
      \ifx\aa@topics\empty
      \addvspace{1.8\aa@beforeletternameskip}
      \else
      \addvspace{\aa@beforeletternameskip}
      \fi
      {\aa@headletternamefont \lettername\par}
      \addvspace{\aa@afterletternameskip}
    \fi
    \ifx\aa@headnote\empty\else
      \ifx\aa@headnote\aa@defaultheadnote\else
        {\aa@headnotefont \aa@headnote\par}
        \addvspace{\aa@afterheadnoteskip}
      \fi
    \fi
    {\aa@titlefont \aa@title\par}
    \ifx\aa@subtitle\empty%
      \ifx\aa@headnote\aa@defaultheadnote
        {\aa@defaultheadnotefont (\aa@headnote)\par}
      \fi
    \addvspace{\aa@aftertitleskip}
    \else
    \addvspace{\aa@aftertitleskip}
      {\aa@subtitlefont \aa@subtitle\par}
      \ifx\aa@headnote\aa@defaultheadnote
        \addvspace{\aa@beforedefaultheadnoteskip}
        {\aa@defaultheadnotefont (\aa@headnote)\par}
      \fi
      \addvspace{\aa@aftersubtitleskip}
    \fi
    \ifnum\aa@nbauthors=\tw@
      \let\and\andtwoname
    \else
      \renewcommand\and{%
        \advance\aa@cnta \@ne
        \ifnum\aa@cnta=\aa@nbauthors
          \andlastname
        \else
          \andname
        \fi
      }
    \fi
    \aa@cnta \@ne
    {\aa@authorfont \aa@author\par}
    \addvspace{\aa@afterauthorskip}
  \endgroup
  \begingroup
    \ifaa@longauth\addvspace{-10pt}
      {\noindent\aa@longauthfont%
       \aa@longauthopen\aa@longauthname\aa@longauthclose\par}
    \else
    \leftskip \aa@headingleftmargin
    \rightskip \aa@headingrightmargin
      \aainstitutename
    \fi
    \vspace{\aa@afterinstituteskip}
    {\noindent\aa@datefont\ifaa@longauth\centering\fi \aa@date\par}
    \ifx\aa@dedication\empty
      \vspace{\aa@afterdateskip}
    \else
    \addvspace\baselineskip
      {\noindent\aa@dedicationfont \aa@dedication\par}
      \vspace{\aa@afterdedicationskip}
    \fi
    {\centering\aa@abstractnamefont \abstractname\par}
    \vspace{\aa@afterabstractnameskip}
    \ifaa@structabstract
      \ifx\aa@abstract@i\longempty\else
        {\noindent\aa@subabstractfont \aa@subabstractiname. }%
        {\aa@abstractfont \ignorespaces\aa@abstract@i\par}
        \vspace{\aa@aftersubabstractskip}
      \fi
      {\noindent\aa@subabstractfont \aa@subabstractiiname. }%
      {\aa@abstractfont \ignorespaces\aa@abstract@ii\par}
      \vspace{\aa@aftersubabstractskip}
      {\noindent\aa@subabstractfont \aa@subabstractiiiname. }%
      {\aa@abstractfont \ignorespaces\aa@abstract@iii\par}
      \vspace{\aa@aftersubabstractskip}
      {\noindent\aa@subabstractfont \aa@subabstractivname. }%
      {\aa@abstractfont \ignorespaces\aa@abstract@iv\par}
      \vspace{\aa@aftersubabstractskip}
      \ifx\aa@abstract@v\longempty\else
        {\noindent\aa@subabstractfont \aa@subabstractvname. }%
        {\aa@abstractfont \ignorespaces\aa@abstract@v\par}
      \fi
    \else
      {\noindent\aa@abstractfont \ignorespaces\aa@abstract@i\par}
    \fi
    \vspace{\aa@afterabstractskip}
    \ifx\aa@keywords\empty\else
      \noindent
      {\aa@keywordnamefont \aa@keywordname. }
      {\aa@keywordfont \ignorespaces\aa@keywords\par}
      \vspace{\aa@afterkeywordsskip}
    \fi
  \endgroup
  \renewcommand*\thanks[1]{}%
  \renewcommand*\inst[1]{}%
  \renewcommand*\fnmsep{}%
  \renewcommand*\and{, }%
  \renewcommand*\\{ }%
  \renewcommand*\vspace[1]{}%
  \renewcommand*\offprints[1]{}%
  \protected@xdef\aa@PDFauthor{\if!\the\authorlink!\aa@author\else\the\authorlink\fi}%
  \protected@xdef\aa@PDFtitle{\noexpand\protect\noexpand\aa@title}%
  \def\object##1{##1}%
  \ifx\aa@authorrunning\empty
    \let\and\relax
    \renewcommand*\\{\unskip\space\ignorespaces}%
    \ifnum\aa@nbauthors=1
      \expandafter\aa@extractoneauthor\aa@author\and\aa@nil
    \else\ifnum\aa@nbauthors=2
      \expandafter\aa@extracttwoauthor\aa@author\and\aa@nil
    \else
      \expandafter\aa@extractoneauthor\aa@author\and\aa@nil
      \protected@xdef\aa@headings{\aa@headings\unskip\space\etalname\space}%
    \fi\fi
  \else
    \xdef\aa@headings{\aa@authorrunning}%
  \fi
  \ifx\aa@titlerunning\empty
    \protected@xdef\aa@headings{\aa@headings\unskip: \aa@title}%
  \else
    \protected@xdef\aa@headings{\aa@headings\unskip: \aa@titlerunning}%
  \fi
  \setbox0=\hbox{%
    \hspace*{\aa@securityheadmargin}%
    \aa@headfont \aa@headings
  }%
  \ifdim\wd0>\aa@fullwidth
    \typeout{^^J^^JIMPORTANT}
    \ClassWarning{aa}
               {The running head built automatically
                 from \string\author\space and \string\title
                 exceeds the pagewidth, please supply a shorter form
                 with \string\authorrunning\space and/or
                 \string\titlerunning\space prior to \string\maketitle}
               {Use \string\titlerunning\space to supply a shorter
                 title and/or \string\authorrunning\space to supply a
                 shorter list of authors.}%
    \typeout{IMPORTANT
    ^^J}
  \fi
  \fi\fi
}
\def\aa@extractoneauthor#1\and#2\aa@nil{%
  \gdef\aa@headings{#1}%
}
\def\aa@extracttwoauthor#1\and#2\and#3\aa@nil{%
  \gdef\aa@headings{#1\andtwoname#2}%
}
\def\sectcounterend{.}
\def\@seccntformat#1{\csname the#1\endcsname\sectcounterend
\hskip\betweenumberspace}
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
          \raggedright
          \hyphenpenalty \@M%
          \interlinepenalty \@M #8\@@par}%
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
      \@svsec #8\hskip\aftertext}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}
\ifaa@contentsAut
  \renewcommand\section{%
    \@startsection{section}
                  {1}
                  {\z@}
                  {-16dd \@plus-4\p@ \@minus-4\p@}
                  {6dd}
                  {\large\bfseries\boldmath\MakeUppercase}%
  }
\else\ifaa@contents
  \renewcommand\section{%
    \@startsection{section}
                  {1}
                  {\z@}
                  {-16dd \@plus-4\p@ \@minus-4\p@}
                  {6dd}
                  {\large\bfseries\boldmath\MakeUppercase}%
  }
\else
  \renewcommand\section{%
    \@startsection{section}
                  {1}
                  {\z@}
                  {-16dd\@plus-4\p@ \@minus-4\p@}
                  {6dd}
                  {\large\sffamily\bfseries\boldmath}%
  }
\fi\fi
\renewcommand\subsection{%
  \@startsection{subsection}
                {2}
                {\z@}
                {-16dd}
                {6dd}
                {\normalsize\sffamily\itshape}%
}
\renewcommand\subsubsection{%
  \@startsection{subsubsection}
                {3}
                {\z@}
                {-16dd}
                {6dd}
                {\normalsize\sffamily}%
}
\renewcommand\paragraph{%
  \@startsection{paragraph}
                {4}
                {\z@}
                {-16dd}
                {-\fontdimen2\font
                  \@plus-\fontdimen3\font
                  \@minus-\fontdimen4\font}
                {\normalsize\sffamily}%
}
\renewcommand\subparagraph{%
  \ClassError{aa}{Command \string\subparagrap\space isn't allowed}
                          {Only commands from \string\section\space to
                            \string\paragraph\space are allowed.}
}
\renewcommand\part{%
  \ClassError{aa}{Command \string\part\space isn't allowed}
                          {Only commands from \string\section\space to
                            \string\paragraph\space are allowed.}
}
\frenchspacing
\clubpenalty=10000
\widowpenalty=10000
\flushbottom
\ifaa@twocolumn
  \sloppy
\fi
\newenvironment{nopar}
  {\let\par\relax%
    \aa@acknowledgementsfont
      {%
        \noindent\aa@acknowledgementsnamefont
        \aa@acknowledgementsname.%
      }%
  \hskip.5em
  }{}
\newenvironment{acknowledgements}
  {%
   \par\addvspace{\aa@beforeacknowledgementsskip}
   \begin{spacing}{0.81}
   \begin{nopar}
  }
  {%
   \end{nopar}
   \end{spacing}
   \addvspace{\aa@afteracknowledgementsskip}
  }
\let\acknowledgement\acknowledgements
\let\endacknowledgement\endacknowledgements
\renewcommand\@makefntext[1]{%
    \noindent
    \hb@xt@\bibindent{\hss\@makefnmark\enspace}#1}
\newcommand\listofobjects{\if@myobjects%
  \section*{\listobjectname}
  \@starttoc{obj}\fi%
}
\newcommand*\l@obj[2]{\noindent`#1' on page #2\par}
\DeclareRobustCommand{\object}[1]{%
   #1\global\let\if@myobjects\iftrue%
   \addcontentsline{obj}{obj}{#1}%
}
\newcommand\longtab[2][]{%
\if!\@aatolm!
  \aalongfloat{table}{#2}%
\else
  \aalongfloat[#1]{section}{#2}%
\fi
}
\newcommand\aalongfloat[3][]{%
  \AtEndDocument{%
    \aa@savecounter \csname c@#2\endcsname \relax
    \csname c@#2\endcsname
  }%
  \edef\aa@nombre{\number\csname c@#2\endcsname}%
\if!\@aatolm!
  \expandafter\AtEndDocument\expandafter{%
    \aa@nombre\relax
    \setlength{\LTcapwidth}{\textwidth}
    \onecolumn
    #3%
    \twocolumn
    \csname c@#2\endcsname \aa@savecounter \relax
  }%
  \stepcounter{#2}%
\else
  \expandafter\AtEndDocument\expandafter{%
    \aa@nombre\relax
    \setlength{\LTcapwidth}{\textwidth}
    \onecolumn
\setcounter{table}{#1}\addtocounter{table}{-1}
\renewcommand\thesection{\Alph{section}}%
\renewcommand\thetable{\thesection.\arabic{table}}
    #3%
    \twocolumn
    \csname c@#2\endcsname \aa@savecounter \relax
  }%
  \setcounter{table}{#1}
\fi
}
\renewcommand*\enddocument{%
   \let\AtEndDocument\@firstofone
   \@enddocumenthook
   \clearpage
   \advance\c@page \m@ne
   \immediate\write\@auxout{\string\newlabel{LastPage}{{\thepage}{\thepage}{}{}{}}}
   \@checkend{document}%
   \clearpage
   \begingroup
     \if@filesw
       \immediate\closeout\@mainaux
       \let\@setckpt\@gobbletwo
       \let\@newl@bel\@testdef
       \@tempswafalse
       \makeatletter \input\jobname.aux
     \fi
     \@dofilelist
     \ifdim \font@submax >\fontsubfuzz\relax
       \@font@warning{Size substitutions with differences\MessageBreak
                  up to \font@submax\space have occurred.\@gobbletwo}%
     \fi
     \@defaultsubs
     \@refundefined
     \if@filesw
       \ifx \@multiplelabels \relax
         \if@tempswa
           \@latex@warning@no@line{Label(s) may have changed.
               Rerun to get cross-references right}%
         \fi
       \else
         \@multiplelabels
       \fi
     \fi
   \endgroup
   \deadcycles\z@\@@end}
\newcommand*\aa@defineentry{%
  \newcommand{\entry}[3]{%
    \par\noindent
    {\aa@entryiifont ##2}%
    \par
    {%
      \noindent
      \raggedright
      {\aa@entryifont ##1}%
      \hfill
      \rlap{\makebox[\aa@entryboxwidth][l]{\hfill##3}}%
      \par
    }%
    \vspace{\aftertext}
  }
}
\ifaa@contents \aa@defineentry \fi
\ifaa@contentsAut \aa@defineentry \fi
\renewcommand\@listi{%
  \leftmargin \leftmargini
  \parsep \z@ \@plus\p@ \@minus\p@
  \topsep 8\p@ \@plus2\p@ \@minus4\p@
  \itemsep \parsep
}
\let\listI\@listi
\@listi
\AtBeginDocument{%
  \let\@listi\listI
  \@listi
}
\renewcommand\@listii{%
  \leftmargin \leftmarginii
  \labelwidth \leftmarginii
  \advance\labelwidth -\labelsep
  \topsep \z@
  \parsep \z@
  \itemsep \z@
}
\renewcommand\@listiii{%
  \leftmargin \leftmarginiii
  \labelwidth \leftmarginiii
  \advance\labelwidth -\labelsep
  \topsep \z@
  \parsep \z@
  \itemsep \z@
  \partopsep \p@ \@minus\p@
 }
\newcommand*\aa@setitemindent[1]{%
  \settowidth{\labelwidth}{\aa@descriptionlabelfont #1}%
  \leftmargini \labelwidth
  \advance\leftmargini \labelsep
  \renewcommand*\@listi{%
    \leftmargin \leftmargini
    \labelwidth \leftmargini
    \advance\labelwidth -\labelsep
    \parsep \parskip
    \topsep \medskipamount
    \itemsep \parskip
    \advance\itemsep -\parsep
  }%
}
\newcommand*\aa@setitemitemindent[1]{%
  \settowidth{\labelwidth}{\aa@descriptionlabelfont #1}%
  \leftmarginii \labelwidth
  \advance\leftmarginii \labelsep
  \renewcommand*\@listii{%
    \leftmargin \leftmarginii
    \labelwidth \leftmarginii
    \advance\labelwidth -\labelsep
    \parsep \parskip
    \topsep \z@
    \itemsep \parskip
    \advance\itemsep -\parsep
  }%
}
\renewcommand*\descriptionlabel[1]{\hspace\labelsep #1\hfil}
\renewcommand*\description{%
  \@ifnextchar[{\aa@description} % ] (pour Emacs)
               {%
                 \list{}{%
                   \labelwidth\z@
                   \itemindent -\leftmargin
                   \let\makelabel\descriptionlabel
                 }%
               }%
}
\newcommand\aa@description[1][]{%
  \relax\ifnum\@listdepth=0
    \aa@setitemindent{#1}%
  \else\ifnum\@listdepth=1
    \aa@setitemitemindent{#1}%
  \fi\fi
  \list{--}{\let\makelabel\aa@descriptionlabel}%
}
\newcommand*\aa@descriptionlabel[1]{%
  \aa@descriptionlabelfont
  #1%
  \hfil
}
\newcommand*\aa@theorem@aftercount{.}
\newcommand*\aa@theorem@aftertitle{.}
\renewcommand*\newtheorem{%
  \@ifstar{\aa@theorem@withstar}{\aa@theorem@nostar}%
}
\newcommand*\aa@theorem@nostar{%
  \aa@numberedtheoremtrue
  \aa@newtheorem
}
\newcommand*\aa@theorem@withstar{%
  \aa@numberedtheoremfalse
  \aa@newtheorem
}
\newcommand*\aa@newtheorem[1]{%
  \def\aa@theorem@name{#1}%
  \expandafter\newif\csname ifaa@theorem@#1@numbered\endcsname
  \ifaa@numberedtheorem
    \csname aa@theorem@#1@numberedtrue\endcsname
  \fi
  \aa@@newtheorem
}
\newif\ifaa@newcounter
\newcommand*\aa@@newtheorem[2][]{%
  \def\aa@arg{#1}%
  \ifx\aa@arg\empty
    \expandafter\edef
      \csname aa@theorem@\aa@theorem@name @counter\endcsname
      {AA\aa@theorem@name}%
    \aa@newcountertrue
  \else
    \expandafter\ifx
        \csname aa@theorem@#1@title\endcsname
        \relax
      \ClassError{aa}
                 {#1 isn't a theorem-like environment}
                 {You have to create #1 with \string\newtheorem\space
                   before use its counter.}%
      \expandafter\edef
        \csname aa@theorem@\aa@theorem@name @counter\endcsname
        {AA\aa@theorem@name}%
      \aa@newcountertrue
    \else
      \expandafter\def
        \csname aa@theorem@\aa@theorem@name @counter\endcsname
        {AA#1}%
      \aa@newcounterfalse
    \fi
  \fi
  \expandafter\newcommand
    \csname aa@theorem@\aa@theorem@name @title\endcsname
    {#2}%
  \aa@@@newtheorem
}
\newcommand*\aa@@@newtheorem[3][]{%
  \ifaa@newcounter
    \def\aa@arg{#1}%
    \ifx\aa@arg\empty
      \newcounter{\csname aa@theorem@\aa@theorem@name @counter\endcsname}%
    \else
      \newcounter{\csname aa@theorem@\aa@theorem@name @counter\endcsname}[#1]%
      \expandafter\xdef
        \csname theaa@theorem@\aa@theorem@name @counter\endcsname{%
          \expandafter\noexpand\csname the#1\endcsname.%
          \noexpand\arabic{aa@theorem@\aa@theorem@name @counter}%
      }%
    \fi
  \fi
  \expandafter\newcommand
    \csname aa@theorem@\aa@theorem@name @titlefont\endcsname{#2}%
  \expandafter\newcommand
    \csname aa@theorem@\aa@theorem@name @bodyfont\endcsname{#3}%
  \newenvironment{\aa@theorem@name}[1][]{%
    \def\aa@theorem@optionnal{##1}%
    \par\addvspace{\VSpaceBeforeTheorem}%
    \noindent\hspace*{\HSpaceIndentTheorem}%
    \begingroup
      \csname aa@theorem@\@currenvir @titlefont\endcsname
      \csname aa@theorem@\@currenvir @title\endcsname
      \csname ifaa@theorem@\@currenvir @numbered\endcsname
        \refstepcounter{\csname aa@theorem@\@currenvir @counter\endcsname}%
        \space
        \csname the\csname aa@theorem@\@currenvir @counter\endcsname\endcsname
      \fi
      \ifx\aa@theorem@optionnal\empty
         \aa@theorem@aftercount
      \else
        \space(\aa@theorem@optionnal)%
        \aa@theorem@aftertitle
      \fi
    \endgroup
    \ifdim\HVSpaceSepTheorem<\z@
      \par\nobreak
      \vspace{-\HVSpaceSepTheorem}%
      \noindent\hspace*{\HSpaceIndentTheorem}%
    \else
      \hspace*{\HVSpaceSepTheorem}%
    \fi
    \csname aa@theorem@\@currenvir @bodyfont\endcsname
  }%
  {%
    \par\addvspace{\VSpaceAfterTheorem}
  }%
}
\newtheorem*{proof}{Proof}{\itshape}{\rmfamily}
\newtheorem{theorem}{Theorem}[section]{\bfseries}{\itshape}
\newtheorem{corollary}[theorem]{Corollary}{\bfseries}{\itshape}
\newtheorem{definition}[theorem]{Definition}{\bfseries}{\rmfamily}
\newtheorem{example}[theorem]{Example}{\itshape}{\rmfamily}
\newtheorem{exercise}[theorem]{Exercise}{\itshape}{\rmfamily}
\newtheorem{lemma}[theorem]{Lemma}{\bfseries}{\itshape}
\newtheorem{note}[theorem]{Note}{\itshape}{\rmfamily}
\newtheorem{problem}[theorem]{Problem}{\itshape}{\rmfamily}
\newtheorem{proposition}[theorem]{Proposition}{\bfseries}{\itshape}
\newtheorem{question}[theorem]{Question}{\itshape}{\rmfamily}
\newtheorem{remark}[theorem]{Remark}{\itshape}{\rmfamily}
\newtheorem{solution}[theorem]{Solution}{\itshape}{\rmfamily}
\setlength{\textfloatsep}{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength{\dbltextfloatsep}{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength{\@fptop}{\z@}
\setlength{\@fpsep}{12\p@}
\setlength{\@dblfptop}{\z@}
\setlength{\@dblfpsep}{12\p@}
\setcounter{topnumber}{4}
\renewcommand*\topfraction{.9}
\setcounter{bottomnumber}{2}
\renewcommand*\bottomfraction{.5}
\setcounter{totalnumber}{6}
\renewcommand*\textfraction{.05}
\renewcommand*\floatpagefraction{.95}
\setcounter{dbltopnumber}{3}
\renewcommand*\dbltopfraction{.9}
\renewcommand*\dblfloatpagefraction{.95}
\renewcommand*\fps@figure{tbp}
\renewcommand*\fps@table{tbp}
\def\sidecaption#1\caption{%
  \setbox0=\hbox{\ignorespaces#1\unskip}%
  \aa@sidecaptionwidth=\linewidth
  \advance\aa@sidecaptionwidth -\wd0
  \advance\aa@sidecaptionwidth -\aa@figgap
  \ifdim\aa@sidecaptionwidth<\aa@sidecaptionminwidth
    \ClassWarning{aa}{%
      No sufficient room for the legend;
      using normal \string\caption.%
    }%
    \let\aa@caption\@caption
    \unhbox0
  \else
    \let\aa@caption\aa@sidecaption
    \dimen0 \dp0
    \advance\dimen0 \ht0
    \dp0 \dp\strutbox
    \advance\dimen0-\dp0
    \ht0 \dimen0
    \unhbox0\hspace{\aa@figgap}%
  \fi
  \refstepcounter\@captype
  \@dblarg{\aa@caption\@captype}%
}
\long\def\aa@sidecaption#1[#2]#3{%
  \addcontentsline
    {\csname ext@#1\endcsname}
    {#1}
    {\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
  \begingroup
    \@parboxrestore
    \aa@makesidecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}%
    \par
  \endgroup
}
\newcommand*\aa@capstrut{%
  \leavevmode
  \vrule \@width\z@ \@height\ht\strutbox \@depth\z@
}
\def\aa@@makecaption#1{%
    {\aa@headcaptionfont #1\aa@floatcounterend}%
}
\long\def\aa@makesidecaption#1#2{%
  \parbox[b]{\aa@sidecaptionwidth}{%
    \aa@captionfont
    \aa@@makecaption{#1}
    \aa@capstrut
    #2%
  }%
}
\renewcommand\@makecaption[2]{%
  \aa@captionfont
  \ifx\@captype\aa@figtype
    \par\vskip\VSpaceBeforeFigCaption
    \aa@@makecaption{#1}%
    \aa@capstrut%
    #2%
    \par
  \else
    {\aa@headcaptionfont #1\aa@floatcounterend}%
    \aa@capstrut
    #2%
    \par
    \vskip\VSpaceAfterTabCaption
  \fi
}
\newcommand\tablebib[1]{\VSpaceBeforeTabFoot=1ex%
  \par\vspace{\VSpaceBeforeTabBib}
  \noindent
  \begin{minipage}{\linewidth}
    {\aa@tablebibnamefont\aa@tablebibname.}~%
    \aa@tablebibfont
    \ignorespaces
    #1%
  \end{minipage}%
}
\newcommand\tablefoot[1]{\VSpaceBeforeTabBib=1ex%
  \par\vspace{\VSpaceBeforeTabFoot}
  \noindent
  \begin{minipage}{\linewidth}
    {\aa@tablefootnamefont\aa@tablefootname.}~%
    \aa@tablefootfont
    \ignorespaces
    #1%
  \end{minipage}%
}
\newcommand*\tablefootmark[1]{%
  \unskip
  \hbox{\@textsuperscript{\normalfont\itshape\ignorespaces#1}}%
  \,%
  \ignorespaces
}
\newcommand\tablefoottext[2]{%
  \hbox{\@textsuperscript{\normalfont({\itshape\ignorespaces#1})}}%
  ~%
  \ignorespaces
  #2\ \ignorespaces%
}
\setlength\bibindent{1em}
\renewenvironment{thebibliography}[1]{%
  \begingroup
  \section*{\refname}%
  \@mkboth{\refname}{\refname}%
  \aa@bibliographyfont
  \list{\@biblabel{\@arabic\c@enumiv}}%
  {%
    \settowidth\labelwidth{\@biblabel{#1}}%
    \ifaa@bibnumber
    \leftmargin\labelwidth
    \advance\leftmargin\labelsep
    \else
    \leftmargin\bibindent
    \setlength{\itemindent}{-\leftmargin}
    \fi
    \@openbib@code
    \usecounter{enumiv}%
    \let\p@enumiv\@empty
    \renewcommand\theenumiv{\@arabic\c@enumiv}%
  }%
  \sloppy
  \clubpenalty4000
  \@clubpenalty \clubpenalty
  \widowpenalty4000%
  \sfcode`\.\@m
}{%
  \endlist\par\endgroup%
  \ifaa@longauth\kern6pt\hrule\kern6pt\aainstitutename\aa@longauthfalse\fi%
  \ifaa@referee\if!\csname appendix\endcsname!\relax\else\listofobjects\fi\fi
  \def\@noitemerr{\ClassWarning{aa}{Empty `thebibliography' environment}}%
}
\def \@floatboxreset {%
     \reset@font
     \ifx\@captype\aa@tabtype
       \small
       \spacetab
     \fi
     \@setnobreak
     \@setminipage
}
\def\spacetab{%
\abovedisplayshortskip-10pt%
\belowdisplayshortskip0pt%
}
\DeclareRobustCommand{\ion}[2]{\textup{#1\,\textsc{\lowercase{#2}}}}
\newcommand*\element[1][]{%
  \def\aa@element@tr{#1}%
  \aa@element
}
\newcommand*\aa@element[1][]{%
  \def\aa@element@tl{#1}%
  \aa@@element
}
\newcommand*\aa@@element[1][]{%
  \def\aa@element@bl{#1}%
  \aa@@@element
}
\newcommand*\aa@@@element[2][]{%
  \def\aa@element@br{#1}%
  \ensuremath{%
    \vphantom{x}^{\aa@element@tl}_{\aa@element@bl}
    \mathrm{#2}
    \vphantom{x}^{\aa@element@tr}_{\aa@element@br}
  }%
}
\renewcommand*\vec[1]{\ensuremath{\boldsymbol{#1}}}
\newcommand*\tens[1]{\ensuremath{\mathsf{#1}}}
\newcommand*\sq{%
  \ifmmode
    \Box
  \else
    \begingroup
      \unskip\nobreak
      \hfil
      \penalty50
      \hskip1em
      \null\nobreak\hfil$\Box$%
      \parfillskip=0pt
      \finalhyphendemerits=0
      \endgraf
    \endgroup
  \fi
}
\newcommand*\sun{\ensuremath{\odot}}
\newcommand*\degr{\ensuremath{^\circ}}
\newcommand*\arcmin{\ensuremath{^\prime}}
\newcommand*\arcsec{\ensuremath{^{\prime\prime}}}
\def\utw{\ensuremath{\smash{\rlap{\lower5pt\hbox{$\sim$}}}}}
\def\udtw{\ensuremath{\smash{\rlap{\lower6pt\hbox{$\approx$}}}}}
\newcommand*\fd{\ensuremath{\overset{\text{d}}{.}}}
\newcommand*\fh{\ensuremath{\overset{\text{h}}{.}}}
\newcommand*\fm{\ensuremath{\overset{\text{m}}{.}}}
\newcommand*\fs{\ensuremath{\overset{\text{s}}{.}}}
\newcommand*\fdg{\ensuremath{\overset{\circ}{.}}}
\newcommand*\farcm{\ensuremath{\overset{\prime}{.}}}
\newcommand*\farcs{\ensuremath{\overset{\prime\prime}{.}}}
\newcommand*\fp{\ensuremath{\overset{\scriptscriptstyle\text{p}}{.}}}
\newcommand*\cor{\ensuremath{\mathrel{\widehat{=}}}}
\newcommand*\aa@centerstack[2]{%
  \ensuremath{%
    \mathrel{%
      \mathchoice{%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\displaystyle##$\hfil\cr#1\cr#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\textstyle##$\hfil\cr#1\cr#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\scriptstyle##$\hfil\cr#1\cr#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\scriptscriptstyle##$\hfil\cr#1\cr#2\cr}%
        }%
      }%
    }%
  }%
}
\newcommand*\la{\aa@centerstack{<}{\sim}}
\newcommand*\ga{\aa@centerstack{>}{\sim}}
\newcommand*\sol{\aa@centerstack{\sim}{<}}
\newcommand*\sog{\aa@centerstack{\sim}{>}}
\newcommand*\lse{\aa@centerstack{<}{\simeq}}
\newcommand*\gse{\aa@centerstack{>}{\simeq}}
\newcommand*\aa@centerstackskip[6]{%
  \ensuremath{%
    \mathrel{%
      \mathchoice{%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\displaystyle##$\hfil\cr
            #1\cr\noalign{\vskip#3}#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\textstyle##$\hfil\cr
            #1\cr\noalign{\vskip#4}#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\scriptstyle##$\hfil\cr
            #1\cr\noalign{\vskip#5}#2\cr}%
        }%
      }
      {%
        \vcenter{%
          \offinterlineskip
          \halign{\hfil$\scriptscriptstyle##$\hfil\cr
            #1\cr\noalign{\vskip#6}#2\cr}%
        }%
      }%
    }%
  }%
}
\newcommand*\grole{\aa@centerstackskip{>}{<}{-1.5\p@}{-1.5\p@}{-1\p@}{-.5\p@}}
\newcommand*\leogr{\aa@centerstackskip{<}{>}{-1.5\p@}{-1.5\p@}{-1\p@}{-.5\p@}}
\newcommand*\loa{\aa@centerstack{<}{\approx}}
\newcommand*\goa{\aa@centerstack{>}{\approx}}
\newcommand*\diameter{%
  \ensuremath{%
    \mathchoice{%
      \ooalign{%
        \hfil\hbox{$\displaystyle/$}\hfil\crcr
        {\hbox{$\displaystyle\mathchar"20D$}}%
      }%
    }
    {%
      \ooalign{%
        \hfil\hbox{$\textstyle/$}\hfil\crcr
        {\hbox{$\textstyle\mathchar"20D$}}%
      }%
    }
    {%
      \ooalign{%
        \hfil\hbox{$\scriptstyle/$}\hfil\crcr
        {\hbox{$\scriptstyle\mathchar"20D$}}%
      }%
    }
    {%
      \ooalign{%
        \hfil\hbox{$\scriptscriptstyle/$}\hfil\crcr
        {\hbox{$\scriptscriptstyle\mathchar"20D$}}%
      }%
    }%
  }%
}
\let\olddiameter\diameter
\newcommand*\getsto{\aa@centerstack{\gets}{\to}}
\newcommand*\lid{\aa@centerstackskip{<}{=}{1.2\p@}{1.2\p@}{1\p@}{0.9\p@}}
\newcommand*\gid{\aa@centerstackskip{>}{=}{1.2\p@}{1.2\p@}{1\p@}{0.9\p@}}
\newcommand*\bbbone{%
  \ensuremath{%
    \mathchoice{\mathrm 1\mskip-4mu l}
               {\mathrm 1\mskip-4mu l}
               {\mathrm 1\mskip-4.5mu l}
               {\mathrm 1\mskip-5mu l}%
  }%
}
\newcommand*\aap{A\&A}
\let\astap=\aap
\newcommand*\aapr{A\&A~Rev.}
\newcommand*\aaps{A\&AS}
\newcommand*\actaa{Acta Astron.}
\newcommand*\aj{AJ}
\newcommand*\ao{Appl.~Opt.}
\let\applopt\ao
\newcommand*\apj{ApJ}
\newcommand*\apjl{ApJ}
\let\apjlett\apjl
\newcommand*\apjs{ApJS}
\let\apjsupp\apjs
\newcommand*\aplett{Astrophys.~Lett.}
\newcommand*\apspr{Astrophys.~Space~Phys.~Res.}
\newcommand*\apss{Ap\&SS}
\newcommand*\araa{ARA\&A}
\newcommand*\azh{AZh}
\newcommand*\baas{BAAS}
\newcommand*\bac{Bull. astr. Inst. Czechosl.}
\newcommand*\bain{Bull.~Astron.~Inst.~Netherlands}
\newcommand*\caa{Chinese Astron. Astrophys.}
\newcommand*\cjaa{Chinese J. Astron. Astrophys.}
\newcommand*\fcp{Fund.~Cosmic~Phys.}
\newcommand*\gca{Geochim.~Cosmochim.~Acta}
\newcommand*\grl{Geophys.~Res.~Lett.}
\newcommand*\iaucirc{IAU~Circ.}
\newcommand*\icarus{Icarus}
\newcommand*\jcap{J. Cosmology Astropart. Phys.}
\newcommand*\jcp{J.~Chem.~Phys.}
\newcommand*\jgr{J.~Geophys.~Res.}
\newcommand*\jqsrt{J.~Quant.~Spectr.~Rad.~Transf.}
\newcommand*\jrasc{JRASC}
\newcommand*\memras{MmRAS}
\newcommand*\memsai{Mem.~Soc.~Astron.~Italiana}
\newcommand*\mnras{MNRAS}
\newcommand*\na{New A}
\newcommand*\nar{New A Rev.}
\newcommand*\nat{Nature}
\newcommand*\nphysa{Nucl.~Phys.~A}
\newcommand*\pasa{PASA}
\newcommand*\pasj{PASJ}
\newcommand*\pasp{PASP}
\newcommand*\physrep{Phys.~Rep.}
\newcommand*\physscr{Phys.~Scr}
\newcommand*\planss{Planet.~Space~Sci.}
\newcommand*\pra{Phys.~Rev.~A}
\newcommand*\prb{Phys.~Rev.~B}
\newcommand*\prc{Phys.~Rev.~C}
\newcommand*\prd{Phys.~Rev.~D}
\newcommand*\pre{Phys.~Rev.~E}
\newcommand*\prl{Phys.~Rev.~Lett.}
\newcommand*\procspie{Proc.~SPIE}
\newcommand*\qjras{QJRAS}
\newcommand*\rmxaa{Rev. Mexicana Astron. Astrofis.}
\newcommand*\skytel{S\&T}
\newcommand*\solphys{Sol.~Phys.}
\newcommand*\sovast{Soviet~Ast.}
\newcommand*\ssr{Space~Sci.~Rev.}
\newcommand*\zap{ZAp}
\DeclareMathAlphabet{\mathsfsl}{OT1}{txss}{m}{sl}
\g@addto@macro{\endlongtable}{\addvspace{-3.5mm}}
\AtBeginDocument{%
\@ifpackageloaded{lscape}{%
\g@addto@macro{\landscape}{\centering\setlength{\LTcapwidth}{\aa@fullheight}}}{}
\@ifpackageloaded{pdflscape}{%
\g@addto@macro{\landscape}{\centering\setlength{\LTcapwidth}{\aa@fullheight}}}{}
\@ifpackageloaded{aalongtable}%
{\def\LS@rot{%
  \setbox\@outputbox\vbox
  to\aa@fullheight{\vfill\centerline{\rotatebox{90}{\box\@outputbox}}\vfill}}}{}
\@ifpackageloaded{subfigure}%
{\renewcommand*{\thesubfigure}{\alph{subfigure}}
\renewcommand*{\@thesubfigure}{(\thesubfigure)\hskip\subfiglabelskip}
\renewcommand*{\@@thesubfigure}{(\thesubfigure)}%
\renewcommand*{\thesubtable}{\alph{subtable}}
\renewcommand*{\@thesubtable}{(\thesubtable)\hskip\subfiglabelskip}
\renewcommand*{\@@thesubtable}{(\thesubtable)}}{}
\@ifpackageloaded{wasysym}%
{\let\diameter\olddiameter}{}
\@ifpackageloaded{txfonts}%
{%
  \newcommand*\bbba{\ensuremath{\varmathbb{A}}}
  \newcommand*\bbbb{\ensuremath{\varmathbb{B}}}
  \newcommand*\bbbc{\ensuremath{\varmathbb{C}}}
  \newcommand*\bbbd{\ensuremath{\varmathbb{D}}}
  \newcommand*\bbbe{\ensuremath{\varmathbb{E}}}
  \newcommand*\bbbf{\ensuremath{\varmathbb{F}}}
  \newcommand*\bbbg{\ensuremath{\varmathbb{G}}}
  \newcommand*\bbbh{\ensuremath{\varmathbb{H}}}
  \newcommand*\bbbi{\ensuremath{\varmathbb{I}}}
  \newcommand*\bbbj{\ensuremath{\varmathbb{J}}}
  \newcommand*\bbbk{\ensuremath{\varmathbb{K}}}
  \newcommand*\bbbl{\ensuremath{\varmathbb{L}}}
  \newcommand*\bbbm{\ensuremath{\varmathbb{M}}}
  \newcommand*\bbbn{\ensuremath{\varmathbb{N}}}
  \newcommand*\bbbo{\ensuremath{\varmathbb{O}}}
  \newcommand*\bbbp{\ensuremath{\varmathbb{P}}}
  \newcommand*\bbbq{\ensuremath{\varmathbb{Q}}}
  \newcommand*\bbbr{\ensuremath{\varmathbb{R}}}
  \newcommand*\bbbs{\ensuremath{\varmathbb{S}}}
  \newcommand*\bbbt{\ensuremath{\varmathbb{T}}}
  \newcommand*\bbbu{\ensuremath{\varmathbb{U}}}
  \newcommand*\bbbv{\ensuremath{\varmathbb{V}}}
  \newcommand*\bbbw{\ensuremath{\varmathbb{W}}}
  \newcommand*\bbbx{\ensuremath{\varmathbb{X}}}
  \newcommand*\bbby{\ensuremath{\varmathbb{Y}}}
  \newcommand*\bbbz{\ensuremath{\varmathbb{Z}}}}%
{%
  \newcommand*\bbba{\ensuremath{\mathbf{A}}}
  \newcommand*\bbbb{\ensuremath{\mathbf{B}}}
  \newcommand*\bbbc{\ensuremath{\mathbf{C}}}
  \newcommand*\bbbd{\ensuremath{\mathbf{D}}}
  \newcommand*\bbbe{\ensuremath{\mathbf{E}}}
  \newcommand*\bbbf{\ensuremath{\mathbf{F}}}
  \newcommand*\bbbg{\ensuremath{\mathbf{G}}}
  \newcommand*\bbbh{\ensuremath{\mathbf{H}}}
  \newcommand*\bbbi{\ensuremath{\mathbf{I}}}
  \newcommand*\bbbj{\ensuremath{\mathbf{J}}}
  \newcommand*\bbbk{\ensuremath{\mathbf{K}}}
  \newcommand*\bbbl{\ensuremath{\mathbf{L}}}
  \newcommand*\bbbm{\ensuremath{\mathbf{M}}}
  \newcommand*\bbbn{\ensuremath{\mathbf{N}}}
  \newcommand*\bbbo{\ensuremath{\mathbf{O}}}
  \newcommand*\bbbp{\ensuremath{\mathbf{P}}}
  \newcommand*\bbbq{\ensuremath{\mathbf{Q}}}
  \newcommand*\bbbr{\ensuremath{\mathbf{R}}}
  \newcommand*\bbbs{\ensuremath{\mathbf{S}}}
  \newcommand*\bbbt{\ensuremath{\mathbf{T}}}
  \newcommand*\bbbu{\ensuremath{\mathbf{U}}}
  \newcommand*\bbbv{\ensuremath{\mathbf{V}}}
  \newcommand*\bbbw{\ensuremath{\mathbf{W}}}
  \newcommand*\bbbx{\ensuremath{\mathbf{X}}}
  \newcommand*\bbby{\ensuremath{\mathbf{Y}}}
  \newcommand*\bbbz{\ensuremath{\mathbf{Z}}}}
\@ifpackageloaded{hyperref}%
{}{\RequirePackage[colorlinks=true,
    linkcolor=blue, citecolor=blue, filecolor=blue, urlcolor=blue]{hyperref}}%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Macro ORCID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\orcid}[1]{\unskip\protect\href{https://orcid.org/#1}{\protect\includegraphics[width=8pt,clip]{logo_orcid}}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}
\setlength\leftmargini  {\parindent}
\leftmargin  \leftmargini
\setlength\leftmarginii {\parindent}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {.5em}
\setlength\leftmarginvi {.5em}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\renewcommand\labelenumi{\theenumi.}
\renewcommand\labelenumii{(\theenumii)}
\renewcommand\labelenumiii{\theenumiii.}
\renewcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\renewcommand\labelitemi{\normalfont\bfseries --}
\renewcommand\labelitemii{\normalfont\bfseries --}
\renewcommand\labelitemiii{$\m@th\bullet$}
\renewcommand\labelitemiv{$\m@th\cdot$}
\def\arxivprefixe{arXiv}
\def\arxivprefixesep{:}
\def\arxiv#1{%
{\tt \arxivprefixe\arxivprefixesep\ignorespaces#1}%
}
\newcommand{\eprint}[2][]{%
{\tt\if!#1!#2\else#1\arxivprefixesep\ignorespaces#2\fi}%
}
% \corr : Mute tag used by LE
\newcommand{\corr}[1]{#1}
\endinput
%%
%% End of file `aa.cls'.
