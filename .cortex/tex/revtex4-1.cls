%%
%% This is file `revtex4-1.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% revtex4-1.dtx  (with options: `kernel')
%% ltxutil.dtx  (with options: `kernel')
%% ltxfront.dtx  (with options: `kernel')
%% ltxgrid.dtx  (with options: `kernel')
%% revtex4-1.dtx  (with options: `options')
%% ltxutil.dtx  (with options: `options')
%% ltxfront.dtx  (with options: `options')
%% ltxgrid.dtx  (with options: `options')
%% revtex4-1.dtx  (with options: `package')
%%
%% This file is part of the APS files in the REVTeX 4 distribution.
%% For the version number, search on the string
%% Original version by David Carlisle
%% Modified by Arthur Ogawa (mailto:arthur_ogawa at sbcglobal dot net)
%%
%% Copyright (c) 2009 The American Physical Society.
%% http://publish.aps.org/revtex4/
%% mailto:revtex4@aps.org
%%
%% See the REVTeX 4 README file for restrictions and more information.
%%
\NeedsTeXFormat{LaTeX2e}[1996/12/01]%
\ProvidesClass{revtex4-1}
 [2010/07/25/20:33:00 4.1r (http://publish.aps.org/revtex4/ for documentation)]% \fileversion
\let\class@name\@gtempa
\GenericInfo{}{\space
 Copyright (c) 2009 The American Physical Society.^^J
 mailto:revtex@aps.org^^J
 Licensed under the LPPL:^^Jhttp://www.ctan.org/tex-archive/macros/latex/base/lppl.txt^^J
 Arthur Ogawa <arthur_ogawa at sbcglobal dot net>^^J
 Based on work by David Carlisle <david at dcarlisle.demon.co.uk>%
 \@gobble
}%
\if@compatibility
  \edef\reserved@a{\errhelp{%
   Change your \string\documentstyle\space statement to
   \string\documentclass\space and rerun.
  }}\reserved@a
  \errmessage{You cannot run \class@name\space in compatability mode}%
  \expandafter\@@end
\fi
\typeout{%
ltxutil%
 [2010/07/25/20:33:00 4.1r utilities package (portions licensed from W. E. Baxter web at superscript.com)]% \fileversion
}%
\def\class@err#1{\ClassError{\class@name}{#1}\@eha}%
\def\class@warn#1{\ClassWarningNoLine{\class@name}{#1}}%
\def\class@info#1{\ClassInfo{\class@name}{#1}}%
\def\obsolete@command#1{%
 \class@warn@end{Command \string#1\space is obsolete.^^JPlease remove from your document}%
 \global\let#1\@empty
 #1%
}%
\def\replace@command#1#2{%
 \class@warn@end{Command \string#1\space is obsolete;^^JUse \string#2\space instead}%
 \global\let#1#2%
 #1%
}%
\def\replace@environment#1#2{%
 \class@warn@end{Environment #1 is obsolete;^^JUse #2 instead}%
 \glet@environment{#1}{#2}%
 \@nameuse{#1}%
}%
\def\incompatible@package#1{%
 \@ifpackageloaded{#1}{%
  \def\@tempa{I cannot continue. You must remove the \string\usepackage\ statement that caused that package to be loaded.}%
  \ClassError{\class@name}{The #1 package cannot be used with \class@name}%
  \@tempa\stop
 }{%
  \class@info{#1 was not loaded (OK!)}%
 }%
}%
\def\class@warn@end#1{%
 \gappdef\class@enddocumenthook{\class@warn{#1}}%
}%
\ifx\undefined\class@name
 \def\class@name{ltxutil}%
 \class@warn{You should define the class name before reading in this package. Using default}%
\fi
\def\t@{to}%
\dimendef\dimen@iii\thr@@
\def\halignt@{\halign\t@}%
\chardef\f@ur=4\relax
\chardef\cat@letter=11\relax
\chardef\other=12\relax
\def\let@environment#1#2{%
 \expandafter\let
 \csname#1\expandafter\endcsname\csname#2\endcsname
 \expandafter\let
 \csname end#1\expandafter\endcsname\csname end#2\endcsname
}%
\def\glet@environment#1#2{%
 \global\expandafter\let
 \csname#1\expandafter\endcsname\csname#2\endcsname
 \global\expandafter\let
 \csname end#1\expandafter\endcsname\csname end#2\endcsname
}%
\newcommand\tracingplain{%
 \tracingonline\z@\tracingcommands\z@\tracingstats\z@
 \tracingpages\z@\tracingoutput\z@\tracinglostchars\@ne
 \tracingmacros\z@\tracingparagraphs\z@\tracingrestores\z@
 \showboxbreadth5\showboxdepth3\relax %\errorstopmode
 }%
\newcommand\traceoutput{%
 \appdef\@resetactivechars{\showoutput}%
}%
\newcommand\say[1]{\typeout{<\noexpand#1=\meaning#1>}}%
\newcommand\saythe[1]{\typeout{<\noexpand#1=\the#1>}}%
\def\fullinterlineskip{\prevdepth\z@}%
\countdef\count@i\@ne
\countdef\count@ii\tw@
\long\def\prepdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@ii\the\toks@}%
}%
\long\def\appdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@\the\toks@ii}%
}%
\long\def\gappdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \global\edef#1{\the\toks@\the\toks@ii}%
}%
\long\def\appdef@val#1#2{%
 \appdef#1{{#2}}%
}%
\long\def\appdef@e#1#2{%
 \expandafter\appdef
 \expandafter#1%
 \expandafter{#2}%
}%
\long\def\appdef@eval#1#2{%
 \expandafter\appdef@val
 \expandafter#1%
 \expandafter{#2}%
}%
\toksdef\toks@ii=\tw@
\long\def\@ifxundefined#1{\@ifx{\undefined#1}}%
\long\def\@ifnotrelax#1#2#3{\@ifx{\relax#1}{#3}{#2}}%
\long\def\@argswap#1#2{#2#1}%
\long\def\@argswap@val#1#2{#2{#1}}%
\def\@ifxundefined@cs#1{\expandafter\@ifx\expandafter{\csname#1\endcsname\relax}}%
\def\@boolean#1#2{%
  \long\def#1{%
    #2% \if<something>
      \expandafter\true@sw
    \else
      \expandafter\false@sw
    \fi
  }%
}%
\def\@boole@def#1#{\@boolean{#1}}% Implicit #2
\def\@booleantrue#1{\let#1\true@sw}%
\def\@booleanfalse#1{\let#1\false@sw}%
\@boole@def\@ifx#1{\ifx#1}%
\@boole@def\@ifx@empty#1{\ifx\@empty#1}%
\@boole@def\@if@empty#1{\if!#1!}%
\def\@if@sw#1#2{#1\expandafter\true@sw\else\expandafter\false@sw#2}%
\@boole@def\@ifdim#1{\ifdim#1}%
\@boole@def\@ifeof#1{\ifeof#1}%
\@boole@def\@ifhbox#1{\ifhbox#1}%
\@boole@def\@ifhmode{\ifhmode}%
\@boole@def\@ifinner{\ifinner}%
\@boole@def\@ifmmode{\ifmmode}%
\@boole@def\@ifnum#1{\ifnum#1}%
\@boole@def\@ifodd#1{\ifodd#1}%
\@boole@def\@ifvbox#1{\ifvbox#1}%
\@boole@def\@ifvmode{\ifvmode}%
\@boole@def\@ifvoid#1{\ifvoid#1}%
\long\def\true@sw#1#2{#1}%
\long\def\false@sw#1#2{#2}%
\long\def\loopuntil#1{#1{}{\loopuntil{#1}}}%
\long\def\loopwhile#1{#1{\loopwhile{#1}}{}}%
\def\@provide#1{%
 \@ifx{\undefined#1}{\true@sw}{\@ifx{\relax#1}{\true@sw}{\false@sw}}%
 {\def#1}{\def\j@nk}%
}%
\prepdef\document{%
 \endgroup
 \document@inithook
 \true@sw{}%
}%
\let\document@inithook\@empty
\appdef\document@inithook{%
 \AtBeginDocument{\class@documenthook}%
}%
\AtEndDocument{%
 \class@enddocumenthook
}%
\let\class@documenthook\@empty
\let\class@enddocumenthook\@empty
\def\enddocument{%
 \let\AtEndDocument\@firstofone
 \@enddocumenthook
 \@checkend{document}%
 \clear@document
 \check@aux
 \deadcycles\z@
 \@@end
}%
\def\check@aux{\do@check@aux}%
\def\do@check@aux{%
 \@if@sw\if@filesw\fi{%
  \immediate\closeout\@mainaux
  \let\@setckpt\@gobbletwo
  \let\@newl@bel\@testdef
  \@tempswafalse
  \makeatletter
  \input\jobname.aux\relax
 }{}%
 \@dofilelist
 \@ifdim{\font@submax >\fontsubfuzz\relax}{%
  \@font@warning{%
   Size substitutions with differences\MessageBreak
   up to \font@submax\space have occured.\@gobbletwo
  }%
 }{}%
 \@defaultsubs
 \@refundefined
 \@if@sw\if@filesw\fi{%
  \@ifx{\@multiplelabels\relax}{%
   \@if@sw\if@tempswa\fi{%
    \@latex@warning@no@line{%
     Label(s) may have changed.
     Rerun to get cross-references right
    }%
   }{}%
  }{%
    \@multiplelabels
  }%
 }{}%
}%
\def\clear@document{%
 \clearpage
 \do@output@cclv{%
  \Call@AfterLastShipout
 }%
}%
\appdef\class@documenthook{%
 \providecommand\Call@AfterLastShipout{}%
}%
\def\class@extension#1#2{%
 \IfFileExists{#1.#2}{%
  \expandafter\class@extensionfile\csname ver@\@currname.\@currext\endcsname{#1}#2%
 }{%
  \csname rtx@#1\endcsname
 }%
}%
\def\class@extensionfile#1#2#3{%
 \@pass@ptions#3\@unusedoptionlist{#2}%
 \global\let\@unusedoptionlist\@empty
 \expandafter\class@ext@hook\csname#2.#3-h@@k\endcsname#1{#2}#3%
}%
\def\class@ext@hook#1#2#3#4{%
 \@pushfilename@ltx
 \makeatletter
 \let\CurrentOption\@empty
 \@reset@ptions
 \let#1\@empty
 \xdef\@currname{#3}%
 \global\let\@currext#4%
 \global\let\@clsextension\@currext
 \input{#3.#4}%
 \@ifl@ter#4{#3}#2{%
  \class@info{Class extension later than: #2}%
 }{%
  \class@info{Class extension earlier: #2}%
  \@@end
 }%
 #1%
 \let#1\@undefined
 \expandafter\@p@pfilename@ltx\@currnamestack@ltx\@nil
 \@reset@ptions
}%
\def\@pushfilename@ltx{%
 \xdef\@currnamestack@ltx{%
  {\@currname}%
  {\@currext}%
  {\@clsextension}%
  {\the\catcode`\@}%
  \@currnamestack@ltx
 }%
}%
\def\@p@pfilename@ltx#1#2#3#4#5\@nil{%
 \gdef\@currname{#1}%
 \gdef\@currext{#2}%
 \gdef\@clsextension{#3}%
 \catcode`\@#4\relax
 \gdef\@currnamestack@ltx{#5}%
}%
\global\let\@currnamestack@ltx\@empty
\def\flushing{%
  \let\\\@normalcr
  \leftskip\z@skip
  \rightskip\z@skip
  \@rightskip\z@skip
  \parfillskip\@flushglue
}%
\expandafter\DeclareRobustCommand\expandafter\@centercr\expandafter{\@centercr}%
\def\eqnarray@LaTeX{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip \tw@\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr
}
\long\def\eqnarray@fleqn@fixed{%
 \stepcounter{equation}\def\@currentlabel{\p@equation\theequation}%
 \global\@eqnswtrue\m@th\global\@eqcnt\z@
 \tabskip\ltx@mathindent
 \let\\=\@eqncr
 \setlength\abovedisplayskip{\topsep}%
 \ifvmode\addtolength\abovedisplayskip{\partopsep}\fi
 \addtolength\abovedisplayskip{\parskip}%
 \setlength\belowdisplayskip{\abovedisplayskip}%
 \setlength\belowdisplayshortskip{\abovedisplayskip}%
 \setlength\abovedisplayshortskip{\abovedisplayskip}%
 $$%
 \everycr{}%
 \halignt@\linewidth\bgroup
  \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
  &\global\@eqcnt\@ne
   \hskip\tw@\eqncolsep
   \hfil${{}##{}}$\hfil
  &\global\@eqcnt\tw@
   \hskip\tw@\eqncolsep
   $\displaystyle{##}$\hfil\tabskip\@centering
  &\global\@eqcnt\thr@@\hb@xt@\z@\bgroup\hss##\egroup
   \tabskip\z@skip
  \cr
}%
\@ifx{\eqnarray\eqnarray@LaTeX}{%
 \class@info{Repairing broken LaTeX eqnarray}%
 \let\eqnarray\eqnarray@fleqn@fixed
 \newlength\eqncolsep
 \setlength\eqncolsep\z@
 \let\eqnarray@LaTeX\relax
 \let\eqnarray@fleqn@fixed\relax
}{}%
\def\ltx@mathindent{\@centering}%
\def\set@eqnarray@skips{}%
\def\prep@math{%
 \@ifvmode{\everypar{{\setbox\z@\lastbox}}}{}%
}%
\def\prep@math@patch{%
 \prepdef\equation{\prep@math}%
 \prepdef\eqnarray{\prep@math}%
}%
\def\footnote{\@ifnextchar[\ltx@xfootnote\ltx@yfootnote}%
\def\ltx@xfootnote[#1]{%
 \ltx@def@footproc\ltx@footmark[#1]%
 \expandafter\ltx@foottext\expandafter{\the\csname c@\@mpfn\endcsname}%
}%
\def\ltx@yfootnote{%
 \ltx@stp@footproc\ltx@footmark
 \expandafter\ltx@foottext\expandafter{\the\csname c@\@mpfn\endcsname}%
}%
\def\footnotemark{\@ifnextchar[\ltx@xfootmark\ltx@yfootmark}%
\def\ltx@xfootmark{\ltx@def@footproc\ltx@footmark}%
\def\ltx@yfootmark{\ltx@stp@footproc\ltx@footmark}%
\def\ltx@footmark#1{%
 \leavevmode
 \ifhmode\edef\@x@sf{\the\spacefactor}\nobreak\fi
 \begingroup
  \expandafter\ltx@make@current@footnote\expandafter{\@mpfn}{#1}%
  \expandafter\@argswap@val\expandafter{\Hy@footnote@currentHref}{\hyper@linkstart {link}}%
   \@makefnmark
  \hyper@linkend
 \endgroup
 \ifhmode\spacefactor\@x@sf\fi
 \relax
}%
\def\footnotetext{\@ifnextchar[\ltx@xfoottext\ltx@yfoottext}%
\def\ltx@xfoottext{\ltx@def@footproc\ltx@foottext}%
\def\ltx@yfoottext{\ltx@stp@footproc\ltx@foottext}%
\long\def\ltx@foottext#1#2{%
 \begingroup
  \expandafter\ltx@make@current@footnote\expandafter{\@mpfn}{#1}%
  \@footnotetext{#2}%
 \endgroup
}%
\def\ltx@def@footproc#1[#2]{%
 \begingroup
   \csname c@\@mpfn\endcsname #2\relax
   \unrestored@protected@xdef\@thefnmark{\thempfn}%
 \expandafter\endgroup
 \expandafter#1%
 \expandafter{\the\csname c@\@mpfn\endcsname}%
}%
\def\ltx@stp@footproc#1{%
 \expandafter\stepcounter\expandafter{\@mpfn}%
 \protected@xdef\@thefnmark{\thempfn}%
 \expandafter#1%
 \expandafter{\the\csname c@\@mpfn\endcsname}%
}%
\appdef\class@documenthook{%
 \let\footnote@latex\footnote
 \@ifpackageloaded{hyperref}{}{%
  \let\H@@footnotetext\@footnotetext
  \def\@footnotetext{\H@@footnotetext}%
  \let\H@@mpfootnotetext\@mpfootnotetext
  \def\@mpfootnotetext{\H@@mpfootnotetext}%
 }%
}%
\def\ltx@make@current@footnote#1#2{%
  \csname c@#1\endcsname#2\relax
  \protected@edef\Hy@footnote@currentHref{\@currentHref-#1.\csname the#1\endcsname}%
}%
\def\thempfootnote@latex{{\itshape \@alph \c@mpfootnote }}%
\def\ltx@thempfootnote{\@alph\c@mpfootnote}%
\@ifx{\thempfootnote\thempfootnote@latex}{%
 \class@info{Repairing hyperref-unfriendly LaTeX definition of \string\mpfootnote}%
 \let\thempfootnote\ltx@thempfootnote
}{}%
\def\@makefnmark{%
 \hbox{%
  \@textsuperscript{%
   \normalfont\itshape\@thefnmark
  }%
 }%
}%
\long\def\@footnotetext{%
 \insert\footins\bgroup
  \make@footnotetext
}%
\long\def\@mpfootnotetext{%
 \minipagefootnote@pick
  \make@footnotetext
}%
\long\def\make@footnotetext#1{%
  \set@footnotefont
  \set@footnotewidth
  \@parboxrestore
  \protected@edef\@currentlabel{%
   \csname p@\@mpfn\endcsname\@thefnmark
  }%
  \color@begingroup
   \@makefntext{%
    \rule\z@\footnotesep\ignorespaces#1%
    \@finalstrut\strutbox\vadjust{\vskip\z@skip}%
   }%
  \color@endgroup
 \minipagefootnote@drop
}%
\def\set@footnotefont{%
  \reset@font\footnotesize
  \interlinepenalty\interfootnotelinepenalty
  \splittopskip\footnotesep
  \splitmaxdepth\dp\strutbox
}%
\def\set@footnotewidth{\set@footnotewidth@one}%
\def\robustify@contents{%
 \let \label \@gobble
 \let \index \@gobble
 \let \glossary \@gobble
 \let\footnote \@gobble
 \def\({\string\(}%
 \def\){\string\)}%
 \def\\{\string\\}%
}%
\long\def\addtocontents#1#2{%
 \protected@write\@auxout{\robustify@contents}{\string \@writefile {#1}{#2}}%
}%
\def\addcontentsline#1#2#3{%
 \addtocontents{#1}{%
  \protect\contentsline{#2}{#3}{\thepage}{}%
 }%
}%
\def\label#1{%
 \@bsphack
  \protected@write\@auxout{}{%
   \string\newlabel{#1}{{\@currentlabel}{\thepage}{}{}{}}%
  }%
 \@esphack
}%
\def\ltx@contentsline#1{%
 \expandafter\@ifnotrelax\csname l@#1\endcsname{}{%
  \expandafter\let\csname l@#1\endcsname\@gobbletwo
 }%
 \contentsline@latex{#1}%
}%
\appdef\document@inithook{%
 \let\contentsline@latex\contentsline
 \let\contentsline\ltx@contentsline
}%
\appdef\class@documenthook{%
 \prepdef\caption{\minipagefootnote@here}%
}%
\def\minipagefootnote@init{%
 \setbox\@mpfootins\box\voidb@x
}%
\def\minipagefootnote@pick{%
 \global\setbox\@mpfootins\vbox\bgroup
  \unvbox\@mpfootins
}%
\def\minipagefootnote@drop{%
 \egroup
}%
\def\minipagefootnote@here{%
    \par
    \@ifvoid\@mpfootins{}{%
      \vskip\skip\@mpfootins
      \fullinterlineskip
      \@ifinner{%
       \vtop{\unvcopy\@mpfootins}%
       {\setbox\z@\lastbox}%
      }{}%
      \unvbox\@mpfootins
    }%
}%
\def\minipagefootnote@foot{%
 \@ifvoid\@mpfootins{}{%
  \insert\footins\bgroup\unvbox\@mpfootins\egroup
 }%
}%
\def\endminipage{%
    \par
    \unskip
    \minipagefootnote@here
    \@minipagefalse   %% added 24 May 89
  \color@endgroup
  \egroup
  \expandafter\@iiiparbox\@mpargs{\unvbox\@tempboxa}%
}%
\@booleantrue\floats@sw
\let\@xfloat@LaTeX\@xfloat
\def\@xfloat#1[#2]{%
  \@xfloat@prep
  \@nameuse{fp@proc@#2}%
  \floats@sw{\@xfloat@LaTeX{#1}[#2]}{\@xfloat@anchored{#1}[]}%
}%
\def\@xfloat@prep{%
  \ltx@footnote@pop
  \def\@mpfn{mpfootnote}%
  \def\thempfn{\thempfootnote}%
  \c@mpfootnote\z@
  \let\H@@footnotetext\H@@mpfootnotetext
}%
\let\ltx@footnote@pop\@empty
\def\@xfloat@anchored#1[#2]{%
 \def\@captype{#1}%
 \begin@float@pagebreak
  \let\end@float\end@float@anchored
  \let\end@dblfloat\end@float@anchored
        \hsize\columnwidth
        \@parboxrestore
        \@floatboxreset
  \minipagefootnote@init
}%
\def\end@float@anchored{%
  \minipagefootnote@here
  \par\vskip\z@skip
 \par
 \end@float@pagebreak
}%
\def\begin@float@pagebreak{\par\addvspace\intextsep}%
\def\end@float@pagebreak{\par\addvspace\intextsep}%
\def\@mpmakefntext#1{%
 \parindent=1em
 \noindent
 \hb@xt@1em{\hss\@makefnmark}%
 #1%
}%
\def\do@if@floats#1#2{%
 \floats@sw{}{%
  \expandafter\newwrite
              \csname#1write\endcsname
  \expandafter\def
              \csname#1@stream\endcsname{\jobname#2}%
  \expandafter\immediate
  \expandafter\openout
              \csname#1write\endcsname
              \csname#1@stream\endcsname\relax
  \@ifxundefined\@float@LaTeX{%
   \let\@float@LaTeX\@float
   \let\@dblfloat@LaTeX\@dblfloat
   \let\@float\write@float
   \let\@dblfloat\write@floats
  }{}%
  \let@environment{#1@float}{#1}%
  \let@environment{#1@floats}{#1*}%
  \@ifxundefined@cs{#1@write}{}{%
   \let@environment{#1}{#1@write}%
  }%
 }%
}%
\def\triggerpar{\leavevmode\@@par}%
\def\oneapage{\def\begin@float@pagebreak{\newpage}\def\end@float@pagebreak{\newpage}}%
\def\print@float#1#2{%
 \lengthcheck@sw{%
  \total@float{#1}%
 }{}%
 \@ifxundefined@cs{#1write}{}{%
  \begingroup
   \@booleanfalse\floats@sw
   #2%
   \raggedbottom
   \def\array@default{v}% floats must
   \let\@float\@float@LaTeX
   \let\@dblfloat\@dblfloat@LaTeX
   \let\trigger@float@par\triggerpar
   \let@environment{#1}{#1@float}%
   \let@environment{#1*}{#1@floats}%
   \expandafter\prepdef\csname#1\endcsname{\trigger@float@par}%
   \expandafter\prepdef\csname#1*\endcsname{\trigger@float@par}%
   \@namedef{fps@#1}{h!}%
   \expandafter\immediate
   \expandafter\closeout
               \csname#1write\endcsname
   \everypar{%
    \global\let\trigger@float@par\relax
    \global\everypar{}\setbox\z@\lastbox
    \@ifxundefined@cs{#1sname}{}{%
     \begin@float@pagebreak
     \expandafter\section
     \expandafter*%
     \expandafter{%
                  \csname#1sname\endcsname
                 }%
    }%
   }%
   \input{\csname#1@stream\endcsname}%
  \endgroup
  \global\expandafter\let\csname#1write\endcsname\relax
 }%
}%
\chardef\@xvi=16\relax
\mathchardef\@twopowerfourteen="4000
\mathchardef\@twopowertwo="4
\def\tally@float#1{%
 \begingroup
  \@tempcnta\count\@currbox
  \divide\@tempcnta\@xxxii
  \multiply\@tempcnta\@xxxii
  \advance\count\@currbox-\@tempcnta
  \divide\@tempcnta\@xxxii
  \@ifnum{\count\@currbox>\@xvi}{%
   \advance\count\@currbox-\@xvi\@booleantrue\@temp@sw
  }{%
   \@booleanfalse\@temp@sw
  }%
  \show@box@size@sw{%
   \class@info{Float #1
    (\the\@tempcnta)[\@temp@sw{16+}{}\the\count\@currbox]^^J%
    (\the\ht\@currbox+\the\dp\@currbox)X\the\wd\@currbox
   }%
  }{}%
 \endgroup
 \expandafter\let
 \expandafter\@tempa
             \csname fbox@\csname ftype@#1\endcsname\endcsname
 \@ifnotrelax\@tempa{%
  \@ifhbox\@tempa{%
   \setbox\@tempboxa\vbox{\unvcopy\@currbox\hrule}%
   \dimen@\ht\@tempboxa
   \divide\dimen@\@twopowerfourteen
   \@ifdim{\wd\@tempboxa<\textwidth}{%
    \advance\dimen@\ht\@tempa
    \global\ht\@tempa\dimen@
   }{%
    \advance\dimen@\dp\@tempa
    \global\dp\@tempa\dimen@
   }%
  }{}%
 }{}%
}%
\def\total@float#1{%
 \expandafter\let
 \expandafter\@tempa
             \csname fbox@\csname ftype@#1\endcsname\endcsname
 \@ifnotrelax\@tempa{%
  \@ifhbox\@tempa{%
   \@tempdima\the\ht\@tempa\divide\@tempdima\@twopowertwo\@tempcnta\@tempdima
   \@tempdimb\the\dp\@tempa\divide\@tempdimb\@twopowertwo\@tempcntb\@tempdimb
   \class@info{Total #1: Column(\the\@tempcnta pt), Page(\the\@tempcnta pt)}%
  }{}%
 }{}%
}%
\def\write@float#1{\write@@float{#1}{#1}}%
\def\endwrite@float{\@Esphack}%
\def\write@floats#1{\write@@float{#1*}{#1}}%
\def\endwrite@floats{\@Esphack}%
\def\write@@float#1#2{%
  \ifhmode
     \@bsphack
  \fi
  \chardef\@tempc\csname#2write\endcsname
  \toks@{\begin{#1}}%
  \def\@tempb{#1}%
  \expandafter\let\csname end#1\endcsname\endwrite@float
  \catcode`\^^M\active
  \@makeother\{\@makeother\}\@makeother\%
  \write@floatline
}%
\begingroup
 \catcode`\[\the\catcode`\{\catcode`\]\the\catcode`\}\@makeother\{\@makeother\}%
 \gdef\float@end@tag#1\end{#2}#3\@nul[%
  \def\@tempa[#2]%
  \@ifx[\@tempa\@tempb][\end[#2]][\write@floatline]%
 ]%
 \obeylines%
 \gdef\write@floatline#1^^M[%
  \begingroup%
   \newlinechar`\^^M%
   \toks@\expandafter[\the\toks@#1]\immediate\write\@tempc[\the\toks@]%
  \endgroup%
  \toks@[]%
  \float@end@tag#1\end{}\@nul%
 ]%
\endgroup
\def\@alph#1{\ifcase#1\or a\or b\or c\or d\else\@ialph{#1}\fi}
\def\@ialph#1{\ifcase#1\or \or \or \or \or e\or f\or g\or h\or i\or j\or
  k\or l\or m\or n\or o\or p\or q\or r\or s\or t\or u\or v\or w\or x\or
  y\or z\or aa\or bb\or cc\or dd\or ee\or ff\or gg\or hh\or ii\or jj\or
  kk\or ll\or mm\or nn\or oo\or pp\or qq\or rr\or ss\or tt\or uu\or
  vv\or ww\or xx\or yy\or zz\else\@ctrerr\fi}
\def\@startsection#1#2#3#4#5#6{%
 \@startsection@hook
 \if@noskipsec \leavevmode \fi
 \par
 \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@
  \@tempskipa -\@tempskipa \@afterindentfalse
 \fi
 \if@nobreak
  \everypar{}%
 \else
  \addpenalty\@secpenalty\addvspace\@tempskipa
 \fi
 \@ifstar
  {\@dblarg{\@ssect@ltx{#1}{#2}{#3}{#4}{#5}{#6}}}%
  {\@dblarg{\@sect@ltx {#1}{#2}{#3}{#4}{#5}{#6}}}%
}%
\def\@startsection@hook{}%
\class@info{Repairing broken LateX \string\@sect}%
\def\@sect@ltx#1#2#3#4#5#6[#7]#8{%
  \@ifnum{#2>\c@secnumdepth}{%
    \def\H@svsec{\phantomsection}%
    \let\@svsec\@empty
  }{%
    \H@refstepcounter{#1}%
    \def\H@svsec{%
     \phantomsection
    }%
    \protected@edef\@svsec{{#1}}%
    \@ifundefined{@#1cntformat}{%
     \prepdef\@svsec\@seccntformat
    }{%
     \expandafter\prepdef
     \expandafter\@svsec
                 \csname @#1cntformat\endcsname
    }%
  }%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfrom@#1}{\@hang@from}{\csname @hangfrom@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{\@svsec}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{%
      \@ifnum{#2>\c@secnumdepth}{%
       \protect\numberline{}%
      }{%
       \protect\numberline{\csname the#1\endcsname}%
      }%
      #8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@to@#1}{\@runin@to}{\csname @runin@to@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{\@svsec}{#8}%
      }%
      \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{%
        \@ifnum{#2>\c@secnumdepth}{%
         \protect\numberline{}%
        }{%
         \protect\numberline{\csname the#1\endcsname}%
        }%
        #8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\@hang@from#1#2#3{\@hangfrom{#1#2}#3}%
\def\@runin@to #1#2#3{#1#2#3}%
\def\@ssect@ltx#1#2#3#4#5#6[#7]#8{%
  \def\H@svsec{\phantomsection}%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfroms@#1}{\@hang@froms}{\csname @hangfroms@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{\protect\numberline{}#8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@tos@#1}{\@runin@tos}{\csname @runin@tos@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{#8}%
      }%
      \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{\protect\numberline{}#8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\@hang@froms#1#2{#1#2}%
\def\@runin@tos #1#2{#1#2}%
\def\init@hyperref{%
 \providecommand\phantomsection{}%
 \providecommand\hyper@makecurrent[1]{}%
 \providecommand\Hy@raisedlink[1]{}%
 \providecommand\hyper@anchorstart[1]{}%
 \providecommand\hyper@anchorend{}%
 \providecommand\hyper@linkstart[2]{}%
 \providecommand\hyper@linkend{}%
 \providecommand\@currentHref{}%
}%
\let\H@refstepcounter\refstepcounter
\appdef\document@inithook{%
 \init@hyperref
}%
\def\sec@upcase#1{\relax{#1}}%
\appdef\document@inithook{%
 \@ifpackageloaded{array}{\switch@array}{\switch@tabular}%
 \prepdef\endtabular{\endtabular@hook}%
 \@provide\endtabular@hook{}%
 \prepdef\endarray{\endarray@hook}%
 \@provide\endarray@hook{}%
 \providecommand\array@hook{}%
 \prepdef\@tabular{\tabular@hook}%
 \@provide\tabular@hook{}%
}%
\def\switch@tabular{%
 \let\@array@sw\@array@sw@array
 \@ifx{\@array\@array@LaTeX}{%
  \@ifx{\multicolumn\multicolumn@LaTeX}{%
   \@ifx{\@tabular\@tabular@LaTeX}{%
    \@ifx{\@tabarray\@tabarray@LaTeX}{%
     \@ifx{\array\array@LaTeX}{%
      \@ifx{\endarray\endarray@LaTeX}{%
       \@ifx{\endtabular\endtabular@LaTeX}{%
        \@ifx{\@mkpream\@mkpream@LaTeX}{%
         \@ifx{\@addamp\@addamp@LaTeX}{%
          \@ifx{\@arrayacol\@arrayacol@LaTeX}{%
           \@ifx{\@tabacol\@tabacol@LaTeX}{%
            \@ifx{\@arrayclassz\@arrayclassz@LaTeX}{%
             \@ifx{\@tabclassiv\@tabclassiv@LaTeX}{%
              \@ifx{\@arrayclassiv\@arrayclassiv@LaTeX}{%
               \@ifx{\@tabclassz\@tabclassz@LaTeX}{%
                \@ifx{\@classv\@classv@LaTeX}{%
                 \@ifx{\hline\hline@LaTeX}{%
                  \@ifx{\@tabularcr\@tabularcr@LaTeX}{%
                   \@ifx{\@xtabularcr\@xtabularcr@LaTeX}{%
                    \@ifx{\@xargarraycr\@xargarraycr@LaTeX}{%
                     \@ifx{\@yargarraycr\@yargarraycr@LaTeX}{%
                      \true@sw
                     }{%
                      \false@sw
                     }%
                    }{%
                     \false@sw
                    }%
                   }{%
                    \false@sw
                   }%
                  }{%
                   \false@sw
                  }%
                 }{%
                  \false@sw
                 }%
                }{%
                 \false@sw
                }%
               }{%
                \false@sw
               }%
              }{%
               \false@sw
              }%
             }{%
              \false@sw
             }%
            }{%
             \false@sw
            }%
           }{%
            \false@sw
           }%
          }{%
           \false@sw
          }%
         }{%
          \false@sw
         }%
        }{%
         \false@sw
        }%
       }{%
        \false@sw
       }%
      }{%
       \false@sw
      }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }{%
    \false@sw
   }%
  }{%
   \false@sw
  }%
 }{%
  \false@sw
 }%
 {%
  \class@info{Patching LaTeX tabular.}%
 }{%
  \class@info{Unrecognized LaTeX tabular. Please update this document class! (Proceeding with fingers crossed.)}%
 }%
 \let\@array\@array@ltx
 \let\multicolumn\multicolumn@ltx
 \let\@tabular\@tabular@ltx
 \let\@tabarray\@tabarray@ltx
 \let\array\array@ltx
 \let\endarray\endarray@ltx
 \let\endtabular\endtabular@ltx
 \let\@mkpream\@mkpream@ltx
 \let\@addamp\@addamp@ltx
 \let\@arrayacol\@arrayacol@ltx
 \let\@tabacol\@tabacol@ltx
 \let\@arrayclassz\@arrayclassz@ltx
 \let\@tabclassiv\@tabclassiv@ltx
 \let\@arrayclassiv\@arrayclassiv@ltx
 \let\@tabclassz\@tabclassz@ltx
 \let\@classv\@classv@ltx
 \let\hline\hline@ltx
 \let\@tabularcr\@tabularcr@ltx
 \let\@xtabularcr\@xtabularcr@ltx
 \let\@xargarraycr\@xargarraycr@ltx
 \let\@yargarraycr\@yargarraycr@ltx
}%
\def\switch@array{%
 \@ifpackageloaded{colortbl}{\let\switch@array@info\colortbl@message}{\let\switch@array@info\array@message}%
 \let\@array@sw\@array@sw@LaTeX
 \@ifx{\@array\@array@array}{%
  \@ifx{\@tabular\@tabular@array}{%
   \@ifx{\@tabarray\@tabarray@array}{%
    \@ifx{\array\array@array}{%
     \@ifx{\endarray\endarray@array}{%
      \@ifx{\endtabular\endtabular@array}{%
       \@ifx{\@mkpream\@mkpream@array}{%
        \@ifx{\@classx\@classx@array}{%
         \@ifx{\insert@column\insert@column@array}{%
          \@ifx{\@arraycr\@arraycr@array}{%
           \@ifx{\@xarraycr\@xarraycr@array}{%
            \@ifx{\@xargarraycr\@xargarraycr@array}{%
             \@ifx{\@yargarraycr\@yargarraycr@array}{%
              \true@sw
             }{%
              \false@sw
             }%
            }{%
             \false@sw
            }%
           }{%
            \false@sw
           }%
          }{%
           \false@sw
          }%
         }{%
          \false@sw
         }%
        }{%
         \false@sw
        }%
       }{%
        \false@sw
       }%
      }{%
       \false@sw
      }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }{%
    \false@sw
   }%
  }{%
   \false@sw
  }%
 }{%
  \false@sw
 }{%
  \class@info{Patching array package.}%
 }{%
  \switch@array@info
 }%
 \let\@array    \@array@array@new
 \let\@@array   \@array % Cosi fan tutti
 \let\@tabular  \@tabular@array@new
 \let\@tabarray \@tabarray@array@new
 \let\array     \array@array@new
 \let\endarray  \endarray@array@new
 \let\endtabular\endtabular@array@new
 \let\@mkpream  \@mkpream@array@new
 \let\@classx   \@classx@array@new
 \let\@arrayacol\@arrayacol@ltx
 \let\@tabacol  \@tabacol@ltx
 \let\insert@column\insert@column@array@new
 \expandafter\let\csname endtabular*\endcsname\endtabular % Cosi fan tutti
 \let\@arraycr  \@arraycr@new
 \let\@xarraycr \@xarraycr@new
 \let\@xargarraycr\@xargarraycr@new
 \let\@yargarraycr\@yargarraycr@new
}%
\def\array@message{%
 \class@info{Unrecognized array package. Please update this document class! (Proceeding with fingers crossed.)}%
}%
\def\colortbl@message{%
 \class@info{colortbl package is loaded. (Proceeding with fingers crossed.)}%
}%
\def\@array@sw@LaTeX{\@ifx{\\\@tabularcr}}%
\def\@array@sw@array{\@ifx{\d@llarbegin\begingroup}}%
\def\@tabular@LaTeX{%
 \leavevmode
 \hbox\bgroup$%
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \let\\\@tabularcr
  \@tabarray
}%
\def\@tabular@ltx{%
  \let\@acoll\@tabacoll
  \let\@acolr\@tabacolr
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \let\\\@tabularcr
  \@tabarray
}%
\def\@tabular@array{%
 \leavevmode
 \hbox\bgroup$%
  \col@sep\tabcolsep
  \let\d@llarbegin\begingroup
  \let\d@llarend\endgroup
  \@tabarray
}%
\def\@tabular@array@new{%
  \let\@acoll\@tabacoll
  \let\@acolr\@tabacolr
  \let\@acol\@tabacol
  \let\d@llarbegin\begingroup
  \let\d@llarend\endgroup
  \@tabarray
}%
\def\@tabarray@LaTeX{%
 \m@th\@ifnextchar[\@array{\@array[c]}%
}%
\def\@tabarray@ltx{%
 \m@th\@ifnextchar[\@array{\expandafter\@array\expandafter[\array@default]}%
}%
\def\@tabarray@array{%
 \@ifnextchar[{\@@array}{\@@array[c]}%
}%
\def\@tabarray@array@new{%
 \@ifnextchar[{\@@array}{\expandafter\@@array\expandafter[\array@default]}%
}%
\newcount\intertabularlinepenalty
\intertabularlinepenalty=100
\newcount\@tbpen
\appdef\samepage{\intertabularlinepenalty\@M}%
\def\@tabularcr@LaTeX{{\ifnum 0=`}\fi \@ifstar \@xtabularcr \@xtabularcr}%
\def\@tabularcr@ltx{{\ifnum 0=`}\fi \@ifstar {\global \@tbpen \@M \@xtabularcr }{\global \@tbpen \intertabularlinepenalty \@xtabularcr }}%
\def\@xtabularcr@LaTeX{\@ifnextchar [\@argtabularcr {\ifnum 0=`{\fi }\cr }}%
\def\@xtabularcr@ltx{\@ifnextchar [\@argtabularcr {\ifnum 0=`{\fi }\cr \noalign {\penalty \@tbpen }}}%
\def\@xargarraycr@LaTeX#1{\@tempdima #1\advance \@tempdima \dp \@arstrutbox \vrule \@height \z@ \@depth \@tempdima \@width \z@ \cr}%
\def\@xargarraycr@ltx#1{\@tempdima #1\advance \@tempdima \dp \@arstrutbox \vrule \@height \z@ \@depth \@tempdima \@width \z@ \cr \noalign {\penalty \@tbpen }}%
\def\@yargarraycr@LaTeX#1{\cr \noalign {\vskip #1}}%
\def\@yargarraycr@ltx#1{\cr \noalign {\penalty \@tbpen \vskip #1}}%
\def\@arraycr@array{%
 \relax
 \iffalse{\fi\ifnum 0=`}\fi
 \@ifstar \@xarraycr \@xarraycr
}%
\def\@arraycr@new{%
 \relax
 \iffalse{\fi\ifnum 0=`}\fi
 \@ifstar {\global \@tbpen \@M \@xarraycr }{\global \@tbpen \intertabularlinepenalty \@xarraycr }%
}%
\def\@xarraycr@array{%
 \@ifnextchar [%]
 \@argarraycr {\ifnum 0=`{}\fi\cr}%
}%
\def\@xarraycr@new{%
 \@ifnextchar [%]
 \@argarraycr {\ifnum 0=`{}\fi\cr \noalign {\penalty \@tbpen }}%
}%
\def\@xargarraycr@array#1{%
 \unskip
 \@tempdima #1\advance\@tempdima \dp\@arstrutbox
 \vrule \@depth\@tempdima \@width\z@
 \cr
}%
\def\@xargarraycr@new#1{%
 \unskip
 \@tempdima #1\advance\@tempdima \dp\@arstrutbox
 \vrule \@depth\@tempdima \@width\z@
 \cr
 \noalign {\penalty \@tbpen }%
}%
\def\@yargarraycr@array#1{%
 \cr
 \noalign{\vskip #1}%
}%
\def\@yargarraycr@new#1{%
 \cr
 \noalign{\penalty \@tbpen \vskip #1}%
}%
\def\array@LaTeX{%
 \let\@acol\@arrayacol
 \let\@classz\@arrayclassz
 \let\@classiv\@arrayclassiv
 \let\\\@arraycr
 \let\@halignto\@empty
 \@tabarray
}%
\def\array@ltx{%
 \@ifmmode{}{\@badmath$}%
 \let\@acoll\@arrayacol
 \let\@acolr\@arrayacol
 \let\@acol\@arrayacol
 \let\@classz\@arrayclassz
 \let\@classiv\@arrayclassiv
 \let\\\@arraycr
 \let\@halignto\@empty
 \@tabarray
}%
\def\array@array{%
 \col@sep\arraycolsep
 \def\d@llarbegin{$}\let\d@llarend\d@llarbegin\gdef\@halignto{}%
 \@tabarray
}
\def\array@array@new{%
 \@ifmmode{}{\@badmath$}%
 \let\@acoll\@arrayacol
 \let\@acolr\@arrayacol
 \let\@acol\@arrayacol
 \def\d@llarbegin{$}%
 \let\d@llarend\d@llarbegin
 \gdef\@halignto{}%
 \@tabarray
}%
\def\@array@LaTeX[#1]#2{%
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi\fi
  \bgroup
  \setbox\@arstrutbox\hbox{%
    \vrule \@height\arraystretch\ht\strutbox
           \@depth\arraystretch \dp\strutbox
           \@width\z@}%
  \@mkpream{#2}%
  \edef\@preamble{%
    \ialign \noexpand\@halignto
      \bgroup \@arstrut \@preamble \tabskip\z@skip \cr}%
  \let\@startpbox\@@startpbox \let\@endpbox\@@endpbox
  \let\tabularnewline\\%
    \let\par\@empty
    \let\@sharp##%
    \set@typeset@protect
    \lineskip\z@skip\baselineskip\z@skip
    \ifhmode \@preamerr\z@ \@@par\fi
    \@preamble
}%
\def\@array@ltx[#1]#2{%
 \@nameuse{@array@align@#1}%
  \set@arstrutbox
  \@mkpream{#2}%
  \prepdef\@preamble{%
    \tabskip\tabmid@skip
    \@arstrut
  }%
  \appdef\@preamble{%
    \tabskip\tabright@skip
    \cr
    \array@row@pre
  }%
  \let\tabularnewline\\%
  \let\par\@empty
  \let\@sharp##%
  \set@typeset@protect
  \lineskip\z@skip\baselineskip\z@skip
  \tabskip\tableft@skip\relax
  \ifhmode \@preamerr\z@ \@@par\fi
  \everycr{}%
  \expandafter\halign\expandafter\@halignto\expandafter\bgroup\@preamble
}%
\def\set@arstrutbox{%
  \setbox\@arstrutbox\hbox{%
    \vrule \@height\arraystretch\ht\strutbox
           \@depth\arraystretch \dp\strutbox
           \@width\z@
  }%
}%
\def\@array@array[#1]#2{%
  \@tempdima \ht \strutbox
  \advance \@tempdima by\extrarowheight
  \setbox \@arstrutbox \hbox{\vrule
             \@height \arraystretch \@tempdima
             \@depth \arraystretch \dp \strutbox
             \@width \z@}%
  \begingroup
  \@mkpream{#2}%
  \xdef\@preamble{\noexpand \ialign \@halignto
                  \bgroup \@arstrut \@preamble
                          \tabskip \z@ \cr}%
  \endgroup
  \@arrayleft
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi \fi
  \bgroup
  \let \@sharp ##\let \protect \relax
  \lineskip \z@
  \baselineskip \z@
  \m@th
  \let\\\@arraycr \let\tabularnewline\\\let\par\@empty \@preamble
}%
\def\@array@array@new[#1]#2{%
  \@tempdima\ht\strutbox
  \advance\@tempdima by\extrarowheight
  \setbox\@arstrutbox\hbox{%
   \vrule \@height\arraystretch\@tempdima
          \@depth \arraystretch\dp\strutbox
          \@width \z@
  }%
  \begingroup
   \@mkpream{#2}%
   \xdef\@preamble{\@preamble}%
  \endgroup
  \prepdef\@preamble{%
   \tabskip\tabmid@skip
    \@arstrut
  }%
  \appdef\@preamble{%
   \tabskip\tabright@skip
   \cr
   \array@row@pre
  }%
  \@arrayleft
  \@nameuse{@array@align@#1}%
  \m@th
  \let\\\@arraycr
  \let\tabularnewline\\%
  \let\par\@empty
  \let\@sharp##%
  \set@typeset@protect
  \lineskip\z@\baselineskip\z@
  \tabskip\tableft@skip
  \everycr{}%
  \expandafter\halign\expandafter\@halignto\expandafter\bgroup\@preamble
}%
\def\endarray@LaTeX{%
 \crcr\egroup\egroup
}%
\def\endarray@ltx{%
 \crcr\array@row@pst\egroup\egroup
}%
\def\endarray@array{%
 \crcr \egroup \egroup \@arrayright \gdef\@preamble{}%
}%
\def\endarray@array@new{%
 \crcr\array@row@pst\egroup\egroup % Same as \endarray@ltx
 \@arrayright
 \global\let\@preamble\@empty
}%
\def\endtabular@LaTeX{%
 \crcr\egroup\egroup $\egroup
}%
\def\endtabular@ltx{%
 \endarray
}%
\def\endtabular@array{%
 \endarray $\egroup
}%
\def\endtabular@array@new{%
 \endarray
}%
\@namedef{endtabular*}{\endtabular}%
\long\def\multicolumn@LaTeX#1#2#3{%
 \multispan{#1}\begingroup
  \@mkpream{#2}%
  \def\@sharp{#3}\set@typeset@protect
  \let\@startpbox\@@startpbox\let\@endpbox\@@endpbox
  \@arstrut \@preamble\hbox{}\endgroup\ignorespaces
}%
\long\def\multicolumn@ltx#1#2#3{%
 \multispan{#1}%
 \begingroup
  \@mkpream{#2}%
  \def\@sharp{#3}%
  \set@typeset@protect
 %\let\@startpbox\@@startpbox\let\@endpbox\@@endpbox
  \@arstrut
  \@preamble
  \hbox{}%
 \endgroup
 \ignorespaces
}%
\def\@array@align@t{\leavevmode\vtop\bgroup}%
\def\@array@align@b{\leavevmode\vbox\bgroup}%
\def\@array@align@c{\leavevmode\@ifmmode{\vcenter\bgroup}{$\vcenter\bgroup\aftergroup$\aftergroup\relax}}%
\def\@array@align@v{%
 \@ifmmode{%
  \@badmath
  \vcenter\bgroup
 }{%
  \@ifinner{%
   $\vcenter\bgroup\aftergroup$
  }{%
   \@@par\bgroup
  }%
 }%
}%
\def\array@default{c}%
\def\array@row@rst{%
 \let\@array@align@v\@array@align@c
}%
\def\array@row@pre{}%
\def\array@row@pst{}%
\newcommand\toprule{\tab@rule{\column@font}{\column@fil}{\frstrut}}%
\newcommand\colrule{\unskip\lrstrut\\\tab@rule{\body@font}{}{\frstrut}}%
\newcommand\botrule{\unskip\lrstrut\\\noalign{\hline@rule}{}}%
\def\hline@LaTeX{%
 \noalign{\ifnum0=`}\fi\hrule \@height \arrayrulewidth \futurelet
   \reserved@a\@xhline
}%
\def\hline@ltx{%
 \noalign{%
  \ifnum0=`}\fi
  \hline@rule
  \futurelet\reserved@a\@xhline
 % \noalign ended in \@xhline
}%
\def\@xhline@unneeded{%
 \say\reserved@a
 \ifx\reserved@a\hline
  \vskip\doublerulesep
  \vskip-\arrayrulewidth
 \fi
 \ifnum0=`{\fi}%
}%
\def\tab@rule#1#2#3{%
 \crcr
 \noalign{%
  \hline@rule
  \gdef\@arstrut@hook{%
   \global\let\@arstrut@hook\@empty
   #3%
  }%
  \gdef\cell@font{#1}%
  \gdef\cell@fil{#2}%
 }%
}%
\def\column@font{}%
\def\column@fil{}%
\def\body@font{}%
\def\cell@font{}%
\def\frstrut{}%
\def\lrstrut{}%
\def\@arstrut@hline{%
 \relax
 \@ifmmode{\copy}{\unhcopy}\@arstrutbox@hline
 \@arstrut@hook
}%
\let\@arstrut@org\@arstrut
\def\@arstrut@hook{%
 \global\let\@arstrut\@arstrut@org
}%
\newbox\@arstrutbox@hline
\appdef\set@arstrutbox{%
  \setbox\@arstrutbox@hline\hbox{%
    \setbox\z@\hbox{$0^{0}_{}$}%
    \dimen@\ht\z@\advance\dimen@\@arstrut@hline@clnc
    \@ifdim{\dimen@<\arraystretch\ht\strutbox}{\dimen@=\arraystretch\ht\strutbox}{}%
    \vrule \@height\dimen@
           \@depth\arraystretch \dp\strutbox
           \@width\z@
  }%
}%
\def\hline@rule{%
 \hrule \@height \arrayrulewidth
 \global\let\@arstrut\@arstrut@hline
}%
\def\@arstrut@hline@clnc{2\p@}% % Klootch: magic number
\def\tableft@skip{\z@skip}%
\def\tabmid@skip{\z@skip}%\@flushglue
\def\tabright@skip{\z@skip}%
\def\tableftsep{\tabcolsep}%
\def\tabmidsep{\tabcolsep}%
\def\tabrightsep{\tabcolsep}%
\def\cell@fil{}%
\def\pbox@hook{}%
\appdef\@arstrut{\@arstrut@hook}%
\let\@arstrut@hook\@empty
\def\@addtopreamble{\appdef\@preamble}%
\def\@mkpream@LaTeX#1{%
  \@firstamptrue\@lastchclass6
  \let\@preamble\@empty
  \let\protect\@unexpandable@protect
  \let\@sharp\relax
  \let\@startpbox\relax\let\@endpbox\relax
  \@expast{#1}%
  \expandafter\@tfor \expandafter
    \@nextchar \expandafter:\expandafter=\reserved@a\do
       {\@testpach\@nextchar
    \ifcase \@chclass \@classz \or \@classi \or \@classii \or \@classiii
      \or \@classiv \or\@classv \fi\@lastchclass\@chclass}%
  \ifcase \@lastchclass \@acol
      \or \or \@preamerr \@ne\or \@preamerr \tw@\or \or \@acol \fi
}%
\def\@mkpream@ltx#1{%
 \@firstamptrue
 \@lastchclass6
 \let\@preamble\@empty
 \let\protect\@unexpandable@protect
 \let\@sharp\relax
 \@expast{#1}%
 \expandafter\@tfor\expandafter\@nextchar\expandafter:\expandafter=\reserved@a
 \do{%
  \expandafter\@testpach\expandafter{\@nextchar}%
  \ifcase\@chclass
   \@classz
  \or
   \@classi
  \or
   \@classii
  \or
   \@classiii
  \or
   \@classiv
  \or
   \@classv
  \fi
  \@lastchclass\@chclass
 }%
 \ifcase\@lastchclass
  \@acolr % right-hand column
 \or
 \or
  \@preamerr\@ne
 \or
  \@preamerr\tw@
 \or
 \or
  \@acolr % right-hand column
 \fi
}%
\def\insert@column@array{%
   \the@toks \the \@tempcnta
   \ignorespaces \@sharp \unskip
   \the@toks \the \count@ \relax
}%
\def\insert@column@array@new{%
 \the@toks\the\@tempcnta
 \array@row@rst\cell@font
 \ignorespaces\@sharp\unskip
 \the@toks\the\count@
 \relax
}%
\def\@mkpream@relax{%
 \let\tableftsep   \relax
 \let\tabmidsep    \relax
 \let\tabrightsep  \relax
 \let\array@row@rst\relax
 \let\cell@font    \relax
 \let\@startpbox   \relax
}%
\def\@mkpream@array#1{%
   \gdef\@preamble{}\@lastchclass 4 \@firstamptrue
   \let\@sharp\relax \let\@startpbox\relax \let\@endpbox\relax
   \@temptokena{#1}\@tempswatrue
   \@whilesw\if@tempswa\fi{\@tempswafalse\the\NC@list}%
   \count@\m@ne
   \let\the@toks\relax
   \prepnext@tok
   \expandafter \@tfor \expandafter \@nextchar
    \expandafter :\expandafter =\the\@temptokena \do
   {\@testpach
   \ifcase \@chclass \@classz \or \@classi \or \@classii
     \or \save@decl \or \or \@classv \or \@classvi
     \or \@classvii \or \@classviii
     \or \@classx
     \or \@classx \fi
   \@lastchclass\@chclass}%
   \ifcase\@lastchclass
   \@acol \or
   \or
   \@acol \or
   \@preamerr \thr@@ \or
   \@preamerr \tw@ \@addtopreamble\@sharp \or
   \or
   \else  \@preamerr \@ne \fi
   \def\the@toks{\the\toks}%
}%
\def\@mkpream@array@new#1{%
 \gdef\@preamble{}%
 \@lastchclass\f@ur
 \@firstamptrue
 \let\@sharp\relax
 \@mkpream@relax
 \@temptokena{#1}\@tempswatrue
 \@whilesw\if@tempswa\fi{\@tempswafalse\the\NC@list}%
 \count@\m@ne
 \let\the@toks\relax
 \prepnext@tok
 \expandafter\@tfor\expandafter\@nextchar\expandafter:\expandafter=\the\@temptokena
 \do{%
  \@testpach
  \ifcase\@chclass
   \@classz
  \or
   \@classi
  \or
   \@classii
  \or
   \save@decl
  \or
  \or
   \@classv
  \or
   \@classvi
  \or
   \@classvii
  \or
   \@classviii
  \or
   \@classx
  \or
   \@classx
  \fi
  \@lastchclass\@chclass
 }%
 \ifcase\@lastchclass
  \@acolr % right-hand column
 \or
 \or
  \@acolr % right-hand column
 \or
  \@preamerr\thr@@
 \or
  \@preamerr\tw@\@addtopreamble\@sharp
 \or
 \or
 \else
  \@preamerr\@ne
 \fi
 \def\the@toks{\the\toks}%
}%
\appdef\@mkpream@relax{%
 \let\CT@setup       \relax
 \let\CT@color       \relax
 \let\CT@do@color    \relax
 \let\color          \relax
 \let\CT@column@color\relax
 \let\CT@row@color   \relax
 \let\CT@cell@color  \relax
}%
\def\@addamp@LaTeX{%
  \if@firstamp\@firstampfalse\else\edef\@preamble{\@preamble &}\fi
}%
\def\@addamp@ltx{%
 \if@firstamp\@firstampfalse\else\@addtopreamble{&}\fi
}%
\def\@arrayacol@LaTeX{%
 \edef\@preamble{\@preamble \hskip \arraycolsep}%
}%
\def\@arrayacol@ltx{%
 \@addtopreamble{\hskip\arraycolsep}%
}%
\def\@tabacoll{%
 \@addtopreamble{\hskip\tableftsep\relax}%
}%
\def\@tabacol@LaTeX{%
 \edef\@preamble{\@preamble \hskip \tabcolsep}%
}%
\def\@tabacol@ltx{%
 \@addtopreamble{\hskip\tabmidsep\relax}%
}%
\def\@tabacolr{%
 \@addtopreamble{\hskip\tabrightsep\relax}%
}%
\def\@arrayclassz@LaTeX{%
 \ifcase \@lastchclass \@acolampacol \or \@ampacol \or
   \or \or \@addamp \or
   \@acolampacol \or \@firstampfalse \@acol \fi
 \edef\@preamble{\@preamble
  \ifcase \@chnum
     \hfil$\relax\@sharp$\hfil \or $\relax\@sharp$\hfil
    \or \hfil$\relax\@sharp$\fi}%
}%
\def\@arrayclassz@ltx{%
 \ifcase\@lastchclass
  \@acolampacol
 \or
  \@ampacol
 \or
 \or
 \or
  \@addamp
 \or
  \@acolampacol
 \or
  \@firstampfalse\@acoll
 \fi
 \ifcase\@chnum
  \@addtopreamble{%
   \hfil\array@row@rst$\relax\@sharp$\hfil
  }%
 \or
  \@addtopreamble{%
   \array@row@rst$\relax\@sharp$\hfil
  }%
 \or
  \@addtopreamble{%
   \hfil\array@row@rst$\relax\@sharp$%
  }%
 \fi
}%
\def\@tabclassz@LaTeX{%
  \ifcase\@lastchclass
    \@acolampacol
  \or
    \@ampacol
  \or
  \or
  \or
    \@addamp
  \or
    \@acolampacol
  \or
    \@firstampfalse\@acol
  \fi
  \edef\@preamble{%
    \@preamble{%
      \ifcase\@chnum
        \hfil\ignorespaces\@sharp\unskip\hfil
      \or
        \hskip1sp\ignorespaces\@sharp\unskip\hfil
      \or
        \hfil\hskip1sp\ignorespaces\@sharp\unskip
      \fi}}%
}%
\def\@tabclassz@ltx{%
 \ifcase\@lastchclass
  \@acolampacol
 \or
  \@ampacol
 \or
 \or
 \or
  \@addamp
 \or
  \@acolampacol
 \or
  \@firstampfalse\@acoll
 \fi
 \ifcase\@chnum
  \@addtopreamble{%
   {\hfil\array@row@rst\cell@font\ignorespaces\@sharp\unskip\hfil}%
  }%
 \or
  \@addtopreamble{%
   {\cell@fil\hskip1sp\array@row@rst\cell@font\ignorespaces\@sharp\unskip\hfil}%
  }%
 \or
  \@addtopreamble{%
   {\hfil\hskip1sp\array@row@rst\cell@font\ignorespaces\@sharp\unskip\cell@fil}%
  }%
 \fi
}%
\def\@tabclassiv@LaTeX{%
 \@addtopreamble\@nextchar
}%
\def\@tabclassiv@ltx{%
 \expandafter\@addtopreamble\expandafter{\@nextchar}%
}%
\def\@arrayclassiv@LaTeX{%
 \@addtopreamble{$\@nextchar$}%
}%
\def\@arrayclassiv@ltx{%
 \expandafter\@addtopreamble\expandafter{\expandafter$\@nextchar$}%
}%
\def\@classv@LaTeX{%
 \@addtopreamble{\@startpbox{\@nextchar}\ignorespaces
 \@sharp\@endpbox}%
}%
\def\@classv@ltx{%
 \expandafter\@addtopreamble
 \expandafter{%
 \expandafter \@startpbox
 \expandafter {\@nextchar}%
 \pbox@hook\array@row@rst\cell@font\ignorespaces\@sharp\@endpbox
 }%
}%
\def\@classx@array{%
  \ifcase \@lastchclass
  \@acolampacol \or
  \@addamp \@acol \or
  \@acolampacol \or
  \or
  \@acol \@firstampfalse \or
  \@addamp
  \fi
}%
\def\@classx@array@new{%
 \ifcase \@lastchclass
  \@acolampacol
 \or
  \@addamp \@acol
 \or
  \@acolampacol
 \or
 \or
  \@firstampfalse\@acoll
 \or
  \@addamp
 \fi
}%
\def\@xbitor@LaTeX #1{\@tempcntb \count#1
   \ifnum \@tempcnta =\z@
   \else
     \divide\@tempcntb\@tempcnta
     \ifodd\@tempcntb \@testtrue\fi
   \fi}%
\def\@xbitor@ltx#1{%
 \@tempcntb\count#1\relax
 \@ifnum{\@tempcnta=\z@}{}{%
  \divide\@tempcntb\@tempcnta
  \@ifodd\@tempcntb{\@testtrue}{}%
 }%
}%
\@ifx{\@xbitor\@xbitor@LaTeX}{%
  \class@info{Repairing broken LaTeX \string\@xbitor}%
}{%
  \class@info{Unrecognized LaTeX \string\@xbitor. Please update this document class! (Proceeding with fingers crossed.)}%
}%
\let\@xbitor\@xbitor@ltx
\newcommand*\@gobble@opt@one[2][]{}%
\def\@starttoc#1{%
  \begingroup
    \toc@pre
    \makeatletter
    \@input{\jobname.#1}%
    \if@filesw
      \expandafter\newwrite\csname tf@#1\endcsname
      \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax
    \fi
    \@nobreakfalse
    \toc@post
  \endgroup
}%
\def\toc@pre{}%
\def\toc@post{}%
\def\toc@@font{}%
\def\ltxu@dotsep{\z@}%
\let\tocdim@section       \leftmargini
\let\tocdim@subsection    \leftmarginii
\let\tocdim@subsubsection \leftmarginiii
\let\tocdim@paragraph     \leftmarginiv
\let\tocdim@appendix      \leftmarginv
\let\tocdim@pagenum       \leftmarginvi
\def\toc@pre@auto{%
  \toc@@font
  \@tempdima\z@
  \toc@setindent\@tempdima{section}%
  \toc@setindent\@tempdima{subsection}%
  \toc@setindent\@tempdima{subsubsection}%
  \toc@setindent\@tempdima{paragraph}%
  \toc@letdimen{appendix}%
  \toc@letdimen{pagenum}%
}%
\def\toc@post@auto{%
  \if@filesw
   \begingroup
    \toc@writedimen{section}%
    \toc@writedimen{subsection}%
    \toc@writedimen{subsubsection}%
    \toc@writedimen{paragraph}%
    \toc@writedimen{appendix}%
    \toc@writedimen{pagenum}%
   \endgroup
  \fi
}%
\def\toc@setindent#1#2{%
 \csname tocdim@#2\endcsname\tocdim@min\relax
 \@ifundefined{tocmax@#2}{\@namedef{tocmax@#2}{\z@}}{}%
 \advance#1\@nameuse{tocmax@#2}\relax
 \expandafter\edef\csname tocleft@#2\endcsname{\the#1}%
}%
\def\toc@letdimen#1{%
 \csname tocdim@#1\endcsname\tocdim@min\relax
 \@ifundefined{tocmax@#1}{\@namedef{tocmax@#1}{\z@}}{}%
 \expandafter\let\csname tocleft@#1\expandafter\endcsname\csname tocmax@#1\endcsname
}%
\def\toc@writedimen#1{%
 \immediate\write\@auxout{%
  \gdef\expandafter\string\csname tocmax@#1\endcsname{%
   \expandafter\the\csname tocdim@#1\endcsname
  }%
 }%
}%
\def\l@@sections#1#2#3#4{%
 \begingroup
  \everypar{}%
  \set@tocdim@pagenum\@tempboxa{#4}%
  \global\@tempdima\csname tocdim@#2\endcsname
  \leftskip\csname tocleft@#2\endcsname\relax
  \dimen@\csname tocleft@#1\endcsname\relax
  \parindent-\leftskip\advance\parindent\dimen@
  \rightskip\tocleft@pagenum plus 1fil\relax
  \skip@\parfillskip\parfillskip\z@
  \let\numberline\numberline@@sections
  \@nameuse{l@f@#2}%
  \ignorespaces#3\unskip\nobreak\hskip\skip@
  \hb@xt@\rightskip{\hfil\unhbox\@tempboxa}\hskip-\rightskip\hskip\z@skip
  \expandafter\par
  \expandafter\aftergroup\csname tocdim@#2%
  \expandafter\endcsname
  \expandafter\endgroup
              \the\@tempdima\relax
}%
\def\set@tocdim@pagenum#1#2{%
 \setbox#1\hbox{\ignorespaces#2}%
 \@ifdim{\tocdim@pagenum<\wd#1}{\global\tocdim@pagenum\wd#1}{}%
}%
\def\numberline@@sections#1{%
 \leavevmode\hb@xt@-\parindent{%
  \hfil
  \@if@empty{#1}{}{%
   \setbox\z@\hbox{#1.\kern\ltxu@dotsep}%
   \@ifdim{\@tempdima<\wd\z@}{\global\@tempdima\wd\z@}{}%
   \unhbox\z@
  }%
 }%
 \ignorespaces
}%
\def\tocdim@min{\z@}%
\def\list#1#2{%
  \ifnum \@listdepth >5\relax
    \@toodeep
  \else
    \global\advance\@listdepth\@ne
  \fi
  \rightmargin\z@
  \listparindent\z@
  \itemindent\z@
  \csname @list\romannumeral\the\@listdepth\endcsname
  \def\@itemlabel{#1}%
  \let\makelabel\@mklab
  \@nmbrlistfalse
  #2\relax
  \@trivlist
  \parskip\parsep
  \set@listindent
  \ignorespaces
}%
\def\set@listindent@parshape{%
 \parindent\listparindent
 \advance\@totalleftmargin\leftmargin
 \advance\linewidth-\rightmargin
 \advance\linewidth-\leftmargin
 \parshape\@ne\@totalleftmargin\linewidth
}%
\def\set@listindent@{%
 \parindent\listparindent
 \advance\@totalleftmargin\leftmargin
 \advance\rightskip\rightmargin
 \advance\leftskip\@totalleftmargin
}%
\let\set@listindent\set@listindent@parshape
\providecommand\href[0]{\begingroup\@sanitize@url\@href}%
\def\@href#1{\@@startlink{#1}\endgroup\@@href}%
\def\@@href#1{#1\@@endlink}%
\providecommand \url  [0]{\begingroup\@sanitize@url \@url }%
\def \@url #1{\endgroup\@href {#1}{\URL@prefix#1}}%
\providecommand \URL@prefix [0]{URL }%
\providecommand\doi[0]{\begingroup\@sanitize@url\@doi}%
\def\@doi#1{\endgroup\@@startlink{\doibase#1}doi:\discretionary {}{}{}#1\@@endlink }%
\providecommand \doibase [0]{http://dx.doi.org/}%
\providecommand \@sanitize@url[0]{\chardef\cat@space\the\catcode`\ \@sanitize\catcode`\ \cat@space}%
\def\@@startlink#1{}%
\def\@@endlink{}%
\@ifxundefined \pdfoutput {\true@sw}{\@ifnum{\z@=\pdfoutput}{\true@sw}{\false@sw}}%
{%
 \def\@@startlink@hypertext#1{\leavevmode\special{html:<a href="#1">}}%
 \def\@@endlink@hypertext{\special{html:</a>}}%
}{%
 \def\@@startlink@hypertext#1{%
  \leavevmode
  \pdfstartlink\pdfstartlink@attr
   user{/Subtype/Link/A<</Type/Action/S/URI/URI(#1)>>}%
  \relax
 }%
 \def\@@endlink@hypertext{\pdfendlink}%
 \def\pdfstartlink@attr{attr{/Border[0 0 1 ]/H/I/C[0 1 1]}}%
}%
\def\hypertext@enable@ltx{%
 \let\@@startlink\@@startlink@hypertext
 \let\@@endlink\@@endlink@hypertext
}%
\def\href@Hy{\hyper@normalise \href@ }%
\def\href@Hy@ltx{\@ifnextchar\bgroup\Hy@href{\hyper@normalise\href@}}%
\def\Hy@href#{\hyper@normalise\href@}%
\begingroup
  \endlinechar=-1 %
  \catcode`\^^A=14 %
  \catcode`\^^M\active
  \catcode`\%\active
  \catcode`\#\active
  \catcode`\_\active
  \catcode`\$\active
  \catcode`\&\active
  \gdef\hyper@normalise@ltx{^^A
    \begingroup
    \catcode`\^^M\active
    \def^^M{ }^^A
    \catcode`\%\active
    \let%\@percentchar
    \let\%\@percentchar
    \catcode`\#\active
    \def#{\hyper@hash}^^A
    \def\#{\hyper@hash}^^A
    \@makeother\&^^A
    \edef&{\string&}^^A
    \edef\&{\string&}^^A
    \edef\textunderscore{\string_}^^A
    \let\_\textunderscore
    \catcode`\_\active
    \let_\textunderscore
    \let~\hyper@tilde
    \let\~\hyper@tilde
    \let\textasciitilde\hyper@tilde
    \let\\\@backslashchar
    \edef${\string$}^^A
    \Hy@safe@activestrue
    \hyper@n@rmalise
  }^^A
  \catcode`\#=6 ^^A
  \gdef\Hy@ActiveCarriageReturn@ltx{^^M}^^A
  \gdef\hyper@n@rmalise@ltx#1#2{^^A
    \def\Hy@tempa{#2}^^A
    \ifx\Hy@tempa\Hy@ActiveCarriageReturn
      \Hy@ReturnAfterElseFi{^^A
        \hyper@@normalise{#1}^^A
      }^^A
    \else
      \Hy@ReturnAfterFi{^^A
        \hyper@@normalise{#1}{#2}^^A
      }^^A
    \fi
  }^^A
  \gdef\hyper@@normalise@ltx#1#2{^^A
    \edef\Hy@tempa{^^A
      \endgroup
      \noexpand#1{\Hy@RemovePercentCr#2%^^M\@nil}^^A
    }^^A
    \Hy@tempa
  }^^A
  \gdef\Hy@RemovePercentCr@ltx#1%^^M#2\@nil{^^A
    #1^^A
    \ifx\limits#2\limits
    \else
      \Hy@ReturnAfterFi{^^A
        \Hy@RemovePercentCr #2\@nil
      }^^A
    \fi
  }^^A
\endgroup
\def\switch@hyperref@href{%
 \expandafter\@ifx\expandafter{\csname href \endcsname\href@Hy}{
  \class@info{Repairing hyperref 6.75r \string\href}%
  \let\hyper@normalise\hyper@normalise@ltx
  \let\hyper@@normalise\hyper@@normalise@ltx
  \let\hyper@n@rmalise\hyper@n@rmalise@ltx
  \let\Hy@ActiveCarriageReturn\Hy@ActiveCarriageReturn@ltx
  \let\Hy@RemovePercentCr\Hy@RemovePercentCr@ltx
  \let\href\href@Hy@ltx
 }{}%
}%
\appdef\document@inithook{\switch@hyperref@href}%
\def\typeout@org#1{%
 \begingroup
  \set@display@protect
  \immediate\write\@unused{#1}%
 \endgroup
}%
\long\def\typeout@ltx#1{%
 \begingroup
  \set@display@protect
  \immediate\write\@unused{#1}%
 \endgroup
}%
\@ifx{\typeout\typeout@org}{%
 \class@info{Making \string\typeout\space \string\long}%
 \let\typeout\typeout@ltx
}{}%
\typeout{%
ltxfront%
 [2010/07/25/20:33:00 4.1r frontmatter package (AO,DPC)]% \fileversion
}%
\appdef\class@documenthook{\frontmatter@init}%
\let\frontmatter@init\@empty
\newcommand\frontmatter@title[2][]{%
 \def\@title{#2}%
 \def\@shorttitle{#1}%
 \let\@AF@join\@title@join
}%
\appdef\frontmatter@init{%
 \def\@title{\class@warn{No title}}%
 \let\@shorttitle\@empty
 \let\@title@aux\@title@aux@cleared
}%
\def\@title@join{\expandafter\@title@join@\@title@aux}%
\def\@title@join@#1#2{%
 \def\@title@aux{{\@join{\@separator}{#1}{#2}}}%
}%
\def\@title@aux@cleared{{}}%
\newcounter{affil}%
\newcounter{collab}%
\appdef\frontmatter@init{%
 \c@affil\z@
 \c@collab\z@
}%
\newcommand\frontmatter@author{% implicit #1
 \@author@def{}% implicit #2
}%
\def\collaboration{% implicit #1
 \@author@def{\@booleantrue\collaboration@sw}% implicit #2
}%
\appdef\frontmatter@init{%
 \@booleanfalse\collaboration@sw
}%
\def\@author@cleared{{}{}{}}%
\def\@author@gobble#1#2#3{}%
\def\@author@init{%
 \let\@author\@author@cleared
 \@booleanfalse\collaboration@sw
}%
\def\@authorclear@sw{\@ifx{\@author\@author@cleared}}%
\appdef\frontmatter@init{%
 \@author@init
}%
\def\@author@def#1#2{%
 \frontmatterverbose@sw{\typeout{\string\author\space\string\collaboration}}{}%
 \move@AU\move@AF\move@AUAF
 \let\@AF@join\@author@join
 #1%
 \def\@author{{#2}{}}%
}%
\def\@author@join@#1#2#3{%
  \def\@author{{#1}{\@join{\@separator}{#2}{#3}}}%
}%
\def\@author@join{\expandafter\@author@join@\@author}%
\def\move@AU{%
 \@authorclear@sw{}{%
  \collaboration@sw{%
   \advance\c@collab\@ne
   \@argswap{\CO@grp\CO@opr}%
  }{%
   \@argswap{\AU@grp\AU@opr}%
  }%
   {%
    \expandafter\@argswap@val
    \expandafter{\@author}%
     {\expandafter\@argswap@val\expandafter{\the\c@collab}{\add@AUCO@grp}}%
   }%
 }%
 \@author@init
}%
\def\add@AUCO@grp#1#2#3#4{%
 \appdef#3{#4{#1}#2}%
 \frontmatterverbose@sw{\say#3}{}%
}%
\def\@author@finish{%
 \frontmatterverbose@sw{\typeout{\string\@author@finish}}{}%
 \move@AU\move@AF
 \@ifx{\AU@grp\@empty}{%
  \@ifx{\CO@grp\@empty}%
 }{%
  \false@sw
 }%
 {}{%
  \@ifx{\AF@grp\@empty}{%
   \begingroup
    \let\href\@secondoftwo
    \let\AU@opr\@secondofthree
    \let\CO@opr\@secondofthree
    \let\footnote\@gobble
    \@ifx{\CO@grp\@empty}{%
     \class@warn{Assuming \string\noaffiliation\space for authors}%
     \frontmatterverbose@sw{\say\AU@grp}%
    }{%
     \class@warn{Assuming \string\noaffiliation\space for collaboration}%
     \frontmatterverbose@sw{\say\CO@grp}{}%
    }%
   \endgroup
   \@affil@none\move@AF
  }{}%
 }%
 \move@AUAF
}%
\def\@secondofthree#1#2#3{#2}%
\def\@join#1#2#3{%
  \@if@empty{#2}{#3}{#2#1#3}%
}%
\def\@separator{;\space}%
\let\surname\@firstofone
\let\firstname\@firstofone
\newcommand\frontmatter@and{\class@err{\protect\and\space is not supported}}
\def\cat@comma@active{\catcode`\,\active}%
{\cat@comma@active\gdef,{\active@comma}}%
\def\active@comma{,\penalty-300\relax}%
\newcommand\affiliation{%
 \frontmatterverbose@sw{\typeout{\string\affiliation}}{}%
 \move@AU\move@AF
 \begingroup
  \cat@comma@active
  \@affiliation
}%
\def\@affiliation#1{%
 \endgroup
 \let\@AF@join\@affil@join
 \@affil@def{#1}%
}%
\newcommand\frontmatter@noaffiliation{%
 \frontmatterverbose@sw{\typeout{\string\noaffiliation}}{}%
 \move@AU\move@AF
 \@affil@none\move@AF
 \move@AUAF
}%
\def\blankaffiliation{{}}%
\def\@affil@cleared{{{}}{}}%
\def\@affil@nil{{\relax}{}}%
\appdef\frontmatter@init{%
 \@affil@init
}%
\def\@affil@none{%
 \let\@affil\@affil@nil
}%
\def\@affil@init{%
 \let\@affil\@affil@cleared
}%
\def\@affilclear@sw{\@ifx{\@affil\@affil@cleared}}%
\def\@affil@def#1{%
 \def\@affil{{#1}{}}%
}%
\def\@affil@join@#1#2#3{%
  \def\@affil{{#1}{\@join{\@separator}{#2}{#3}}}%
}%
\def\@affil@join{\expandafter\@affil@join@\@affil}%
\def\move@AF{%
 \@affilclear@sw{}{%
  \@booleanfalse\temp@sw
  \let\@tempd\@empty
  \@affils@sw{%
    \expandafter\@affil@addr@def\expandafter\@tempa\@affil
    \def\AFF@opr{\@affil@match\@tempa}%
    \@AFF@list
  }{}\temp@sw
  {%
   \expandafter\@affil@aux@def\expandafter\@tempb\@affil
   \@ifx{\@tempb\@empty}{}{%
    \@ifx{\@tempb\@tempd}{}{%
     \class@warn{%
      Ancillary information for \@tempa\space must not be different!
      Please put all of it on the first instance%
     }%
    }%
   }%
  }%
  {%
   \@ifx{\@affil\@affil@nil}{%
    \def\@tempc{0}%
    \@argswap@val{0}%
   }{%
    \advance\c@affil\@ne
    \expandafter\def\expandafter\@tempc\expandafter{\the\c@affil}%
    \expandafter\@argswap@val\expandafter{\the\c@affil}%
   }%
   {%
    \expandafter\@argswap@val\expandafter{\the\c@collab}{%
     \expandafter\@argswap@val\expandafter{\@affil}{%
      \add@list@val@val@val\@AFF@list\AFF@opr
     }%
    }%
   }%
  }%
  \appdef@eval\AF@grp\@tempc
  \frontmatterverbose@sw{\say\AF@grp}{}%
  \@affil@init
 }%
}%
\def\@affil@addr@def#1#2#3{%
 \def#1{#2}%
}%
\def\@affil@aux@def#1#2#3{%
 \def#1{#3}%
}%
\def\add@list@val@val@val#1#2#3#4#5{%
 \appdef#1{#2{#5}{#4}#3}%
 \frontmatterverbose@sw{\say#1}{}%
}%
\def\@affil@match#1#2#3#4#5{%
 \temp@sw{}{%
  \def\@tempifx{#4}%
  \@ifx{\@tempifx#1}{%
   \groupauthors@sw{%
    \@ifnum{#3=\c@collab}{%
     \true@sw
    }{%
     \false@sw
    }%
   }{%
    \true@sw
   }%
  }{%
   \false@sw
  }%
  {%
   \@booleantrue\temp@sw
   \def\@tempc{#2}%
   \def\@tempd{#5}%
  }{%
  }%
 }%
}%
\def\move@AUAF{%
 \frontmatterverbose@sw{\say\AU@grp\say\AF@grp\say\CO@grp}{}%
 \@ifx{\AF@grp\@empty}{%
    \@ifx{\@empty\CO@grp}{%
    }{%
     \appdef     \@AAC@list{\AF@opr{{0}}}%
     \appdef@e   \@AAC@list{\CO@grp}%
     \appdef@e   \@AFG@list{\CO@grp}%
     \let\CO@grp\@empty
    }%
 }{%
      \appdef     \@AAC@list{\AF@opr}%
      \appdef@eval\@AAC@list{\AF@grp}%
      \appdef@e   \@AAC@list{\AU@grp}%
    \@ifx{\@empty\AU@grp}{%
     \@ifx{\@empty\CO@grp}%
    }{%
     \false@sw
    }%
    {%
    }{%
      \@booleanfalse\temp@sw
      \def\AFG@opr{\x@match\AF@grp}%
      \let\CO@opr\@author@gobble
      \@AFG@list
      \temp@sw{}{%
        \appdef     \@AFG@list{\AFG@opr}%
        \appdef@eval\@AFG@list{\AF@grp}%
      }%
     \@ifx{\@empty\CO@grp}{}{%
       \appdef@e   \@AAC@list{\CO@grp}%
       \appdef@e   \@AFG@list{\CO@grp}%
       \let\CO@grp\@empty
     }%
    }%
    \let\CO@grp\@empty
    \let\AU@grp\@empty
    \let\AF@grp\@empty
 }%
 \frontmatterverbose@sw{\say\@AAC@list\say\@AFG@list}{}%
}%
\appdef\frontmatter@init{%
 \let\AU@grp\@empty
 \let\CO@grp\@empty
 \let\AF@grp\@empty
 \let\@AAC@list\@empty
 \let\@AFG@list\@empty
 \let\@AFF@list\@empty
}%
\appdef\frontmatter@init{%
 \let\@AF@join\@AF@join@error
}%
\def\@AF@join@error#1{%
 \class@warn{%
  \string\email, \string\homepage, \string\thanks, or \string\altaffiliation\space
  appears in wrong context.
 }%
}%
\def\sanitize@url{%
 \@makeother\%%
 \@makeother\~%
 \@makeother\_%
}%
\newcommand*\email[1][]{\begingroup\sanitize@url\@email{#1}}%
\def\@email#1#2{%
 \endgroup
 \@AF@join{#1\href{mailto:#2}{#2}}%
}%
\newcommand*\homepage[1][]{\begingroup\sanitize@url\@homepage{#1}}%
\def\@homepage#1#2{%
 \endgroup
 \@AF@join{#1\href{#2}{#2}}%
}%
\appdef\class@documenthook{%
 \providecommand\href[1]{}%
}%
\def\frontmatter@thanks{% implicit #1
  \@AF@join
}%
\newcommand*\altaffiliation[2][]{%
  \@AF@join{#1#2}%
}%
\def\set@listcomma@list#1{%
  \expandafter\@reset@ac\expandafter#1#1{0}\@reset@ac{%
   \let\@listcomma\relax
  }{%
   \let\@listcomma\@listcomma@comma
  }%
}%
\def\set@listcomma@count#1{%
  \@ifnum{#1=\tw@}{%
    \let\@listcomma\relax
  }{%
    \let\@listcomma\@listcomma@comma
  }%
}%
\def\@reset@ac#1#2#3\@reset@ac{%
  \def#1{#3}%
  \@tempcnta#2\relax
  \@ifnum{#2=\tw@}%
}%
\def\@listand{\@ifnum{\@tempcnta=\tw@}{\andname\space}{}}%
\def\@listcomma@comma{\@ifnum{\@tempcnta>\@ne}{,}{}}%
\def\@listcomma@comma@UK{\@ifnum{\@tempcnta>\tw@}{,}{}}%
\def\@collaboration@gobble#1#2#3{}%
\def\doauthor#1#2#3{%
  \ignorespaces#1\unskip\@listcomma
  \begingroup
   #3%
  \@if@empty{#2}{\endgroup{}{}}{\endgroup{\comma@space}{}\frontmatter@footnote{#2}}%
  \space \@listand
}%
\def\x@match#1#2{%
 \temp@sw{}{%
  \def\@tempifx{#2}%
  \@ifx{\@tempifx#1}{%
    \@booleantrue\temp@sw
  }{%
  }%
 }%
}%
\def\y@match#1#2#3{%
 \temp@sw{}{%
  \def\@tempifx{#3}%
  \@ifx{\@tempifx#1}{%
    \@booleantrue\temp@sw
    \def\@tempb{#2}%
  }{%
  }%
 }%
}%
\def\frontmatter@footnote#1{%
 \begingroup
  \@booleanfalse\temp@sw
  \def\@tempa{#1}%
  \let\@tempb\@empty
  \def\@TBN@opr{\y@match\@tempa}%
  \@FMN@list
  \temp@sw{%
   \expandafter\frontmatter@footnotemark
   \expandafter{\@tempb}%
  }{%
   \stepcounter\@mpfn
   \expandafter\expandafter
   \expandafter\frontmatter@foot@mark
   \expandafter\expandafter
   \expandafter{%
   \expandafter \the\csname c@\@mpfn\endcsname
               }{#1}%
  }%
 \endgroup
}%
\def\frontmatter@foot@mark#1#2{%
 \frontmatter@footnotemark{#1}%
 \g@addto@macro\@FMN@list{\@TBN@opr{#1}{#2}}%
}%
\appdef\frontmatter@init{%
 \global\let\@FMN@list\@empty
}%
\def\frontmatter@footnotemark#1{%
 \leavevmode
 \ifhmode\edef\@x@sf{\the\spacefactor}\nobreak\fi
  \begingroup
   \hyper@linkstart {link}{frontmatter.#1}%
    \csname c@\@mpfn\endcsname#1\relax
    \def\@thefnmark{\frontmatter@thefootnote}%
    \@makefnmark
   \hyper@linkend
  \endgroup
 \ifhmode\spacefactor\@x@sf\fi
 \relax
}%
\def\keywords#1{%
  \aftermaketitle@chk{\keywords}%
  \gdef\@keywords{#1}%
}%
\appdef\frontmatter@init{%
 \let\@keywords\@empty
}%
\def\features#1{%
  \aftermaketitle@chk{\features}%
  \gdef\@features{#1}%
}%
\appdef\frontmatter@init{%
 \let\@features\@empty
}%
\newcommand*\frontmatter@date[2][\Dated@name]{\def\@date{#1#2}}%
\def\@date{}%
\newcommand*\received[2][\Received@name]{\def\@received{#1#2}}%
\def\@received{}%
\newcommand*\revised[2][\Revised@name]{\def\@revised{#1#2}}%
\def\@revised{}%
\newcommand*\accepted[2][\Accepted@name]{\def\@accepted{#1#2}}%
\def\@accepted{}%
\newcommand*\published[2][\Published@name]{\def\@published{#1#2}}%
\def\@published{}%
\def\pacs#1{%
  \aftermaketitle@chk{\pacs}%
  \gdef\@pacs{#1}%
}%
\appdef\frontmatter@init{%
 \let\@pacs\@empty
}%
\def\preprint#1{\gappdef\@preprint{\preprint{#1}}}%
\appdef\frontmatter@init{%
 \let\@preprint\@empty
}%
\newbox\absbox
\def\toclevel@abstract{1}%
\def\addcontents@abstract{%
 \phantomsection
 \expandafter\def\csname Parent0\endcsname{section*.2}%
 \expandafter\@argswap@val\expandafter{\abstractname}{\addcontentsline{toc}{abstract}}%
}%
\newenvironment{frontmatter@abstract}{%
  \aftermaketitle@chk{\begin{abstract}}%
  \global\setbox\absbox\vbox\bgroup
   \color@begingroup
   \columnwidth\textwidth
   \hsize\columnwidth
   \@parboxrestore
   \def\@mpfn{mpfootnote}\def\thempfn{\thempfootnote}\c@mpfootnote\z@
   \let\@footnotetext\frontmatter@footnotetext
   \minipagefootnote@init
   \let\set@listindent\set@listindent@
   \let\@listdepth\@mplistdepth \@mplistdepth\z@
   \let@environment{description}{frontmatter@description}%
   \@minipagerestore
   \@setminipage
    \frontmatter@abstractheading
    \frontmatter@abstractfont
    \let\footnote\mini@note
    \expandafter\everypar\expandafter{\the\everypar\addcontents@abstract\everypar{}}%
}{%
    \par
    \unskip
    \minipagefootnote@here
    \@minipagefalse   %% added 24 May 89
    \color@endgroup
  \egroup
}%
\long\def\frontmatter@footnotetext#1{%
  \minipagefootnote@pick
    \set@footnotefont
    \set@footnotewidth
    \@parboxrestore
    \protected@edef\@currentlabel{\csname p@\@mpfn\endcsname\@thefnmark}%
    \color@begingroup
      \frontmatter@makefntext{%
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox\vadjust{\vskip\z@skip}%
      }%
    \color@endgroup
  \minipagefootnote@drop
}%
\def\ltx@no@footnote{%
 \let\ltx@xfootnote\ltx@no@xfootnote\let\ltx@yfootnote\ltx@no@yfootnote
 \let\ltx@xfootmark\ltx@no@xfootmark\let\ltx@yfootmark\ltx@no@yfootmark
 \let\ltx@xfoottext\ltx@no@xfoottext\let\ltx@yfoottext\ltx@no@yfoottext
}%
\def\ltx@no@xfootnote[#1]#2{\ltx@no@footwarn\footnote}%
\def\ltx@no@yfootnote#1{\ltx@no@footwarn\footnote}%
\def\ltx@no@xfootmark[#1]{\ltx@no@footwarn\footnotemark}%
\def\ltx@no@yfootmark{\ltx@no@footwarn\footnotemark}%
\def\ltx@no@xfoottext[#1]#2{\ltx@no@footwarn\footnotetext}%
\def\ltx@no@yfoottext#1{\ltx@no@footwarn\footnotetext}%
\def\ltx@no@footwarn#1{%
 \class@warn{%
  The \string#1\space command is not legal on the title page;
  using \string\thanks\space instead might suit you: consult the manual for details%
 }%
}%
\def\frontmatter@abstractheading{%
 \begingroup
  \centering\large
  \abstractname
  \par
 \endgroup
}%
\def\frontmatter@abstractfont{}%
\newenvironment{frontmatter@description}{%
 \list{}{%
  \leftmargin\z@
  \labelwidth\z@
  \itemindent\z@
  \let\makelabel\frontmatter@descriptionlabel
 }%
}{%
 \endlist
}%
\def\frontmatter@descriptionlabel#1{%
 \hspace\labelsep
 \normalfont\bfseries
 #1:%
}%
\def\frontmatter@abstractwidth{\textwidth}
\def\frontmatter@abstract@produce{%
  \par
  \preprintsty@sw{%
   \do@output@MVL{%
    \vskip\frontmatter@preabstractspace
    \vskip200\p@\@plus1fil
    \penalty-200\relax
    \vskip-200\p@\@plus-1fil
   }%
  }{%
   \addvspace{\frontmatter@preabstractspace}%
  }%
   \begingroup
    \dimen@\baselineskip
    \setbox\z@\vtop{\unvcopy\absbox}%
    \advance\dimen@-\ht\z@\advance\dimen@-\prevdepth
    \@ifdim{\dimen@>\z@}{\vskip\dimen@}{}%
   \endgroup
   \begingroup
    \prep@absbox
    \unvbox\absbox
    \post@absbox
   \endgroup
  \@ifx{\@empty\mini@notes}{}{\mini@notes\par}%
  \addvspace\frontmatter@postabstractspace
}%
\appdef\frontmatter@init{\let\mini@notes\@empty}%
\let\prep@absbox\@empty
\let\post@absbox\@empty
\def\frontmatter@preabstractspace{.5\baselineskip}
\def\frontmatter@postabstractspace{.5\baselineskip}
\newenvironment{frontmatter@titlepage}{%
      \twocolumn@sw{\onecolumngrid}{\newpage}%
      \thispagestyle{titlepage}%
      \setcounter{page}\@ne
}{%
     \twocolumn@sw{\twocolumngrid}{\newpage}%
     \twoside@sw{}{%
        \setcounter{page}\@ne
     }%
}%
\def\frontmatter@maketitle{%
  \@author@finish
  \title@column\titleblock@produce
  \suppressfloats[t]%
  \let\and\relax
  \let\affiliation\@gobble
  \let\author\@gobble
  \let\@AAC@list\@empty
  \let\@AFF@list\@empty
  \let\@AFG@list\@empty
  \let\@AF@join\@AF@join@error
  \let\email\@gobble
  \let\@address\@empty
  \let\maketitle\relax
  \let\thanks\@gobble
  \let\abstract\@undefined\let\endabstract\@undefined
  \titlepage@sw{%
   \vfil
   \clearpage
  }{}%
}%
\def\maketitle@Hy{%
  \let\Hy@saved@footnotemark\@footnotemark
  \let\Hy@saved@footnotetext\@footnotetext
  \let\@footnotemark\H@@footnotemark
  \let\@footnotetext\H@@footnotetext
  \@ifnextchar[%]
   \Hy@maketitle@optarg
   {%
    \HyOrg@maketitle
    \Hy@maketitle@end
   }%
}%
\appdef\class@documenthook{%
  \@ifx{\maketitle\maketitle@Hy}{%
   \class@info{Taking \string\maketitle\space back from hyperref}%
   \let\maketitle\frontmatter@maketitle
  }{%
  }%
}%
\def\titleblock@produce{%
 \begingroup
  \ltx@footnote@pop
  \def\@mpfn{mpfootnote}%
  \def\thempfn{\thempfootnote}%
  \c@mpfootnote\z@
  \let\@makefnmark\frontmatter@makefnmark
  \frontmatter@setup
  \thispagestyle{titlepage}\label{FirstPage}%
  \frontmatter@title@produce
  \groupauthors@sw{%
   \frontmatter@author@produce@group
  }{%
   \frontmatter@author@produce@script
  }%
  \frontmatter@RRAPformat{%
   \expandafter\produce@RRAP\expandafter{\@date}%
   \expandafter\produce@RRAP\expandafter{\@received}%
   \expandafter\produce@RRAP\expandafter{\@revised}%
   \expandafter\produce@RRAP\expandafter{\@accepted}%
   \expandafter\produce@RRAP\expandafter{\@published}%
  }%
  \frontmatter@abstract@produce
  \@ifx@empty\@pacs{}{%
   \@pacs@produce\@pacs
  }%
  \@ifx@empty\@keywords{}{%
   \@keywords@produce\@keywords
  }%
  \@ifx@empty\@features{}{%
   \@features@produce\@features
  }%
  \par
  \frontmatter@finalspace
 \endgroup
}%
\def\toclevel@title{0}%
\def\frontmatter@title@produce{%
 \begingroup
  \frontmatter@title@above
  \frontmatter@title@format
  \@title
  \unskip
  \phantomsection\expandafter\@argswap@val\expandafter{\@title}{\addcontentsline{toc}{title}}%
  \@ifx{\@title@aux\@title@aux@cleared}{}{%
   \expandafter\frontmatter@footnote\expandafter{\@title@aux}%
  }%
  \par
  \frontmatter@title@below
 \endgroup
}%
\appdef\let@mark{\let\\\relax}%
\def\frontmatter@title@above{}%
\def\frontmatter@title@format{}%
\def\frontmatter@title@below{\addvspace{\baselineskip}}%
\def\frontmatter@author@produce@script{%
  \begingroup
    \let\@author@present\@author@present@script
    \frontmatterverbose@sw{\typeout{\string\frontmatter@author@produce@script:}\say\@AAC@list\say\@AFF@list\say\@AFG@list}{}%
    \let\AU@temp\@empty
    \@tempcnta\z@
    \let\AF@opr \@gobble
    \def\AU@opr{\@author@count\@tempcnta}%
    \def\CO@opr{\@collaboration@count\AU@temp\@tempcnta}%
    \@AAC@list
    \expandafter\CO@opr\@author@cleared
    \begingroup
     \frontmatter@authorformat
     \let\AF@opr \@affilID@def
     \let\AU@opr \@author@present
     \def\CO@opr{\@collaboration@present\AU@temp}%
     \set@listcomma@list\AU@temp
     \@AAC@list
     \unskip\unskip
     \par
    \endgroup
    \begingroup
     \frontmatter@above@affiliation@script
     \let\AFF@opr \@affil@script
     \@AFF@list
     \frontmatter@footnote@produce
     \par
    \endgroup
  \endgroup
}%
\def\@author@count#1{%
 \advance#1\@ne
 \@author@gobble
}%
\def\@collaboration@present#1#2#3#4{%
 \par
 \begingroup
  \frontmatter@collaboration@above
  \@affilID@def{}%
  \@tempcnta\z@
  \@author@present{}{(\ignorespaces#3\unskip)}{#4}%
  \par
 \endgroup
 \set@listcomma@list#1%
}%
\def\frontmatter@collaboration@above{}%
\def\@collaboration@count#1#2{%
 \appdef@eval#1{\the#2}#2\z@
 \@author@gobble
}%
\def\@affilID@def{\def\@affilID@temp}%
\let\@affilID@temp\@empty
\def\affil@script#1#2#3{%
 \def\@tempifx{#1}\@ifx{\@tempifx\@tempa}{%
  \@if@empty{#2}{}{%
   \par
   \begingroup
    \def\@thefnmark{#1}\@makefnmark\ignorespaces
    #2%
    \@if@empty{#3}{}{\frontmatter@footnote{#3}}%
    \par
   \endgroup
  }%
 }{}%
}%
\def\@affil@script#1#2#3#4{%
 \@ifnum{#1=\z@}{}{%
  \par
  \begingroup
   \frontmatter@affiliationfont
   \@ifnum{\c@affil<\affil@cutoff}{}{%
    \def\@thefnmark{#1}\@makefnmark
   }%
   \ignorespaces#3%
   \@if@empty{#4}{}{\frontmatter@footnote{#4}}%
   \par
  \endgroup
 }%
}%
\let\affil@cutoff\@ne
\def\@author@present@script#1#2#3{%
 \begingroup
  \gdef\comma@space{\textsuperscript{,\,}}%
  \doauthor{#2}{#3}{\@affil@present@script}%
 \endgroup
 \advance\@tempcnta\m@ne
}%
\def\@affilcomma#1#2{%
 \@ifx{\z@#1}{%
  \@ifx{\relax#2}{}{%
   \@affilcomma{#2}%
  }%
 }{%
  #1%
  \@ifx{\relax#2}{}{%
   \@ifx{\z@#2}{%
    \@affilcomma
   }{%
    ,\,\@affilcomma{#2}%
   }%
  }%
 }%
}%
\def\@affil@present@script{%
 \let\@tempa\@empty
 \expandafter\@affil@present@script@\@affilID@temp\relax
}%
\def\@affil@present@script@#1{%
 \@ifx{\relax#1}{%
  \@ifx{\@tempa\@empty}{%
   \aftergroup\false@sw
  }{%
   \textsuperscript{\expandafter\@affilcomma\@tempa\relax\relax}%
   \aftergroup\true@sw
  }%
 }{%
  \@ifnum{#1=\z@}{}{\appdef\@tempa{{#1}}}%
  \@affil@present@script@
 }%
}%
\@provide\@author@parskip{\z@skip}%
\def\frontmatter@author@produce@group{%
  \begingroup
    \let\@author@present\@author@present@group
    \frontmatter@authorformat
    \frontmatterverbose@sw{\typeout{\string\frontmatter@author@produce@group:}\say\@AAC@list\say\@AFF@list\say\@AFG@list}{}%
    \let\AU@temp\@empty
    \set@listcomma@list\AU@temp
    \def\CO@opr{\@collaboration@present\AU@temp}%
    \let\AFG@opr \affils@present@group
    \let\@listcomma\relax
    \@AFG@list
    \frontmatter@footnote@produce
    \par
  \endgroup
  \frontmatter@authorbelow
}%
\@provide\frontmatter@authorbelow{}%
\def\affils@present@group#1{%
 \begingroup
   \def\AF@temp{#1}%
   \@tempcnta\z@
   \let\AU@opr \@undefined
   \let\CO@opr \@undefined
   \def\AF@opr{\@affilID@count\AF@temp\@tempcnta}%
   \@AAC@list
   \@ifnum{\@tempcnta=\z@}{}{%
    \begingroup
     \frontmatter@above@affilgroup
     \set@listcomma@count\@tempcnta
     \let\AU@opr \@undefined
     \let\CO@opr \@undefined
     \def\AF@opr{\@affilID@match\AF@temp}%
     \@AAC@list
    \endgroup
    \begingroup
     \par
     \frontmatter@above@affiliation
     \frontmatter@affiliationfont
     \let\\\frontmatter@addressnewline
     \@tempcnta\z@
     \@tfor\AF@temp:=#1\do{%
      \expandafter\@ifx\expandafter{\expandafter\z@\AF@temp}{}{%
       \advance\@tempcnta\@ne
      }%
     }%
     \@ifnum{\@tempcnta=\tw@}{%
      \let\@listcomma\relax
     }{}%
     \def@after@address
     \runinaddress@sw{%
     }{%
      \tightenlines@sw{}{%
       \parskip\z@
      }%
      \appdef\after@address\par
     }%
     \let\AFF@opr \@affil@group
     \do@affil@fromgroup\@AFF@list#1\relax
    \endgroup
   }%
   \par
 \endgroup
}%
\def\def@after@address{\def\after@address{\@listcomma\ \@listand}}%
\def\def@after@address@empty{\let\after@address\@empty}%
\def\@affilID@count#1#2#3{%
  \def\@tempifx{#3}%
  \@ifx{\@tempifx#1}{%
    \def\AU@opr{\@author@count#2}%
  }{%
    \let\AU@opr \@author@gobble
  }%
  \let\CO@opr \@collaboration@gobble
}%
\def\@affilID@match#1#2{%
 \def\@tempifx{#2}%
 \@ifx{\@tempifx#1}{%
   \let\AU@opr \@author@present
 }{%
   \let\AU@opr \@author@gobble
 }%
  \let\CO@opr \@collaboration@gobble
}%
\def\do@affil@fromgroup#1#2{%
  \@ifx{\relax#2}{}{%
    \count@#2\relax
    \@ifnum{\z@=\count@}{}{#1}%
    \do@affil@fromgroup#1%
  }%
}%
\def\@affil@group#1#2#3#4{%
  \@ifnum{#1=\count@}{%
   \def\@tempa{#3}%
   \@ifx{\@tempa\blankaffiliation}{}{%
    #3%
    \@if@empty{#4}{}{%
     \frontmatter@footnote{#4}%
    }%
    \after@address
   }%
   \advance\@tempcnta\m@ne
  }{}%
}%
\def\@author@present@group#1#2#3{%
  \gdef\comma@space{\gdef\comma@space{\textsuperscript{,\,}}}%
  \doauthor{#2}{#3}{\@affil@present@group}%
  \advance\@tempcnta\m@ne
}%
\def\@affil@present@group{%
 \aftergroup\false@sw
}%
\def\@pacs@produce#1{%
 \showPACS@sw{%
  \begingroup
   \frontmatter@PACS@format
   \@pacs@name#1\par
  \endgroup
 }{%
  \@if@empty{#1}{}{%
   \class@warn{\PACS@warn}%
  }%
 }%
}%
\def\PACS@warn{If you want your PACS to appear in your output, use document class option showpacs}%
\def\@keywords@produce#1{%
 \showKEYS@sw{%
  \begingroup
   \frontmatter@keys@format
   \@keys@name#1\par
  \endgroup
 }{%
  \@if@empty{#1}{}{%
   \class@warn{If you want your keywords to appear in your output, use document class option showkeys}%
  }%
 }%
}%
\def\@features@produce#1{%
 \showFEATS@sw{%
  \begingroup
   \frontmatter@feats@format
   \@feats@name#1\par
  \endgroup
 }{%
  \@if@empty{#1}{}{%
   \class@warn{If you want your features to appear in your output, use document class option showfeats}%
  }%
 }%
}%
\def\frontmatter@footnote@produce@footnote{%
 \let\@TBN@opr\present@FM@footnote
 \@FMN@list
 \global\let\@FMN@list\@empty
}%
\def\present@FM@footnote#1#2{%
 \begingroup
  \csname c@\@mpfn\endcsname#1\relax
  \def\@thefnmark{\frontmatter@thefootnote}%
  \frontmatter@footnotetext{#2}%
 \endgroup
}%
\def\frontmatter@footnote@produce@endnote{%
}%
\appdef\frontmatter@init{%
 \@ifxundefined\title@column         {\let\title@column\@empty}{}%
 \@ifxundefined\preprintsty@sw       {\@booleanfalse\preprintsty@sw}{}%
 \@ifxundefined\frontmatter@footnote@produce{\let\frontmatter@footnote@produce\frontmatter@footnote@produce@footnote}{}%
 \@ifxundefined\do@output@MVL        {\let\do@output@MVL\@firstofone}{}%
 \@ifxundefined\comma@space          {\let\comma@space\@empty}{}%
}%
\def\frontmatter@thefootnote{%
 \altaffilletter@sw{\@alph}{\@fnsymbol}{\csname c@\@mpfn\endcsname}%
}%
\@ifx{\altaffilletter@sw\@undefined}{\@booleantrue\altaffilletter@sw}{}%
\def\frontmatter@makefnmark{%
 \@textsuperscript{%
  \normalfont\@thefnmark
 }%
}%
\long\def\frontmatter@makefntext#1{%
 \parindent 1em
 \noindent
 \Hy@raisedlink{\hyper@anchorstart{frontmatter.\expandafter\the\csname c@\@mpfn\endcsname}\hyper@anchorend}%
 \@makefnmark
 #1%
}%
\def\frontmatter@setup{}%
\def\frontmatter@RRAPformat#1{%
 \removelastskip
 \begingroup
  \frontmatter@RRAP@format
  #1\par
 \endgroup
}%
\def\punct@RRAP{; }%
\def\produce@RRAP#1{%
  \@if@empty{#1}{}{%
   \@ifvmode{\leavevmode}{\unskip\punct@RRAP\ignorespaces}%
   #1%
  }%
}%
\def\frontmatter@authorformat{}%
\def\frontmatter@above@affilgroup{}%
\def\frontmatter@above@affiliation{}%
\def\frontmatter@above@affiliation@script{}%
\def\frontmatter@affiliationfont{\itshape\selectfont}%
\def\frontmatter@RRAP@format{}%
\def\frontmatter@PACS@format{}%
\def\frontmatter@keys@format{}%
\def\frontmatter@feats@format{}%
\def\frontmatter@finalspace{\addvspace{18\p@}}
\def\frontmatter@addressnewline{%
  \@ifhmode{\skip@\lastskip\unskip\unpenalty\break\hskip\skip@}{}%
  % was: \vskip-.5ex
}%
\def\frontmatter@preabstractspace{5.5\p@}
\def\frontmatter@postabstractspace{6.5\p@}
\def\aftermaketitle@chk#1{%
  \@ifx{\maketitle\relax}{%
    \class@err{\protect#1 must be used before \protect\maketitle}%
  }{}%
}%
\def\ps@titlepage{\ps@empty}%
\def\volumeyear#1{\gdef\@volumeyear{#1}}%
\def\@volumeyear{}%
\def\volumenumber#1{\gdef\@volumenumber{#1}}%
\def\@volumenumber{}%
\def\issuenumber#1{\gdef\@issuenumber{#1}}%
\def\@issuenumber{}%
\def\eid#1{\gdef\@eid{#1}}%
\def\@eid{}%
\def\startpage#1{\gdef\@startpage{#1}\c@page#1\relax}%
\def\@startpage{\pageref{FirstPage}}%
\def\endpage#1{\gdef\@endpage{#1}}%
\def\@endpage{\pageref{LastPage}}%
\def\print@toc#1{%
 \begingroup
  \expandafter\section
  \expandafter*%
  \expandafter{%
              \csname#1name\endcsname
              }%
  \let\appendix\appendix@toc
  \@starttoc{#1}%
 \endgroup
}%
\def\appendix@toc{}%
\def\Dated@name{Dated }%
\def\Received@name{Received }%
\def\Revised@name{Revised }%
\def\Accepted@name{Accepted }%
\def\Published@name{Published }%
\appdef\robustify@contents{%
 \let\thanks\@gobble\let\class@warn\@gobble
 \def\begin{\string\begin}\def\end{\string\end}%
}%
\@ifxundefined\frontmatter@syntax@sw{\@booleantrue\frontmatter@syntax@sw}{}%
\frontmatter@syntax@sw{%
 \let\title         \frontmatter@title
 \let\author        \frontmatter@author
 \let\date          \frontmatter@date
 \@ifxundefined\@maketitle{%
  \let\maketitle    \frontmatter@maketitle
  \@booleantrue     \titlepage@sw
 }{%
  \let\@maketitle   \frontmatter@maketitle
  \prepdef\maketitle\@author@finish
 }%
 \let\noaffiliation \frontmatter@noaffiliation
 \let\thanks@latex  \thanks
 \let\thanks        \frontmatter@thanks
 \let\and@latex     \and
 \let\and           \frontmatter@and
 \let@environment{titlepage}{frontmatter@titlepage}%
 \let@environment{abstract}{frontmatter@abstract}%
}{%
 \let\noaffiliation\@empty
}%
\typeout{%
ltxgrid%
 [2010/07/25/20:33:00 4.1r page grid package (portions licensed from W. E. Baxter web at superscript.com)]% \fileversion
}%
\newcounter{linecount}
\def\loop@line#1#2{%
 \par
 \hb@xt@\hsize{%
  \global\advance#1\@ne
  \edef\@tempa{\@ifnum{100>#1}{0}{}\@ifnum{10>#1}{0}{}\number#1}%
  \@tempa\edef\@tempa{\special{line:\@tempa}}\@tempa
  \vrule depth2.5\p@#2\leaders\hrule\hfil
 }%
}%
\def\lineloop#1{%
 \loopwhile{\loop@line\c@linecount{}\@ifnum{#1>\c@linecount}}%
}%
\def\linefoot#1{%
 \loop@line\c@linecount{%
  \footnote{%
   #1\special{foot:#1}\vrule depth2.5\p@\leaders\hrule\hfill
  }%
 }%
}%
\let\@@mark\mark
\let\@@topmark\topmark
\let\@@firstmark\firstmark
\let\@@botmark\botmark
\let\@@splitfirstmark\splitfirstmark
\let\@@splitbotmark\splitbotmark
\def\@themark{{}{}{}{}}%
\def\nul@mark{{}{}{}{}\@@nul}%
\def\set@mark@netw@#1#2#3#4#5#6#7{\gdef#1{{#6}{#7}{#4}{#5}}\do@mark}%
\def\set@marktw@#1#2#3#4#5#6{\gdef#1{{#2}{#6}{#4}{#5}}\do@mark}%
\def\set@markthr@@#1#2#3#4#5#6{\gdef#1{{#2}{#3}{#6}{#5}}\do@mark}%
\def\get@mark@@ne#1#2#3#4#5\@@nul{#1}%
\def\get@mark@tw@#1#2#3#4#5\@@nul{#2}%
\def\get@mark@thr@@#1#2#3#4#5\@@nul{#3}%
\def\get@mark@f@ur#1#2#3#4#5\@@nul{#4}%
\def\mark@netw@{\expandafter\set@mark@netw@\expandafter\@themark\@themark}%
\def\marktw@{\expandafter\set@marktw@\expandafter\@themark\@themark}%
\def\markthr@@{\expandafter\set@markthr@@\expandafter\@themark\@themark}%
\def\do@mark{\do@@mark\@themark\nobreak@mark}%
\def\do@@mark#1{%
 \begingroup
  \let@mark
  \@@mark{#1}%
 \endgroup
}%
\def\let@mark{%
 \let\protect\@unexpandable@protect
 \let\label\relax
 \let\index\relax
 \let\glossary\relax
}%
\def\nobreak@mark{%
 \@if@sw\if@nobreak\fi{\@ifvmode{\nobreak}{}}{}%
}%
\def\mark@envir{\markthr@@}%
\def\bot@envir{%
 \expandafter\expandafter
 \expandafter\get@mark@thr@@
 \expandafter\@@botmark
             \nul@mark
}%
\def\markboth{\mark@netw@}%
\def\markright{\marktw@}%
\def\leftmark{%
 \expandafter\expandafter
 \expandafter\get@mark@@ne
 \expandafter\saved@@botmark
             \nul@mark
}%
\def\rightmark{%
 \expandafter\expandafter
 \expandafter\get@mark@tw@
 \expandafter\saved@@firstmark
             \nul@mark
}%
\let\primitive@output\output
\long\def\@tempa#1\@@nil{#1}%
            \toks@
\expandafter\expandafter
\expandafter{%
\expandafter \@tempa
             \the\primitive@output
             \@@nil
             }%
\newtoks\output@latex
\output@latex\expandafter{\the\toks@}%
\let\output\output@latex
\primitive@output{\dispatch@output}%
\def\dispatch@output{%
 \let\par\@@par
 \expandafter\let\expandafter\output@procedure\csname output@\the\outputpenalty\endcsname
 \@ifnotrelax\output@procedure{}{%
  \expandafter\def\expandafter\output@procedure\expandafter{\the\output@latex}%
 }%
 \expandafter\@ifx\expandafter{\csname output@-\the\execute@message@pen\endcsname\output@procedure}{%
  \let\output@procedure\@message@saved
 }{}%
 \ltxgrid@info@sw{\class@info{\string\dispatch@output}\say\output@procedure\saythe\holdinginserts}{}%
 \outputdebug@sw{\output@debug}{}%
 \output@procedure
}%
\def\set@output@procedure#1#2{%
 \count@\outputpenalty\advance\count@-#2%
 \expandafter\let\expandafter#1\csname output@\the\count@\endcsname
}%
\def\output@debug{%
 \def\@tempa{\save@message}%
 \@ifx{\output@procedure\@tempa}{%
  \true@sw
 }{%
  \@ifnum{\outputpenalty=-\save@column@insert@pen}{%
   \@ifnum{\holdinginserts>\z@}%
  }{%
   \false@sw
  }%
 }%
 {}{\output@debug@}%
}%
\def\output@debug@{%
  \saythe\outputpenalty
  \saythe\interlinepenalty
  \saythe\brokenpenalty
  \saythe\clubpenalty
  \saythe\widowpenalty
  \saythe\displaywidowpenalty
  \saythe\predisplaypenalty
  \saythe\interdisplaylinepenalty
  \saythe\postdisplaypenalty
  \saythe\badness
  \say\thepagegrid
  \saythe\pagegrid@col
  \saythe\pagegrid@cur
  \saythe\insertpenalties
  \say\@@botmark
  \saythe\pagegoal
  \saythe\pagetotal
  \saythe{\badness\@cclv}%
  \say\@toplist
  \say\@botlist
  \say\@dbltoplist
  \say\@deferlist
  \trace@scroll{%
  \showbox\@cclv
  \showbox\@cclv@saved
  \showbox\pagesofar
  \showbox\csname col@1\endcsname
  \showbox\footsofar
  \showbox\footins
  \showbox\footins@saved
  \showlists
  }%
}%
\@ifxundefined{\outputdebug@sw}{%
 \@booleanfalse\outputdebug@sw
}{}%
\def\trace@scroll#1{\begingroup\showboxbreadth\maxdimen\showboxdepth\maxdimen\scrollmode#1\endgroup}%
\def\trace@box#1{\trace@scroll{\showbox#1}}%
\prepdef\@outputpage{\@outputpage@head}%
\let\@outputpage@head\@empty
\appdef\@outputpage{\@outputpage@tail}%
\let\@outputpage@tail\@empty
\def\show@box@size#1#2{%
 \show@box@size@sw{%
  \begingroup
   \setbox\z@\vbox{\unvcopy#2\hrule}%
   \class@info{Show box size: #1^^J%
    (\the\ht\z@\space X \the\wd\z@)
    \the\c@page\space\space\the\pagegrid@cur\space\the\pagegrid@col
   }%
  \endgroup
 }{}%
}%
\def\show@text@box@size{%
 \show@box@size{Text column}\@outputbox
 \tally@box@size@sw{%
  \@ifdim{\wd\@outputbox>\z@}{%
   \dimen@\ht\@outputbox\divide\dimen@\@twopowerfourteen
   \advance\dimen@-\dp\csname box@size@\the\pagegrid@col\endcsname
   \@ifdim{\dimen@>\z@}{%
    \advance\dimen@ \ht\csname box@size@\the\pagegrid@col\endcsname
    \global\ht\csname box@size@\the\pagegrid@col\endcsname\dimen@
    \show@box@size@sw{%
     \class@info{Column: \the\dimen@}%
    }{}%
   }{}%
  }{}%
  \global\dp\csname box@size@\the\pagegrid@col\endcsname\z@
 }{}%
}%
\def\show@pagesofar@size{%
 \show@box@size{Page so far}\pagesofar
 \dimen@\ht\pagesofar\divide\dimen@\@twopowerfourteen
 \global\dp\csname box@size@1\endcsname\dimen@
 \show@box@size@sw{%
  \class@info{Pagesofar: \the\dimen@}%
 }{}%
}%
\@booleanfalse\tally@box@size@sw
\@booleanfalse\show@box@size@sw
\expandafter\newbox\csname box@size@1\endcsname
\expandafter\setbox\csname box@size@1\endcsname\hbox{}%
\expandafter\newbox\csname box@size@2\endcsname
\expandafter\setbox\csname box@size@2\endcsname\hbox{}%
\def\total@text{%
 \@tempdima\the\ht\csname box@size@2\endcsname\divide\@tempdima\@twopowertwo\@tempcnta\@tempdima
 \@tempdimb\the\ht\csname box@size@1\endcsname\divide\@tempdimb\@twopowertwo\@tempcntb\@tempdimb
 \class@info{Total text: Column(\the\@tempcnta pt), Page(\the\@tempcntb pt)}%
}%
\def\natural@output{\toggle@insert{\output@holding}{\output@moving}}%
\output@latex{\natural@output}%
\def\output@holding{%
 \csname output@init@\bot@envir\endcsname
 \@if@exceed@pagegoal{\unvcopy\@cclv}{%
  \setbox\z@\vbox{\unvcopy\@cclv}%
  \outputdebug@sw{\trace@box\z@}{}%
  \dimen@\ht\@cclv\advance\dimen@-\ht\z@
  \dead@cycle@repair\dimen@
 }{%
  \dead@cycle
 }%
}%
\def\@if@exceed@pagegoal#1{%
 \begingroup
  \setbox\z@\vbox{#1}%
  \dimen@\ht\z@\advance\dimen@\dp\z@
  \outputdebug@sw{\saythe\dimen@}{}%
  \@ifdim{\dimen@>\pagegoal}{%
   \setbox\z@\vbox{\@@mark{}\unvbox\z@}%
   \splittopskip\topskip
   \splitmaxdepth\maxdepth
   \vbadness\@M
   \vfuzz\maxdimen
   \setbox\tw@\vsplit\z@ to\pagegoal
   \outputdebug@sw{\trace@scroll{\showbox\tw@\showbox\z@}}{}%
   \setbox\tw@\vbox{\unvbox\tw@}%
   \@ifdim{\ht\tw@=\z@}{%
    \ltxgrid@info{Found overly large chunk while preparing to move insertions. Attempting repairs}%
    \aftergroup\true@sw
   }{%
    \aftergroup\false@sw
   }%
 }{%
  \aftergroup\false@sw
 }%
 \endgroup
}%
\def\output@moving{%
 \set@top@firstmark
 \@ifnum{\outputpenalty=\do@newpage@pen}{%
  \setbox\@cclv\vbox{%
   \unvbox\@cclv
   \remove@lastbox
   \@ifdim{\ht\z@=\ht\@protection@box}{\box\lastbox}{\unskip}%
  }%
 }{}%
 \@cclv@nontrivial@sw{%
  \expandafter\output@do@prep\csname output@prep@\bot@envir \endcsname
  \@makecolumn\true@sw
  \expandafter\output@column@do\csname output@column@\thepagegrid\endcsname
  \protect@penalty\do@startcolumn@pen
  \clearpage@sw{%
   \protect@penalty\do@endpage@pen
  }{}%
  \expandafter\let\expandafter\output@post@\csname output@post@\bot@envir \endcsname
  \outputdebug@sw{\say\output@post@}{}%
  \@ifx{\output@post@\relax}{\output@post@document}{\output@post@}%
 }{%
  \void@cclv
 }%
 \set@colht
 \global\@mparbottom\z@
 \global\@textfloatsheight\z@
}%
\def\output@do@prep#1{%
 \outputdebug@sw{\class@info{Prep: \string#1}}{}%
 \@ifx{#1\relax}{\output@prep@document}{#1}%
}%
\def\output@column@do#1{%
  \outputdebug@sw{\class@info{Output column: \string#1}}{}%
  \@ifx{#1\relax}{\output@column@one}{#1}%
}%
\def\void@cclv{\begingroup\setbox\z@\box\@cclv\endgroup}%
\def\remove@lastbox{\setbox\z@\lastbox}%
\def\@cclv@nontrivial@sw{%
 \@ifx@empty\@toplist{%
  \@ifx@empty\@botlist{%
   \@ifvoid\footins{%
    \@ifvoid\@cclv{%
     \false@sw
    }{%
     \setbox\z@\vbox{\unvcopy\@cclv}%
     \@ifdim{\ht\z@=\topskip}{%
      \setbox\z@\vbox\bgroup
       \unvbox\z@
       \remove@lastbox
       \dimen@\lastskip\unskip
       \@ifdim{\ht\z@=\ht\@protection@box}{%
        \advance\dimen@\ht\z@
        \@ifdim{\dimen@=\topskip}{%
         \aftergroup\true@sw
        }{%
         \aftergroup\false@sw
        }%
       }{%
        \aftergroup\false@sw
       }%
      \egroup
      {%
       \false@sw
      }{%
       \true@sw
      }%
     }{%
      \@ifdim{\ht\z@=\z@}{%
       \ltxgrid@info{Found trivial column. Discarding it}%
       \outputdebug@sw{\trace@box\@cclv}{}%
       \false@sw
      }{%
       \true@sw
      }%
     }%
    }%
   }{%
    \true@sw
   }%
  }{%
   \true@sw
  }%
 }{%
  \true@sw
 }%
}%
\def\protect@penalty#1{\protection@box\penalty-#1\relax}%
\newbox\@protection@box
\setbox\@protection@box\vbox to1986sp{\vfil}%
\def\protection@box{\nointerlineskip\copy\@protection@box}%
\def\dead@cycle@repair#1{%
 \expandafter\do@@mark
 \expandafter{%
              \@@botmark
             }%
 \unvbox\@cclv
 \nointerlineskip
 \vbox to#1{\vss}%
 \@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
}%
\def\dead@cycle@repair@protected#1{%
 \expandafter\do@@mark
 \expandafter{%
              \@@botmark
             }%
 \begingroup
  \unvbox\@cclv
  \remove@lastbox
  \nointerlineskip
  \advance#1-\ht\@protection@box
  \vbox to#1{\vss}%
  \protection@box % Reinsert protection box
  \@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
 \endgroup
}%
\def\dead@cycle{%
 \expandafter\do@@mark
 \expandafter{%
              \@@botmark
             }%
 \unvbox\@cclv
 \@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
}%
\def\output@init@document{%
 \ltxgrid@info@sw{\class@info{\string\output@init@document}}{}%
 \global\vsize\vsize
}%
\def\output@prep@document{%
 \ltxgrid@foot@info@sw{\class@info{\string\output@prep@document}\trace@scroll{\showbox\footins\showbox\footsofar}}{}%
 \@ifvoid\footsofar{%
 }{%
  \global\setbox\footins\vbox\bgroup
   \unvbox\footsofar
   \@ifvoid\footins{}{%
    \marry@baselines
    \unvbox\footins
   }%
  \egroup
  \ltxgrid@foot@info@sw{\trace@box\footins}{}%
 }%
}%
\def\output@post@document{}%
\let\@opcol\@undefined
\def\@makecolumn#1{%
 \ltxgrid@foot@info@sw{\class@info{\string\@makecolumn\string#1}}{}%
 \setbox\@outputbox\vbox\bgroup
  \boxmaxdepth\@maxdepth
  \@tempdima\dp\@cclv
  \unvbox\@cclv
  \vskip-\@tempdima
 \egroup
 \xdef\@freelist{\@freelist\@midlist}\global\let\@midlist\@empty
 \show@text@box@size
 \@combinefloats
 #1{%
  \@combineinserts\@outputbox\footins
 }{%
  \combine@foot@inserts\footsofar\footins
 }%
 \set@adj@colht\dimen@
 \count@\vbadness
 \vbadness\@M
 \setbox\@outputbox\vbox to\dimen@\bgroup
  \@texttop
  \dimen@\dp\@outputbox
  \unvbox\@outputbox
  \vskip-\dimen@
  \@textbottom
 \egroup
 \vbadness\count@
 \global\maxdepth\@maxdepth
}%
\let\@makespecialcolbox\@undefined
\def\@combineinserts#1#2{%
 \ltxgrid@foot@info@sw{\class@info{\string\@combineinserts\string#1\string#2}\trace@box#2}{}%
 \setbox#1\vbox\bgroup
  \unvbox#1%
  \@ifvoid{#2}{}{%
   \dimen@\ht#2\advance\dimen@\dp#2\advance\dimen@\skip#2%
   \show@box@size{Combining inserts}#2%
   \vskip\skip#2%
   \setbox\z@\vbox{\footnoterule}\dimen@i\ht\z@
   \color@begingroup
   \normalcolor
   \cleaders\box\z@\vskip\dimen@i\kern-\dimen@i
   \csname combine@insert@\the\pagegrid@col\endcsname#2%
   \color@endgroup
   \kern-\dimen@\kern\dimen@
  }%
 \egroup
 \ltxgrid@foot@info@sw{\trace@box#1}{}%
}%
\def\combine@insert@tw@#1{%
 \compose@footnotes@two#1\@ifvbox{#1}{\unvbox}{\box}#1%
}%
\def\combine@insert@@ne#1{%
 \compose@footnotes@one#1\@ifvbox{#1}{\unvbox}{\box}#1%
}%
\def\twocolumn@grid@setup{%
 \expandafter\let\csname combine@insert@1\endcsname\combine@insert@tw@
 \expandafter\let\csname combine@insert@2\endcsname\combine@insert@@ne
}%
\def\onecolumn@grid@setup{%
 \expandafter\let\csname combine@insert@1\endcsname\combine@insert@@ne
 \expandafter\let\csname combine@insert@2\endcsname\combine@insert@@ne
}%
\let\columngrid@setup\onecolumn@grid@setup
\columngrid@setup
\appdef\@floatplacement{%
 \global\@fpmin\@fpmin
}%
\mathchardef\pagebreak@pen=\@M
\expandafter\let\csname output@-\the\pagebreak@pen\endcsname\relax
\mathchardef\do@startcolumn@pen=10005
\@namedef{output@-\the\do@startcolumn@pen}{\do@startcolumn}%
\def\do@startcolumn{%
 \setbox\@cclv\vbox{\unvbox\@cclv\remove@lastbox\unskip}%
 \clearpage@sw{\@clearfloatplacement}{\@floatplacement}%
 \set@colht
 \@booleanfalse\pfloat@avail@sw
 \begingroup
  \@colht\@colroom
  \@booleanfalse\float@avail@sw
  \@tryfcolumn\test@colfloat
  \float@avail@sw{\aftergroup\@booleantrue\aftergroup\pfloat@avail@sw}{}%
 \endgroup
 \fcolmade@sw{%
  \setbox\@cclv\vbox{\unvbox\@outputbox\unvbox\@cclv}%
  \outputpenalty-\pagebreak@pen
  \dead@cycle
 }{%
  \begingroup
   \let\@elt\@scolelt
   \let\reserved@b\@deferlist\global\let\@deferlist\@empty\reserved@b
  \endgroup
  \clearpage@sw{%
   \outputpenalty\@M
  }{%
   \outputpenalty\do@newpage@pen
  }%
  \dead@cycle
 }%
 \check@deferlist@stuck\do@startcolumn
 \set@vsize
}%
\def\@scolelt#1{\def\@currbox{#1}\@addtonextcol}%
\def\test@colfloat#1{%
 \csname @floatselect@sw@\thepagegrid\endcsname#1{}{\@testtrue}%
 \@if@sw\if@test\fi{}{\aftergroup\@booleantrue\aftergroup\float@avail@sw}%
}%
\def\@addtonextcol{%
 \begingroup
  \@insertfalse
  \@setfloattypecounts
  \csname @floatselect@sw@\thepagegrid\endcsname\@currbox{%
   \@ifnum{\@fpstype=8 }{}{%
     \@ifnum{\@fpstype=24 }{}{%
       \@flsettextmin
       \@reqcolroom \ht\@currbox
       \advance \@reqcolroom \@textmin
       \advance \@reqcolroom \vsize % take into account split insertions
       \advance \@reqcolroom -\pagegoal
       \@ifdim{\@colroom>\@reqcolroom}{%
         \@flsetnum \@colnum
         \@ifnum{\@colnum>\z@}{%
            \@bitor\@currtype\@deferlist
            \@if@sw\if@test\fi{}{%
              \@addtotoporbot
            }%
         }{}%
       }{}%
     }%
   }%
  }{}%
  \@if@sw\if@insert\fi{}{%
    \@cons\@deferlist\@currbox
  }%
 \endgroup
}%
\mathchardef\do@startpage@pen=10006
\@namedef{output@-\the\do@startpage@pen}{\do@startpage}%
\def\do@startpage{%
 \setbox\@cclv\vbox{\unvbox\@cclv\remove@lastbox\unskip}%
 \clearpage@sw{\@clearfloatplacement}{\@dblfloatplacement}%
 \set@colht
 \@booleanfalse\pfloat@avail@sw
 \begingroup
  \@booleanfalse\float@avail@sw
  \@tryfcolumn\test@dblfloat
  \float@avail@sw{\aftergroup\@booleantrue\aftergroup\pfloat@avail@sw}{}%
 \endgroup
 \fcolmade@sw{%
  \global\setbox\pagesofar\vbox{\unvbox\pagesofar\unvbox\@outputbox}%
  \@output@combined@page
 }{%
  \begingroup
   \@booleanfalse\float@avail@sw
   \let\@elt\@sdblcolelt
   \let\reserved@b\@deferlist\global\let\@deferlist\@empty\reserved@b
  \endgroup
  \@ifdim{\@colht=\textheight}{% No luck...
   \pfloat@avail@sw{% ...but a float *was* available!
    \forcefloats@sw{%
     \ltxgrid@warn{Forced dequeueing of floats stalled}%
    }{%
     \ltxgrid@warn{Dequeueing of floats stalled}%
    }%
   }{}%
  }{}%
  \outputpenalty\@M
  \dead@cycle
 }%
 \check@deferlist@stuck\do@startpage
 \set@colht
}%
\def\@output@combined@page{%
 \@combinepage\true@sw
 \@combinedblfloats
 \@outputpage
 \global\pagegrid@cur\@ne
 \protect@penalty\do@startpage@pen
}%
\def\@sdblcolelt#1{\def\@currbox{#1}\@addtodblcol}%
\def\test@dblfloat#1{%
 \@if@notdblfloat{#1}{\@testtrue}{}%
 \@if@sw\if@test\fi{}{\aftergroup\@booleantrue\aftergroup\float@avail@sw}%
}%
\def\@if@notdblfloat#1{\@ifdim{\wd#1<\textwidth}}%
\@booleanfalse\forcefloats@sw
\def\@addtodblcol{%
 \begingroup
  \@if@notdblfloat{\@currbox}{%
   \false@sw
  }{%
   \@setfloattypecounts
   \@getfpsbit \tw@
   \@bitor \@currtype \@deferlist
   \@if@sw\if@test\fi{%
    \false@sw
   }{%
    \@ifodd\@tempcnta{%
     \aftergroup\@booleantrue\aftergroup\float@avail@sw
     \@flsetnum \@dbltopnum
     \@ifnum{\@dbltopnum>\z@}{%
       \@ifdim{\@dbltoproom>\ht\@currbox}{%
        \true@sw
       }{%
        \@ifnum{\@fpstype<\sixt@@n}{%
         \begingroup
          \advance \@dbltoproom \@textmin
          \@ifdim{\@dbltoproom>\ht\@currbox}{%
           \endgroup\true@sw
          }{%
           \endgroup\false@sw
          }%
        }{%
         \false@sw
        }%
       }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }%
  }%
  {%
   \@tempdima -\ht\@currbox
   \advance\@tempdima
    -\@ifx{\@dbltoplist\@empty}{\dbltextfloatsep}{\dblfloatsep}%
   \global \advance \@dbltoproom \@tempdima
   \global \advance \@colht \@tempdima
   \global \advance \@dbltopnum \m@ne
   \@cons \@dbltoplist \@currbox
  }{%
   \@cons \@deferlist \@currbox
  }%
 \endgroup
}%
\def\@tryfcolumn#1{%
  \global\@booleanfalse\fcolmade@sw
  \@ifx@empty\@deferlist{}{%
    \global\let\@trylist\@deferlist
    \global\let\@failedlist\@empty
    \begingroup
      \dimen@\vsize\advance\dimen@-\pagegoal\@ifdim{\dimen@>\z@}{%
       \advance\@fpmin-\dimen@
      }{}%
      \def\@elt{\@xtryfc#1}\@trylist
    \endgroup
    \fcolmade@sw{%
      \global\setbox\@outputbox\vbox{\vskip \@fptop}%
      \let \@elt \@wtryfc \@flsucceed
      \global\setbox\@outputbox\vbox{\unvbox\@outputbox
        \unskip \vskip \@fpbot
      }%
      \let \@elt \relax
      \xdef\@deferlist{\@failedlist\@flfail}%
      \xdef\@freelist{\@freelist\@flsucceed}%
    }{}%
  }%
}%
\def\@wtryfc #1{%
  \global\setbox\@outputbox\vbox{\unvbox\@outputbox
    \box #1\vskip\@fpsep
  }%
}%
\def\@xtryfc#1#2{%
  \@next\reserved@a\@trylist{}{}% trim \@trylist. Ugly!
  \@currtype \count #2%
  \divide\@currtype\@xxxii\multiply\@currtype\@xxxii
  \@bitor \@currtype \@failedlist
  \@testfp #2%
  #1#2%
  \@ifdim{\ht #2>\@colht   }{\@testtrue}{}%
  \@if@sw\if@test\fi{%
   \@cons\@failedlist #2%
  }{%
   \begingroup
     \gdef\@flsucceed{\@elt #2}%
     \global\let\@flfail\@empty
     \@tempdima\ht #2%
     \def \@elt {\@ztryfc#1}\@trylist
     \@ifdim{\@tempdima >\@fpmin}{%
       \global\@booleantrue\fcolmade@sw
     }{%
       \@cons\@failedlist #2%
     }%
   \endgroup
   \fcolmade@sw{%
     \let \@elt \@gobble
   }{}%
  }%
}%
\def\@ztryfc #1#2{%
  \@tempcnta \count#2%
  \divide\@tempcnta\@xxxii\multiply\@tempcnta\@xxxii
  \@bitor \@tempcnta {\@failedlist \@flfail}%
  \@testfp #2%
  #1#2%
  \@tempdimb\@tempdima
  \advance\@tempdimb \ht#2\advance\@tempdimb\@fpsep
  \@ifdim{\@tempdimb >\@colht}{%
    \@testtrue
  }{}%
  \@if@sw\if@test\fi{%
    \@cons\@flfail #2%
  }{%
    \@cons\@flsucceed #2%
    \@tempdima\@tempdimb
  }%
}%
\def\newpage@prep{%
  \if@noskipsec
    \ifx \@nodocument\relax
      \leavevmode
      \global \@noskipsecfalse
    \fi
  \fi
  \if@inlabel
    \leavevmode
    \global \@inlabelfalse
  \fi
  \if@nobreak \@nobreakfalse \everypar{}\fi
  \par
}%
\def \newpage {%
 \newpage@prep
 \do@output@MVL{%
  \vfil
  \penalty-\pagebreak@pen
 }%
}%
\def\clearpage{%
 \newpage@prep
 \do@output@MVL{%
  \vfil
  \penalty-\pagebreak@pen
  \global\@booleantrue\clearpage@sw
  \protect@penalty\do@startcolumn@pen
  \protect@penalty\do@endpage@pen
 }%
 \do@output@MVL{%
  \global\@booleanfalse\clearpage@sw
 }%
}%
\def\cleardoublepage{%
 \clearpage
 \@if@sw\if@twoside\fi{%
  \@ifodd\c@page{}{%
   \null\clearpage
  }%
 }{}%
}%
\@booleanfalse\clearpage@sw
\mathchardef\do@endpage@pen=10007
\@namedef{output@-\the\do@endpage@pen}{\csname end@column@\thepagegrid\endcsname}%
\mathchardef\do@newpage@pen=10001
\expandafter\let\csname output@-\the\do@newpage@pen\endcsname\relax
\def\@clearfloatplacement{%
 \global\@topnum     \maxdimen
 \global\@toproom    \maxdimen
 \global\@botnum     \maxdimen
 \global\@botroom    \maxdimen
 \global\@colnum     \maxdimen
 \global\@dbltopnum  \maxdimen
 \global\@dbltoproom \maxdimen
 \global\@textmin    \z@
 \global\@fpmin      \z@
 \let\@testfp\@gobble
 \appdef\@setfloattypecounts{\@fpstype16\advance\@fpstype\m@ne}%
}%
\let\@doclearpage\@undefined
\let\@makefcolumn\@undefined
\let\@makecol\@undefined
\def\clr@top@firstmark{%
 \global\let\saved@@topmark\@undefined
 \global\let\saved@@firstmark\@empty
 \global\let\saved@@botmark\@empty
}%
\clr@top@firstmark
\def\set@top@firstmark{%
 \@ifxundefined\saved@@topmark{\expandafter\gdef\expandafter\saved@@topmark\expandafter{\@@topmark}}{}%
 \@if@empty\saved@@firstmark{\expandafter\gdef\expandafter\saved@@firstmark\expandafter{\@@firstmark}}{}%
 \@if@empty\@@botmark{}{\expandafter\gdef\expandafter\saved@@botmark\expandafter{\@@botmark}}%
}%
\appdef\@outputpage@tail{%
 \clr@top@firstmark
}%
\def\@float#1{%
 \@ifnextchar[{%
  \@yfloat\width@float{#1}%
 }{%
  \@ifxundefined@cs{fps@#1}{}{\expandafter\let\expandafter\fps@\csname fps@#1\endcsname}%
  \expandafter\@argswap\expandafter{\expandafter[\fps@]}{\@yfloat\width@float{#1}}%
 }%
}%
\def\@dblfloat#1{%
 \@ifnum{\pagegrid@col=\@ne}{%
  \@float{#1}%
 }{%
  \@ifnextchar[{%
   \@yfloat\widthd@float{#1}%
  }{%
   \@ifxundefined@cs{fpsd@#1}{}{\expandafter\let\expandafter\fpsd@\csname fpsd@#1\endcsname}%
   \expandafter\@argswap\expandafter{\expandafter[\fpsd@]}{\@yfloat\widthd@float{#1}}%
  }%
 }%
}%
\def\@yfloat#1#2[#3]{%
 \@xfloat{#2}[#3]%
 \hsize#1\linewidth\hsize
 \let\set@footnotewidth\@empty
 \minipagefootnote@init
}%
\def\fps@{tbp}%
\def\fpsd@{tp}%
\def\width@float{\columnwidth}%
\def\widthd@float{\textwidth}%
\def\end@float{%
 \end@@float{%
  \check@currbox@count
 }%
}%
\def\end@dblfloat{%
 \@ifnum{\pagegrid@col=\@ne}{%
  \end@float
 }{%
  \end@@float{%
   \@iffpsbit\@ne{\global\advance\count\@currbox\m@ne}{}%
   \@iffpsbit\f@ur{\global\advance\count\@currbox-4\relax}{}%
   \global\wd\@currbox\textwidth % Klootch
   \check@currbox@count
  }%
 }%
}%
\def\end@@float#1{%
 \minipagefootnote@here
 \@endfloatbox
 #1%
 \@ifnum{\@floatpenalty <\z@}{%
  \@largefloatcheck
  \@cons\@currlist\@currbox
  \@ifnum{\@floatpenalty <-\@Mii}{%
   \do@output@cclv{\@add@float}%
  }{%
   \vadjust{\do@output@cclv{\@add@float}}%
   \@Esphack
  }%
 }{}%
}%
\newcommand\float@end@float{%
 \@endfloatbox
 \global\setbox\@currbox\float@makebox\columnwidth
 \let\@endfloatbox\relax
 \end@float
}%
\newcommand\float@end@ltx{%
 \end@@float{%
  \global\setbox\@currbox\float@makebox\columnwidth
  \check@currbox@count
 }%
}%
\newcommand\newfloat@float[3]{%
 \@namedef{ext@#1}{#3} %!
 \let\float@do=\relax
 \xdef\@tempa{\noexpand\float@exts{\the\float@exts \float@do{#3}}}%
 \@tempa
 \floatplacement{#1}{#2}%
 \@ifundefined{fname@#1}{\floatname{#1}{#1}}{} %!
 \expandafter\edef\csname ftype@#1\endcsname{\value{float@type}}%
 \addtocounter{float@type}{\value{float@type}} %!
 \restylefloat{#1}%
 \expandafter\edef\csname fnum@#1\endcsname{%
  \expandafter\noexpand\csname fname@#1\endcsname{} %!
  \expandafter\noexpand\csname the#1\endcsname
 }
 \@ifnextchar[%]
  {%
   \float@newx{#1}%
  }{%
   \@ifundefined{c@#1}{\newcounter{#1}\@namedef{the#1}{\arabic{#1}}}{}%
  }%
}%
\newcommand\newfloat@ltx[3]{%
 \@namedef{ext@#1}{#3}%
 \let\float@do=\relax
 \xdef\@tempa{\noexpand\float@exts{\the\float@exts \float@do{#3}}}%
 \@tempa
 \floatplacement{#1}{#2}%
 \@ifundefined{fname@#1}{\floatname{#1}{#1}}{}%
 \expandafter\edef\csname ftype@#1\expandafter\endcsname\expandafter{\the\c@float@type}%
 \addtocounter{float@type}{\value{float@type}}%
 \restylefloat{#1}%
 \expandafter\edef\csname fnum@#1\endcsname{%
  \expandafter\noexpand\csname fname@#1\endcsname{}%
  \expandafter\noexpand\csname the#1\endcsname
 }
 \@ifnextchar[%]
  {%
   \float@newx{#1}%
  }{%
   \@ifundefined{c@#1}{\newcounter{#1}\@namedef{the#1}{\arabic{#1}}}{}%
  }%
}%
\appdef\document@inithook{%
 \@ifxundefined\newfloat{}{%
  \@ifx{\float@end\float@end@float}{%
   \@ifx{\newfloat\newfloat@float}{\true@sw}{\false@sw}%
   }{\false@sw}%
   {%
    \class@warn{Repair the float package}%
    \let\float@end\float@end@ltx
    \let\newfloat\newfloat@ltx
   }{%
    \class@warn{Failed to patch the float package}%
   }%
 }%
}%
\def\@iffpsbit#1{%
 \begingroup
  \@tempcnta\count\@currbox
  \divide\@tempcnta#1\relax
  \@ifodd\@tempcnta{\aftergroup\true@sw}{\aftergroup\false@sw}%
 \endgroup
}%
\def\check@currbox@count{%
 \@ifnum{\count\@currbox>\z@}{%
  \count@\count\@currbox\divide\count@\sixt@@n\multiply\count@\sixt@@n
  \@tempcnta\count\@currbox\advance\@tempcnta-\count@
  \@ifnum{\@tempcnta=\z@}{%
   \ltxgrid@warn{Float cannot be placed}%
  }{}%
  \expandafter\tally@float\expandafter{\@captype}%
 }{%
 }%
}%
\providecommand\minipagefootnote@init{}%
\providecommand\minipagefootnote@here{}%
\providecommand\tally@float[1]{}%
\let\@specialoutput\@undefined
\def\@add@float{%
 \@pageht\ht\@cclv\@pagedp\dp\@cclv
 \unvbox\@cclv
 \@next\@currbox\@currlist{%
  \csname @floatselect@sw@\thepagegrid\endcsname\@currbox{%
   \@ifnum{\count\@currbox>\z@}{%
    \advance \@pageht \@pagedp
    \advance \@pageht \vsize \advance \@pageht -\pagegoal
    \@addtocurcol
   }{%
    \@addmarginpar
   }%
  }{%
   \@resethfps
   \@cons\@deferlist\@currbox
  }%
 }{\@latexbug}%
 \@ifnum{\outputpenalty<\z@}{%
  \@if@sw\if@nobreak\fi{%
   \nobreak
  }{%
   \addpenalty \interlinepenalty
  }%
 }{}%
 \set@vsize
}%
\let\@reinserts\@undefined
\def \@addtocurcol {%
   \@insertfalse
   \@setfloattypecounts
   \ifnum \@fpstype=8
   \else
     \ifnum \@fpstype=24
     \else
       \@flsettextmin
       \advance \@textmin \@textfloatsheight
       \@reqcolroom \@pageht
       \ifdim \@textmin>\@reqcolroom
         \@reqcolroom \@textmin
       \fi
       \advance \@reqcolroom \ht\@currbox
       \ifdim \@colroom>\@reqcolroom
         \@flsetnum \@colnum
         \ifnum \@colnum>\z@
           \@bitor\@currtype\@deferlist
           \if@test
           \else
             \@bitor\@currtype\@botlist
             \if@test
               \@addtobot
             \else
               \ifodd \count\@currbox
                 \advance \@reqcolroom \intextsep
                 \ifdim \@colroom>\@reqcolroom
                   \global \advance \@colnum \m@ne
                   \global \advance \@textfloatsheight \ht\@currbox
                   \global \advance \@textfloatsheight 2\intextsep
                   \@cons \@midlist \@currbox
                   \if@nobreak
                     \nobreak
                     \@nobreakfalse
                     \everypar{}%
                   \else
                     \addpenalty \interlinepenalty
                   \fi
                   \vskip \intextsep
                   \unvbox\@currbox %AO
                   \penalty\interlinepenalty
                   \vskip\intextsep
                   \ifnum\outputpenalty <-\@Mii \vskip -\parskip\fi
                   \outputpenalty \z@
                   \@inserttrue
                 \fi
               \fi
               \if@insert
               \else
                 \@addtotoporbot
               \fi
             \fi
           \fi
         \fi
       \fi
     \fi
   \fi
   \if@insert
   \else
     \@resethfps
     \@cons\@deferlist\@currbox
   \fi
}%
\@twocolumnfalse
\let\@twocolumntrue\@twocolumnfalse
\def\@addmarginpar{%
 \@next\@marbox\@currlist{%
  \@cons\@freelist\@marbox\@cons\@freelist\@currbox
 }\@latexbug
 \setbox\@marbox\hb@xt@\columnwidth{%
  \csname @addmarginpar@\thepagegrid\endcsname{%
   \hskip-\marginparsep\hskip-\marginparwidth
   \box\@currbox
  }{%
   \hskip\columnwidth\hskip\marginparsep
   \box\@marbox
  }%
  \hss
 }%
 \setbox\z@\box\@currbox
    \@tempdima\@mparbottom
    \advance\@tempdima -\@pageht
    \advance\@tempdima\ht\@marbox
 \@ifdim{\@tempdima >\z@}{%
   \@latex@warning@no@line {Marginpar on page \thepage\space moved}%
 }{%
   \@tempdima\z@
 }%
    \global\@mparbottom\@pageht
    \global\advance\@mparbottom\@tempdima
    \global\advance\@mparbottom\dp\@marbox
    \global\advance\@mparbottom\marginparpush
    \advance\@tempdima -\ht\@marbox
    \global\setbox \@marbox
                   \vbox {\vskip \@tempdima
                          \box \@marbox}%
    \global \ht\@marbox \z@
    \global \dp\@marbox \z@
    \kern -\@pagedp
    \nointerlineskip
  \box\@marbox
    \nointerlineskip
    \hbox{\vrule \@height\z@ \@width\z@ \@depth\@pagedp}%
}%
\newenvironment{turnpage}{%
 \def\width@float{\textheight}%
 \def\widthd@float{\textheight}%
 \appdef\@endfloatbox{%
  \@ifxundefined\@currbox{%
   \ltxgrid@warn{Cannot rotate! Not a float}%
  }{%
   \setbox\@currbox\vbox to\textwidth{\vfil\unvbox\@currbox\vfil}%
   \global\setbox\@currbox\vbox{\rotatebox{90}{\box\@currbox}}%
  }%
 }%
}{%
}%
\def\rotatebox@dummy#1#2{%
 \ltxgrid@warn{You must load the graphics or graphicx package in order to use the turnpage environment}%
 #2%
}%
\appdef\document@inithook{%
 \@ifxundefined\rotatebox{\let\rotatebox\rotatebox@dummy}{}%
}%
\@namedef{output@-1073741824}{%
 \deadcycles\z@
 \void@cclv
}%
\mathchardef\save@column@pen=10016
\@namedef{output@-\the\save@column@pen}{\save@column}%
\let \@cclv@saved \@holdpg
\let \@holdpg \@undefined
\def\save@column{%
 \@ifvoid\@cclv@saved{%
  \set@top@firstmark
  \global\@topmark@saved\expandafter{\@@topmark}%
 }{}%
 \global\setbox\@cclv@saved\vbox{%
  \@ifvoid\@cclv@saved{}{%
   \unvbox\@cclv@saved
   \marry@baselines
  }%
  \unvbox\@cclv
  \lose@breaks
  \remove@lastbox
 }%
}%
\newtoks\@topmark@saved
\def\prep@cclv{%
 \void@cclv
 \setbox\@cclv\box\@cclv@saved
 \vbadness\@M
}%
\mathchardef\save@column@insert@pen=10017
\@namedef{output@-\the\save@column@insert@pen}{\toggle@insert{\savecolumn@holding}{\savecolumn@moving}}%
\def\savecolumn@holding{%
 \@if@exceed@pagegoal{\unvcopy\@cclv\remove@lastbox}{%
  \setbox\z@\vbox{\unvcopy\@cclv\remove@lastbox}%
  \outputdebug@sw{\trace@box\z@}{}%
  \dimen@\ht\@cclv\advance\dimen@-\ht\z@
  \dead@cycle@repair@protected\dimen@
 }{%
  \dead@cycle
 }%
}%
\def\savecolumn@moving{%
 \ltxgrid@info@sw{\class@info{\string\savecolumn@moving}}{}%
 \@cclv@nontrivial@sw{%
  \save@column
 }{%
  \void@cclv
 }%
 \@ifvoid\footins{}{%
  \ltxgrid@foot@info@sw{\class@info{\string\savecolumn@moving}\trace@scroll{\showbox\footins@saved\showbox\footins}}{}%
  \@ifvoid\footins@saved{%
   \global\setbox\footins@saved\box\footins
  }{%
   \global\setbox\footins@saved\vbox\bgroup
    \unvbox\footins@saved
    \marry@baselines
    \unvbox\footins
   \egroup
  }%
  \ltxgrid@foot@info@sw{\trace@box\footins@saved}{}%
  \protect@penalty\save@column@insert@pen
 }%
}%
\newbox\footins@saved
\newbox\footins@recovered
\newbox\column@recovered
\mathchardef\save@message@pen=10018
\@namedef{output@-\the\save@message@pen}{\save@message}%
\def\save@message{%
 \void@cclv
 \toks@\expandafter{\@@firstmark}%
 \expandafter\gdef\expandafter\@message@saved\expandafter{\the\toks@}%
 \expandafter\do@@mark\expandafter{\the\@topmark@saved}%
}%
\gdef\@message@saved{}%
\mathchardef\execute@message@pen=10019
\@namedef{output@-\the\execute@message@pen}{\@message@saved}%
\def\execute@message{%
 \@execute@message\save@column@pen
}%
\def\execute@message@insert#1{%
 \@execute@message\save@column@insert@pen{%
  \setbox \footins \box \footins@saved
  \ltxgrid@foot@info@sw{\class@info{\string\execute@message@insert}\trace@box\footins}{}%
  #1%
 }%
}%
\long\def\@execute@message#1#2{%
 \begingroup
  \dimen@\prevdepth\@ifdim{\dimen@<\z@}{\dimen@\z@}{}%
  \setbox\z@\vbox{%
   \protect@penalty#1%
   \protection@box
   \toks@{\prep@cclv#2}%
   \@@mark{\the\toks@}%
   \penalty-\save@message@pen
   \setbox\z@\null\dp\z@\dimen@\ht\z@-\dimen@
   \nointerlineskip\box\z@
   \penalty-\execute@message@pen
  }\unvbox\z@
 \endgroup
}%
\def\do@output@cclv{\execute@message}%
\def\do@output@MVL#1{%
 \@ifvmode{%
  \begingroup\execute@message{\unvbox\@cclv#1}\endgroup
 }{%
  \@ifhmode{%
   \vadjust{\execute@message{\unvbox\@cclv#1}}%
  }{%
   \@latexerr{\string\do@output@MVL\space cannot be executed in this mode!}\@eha
  }%
 }%
}%
\def\lose@breaks{%
 \loopwhile{%
  \count@\lastpenalty
  \@ifnum{\count@=\@M}{%
   \unpenalty\true@sw
  }{%
   \false@sw
  }%
 }%
}%
\def\removestuff{\do@output@MVL{\unskip\unpenalty}}%
\def\removephantombox{%
 \vadjust{%
  \execute@message{%
   \unvbox\@cclv
   \remove@lastbox
   \unskip
   \unskip
   \unpenalty
   \penalty\predisplaypenalty
   \vskip\abovedisplayskip
  }%
 }%
}%
\def\addstuff#1#2{\edef\@tempa{\noexpand\do@output@MVL{\noexpand\@addstuff{#1}{#2}}}\@tempa}%
\def\@addstuff#1#2{%
 \skip@\lastskip\unskip
 \count@\lastpenalty\unpenalty
 \@if@empty{#1}{}{\penalty#1\relax}%
 \@ifnum{\count@=\z@}{}{\penalty\count@}%
 \vskip\skip@
 \@if@empty{#2}{}{\vskip#2\relax}%
}%
\def\replacestuff#1#2{\edef\@tempa{\noexpand\do@output@MVL{\noexpand\@replacestuff{#1}{#2}}}\@tempa}%
\def\@replacestuff#1#2{%
 \skip@\lastskip\unskip
 \count@\lastpenalty\unpenalty
 \@if@empty{#1}{}{%
 \@ifnum{\count@>\@M}{}{%
   \@ifnum{\count@=\z@}{\count@=#1\relax}{%
    \@ifnum{\count@<#1\relax}{}{%
     \count@=#1\relax
    }%
   }%
 }%
 }%
 \@ifnum{\count@=\z@}{}{\penalty\count@}%
 \@if@empty{#2}{}{%
  \@tempskipa#2\relax
  \@ifdim{\z@>\@tempskipa}{%
   \advance\skip@-\@tempskipa
  }{%
   \@ifdim{\skip@>\@tempskipa}{}{%
    \skip@\@tempskipa
   }%
  }%
 }%
 \vskip\skip@
}%
\def\move@insertions{\global\holdinginserts\z@}%
\def\hold@insertions{\global\holdinginserts\@ne}%
\hold@insertions
\def\toggle@insert#1#2{%
 \@ifnum{\holdinginserts>\z@}{\move@insertions#1}{\hold@insertions#2}%
}%
\def\do@columngrid#1#2{%
 \par
 \expandafter\let\expandafter\@tempa\csname open@column@#1\endcsname
 \@ifx{\relax\@tempa}{%
  \ltxgrid@warn{Unknown page grid #1. No action taken}%
 }{%
  \do@output@MVL{\start@column{#1}{#2}}%
 }%
}%
\def\start@column#1#2{%
 \def\@tempa{#1}\@ifx{\@tempa\thepagegrid}{%
  \ltxgrid@info{Already in page grid \thepagegrid. No action taken}%
 }{%
  \expandafter\execute@message@insert
  \expandafter{%
               \csname shut@column@\thepagegrid\expandafter\endcsname
               \csname open@column@#1\endcsname{#2}%
               \set@vsize
             }%
 }%
}%
\def\thepagegrid{one}%
\newbox\pagesofar
\newbox\footsofar
\def\combine@foot@inserts#1#2{%
  \ltxgrid@info@sw{\class@info{\string\combine@foot@inserts\string#1\string#2}}{}%
  \@ifvoid#1{%
    \ltxgrid@foot@info@sw{\trace@box#2}{}\global\setbox#1\box#2%
  }{%
   \global\setbox#1\vbox\bgroup
    \ltxgrid@foot@info@sw{\trace@box#1}{}\unvbox#1%
    \@ifvoid#2{}{%
     \marry@baselines
     \ltxgrid@foot@info@sw{\trace@box#2}{}\unvbox#2%
    }%
   \egroup
  }%
  \ltxgrid@foot@info@sw{\trace@scroll{\showbox#1\showbox#2}}{}%
}%
\newcommand\onecolumngrid{\do@columngrid{one}{\@ne}}%
\let\onecolumn\@undefined
\def\open@column@one#1{%
 \ltxgrid@info@sw{\class@info{\string\open@column@one\string#1}}{}%
 \unvbox\pagesofar
 \@ifvoid{\footsofar}{}{%
  \insert\footins\bgroup\unvbox\footsofar\egroup
  \penalty\z@
 }%
 \gdef\thepagegrid{one}%
 \global\pagegrid@col#1%
 \global\pagegrid@cur\@ne
 \global\count\footins\@m
 \global\divide\count\footins\tw@
 \set@column@hsize\pagegrid@col
 \set@colht
}%
\def\shut@column@one{%
 \ltxgrid@info@sw{\class@info{\string\shut@column@one}}{}%
 \@makecolumn\false@sw
 \global\setbox\pagesofar\vbox\bgroup
  \recover@column\@outputbox\footsofar\column@recovered\footins@recovered
 \egroup
 \begingroup\setbox\z@\box\@outputbox\endgroup
 \combine@foot@inserts\footsofar\footins
 \set@colht
}%
\def\float@column@one{%
 \@makecolumn\true@sw
 \@outputpage
}%
\def\end@column@one{%
 \unvbox\@cclv\remove@lastbox
 \protect@penalty\do@newpage@pen
}%
\def\output@column@one{%
 \@outputpage
}%
\def\@addmarginpar@one{%
 \@if@sw\if@mparswitch\fi{%
  \@ifodd\c@page{\false@sw}{\true@sw}%
 }{\false@sw}{%
  \@if@sw\if@reversemargin\fi{\false@sw}{\true@sw}%
 }{%
  \@if@sw\if@reversemargin\fi{\true@sw}{\false@sw}%
 }%
}%
\def\@floatselect@sw@one#1{\true@sw}%
\def\onecolumngrid@push{%
 \do@output@MVL{%
  \@ifnum{\pagegrid@col=\@ne}{%
   \global\let\restorecolumngrid\@empty
  }{%
   \xdef\restorecolumngrid{%
    \noexpand\start@column{\thepagegrid}{\the\pagegrid@col}%
   }%
   \start@column{one}{\@ne}%
  }%
 }%
}%
\def\onecolumngrid@pop{%
 \do@output@MVL{\restorecolumngrid}%
}%
\newcommand\twocolumngrid{\do@columngrid{mlt}{\tw@}}%
\let\twocolumn\@undefined
\let\@topnewpage\@undefined
\def\open@column@mlt#1{%
 \ltxgrid@info@sw{\class@info{\string\open@column@mlt\string#1}}{}%
 \@ifvoid{\footsofar}{}{%
  \insert\footins\bgroup\unvbox\footsofar\egroup
 }%
 \gdef\thepagegrid{mlt}%
 \global\pagegrid@col#1%
 \global\pagegrid@cur\@ne
 \global\count\footins\@m
 \set@column@hsize\pagegrid@col
 \set@colht
}%
\def\shut@column@mlt{%
 \ltxgrid@info@sw{\class@info{\string\shut@column@mlt}}{}%
 \@cclv@nontrivial@sw{%
  \@makecolumn\false@sw
  \@ifnum{\pagegrid@cur<\pagegrid@col}{%
   \expandafter\global\expandafter\setbox\csname col@\the\pagegrid@cur\endcsname\box\@outputbox
   \global\advance\pagegrid@cur\@ne
  }{}%
 }{%
  \void@cclv
 }%
 \@ifnum{\pagegrid@cur>\@ne}{%
  \csname balance@\the\pagegrid@col\endcsname
  \grid@column\@outputbox{}%
  \@combinepage\false@sw
  \@combinedblfloats
  \global\setbox\pagesofar\box\@outputbox
  \show@pagesofar@size
 }{}%
 \set@colht
}%
\def\float@column@mlt{%
  \@output@combined@page
}%
\def\end@column@mlt{%
 \@ifx@empty\@toplist{%
  \@ifx@empty\@botlist{%
   \@ifx@empty\@dbltoplist{%
    \@ifx@empty\@deferlist{%
     \@ifnum{\pagegrid@cur=\@ne}{%
      \false@sw
     }{%
      \true@sw
     }%
    }{%
     \true@sw
    }%
   }{%
    \true@sw
   }%
  }{%
   \true@sw
  }%
 }{%
  \true@sw
 }%
 % true = kick out a column and try again
 {%
  \@cclv@nontrivial@sw{%
   \unvbox\@cclv\remove@lastbox
  }{%
   \unvbox\@cclv\remove@lastbox\unskip\null
  }%
  \protect@penalty\do@newpage@pen
  \protect@penalty\do@endpage@pen
 }{%
  \unvbox\@cclv\remove@lastbox
 }%
}%
\def\output@column@mlt{%
 \@ifnum{\pagegrid@cur<\pagegrid@col}{%
  \expandafter\global\expandafter\setbox\csname col@\the\pagegrid@cur\endcsname\box\@outputbox
  \global\advance\pagegrid@cur\@ne
 }{%
  \set@adj@colht\dimen@
  \grid@column\@outputbox{}%
  \@output@combined@page
 }%
}%
\let\@outputdblcol\@undefined
\def\@floatselect@sw@mlt#1{\@if@notdblfloat{#1}}%
\def\@addmarginpar@mlt{% emits a boolean
 \@ifnum{\pagegrid@cur=\@ne}%
}%
\def\set@footnotewidth@one{%
 \hsize\columnwidth
 \linewidth\hsize
}%
\def\set@footnotewidth@two{\set@footnotewidth@mlt\tw@}%
\def\set@footnotewidth@mlt#1{%
 \hsize\textwidth
 \advance\hsize\columnsep
 \divide\hsize#1%
 \advance\hsize-\columnsep
 \linewidth\hsize
}%
\def\compose@footnotes@one#1{%
 \ltxgrid@foot@info@sw{\class@info{\string\compose@footnotes@one\string#1}\trace@box#1}{}%
}%
\let\compose@footnotes\compose@footnotes@one
\def\compose@footnotes@two#1{%
 \ltxgrid@foot@info@sw{\class@info{\string\compose@footnotes@two\string#1}\trace@box#1}{}%
 \setbox\z@\box\@tempboxa
 \let\recover@column\recover@column@null
 \let\marry@baselines\@empty
 \balance@two#1\@tempboxa
 \global\setbox#1\hbox to\textwidth{\box#1\hfil\box\@tempboxa}%
 \ltxgrid@foot@info@sw{\trace@box#1}{}%
}%
\let\pagegrid@cur\col@number
\let\col@number\@undefined
\newcount\pagegrid@col
\pagegrid@cur\@ne
\expandafter\let\csname col@\the\pagegrid@cur\endcsname\@leftcolumn
\let\@leftcolumn\@undefined
\pagegrid@col\tw@
\def\pagegrid@init{%
 \advance\pagegrid@cur\@ne
 \@ifnum{\pagegrid@cur<\pagegrid@col}{%
  \csname newbox\expandafter\endcsname\csname  col@\the\pagegrid@cur\endcsname
  \pagegrid@init
 }{%
 }%
}%
\appdef\class@documenthook{%
 \pagegrid@init
}%
\def\grid@column#1#2{%
 \ltxgrid@info@sw{\class@info{\string\grid@column\string#1}}{}%
 \global\setbox#1\vbox\bgroup
  \hb@xt@\textwidth\bgroup
   \vrule\@height\z@\@width\z@\@if@empty{#2}{}{\@depth#2}%
   \pagegrid@cur\@ne
   \@ifnum{\pagegrid@cur<\pagegrid@col}{\loopwhile{\append@column@\pagegrid@cur\pagegrid@col}}{}%
   \box@column#1%
  \egroup
  \vskip\z@skip
 \egroup
}%
\def\append@column@#1#2{%
 \expandafter\box@column\csname col@\the#1\endcsname
 \hfil\vrule\@width\columnseprule\hfil
 \advance#1\@ne
 \@ifnum{#1<#2}%
}%
\def\box@column#1{%
 \ltxgrid@info@sw{\class@info{\string\box@column\string#1}}{}%
 \raise\topskip
 \hb@xt@\columnwidth\bgroup
  \dimen@\ht#1\@ifdim{\dimen@>\@colht}{\dimen@\@colht}{}%
  \count@\vbadness\vbadness\@M
  \dimen@ii\vfuzz\vfuzz\maxdimen
  \ltxgrid@info@sw{\saythe\@colht\saythe\dimen@}{}%
  \vtop to\dimen@\bgroup
   \hrule\@height\z@
   \unvbox#1%
   \raggedcolumn@skip
  \egroup
  \vfuzz\dimen@ii
  \vbadness\count@
  \hss
 \egroup
}%
\def\marry@baselines{%
 \begingroup
  \setbox\z@\lastbox
  \@ifvoid{\z@}{%
   \endgroup
  }{%
   \aftergroup\kern
   \aftergroup-%
   \expandafter\box\expandafter\z@\expandafter\endgroup\the\dp\z@\relax
  }%
 \vskip\marry@skip\relax
}%
\gdef\marry@skip{\z@skip}%
\def\set@marry@skip{%
 \begingroup
  \skip@\baselineskip\advance\skip@-\topskip
  \@ifdim{\skip@>\z@}{%
   \xdef\marry@skip{\the\skip@}%
  }{}%
 \endgroup
}%
\appdef\document@inithook{%
 \@ifxundefined\raggedcolumn@sw{\@booleanfalse\raggedcolumn@sw}{}%
}%
\def\raggedcolumn@skip{%
 \vskip\z@\raggedcolumn@sw{\@plus.0001fil\@minus.0001fil}{}\relax
}%
\def\@combinepage#1{%
 \ltxgrid@foot@info@sw{\class@info{\string\@combinepage\string#1}}{}%
 \@ifvoid\pagesofar{}{%
  \setbox\@outputbox\vbox{%
   \unvbox\pagesofar
   \marry@baselines
   \unvbox\@outputbox
  }%
 }%
 #1{%
  \@ifvoid\footsofar{}{%
   \show@box@size{Combining page footnotes}\footsofar
   \setbox\footins\box\footsofar
   \compose@footnotes
   \@combineinserts\@outputbox\footins
  }%
 }{%
 }%
}%
\def \@cflt{%
 \let \@elt \@comflelt
 \setbox\@tempboxa \vbox{}%
 \@toplist
 \setbox\@outputbox \vbox{%
  \boxmaxdepth \maxdepth
  \unvbox\@tempboxa\unskip
  \topfigrule\vskip \textfloatsep
  \unvbox\@outputbox
 }%
 \let\@elt\relax
 \xdef\@freelist{\@freelist\@toplist}%
 \global\let\@toplist\@empty
}%
\def \@cflb {%
 \let\@elt\@comflelt
 \setbox\@tempboxa \vbox{}%
 \@botlist
 \setbox\@outputbox \vbox{%
  \unvbox\@outputbox
  \vskip \textfloatsep\botfigrule
  \unvbox\@tempboxa\unskip
 }%
 \let\@elt\relax
 \xdef\@freelist{\@freelist\@botlist}%
 \global \let \@botlist\@empty
}%
\def\@combinedblfloats{%
 \@ifx@empty\@dbltoplist{}{%
  \setbox\@tempboxa\vbox{}%
  \let\@elt\@comdblflelt\@dbltoplist
  \let\@elt\relax\xdef\@freelist{\@freelist\@dbltoplist}%
  \global\let\@dbltoplist\@empty
  \setbox\@outputbox\vbox{%
   %\boxmaxdepth\maxdepth   %% probably not needed, CAR
   \unvbox\@tempboxa\unskip
   \@ifnum{\@dbltopnum>\m@ne}{\dblfigrule}{}%FIXME: how is \@dbltopnum maintained?
   \vskip\dbltextfloatsep
   \unvbox\@outputbox
  }%
 }%
}%
\def\set@column@hsize#1{%
 \pagegrid@col#1%
 \global\columnwidth\textwidth
 \global\advance\columnwidth\columnsep
 \global\divide\columnwidth\pagegrid@col
 \global\advance\columnwidth-\columnsep
 \global\hsize\columnwidth
 \global\linewidth\columnwidth
 \skip@\baselineskip\advance\skip@-\topskip
 \@ifnum{\pagegrid@col>\@ne}{\set@marry@skip}{}%
}%
\def\set@colht{%
 \set@adj@textheight\@colht
 \global\let\enlarge@colroom\@empty
 \set@colroom
}%
\def\set@adj@textheight#1{%
 \ltxgrid@info@sw{\class@info{\string\set@adj@textheight\string#1}\saythe\textheight}{}%
 #1\textheight
 \def\@elt{\adj@page#1}%
 \@booleantrue\firsttime@sw\@dbltoplist
 \let\@elt\relax
 \global#1#1\relax
 \ltxgrid@info@sw{\saythe#1}{}%
}%
\def\set@colroom{%
 \ltxgrid@info@sw{\class@info{\string\set@colroom}}{}%
 \set@adj@colht\@colroom
 \@if@empty\enlarge@colroom{}{%
  \global\advance\@colroom\enlarge@colroom\relax
  \ltxgrid@info@sw{\saythe\@colroom}{}%
 }%
 \@ifdim{\@colroom>\topskip}{}{%
  \ltxgrid@info{Not enough room: \string\@colroom=\the\@colroom; increasing to \the\topskip}%
  \@colroom\topskip
 }%
 \global\@colroom\@colroom
 \set@vsize
}%
\def\set@vsize{%
 \global\vsize\@colroom
 \ltxgrid@info@sw{\class@info{\string\set@vsize\string\vsize=\string\colroom}\saythe\vsize}{}%
}%
\def\set@adj@colht#1{%
 #1\@colht
 \ltxgrid@info@sw{\class@info{\string\set@adj@colht\string#1-\string\pagesofar}\saythe#1}{}%
 \@ifvoid\pagesofar{}{%
  \advance#1-\ht\pagesofar\advance#1-\dp\pagesofar
  \ltxgrid@info@sw{\class@info{\string\pagesofar}\saythe#1}{}%
 }%
 \def\@elt{\adj@column#1}%
 \@booleantrue\firsttime@sw\@toplist
 \@booleantrue\firsttime@sw\@botlist
 \let\@elt\relax
}%
\def\adj@column#1#2{%
 \advance#1-\ht#2%
 \advance#1-\firsttime@sw{\textfloatsep\@booleanfalse\firsttime@sw}{\floatsep}%
 \ltxgrid@info@sw{\class@info{\string\adj@column\string#1-\string#2}\saythe#1}{}%
}%
\def\adj@page#1#2{%
 \advance#1-\ht#2%
 \advance#1-\firsttime@sw{\dbltextfloatsep\@booleanfalse\firsttime@sw}{\dblfloatsep}%
 \ltxgrid@info@sw{\class@info{\string\adj@page\string#1-\string#2}\saythe#1}{}%
}%
\def\set@adj@box#1#2{%
 \@ifvoid#2{}{%
  \advance#1-\ht#2\advance#1-\dp#2%
  \@booleantrue\temp@sw
  \ltxgrid@foot@info@sw{\class@info{\string\set@adj@box\string#2}\saythe#1}{}%
 }%
}%
\appdef\@outputpage@tail{%
 \set@colht          % FIXME: needed?
 \@floatplacement    % FIXME: needed?
 \@dblfloatplacement % FIXME: needed?
}%
\begingroup
 \catcode`\1=\cat@letter
 \catcode`\2=\cat@letter
 \toks@{%
  \setbox\footins\box\footsofar
  \balance@two\col@1\@outputbox
  \global\setbox\col@1\box\col@1
  \global\setbox\@outputbox\box\@outputbox
  \combine@foot@inserts\footsofar\footins
 }%
 \aftergroup\def\aftergroup\balance@2\expandafter
\endgroup\expandafter{\the\toks@}%
\def\balance@two#1#2{%
 \ltxgrid@info@sw{\class@info{\string\balance@two\string#1\string#2}}{}%
 \outputdebug@sw{\trace@scroll{\showbox#1\showbox#2}}{}%
 \setbox\thr@@\copy\footsofar
 \setbox\@ne\vbox\bgroup
  \@ifvoid{#1}{}{%
   \recover@column#1\footsofar\column@recovered\footins@recovered
   \@ifvoid{#2}{}{\marry@baselines}%
  }%
  \@ifvoid{#2}{}{%
   \recover@column#2\footsofar\column@recovered\footins@recovered
  }%
 \egroup
 \outputdebug@sw{\trace@scroll{\showbox\@ne}}{}%
 \ltxgrid@foot@info@sw{\trace@scroll{\showbox\footsofar}}{}%
 \dimen@\ht\@ne\divide\dimen@\tw@
 \dimen@i\dimen@
 \vbadness\@M
 \vfuzz\maxdimen
 \splittopskip\topskip
 \loopwhile{%
  \setbox\z@\copy\@ne\setbox\tw@\vsplit\z@ to\dimen@
  \remove@depth\z@\remove@depth\tw@
  \dimen@ii\ht\tw@\advance\dimen@ii-\ht\z@
  \dimen@i=.5\dimen@i
  \ltxgrid@info@sw{\saythe\dimen@\saythe\dimen@i\saythe\dimen@ii}{}%
  \@ifdim{\dimen@ii<.5\p@}{%
   \@ifdim{\dimen@ii>-.5\p@}%
  }{%
   \false@sw
  }%
  {%
   \true@sw
  }{%
   \@ifdim{\dimen@i<.5\p@}%
  }%
  {%
   \false@sw
  }%
  {%
   \advance\dimen@\@ifdim{\dimen@ii<\z@}{}{-}\dimen@i
   \true@sw
  }%
 }%
 \ltxgrid@info@sw{\saythe\dimen@\saythe\dimen@i\saythe\dimen@ii}{}%
 \@ifdim{\ht\z@=\z@}{%
  \@ifdim{\ht\tw@=\z@}%
 }{%
  \true@sw
 }%
 {%
 }{%
  \ltxgrid@info{Unsatifactorily balanced columns: giving up}%
  \setbox\tw@\box#1%
  \setbox\z@ \box#2%
  \global\setbox\footsofar\box\thr@@
 }%
 \setbox\tw@\vbox{\unvbox\tw@\vskip\z@skip}%
 \setbox\z@ \vbox{\unvbox\z@ \vskip\z@skip}%
 \set@colht
 \dimen@\ht\z@\@ifdim{\dimen@<\ht\tw@}{\dimen@\ht\tw@}{}%
 \@ifdim{\dimen@>\@colroom}{\dimen@\@colroom}{}%
 \ltxgrid@info@sw{\saythe{\ht\z@}\saythe{\ht\tw@}\saythe\@colroom\saythe\dimen@}{}%
 \setbox#1\vbox to\dimen@{\unvbox\tw@\unskip\raggedcolumn@skip}%
 \setbox#2\vbox to\dimen@{\unvbox\z@ \unskip\raggedcolumn@skip}%
 \outputdebug@sw{\trace@scroll{\showbox#1\showbox#2}}{}%
}%
\def\remove@depth#1{%
  \setbox#1\vbox\bgroup
   \unvcopy#1%
   \setbox\z@\vbox\bgroup
    \unvbox#1%
    \setbox\z@\lastbox
    \aftergroup\kern\aftergroup-\expandafter
   \egroup
   \the\dp\z@\relax
  \egroup
}%
\def\recover@column#1#2#3#4{%
 \ltxgrid@info@sw{\class@info{\string\recover@column\string#1\string#2\string#3\string#4}}{}%
 \setbox#4\vbox{\unvcopy#1}%
 \ltxgrid@foot@info@sw{\trace@scroll{\showbox#4}}{}%
 \dimen@\ht#4%
 \ltxgrid@foot@info@sw{\saythe\dimen@}{}%
 \setbox#4\vbox\bgroup
  \unvbox#4\unskip
  \dimen@i\lastkern\unkern\advance\dimen@i\lastkern
  \@ifdim{\dimen@i=\z@}{%
   \dimen@i\lastkern\unkern
   \ltxgrid@foot@info@sw{\saythe\dimen@i}{}%
   \aftergroup\dimen@i
   \expandafter\egroup\the\dimen@i\relax
  }{%
   \egroup
  }%
 \@ifdim{\dimen@i<\z@}{%
  \advance\dimen@\dimen@i
  \ltxgrid@foot@info@sw{\saythe\dimen@i\saythe\dimen@}{}%
  \splittopskip\z@skip
  \global\setbox#3\vsplit#4 to\dimen@
  \global\setbox#4\vbox{\unvbox#4}%
  \ltxgrid@foot@info@sw{\trace@scroll{\showbox#1\showbox#2\showbox#3\showbox#4}}{}%
  \global\setbox#2\vbox\bgroup\unvbox#2\vskip\z@skip\unvbox#4\egroup
 }{%
  \setbox#3\box#4%
  \ltxgrid@foot@info@sw{\trace@scroll{\showbox#1\showbox#2\showbox#3\showbox#4}}{}%
 }%
 \unvbox#3%
 \loopwhile{\dimen@\lastskip\@ifdim{\dimen@>\z@}{\unskip\true@sw}{\false@sw}}%
}%
\def\recover@column@null#1#2#3#4{%
 \unvcopy#1%
}%
\prepdef\@begindocumenthook{%
 \open@column@one\@ne
 \set@colht
 \@floatplacement
 \@dblfloatplacement
}%
\def\longtable@longtable{%
 \par
 \ifx\multicols\@undefined\else\ifnum\col@number>\@ne\@twocolumntrue\fi\fi
 \if@twocolumn\LT@err{longtable not in 1-column mode}\@ehc\fi
 \begingroup
 \@ifnextchar[\LT@array{\LT@array[x]}%
}%
\def\longtable@new{%
 \par
  \@ifnextchar[\LT@array{\LT@array[x]}%
}%
\def\endlongtable@longtable{%
  \crcr
  \noalign{%
    \let\LT@entry\LT@entry@chop
    \xdef\LT@save@row{\LT@save@row}}%
  \LT@echunk
  \LT@start
  \unvbox\z@
  \LT@get@widths
  \if@filesw
    {\let\LT@entry\LT@entry@write\immediate\write\@auxout{%
      \gdef\expandafter\noexpand
        \csname LT@\romannumeral\c@LT@tables\endcsname
          {\LT@save@row}}}%
  \fi
  \ifx\LT@save@row\LT@@save@row
  \else
    \LT@warn{Column \@width s have changed\MessageBreak
             in table \thetable}%
    \LT@final@warn
  \fi
  \endgraf\penalty -\LT@end@pen
  \endgroup
  \global\@mparbottom\z@
  \pagegoal\vsize
  \endgraf\penalty\z@\addvspace\LTpost
  \ifvoid\footins\else\insert\footins{}\fi
}%
\def\endlongtable@new{%
  \crcr
  \noalign{%
   \let\LT@entry\LT@entry@chop
   \xdef\LT@save@row{\LT@save@row}%
  }%
  \LT@echunk
  \LT@start
  \unvbox\z@
  \LT@get@widths
  \@if@sw\if@filesw\fi{%
   {%
    \let\LT@entry\LT@entry@write
    \immediate\write\@auxout{%
     \gdef\expandafter\noexpand\csname LT@\romannumeral\c@LT@tables\endcsname
     {\LT@save@row}%
    }%
   }%
  }{}%
  \@ifx{\LT@save@row\LT@@save@row}{}{%
   \LT@warn{%
    Column \@width s have changed\MessageBreak in table \thetable
   }\LT@final@warn
  }%
  \endgraf
  \nobreak
  \box\@ifvoid\LT@lastfoot{\LT@foot}{\LT@lastfoot}%
 \global\@mparbottom\z@
 \endgraf
 \LT@post
}%
\def\LT@start@longtable{%
  \let\LT@start\endgraf
  \endgraf\penalty\z@\vskip\LTpre
  \dimen@\pagetotal
  \advance\dimen@ \ht\ifvoid\LT@firsthead\LT@head\else\LT@firsthead\fi
  \advance\dimen@ \dp\ifvoid\LT@firsthead\LT@head\else\LT@firsthead\fi
  \advance\dimen@ \ht\LT@foot
  \dimen@ii\vfuzz
  \vfuzz\maxdimen
    \setbox\tw@\copy\z@
    \setbox\tw@\vsplit\tw@ to \ht\@arstrutbox
    \setbox\tw@\vbox{\unvbox\tw@}%
  \vfuzz\dimen@ii
  \advance\dimen@ \ht
        \ifdim\ht\@arstrutbox>\ht\tw@\@arstrutbox\else\tw@\fi
  \advance\dimen@\dp
        \ifdim\dp\@arstrutbox>\dp\tw@\@arstrutbox\else\tw@\fi
  \advance\dimen@ -\pagegoal
  \ifdim \dimen@>\z@\vfil\break\fi
      \global\@colroom\@colht
  \ifvoid\LT@foot\else
    \advance\vsize-\ht\LT@foot
    \global\advance\@colroom-\ht\LT@foot
    \dimen@\pagegoal\advance\dimen@-\ht\LT@foot\pagegoal\dimen@
    \maxdepth\z@
  \fi
  \ifvoid\LT@firsthead\copy\LT@head\else\box\LT@firsthead\fi
\nobreak
  \output{\LT@output}%
}%
\def\LT@start@new{%
 \let\LT@start\endgraf
 \endgraf
 \markthr@@{}%
 \LT@pre
 \@ifvoid\LT@firsthead{\LT@top}{\box\LT@firsthead\nobreak}%
 \mark@envir{longtable}%
}%
\def\LT@end@hd@ft@longtable#1{%
 \LT@echunk
 \ifx\LT@start\endgraf
  \LT@err{Longtable head or foot not at start of table}{Increase LTchunksize}%
 \fi
 \setbox#1\box\z@
 \LT@get@widths\LT@bchunk
}%
\def\LT@end@hd@ft@new#1{%
 \LT@echunk
 \@ifx{\LT@start\endgraf}{%
  \LT@err{Longtable head or foot not at start of table}{Increase LTchunksize}%
 }%
 \global\setbox#1\box\z@
 \LT@get@widths
 \LT@bchunk
}%
\def\LT@array@longtable[#1]#2{%
  \refstepcounter{table}\stepcounter{LT@tables}%
  \if l#1%
    \LTleft\z@ \LTright\fill
  \else\if r#1%
    \LTleft\fill \LTright\z@
  \else\if c#1%
    \LTleft\fill \LTright\fill
  \fi\fi\fi
  \let\LT@mcol\multicolumn
  \let\LT@@tabarray\@tabarray
  \let\LT@@hl\hline
  \def\@tabarray{%
    \let\hline\LT@@hl
    \LT@@tabarray}%
  \let\\\LT@tabularcr\let\tabularnewline\\%
  \def\newpage{\noalign{\break}}%
  \def\pagebreak{\noalign{\ifnum`}=0\fi\@testopt{\LT@no@pgbk-}4}%
  \def\nopagebreak{\noalign{\ifnum`}=0\fi\@testopt\LT@no@pgbk4}%
  \let\hline\LT@hline \let\kill\LT@kill\let\caption\LT@caption
  \@tempdima\ht\strutbox
  \let\@endpbox\LT@endpbox
  \ifx\extrarowheight\@undefined
    \let\@acol\@tabacol
    \let\@classz\@tabclassz \let\@classiv\@tabclassiv
    \def\@startpbox{\vtop\LT@startpbox}%
    \let\@@startpbox\@startpbox
    \let\@@endpbox\@endpbox
    \let\LT@LL@FM@cr\@tabularcr
  \else
    \advance\@tempdima\extrarowheight
    \col@sep\tabcolsep
    \let\@startpbox\LT@startpbox\let\LT@LL@FM@cr\@arraycr
  \fi
  \setbox\@arstrutbox\hbox{\vrule
    \@height \arraystretch \@tempdima
    \@depth \arraystretch \dp \strutbox
    \@width \z@}%
  \let\@sharp##\let\protect\relax
   \begingroup
    \@mkpream{#2}%
    \xdef\LT@bchunk{%
       \global\advance\c@LT@chunks\@ne
       \global\LT@rows\z@\setbox\z@\vbox\bgroup
       \LT@setprevdepth
       \tabskip\LTleft \noexpand\halign to\hsize\bgroup
      \tabskip\z@ \@arstrut \@preamble \tabskip\LTright \cr}%
  \endgroup
  \expandafter\LT@nofcols\LT@bchunk&\LT@nofcols
  \LT@make@row
  \m@th\let\par\@empty
  \everycr{}\lineskip\z@\baselineskip\z@
  \LT@bchunk}%
\def\LT@LR@l{\LTleft\z@   \LTright\fill}%
\def\LT@LR@r{\LTleft\fill \LTright\z@  }%
\def\LT@LR@c{\LTleft\fill \LTright\fill}%
\def\LT@array@new[#1]#2{%
 \refstepcounter{table}\stepcounter{LT@tables}%
 \table@hook
 \LTleft\fill \LTright\fill
 \csname LT@LR@#1\endcsname
 \let\LT@mcol\multicolumn
 \let\LT@@hl\hline
 \prepdef\@tabarray{\let\hline\LT@@hl}%
 \let\\\LT@tabularcr
 \let\tabularnewline\\%
 \def\newpage{\noalign{\break}}%
 \def\pagebreak{\noalign{\ifnum`}=0\fi\@testopt{\LT@no@pgbk-}4}%
 \def\nopagebreak{\noalign{\ifnum`}=0\fi\@testopt\LT@no@pgbk4}%
 \let\hline\LT@hline
 \let\kill\LT@kill
 \let\caption\LT@caption
 \@tempdima\ht\strutbox
 \let\@endpbox\LT@endpbox
 \@ifxundefined\extrarowheight{%
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \def\@startpbox{\vtop\LT@startpbox}%
  \let\@@startpbox\@startpbox
  \let\@@endpbox\@endpbox
  \let\LT@LL@FM@cr\@tabularcr@LaTeX
  \let\@xtabularcr\@xtabularcr@LaTeX
 }{%
  \advance\@tempdima\extrarowheight
  \col@sep\tabcolsep
  \let\@startpbox\LT@startpbox
  \let\LT@LL@FM@cr\@arraycr@array
 }%
 \let\@acoll\@tabacoll
 \let\@acolr\@tabacolr
 \let\@acol\@tabacol
 \setbox\@arstrutbox\hbox{%
  \vrule
  \@height \arraystretch \@tempdima
  \@depth \arraystretch \dp \strutbox
  \@width \z@
 }%
 \let\@sharp##%
 \let\protect\relax
 \begingroup
  \@mkpream{#2}%
  \@mkpream@relax
  \edef\@preamble{\@preamble}%
  \prepdef\@preamble{%
   \global\advance\c@LT@chunks\@ne
   \global\LT@rows\z@
   \setbox\z@\vbox\bgroup
    \LT@setprevdepth
    \tabskip\LTleft
    \halign to\hsize\bgroup
     \tabskip\z@
     \@arstrut
  }%
  \appdef\@preamble{%
     \tabskip\LTright
     \cr
  }%
  \global\let\LT@bchunk\@preamble
 \endgroup
 \expandafter\LT@nofcols\LT@bchunk&\LT@nofcols
 \LT@make@row
 \m@th
 \let\par\@empty
 \everycr{}%
 \lineskip\z@
 \baselineskip\z@
 \LT@bchunk
}%
\appdef\table@hook{}%
\def\switch@longtable{%
 \@ifpackageloaded{longtable}{%
  \@ifx{\longtable\longtable@longtable}{%
   \@ifx{\endlongtable\endlongtable@longtable}{%
    \@ifx{\LT@start\LT@start@longtable}{%
     \@ifx{\LT@end@hd@ft\LT@end@hd@ft@longtable}{%
      \@ifx{\LT@array\LT@array@longtable}{%
       \true@sw
      }{\false@sw}%
     }{\false@sw}%
    }{\false@sw}%
   }{\false@sw}%
  }{\false@sw}%
  {%
   \class@info{Patching longtable package}%
  }{%
   \class@info{Patching unrecognized longtable package. (Proceeding with fingers crossed)}%
  }%
  \let\longtable\longtable@new
  \let\endlongtable\endlongtable@new
  \let\LT@start\LT@start@new
  \let\LT@end@hd@ft\LT@end@hd@ft@new
  \let\LT@array\LT@array@new
  \newenvironment{longtable*}{%
   \onecolumngrid@push
   \longtable
  }{%
   \endlongtable
   \onecolumngrid@pop
  }%
 }{}%
}%
\def\LT@pre{\penalty\z@\vskip\LTpre}%
\def\LT@bot{\nobreak\copy\LT@foot\vfil}%
\def\LT@top{\copy\LT@head\nobreak}%
\def\LT@post{\penalty\z@\addvspace\LTpost\mark@envir{\curr@envir}}%
\def\LT@adj{%
 \setbox\z@\vbox{\null}\dimen@-\ht\z@
 \setbox\z@\vbox{\unvbox\z@\LT@bot}\advance\dimen@\ht\z@
 \global\advance\vsize-\dimen@
}%
\def\output@init@longtable{\LT@adj}%
\def\output@prep@longtable{\setbox\@cclv\vbox{\unvbox\@cclv\LT@bot}}%
\def\output@post@longtable{\LT@top}%
\let\output@init@theindex\@empty
\let\output@prep@theindex\@empty
\def\output@post@theindex{%
 \@ifodd\c@page{}{%
  \@ifnum{\pagegrid@cur=\@ne}{%
  }%
 }%
}%
\def\check@aux{\do@output@MVL{\do@check@aux}}%
\def\check@deferlist@stuck#1{%
 \@ifx{\@deferlist@postshipout\@empty}{}{%
  \@ifx{\@deferlist@postshipout\@deferlist}{%
   \@fltstk
   \clearpage@sw{%
    \ltxgrid@warn{Deferred float stuck during \string\clearpage\space processing}%
   }{%
    \force@deferlist@stuck#1%
   }%
  }{%
  }%
  \global\let\@deferlist@postshipout\@empty
 }%
}%
\def\@fltstk{%
 \@latex@warning{A float is stuck (cannot be placed without \string\clearpage)}%
}%
\appdef\@outputpage@tail{%
 \global\let\@deferlist@postshipout\@deferlist
}%
\def\@next#1#2{%
 \@ifx{#2\@empty}{\false@sw}{%
  \expandafter\@xnext#2\@@#1#2%
  \true@sw
 }%
}%
\def\@xnext\@elt#1#2\@@#3#4{%
 \def#3{#1}%
 \gdef#4{#2}%
 \def\@tempa{#4}\def\@tempb{\@freelist}%
 \@ifx{\@tempa\@tempb}{%
  \@ifx{#4\@empty}{%
   \force@deferlist@empty%{Float register pool exhausted}%
  }{}%
 }{}%
}%
\def\force@deferlist@stuck#1{%
 \force@deferlist@sw{%
  \@booleantrue\clearpage@sw
  \@booleantrue\forcefloats@sw
  #1%
 }{%
 }%
}%
\def\force@deferlist@empty{%
 \force@deferlist@sw{%
  \penalty-\pagebreak@pen
  \protect@penalty\do@forcecolumn@pen
 }{%
 }%
}%
\@booleanfalse\force@deferlist@sw
\mathchardef\do@forcecolumn@pen=10009
\@namedef{output@-\the\do@forcecolumn@pen}{\do@forcecolumn}%
\def\do@forcecolumn{%
 \@booleantrue\clearpage@sw
 \@booleantrue\forcefloats@sw
 \do@startcolumn
}%
\def\enlargethispage{%
 \@ifstar{%
  \@enlargethispage{}%
 }{%
  \@enlargethispage{}%
 }%
}%
\def\@enlargethispage#1#2{%
 \begingroup
  \dimen@#2\relax
  \edef\@tempa{#1}%
  \edef\@tempa{\noexpand\@@enlargethispage{\@tempa}{\the\dimen@}}%
  \expandafter\do@output@MVL\expandafter{\@tempa}%
 \endgroup
}%
\def\@@enlargethispage#1#2{%
 \def\@tempa{one}%
 \@ifx{\thepagegrid\@tempa}{%
  \true@sw
 }{%
  \def\@tempa{mlt}%
  \@ifx{\thepagegrid\@tempa}{%
   \@ifnum{\pagegrid@cur=\@ne}{%
    \gdef\enlarge@colroom{#2}%
    \true@sw
   }{%
    \ltxgrid@warn{Too late to enlarge this page; move the command to the first column.}%
    \false@sw
   }%
  }{%
   \ltxgrid@warn{Unable to enlarge a page of this kind.}%
   \false@sw
  }%
 }%
 {%
  \class@info{Enlarging page \thepage\space by #2}%
  \global\advance\@colroom#2\relax
  \set@vsize
 }{%
 }%
}%
\let\enlarge@colroom\@empty
\let\@kludgeins\@undefined
\@booleantrue\textheight@sw
\prepdef\@outputpage@head{%
 \textheight@sw{%
  \count@\vbadness\vbadness\@M
  \dimen@\vfuzz\vfuzz\maxdimen
  \setbox\@outputbox\vbox to\textheight{\unvbox\@outputbox}%
  \vfuzz\dimen@
  \vbadness\count@
 }{}%
}%
\appdef\@outputpage@head{%
 \@ifx{\LS@rot\@undefined}{}{\LS@rot}%
}%
\def\ltxgrid@info{%
 \ltxgrid@info@sw{\class@info}{\@gobble}%
}%
\@booleanfalse\ltxgrid@info@sw
\def\ltxgrid@warn{%
 \ltxgrid@warn@sw{\class@warn}{\@gobble}%
}%
\@booleantrue\ltxgrid@warn@sw
\@booleanfalse\ltxgrid@foot@info@sw
\def\def@next@handler#1#2#3{%
 \advance#1\@ne\mathchardef#2\the#1%
 \expandafter\def\csname output@-\the#1\endcsname{#3}%
}%
\def\def@line@handler#1#2{%
 \begingroup
  \@tempcnta\int@parpenalty
  \advance\@tempcnta-#1%
  \aftergroup\def
  \expandafter\aftergroup\csname output@-\the\@tempcnta\endcsname
 \endgroup{#2}%
}%
\mathchardef\int@parpenalty11012
\def@line@handler\z@{\@handle@line@ltx{}{}{}}%
\def@line@handler\@ne{\@handle@line@ltx{}{}{\brokenpenalty@ltx}}%
\def@line@handler\tw@{\@handle@line@ltx{}{\clubpenalty@ltx}{}}%
\def@line@handler\thr@@{\@handle@line@ltx{\clubpenalty@ltx}{}{\brokenpenalty@ltx}}%
\def@line@handler\f@ur{\@handle@line@ltx{\widowpenalty@ltx}{}{}}%
\def@line@handler{5}{\@handle@line@ltx{\widowpenalty@ltx}{}{\brokenpenalty@ltx}}%
\def@line@handler{6}{\@handle@line@ltx{\widowpenalty@ltx}{\clubpenalty@ltx}{}}%
\def@line@handler{7}{\@handle@line@ltx{\widowpenalty@ltx}{\clubpenalty@ltx}{\brokenpenalty@ltx}}%
\def@line@handler{8}{\@handle@line@ltx{\displaywidowpenalty@ltx}{}{}}%
\def@line@handler{9}{\@handle@line@ltx{\displaywidowpenalty@ltx}{}{\brokenpenalty@ltx}}%
\def@line@handler{10}{\@handle@line@ltx{\displaywidowpenalty@ltx}{\clubpenalty@ltx}{}}%
\def@line@handler{11}{\@handle@line@ltx{\displaywidowpenalty@ltx}{\clubpenalty@ltx}{\brokenpenalty@ltx}}%
\def\@handle@line@ltx#1#2#3{%
 \@@handle@line@ltx
 \@tempcnta\lastpenalty
 \@tempcntb\interlinepenalty@ltx\relax
 \@if@empty{#1}{}{\advance\@tempcntb#1\relax}%
 \@if@empty{#2}{}{\advance\@tempcntb#2\relax}%
 \@if@empty{#3}{}{\advance\@tempcntb#3\relax}%
 \penalty\@ifnum{\@tempcnta<\@tempcntb}{\@tempcntb}{\@tempcnta}%
}%
\let\@@handle@line@ltx\@empty
\@tempcnta\int@parpenalty
\def@next@handler\@tempcnta\int@postparpenalty{\reset@queues@ltx\handle@par@ltx}%
\def@next@handler\@tempcnta\int@vadjustpenalty{\handle@vadjust@ltx}%
\def@next@handler\@tempcnta\int@whatsitpenalty{\handle@whatsit@ltx}%
\def@next@handler\@tempcnta\int@predisplaypenalty{\reset@queues@ltx\@handle@display@ltx{\predisplaypenalty@ltx}}%
\def@next@handler\@tempcnta\int@interdisplaylinepenalty{\@handle@display@ltx{\interdisplaylinepenalty@ltx}}%
\def@next@handler\@tempcnta\int@postdisplaypenalty{\@handle@display@ltx{\postdisplaypenalty@ltx}}%
\def\@handle@display@ltx#1{%
 \@@handle@display@ltx
 \@tempcnta\lastpenalty
 \@tempcntb#1%
 \penalty\@ifnum{\@tempcnta<\@tempcntb}{\@tempcntb}{\@tempcnta}%
}%
\let\@@handle@display@ltx\@empty
\def\handle@par@ltx{}%
\def\set@linepenalties{%
 \expandafter\def\expandafter\interlinepenalty@ltx\expandafter{\the\interlinepenalty}%
 \interlinepenalty-\int@parpenalty
 \expandafter\def\expandafter\brokenpenalty@ltx\expandafter{\the\brokenpenalty}%
 \brokenpenalty\@ne
 \expandafter\def\expandafter\clubpenalty@ltx\expandafter{\the\clubpenalty}%
 \clubpenalty\tw@
 \expandafter\def\expandafter\widowpenalty@ltx\expandafter{\the\widowpenalty}%
 \widowpenalty\f@ur
 \expandafter\def\expandafter\displaywidowpenalty@ltx\expandafter{\the\displaywidowpenalty}%
 \displaywidowpenalty8\relax
}%
\def\restore@linepenalties{%
 \interlinepenalty\interlinepenalty@ltx
 \brokenpenalty\brokenpenalty@ltx
 \clubpenalty\clubpenalty@ltx
 \widowpenalty\widowpenalty@ltx
 \displaywidowpenalty\displaywidowpenalty@ltx
 \relax
}%
\def\set@displaypenalties#1{%
 \expandafter\def\expandafter\predisplaypenalty@ltx\expandafter{\the\predisplaypenalty}%
 \expandafter\def\expandafter\interdisplaylinepenalty@ltx\expandafter{\the\interdisplaylinepenalty}%
 \expandafter\def\expandafter\postdisplaypenalty@ltx\expandafter{\the\postdisplaypenalty}%
 \@ifhmode{\predisplaypenalty-\int@predisplaypenalty\relax}{}%
 #1{\interdisplaylinepenalty-\int@interdisplaylinepenalty\relax}{}%
 #1{\postdisplaypenalty-\int@postdisplaypenalty\relax}{}%
}%
\def\enqueue@whatsit@ltx#1{%
 \gappdef\g@whatsit@queue{{#1}}%
 \vadjust{\penalty-\int@whatsitpenalty}%
}%
\def\handle@whatsit@ltx{%
 \unvbox\@cclv
 \g@pop@ltx\g@whatsit@queue\@tempa
 \expandafter\do@whatsit\expandafter{\@tempa}%
}%
\def\do@whatsit#1{}%
\def\g@pop@ltx#1#2{%
 \expandafter\@g@pop@ltx#1{}{}\@@#1#2%
}%
\def\@g@pop@ltx#1#2\@@#3#4{%
 \gdef#3{#2}%
 \def#4{#1}%
}%
\let\vspace@ltx\vspace
\let\pagebreak@ltx\pagebreak
\let\nopagebreak@ltx\nopagebreak
\let\endline@ltx\\
\let\@arrayparboxrestore@ltx\@arrayparboxrestore
\def\@vspace@org #1{%
  \ifvmode
    \vskip #1
    \vskip\z@skip
   \else
     \@bsphack
     \vadjust{\@restorepar
              \vskip #1
              \vskip\z@skip
              }%
     \@esphack
   \fi
}%
\def\@vspace@ltx#1{%
 \@ifvmode{%
  \vskip#1\vskip\z@skip
 }{%
  \@bsphack
  \ex@vadjust@ltx{%
   \@restorepar
   \nobreak
   \vskip#1\vskip\z@skip
  }%
  \@esphack
 }%
}%
\def\@vspacer@org#1{%
  \ifvmode
    \dimen@\prevdepth
    \hrule \@height\z@
    \nobreak
    \vskip #1
    \vskip\z@skip
    \prevdepth\dimen@
  \else
    \@bsphack
    \vadjust{\@restorepar
             \hrule \@height\z@
             \nobreak
             \vskip #1
             \vskip\z@skip}%
    \@esphack
\fi
}%
\def\@vspacer@ltx#1{%
 \@ifvmode{%
  \dimen@\prevdepth
  \hrule\@height\z@
  \nobreak
  \vskip#1\vskip\z@skip
  \prevdepth\dimen@
 }{%
  \@bsphack
  \ex@vadjust@ltx{%
   \@restorepar
   \hrule\@height\z@
   \nobreak
   \vskip#1\vskip\z@skip
  }%
  \@esphack
 }%
}%
\def\@no@pgbk@org #1[#2]{%
  \ifvmode
    \penalty #1\@getpen{#2}%
  \else
    \@bsphack
    \vadjust{\penalty #1\@getpen{#2}}%
    \@esphack
  \fi
}%
\def\@no@pgbk@ltx#1[#2]{%
 \@ifvmode{%
  \penalty#1\@getpen{#2}%
 }{%
  \@bsphack
  \ex@vadjust@ltx{%
   \penalty#1\@getpen{#2}%
  }%
  \@esphack
 }%
}%
\long\def\end@line@org{%
 \let\reserved@e\relax
 \let\reserved@f\relax
 \@ifstar{%
  \let\reserved@e\vadjust
  \let\reserved@f\nobreak
  \@xnewline
 }%
 \@xnewline
}%
\long\def\end@line@ltx{%
 \let\reserved@e\relax
 \let\reserved@f\relax
 \@ifstar{%
  \let\reserved@e\ex@vadjust@ltx
  \let\reserved@f\nobreak
  \@xnewline
 }{%
  \@xnewline
 }%
}%
\def\@newline@org[#1]{%
 \let\reserved@e\vadjust
 \@gnewline{\vskip#1}%
}%
\def\@newline@ltx[#1]{%
 \let\reserved@e\ex@vadjust@ltx
 \@gnewline{\vskip#1}%
}%
 \@ifx{\@vspace\@vspace@org}{%
  \@ifx{\@vspacer\@vspacer@org}{%
   \@ifx{\@no@pgbk\@no@pgbk@org}{%
    \@ifx{\@newline\@newline@org}{%
     \expandafter\@ifx\expandafter{\csname\expandafter\@gobble\string\\ \endcsname\end@line@org}{%
      \true@sw
     }{\false@sw}%
    }{\false@sw}%
   }{\false@sw}%
  }{\false@sw}%
 }{\false@sw}%
 {%
  \class@info{Overriding \string\@vspace, \string\@vspacer, \string\@no@pgbk, \string\@newline, and \string\\ }%
  \let\@normalcr\end@line@ltx
  \expandafter\let\csname\expandafter\@gobble\string\\ \endcsname\@normalcr
  \let\@newline\@newline@ltx
  \let\@vspace\@vspace@ltx
  \let\@vspacer\@vspacer@ltx
  \let\@no@pgbk\@no@pgbk@ltx
 }{%
  \class@warn{%
   Failed to recognize \string\@vspace, \string\@vspacer, \string\@no@pgbk, \string\@newline, and \string\\;
   no patches applied. Please get a more up-to-date class,
  }%
 }%
\let\ex@vadjust@ltx\vadjust
\def\enqueue@vadjust@ltx#1{%
 \gappdef\g@vadjust@queue{{#1}}%
 \vadjust{\penalty-\int@vadjustpenalty}%
}%
\def\handle@vadjust@ltx{%
 \unvbox\@cclv
 \g@pop@ltx\g@vadjust@queue\@tempa
 \expandafter\gappdef\expandafter\g@vadjust@line\expandafter{\@tempa}%
}%
\let\g@vadjust@line\@empty
\def\reset@queues@ltx{%
 \global\let\g@whatsit@queue\@empty
 \global\let\g@vadjust@queue\@empty
}%
\newcommand\linenomathWithnumbers@LN{%
  \ifLineNumbers
    \ifnum\interlinepenalty>-\linenopenaltypar
      \global\holdinginserts\thr@@
      \advance\interlinepenalty \linenopenalty
     \ifhmode
      \advance\predisplaypenalty \linenopenalty
     \fi
      \advance\postdisplaypenalty \linenopenalty
      \advance\interdisplaylinepenalty \linenopenalty
    \fi
  \fi
  \ignorespaces
}%
\newcommand\linenomathNonumbers@LN{%
  \ifLineNumbers
    \ifnum\interlinepenalty>-\linenopenaltypar
      \global\holdinginserts\thr@@
      \advance\interlinepenalty \linenopenalty
     \ifhmode
      \advance\predisplaypenalty \linenopenalty
     \fi
    \fi
  \fi
  \ignorespaces
}%
\def\endlinenomath@LN{%
  \ifLineNumbers
   \global\holdinginserts\@LN@outer@holdins
  \fi
  \global\@ignoretrue
}
\def\linenumberpar@LN{%
  \ifvmode \@@@par \else
    \ifinner \@@@par \else
      \xdef\@LN@outer@holdins{\the\holdinginserts}%
      \advance \interlinepenalty \linenopenalty
      \linenoprevgraf \prevgraf
      \global \holdinginserts \thr@@
      \@@@par
      \ifnum\prevgraf>\linenoprevgraf
        \penalty-\linenopenaltypar
      \fi
      \@LN@parpgbrk
      \global\holdinginserts\@LN@outer@holdins
      \advance\interlinepenalty -\linenopenalty
    \fi
  \fi
}%
\appdef\class@documenthook{%
 \@ifpackageloaded{lineno}{%
  \@ifx{\linenomathWithnumbers\linenomathWithnumbers@LN}{%
   \@ifx{\linenomathNonumbers\linenomathNonumbers@LN}{%
    \@ifx{\endlinenomath\endlinenomath@LN}{%
     \@ifx{\linenumberpar\linenumberpar@LN}{%
      \true@sw
     }{\false@sw}%
    }{\false@sw}%
   }{\false@sw}%
  }{\false@sw}%
  {%
   \class@info{Overriding lineo.sty, restoring output routine,}%
   \let\linenumberpar\linenumberpar@ltx
   \let\endlinenomath\endlinenomath@ltx
   \expandafter\let\csname endlinenomath*\endcsname\endlinenomath@ltx
   \let\linenomathWithnumbers\linenomathWithnumbers@ltx
   \let\linenomathNonumbers\linenomathNonumbers@ltx
   \let\ex@vadjust@ltx\ex@vadjust@line
   \let\@LN@postlabel\enqueue@whatsit@ltx
   \let\do@whatsit\write@linelabel
   \let\handle@par@ltx\handle@par@LN
   \let\@@handle@line@ltx\Make@LineNo@ltx
   \let\@@handle@display@ltx\Make@LineNo@ltx
   \output@latex{\natural@output}%
   \let\vspace\vspace@ltx
   \let\pagebreak\pagebreak@ltx
   \let\nopagebreak\nopagebreak@ltx
   \let\@arrayparboxrestore\@arrayparboxrestore@ltx
   \let\\\endline@ltx
   \appdef\set@footnotefont{%
    \let\par\@@@par
    \let\@@par\@@@par
   }%
   \@if@sw\ifLineNumbers\fi{%
    \class@info{Reinvoke \string\linenumbers}%
    \let\@@par\linenumberpar
    \@ifx{\@par\linenumberpar@LN}{\let\@par\linenumberpar}{}%
    \@ifx{\par\linenumberpar@LN}{\let\par\linenumberpar}{}%
   }{%
    \class@info{Line numbering not turned on yet}%
   }%
  }{%
   \class@warn{Failed to recognize lineno.sty procedures; no patches applied. Please get a more up-to-date class.}%
  }%
 }{%
 }%
}%
\def\linenumberpar@ltx{\@ifvmode{\@@@par}{\@linenumberpar}}%
\def\@linenumberpar{%
   \linenoprevgraf\prevgraf
   \set@linepenalties
   \@@@par
   \@ifnum{\prevgraf>\linenoprevgraf}{
    \penalty-\int@postparpenalty
   }{}%
   \@LN@parpgbrk
   \restore@linepenalties
}%
\newcommand\linenomathWithnumbers@ltx{\@linenomathnumbers@ltx\true@sw}%
\newcommand\linenomathNonumbers@ltx{\@linenomathnumbers@ltx\false@sw}%
\def\@linenomathnumbers@ltx#1{%
 \@if@sw\ifLineNumbers\fi{%
   \set@linepenalties
   \set@displaypenalties#1%
 }{}%
 \ignorespaces
}%
\def\endlinenomath@ltx{%
 \global\@ignoretrue
}%
\def\handle@par@LN{%
 \Make@LineNo@ltx
 \@tempcnta\lastpenalty
 \@ifnum{\@tempcnta=\z@}{}{%
  \expandafter\gdef
  \expandafter\@LN@parpgbrk
  \expandafter{%
   \expandafter\penalty
               \the\@tempcnta
   \global\let\@LN@parpgbrk\@LN@screenoff@pen
  }%
 }%
}%
\def\Make@LineNo@ltx{%
 \@LN@maybe@normalLineNumber
 \boxmaxdepth\maxdimen\setbox\z@\vbox{\unvbox\@cclv}%
 \@tempdima\dp\z@
 \unvbox\z@
 \sbox\@tempboxa{\hb@xt@\z@{\makeLineNumber}}%
 \ht\@tempboxa\z@
 \@LN@depthbox
 \stepLineNumber
 \g@vadjust@line
 \global\let\g@vadjust@line\@empty
}%
\def\write@linelabel#1{%
 \protected@write\@auxout{}{%
  \string\newlabel{#1}{{\theLineNumber}{\thepage}{}{}{}}%
 }%
}%
\def\ex@vadjust@line{%
 \@if@sw\ifLineNumbers\fi{\enqueue@vadjust@ltx}{\vadjust}%
}%
\let\setup@hook\@empty
\appdef\document@inithook{%
 \switch@longtable
 \let\LT@makecaption\LT@makecaption@rtx
}%
\def\LT@makecaption@rtx#1#2#3{%
  \LT@mcol\LT@cols c{%
    \hbox to\z@{%
     \hss
     \parbox[t]\LTcapwidth{%
      \sbox\@tempboxa{#1{#2: }#3\unskip\nobreak\vrule\@width\z@\@height\z@\@depth .5\baselineskip}%
      \ifdim\wd\@tempboxa>\hsize
        #1{#2: }#3\unskip\nobreak\vrule\@width\z@\@height\z@\@depth .5\baselineskip
      \else
        \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
      \fi
      \endgraf
     }%
     \hss
    }%
  }%
}%
\def\protectdef@boldmath{%
 \expandafter\@ifnotrelax\csname boldmath \endcsname{}{%
  \class@info{Robustifying \string\LaTeX's \string\boldmath\space command}%
  \expandafter\let\csname boldmath \endcsname\boldmath
  \expandafter\def\expandafter\boldmath\expandafter{\expandafter\protect\csname boldmath \endcsname}%
 }%
}%
\appdef\document@inithook{%
 \protectdef@boldmath
}%
\DeclareOption{checkin}{%
  \@booleantrue\dateinRH@sw
  \@booleantrue\preprintsty@sw
  \def\@pointsize{12}%
  \@booleantrue\showPACS@sw
  \@booleantrue\showKEYS@sw
  \@booleantrue\showFEATS@sw
  \def\fp@proc@h{\allow@breaking@tables}%
  \def\fp@proc@H{\allow@breaking@tables}%
}%
\@booleanfalse\dateinRH@sw
\def\checkindate{\dateinRH@sw{{\tiny(\today)}}{}}%
\def\allow@breaking@tables{%
 \def\array@default{v}% tabular can break over pages
 \@booleanfalse\floats@sw % table can break over pages
}%
\DeclareOption{preprint}{%
 \@booleantrue\preprintsty@sw
 \ExecuteOptions{12pt}%
}%
\DeclareOption{reprint}{%
 \@booleanfalse\preprintsty@sw
 \@booleantrue\twocolumn@sw
 \ExecuteOptions{10pt}%
}%
\DeclareOption{manuscript}{%
 \class@warn{Document class option manuscript is obsolete; use preprint instead}%
 \ExecuteOptions{preprint}%
}%
\@booleanfalse\preprintsty@sw
\DeclareOption{showpacs}{%
  \@booleantrue\showPACS@sw
}%
\DeclareOption{noshowpacs}{%
  \@booleanfalse\showPACS@sw
}%
\DeclareOption{showkeys}{%
  \@booleantrue\showKEYS@sw
}%
\DeclareOption{noshowkeys}{%
  \@booleanfalse\showKEYS@sw
}%
\DeclareOption{showfeats}{%
  \@booleantrue\showFEATS@sw
}%
\DeclareOption{noshowfeats}{%
  \@booleanfalse\showFEATS@sw
}%
\@booleanfalse\showPACS@sw
\@booleanfalse\showKEYS@sw
\@booleanfalse\showFEATS@sw
\DeclareOption{balancelastpage}{%
  \@booleantrue\balancelastpage@sw
}%
\DeclareOption{nobalancelastpage}{%
  \@booleanfalse\balancelastpage@sw
}%
\@booleantrue\balancelastpage@sw
\DeclareOption{nopreprintnumbers}{\@booleanfalse\preprint@sw}%
\DeclareOption{preprintnumbers}{\@booleantrue\preprint@sw}%
\appdef\setup@hook{%
 \@ifxundefined\preprint@sw{\let\preprint@sw\preprintsty@sw}{}%
}%
\DeclareOption{hyperref}{%
 \class@warn{Class option "hyperref" is no longer supported.^^JEmploy \string\usepackage{hyperref} instead}%
}%
\DeclareOption{10pt}{\def\@pointsize{10}}%
\DeclareOption{11pt}{\def\@pointsize{11}}%
\DeclareOption{12pt}{\def\@pointsize{12}}%
\let\@pointsize\@undefined
\DeclareOption{a4paper}{%
    \setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}%
}%
\DeclareOption{a5paper}{%
    \setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}%
}%
\DeclareOption{b5paper}{%
    \setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}%
}%
\DeclareOption{letterpaper}{%
    \setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}%
}%
\DeclareOption{legalpaper}{%
    \setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}%
}%
\DeclareOption{executivepaper}{%
    \setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}%
}%
\DeclareOption{landscape}{%
    \setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}%
}%
\ds@letterpaper
\DeclareOption{bibnotes}{\let\frontmatter@footnote@produce\frontmatter@footnote@produce@endnote}%
\DeclareOption{nobibnotes}{\let\frontmatter@footnote@produce\frontmatter@footnote@produce@footnote}%
\let\frontmatter@footnote@produce\frontmatter@footnote@produce@footnote
\appdef\class@enddocumenthook{\auto@bib}%
\DeclareOption{footinbib}{\@booleantrue\footinbib@sw}
\DeclareOption{nofootinbib}{\@booleanfalse\footinbib@sw}
\@booleanfalse\footinbib@sw
\DeclareOption{altaffilletter}{\@booleantrue\altaffilletter@sw}%
\DeclareOption{altaffilsymbol}{\@booleanfalse\altaffilletter@sw}%
\@booleanfalse\altaffilletter@sw
\DeclareOption{superbib}{%
 \let\place@bibnumber\place@bibnumber@sup
}%
\def\place@bibnumber{\NATx@bibnumfmt}%
\def\place@bibnumber@sup#1{\textsuperscript{#1}}%
\def\place@bibnumber@inl#1{[#1]}%
\DeclareOption{citeautoscript}{\@booleantrue\citeautoscript@sw}%
\@booleanfalse\citeautoscript@sw
\DeclareOption{longbibliography}{\@booleantrue\longbibliography@sw}%
\DeclareOption{nolongbibliography}{\@booleanfalse\longbibliography@sw}%
\@booleanfalse\longbibliography@sw
\DeclareOption{eprint}{\@booleantrue\eprint@enable@sw}%
\DeclareOption{noeprint}{\@booleanfalse\eprint@enable@sw}%
\@booleantrue\eprint@enable@sw
\@booleanfalse\twoside@sw
\appdef\document@inithook{%
 \twoside@sw{\@twosidetrue}{\@twosidefalse}%
}%
\DeclareOption{twoside}{\@booleantrue \twoside@sw\@mparswitchfalse}%
\DeclareOption{oneside}{\@booleanfalse\twoside@sw\@mparswitchtrue}%
\DeclareOption{onecolumn}{\@booleanfalse\twocolumn@sw}%
\DeclareOption{twocolumn}{\@booleantrue \twocolumn@sw}%
\@booleanfalse\twocolumn@sw
\def\select@column@grid{%
 \twocolumn@sw{%
  \twocolumn@grid@setup
  \open@twocolumn
 }{%
  \onecolumn@grid@setup
 }%
}%
\appdef\class@documenthook{%
 \select@column@grid
}%
\appdef\setup@hook{%
 \let\clearpage@ltx\clearpage
 \prepdef\clear@document{\let\clearpage\clearpage@ltx\let\clear@document\@empty\close@column}%
 \appdef\class@documenthook{%
  \appdef\class@enddocumenthook{%
   \let\clearpage\clear@document
  }%
 }%
}%
\DeclareOption{author-year}{\@booleantrue\authoryear@sw}%
\DeclareOption{numerical}{\@booleanfalse\authoryear@sw}%
\@booleanfalse\authoryear@sw
\DeclareOption{galley}{%
  \ExecuteOptions{onecolumn}%
  \@booleantrue\galley@sw
  \@booleanfalse\preprintsty@sw
  \appdef\setup@hook{%
    \advance\textwidth-\columnsep
    \textwidth.5\textwidth
  }%
}%
\@booleanfalse\galley@sw
\DeclareOption{raggedbottom}{\@booleantrue\raggedcolumn@sw}
\DeclareOption{flushbottom}{\@booleanfalse\raggedcolumn@sw}
\@booleanfalse\raggedcolumn@sw
\appdef\setup@hook{%
 \raggedcolumn@sw{\raggedbottom}{\flushbottom}%
}%
\DeclareOption{tightenlines}{\@booleantrue\tightenlines@sw}
\@booleanfalse\tightenlines@sw
\@booleanfalse\lengthcheck@sw
\DeclareOption{lengthcheck}{%
 \@booleantrue\lengthcheck@sw
 \ExecuteOptions{reprint}%
}%
\appdef\setup@hook{%
 \lengthcheck@sw{\@booleantrue\tally@box@size@sw}{}%
}%
\appdef\setup@hook{%
 \draft@sw{\overfullrule 5\p@}{\overfullrule\z@}%
}%
\DeclareOption{draft}{\@booleantrue\draft@sw}%
\DeclareOption{final}{\@booleanfalse\draft@sw}%
\@booleanfalse\draft@sw
\appdef\setup@hook{%
 \eqsecnum@sw{%
  \@addtoreset{equation}{section}%
  \def\theequation@prefix{\arabic{section}.}%
 }{}%
}%
\DeclareOption{eqsecnum}{\@booleantrue\eqsecnum@sw}%
\@booleanfalse\eqsecnum@sw
\appdef\setup@hook{%
 \setup@secnums
}%
\DeclareOption{secnumarabic}{%
 \def\setup@secnums{\secnums@arabic}%
}%
\def\setup@secnums{\secnums@rtx}%
\DeclareOption{fleqn}{%
 \input{fleqn.clo}%
}%
\DeclareOption{floats}{\@booleantrue\floats@sw\@booleanfalse\floatp@sw}
\DeclareOption{endfloats}{\@booleanfalse\floats@sw\@booleanfalse\floatp@sw}
\DeclareOption{endfloats*}{\@booleanfalse\floats@sw\@booleantrue\floatp@sw}
\@booleantrue\floats@sw
\@booleantrue\floatp@sw
\DeclareOption{titlepage}{\@booleantrue\titlepage@sw}
\DeclareOption{notitlepage}{\@booleanfalse\titlepage@sw}
\@booleanfalse\titlepage@sw
\def\change@society#1{%
 \def\@tempa{#1}%
 \@ifxundefined\@society{%
  \class@info{Selecting society \@tempa}%
  \let\@society\@tempa
 }{%
  \@ifx{\@tempa\@society}{}{%
   \class@warn{Conflicting society \@tempa<>\@society; not selected}%
  }%
 }%
}%
\def\change@journal#1{%
 \def\@tempa{#1}%
 \@ifxundefined\@journal{%
  \class@info{Selecting journal \@tempa}%
  \let\@journal\@tempa
 }{%
  \@ifx{\@tempa\@journal}{}{%
   \class@warn{Conflicting journal \@tempa<>\@journal; not selected}%
  }%
 }%
}%
\DeclareOption{osa}{\change@society{osa}\let\@journal\@undefined}%
\DeclareOption{osameet}{\change@society{osa}\def\@journal{osameet}}%
\DeclareOption{opex}{\change@society{osa}\def\@journal{opex}}%
\DeclareOption{tops}{\change@society{osa}\def\@journal{tops}}%
\DeclareOption{josa}{\change@society{osa}\def\@journal{josa}}%
\let\rtx@require@packages\@empty
\appdef\rtx@require@packages{%
 \RequirePackage[overload]{textcase}%
}%
\DeclareOption{amsfonts}{%
  \def\class@amsfonts{\RequirePackage{amsfonts}}%
}%
\DeclareOption{noamsfonts}{%
  \let\class@amsfonts\@empty
}%
\appdef\rtx@require@packages{%
 \@ifxundefined\class@amsfonts{}{\class@amsfonts}%
}%
\DeclareOption{amssymb}{%
  \def\class@amssymb{\RequirePackage{amssymb}}%
}%
\DeclareOption{noamssymb}{%
  \let\class@amssymb\@empty
}%
\appdef\rtx@require@packages{%
 \@ifxundefined\class@amssymb{}{\class@amssymb}%
}%
\DeclareOption{amsmath}{%
  \def\class@amsmath{\RequirePackage{amsmath}[\ver@amsmath@prefer]}%
}%
\DeclareOption{noamsmath}{%
  \let\class@amsmath\@empty
}%
\appdef\rtx@require@packages{%
 \preserve@LaTeX
 \@ifxundefined\class@amsmath{}{\class@amsmath}%
 \appdef\class@enddocumenthook{\test@amsmath@ver}%
}%
\appdef\preserve@LaTeX{%
 \let\@ifl@aded@LaTeX\@ifl@aded
 \let\@ifpackageloaded@LaTeX\@ifpackageloaded
 \let\@pkgextension@LaTeX\@pkgextension
 \let\@ifpackagelater@LaTeX\@ifpackagelater
 \let\@ifl@ter@LaTeX\@ifl@ter
 \let\@ifl@t@r@LaTeX\@ifl@t@r
 \let\@parse@version@LaTeX\@parse@version
}%
\appdef\restore@LaTeX{%
 \let\@ifl@aded\@ifl@aded@LaTeX
 \let\@ifpackageloaded\@ifpackageloaded@LaTeX
 \let\@pkgextension\@pkgextension@LaTeX
 \let\@ifpackagelater\@ifpackagelater@LaTeX
 \let\@ifl@ter\@ifl@ter@LaTeX
 \let\@ifl@t@r\@ifl@t@r@LaTeX
 \let\@parse@version\@parse@version@LaTeX
}%
\def\test@amsmath@ver{%
 \begingroup
  \restore@LaTeX
  \@ifpackageloaded{amsmath}{%
   \@ifpackagelater{amsmath}{\ver@amsmath@prefer}{}{%
     \class@warn{%
      You have loaded amsmath, version "\csname ver@amsmath.sty\endcsname",\MessageBreak
      but this class requires version "\ver@amsmath@prefer", or later.\MessageBreak
      Please update your LaTeX installation.
     }%
   }%
  }{%
  }%
 \endgroup
}%
\def\ver@amsmath@prefer{2000/01/15 v2.05 AMS math features}%
\DeclareOption{byrevtex}{\@booleantrue\byrevtex@sw}%
\@booleanfalse\byrevtex@sw
\DeclareOption{floatfix}{\@booleantrue\force@deferlist@sw}%
\DeclareOption{nofloatfix}{\@booleanfalse\force@deferlist@sw}%
\@booleanfalse\force@deferlist@sw
\gdef\@fltovf{%
 \@latex@error{%
  Too many unprocessed floats%
  \force@deferlist@sw{}{; try class option [floatfix]}%
 }\@ehb
}%
\def\@fltstk{%
 \@latex@warning{%
  A float is stuck (cannot be placed)%
  \force@deferlist@sw{}{; try class option [floatfix]}%
 }%
}%
\DeclareOption{ltxgridinfo}{%
 \@booleantrue\ltxgrid@info@sw
}%
\DeclareOption{outputdebug}{%
 \@booleantrue\outputdebug@sw
 \@booleantrue\ltxgrid@info@sw
 \@booleantrue\ltxgrid@foot@info@sw
 \traceoutput
}%
\DeclareOption{raggedfooter}{\@booleanfalse\textheight@sw}%
\DeclareOption{noraggedfooter}{\@booleantrue\textheight@sw}%
\DeclareOption{frontmatterverbose}{\@booleantrue\frontmatterverbose@sw}%
\@booleanfalse\frontmatterverbose@sw
\DeclareOption{linenumbers}{%
 \appdef
 \class@documenthook{%
  \RequirePackage{lineno}[2005/11/02 v4.41]%
  \linenumbersep4pt\relax
  \linenumbers\relax
 }%
}%
\DeclareOption{nomerge}{%
 \appdef\setup@hook{%
  \@ifnum{\NAT@merge>\z@}{\let\NAT@merge\z@}{}%
 }%
}%
\def\@parse@class@options@society{%
 \edef\@tempa{\@ptionlist{\@currname.\@currext}}%
 \expandafter\@for\expandafter\CurrentOption\expandafter:\expandafter=\@tempa\do{%
  \expandafter\@ifnotrelax\csname ds@\CurrentOption\endcsname{}{%
   \IfFileExists{\CurrentOption\substyle@post.\substyle@ext}{%
    \expandafter\change@society\expandafter{\CurrentOption}%
    \expandafter\let\csname ds@\CurrentOption\endcsname\@empty
   }{}%
  }%
 }%
}%
\def\@parse@class@options@#1{%
 \edef\@tempa{\@ptionlist{\@currname.\@currext}}%
 \expandafter\@for\expandafter\CurrentOption\expandafter:\expandafter=\@tempa\do{%
  \expandafter\@ifnotrelax\csname ds@\CurrentOption\endcsname{%
   \begingroup\csname ds@\CurrentOption\endcsname
    \@ifxundefined#1{%
     \endgroup
    }{%
     \expandafter\endgroup\expandafter\def\expandafter#1\expandafter{#1}%
    }%
  }{}%
 }%
}%
\def\@parse@class@options@journal{%
 \edef\@tempa{\@ptionlist{\@currname.\@currext}}%
 \expandafter\@for\expandafter\CurrentOption\expandafter:\expandafter=\@tempa\do{%
  \expandafter\@ifnotrelax\csname ds@\CurrentOption\endcsname{%
   \begingroup
    \csname ds@\CurrentOption\endcsname
    \@ifxundefined\@journal{%
     \endgroup
    }{%
     \expandafter\endgroup\expandafter\def\expandafter\@journal\expandafter{\@journal}%
    }%
  }{}%
 }%
}%
\def\@parse@class@options{%
 \edef\@tempa{\@ptionlist{\@currname.\@currext}}%
 \expandafter\@for\expandafter\CurrentOption\expandafter:\expandafter=\@tempa\do{%
  \expandafter\@ifnotrelax\csname ds@\CurrentOption\endcsname{%
   \begingroup
    \csname ds@\CurrentOption\endcsname
    \@ifxundefined\@pointsize{%
     \endgroup
    }{%
     \expandafter\endgroup\expandafter\def\expandafter\@pointsize\expandafter{\@pointsize}%
    }%
  }{%
   \IfFileExists{\CurrentOption\substyle@post.\substyle@ext}{%
    \expandafter\change@society\expandafter{\CurrentOption}%
    \expandafter\let\csname ds@\CurrentOption\endcsname\@empty
   }{}%
  }%
 }%
}%
\DeclareOption{hypertext}{\hypertext@enable@ltx}%
\appdef\document@inithook{\@ifpackageloaded{hyperref}{\hypertext@enable@ltx}{}}%
\DeclareOption{frontmatterverbose}{\@booleantrue\frontmatterverbose@sw}%
\@booleanfalse\frontmatterverbose@sw
\DeclareOption{inactive}{\@booleanfalse\frontmatter@syntax@sw}%
\@booleantrue\frontmatter@syntax@sw
\@booleanfalse\runinaddress@sw
\@booleantrue\@affils@sw
\@booleanfalse\groupauthors@sw
\DeclareOption{groupedaddress}{\clo@groupedaddress}%
\def\clo@groupedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{unsortedaddress}{\clo@unsortedaddress}%
\def\clo@unsortedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleanfalse\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{runinaddress}{\clo@runinaddress}%
\def\clo@runinaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleantrue\runinaddress@sw
}%
\DeclareOption{superscriptaddress}{\clo@superscriptaddress}%
\def\clo@superscriptaddress{%
 \@booleanfalse\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
%%%  @LaTeX-file{
%%%     filename        = "revtex4-1.dtx",
%%%     version         = "4.1r",
%%%     date            = "2010/07/25",
%%%     time            = "20:33:00 UT-8",
%%%     checksum        = "7177",
%%%     author          = "Arthur Ogawa (mailto:arthur_ogawa at sbcglobal.net),
%%%                        commissioned by the American Physical Society.
%%%                        ",
%%%     copyright       = "Copyright (C) 1999, 2009 Arthur Ogawa,
%%%                        distributed under the terms of the
%%%                        LaTeX Project Public License, see
%%%                        ftp://ctan.tug.org/macros/latex/base/lppl.txt
%%%                        ",
%%%     address         = "Arthur Ogawa,
%%%                        USA",
%%%     telephone       = "",
%%%     FAX             = "",
%%%     email           = "mailto colon arthur_ogawa at sbcglobal.net",
%%%     codetable       = "ISO/ASCII",
%%%     keywords        = "latex, page grid, main vertical list",
%%%     supported       = "yes",
%%%     abstract        = "package to change page grid, MVL",
%%%     docstring       = "The checksum field above generated by ltxdoc",
%%%  }
\def\substyle@post{4-1}%
\def\substyle@ext{rtx}%
\DeclareOption*{\OptionNotUsed}%
\def\@process@society#1{%
 \@ifxundefined\@society{%
  \class@warn{No Society specified, using default society #1}%
  \def\@society{#1}\let\@journal\@undefined
 }{}%
 \expandafter\input\expandafter{\@society\substyle@post.\substyle@ext}%
}%
\def\@process@journal#1{%
 \@ifxundefined\@journal{%
  \class@warn{No journal specified, using default #1}%
  \def\@journal{#1}%
 }{}%
 \expandafter\expandafter
 \expandafter\rtx@do@substyle
 \expandafter\expandafter
 \expandafter{\expandafter\@society\@journal}%
}%
\def\rtx@do@substyle#1{%
 \InputIfFileExists{#1\substyle@post.\substyle@ext}{}{\csname rtx@#1\endcsname}%
}%
\def\@process@pointsize#1{%
 \@ifxundefined\@pointsize{%
  \def\@pointsize{#1}%
  \class@warn{No type size specified, using default \@pointsize}%
 }{}%
 \expandafter\expandafter
 \expandafter\rtx@do@substyle
 \expandafter\expandafter
 \expandafter{\expandafter\@society\@pointsize pt}%
}%
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\slshape\leftmark}%
      \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\MakeUppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\quad
        \fi
        ##1}}}%
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }%
\def\ps@article{%
    \def\@evenhead{\let\\\heading@cr\thepage\quad\checkindate\hfil{\leftmark}}%
    \def\@oddhead{\let\\\heading@cr{\rightmark}\hfil\checkindate\quad\thepage}%
    \def\@oddfoot{}%
    \def\@evenfoot{}%
    \let\@mkboth\markboth
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\def\ps@article@final{%
    \def\@evenhead{\let\\\heading@cr\thepage\quad\checkindate\hfil{\leftmark}}%
    \def\@oddhead{\let\\\heading@cr{\rightmark}\hfil\checkindate\quad\thepage}%
    \def\@oddfoot{}%
    \def\@evenfoot{}%
    \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth{%
       \MakeTextUppercase{%
        \@ifnum{\c@secnumdepth >\z@}{\thesection\hskip 1em\relax}{}%
         ##1%
       }%
       }{}%
    }%
    \def\subsectionmark##1{%
      \markright {%
        \@ifnum{\c@secnumdepth >\@ne}{\thesubsection\hskip 1em\relax}{}%
         ##1%
      }%
    }%
}%
\def\heading@cr{\unskip\space\ignorespaces}%
\def\ps@preprint{%
  \def\@oddfoot{\hfil\thepage\quad\checkindate\hfil}%
  \def\@evenfoot{\hfil\thepage\quad\checkindate\hfil}%
  \def\@oddhead{}%
  \def\@evenhead{}%
  \let\@mkboth\@gobbletwo
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\let\@oddhead\@empty
\let\@evenhead\@empty
\let\@oddfoot\@empty
\let\@evenfoot\@empty
\def\lastpage@putlabel{%
 \if@filesw
  \begingroup
    \advance\c@page\m@ne
    \immediate\write\@auxout{\string\newlabel{LastPage}{{}{\thepage}{}{}{}}}%
  \endgroup
 \fi
}%
\appdef\clear@document{%
 \do@output@cclv{%
  \lastpage@putlabel
  \tally@box@size@sw{\total@text}{}%
 }%
}%
\providecommand\write@column@totals{}%
\appdef\rtx@require@packages{%
 \RequirePackage{url}%
}%
\appdef\document@inithook{%
 \incompatible@package{cite}%
 \incompatible@package{mcite}%
 \incompatible@package{multicol}%
}%
\def\labelenumi{\theenumi.}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{(\theenumii)}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv.}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}
\def\labelitemi{\textbullet}
\def\labelitemii{\normalfont\bfseries\textendash}
\def\labelitemiii{\textasteriskcentered}
\def\labelitemiv{\textperiodcentered}
\pagenumbering{arabic}
\setcounter{topnumber}{2}
\def\topfraction{.9}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.9}
\setcounter{totalnumber}{3}
\def\textfraction{.1}
\def\floatpagefraction{.9}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.9}
\def\dblfloatpagefraction{.9}
\newenvironment{verse}{%
  \let\\=\@centercr
  \list{}{%
    \itemsep\z@ \itemindent -1.5em\listparindent \itemindent
    \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]%
}{%
 \endlist
}%
\newenvironment{quotation}{%
  \list{}{%
    \listparindent 1.5em
    \itemindent\listparindent
    \rightmargin\leftmargin \parsep \z@ \@plus\p@}\item[]%
}{%
 \endlist
}%
\newenvironment{quote}{%
  \list{}{%
  \rightmargin\leftmargin}\item[]%
}{%
 \endlist
}%
\def\descriptionlabel#1{%
  \hspace\labelsep \normalfont\bfseries #1\unskip:%
}%
\newenvironment{description}{%
 \list{}{%
    \labelwidth\z@ \itemindent-\leftmargin
    \let\makelabel\descriptionlabel
 }%
}{%
 \endlist
}%
\newcounter{part}%
\let\thepart\@undefined
\newcounter{section}%
\let\thesection\@undefined
\newcounter{subsection}[section]%
\let\thesubsection\@undefined
\newcounter{subsubsection}[subsection]%
\let\thesubsubsection\@undefined
\newcounter{paragraph}[subsubsection]%
\let\theparagraph\@undefined
\newcounter{subparagraph}[paragraph]%
\let\thesubparagraph\@undefined
\def\secnums@rtx{%
 \@ifxundefined\thepart{%
  \def\thepart{\Roman{part}}%
 }{}%
 \@ifxundefined\thesection{%
  \def\thesection       {\Roman{section}}%
  \def\p@section        {}%
 }{}%
 \@ifxundefined\thesubsection{%
  \def\thesubsection    {\Alph{subsection}}%
  \def\p@subsection     {\thesection\,}%
 }{}%
 \@ifxundefined\thesubsubsection{%
  \def\thesubsubsection {\arabic{subsubsection}}%
  \def\p@subsubsection  {\thesection\,\thesubsection\,}%
 }{}%
 \@ifxundefined\theparagraph{%
  \def\theparagraph     {\alph{paragraph}}%
  \def\p@paragraph      {\thesection\,\thesubsection\,\thesubsubsection\,}%
 }{}%
 \@ifxundefined\thesubparagraph{%
  \def\thesubparagraph  {\arabic{subparagraph}}%
  \def\p@subparagraph   {\thesection\,\thesubsection\,\thesubsubsection\,\theparagraph\,}%
 }{}%
}%
\def\secnums@arabic{%
 \@ifxundefined\thepart{%
  \def\thepart          {\Roman{part}}%
 }{}%
 \@ifxundefined\thesection{%
  \def\thesection       {\Roman{section}}%
  \def\p@section        {}%
 }{}%
 \@ifxundefined\thesubsection{%
  \def\thesubsection    {\thesection.\arabic{subsection}}%
  \def\p@subsection     {}%
 }{}%
 \@ifxundefined\thesubsubsection{%
  \def\thesubsubsection {\thesubsection.\arabic{subsubsection}}%
  \def\p@subsubsection  {}%
 }{}%
 \@ifxundefined\theparagraph{%
  \def\theparagraph     {\thesubsubsection.\arabic{paragraph}}%
  \def\p@paragraph      {}%
 }{}%
 \@ifxundefined\thesubparagraph{%
  \def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}%
  \def\p@subparagraph   {}%
 }{}%
}%
\newenvironment{acknowledgments}{%
 \acknowledgments@sw{%
  \expandafter\section\expandafter*\expandafter{\acknowledgmentsname}%
 }{%
  \par
  \phantomsection
  \addcontentsline{toc}{section}{\protect\numberline{}\acknowledgmentsname}%
 }%
}{%
 \par
}%
\@booleantrue\acknowledgments@sw
\newenvironment{acknowledgements}{%
 \replace@environment{acknowledgements}{acknowledgments}%
}{%
 \endacknowledgments
}%
\def\part{\par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}%
\def\@part[#1]#2{%
 \@ifnum{\c@secnumdepth >\m@ne}{%
        \refstepcounter{part}%
        \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
 }{%
      \addcontentsline{toc}{part}{#1}%
 }%
 \begingroup
    \parindent \z@ \raggedright
    \interlinepenalty\@M
    \@ifnum{\c@secnumdepth >\m@ne}{%
      \Large \bf \partname~\thepart%
      \par\nobreak
    }{}%
    \huge \bf
    #2%
    \markboth{}{}\par
 \endgroup
   \nobreak
   \vskip 3ex
   \@afterheading
}%
\def\@spart#1{{\parindent \z@ \raggedright
    \interlinepenalty\@M
    \huge \bf
    #1\par}
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\section{%
  \@startsection
    {section}%
    {1}%
    {\z@}%
    {0.8cm \@plus1ex \@minus .2ex}%
    {0.5cm}%
    {\normalfont\small\bfseries}%
}%
\def\subsection{%
  \@startsection
    {subsection}%
    {2}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm}%
    {\normalfont\small\bfseries}%
}%
\def\subsubsection{%
  \@startsection
    {subsubsection}%
    {3}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm}%
    {\normalfont\small\itshape}%
}%
\def\paragraph{%
  \@startsection
    {paragraph}%
    {4}%
    {\parindent}%
    {\z@}%
    {-1em}%
    {\normalfont\normalsize\itshape}%
}%
\def\subparagraph{%
  \@startsection
    {subparagraph}%
    {5}%
    {\parindent}%
    {3.25ex \@plus1ex \@minus .2ex}%
    {-1em}%
    {\normalfont\normalsize\bfseries}%
}%
\def\theequation{%
 \theequation@prefix\arabic{equation}%
}%
\def\theequation@prefix{}%
\setcounter{secnumdepth}{4}
\lineskip 1pt
\normallineskip 1pt
\def\baselinestretch{1}%
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\arraycolsep 3pt
\tabcolsep 2pt
\arrayrulewidth .4pt
\doublerulesep 2pt
\skip\@mpfootins = 0pt
\fboxsep  = 3.0pt
\fboxrule = 0.4pt
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\def\listoffigures{\print@toc{lof}}%
\def\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{2\p@}
\long\def\@makecaption#1#2{%
  \par
  \vskip\abovecaptionskip
  \begingroup
   \small\rmfamily
   \sbox\@tempboxa{%
    \let\\\heading@cr
    \@make@capt@title{#1}{#2}%
   }%
   \@ifdim{\wd\@tempboxa >\hsize}{%
    \begingroup
     \samepage
     \flushing
     \let\footnote\@footnotemark@gobble
     \@make@capt@title{#1}{#2}\par
    \endgroup
   }{%
     \global \@minipagefalse
     \hb@xt@\hsize{\hfil\unhbox\@tempboxa\hfil}%
   }%
  \endgroup
  \vskip\belowcaptionskip
}%
\def\@make@capt@title#1#2{%
 \@ifx@empty\float@link{\@firstofone}{\expandafter\href\expandafter{\float@link}}%
  {#1}\@caption@fignum@sep#2%
}%
\def\@footnotemark@gobble{%
 \@footnotemark
 \@ifnextchar[{\@gobble@opt@i}{\@gobble}%
}%
\def\@gobble@opt@i[#1]#2{}%
\def\@mpmakefntext#1{%
 \flushing
 \parindent=1em
 \noindent
 \hb@xt@1em{\hss\@makefnmark}%
 #1%
}%
\def\@caption@fignum@sep{. }%
\def\setfloatlink{\def\float@link}%
\let\float@link\@empty
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\expandafter\newbox\csname fbox@\ftype@figure\endcsname
\expandafter\setbox\csname fbox@\ftype@figure\endcsname\hbox{}%
\appdef\class@documenthook{%
 \do@if@floats{figure}{.fgx}%
}%
\appdef\class@enddocumenthook{%
 \printfigures\relax
}%
\newcommand\printfigures{%
 \@ifstar{\true@sw}{\floatp@sw{\true@sw}{\false@sw}}%
 {%
  \print@float{figure}{\oneapage}%
 }{%
  \print@float{figure}{}%
 }%
}%
\appdef\@xfloat@prep{%
 \appdef\@parboxrestore{\centering}%
}%
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@Roman\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\expandafter\newbox\csname fbox@\ftype@table\endcsname
\expandafter\setbox\csname fbox@\ftype@table\endcsname\hbox{}%
\def\listoftables{\print@toc{lot}}%
\let\l@table\l@figure
\def\table@hook{\small}%
\def\squeezetable{\def\table@hook{\scriptsize}}%
\appdef\@floatboxreset{\table@hook}%
\def\set@table@environments{%
 \floats@sw{}{%
  \let@environment{longtable@float}{longtable}%
  \let@environment{longtable}{longtable@write}%
  \let@environment{longtable*@float}{longtable*}%
  \let@environment{longtable*}{longtable*@write}%
  \let@environment{turnpage@float}{turnpage}%
  \let@environment{turnpage}{turnpage@write}%
 }%
 \do@if@floats{table}{.tbx}%
}%
\appdef\document@inithook{%
 \set@table@environments
}%
\appdef\class@enddocumenthook{%
 \printtables\relax
}%
\newenvironment{longtable@write}{%
 \write@@float{longtable}{table}%
}{%
 \endwrite@float
}%
\newenvironment{longtable*@write}{%
 \write@@float{longtable*}{table}%
}{%
 \endwrite@float
}%
\newenvironment{turnpage@write}{%
 \immediate\write\tablewrite{\string\begin{turnpage}}%
}{%
 \immediate\write\tablewrite{\string\end{turnpage}}%
}%
\newcommand\printtables{%
 \begingroup
  \let@environment{longtable}{longtable@float}%
  \let@environment{longtable*}{longtable*@float}%
  \let@environment{turnpage}{turnpage@anchored}%
  \prepdef\longtable{\trigger@float@par}%
  \expandafter\prepdef\csname longtable*\endcsname{\trigger@float@par}%
  \expandafter\prepdef\csname table@floats\endcsname{%
   \onecolumngrid@push
  }%
  \expandafter\appdef\csname endtable@floats\endcsname{%
   \onecolumngrid@pop
  }%
  \@ifstar{\true@sw}{\floatp@sw{\true@sw}{\false@sw}}%
  {%
   \print@float{table}{\oneapage}%
  }{%
   \print@float{table}{}%
  }%
 \endgroup
}%
\newenvironment{turnpage@anchored}{%
 \onecolumngrid@push
 \setbox\z@\vbox to\textwidth\bgroup
  \columnwidth\textheight
}{%
  \vfil
 \egroup
 \rotatebox{90}{\box\z@}%
 \onecolumngrid@pop
}%
\newenvironment{video}
 {\@float{video}}
 {\end@float}%
\newenvironment{video*}
 {\@dblfloat{video}}
 {\end@dblfloat}%
\newcounter{video}
\renewcommand \thevideo {\@arabic\c@video}
\def\ext@video{lov}%
\def\fname@video{Video}%
\def\lovname{List of Videos}%
\def\fps@video{tbp}%
\def\ftype@video{4}%
\def\fnum@video{\fname@video~\thevideo}%
\appdef\document@inithook{%
 \@ifxundefined\c@float@type{}{%
  \global\setcounter{float@type}{8}%
 }%
}%
\expandafter\newbox\csname fbox@\ftype@video\endcsname
\expandafter\setbox\csname fbox@\ftype@video\endcsname\hbox{}%
\let\theHvideo\thevideo
\def\listofvideos{\print@toc{lov}}%
\let\l@video\l@figure
\appdef\class@documenthook{%
 \do@if@floats{video}{.vdx}%
}%
\appdef\class@enddocumenthook{%
 \printvideos\relax
}%
\newcommand\printvideos{%
 \@ifstar{\true@sw}{\floatp@sw{\true@sw}{\false@sw}}%
 {%
  \print@float{video}{\oneapage}%
 }{%
  \print@float{video}{}%
 }%
}%
\def\endtabular@hook{}%
\appdef\document@inithook{%
 \@ifpackageloaded{dcolumn}{%
  \expandafter\@ifnotrelax\csname NC@find@d\endcsname{}{%
   \newcolumntype{d}{D{.}{.}{-1}}%
  }%
 }{}%
}%
\def\toprule{\hline\hline}%
\def\colrule{\hline}%
\def\botrule{\hline\hline}%
\newenvironment{ruledtabular}{%
 \def\array@default{v}%
 \appdef\tabular@hook{\def\@halignto{to\hsize}}%
 \let\tableft@skip@default\tableft@skip
 \let\tableft@skip\tableft@skip@float
 \let\tabmid@skip@default\tabmid@skip
 \let\tabmid@skip\tabmid@skip@float
 \let\tabright@skip@default\tabright@skip
 \let\tabright@skip\tabright@skip@float
 \let\array@row@pre@default\array@row@pre
 \let\array@row@pre\array@row@pre@float
 \let\array@row@pst@default\array@row@pst
 \let\array@row@pst\array@row@pst@float
 \appdef\array@row@rst{%
  \let\array@row@pre\array@row@pre@default
  \let\array@row@pst\array@row@pst@default
  \let\tableft@skip\tableft@skip@default
  \let\tabmid@skip\tabmid@skip@default
  \let\tabright@skip\tabright@skip@default
  \appdef\tabular@hook{\let\@halignto\@empty}%
 }%
}{%
}%
\def\@makefntext#1{%
  \def\baselinestretch{1}%
  \parindent1em%
  \noindent
  \hb@xt@1.8em{%
   \hss\@makefnmark
  }%
  #1%
  \par
}%
\def\@makefnmark{%
 \hbox{%
  \@textsuperscript{%
   \normalfont\@thefnmark
  }%
 }%
}%
\expandafter\DeclareRobustCommand
\expandafter\rev@citet
\expandafter{%
 \expandafter\begingroup
  \expandafter\rtx@swap@citea
  \expandafter\g@bblefirsttoken
              \csname citet \endcsname
}%
\expandafter\DeclareRobustCommand
\expandafter\rev@citealp
\expandafter{%
 \expandafter\begingroup
  \expandafter\rtx@swap@citea
  \expandafter\g@bblefirsttoken
              \csname citealp \endcsname
}%
\expandafter\DeclareRobustCommand
\expandafter\rev@citealpnum
\expandafter{%
 \expandafter\begingroup
  \expandafter\rtx@swap@citenum
  \expandafter\g@bblefirsttoken
              \csname citealp \endcsname
}%
\def\rtx@swap@citenum{%
  \rtx@swap@citea
  \let\@cite\NAT@citenum
  \let\NAT@mbox\mbox
  \let\citeyear\NAT@citeyear
  \let\NAT@space\NAT@spacechar
}%
\def\g@bblefirsttoken{%
 \expandafter\true@sw
 \expandafter\@empty
}%
\newcommand\rtx@citesuper[3]{%
 \ifNAT@swa
  \leavevmode
  \unskip
  \textsuperscript{\normalfont#1}%
  \if*#3*\else\ (#3)\fi
 \else
  #1%
 \fi
 \endgroup
}%
\def\@makefnmark@cite{\begingroup\NAT@swatrue\@cite{{\@thefnmark}}{}{}}%
\def\rtx@bibsection{%
 \@ifx@empty\refname{%
  \par
 }{%
  \expandafter\section\expandafter*\expandafter{\refname}%
  \@nobreaktrue
 }%
}%
\def\rtx@swap@citea{%
 \let\NAT@def@citea\rtx@def@citea
 \let\NAT@def@citea@close\rtx@def@citea@close
 \let\NAT@def@citea@box\rtx@def@citea@box
}%
\def\rtx@def@citea{%
 \def\@citea{\NAT@separator\NAT@space}%
 \advance\c@NAT@ctr\@ne
 \@ifnum{\count@>\tw@}{%
  \@ifnum{\c@NAT@ctr=\count@}{\appdef\@citea{\NAT@conj\NAT@space}}{}%
 }{%
  \def\@citea{\NAT@space\NAT@conj\NAT@space}%
 }%
}%
\def\rtx@def@citea@close{%
 \rtx@def@citea
 \prepdef\@citea{\NAT@@close}%
}%
\def\rtx@def@citea@box{%
 \rtx@def@citea@close
 \expandafter\def\expandafter\@citea\expandafter{\expandafter\NAT@mbox\expandafter{\@citea}}%
}%
\def\NAT@conj{and}%
\def\NAT@BibitemShut#1{%
 \def\@bibstop{#1}%
 \let\bibitem@Stop\bibitemStop
 \let\bibitem@NoStop\bibitemNoStop
 \@ifx{\bibitemShut\relax}{\let\@bibitemShut\@empty}{%
  \expandafter\def\expandafter\@bibitemShut\expandafter{\bibitemShut}%
 }%
}%
\def\BibitemShut@ltx#1{%
 \unskip
 \def\@bibstop{#1}%
 \let\bibitem@Stop\bibitemStop
 \let\bibitem@NoStop\bibitemNoStop
 \@ifx{\bibitemShut\relax}{\let\@bibitemShut\@empty}{%
  \expandafter\def\expandafter\@bibitemShut\expandafter{\bibitemShut}%
 }%
}%
\newenvironment{thebibliography}{}{}%
\let\@listi\@empty
\appdef\rtx@require@packages{%
 \RequirePackage[sort&compress]{natbib}[2009/11/07 8.31a (PWD, AO)]%
 \let@environment{NAT@thebibliography}{thebibliography}%
 \let@environment{thebibliography}{rtx@thebibliography}%
 \let\bibliographystyle@latex\bibliographystyle
 \let\NAT@citesuper\rtx@citesuper
\let\bibsection\rtx@bibsection
\let\NATx@bibsetnum\NAT@bibsetnum
\def\NAT@bibsetnum#1{%
 \setlength{\topsep}{\z@}%
 \NATx@bibsetnum{\ref{LastBibItem}}%
}%
\let\NATx@bibsetup\NAT@bibsetup
\def\NAT@bibsetup{%
 \setlength{\labelwidth}{\z@}%
 \setlength{\labelsep}{\z@}%
 \setlength{\itemindent}{\z@}%
 \setlength{\listparindent}{\z@}%
 \setlength{\topsep}{\z@}%
 \setlength{\parsep}{\z@}%
 \NATx@bibsetup
}%
\let\bibpreamble\@empty
\def\newblock{\ }%
\let\NATx@bibnumfmt\bibnumfmt
\def\bibnumfmt{\place@bibnumber}%
\let\NAT@merge\thr@@
\let\NAT@citeyear\citeyear
\let\onlinecite\rev@citealp
\let\textcite\rev@citet
\@ifx{\BibitemShut\NAT@BibitemShut}{%
 \class@info{Repairing natbib's \string\BibitemShut}%
 \let\BibitemShut\BibitemShut@ltx
}{}%
\let\bibliographystyle@latex\bibliographystyle
\def\bibliographystyle{\@booleantrue\bibliographystyle@sw\def\@bibstyle}%
\@booleanfalse\bibliographystyle@sw
\def\NAT@bibitem@cont{%
 \let\bibitem@Stop\bibitemContinue@Stop
 \let\bibitem@NoStop\bibitemContinue
}%
\def\bibitemNoStop{%
 \@ifx@empty\@bibitemShut{.\spacefactor\@mmm\space}{\@bibitemShut}%
}%
\def\bibitemContinue{%
 \@ifx@empty\@bibitemShut{;\spacefactor\@mmm\space}{\@bibitemShut}%
}%
\def\bibitemContinue@Stop{%
 \@ifx@empty\@bibitemShut{\spacefactor\@mmm\space}{\@bibitemShut}%
}%
}%
\DeclareRobustCommand\onlinecite{\@onlinecite}%
\DeclareRobustCommand\textcite{\@textcite}%
\let\bibliography@latex\bibliography
\def\bibliography#1{%
 \auto@bib@empty
 \begingroup
  \let\auto@bib@innerbib\@empty
  \@ifx@empty{\pre@bibdata}{%
   \bibliography@latex{#1}%
  }{%
   \@if@empty{#1}{%
    \expandafter\bibliography@latex\expandafter{\pre@bibdata}%
   }{%
    \expandafter\bibliography@latex\expandafter{\pre@bibdata,#1}%
   }%
  }%
 \endgroup
}%
\let\pre@bibdata\@empty
\newenvironment{rtx@thebibliography}[1]{%
 \NAT@thebibliography{#1}%
 \let\@TBN@opr\present@bibnote
 \@FMN@list
}{%
 \auto@bib@innerbib
 \edef\@currentlabel{\arabic{NAT@ctr}}%
 \label{LastBibItem}%
 \endNAT@thebibliography
 \aftergroup\auto@bib@empty
}%
\def\present@bibnote#1#2{%
 \item[%
  \textsuperscript{%
   \normalfont
   \Hy@raisedlink{\hyper@anchorstart{frontmatter.#1}\hyper@anchorend}%
   \begingroup
    \csname c@\@mpfn\endcsname#1\relax
    \frontmatter@thefootnote
   \endgroup
  }%
 ]#2\par
}%
\def\write@bibliographystyle{%
 \@ifxundefined\@bibstyle{}{%
  \expandafter\bibliographystyle@latex\expandafter{\@bibstyle}%
  \bibliographystyle@sw{}{\@bibdataout@rev}%
 }%
 \global\let\write@bibliographystyle\relax
}%
\AtEndDocument{\write@bibliographystyle}%
\def\rtx@@citetp[#1]{\@ifnextchar[{\rtx@citex[#1]}{\rtx@citex[][#1]}}%
\def\rtx@citex[#1][#2]#3{%
 \begingroup
  \def\@tempa{[#1][#2]{#3}}%
  \@ifx{\@cite\NAT@citesuper}{%
   \leavevmode
   \skip@\lastskip
   \unskip
   \super@cite@let
  }{%
   \super@cite@end
  }%
}%
\def\super@cite@let{%
 \futurelet\@let@token\super@cite@check
}%
\def\super@cite@end{%
 \aftergroup\@citex\expandafter\endgroup\@tempa
}%
\def\super@cite@check{%
 \@ifx{\@let@token\@sptoken}{%
  \super@cite@end
 }{%
  \super@cite@swap
 }%
}%
\long\def\super@cite@swap#1{%
 \expandafter\@ifx\expandafter{\csname rtx@automove#1\endcsname\@empty}{%
  #1%
  \super@cite@let
 }{%
  \super@cite@end
  #1%
 }%
}%
\expandafter\let\csname rtx@automove.\endcsname\@empty
\expandafter\let\csname rtx@automove,\endcsname\@empty
\expandafter\let\csname rtx@automove:\endcsname\@empty
\expandafter\let\csname rtx@automove;\endcsname\@empty
\appdef\class@documenthook{%
 \citeautoscript@sw{%
  \@ifx{\@cite\NAT@citesuper}{%
   \let\NAT@@citetp\rtx@@citetp
  }{}%
 }{}%
}%
\def\mini@note{\save@note\mini@notes}%Implicit #2
\def\save@note#1#2{%
  \stepcounter\@mpfn
  \protected@xdef\@thefnmark{\thempfn}%
  \@footnotemark
  \expandafter\g@addto@macro
  \expandafter#1%
  \expandafter{%
  \expandafter \@@footnotetext
  \expandafter {\@thefnmark}{#2}%
              }%
}%
\long\def\@@footnotetext#1{\def\@thefnmark{#1}\@footnotetext}%
\let\mini@notes\@empty
\def\rev@citemark#1{%
 \expandafter\cite\expandafter{\@thefnmark}%
}%
\def\rev@endtext#1{%
 \let\@endnotelabel\@thefnmark
 \@endnotetext
}%
\def\endnote@ext{.end}%
\def\bibdata@app{Notes}%
\def\bibdata@ext{bib}%
\long\def\@endnotetext#1{%
  \begingroup
    \endnote@relax
    \immediate\write\@bibdataout{%
     @FOOTNOTE{%
      \@endnotelabel,%
      key="\@endnotelabel",%
      note="#1"%
     }%
    }%
  \endgroup
}%
\newwrite\@bibdataout
\def\endnote@relax{%
 \let\label\relax \let\index\relax \let\glossary\relax
 \let\cite \relax \let\ref  \relax \let\pageref \relax
 \let\(    \relax \let\)    \relax \let\\       \relax
 \let~\relax
 \let \protect \@unexpandable@protect
 \newlinechar`\^^M%
 \let\begin\relax \let\end\relax
}%
\appdef\class@documenthook{\@bibdataout@init}%
\def\@bibdataout@init{%
 \immediate\openout\@bibdataout\pre@bibdata.\bibdata@ext\relax
}%
\def\@bibdataout@rev{%
 \immediate\write\@bibdataout{%
  @CONTROL{%
   REVTEX41Control%
   \eprint@enable@sw{}{,eprint=""}%
  }%
 }%
 \if@filesw
  \immediate\write\@auxout{\string\citation{REVTEX41Control}}%
 \fi
}%
\def\printendnotes{%
 \class@warn{The \string\printendnotes\space command no longer serves any function. Please remove it from your document.}%
}%
\def\make@footnote@endnote{%
 \footinbib@sw{%
  \authoryear@sw{}{%
   \ltx@footnote@push
   \def\thempfn{Note\thefootnote}%
   \let\ltx@footmark\rev@citemark
   \let\ltx@foottext\rev@endtext
  }%
 }{}%
}%
\def\ltx@footnote@push{%
 \let\ltx@footmark@latex\ltx@footmark
 \let\ltx@foottext@latex\ltx@foottext
 \let\thempfn@latex\thempfn
 \def\ltx@footnote@pop{%
  \let\ltx@footmark\ltx@footmark@latex
  \let\ltx@foottext\ltx@foottext@latex
  \let\thempfn\thempfn@latex
 }%
}%
\appdef\class@documenthook{%
 \make@footnote@endnote
}%
\def\auto@bib{%
 \@ifx@empty\@FMN@list{%
  \footinbib@sw{%
   \@ifnum{\csname c@\@mpfn\endcsname>\z@}{%
    \true@sw
   }{%
    \test@bbl@sw
   }%
  }{%
   \test@bbl@sw
  }%
 }{%
  \true@sw
 }%
 {%
  \bibliography{}%
 }{}%
}%
\def\auto@bib@empty{%
 \let\auto@bib\@empty
}%
\def\test@bbl@sw{%
 \setbox\z@\vbox\bgroup
  \let\providecommand\providecommand@j@nk
  \let\bibfield\@gobbletwo
  \let\bibinfo\@gobbletwo
  \let\translation\@gobble
  \let\BibitemOpen\@empty
  \let\bibitemStop\@empty
  \let\bibitemNoStop\@empty
  \let\EOS\@empty
  \let\BibitemShut\@gobble
  \let\bibAnnoteFile\@gobbletwo
  \let\bibAnnote\@gobblethree
  \let\textbf\@gobble
  \let\emph\@gobble
  \@booleanfalse\bibitem@sw
  \let\bibitem\bibitem@set
  \auto@bib@innerbib
  \bibitem@sw{\aftergroup\true@sw}{\aftergroup\false@sw}%
 \egroup
}%
\newcommand\bibitem@set[1][]{%
 \bibitem@sw{}{%
  \@booleantrue\bibitem@sw
  \aftergroup\@booleantrue\aftergroup\bibitem@sw
 }%
}%
\def\auto@bib@innerbib{%
 \begingroup
  \let@environment{thebibliography}{thebibliography@nogroup}%
  \bibliography{}%
 \endgroup
}%
\def\thebibliography@nogroup#1{%
 \endgroup
 \def\@currenvir{thebibliography}%
}%
\def\endthebibliography@nogroup{\begingroup}%
\long\def \@gobblethree #1#2#3{}%
\def\providecommand@j@nk#1[#2]{%
 \@ifnum{#2=\z@}{\def\j@nk}{%
  \@ifnum{#2=\@ne}{\def\j@nk##1}{%
   \@ifnum{#2=\tw@}{\def\j@nk##1##2}{%
    \@ifnum{#2=\thr@@}{\def\j@nk##1##2##3}{%
    }%
   }%
  }%
 }%
}%
\def\thepage{\@arabic\c@page}%
\appdef\setup@hook{%
 \tabbingsep \labelsep
 \leftmargin\leftmargini
 \labelwidth\leftmargin\advance\labelwidth-\labelsep
 \let\@listi\@listI
 \@listi
}%
\appdef\class@documenthook{%
 \global\c@page\@ne
 \def\curr@envir{document}%
 \mark@envir{\curr@envir}%
}%
\def\open@onecolumn{%
 \open@column@one\@ne
 \set@colht
 \@floatplacement
 \@dblfloatplacement
}%
\def\open@twocolumn{%
 \open@column@mlt\tw@
 \set@colht
 \@floatplacement
 \@dblfloatplacement
 \sloppy
 \let\set@listindent\set@listindent@
}%
\def\appendix{%
 \par
 \setcounter{section}\z@
 \setcounter{subsection}\z@
 \setcounter{subsubsection}\z@
 \def\thesubsection{\arabic{subsection}}%
 \def\thesubsubsection{\alph{subsubsection}}%
 \@addtoreset{equation}{section}%
 \def\theequation@prefix{\thesection}%
 \addtocontents{toc}{\protect\appendix}%
 \@ifstar{%
  \def\thesection{\unskip}%
  \def\theequation@prefix{A.}%
 }{%
  \def\thesection{\Alph{section}}%
 }%
}%
\def\title@column#1{%
 \minipagefootnote@init
 #1%
 \minipagefootnote@foot
}%
\def\close@column{%
 \newpage
}%
\def\galley@outdent{\rightmargin-\columnwidth\advance\rightmargin-\columnsep}%
\let\widetext@outdent\@empty
\newenvironment{widetext@galley}{%
  \list{}{%
    \topsep        \z@skip
    \listparindent \parindent
    \itemindent    \parindent
    \leftmargin    \z@
    \parsep        \z@\@plus\p@
    \widetext@outdent
    \relax
  }%
  \item\relax
}{
  \endlist
}%
\def\title@column@grid#1{%
 \minipagefootnote@init
  \onecolumngrid
  \begingroup
   \let\@footnotetext\frontmatter@footnotetext
   \ltx@no@footnote
   #1%
  \endgroup
  \twocolumngrid
 \minipagefootnote@foot
}%
\def\close@column@grid{%
 \balancelastpage@sw{%
  \onecolumngrid
 }{}%
}%
\newenvironment{widetext@grid}{%
  \par\ignorespaces
  \setbox\widetext@top\vbox{%
   \hb@xt@\hsize{%
    \leaders\hrule\hfil
    \vrule\@height6\p@
   }%
  }%
  \setbox\widetext@bot\hb@xt@\hsize{%
    \vrule\@depth6\p@
    \leaders\hrule\hfil
  }%
  \onecolumngrid
  \vskip10\p@
  \dimen@\ht\widetext@top\advance\dimen@\dp\widetext@top
  \cleaders\box\widetext@top\vskip\dimen@
  \vskip6\p@
  \prep@math@patch
}{%
  \par
  \vskip6\p@
  \setbox\widetext@bot\vbox{%
   \hb@xt@\hsize{\hfil\box\widetext@bot}%
  }%
  \dimen@\ht\widetext@bot\advance\dimen@\dp\widetext@bot
  \cleaders\box\widetext@bot\vskip\dimen@
  \vskip8.5\p@
  \twocolumngrid\global\@ignoretrue
  \@endpetrue
}%
\newbox\widetext@top
\newbox\widetext@bot
\def\set@page@grid{%
 \twocolumn@sw{%
  \let\set@footnotewidth\set@footnotewidth@two
  \let\compose@footnotes\compose@footnotes@two
  \let@environment{widetext}{widetext@grid}%
  \let\title@column\title@column@grid
  \let\close@column\close@column@grid
 }{%
  \let@environment{widetext}{widetext@galley}%
  \preprintsty@sw{%
  }{%
   \galley@sw{%
    \let\widetext@outdent\galley@outdent
   }{}%
  }%
 }%
}%
\appdef\setup@hook{\set@page@grid}%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\def\partname{Part}
\def\tocname{Contents}
\def\lofname{List of Figures}
\def\lotname{List of Tables}
\def\refname{References}
\def\indexname{Index}
\def\figurename{FIG.}
\def\figuresname{Figures}%
\def\tablename{TABLE}
\def\tablesname{Tables}%
\def\abstractname{Abstract}
\def\appendixesname{Appendixes}%
\def\appendixname{Appendix}%
\def\acknowledgmentsname{Acknowledgments}
\def\journalname{??}
\def\copyrightname{??}
\def\andname{and}
\def\@pacs@name{PACS numbers: }%
\def\@keys@name{Keywords: }%
\def\@feats@name{Features: }%
\def\ppname{pp}
\def\numbername{number}
\def\volumename{volume}
\def\Dated@name{Dated: }%
\def\Received@name{Received }%
\def\Revised@name{Revised }%
\def\Accepted@name{Accepted }%
\def\Published@name{Published }%
\def\address{\replace@command\address\affiliation}%
\def\altaddress{\replace@command\altaddress\altaffiliation}%
\newenvironment{references}{%
 \class@warn@end{The references environment is not supported; use thebibliography instead.}
 \gdef\references{\thebibliography{}}\references
}{%
 \endthebibliography
}%
\def\draft{%
 \class@warn@end{Command \string\draft\space is obsolete;^^JInvoke option draft instead.}%
 \@booleantrue\draft@sw
}%
\def\tighten{%
 \class@warn@end{Command \string\tighten\space is obsolete;^^JInvoke option tightenlines instead.}%
 \@booleantrue\tightenlines@sw
}%
\def\tableline{%
 \noalign{%
  \class@warn@end{Command \string\tableline\space is obsolete;^^JUse \string\colrule\space instead.}%
  \global\let\tableline\colrule
 }%
 \tableline
}%
\def\case{\replace@command\case\frac}%
\def\slantfrac{\replace@command\slantfrac\frac}%
\def\tablenote{\replace@command\tablenote\footnote}%
\def\tablenotemark{\replace@command\tablenotemark\footnotemark}%
\def\tablenotetext{\replace@command\tablenotetext\footnotetext}%
\DeclareRobustCommand\REV@text[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }#1}}}%
   {\hbox{{\everymath{\textstyle        }#1}}}%
   {\hbox{{\everymath{\scriptstyle      }\let\f@size\sf@size\selectfont#1}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\let\f@size\ssf@size\selectfont#1}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bbox[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }\boldmath$#1$}}}%
   {\hbox{{\everymath{\textstyle        }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptstyle      }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\boldmath$#1$}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bm[1]{%
 \class@warn@end{To use \string\bm, please load the bm package!}%
 \global\let\bm\relax
}%
\def\FL{\obsolete@command\FL}%
\def\FR{\obsolete@command\FR}%
\def\narrowtext{\obsolete@command\narrowtext}%
\def\mediumtext{\obsolete@command\mediumtext}%
\newenvironment{quasitable}{%
 \let@environment{tabular}{longtable}%
}{%
}%
\let\text\REV@text
\let\bm\REV@bm
\appdef\setup@hook{%
 \providecommand\bibinfo[2]{#2}%
 \providecommand\eprint[2][]{#2}%
}%
\def\bbox#1{%
 \class@warn@end{\string\bbox\space is obsolete,^^Jload the bm package and use \string\bm\space instead.}%
 \global\let\bbox\relax
}%
\newenvironment{mathletters}{%
 \class@warn@end{Environment {mathletters} is obsolete;^^Jload the amsmath package and use {subequations}!}%
 \global\let\mathletters\@empty
}{%
}%
\def\eqnum#1{%
 \class@warn@end{\string\eqnum\space is obsolete, load the amsmath package and use \string\tag!}%
 \global\let\eqnum\@gobble
}%
\appdef\rtx@require@packages{%
 \RequirePackage{revsymb4-1}%
}%
\appdef\class@documenthook{\revsymb@inithook}%
\def\@LN@LLerror@org{%
 \PackageError{lineno}{%
  \string\linelabel\space without \string\linenumbers
 }{%
  Just see documentation. (New feature v4.11)%
 }%
 \@gobble
}%
\def\@LN@LLerror@ltx{%
 \PackageWarning{lineno}{%
  To make the \string\linelabel\space command work, you must issue the \string\linenumbers\ command
 }%
 \@gobble
}%
\appdef\class@documenthook{%
 \@ifx{\@LN@LLerror\@LN@LLerror@org}{%
  \class@info{Overriding \string\@LN@LLerror}%
  \let\@LN@LLerror\@LN@LLerror@ltx
 }{}%
 \@ifpackageloaded{lineno}{%
  \@ifxundefined{\set@linepenalties}{}{%
   \def\prep@absbox{\set@linepenalties}%
   \def\post@absbox{\let\@LN@parpgbrk\@empty\@linenumberpar}%
  }%
 }{}%
}%
\appdef\rtx@require@packages{%
 \InputIfFileExists{\jobname.rty}{}{}%
}%
\@parse@class@options@society
\@process@society{aps}%
\@parse@class@options@\@journal
\expandafter\@process@journal\expandafter{\@journal@default}%
\@parse@class@options@\@pointsize
\expandafter\@process@pointsize\expandafter{\@pointsize@default}%
\@options
\rtx@require@packages
\setup@hook
\endinput
%%
%% End of file `revtex4-1.cls'.
